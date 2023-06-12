import time
import threading

from nomadic.util.logging_config import LoggingFascade
from nomadic.util.metadata import MetadataTableParser
from nomadic.util.dirs import ExperimentDirectories
from nomadic.realtime.watchers import BarcodeWatcher

# from nomadic.realtime.barcode_pipelines import BasicBarcodeInRT, MappingInRT2
from nomadic.realtime.pipelines.barcode import BarcodePipelineRT
from nomadic.realtime.pipelines.experiment import MappingInRTExptPipeline
from nomadic.realtime.dashboard.creator import create_dashboard


WAIT_INTERVAL = 5


def main(expt_name: str, fastq_dir: str, metadata_csv: str, region_bed: str) -> None:
    """
    Run nomadic in realtime

    """

    # PARSE INPUT
    log = LoggingFascade()
    log.info("Input parameters:")
    log.info(f"  Experiment Name: {expt_name}")
    log.info(f"  FASTQ (.fastq): {fastq_dir}")
    log.info(f"  Metadata (.csv): {metadata_csv}")
    log.info(f"  Regions (.bed): {region_bed}")
    log.info("Done.")

    # PREPARE TO RUN
    metadata = MetadataTableParser(metadata_csv)
    expt_dirs = ExperimentDirectories(expt_name, metadata)
    log.info(f"Found {len(metadata.barcodes)} barcodes to track.")

    # INITIALISE WATCHERS
    watchers = [
        BarcodeWatcher(
            barcode_fastq_dir=f"{fastq_dir}/{b}",
            barcode_pipeline=BarcodePipelineRT(barcode_name=b, expt_dirs=expt_dirs, bed_path=region_bed),
        )
        for b in metadata.barcodes
    ]
    expt_pipeline = MappingInRTExptPipeline(metadata, expt_dirs)

    # Initiliase the dashboard
    dashboard = create_dashboard(
        expt_name=expt_name,
        flagstats_csv=f"{expt_dirs.approach_dir}/summary.bam_flagstats.csv",
        bedcov_csv=f"{expt_dirs.approach_dir}/summary.bedcov.csv"
    )
    dashboard_thread = threading.Thread(target=dashboard.run, name="dashboard")
    dashboard_thread.start()

    # BEGIN REALTIME WATCHING
    try:
        while True:
            log.info("Scanning barcodes for new FASTQ files...")
            updated = [watcher.update() for watcher in watchers]
            n_updated = sum(updated)
            if n_updated > 0:
                log.info(f"Have updated {n_updated} barcodes.")
                log.info(f"Running experiment pipeline...")
                expt_pipeline.run()
            else:
                log.info("No barcodes updated. Waiting before scannning again.")
                time.sleep(WAIT_INTERVAL)
    except KeyboardInterrupt:
        log.info("")
        log.info("Program has been interrupted by user. Exiting.")
        # TODO: can run final status checks here

