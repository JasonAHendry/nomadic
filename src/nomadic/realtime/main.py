import time

from nomadic.util.logging_config import LoggingFascade
from nomadic.util.metadata import MetadataTableParser
from nomadic.util.dirs import ExperimentDirectories
from nomadic.util.regions import RegionBEDParser

from nomadic.realtime.watchers import BarcodeWatcher
from nomadic.realtime.pipelines.barcode import BarcodePipelineRT
from nomadic.realtime.pipelines.experiment import MappingInRTExptPipeline
from nomadic.realtime.dashboard.builders import MappingRTDashboard


WAIT_INTERVAL = 5


def main(expt_name: str, fastq_dir: str, metadata_csv: str, region_bed: str, verbose: bool) -> None:
    """
    Run nomadic in realtime

    """

    # PARSE INPUT
    log = LoggingFascade(logger_name="nomadic", verbose=verbose)
    log.info("Input parameters:")
    log.info(f"  Experiment Name: {expt_name}")
    log.info(f"  FASTQ (.fastq): {fastq_dir}")
    log.info(f"  Metadata (.csv): {metadata_csv}")
    log.info(f"  Regions (.bed): {region_bed}")
    log.info("Processing...")

    # PREPARE TO RUN
    metadata = MetadataTableParser(metadata_csv)
    regions = RegionBEDParser(region_bed)
    expt_dirs = ExperimentDirectories(expt_name, metadata, regions)
    log.info(f"  Found {len(metadata.barcodes) - 1} barcodes to track.")
    log.info(f"  Found {regions.n_regions} regions of interest.")
    log.info(f"  Outputs will be written to: {expt_dirs.expt_dir}.")
    log.info("Done.\n")

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
    dashboard = MappingRTDashboard(
        expt_name=expt_name,
        regions=regions,
        fastq_csv=f"{expt_dirs.approach_dir}/summary.fastq.csv",
        flagstats_csv=f"{expt_dirs.approach_dir}/summary.bam_flagstats.csv",
        bedcov_csv=f"{expt_dirs.approach_dir}/summary.bedcov.csv",
        depth_csv=f"{expt_dirs.approach_dir}/summary.depth.csv"
    )
    dashboard.run(in_thread=True)

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

