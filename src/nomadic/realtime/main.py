import time

from nomadic.util.logging_config import LoggingFascade
from nomadic.util.metadata import MetadataTableParser
from nomadic.util.dirs import ExperimentDirectories
from nomadic.util.regions import RegionBEDParser

from nomadic.realtime.factory import PipelineFactory

WAIT_INTERVAL = 5


def main(expt_name: str, fastq_dir: str, metadata_csv: str, region_bed: str, call: bool, verbose: bool) -> None:
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
    log.info(f"  Performing variant calling: {call}")
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
    factory = PipelineFactory(metadata, 
                              regions,
                              expt_dirs,
                              fastq_dir,
                              call)
    
    watchers = factory.get_watchers()
    expt_pipeline = factory.get_expt_pipeline()
    dashboard = factory.get_dashboard()
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

