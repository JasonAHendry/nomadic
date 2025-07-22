import time
import webbrowser
from datetime import datetime

from nomadic.download.references import REFERENCE_COLLECTION
from nomadic.realtime.factory import PipelineFactory
from nomadic.util.dirs import (
    ExperimentDirectories,
)
from nomadic.util.logging_config import LoggingFascade
from nomadic.util.metadata import MetadataTableParser
from nomadic.util.regions import RegionBEDParser
from nomadic.util.settings import (
    ExperimentSettings,
    load_settings,
    save_settings,
    verify_compatible_settings,
)

WAIT_INTERVAL = 5


def main(
    expt_name: str,
    output: str,
    workspace: str,
    fastq_dir: str,
    metadata_csv: str,
    region_bed: str,
    reference_name: str,
    call: bool,
    verbose: bool,
) -> None:
    """
    Run nomadic in realtime

    """

    # PARSE INPUT
    log = LoggingFascade(logger_name="nomadic", verbose=verbose)
    log.info("Input parameters:")
    log.info(f"  Experiment Name: {expt_name}")
    log.info(f"  Workspace: {workspace}")
    log.info(f"  Output dir: {output}")
    log.info(f"  FASTQ (.fastq): {fastq_dir}")
    log.info(f"  Metadata (.csv): {metadata_csv}")
    log.info(f"  Regions (.bed): {region_bed}")
    log.info(f"  Reference genome: {reference_name}")
    REFERENCE_COLLECTION[reference_name].confirm_downloaded()
    log.info(f"  Performing variant calling: {call}")
    log.info("Processing...")

    # PREPARE TO RUN
    metadata = MetadataTableParser(metadata_csv)
    regions = RegionBEDParser(region_bed)
    expt_dirs = ExperimentDirectories(output, metadata, regions)
    log.info(f"  Found {len(metadata.barcodes) - 1} barcodes to track.")
    log.info(f"  Found {regions.n_regions} regions of interest.")
    log.info(f"  Outputs will be written to: {expt_dirs.expt_dir}.")
    log.info("Done.\n")

    # LOAD/STORE EXPERIMENT SETTINGS
    previous_settings = load_settings(expt_dirs.get_settings_file())
    experiment_settings = ExperimentSettings(
        name=expt_name,
        start_date=datetime.now().replace(microsecond=0),
        fastq_dir=fastq_dir,
        metadata_csv=metadata_csv,
        region_bed=region_bed,
        reference_name=reference_name,
        n_barcodes=len(metadata.barcodes) - 1,
        n_regions=regions.n_regions,
        call=call,
    )
    start_time = None
    if previous_settings is None:
        # first run
        save_settings(
            expt_dirs.get_settings_file(), experiment_settings=experiment_settings
        )
    else:
        verify_compatible_settings(
            old_settings=previous_settings, new_settings=experiment_settings
        )
        start_time = previous_settings.start_date

    # INITIALISE WATCHERS
    factory = PipelineFactory(
        expt_name, metadata, regions, expt_dirs, fastq_dir, call, reference_name
    )

    watchers = factory.get_watchers()
    expt_pipeline = factory.get_expt_pipeline()
    dashboard = factory.get_dashboard(start_time=start_time)
    dashboard.run(in_thread=True)

    webbrowser.open("http://127.0.0.1:8050")

    # CATCH UP FROM WORK LOG IF WE RESUME
    catch_up_info = [watcher.catch_up_from_work_log() for watcher in watchers]
    processed = sum(info[0] for info in catch_up_info)
    reprocessed = sum(info[1] for info in catch_up_info)
    if processed > 0 or reprocessed > 0:
        log.info("Recovered from work log:")
        log.info(f"{processed} fastq files already processed")
        log.info(
            f"{reprocessed} fastq files reprocessed because they were not finished"
        )
        log.info("Running experiment pipeline...")
        expt_pipeline.run()
        log.info("Resuming pipeline...")

    # BEGIN REALTIME WATCHING
    try:
        while True:
            log.info("Scanning barcodes for new FASTQ files...")
            updated = [watcher.update() for watcher in watchers]
            n_updated = sum(updated)
            if n_updated > 0:
                log.info(f"Have updated {n_updated} barcodes.")
                log.info("Running experiment pipeline...")
                expt_pipeline.run()
            else:
                log.info("No barcodes updated. Waiting before scannning again.")
                time.sleep(WAIT_INTERVAL)
    except KeyboardInterrupt:
        log.info("")
        log.info("Program has been interrupted by user. Exiting.")
        # TODO: can run final status checks here
