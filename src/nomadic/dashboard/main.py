import os
from nomadic.realtime.dashboard.builders import MappingRTDashboard, CallingRTDashboard
from nomadic.util.experiment import ExperimentDirectories, find_metadata, find_regions
from nomadic.util.settings import load_settings


def variant_calling_performed(expt_dirs: ExperimentDirectories) -> bool:
    """
    Check if the variant calling TSV is present
    """

    return os.path.exists(expt_dirs.get_summary_files().variants)


def main(input_dir: str):
    """
    Main execution code for just running the dashboard

    TODO:
    Some of the logic for looking for specific directories
    should probably come from ExperimentDirectories;
    Otherwise, I am now defining expected path to metadata
    independently in two places.

    TODO:
    - What happens if the experiment directory doesn't exist?
    - Lot's of other verification might be useful here, like if the
    expected summary files exist

    """

    metadata = find_metadata(input_dir)
    expt_dirs = ExperimentDirectories(input_dir, metadata)
    regions = find_regions(input_dir)
    settings = load_settings(expt_dirs.get_settings_file())

    if settings is not None:
        expt_name = settings.name
        start_time = settings.start_date
    else:
        expt_name = os.path.basename(input_dir)
        start_time = None

    print("Input Parameters:")
    print(f"  Input dir: {input_dir}")
    print(f"  Experiment Name: {expt_name}")
    print(f"  Metadata (.csv): {metadata.csv}")
    print(f"  Regions (.bed): {regions.path}")
    print(f"  Found {len(metadata.barcodes) - 1} barcodes in this experiment.")
    print(f"  Found {regions.n_regions} regions of interest.")

    summary_files = expt_dirs.get_summary_files()

    shared_kwargs = {
        "expt_name": expt_name,
        "regions": regions,
        "metadata": metadata,
        "fastq_csv": summary_files.fastqs_processed,
        "read_mapping_csv": summary_files.read_mapping,
        "region_coverage_csv": summary_files.region_coverage,
        "depth_profiles_csv": summary_files.depth_profiles,
        "start_time": start_time,
    }

    if variant_calling_performed(expt_dirs):
        print("  Variant calling: True")
        dashboard = CallingRTDashboard(
            **shared_kwargs,
            variant_csv=summary_files.variants,
            is_realtime=False,
        )
    else:
        print("  Variant calling: False")
        dashboard = MappingRTDashboard(**shared_kwargs, is_realtime=False)
    print("Done.")

    print("")
    print("Launching dashboard (press CNTRL+C to exit):")
    print("")

    dashboard.run(debug=False)
