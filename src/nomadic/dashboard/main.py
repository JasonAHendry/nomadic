import os
import difflib
from nomadic.realtime.dashboard.builders import MappingRTDashboard, CallingRTDashboard
from nomadic.util.metadata import MetadataTableParser
from nomadic.util.dirs import ExperimentDirectories
from nomadic.util.regions import RegionBEDParser


class ExperimentNotFound(Exception):
    """Could not find experiment"""


def verify_experiment_exists(expt_name: str) -> None:
    """
    When running just the dashboard, we want to verify the expected experiment exists
    Otherwise, we raise exceptions

    """

    result_dir = "results"
    if not os.path.exists(result_dir):
        raise ExperimentNotFound(
            "Cannot find 'results' directory. "
            "Make sure there is the 'results' folder is in your current directory."
        )

    expts_found = os.listdir(result_dir)
    if not expts_found:
        raise ExperimentNotFound("No experiments found in 'results' directory.")

    if not expt_name in expts_found:
        exception_str = (
            f"Found {len(expts_found)} experiments, but none matched '{expt_name}'."
        )
        close_match = difflib.get_close_matches(expt_name, expts_found, cutoff=0.8, n=1)
        if close_match:
            exception_str += f" Did you mean: '{close_match[0]}'?"
        raise ExperimentNotFound(exception_str)


def find_metadata(expt_name: str) -> MetadataTableParser:
    """
    Given an experiment directory, search for the metadata CSV file in thee
    expected location

    """

    metadata_dir = f"results/{expt_name}/metadata"
    csvs = [
        f"{metadata_dir}/{file}"
        for file in os.listdir(metadata_dir)
        if file.endswith(".csv")
    ]  # TODO: what about no-suffix files?

    print(csvs)
    if len(csvs) != 1:  # Could alternatively load and LOOK
        raise FileNotFoundError(
            f"Expected one metadata CSV file (*.csv) at {metadata_dir}, but found {len(csvs)}."
        )

    return MetadataTableParser(csvs[0])


def find_regions(expt_name: str) -> RegionBEDParser:
    """
    Given an experiment directory, search for the metadata CSV file in thee
    expected location

    TODO: Bad duplication from above, can write inner function
    """

    metadata_dir = f"results/{expt_name}/metadata"
    beds = [
        f"{metadata_dir}/{file}"
        for file in os.listdir(metadata_dir)
        if file.endswith(".bed")
    ]  # TODO: what about no-suffix files?

    if len(beds) != 1:  # Could alternatively load and LOOK
        raise FileNotFoundError(
            f"Expected one region BED file (*.bed) at {metadata_dir}, but found {len(beds)}."
        )

    return RegionBEDParser(beds[0])


def variant_calling_performed(expt_dirs: ExperimentDirectories) -> bool:
    """
    Check if the variant calling TSV is present
    """

    return os.path.exists(f"{expt_dirs.approach_dir}/summary.variants.csv")


def main(expt_name):
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

    verify_experiment_exists(expt_name)
    metadata = find_metadata(expt_name)
    expt_dirs = ExperimentDirectories(expt_name, metadata)
    regions = find_regions(expt_name)

    print("Input Parameters:")
    print(f"  Experiment Name: {expt_name}")
    print(f"  Metadata (.csv): {metadata.csv}")
    print(f"  Regions (.bed): {regions.path}")
    print(f"  Found {len(metadata.barcodes) - 1} barcodes in this experiment.")
    print(f"  Found {regions.n_regions} regions of interest.")

    shared_kwargs = {
        "expt_name": expt_name,
        "regions": regions,
        "metadata": metadata,
        "fastq_csv": f"{expt_dirs.approach_dir}/summary.fastq.csv",
        "flagstats_csv": f"{expt_dirs.approach_dir}/summary.bam_flagstats.csv",
        "bedcov_csv": f"{expt_dirs.approach_dir}/summary.bedcov.csv",
        "depth_csv": f"{expt_dirs.approach_dir}/summary.depth.csv",
    }

    if variant_calling_performed(expt_dirs):
        print("  Variant calling: True")
        dashboard = CallingRTDashboard(
            **shared_kwargs,
            variant_csv=f"{expt_dirs.approach_dir}/summary.variants.csv",
        )
    else:
        print("  Variant calling: False")
        dashboard = MappingRTDashboard(**shared_kwargs)
    print("Done.")

    print("")
    print("Launching dashboard (press CNTRL+C to exit):")
    print("")
    dashboard.run(debug=False)
