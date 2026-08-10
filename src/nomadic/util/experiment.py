import glob
import os
import json
import shutil
import pandas as pd
from pathlib import Path
from typing import Literal, NamedTuple, Any, Optional
from dataclasses import dataclass

from nomadic.util.dirs import produce_dir
from nomadic.util.exceptions import MetadataFormatError
from nomadic.util.metadata import (
    MetadataTableParser,
    ExtendedMetadataTableParser,
    STANDARD_METADATA_FILENAME,
)
from nomadic.util.regions import RegionBEDParser


# --------------------------------------------------------------------------------
# Handle summary file names: legacy and current
#
# --------------------------------------------------------------------------------


class SummaryFiles(NamedTuple):
    """Define summary file names / paths"""

    fastqs_processed: str
    read_mapping: str
    region_coverage: str
    depth_profiles: str
    variants: str


# Currently used summary file names
DEFAULT_CONFIG_PATH = "config/defaults.json"
SUMMARY_NAMES = SummaryFiles(
    fastqs_processed="summary.fastqs_processed.csv",
    read_mapping="summary.read_mapping.csv",
    region_coverage="summary.region_coverage.csv",
    depth_profiles="summary.depth_profiles.csv",
    variants="summary.variants.csv",
)

# Legacy summary file names for backward compatibility
SUMMARY_NAMES_LEGACY = SummaryFiles(
    fastqs_processed="summary.fastq.csv",
    read_mapping="summary.bam_flagstats.csv",
    region_coverage="summary.bedcov.csv",
    depth_profiles="summary.depth.csv",
    variants="summary.variants.csv",
)


def get_summary_files(expt_dir: Path) -> SummaryFiles:
    """
    Determine whether the summary files are use the legacy or current names,
    and return a SummaryFiles object with the appropriate file names

    """

    if not expt_dir.exists():
        raise FileNotFoundError(f"Experiment path does not exist: {expt_dir}")

    if (expt_dir / SUMMARY_NAMES_LEGACY.read_mapping).exists():
        # Detect legacy format using *one* of the differentiating file names
        format_used = SUMMARY_NAMES_LEGACY
    else:
        format_used = SUMMARY_NAMES
    return SummaryFiles(*[str(expt_dir / field) for field in format_used])


# --------------------------------------------------------------------------------
# Define experiment directories
#
# --------------------------------------------------------------------------------


class ExperimentDirectories:
    """
    Put all the information about experimental
    directory structure in a single place

    Advantage is any future changes to how the directories
    are organised can be managed here

    Best practice is that objects/functions *do not* directly depend on
    this object, better to use its members as arguments to other
    objects or functions

    TODO:
    - For specific analysis dirs, can have two dictionaries that grow
      - They are forced to be in particular relation to the existing dirs
      - But can be added to or reduced depending on pipeline

    """

    def __init__(
        self,
        output_dir: str,
        barcodes: list[str],
        regions_basename: str,
    ):
        self.output_dir = output_dir

        self.metadata_dir = os.path.join(self.output_dir, "metadata")
        self.barcodes_dir = os.path.join(self.output_dir, "barcodes")
        self._barcode_dirs = {b: os.path.join(self.barcodes_dir, b) for b in barcodes}
        self.regions_bed = os.path.join(self.metadata_dir, regions_basename)

    def setup_dirs(self):
        """
        Create all the directories for the experiment

        """

        produce_dir(self.output_dir)
        produce_dir(self.metadata_dir)
        produce_dir(self.barcodes_dir)
        for barcode_dir in self._barcode_dirs.values():
            produce_dir(barcode_dir)

    def get_barcode_dir(self, barcode_name: str):
        """
        Get the path to a particular `barcode_name`

        """

        return self._barcode_dirs[barcode_name]

    def get_settings_file(self) -> str:
        """Get the path to the setting file for the experiment"""
        return os.path.join(self.metadata_dir, "settings.json")

    def get_summary_files(self) -> SummaryFiles:
        return get_summary_files(Path(self.output_dir))

    def setup_metadata_dir(
        self, metadata: MetadataTableParser, regions: RegionBEDParser
    ) -> None:
        """
        Move metadata CSV and regions BED into the metadata directory,
        and store their paths as attributes
        """
        if metadata is not None:
            metadata_csv = f"{self.metadata_dir}/{STANDARD_METADATA_FILENAME}"
            if not os.path.exists(metadata_csv):
                metadata.df.to_csv(metadata_csv, index=False)

        if regions is not None:
            if not os.path.exists(self.regions_bed):
                shutil.copy(regions.path, self.regions_bed)


# --------------------------------------------------------------------------------
# Checks on experiment outputs
#
# --------------------------------------------------------------------------------


@dataclass
class ExperimentOutputs:
    """Store information about outputs in `expt_dir`"""

    expt_dir: str  # TODO: change to Path
    metadata: pd.DataFrame
    regions: RegionBEDParser
    summary_files: SummaryFiles
    settings: dict[str, Any]

    # Variant calling outputs
    caller: str
    reference_name: str
    has_complete_vcf: bool
    has_filtered_vcf: bool


def find_metadata(
    expt_dir: str, Parser: MetadataTableParser = MetadataTableParser
) -> MetadataTableParser:
    """
    Given an experiment directory, search for the metadata CSV file in thee
    expected location and load it

    """
    # first check if the file with standard name exists
    standard_path = os.path.join(expt_dir, "metadata", STANDARD_METADATA_FILENAME)
    if os.path.isfile(standard_path):
        return Parser(standard_path)

    # In legacy cases, it should have the name of the experiment
    csv = os.path.join(expt_dir, "metadata", f"{os.path.basename(expt_dir)}.csv")
    if os.path.exists(csv):
        return Parser(csv)

    # finally, look for any CSV file in the metadata directory
    csv = glob.glob(f"{expt_dir}/metadata/*.csv")
    if len(csv) == 1:
        return Parser(csv[0])

    raise ValueError(
        f"Found {len(csv)} *.csv files in '{expt_dir}/metadata', cannot determine which is metadata."
    )


def find_regions(expt_dir: str) -> RegionBEDParser:
    """
    Given an experiment directory, search for the metadata CSV file in thee
    expected location

    """

    bed = [
        f
        for f in glob.glob(f"{expt_dir}/metadata/*.bed")
        if f.endswith(".bed") and not f.endswith(".lowcomplexity_mask.bed")
    ]

    if len(bed) != 1:
        raise FileNotFoundError(
            f"Expected one region BED file (*.bed) at '{expt_dir}/metadata', but found {len(bed)}."
        )

    return RegionBEDParser(bed[0])


def experiment_outputs(
    expt_dir: str,
    allow_missing_files: Optional[
        list[Literal["fastq", "depth", "region", "variant", "read"]]
    ],
) -> ExperimentOutputs:
    """For a given `expt_dir` return the experiment outputs as an `ExperimentOutputs` object.

    Will raise exceptions if data required for summarising is missing.

    """
    if allow_missing_files is None:
        allow_missing_files = []

    # Existence of directory
    if not os.path.isdir(expt_dir):
        raise FileNotFoundError(f"Experiment directory {expt_dir} does not exist.")

    # Existence of metadata
    try:
        parser = find_metadata(expt_dir, Parser=ExtendedMetadataTableParser)
    except MetadataFormatError as e:
        raise MetadataFormatError(f"Error in metadata for '{expt_dir}': {e}") from e
    metadata = parser.df
    metadata.insert(0, "expt_name", os.path.basename(expt_dir))

    # Existence of regions
    regions = find_regions(expt_dir)

    # Existence of summary Files
    summary_files = get_summary_files(Path(expt_dir))
    for file_name, file_path in summary_files._asdict().items():
        if any(f in file_name for f in allow_missing_files):
            continue
        if not os.path.exists(file_path):
            raise FileNotFoundError(
                f"Missing {file_name} file: '{file_path}' in {expt_dir}."
            )

    # Existence of settings / caller
    settings_path = f"{expt_dir}/metadata/settings.json"
    if not os.path.exists(settings_path):
        settings = None
        caller = "bcftools"  # if no settings, was using bcftools
        reference_name = "Unknown"
    else:
        settings = json.load(open(settings_path, "r"))
        caller = settings["caller"]
        reference_name = str(settings["reference_name"])

    return ExperimentOutputs(
        expt_dir=expt_dir,
        metadata=metadata,
        regions=regions,
        reference_name=reference_name,
        summary_files=summary_files,
        settings=settings,
        caller=caller,
        has_complete_vcf=os.path.exists(f"{expt_dir}/vcfs/summary.variants.vcf.gz"),
        has_filtered_vcf=os.path.exists(
            f"{expt_dir}/vcfs/summary.variants.filtered.annotated.vcf.gz"
        ),
    )
