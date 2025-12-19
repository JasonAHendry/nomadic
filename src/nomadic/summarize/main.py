import glob
import os
import shutil
from typing import Iterable, Optional
from warnings import warn
from enum import StrEnum, auto
from collections import Counter
from pathlib import Path

import pandas as pd
import numpy as np

from nomadic.summarize.compute import (
    calc_amplicons_summary,
    calc_samples_summary,
    compute_variant_prevalence,
    filter_false_positives,
    gene_deletion_prevalence_by,
    gene_deletions,
)
from nomadic.summarize.dashboard.builders import BasicSummaryDashboard
from nomadic.util.panel import get_panel_settings
from nomadic.util.workspace import Workspace
from nomadic.util.dirs import produce_dir
from nomadic.util.regions import RegionBEDParser
from nomadic.util.experiment import (
    get_summary_files,
    check_experiment_outputs,
)
from nomadic.util.logging_config import LoggingFascade
from nomadic.util.summary import Settings, get_master_columns_mapping, load_settings


# --------------------------------------------------------------------------------
# Check for completion and consistency across experiments
#
# --------------------------------------------------------------------------------


def check_regions_consistent(expt_regions: list[RegionBEDParser]) -> None:
    """
    Check that the regions are consistent across all experiment directories

    TODO:
    - Might make sense to *extract* the region that was used and save it;

    """
    if len(expt_regions) == 0:
        # Nothing to check
        return
    base = expt_regions[0]
    for r in expt_regions:
        if not (r.df == base.df).all().all():
            raise ValueError(
                "Different regions used across experiments, this is not supported. Check region BED files are the same."
            )


def check_calling_consistent(expt_callers: list[str]) -> Optional[str]:
    """
    Check that the same variant caller was used across all experiments,
    where `expt_callers` is a list of used variant callers
    """
    if len(expt_callers) == 0:
        return None
    caller_counts = Counter([caller for caller in expt_callers])
    if len(caller_counts) > 1:
        raise ValueError(
            "Found more than one variant caller used across experiments: "
            + f"{', '.join([f'{v} experiment(s) used {c}' for c, v in caller_counts.items()])}."
        )
    return caller_counts.most_common()[0][0]


def get_shared_metadata_columns(
    metadata_dfs: list[pd.DataFrame],
    fixed_columns: list[str] = ["expt_name", "barcode", "sample_id", "sample_type"],
) -> list[str]:
    """Get metadata columns that are shared acrossa all experiments"""

    shared_columns = set(metadata_dfs[0].columns)
    for df in metadata_dfs[1:]:
        shared_columns.intersection_update(df.columns)
    shared_columns.difference_update(fixed_columns)  # why am I doing this?
    return list(shared_columns)


# --------------------------------------------------------------------------------
# Throughput
#
# --------------------------------------------------------------------------------


def compute_throughput(metadata: pd.DataFrame, add_unique: bool = True) -> pd.DataFrame:
    """
    Compute a simple throughput crosstable

    Also add information about uniqueness

    """

    throughput_df = pd.crosstab(
        metadata["sample_type"], metadata["expt_name"], margins=True
    )

    if add_unique:
        um = metadata.drop_duplicates("sample_id")
        throughput_df.loc["field_unique"] = pd.crosstab(
            um["sample_type"], um["expt_name"], margins=True
        ).loc["field"]

    throughput_df.fillna(0, inplace=True)
    throughput_df = throughput_df.astype(int)

    return throughput_df


def get_region_coverage_dataframe(
    expt_dirs: Iterable[str], metadata: pd.DataFrame
) -> pd.DataFrame:
    """
    Here we load a consolidated region coverage dataframe and include information required
    for quality control
    """

    # Load coverage data
    bed_dfs = []
    for expt_dir in expt_dirs:
        bed_csv = get_summary_files(Path(expt_dir)).region_coverage

        bed_df = pd.read_csv(bed_csv)
        bed_df.insert(0, "expt_name", os.path.basename(expt_dir))
        bed_df.query("barcode != 'unclassified'", inplace=True)
        if "sample_id" in bed_df.columns:
            bed_df.drop(columns=["sample_id"], inplace=True)

        # TODO: Do checks
        bed_df = pd.merge(
            left=metadata[["expt_name", "barcode", "sample_id", "sample_type"]],
            right=bed_df,
            on=["expt_name", "barcode"],
            how="inner",  # ensure we only take samples that are in the master metadata
        )
        # TODO: Do checks
        bed_dfs.append(bed_df)
    concat_df = pd.concat(bed_dfs)

    # Get negative control data
    neg_df = (
        concat_df.query("sample_type == 'neg'")
        .groupby(["expt_name", "name"])
        .mean_cov.mean()
        .reset_index()
        .rename({"mean_cov": "mean_cov_neg"}, axis=1)
    )

    # TODO: do checks
    coverage_df = pd.merge(
        left=concat_df[
            ["expt_name", "barcode", "sample_id", "sample_type", "name", "mean_cov"]
        ],  # sample ID, will want it at some point
        right=neg_df,
        on=["expt_name", "name"],
        how="left",
    )
    # TODO: do checks

    return coverage_df


def calc_unknown_samples(
    inventory_metadata: pd.DataFrame, master_metadata: pd.DataFrame
):
    """
    Returns the sample ids that are not in the masterdata file.
    """
    field_samples = inventory_metadata.query("sample_type == 'field'")
    unknown_samples = (
        field_samples.loc[
            ~field_samples["sample_id"].isin(master_metadata["sample_id"]),
            "sample_id",
        ]
        .unique()
        .tolist()
    )
    return unknown_samples


def calc_quality_control_columns(
    df: pd.DataFrame, *, min_coverage: int = 50, max_contam: float = 0.1
) -> None:
    """
    Calculate columns evaluating whether samples have passed quality control
    """

    # Do we have enough coverage?
    df["fail_lowcov"] = df["mean_cov"] < min_coverage

    # Check if coverage of negative control exceeds `max_contam`
    df["fail_contam_rel"] = df["mean_cov_neg"] / (df["mean_cov"] + 0.01) >= max_contam
    df["fail_contam_abs"] = df["mean_cov_neg"] >= min_coverage

    df["fail_contam"] = (
        (df["fail_contam_rel"] & ~df["fail_lowcov"]) | df["fail_contam_abs"]
    )  # If already failed low coverage, don't consider contamination, unless it's absolute threshold is exceeded

    # Finally, define passing
    df["passing"] = ~df["fail_contam"] & ~df["fail_lowcov"]


# --------------------------------------------------------------------------------
# Quality Control
#
# --------------------------------------------------------------------------------


class QcStatus(StrEnum):
    PASS = auto()
    LOWCOV = auto()
    CONTAM = auto()
    DUPLICATE = auto()
    CONTROL = auto()


def _add_qc_status_no_duplicates(df: pd.DataFrame) -> list[str]:
    """
    Adds a status to each replicate/amplicon to see which ones passed QC
    and if they didn't, why.
    """
    status_strs = []
    for _, row in df.iterrows():
        status = []
        if row["sample_type"] in ["pos", "neg"]:
            status.append(QcStatus.CONTROL)
        else:
            if row["fail_contam"]:
                status.append(QcStatus.CONTAM)
            if row["fail_lowcov"]:
                status.append(QcStatus.LOWCOV)
            if not status:
                status.append(QcStatus.PASS)
        status_strs.append(";".join(status))
    df["status"] = status_strs


def _mark_duplicates(df: pd.DataFrame) -> None:
    """
    Mark all field samples as duplicates, if there is a better covered replicate for the same amplicon
    Replicates marked as duplicate will not be used for prevalance evaluation.
    """

    def _update_duplicate(status: str, idx: int, keep_idx: int) -> str:
        if idx == keep_idx:
            return status
        return f"{status};{QcStatus.DUPLICATE}"

    for (_, _), data in df.query("sample_type == 'field'").groupby(
        ["sample_id", "name"]
    ):
        # Select an index to keep, i.e. the best sample that should
        # marked as duplicate
        passing = data["status"] == "pass"
        if passing.sum() == 1:
            keep_idx = passing.idxmax()
        elif passing.sum() > 1:
            keep_idx = data[passing]["mean_cov"].idxmax()  # keep maximum coverage
        else:  # none are passing, arbrarily keep first
            keep_idx = data.index[0]

        # NB: updating in-place
        # must be a view, not a slice hence [,]
        df.loc[data.index, "status"] = [
            _update_duplicate(status, idx, keep_idx)
            for idx, status in data["status"].items()
        ]


def add_quality_control_status_column(df: pd.DataFrame) -> None:
    """
    Add a QC status column in-place

    Note:
    - When we mark duplicates; we do it on an AMPLICON x SAMPLE level;
    not on a per-sample level. So we could take amplicons from separate
    samples to get the best data for that sample.

    """
    _add_qc_status_no_duplicates(df)
    _mark_duplicates(df)


# --------------------------------------------------------------------------------
# Variant analysis
#
# --------------------------------------------------------------------------------


def load_variant_summary_csv(
    variants_csv: str, define_gene: bool = True
) -> pd.DataFrame:
    """
    Load an clean `summary.variants.csv` data produced by `nomadic`

    """

    # Settings
    NUMERIC_COLUMNS = ["dp", "wsaf"]
    UNPHASED_GT_TO_INT = {
        "./.": -1,
        "0/0": 0,
        "0/1": 1,
        "1/1": 2,
    }  # TODO What about multiallelic sites

    # Load
    variants_df = pd.read_csv(variants_csv)

    # Reformat numeric columns to be floats
    for c in NUMERIC_COLUMNS:
        variants_df[c] = [float(v) if v != "." else None for v in variants_df[c]]

    # Reformat unphased genotypes as integers
    variants_df.insert(
        variants_df.columns.tolist().index("gt") + 1,
        "gt_int",
        variants_df["gt"].map(UNPHASED_GT_TO_INT),
    )

    # Optionally reformat amplicon name to gene; assuming like  gene-...
    if define_gene:
        variants_df.insert(
            variants_df.columns.get_loc("amplicon") + 1,
            "gene",
            [a.split("-")[0] for a in variants_df["amplicon"]],
        )
        variants_df.insert(
            variants_df.columns.get_loc("gene") + 1,
            "mutation",
            [
                f"{gene}-{aa_change}"
                for gene, aa_change in zip(
                    variants_df["gene"], variants_df["aa_change"]
                )
            ],
        )

    return variants_df


def load_and_concat_variants(expt_dirs: list[str]) -> pd.DataFrame:
    """
    Load all of the variant calls for a set of experiment dirs

    Note that because we do note have the unfiltered VCF files, we have to do
    some additional work in order to ensure all mutations are represented across
    all experiments;
    """

    # Load data
    variant_dfs = []
    for expt_dir in expt_dirs:
        variant_csv = f"{expt_dir}/summary.variants.csv"
        variant_df = load_variant_summary_csv(variant_csv)
        variant_df.insert(0, "expt_name", os.path.basename(expt_dir))
        variant_df.query("barcode != 'unclassified'", inplace=True)
        variant_dfs.append(variant_df)
    variant_df = pd.concat(variant_dfs)

    # Get all unique mutations
    MUT_COLUMNS = [
        "chrom",
        "pos",
        "ref",
        "alt",
        "strand",
        "aa_change",
        "aa_pos",
        "mut_type",
        "mutation",
        "amplicon",
        "gene",
    ]
    uniq_mutation_df = variant_df[MUT_COLUMNS].drop_duplicates()

    # Now we merge these back in for each barcode
    # - By doing a 'right' merge, we make sure all variants are present for each barcode
    # - The 'NaNs' that are present when the variant doesn't exist for that barcode
    #   get filled with zeros, i.e. we assume homozygous reference
    # - In reality, it could have been EITHER 0/0 or ./. (i.e. filtered), but we handle
    #   this afterwards when we merge with QC data;
    full_variant_dfs = []
    for (expt_name, barcode), bdf in variant_df.groupby(["expt_name", "barcode"]):
        mdf = pd.merge(bdf, uniq_mutation_df, on=MUT_COLUMNS, how="right")
        mdf["expt_name"] = expt_name
        mdf["barcode"] = barcode

        # Filled by default with hom reference
        mdf["gt"] = mdf["gt"].fillna("0/0")
        mdf["gt_int"] = mdf["gt_int"].fillna(0.0)
        mdf["wsaf"] = mdf["wsaf"].fillna(0.0)

        full_variant_dfs.append(mdf)

    return pd.concat(full_variant_dfs)


def replicates_qc(
    coverage_df: pd.DataFrame, REPLICATE_PASSING_THRESHOLD: float
) -> pd.DataFrame:
    """
    Calculates which of the replicates (repeated runs of a sample) passed QC as a whole
    (more than REPLICATE_PASSING_THRESHOLD passed)
    """
    replicates_qc_df = (
        coverage_df.query("sample_type == 'field'")
        .groupby(["expt_name", "barcode", "sample_id"])
        .agg(
            n_amplicons=pd.NamedAgg("name", "count"),
            n_passing=pd.NamedAgg("passing", "sum"),
            n_fail_contam=pd.NamedAgg("fail_contam", "sum"),
            n_fail_lowcov=pd.NamedAgg("fail_lowcov", "sum"),
        )
        .reset_index()
    )
    replicates_qc_df["passing"] = (
        replicates_qc_df["n_passing"] / replicates_qc_df["n_amplicons"]
        >= REPLICATE_PASSING_THRESHOLD
    )

    return replicates_qc_df


def replicates_amplicon_qc(coverage_df):
    return coverage_df.query("sample_type == 'field'")


# --------------------------------------------------------------------------------
# Main
#
# --------------------------------------------------------------------------------


def main(
    *,
    workspace: Workspace,
    output_dir: Path,
    expt_dirs: tuple[str],
    summary_name: str,
    metadata_path: Optional[Path],
    settings_file_path: Path,
    show_dashboard: bool = True,
    prevalence_by: list[str],
    no_master_metadata: bool = False,
    qc_min_coverage: int,
    qc_max_contam: float,
) -> None:
    """
    Define the main function for the summary analysis


    """

    assert (metadata_path is not None) or no_master_metadata

    produce_dir(str(output_dir))

    # PARSE EXPERIMENT DIRECTORIES
    log = LoggingFascade(logger_name="nomadic")
    log.info("Input parameters:")
    log.info(f"  Summary Name: {summary_name}")
    if not no_master_metadata:
        log.info(f"  Master metadata: {metadata_path}")
    else:
        log.info("  No master metadata will be used.")
    log.info(f"  Setting file: {settings_file_path}")
    log.info(f"  Found {len(expt_dirs)} experiment directories.")

    # Check experiments are complete
    expts = [check_experiment_outputs(expt_dir) for expt_dir in expt_dirs]
    log.info("  All experiments are complete.")

    # Check experiments are consistent
    check_regions_consistent([expt.regions for expt in expts])
    log.info("  All experiments use the same regions.")
    if expts:
        panel_name = expts[0].regions.name
    else:
        panel_name = "Unknown"
    log.info(f"  Panel used: {panel_name}")
    caller = check_calling_consistent([expt.caller for expt in expts])
    log.info(f"  All experiments use same variant caller: {caller}")

    panel_settings = get_panel_settings(panel_name)
    log.info(f"  Loaded panel settings for panel '{panel_settings.name}'.")

    settings: Settings = Settings()
    if settings_file_path.exists():
        settings = load_settings(settings_file_path)
        log.info(f"  Loaded summary settings from {settings_file_path}.")

    # CHECK METADATA IS VALID
    FIXED_COLUMNS = ["expt_name", "barcode", "sample_id", "sample_type"]
    shared_columns = get_shared_metadata_columns(
        [expt.metadata for expt in expts], fixed_columns=FIXED_COLUMNS
    )
    log.info("  All metadata tables pass completion checks.")
    log.info(
        f"  Found {len(shared_columns)} non-essential shared columns across all metadata files: {', '.join(shared_columns)}"
    )

    # for now we use the master metadata file
    inventory_metadata = pd.concat([expt.metadata[FIXED_COLUMNS] for expt in expts])
    if metadata_path is not None and not no_master_metadata:
        master_metadata = pd.read_csv(metadata_path).rename(
            columns=get_master_columns_mapping(settings)
        )
    else:
        # create metadata from experiment meta data files
        shared_columns = ["sample_id"] + list(shared_columns)
        master_metadata = pd.concat(
            [
                expt.metadata.query("sample_type == 'field'")[shared_columns]
                for expt in expts
            ]
        )

    master_metadata = master_metadata.astype(
        {"sample_id": "str"}
    )  # ensure sample IDs are strings
    inventory_metadata = inventory_metadata.astype(
        {"sample_id": "str"}
    )  # ensure sample IDs are strings
    # strip whitespaces from sample IDs
    inventory_metadata["sample_id"] = inventory_metadata["sample_id"].str.strip()
    master_metadata["sample_id"] = master_metadata["sample_id"].str.strip()

    unknown_samples = calc_unknown_samples(inventory_metadata, master_metadata)
    if unknown_samples:
        warn(
            f"Samples in experiments that are not in master metadata: {unknown_samples}"
        )

    # Delete unknown samples
    inventory_metadata["status"] = np.where(
        inventory_metadata["sample_id"].isin(unknown_samples), "unknown", "known"
    )
    # set all controls
    inventory_metadata["status"] = np.where(
        inventory_metadata["sample_type"].isin(["pos", "neg"]),
        "control",
        inventory_metadata["status"],
    )
    inventory_metadata.to_csv(f"{output_dir}/inventory.csv", index=False)
    inventory_metadata = inventory_metadata.query("status != 'unknown'")

    # Throughput data
    log.info("Overall sequencing throughput:")
    throughput_df = compute_throughput(inventory_metadata)
    log.info(f"  Positive controls: {throughput_df.loc['pos', 'All']}")
    log.info(f"  Negative controls: {throughput_df.loc['neg', 'All']}")
    log.info(f"  Fields samples sequenced (total): {throughput_df.loc['field', 'All']}")
    log.info(f"  Field samples (unique): {throughput_df.loc['field_unique', 'All']}")
    log.info(f"  Unknown samples (excluded): {len(unknown_samples)}")
    throughput_df.to_csv(f"{output_dir}/summary.throughput.csv", index=True)

    # Now let's evaluate coverage
    coverage_df = get_region_coverage_dataframe(expt_dirs, inventory_metadata)
    calc_quality_control_columns(
        coverage_df, min_coverage=qc_min_coverage, max_contam=qc_max_contam
    )

    log.info("Amplicon-Sample QC Statistics:")
    field_coverage_df = coverage_df.query("sample_type == 'field'")
    n = field_coverage_df.shape[0]
    n_lowcov = field_coverage_df["fail_lowcov"].sum()
    n_contam = field_coverage_df["fail_contam"].sum()
    n_pass = field_coverage_df["passing"].sum()
    log.info(
        f"  Coverage below <{qc_min_coverage}x: {n_lowcov} ({100 * n_lowcov / n:.2f}%)"
    )
    log.info(
        f"  Contamination >{qc_max_contam}: {n_contam} ({100 * n_contam / n:.2f}%)"
    )
    log.info(f"  Passing QC: {n_pass} ({100 * n_pass / n:.2f}%)")
    add_quality_control_status_column(coverage_df)
    log.info(str(coverage_df["status"].value_counts()))
    coverage_df.to_csv(f"{output_dir}/summary.coverage.csv", index=False)

    REPLICATE_PASSING_THRESHOLD = 0.8
    replicates_qc_df = replicates_qc(coverage_df, REPLICATE_PASSING_THRESHOLD)
    replicates_qc_df.to_csv(f"{output_dir}/summary.replicates_qc.csv", index=False)

    samples_summary_df = calc_samples_summary(master_metadata, replicates_qc_df)
    samples_summary_df.to_csv(f"{output_dir}/summary.samples_qc.csv", index=False)

    samples_by_amplicon_summary_df = calc_amplicons_summary(
        master_metadata, replicates_amplicon_qc(coverage_df)
    )
    samples_by_amplicon_summary_df.to_csv(
        f"{output_dir}/summary.samples_amplicons_qc.csv", index=False
    )

    final_df = (
        coverage_df.query("sample_type == 'field'")
        .groupby(["expt_name", "name"])
        .agg(
            mean_cov_field=pd.NamedAgg("mean_cov", "median"),
            mean_cov_neg=pd.NamedAgg("mean_cov_neg", "median"),
            n_field=pd.NamedAgg("barcode", len),
            n_field_passing=pd.NamedAgg("passing", lambda x: x.sum()),
            per_field_contam=pd.NamedAgg("fail_contam", lambda x: 100 * x.mean()),
            per_field_lowcov=pd.NamedAgg("fail_lowcov", lambda x: 100 * x.mean()),
            per_field_passing=pd.NamedAgg("passing", lambda x: 100 * x.mean()),
        )
        .reset_index()
    )
    final_df.to_csv(f"{output_dir}/summary.experiments_qc.csv", index=False)

    # --------------------------------------------------------------------------------
    # Let's move onto to variant calling results
    #
    # --------------------------------------------------------------------------------

    log.info("Loading variants...")
    variant_df = load_and_concat_variants(expt_dirs)

    if "sample_id" in variant_df.columns:
        variant_df.drop(columns=["sample_id"], inplace=True)
    # Merge with the quality control results, then we can subset to the analysis set
    variant_df = pd.merge(
        left=coverage_df.rename({"name": "amplicon"}, axis=1)[
            ["expt_name", "barcode", "sample_id", "sample_type", "amplicon", "status"]
        ],
        right=variant_df,
        on=["expt_name", "barcode", "amplicon"],
    )

    log.info("Filtering to analysis set...")
    remove_amplicons = panel_settings.excluded_amplicons  # noqa: F841 later used in query
    remove_mutations = panel_settings.filtered_mutations  # noqa: F841 later used in query
    analysis_df = (
        variant_df.query("status == 'pass'")
        .query("mut_type == 'missense'")
        .query("amplicon not in @remove_amplicons")
        .query("mutation not in @remove_mutations")
    )

    # Filter out false positives
    analysis_df = filter_false_positives(analysis_df)
    analysis_df.to_csv(f"{output_dir}/summary.variants.analysis_set.csv", index=False)

    # Then we will compute prevalence
    prev_df = compute_variant_prevalence(analysis_df)
    prev_df.to_csv(f"{output_dir}/summary.variants.prevalence.csv", index=False)

    for col in prevalence_by:
        prev_by_col_df = compute_variant_prevalence(analysis_df, master_metadata, [col])
        prev_by_col_df.to_csv(
            f"{output_dir}/summary.variants.prevalence-{col}.csv", index=False
        )

    # --------------------------------------------------------------------------------
    # Gene deletion analysis
    #
    # --------------------------------------------------------------------------------

    if panel_settings.deletion_genes:
        log.info("Calculate gene deletions...")
        gene_deletion_df = gene_deletions(coverage_df, panel_settings.deletion_genes)
        gene_deletion_df.to_csv(f"{output_dir}/summary.gene_deletions.csv", index=False)

        prev_gen_deletions_df = gene_deletion_prevalence_by(
            gene_deletion_df, master_metadata, []
        )
        prev_gen_deletions_df.to_csv(
            f"{output_dir}/summary.gene-deletions.prevalence.csv", index=False
        )

    for col in prevalence_by:
        prev_gen_deletion_by_col_df = gene_deletion_prevalence_by(
            gene_deletion_df, master_metadata, [col]
        )
        prev_gen_deletion_by_col_df.to_csv(
            f"{output_dir}/summary.gene-deletions.prevalence-{col}.csv", index=False
        )

    master_metadata.to_csv(f"{output_dir}/metadata.csv", index=False)

    log.info("Copy relevant files to summary output directory...")

    # Copy relevant files
    if expts:
        shutil.copy(
            expts[0].regions.path,
            os.path.join(output_dir, os.path.basename(expts[0].regions.path)),
        )
    for file in glob.glob(f"{workspace.get_metadata_dir()}/{summary_name}*.geojson"):
        shutil.copy(file, f"{output_dir}/{file.split('-')[-1]}")
    coords_file = f"{workspace.get_metadata_dir()}/{summary_name}.coords.csv"
    if os.path.isfile(coords_file):
        shutil.copy(
            coords_file,
            os.path.join(output_dir, "coords.csv"),
        )
    summary_settings_file = workspace.get_summary_settings_file(summary_name)
    if os.path.isfile(summary_settings_file):
        shutil.copy(
            summary_settings_file,
            os.path.join(output_dir, "settings.yaml"),
        )

    log.info("Summary analysis complete.")

    # --------------------------------------------------------------------------------
    # Dashboard
    #
    # --------------------------------------------------------------------------------

    if show_dashboard:
        view(output_dir, summary_name)


def view(input_dir: Path, summary_name: str) -> None:
    """
    View the summary dashboard for a given summary
    """
    settings: Settings = Settings()
    settings_file = Path(f"{input_dir}/settings.yaml")
    if settings_file.exists():
        print(f"Loading settings from {settings_file}...")
        settings = load_settings(settings_file)

    bed_files = glob.glob(f"{input_dir}/*.bed")
    if bed_files:
        panel_name = os.path.basename(bed_files[0]).split(".")[0]
        print(f"Use panel name from regions BED file: {panel_name}")
        amplicons = RegionBEDParser(bed_files[0]).names
    else:
        raise ValueError("No regions BED file found in summary directory.")

    panel_settings = get_panel_settings(panel_name)
    amplicon_sets = panel_settings.amplicon_sets
    deletion_genes = panel_settings.deletion_genes

    print("Load data...")

    dashboard = BasicSummaryDashboard(
        summary_name,
        throughput_csv=f"{input_dir}/summary.throughput.csv",
        samples_csv=f"{input_dir}/summary.samples_qc.csv",
        samples_amplicons_csv=f"{input_dir}/summary.samples_amplicons_qc.csv",
        coverage_csv=f"{input_dir}/summary.experiments_qc.csv",
        analysis_csv=f"{input_dir}/summary.variants.analysis_set.csv",
        gene_deletions_csv=f"{input_dir}/summary.gene_deletions.csv",
        master_csv=f"{input_dir}/metadata.csv",
        geojson_glob=f"{input_dir}/*.geojson",
        location_coords_csv=f"{input_dir}/coords.csv",
        settings=settings,
        amplicons=amplicons,
        amplicon_sets=amplicon_sets,
        deletion_genes=deletion_genes,
    )

    print("")
    print("Launching dashboard (press CNTRL+C to exit):")
    print("")
    debug = bool(os.getenv("NOMADIC_DEBUG"))
    dashboard.run(debug=debug, auto_open=not debug)
