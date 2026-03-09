import glob
import os
import shutil
from typing import Iterable, Optional
from warnings import warn
from enum import StrEnum, auto
from collections import Counter
from pathlib import Path
import subprocess
import time

import pandas as pd
import numpy as np

from nomadic.download.references import REFERENCE_COLLECTION
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
from nomadic.util.vcf import VariantAnnotator
from nomadic.util.workspace import Workspace
from nomadic.util.dirs import produce_dir
from nomadic.util.regions import RegionBEDParser
from nomadic.util.experiment import (
    get_summary_files,
    check_experiment_outputs,
)
from nomadic.util.logging_config import LoggingFascade
from nomadic.util.summary import Settings, get_master_columns_mapping, load_settings
from nomadic.util.wrappers import bcftools


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
        "0/2": 1,
        "0/3": 1,
        "1/1": 2,
        "2/2": 2,
        "3/3": 2,
    }

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
        # variants_df.insert(
        #     variants_df.columns.get_loc("amplicon") + 1,
        #     "gene",
        #     [a.split("-")[0] for a in variants_df["amplicon"]],
        # )
        # variants_df.insert(
        #     variants_df.columns.get_loc("gene") + 1,
        #     "mutation",
        #     [
        #         f"{gene}-{aa_change}"
        #         for gene, aa_change in zip(
        #             variants_df["gene"], variants_df["aa_change"]
        #         )
        #     ],
        # )
        # --- gene ---
        gene_values = variants_df["amplicon"].str.split("-").str[0]

        variants_df.insert(
            variants_df.columns.get_loc("amplicon") + 1,
            "gene",
            gene_values.where(variants_df["amplicon"].notna()),
        )

        # --- mutation ---
        mutation_values = variants_df["gene"] + "-" + variants_df["aa_change"]

        variants_df.insert(
            variants_df.columns.get_loc("gene") + 1,
            "mutation",
            mutation_values.where(
                variants_df["gene"].notna() & variants_df["aa_change"].notna()
            ),
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


def load_variants_from_vcfs(
    expt_dirs: Iterable[str],
    *,
    caller: str,
    output_dir: Path,
    summary_regions: RegionBEDParser,
    reference_name: str,
) -> pd.DataFrame:
    """
    Load variants directly from VCF files, rather than the summary CSVs

    This is more work, but allows us to get all variants that were called, even those that didn't pass filtering and thus don't appear in the summary CSVs.
    """

    timer = Timer()
    timer.start()
    seperator = "___"

    if any(seperator in expt_dir for expt_dir in expt_dirs):
        raise ValueError(
            f"Experiment directories can not contain the string '{seperator}', as this is used to separate experiment name and barcode when loading from VCFs. Please rename the following directories: {', '.join([d for d in expt_dirs if seperator in d])}."
        )

    temp_dir = output_dir / "temp_vcf_processing"
    temp_dir.mkdir(exist_ok=True)

    # Record all samples for for sanity check after
    # expt_name -> set of samples in that experiment
    experiment_sample_mapping: dict[str, set[str]] = dict()

    timer.time("setup")

    for expt_dir in expt_dirs:
        expt_name = os.path.basename(expt_dir)
        vcf_dir = Path(expt_dir) / "vcfs"
        vcf_file = vcf_dir / "summary.variants.vcf.gz"
        if not vcf_file.exists():
            raise FileNotFoundError(
                f"VCF file not found at {vcf_file}. Cannot load variants from VCFs."
            )

        experiment_samples = (
            subprocess.check_output(["bcftools", "query", "-l", str(vcf_file)])
            .decode()
            .splitlines()
        )
        experiment_sample_set = set(experiment_samples)
        if len(experiment_sample_set) != len(experiment_samples):
            raise ValueError(
                f"Duplicate sample names found in VCF file {vcf_file}. Sample names must be unique to load from VCFs."
            )
        experiment_sample_mapping[expt_name] = experiment_sample_set
        # Make sample names unique by concating expt name and experiment sample name (the barcode)
        unique_samples = [f"{expt_name}{seperator}{s}" for s in experiment_samples]
        unique_samples = [
            s.replace(" ", r"\ ") for s in unique_samples
        ]  # replace space, bcftools treat spaces in samples names special

        ### Reheader vcf file and move
        temp_vcf = temp_dir / f"{Path(expt_dir).name}.variants.temp.vcf.gz"
        subprocess.run(
            [
                "bcftools",
                "reheader",
                "-n",
                ",".join(unique_samples),
                "-o",
                str(temp_vcf),
                str(vcf_file),
            ],
            check=True,
        )
        bcftools.index(temp_vcf)

    timer.time("Loading and reheadering VCF files")

    # Now we can merge all the temp VCFs together
    temp_vcfs = list(temp_dir.glob("*.temp.vcf.gz"))
    merged_vcf = output_dir / "summary.variants.vcf.gz"

    subprocess.run(
        ["bcftools", "merge", "-Fx", "-Oz", "--force-single", "-o", str(merged_vcf)]
        + [str(v) for v in temp_vcfs],
        check=True,
    )

    bcftools.index(merged_vcf)
    timer.time("Merging VCF files")

    # Filtering
    filtered_vcf = output_dir / "summary.variants.filtered.vcf.gz"
    subprocess.run(
        [
            "bcftools",
            "view",
            "--apply-filters",
            "PASS",
            "--types",
            "snps",
            "--min-alleles",
            "2",
            "-Oz",
            "-o",
            str(filtered_vcf),
            str(merged_vcf),
        ],
        check=True,
    )

    timer.time("Filtering VCF file")

    REFERENCE_COLLECTION[reference_name].confirm_downloaded()
    annotator = VariantAnnotator(
        input_vcf=str(filtered_vcf),
        bed_path=summary_regions.path,
        reference=REFERENCE_COLLECTION[reference_name],
        caller=caller,
        output_vcf=str(output_dir / "summary.variants.annotated.vcf.gz"),
    )
    annotator.run()

    timer.time("Annotating variants")

    # TODO saving and loading of this file should be removed
    annotator.convert_to_csv(f"{temp_dir}/summary.variants.merged.csv")

    timer.time("Converting VCF to CSV")

    variant_df = load_variant_summary_csv(f"{temp_dir}/summary.variants.merged.csv")

    timer.time("Loading variants into dataframe")

    # TODO the following code should probably be combined with the anotator
    variant_sample_names = variant_df["barcode"].str.replace(
        r"\ ", " "
    )  # reverse escaping of spaces
    variant_df.insert(0, "expt_name", variant_sample_names.str.split(seperator).str[0])
    variant_df["barcode"] = variant_sample_names.str.split(seperator).str[1]

    # Sanity checks that all worked and we have the same samples as before
    # Bcftools has some magic around sample names, so good to check they are consistent with what we expected
    all_samples_after = variant_df.groupby(["expt_name"])["barcode"].unique().to_dict()
    assert set(experiment_sample_mapping.keys()) == set(all_samples_after.keys()), (
        "Mismatch in experiments after loading from VCFs."
    )
    assert all(
        samples == set(all_samples_after[expt_dir])
        for expt_dir, samples in experiment_sample_mapping.items()
    ), "Mismatch in samples after loading from VCFs."

    timer.time("fixing sample names and sanity checking")

    variant_df = variant_df.query("aa_pos != -1")

    timer.time("Filtering out variants without aa_pos")

    # identifies an aa variant, note this does not include aa_change, as there might be different
    # changes for the same position
    groups = ["expt_name", "barcode", "chrom", "amplicon", "gene", "aa_pos"]
    # wt_df = (
    #     variant_df.groupby(groups)
    #     .agg(
    #         WT=pd.NamedAgg(column="gt_int", aggfunc=lambda x: (x == 0).all()),
    #         ref=("ref", "first"),
    #         alt=("alt", "first"),
    #         pos=("pos", "first"),
    #         qual=("qual", "first"),
    #         dp=("dp", "first"),
    #         wsaf=("wsaf", "first"),
    #         gt=("gt", "first"),
    #         gt_int=("gt_int", "first"),
    #     )
    #     .reset_index()
    #     .query("WT == True")
    # )
    # wt_df.drop(columns=["WT"], inplace=True)
    # wt_df["type"] = "wt"

    # fully vectorized way to get WT variants:
    # Wt is where gt_int is 0 for all variants in the group; if any variant has gt_int != 0, then it's not WT
    is_zero = variant_df["gt_int"].eq(0)
    wt_mask = is_zero.groupby([variant_df[c] for c in groups]).transform("all")
    # TODO, decide how to record gt, wsaf, dp here, currently just taking first
    wt_df = variant_df.loc[wt_mask].groupby(groups, as_index=False).first()
    wt_df["type"] = "wt"

    timer.time("Processing WT variants into final dataframe")

    # filtered_df = (
    #     variant_df.groupby(groups)
    #     .agg(
    #         filtered=pd.NamedAgg(column="gt_int", aggfunc=lambda x: (x == -1).any()),
    #         ref=("ref", "first"),
    #         alt=("alt", "first"),
    #         pos=("pos", "first"),
    #         qual=("qual", "first"),
    #         dp=("dp", "first"),
    #         wsaf=("wsaf", "first"),
    #         gt=("gt", "first"),
    #         gt_int=("gt_int", "first"),
    #     )
    #     .reset_index()
    #     .query("filtered == True")
    # )
    # filtered_df.drop(columns=["filtered"], inplace=True)
    # filtered_df["type"] = "filtered"

    # Fully vectorized way to get filtered variants: filtered is where gt_int is -1 for any variant in the group
    # TODO, decide if we actually want to record this at all
    is_filtered = variant_df["gt_int"].eq(-1)
    filtered_mask = is_filtered.groupby([variant_df[c] for c in groups]).transform(
        "any"
    )
    filtered_df = variant_df.loc[filtered_mask].groupby(groups, as_index=False).first()
    filtered_df["type"] = "filtered"

    timer.time("Processing filtered variants into final dataframe")

    # TODO, decide also here how to record gt, wsaf, dp. Currently just taking the one
    # that is in the row of the vcf where the aa_change was recorded
    mut_df = variant_df.dropna(subset=["aa_change"])
    mut_df["type"] = (
        mut_df.groupby(groups + ["aa_change"])["gt_int"]
        .transform(lambda x: (x == 1).any())
        .map({True: "mixed_mut", False: "mut"})
    )

    timer.time("Processing mutant variants into final dataframe")

    variant_df = pd.concat([wt_df, filtered_df, mut_df]).reset_index(drop=True)

    # drop duplicates if they exist. This should only happen when bcftools record a variant, but we want to filter it
    # drop so that type=filtered stays
    # NOTE this will have to be changed once we allow more than one mut per sample
    pref = variant_df["type"].eq("filtered")
    idx = pref.groupby([variant_df[c] for c in groups]).idxmax()
    variant_df = variant_df.loc[idx].reset_index(drop=True)

    # sanity checks
    assert variant_df.duplicated(subset=groups).sum() == 0, (
        "Duplicate variants found in final dataframe after processing. This should not happen"
    )

    timer.time("Concatenating final dataframe")

    timer.report()

    return variant_df


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
    maps: list[str],
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
    if caller is None:
        raise ValueError("Can only summarize variants if a variant caller was used")
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
    inventory_metadata = pd.concat(
        [expt.metadata[FIXED_COLUMNS] for expt in expts]
    ).reset_index()
    if metadata_path is not None and not no_master_metadata:
        master_metadata = pd.read_csv(metadata_path, dtype={"sample_id": "str"}).rename(
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
        # Note, problematic if same sample ID has different metadata across experiments
        master_metadata.drop_duplicates(subset=["sample_id"], inplace=True)

    master_metadata = master_metadata.astype(
        {"sample_id": "str"}
    )  # ensure sample IDs are strings
    inventory_metadata = inventory_metadata.astype(
        {"sample_id": "str"}
    )  # ensure sample IDs are strings
    # strip whitespaces from sample IDs
    inventory_metadata["sample_id"] = inventory_metadata["sample_id"].str.strip()
    master_metadata["sample_id"] = master_metadata["sample_id"].str.strip()

    # Check no duplicate sample IDs in master metadata
    duplicate_sample_ids = master_metadata[
        master_metadata.duplicated(subset=["sample_id"])
    ]["sample_id"].unique()
    if len(duplicate_sample_ids) > 0:
        raise ValueError(
            f"Duplicate sample IDs found in master metadata: {', '.join(duplicate_sample_ids)}"
        )

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

    if len(inventory_metadata.query("sample_type == 'field'")) == 0:
        log.info("No known field samples, exiting...")
        return

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
    # variant_df = load_and_concat_variants(expt_dirs)

    # if "sample_id" in variant_df.columns:
    #     variant_df.drop(columns=["sample_id"], inplace=True)

    timer = Timer()
    timer.start()
    variant_df = load_variants_from_vcfs(
        expt_dirs,
        caller=caller,
        output_dir=output_dir,
        summary_regions=expts[0].regions,
        reference_name="Pf3D7",
    )
    timer.time("Loading and annotating variants from VCFs")

    variant_df.to_csv(f"{output_dir}/intermediate.summary.variants.csv", index=False)

    # # Merge with the quality control results, then we can subset to the analysis set
    variant_df = pd.merge(
        left=coverage_df.rename({"name": "amplicon"}, axis=1)[
            ["expt_name", "barcode", "sample_id", "sample_type", "amplicon", "status"]
        ],
        right=variant_df,
        on=["expt_name", "barcode", "amplicon"],
    )
    timer.time("Merging variants with coverage data")

    log.info("Filtering to analysis set...")
    remove_amplicons = panel_settings.excluded_amplicons  # noqa: F841 later used in query
    remove_mutations = panel_settings.filtered_mutations  # noqa: F841 later used in query
    analysis_df = (
        variant_df.query("status == 'pass'")
        .query("amplicon not in @remove_amplicons")
        .query("mutation not in @remove_mutations")
    )
    timer.time("Filtering to analysis set")

    # Filter out false positives
    analysis_df = filter_false_positives(analysis_df)
    analysis_df.to_csv(f"{output_dir}/summary.variants.analysis_set.csv", index=False)
    timer.time("Filtering false positives")

    # Then we will compute prevalence
    prev_df = compute_variant_prevalence(analysis_df)
    prev_df.to_csv(f"{output_dir}/summary.variants.prevalence.csv", index=False)
    timer.time("Computing variant prevalence")

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
    for map_name in maps:
        file = Path(workspace.path) / "maps" / f"{map_name}.geojson"
        if file.exists():
            shutil.copy(file, f"{output_dir}/{map_name.split('-')[-1]}.geojson")
    coords_file = f"{workspace.get_metadata_dir()}/{summary_name}.coords.csv"
    if os.path.isfile(coords_file):
        shutil.copy(
            coords_file,
            os.path.join(output_dir, "coords.csv"),
        )
    if os.path.isfile(settings_file_path):
        shutil.copy(
            settings_file_path,
            os.path.join(output_dir, "settings.yaml"),
        )

    log.info("Summary analysis complete.")

    timer.report()

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
    print(f'View summary dashboard for "{summary_name}".')
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

    print(f"Load data from {input_dir}...")

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


class Timer:
    """
    Simple timer class for measuring execution time of code blocks
    """

    def __init__(self):
        self.start_time = None
        self.timings = {}

    def start(self):
        self.start_time = time.time()

    def time(self, name: str):
        if self.start_time is None:
            raise RuntimeError("time can only be called after start")
        end_time = time.time()
        elapsed_time = end_time - self.start_time
        self.timings[name] = elapsed_time
        self.start_time = end_time

    def report(self):
        print("Execution time report:")
        for name, timing in self.timings.items():
            print(f"  {name}: {timing:.2f} seconds")
