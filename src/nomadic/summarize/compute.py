import enum
from typing import Optional

import pandas as pd
from statsmodels.stats.proportion import proportion_confint


# These columns are used to define unique variants
VARIANTS_GROUP_COLUMNS = [
    "gene",
    "chrom",
    "pos",
    # "ref",
    # The alt might differ for multiallelic sites
    # "alt",
    "aa_change",
    "aa_pos",
    "mut_type",
    "mutation",
]


class Status(enum.Enum):
    PASSING = "passing"
    FAILING = "failing"
    NOT_SEQUENCED = "not_sequenced"


def calc_samples_summary(
    master_metadata_df: pd.DataFrame, replicates_qc_df: pd.DataFrame
) -> pd.DataFrame:
    """
    Calculates a summary of which samples have how many replicates that are passing or failing,
    and if it has at least one passing replicate, it is concidered as passing.

    This can be used to get a list of to be resequenced samples.
    """
    samples_summary_df = (
        replicates_qc_df.groupby(["sample_id"])
        .agg(
            n_replicates=pd.NamedAgg("barcode", "count"),
            n_passing=pd.NamedAgg("passing", "sum"),
        )
        .reset_index()
    )
    samples_summary_df = (
        samples_summary_df.merge(
            master_metadata_df[["sample_id"]], how="right", on="sample_id"
        )
        .fillna({"n_replicates": 0, "n_passing": 0})
        .astype({"n_replicates": int, "n_passing": int})
    )
    samples_summary_df["status"] = samples_summary_df.apply(
        lambda row: Status.PASSING.value
        if row["n_passing"] > 0
        else Status.FAILING.value
        if row["n_replicates"] > 0
        else Status.NOT_SEQUENCED.value,
        axis=1,
    )
    samples_summary_df.sort_values(
        by=["n_passing", "n_replicates", "sample_id"],
        inplace=True,
        ascending=[False, False, True],
    )

    return samples_summary_df


def calc_amplicons_summary(master_metadata, replicates_amplicon_qc_df):
    """
    Calculates a summary of which samples have how many replicates per amplicon that are passing or failing,
    and if it has at least one passing replicate over that amplicon, it is concidered as passing.

    This can be used to understand if there are certain amplicons of samples that have no coverage yet
    and make decisions on resampling, that are more fine granular than per sample.
    """
    samples_by_amplicons_summary_df = (
        replicates_amplicon_qc_df.groupby(["sample_id", "name"])
        .agg(
            n_replicates=pd.NamedAgg("barcode", "count"),
            n_passing=pd.NamedAgg("passing", "sum"),
        )
        .reset_index()
    )
    samples_by_amplicons_summary_df = (
        samples_by_amplicons_summary_df.merge(
            master_metadata[["sample_id"]], how="right", on="sample_id"
        )
        .fillna({"n_replicates": 0, "n_passing": 0})
        .astype({"n_replicates": int, "n_passing": int})
    )
    samples_by_amplicons_summary_df["status"] = samples_by_amplicons_summary_df.apply(
        lambda row: Status.PASSING.value
        if row["n_passing"] > 0
        else Status.FAILING.value
        if row["n_replicates"] > 0
        else Status.NOT_SEQUENCED.value,
        axis=1,
    )
    samples_by_amplicons_summary_df.sort_values(
        by=["n_passing", "n_replicates", "sample_id"],
        inplace=True,
        ascending=[False, False, True],
    )

    return samples_by_amplicons_summary_df


def filter_false_positives(
    variants_df: pd.DataFrame, min_obs: int = 1, min_wsaf: float = 0.15
):
    """Filter out likely false positive variant calls.

    Remove variants that have only been observed `min_obs` times across all samples in
    the analysis set, and with a WSAF below `min_wsaf`

    This removes likely false-positive mutations, which are usually rare and at low
    WSAF.

    """

    mut = variants_df[variants_df["gt_int"].isin([1, 2])]
    df = variants_df.merge(
        mut.groupby(VARIANTS_GROUP_COLUMNS).agg(
            n_mut=pd.NamedAgg("gt_int", len), wsaf_max=pd.NamedAgg("wsaf", "max")
        ),
        on=VARIANTS_GROUP_COLUMNS,
    )
    df = df[~(df["n_mut"].le(min_obs) & df["wsaf_max"].lt(min_wsaf))].drop(
        columns=["n_mut", "wsaf_max"]
    )
    return df


def compute_variant_prevalence(
    variants_df: pd.DataFrame,
    master_df: Optional[pd.DataFrame] = None,
    additional_groups: Optional[list[str]] = None,
) -> pd.DataFrame:
    """
    Compute the prevalence of each mutation in `variants_df`
    """
    if additional_groups is None:
        additional_groups = []

    if additional_groups:
        assert master_df is not None, (
            "master_df must be provided if additional_groups are used"
        )
        assert all(group in master_df.columns for group in additional_groups), (
            "all additional_groups must be columns in master_df"
        )
        variants_df = variants_df.merge(
            master_df[["sample_id", *additional_groups]], on="sample_id", how="left"
        )

    prev_df = (
        variants_df.groupby(
            VARIANTS_GROUP_COLUMNS + additional_groups,
        )
        .agg(
            n_samples=pd.NamedAgg("gt_int", len),
            n_passed=pd.NamedAgg("gt_int", lambda x: sum(x != -1)),
            n_wt=pd.NamedAgg("gt_int", lambda x: sum(x == 0)),
            n_mixed=pd.NamedAgg("gt_int", lambda x: sum(x == 1)),
            n_mut=pd.NamedAgg("gt_int", lambda x: sum(x == 2)),
        )
        .reset_index()
    )

    # Compute frequencies
    prev_df["per_wt"] = 100 * prev_df["n_wt"] / prev_df["n_passed"]
    prev_df["per_mixed"] = 100 * prev_df["n_mixed"] / prev_df["n_passed"]
    prev_df["per_mut"] = 100 * prev_df["n_mut"] / prev_df["n_passed"]

    # Compute prevalence
    prev_df["prevalence"] = prev_df["per_mixed"] + prev_df["per_mut"]

    # Compute prevalence 95% confidence intervals
    low, high = proportion_confint(
        prev_df["n_mut"] + prev_df["n_mixed"],
        prev_df["n_passed"],
        alpha=0.05,
        method="beta",
    )
    prev_df["prevalence_lowci"] = 100 * low
    prev_df["prevalence_highci"] = 100 * high

    return prev_df


def gene_deletions(coverage_df: pd.DataFrame, genes: list[str]) -> pd.DataFrame:
    """
    Analyze gene deletions based on coverage data.

    Assumes `coverage_df` contains columns:
    - 'sample_id'
    - 'name'
    - 'mean_cov'

    A gene is considered deleted in a sample if its mean coverage is below a threshold.
    This is a very simple heuristic and might be improved in the future.

    Returns a DataFrame with deletion prevalence per gene.
    """
    # How many of the other amplolicons must be covered to consider the sample for deletion analysis
    AMPLICONS_QC_CUTOFF = 0.8
    # Minimum coverage to consider a gene deleted
    DELETION_COVERAGE_THRESHOLD = 5

    coverage_df = coverage_df[coverage_df["sample_type"] == "field"].copy()
    coverage_df["gene"] = coverage_df["name"].str.split("-").str[0]

    # QC: only analyze samples with sufficient control amplicon coverage
    analysis_samples = (
        coverage_df[~coverage_df["gene"].isin(genes)]
        .groupby(["expt_name", "barcode"])
        .agg(
            n_ctrl_amplicons=pd.NamedAgg("name", "nunique"),
            n_passing_ctrl_amplicons=pd.NamedAgg("passing", "sum"),
        )
    )

    analysis_set = coverage_df.merge(
        analysis_samples, on=["expt_name", "barcode"], how="left"
    )
    analysis_set = analysis_set[
        analysis_set["n_passing_ctrl_amplicons"]
        >= AMPLICONS_QC_CUTOFF * analysis_set["n_ctrl_amplicons"]
    ]

    # QC: don't analyze amplicons that did not pass negative control
    analysis_set = analysis_set[~analysis_set["fail_contam_abs"]]

    # Determine deletions
    analysis_set["is_deleted"] = analysis_set["mean_cov"] < DELETION_COVERAGE_THRESHOLD

    # consider a gene deleted if all amplicons that belong to the gene are deleted
    result = (
        analysis_set.groupby(["expt_name", "barcode", "sample_id", "gene"])
        .agg(
            is_deleted=pd.NamedAgg("is_deleted", "all"),
        )
        .reset_index()
    )

    # consider a gene deleted, if it is deleted for all replicates of a sample
    result = (
        result.groupby(["sample_id", "gene"])
        .agg(
            is_deleted=pd.NamedAgg("is_deleted", "all"),
            n_deleted=pd.NamedAgg("is_deleted", "sum"),
            n_replicates=pd.NamedAgg("barcode", "count"),
        )
        .reset_index()
    )

    return result[result["gene"].isin(genes)]


def gene_deletion_prevalence_by(
    gene_deletions_df: pd.DataFrame, master_df: pd.DataFrame, fields: list[str]
) -> pd.DataFrame:
    """
    Compute the prevalence of gene deletions in `gene_deletions_df`
    stratified by columns in `fields`.
    """
    gene_deletions_df = gene_deletions_df.merge(
        master_df[["sample_id", *fields]], on="sample_id", how="left"
    )

    prev_df = (
        gene_deletions_df.groupby(["gene", *fields])
        .agg(
            n_samples=pd.NamedAgg("is_deleted", len),
            n_passed=pd.NamedAgg("is_deleted", lambda x: sum(x.notnull())),
            n_deleted=pd.NamedAgg("is_deleted", lambda x: sum(x)),
        )
        .reset_index()
    )

    # Compute prevalence
    prev_df["prevalence"] = 100 * prev_df["n_deleted"] / prev_df["n_passed"]

    # Compute prevalence 95% confidence intervals
    low, high = proportion_confint(
        prev_df["n_deleted"],
        prev_df["n_passed"],
        alpha=0.05,
        method="beta",
    )
    prev_df["prevalence_lowci"] = 100 * low
    prev_df["prevalence_highci"] = 100 * high

    return prev_df
