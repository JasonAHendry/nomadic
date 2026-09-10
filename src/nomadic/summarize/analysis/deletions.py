import pandas as pd
from statsmodels.stats.proportion import proportion_confint


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
    AMPLICONS_QC_CUTOFF = 0.7
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
    )[
        [
            "expt_name",
            "barcode",
            "sample_id",
            "sample_type",
            "name",
            "mean_cov",
            "mean_cov_neg",
            "fail_contam_abs",
            "gene",
            "n_ctrl_amplicons",
            "n_passing_ctrl_amplicons",
        ]
    ]

    analysis_set = analysis_set[analysis_set["gene"].isin(genes)]

    # QC: only analyze samples with sufficient control amplicon coverage
    analysis_set["pass_ctrl_amplicons"] = (
        analysis_set["n_passing_ctrl_amplicons"]
        >= AMPLICONS_QC_CUTOFF * analysis_set["n_ctrl_amplicons"]
    )

    # analysis_set.to_csv("analysis_set_gene_deletions.csv", index=False)  # For debugging

    analysis_set = analysis_set[analysis_set["pass_ctrl_amplicons"]]

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

    return result


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
