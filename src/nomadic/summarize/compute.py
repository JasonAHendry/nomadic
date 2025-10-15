import pandas as pd

from statsmodels.stats.proportion import proportion_confint

variants_group_columns = [
    "gene",
    "chrom",
    "pos",
    "ref",
    "alt",
    "aa_change",
    "mut_type",
    "mutation",
]


def filter_false_positives(variants_df: pd.DataFrame):
    mut = variants_df[variants_df["gt_int"].isin([1,2])]
    df = variants_df.merge(mut.groupby(variants_group_columns).agg(n_mut=pd.NamedAgg("gt_int", len), wsaf_max=pd.NamedAgg("wsaf", "max")), on=variants_group_columns)
    df = df[~(df["n_mut"].le(1) & df["wsaf_max"].lt(0.1))].drop(columns=["n_mut", "wsaf_max"])
    return df


def compute_variant_prevalence(variants_df: str) -> pd.DataFrame:
    """
    Compute the prevalence of each mutation in `variants_df`

    Assumes columns several columns exist; compute across all samples
    in data;

    """

    prev_df = (
        variants_df.groupby(
            variants_group_columns,
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


def compute_variant_prevalence_per(
    variants_df, master_df, fields: list[str]
) -> pd.DataFrame:
    """
    Compute the prevalence of each mutation in `variants_df`

    Assumes columns several columns exist; compute across all samples
    in data;

    """
    variants_df = variants_df.merge(
        master_df[["sample_id", *fields]], on="sample_id", how="left"
    )

    prev_df = (
        variants_df.groupby([*variants_group_columns, *fields])
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
