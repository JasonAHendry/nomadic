import enum
from collections.abc import Iterable
from dataclasses import dataclass
from enum import StrEnum, auto
from pathlib import Path

import numpy as np
import pandas as pd

from nomadic.util.experiment import get_summary_files


def create_region_coverage_df(
    expt_dirs: Iterable[Path], inventory_df: pd.DataFrame
) -> pd.DataFrame:
    """
    Here we load a consolidated region coverage dataframe and include information required
    for quality control
    """

    # Load coverage data
    bed_dfs = []
    for expt_dir in expt_dirs:
        bed_csv = get_summary_files(expt_dir).region_coverage

        bed_df = pd.read_csv(bed_csv)
        bed_df.insert(0, "expt_name", expt_dir.name)
        bed_df.query("barcode != 'unclassified'", inplace=True)
        if "sample_id" in bed_df.columns:
            bed_df.drop(columns=["sample_id"], inplace=True)

        bed_df = pd.merge(
            left=inventory_df[["expt_name", "barcode", "sample_id", "sample_type"]],
            right=bed_df,
            on=["expt_name", "barcode"],
            how="inner",  # ensure we only take samples that are in the master metadata
            validate="1:m",
        )
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

    coverage_df = pd.merge(
        left=concat_df[
            [
                "expt_name",
                "barcode",
                "sample_id",
                "sample_type",
                "chrom",
                "name",
                "mean_cov",
            ]
        ],  # sample ID, will want it at some point
        right=neg_df,
        on=["expt_name", "name"],
        how="left",
        validate="m:1",
    )

    # Checkk that all expected columns are present
    expected_columns = [
        "expt_name",
        "barcode",
        "sample_id",
        "sample_type",
        "chrom",
        "name",
        "mean_cov",
        "mean_cov_neg",
    ]
    if not all(col in coverage_df.columns for col in expected_columns):
        raise ValueError(
            f"Coverage dataframe is missing expected columns. Expected: {expected_columns}, got: {coverage_df.columns.tolist()}"
        )

    return coverage_df


def add_quality_control_columns(
    region_coverage_df: pd.DataFrame, *, min_coverage: int = 50, max_contam: float = 0.1
) -> pd.DataFrame:
    """
    Return a dataframe with columns added evaluating whether samples have passed quality control.

    Quality control is based on two criteria:
    1. Minimum coverage: If the mean coverage of a sample is below `min_coverage`, it is considered to have failed quality control.
    2. Contamination: If the mean coverage of the negative control is relative to the sample coverage more than `max_contam`, it is considered to have failed quality control.
    3. Absolute contamination: If the mean coverage of the negative control is more than `min_coverage`, it is considered to have failed quality control.

    Input dataframe must contain the following columns:
    - mean_cov: The mean coverage of the sample
    - mean_cov_neg: The mean coverage of the negative control

    Output dataframe will contain the following columns:
    - fail_lowcov: True if the sample has failed quality control due to low coverage
    - fail_contam_rel: True if the sample has failed quality control due to relative contamination
    - fail_contam_abs: True if the sample has failed quality control due to absolute contamination
    - fail_contam: True if the sample has failed quality control due to either relative or absolute contamination
    - passing: True if the sample has passed quality control
    """

    # Do we have enough coverage?
    fail_lowcov = region_coverage_df["mean_cov"] < min_coverage

    # Check if coverage of negative control exceeds `max_contam`
    fail_contam_rel = (
        region_coverage_df["mean_cov_neg"] / (region_coverage_df["mean_cov"] + 0.01)
    ) >= max_contam
    fail_contam_abs = region_coverage_df["mean_cov_neg"] >= min_coverage

    fail_contam = (
        (fail_contam_rel & ~fail_lowcov) | fail_contam_abs
    )  # If already failed low coverage, don't consider contamination, unless absolute threshold is exceeded

    # Finally, define passing
    passing = ~fail_contam & ~fail_lowcov

    return region_coverage_df.assign(
        fail_lowcov=fail_lowcov,
        fail_contam_rel=fail_contam_rel,
        fail_contam_abs=fail_contam_abs,
        fail_contam=fail_contam,
        passing=passing,
    )


@dataclass
class FieldCoverageSummary:
    n: int
    n_lowcov: int
    n_contam: int
    n_pass: int


def compute_field_coverage_summary(coverage_df: pd.DataFrame) -> FieldCoverageSummary:
    """
    Compute a summary of the field coverage dataframe, including the number of samples,
    the number of samples that failed quality control due to low coverage, the number of samples that failed quality control due to contamination, and the number of samples that passed quality control.
    """

    field_coverage_df = coverage_df.query("sample_type == 'field'")
    n = field_coverage_df.shape[0]
    n_lowcov = field_coverage_df["fail_lowcov"].sum()
    n_contam = field_coverage_df["fail_contam"].sum()
    n_pass = field_coverage_df["passing"].sum()

    return FieldCoverageSummary(
        n=n, n_lowcov=n_lowcov, n_contam=n_contam, n_pass=n_pass
    )


class QcStatus(StrEnum):
    PASS = auto()
    LOWCOV = auto()
    CONTAM = auto()
    DUPLICATE = auto()
    CONTROL = auto()


def add_qc_status(region_qc_df: pd.DataFrame) -> pd.DataFrame:
    """
    Adds a status to each replicate/amplicon to see which ones passed QC
    and if they didn't, why.
    """
    is_control = region_qc_df["sample_type"].isin(["pos", "neg"])
    fail_contam = region_qc_df["fail_contam"]
    fail_lowcov = region_qc_df["fail_lowcov"]

    status = np.select(
        [
            is_control,
            fail_contam & fail_lowcov,
            fail_contam,
            fail_lowcov,
        ],
        [
            QcStatus.CONTROL,
            f"{QcStatus.CONTAM};{QcStatus.LOWCOV}",
            QcStatus.CONTAM,
            QcStatus.LOWCOV,
        ],
        default=QcStatus.PASS,
    )

    return region_qc_df.assign(status=status)


def mark_duplicates(region_qc_df: pd.DataFrame) -> pd.DataFrame:
    """
    Mark all field samples as duplicates, if there is a better covered replicate for the same amplicon
    Replicates marked as duplicate will not be used for prevalance evaluation.
    """
    field_mask = region_qc_df["sample_type"].eq("field")
    field = region_qc_df.loc[field_mask, ["sample_id", "name", "status", "mean_cov"]]

    # Sort by passing/non passing and then mean cov
    ordered = field.assign(_passing=field["status"].eq(QcStatus.PASS)).sort_values(
        ["sample_id", "name", "_passing", "mean_cov"],
        ascending=[True, True, False, False],
    )

    # keep the first
    duplicate_idx = ordered.index[
        ordered.duplicated(["sample_id", "name"], keep="first")
    ]

    status = region_qc_df["status"].copy()
    status.loc[duplicate_idx] = (
        status.loc[duplicate_idx].astype(str) + f";{QcStatus.DUPLICATE}"
    )

    return region_qc_df.assign(status=status)


def add_quality_control_status_column(df: pd.DataFrame) -> pd.DataFrame:
    """
    Add a QC status column and return the result without mutating the input

    Note:
    - When we mark duplicates; we do it on an AMPLICON x SAMPLE level;
    not on a per-sample level. So we could take amplicons from separate
    samples to get the best data for that sample.

    """
    return df.pipe(add_qc_status).pipe(mark_duplicates)


class Status(enum.Enum):
    PASSING = "passing"
    FAILING = "failing"
    NOT_SEQUENCED = "not_sequenced"


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


def samples_qc(
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
        lambda row: (
            Status.PASSING.value
            if row["n_passing"] > 0
            else Status.FAILING.value
            if row["n_replicates"] > 0
            else Status.NOT_SEQUENCED.value
        ),
        axis=1,
    )
    samples_summary_df.sort_values(
        by=["n_passing", "n_replicates", "sample_id"],
        inplace=True,
        ascending=[False, False, True],
    )

    return samples_summary_df


def replicates_amplicon_qc(coverage_df):
    return coverage_df.query("sample_type == 'field'")


def amplicons_qc_summary(master_metadata, replicates_amplicon_qc_df):
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
        lambda row: (
            Status.PASSING.value
            if row["n_passing"] > 0
            else Status.FAILING.value
            if row["n_replicates"] > 0
            else Status.NOT_SEQUENCED.value
        ),
        axis=1,
    )
    samples_by_amplicons_summary_df.sort_values(
        by=["n_passing", "n_replicates", "sample_id"],
        inplace=True,
        ascending=[False, False, True],
    )

    return samples_by_amplicons_summary_df


def experiment_qc_summary(amplicon_qc_df: pd.DataFrame):
    return (
        amplicon_qc_df.groupby(["expt_name", "name"])
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
