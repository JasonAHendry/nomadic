import warnings
from typing import Optional

import pandas as pd

from nomadic.util.summary_settings import Settings, get_master_columns_mapping


def load_master_metadata(metadata_path, *, settings: Settings) -> pd.DataFrame:
    mapping = get_master_columns_mapping(settings)
    sample_id_column = next(
        (
            source
            for source, destination in mapping.items()
            if destination == "sample_id"
        ),
        "sample_id",
    )
    return pd.read_csv(metadata_path, dtype={sample_id_column: "str"}).rename(
        columns=mapping
    )


def master_metadata_from_expts(expts, *, shared_columns: list[str]) -> pd.DataFrame:
    """Create a master metadata dataframe from a list of experiments and shared columns"""
    shared_columns = ["sample_id"] + list(shared_columns)
    master_metadata = pd.concat(
        [
            expt.metadata.query("sample_type == 'field'")[shared_columns]
            for expt in expts
        ]
    )

    # Drop safe duplicates
    master_metadata = master_metadata.drop_duplicates()

    duplicated_sample_ids = master_metadata[
        master_metadata.duplicated(subset=["sample_id"], keep=False)
    ]

    # Log the duplicated rows as a warning
    if not duplicated_sample_ids.empty:
        warnings.warn(
            f"Duplicated sample_ids found with different metadata across experiments, will keep the first occurrence. Duplicates:\n{duplicated_sample_ids}",
            category=UserWarning,
        )

    # Note, problematic if same sample ID has different metadata across experiments
    return master_metadata.drop_duplicates(subset=["sample_id"])


def get_shared_metadata_columns(
    metadata_dfs: list[pd.DataFrame],
    fixed_columns: Optional[list[str]] = None,
) -> list[str]:
    """Get metadata columns that are shared acrossa all experiments"""
    if fixed_columns is None:
        fixed_columns = ["expt_name", "barcode", "sample_id", "sample_type"]

    shared_columns = set(metadata_dfs[0].columns)
    for df in metadata_dfs[1:]:
        shared_columns.intersection_update(df.columns)
    shared_columns.difference_update(fixed_columns)  # why am I doing this?
    return list(shared_columns)


def normalize_sample_id(df: pd.DataFrame) -> pd.DataFrame:
    """Normalize the sample_id column in a dataframe to be str and stripped of whitespace"""
    return df.assign(sample_id=df["sample_id"].astype(str).str.strip())


METADATA_COLUMN_PREFIX = "metadata__"


def prefix_metadata_columns(df: pd.DataFrame) -> pd.DataFrame:
    """Prefix metadata columns to avoid collisions"""
    return df.rename(
        columns={
            col: f"{METADATA_COLUMN_PREFIX}{col}"
            for col in df.columns
            if col != "sample_id"
        }
    )


def remove_metadata_prefix(df: pd.DataFrame) -> pd.DataFrame:
    """Remove metadata prefix from columns"""
    return df.rename(
        columns={
            col: col.replace(METADATA_COLUMN_PREFIX, "")
            for col in df.columns
            if col.startswith(METADATA_COLUMN_PREFIX)
        }
    )


def normalize_metadata(df: pd.DataFrame) -> pd.DataFrame:
    """Normalize the metadata dataframe to have consistent sample_id and column names"""
    return df.pipe(normalize_sample_id).pipe(prefix_metadata_columns)


def validate_metadata(df: pd.DataFrame) -> pd.DataFrame:
    """Validate that the metadata dataframe has the required columns and no duplicate sample_ids"""
    required_columns = ["sample_id"]
    missing_columns = set(required_columns) - set(df.columns)
    if missing_columns:
        raise ValueError(f"Missing required columns in metadata: {missing_columns}")
    if df["sample_id"].duplicated().any():
        raise ValueError("Duplicate sample_ids found in metadata")
    return df
