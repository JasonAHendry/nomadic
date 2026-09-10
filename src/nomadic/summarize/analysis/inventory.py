"""Inventory and throughput calculations for a set of experiments.
Inventory is a dataframe with the following columns:
- expt_name: The name of the experiment
- barcode: The barcode of the sample
- sample_id: The sample ID
- sample_type: The type of the sample (field, pos, neg)
- status: The status of the sample (included, excluded, control)
"""

from dataclasses import dataclass
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd

from nomadic.util.exceptions import UserInputError
from nomadic.util.experiment import ExperimentOutputs

REQUIRED_INVENTORY_COLUMNS = ["expt_name", "barcode", "sample_id", "sample_type"]


def create_inventory_df(expts: list[ExperimentOutputs]) -> pd.DataFrame:
    """Create a dataframe containing the inventory of all experiments.

    The dataframe will contain the following columns:
    - expt_name: The name of the experiment
    - barcode: The barcode of the sample
    - sample_id: The sample ID
    - sample_type: The type of the sample (field, pos, neg)
    """
    inventory_df = pd.concat(
        [expt.metadata[REQUIRED_INVENTORY_COLUMNS] for expt in expts]
    ).reset_index(drop=True)
    inventory_df["sample_id"] = inventory_df["sample_id"].astype(str).str.strip()
    # Check for duplicates
    if inventory_df.duplicated(subset=["expt_name", "barcode"]).any():
        raise UserInputError(
            "Duplicate item found in inventory dataframe. "
            "This should not happen, please check your metadata files."
            f" Duplicates: {inventory_df[inventory_df.duplicated(subset=['expt_name', 'barcode'], keep=False)]}"
        )
    return inventory_df


def add_inventory_status(
    inventory_df: pd.DataFrame, master_metadata: pd.DataFrame
) -> pd.DataFrame:
    """
    Adds a column status with either control (for controls), or included/excluded for field samples
    """
    field_samples = inventory_df.query("sample_type == 'field'")
    excluded_samples = (
        field_samples.loc[
            ~field_samples["sample_id"].isin(master_metadata["sample_id"]),
            "sample_id",
        ]
        .unique()
        .tolist()
    )
    # Mark excluded/included samples
    inventory_df = inventory_df.assign(
        status=np.where(
            inventory_df["sample_id"].isin(excluded_samples), "excluded", "included"
        )
    )
    # set all controls
    inventory_df = inventory_df.assign(
        status=np.where(
            inventory_df["sample_type"].isin(["pos", "neg"]),
            "control",
            inventory_df["status"],
        )
    )
    return inventory_df


def drop_excluded_samples(inventory_df: pd.DataFrame) -> pd.DataFrame:
    """Drop all excluded samples and the full experiment including controls if no sample is included in an experiment"""
    keep_mask = (
        inventory_df["status"]
        .eq("included")
        .groupby(inventory_df["expt_name"])
        .transform("any")
    ) & inventory_df["status"].ne("excluded")
    return inventory_df[keep_mask]


def n_excluded_samples(inventory_df: pd.DataFrame) -> int:
    """Return the number of excluded samples in the inventory dataframe"""
    return inventory_df.loc[inventory_df["status"] == "excluded", "sample_id"].nunique()


def n_field_samples(inventory_df: pd.DataFrame) -> int:
    """Return the number of field samples in the inventory dataframe"""
    return inventory_df.loc[
        inventory_df["sample_type"] == "field", "sample_id"
    ].nunique()


def experiments_in_inventory(
    inventory_df: pd.DataFrame, expt_dirs: Optional[list[Path]]
) -> list[Path]:
    """Return a list of experiments in the inventory dataframe"""
    expt_names_in_inventory = inventory_df["expt_name"].unique().tolist()
    if expt_dirs is None:
        return expt_names_in_inventory

    return [
        expt_dir for expt_dir in expt_dirs if expt_dir.name in expt_names_in_inventory
    ]


@dataclass
class Throughput:
    """Class to store throughput information"""

    n_expts: int
    n_expts_included: int
    n_pos: int
    n_neg: int
    n_field_total: int
    n_field_unique: int


def compute_throughput(
    inventory_df: pd.DataFrame, add_unique: bool = True
) -> tuple[Throughput, pd.DataFrame]:
    """
    Compute throughput crosstable

    Also add information about uniqueness

    """
    inventory_df = inventory_df.assign(
        sample_type=np.where(
            inventory_df["status"] == "control",
            inventory_df["sample_type"],
            inventory_df["sample_type"] + "_" + inventory_df["status"],
        )
    )
    throughput_df = pd.crosstab(
        inventory_df["expt_name"], inventory_df["sample_type"], margins=True
    )

    inventory_included_expts_df = inventory_df.groupby("expt_name").filter(
        lambda x: (x["status"] == "included").sum() > 0
    )

    throughtput_included_expts_df = pd.crosstab(
        inventory_included_expts_df["expt_name"],
        inventory_included_expts_df["sample_type"],
        margins=True,
    )

    throughput_df.loc["included_expts"] = throughtput_included_expts_df.loc["All"]

    if add_unique:
        um = inventory_df.drop_duplicates("sample_id").query("status == 'included'")
        throughput_df["field_unique"] = pd.crosstab(
            um["expt_name"], um["sample_type"], margins=True
        )["field_included"]

    throughput_df.loc["included_expts", "field_unique"] = throughput_df.loc[
        "All", "field_unique"
    ]

    # Ensure all expected rows are present
    throughput_df = (
        throughput_df.reindex(
            columns=["All", "pos", "neg", "field_included", "field_unique"],
            fill_value=0,
        )
        .fillna(0)
        .astype(int)
    )

    return Throughput(
        n_expts=int(throughput_df.shape[0]),
        n_expts_included=int(throughput_df.query("field_included > 0").shape[0] - 2),
        n_pos=int(throughput_df.loc["included_expts", "pos"]),
        n_neg=int(throughput_df.loc["included_expts", "neg"]),
        n_field_total=int(throughput_df.loc["All", "field_included"]),
        n_field_unique=int(throughput_df.loc["All", "field_unique"]),
    ), throughput_df
