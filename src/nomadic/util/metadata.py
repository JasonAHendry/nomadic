import os
import re
import warnings
from typing import List, Optional

import pandas as pd

from .exceptions import MetadataFormatError


def get_csv_delimiter(csv_path: str, delimiters: List[str] = [",", ";", "\t"]):
    """
    Determine which delimiter is being used in a CSV

    I couldn't find a solution just using `pd.read_csv`; using
    a regex sep=[,;] does not work.

    """

    with open(csv_path, "r") as csv:
        for header in csv:
            break

        used = [delimiter for delimiter in delimiters if header.find(delimiter) != -1]

        if not used:
            return None

        if not len(used) == 1:
            raise MetadataFormatError(
                f"Found multiple delimiters ({' '.join(used)}) in header: {header.strip()}."
            )

    return used[0]


def correct_barcode_format(barcode: str, try_to_fix: bool = True) -> str:
    """
    Check that the format of a barcode is as expected, and optionally
    try and fix if it is not

    """

    # Some hard-coded settings
    MAX_BARCODE = 96
    EXPECTED = "barcode[0-9]{2}$"
    EXAMPLE = "barcode01"

    if not isinstance(barcode, str):
        barcode = str(barcode)

    if not re.match(EXPECTED, barcode):
        if not try_to_fix:
            raise MetadataFormatError(
                f"Barcode '{barcode}' has bad format: must conform to '{EXAMPLE}'."
            )

        # Raise a warning
        warnings.warn(
            f"Barcode '{barcode}' has bad format: must conform to '{EXAMPLE}'. Trying to fix..."
        )

        nums = re.findall("[0-9]+", barcode)

        if not nums:
            raise MetadataFormatError(
                f"Barcode '{barcode}' has bad format: must conform to '{EXAMPLE}'."
            )

        if len(nums) > 1:
            raise MetadataFormatError(f"Multiple numbers found in barcode: {barcode}.")

        barcode_int = int(nums[0])
        if barcode_int > 96:
            raise MetadataFormatError(
                f"Barcode '{barcode}' exceeds maximum of {MAX_BARCODE}."
            )

        barcode = f"barcode{barcode_int:02d}"

    return barcode


class MetadataTableParser:
    """
    Parse the `metadata_csv` table, and make sure that it is formatted
    correctly. If not formatted correctly, try and fix it if possible.

    """

    REQUIRED_COLUMNS = ["barcode", "sample_id"]
    UNIQUE_COLUMNS = ["barcode"]

    # If the required columns are not found, try these alternative names, case insensitive
    ALTERNATIVE_NAMES = {
        "barcode": ["barcodes"],
        "sample_id": [
            "sample",
            "sampleid",
            "sample-id",
            "sample_id",
            "sampleids",
            "sample-ids",
            "sample_ids",
            "sample id",
            "sample ids",
        ],
    }

    def __init__(self, metadata_path: str, include_unclassified: bool = True):
        """
        Load and sanity check the metadata table

        """

        self.path = metadata_path
        if os.path.splitext(self.path)[1].lower() == ".xlsx":
            data = pd.read_excel(self.path, sheet_name="nomadic")
            data.dropna(how="all", inplace=True)
            self.df = data
        else:
            self.df = pd.read_csv(self.path, delimiter=get_csv_delimiter(self.path))

        self._correct_columns()
        self._check_entries_unique()
        self._correct_all_barcodes()

        self.barcodes = self.df["barcode"].tolist()
        self.required_metadata = self.df[self.REQUIRED_COLUMNS].set_index("barcode")

        if include_unclassified:
            self.required_metadata.loc["unclassified"] = ["unclassified"]
            self.barcodes.append("unclassified")

    def get_sample_id(self, barcode: str) -> Optional[str]:
        if barcode == "unclassified":
            return barcode
        metadata = self.required_metadata
        if barcode not in metadata.index:
            return None
        return metadata.loc[barcode].get("sample_id", None)

    def _correct_columns(self):
        """
        Check the correct columns are present and if not try to find them under different names

        """

        normalized_column_names = [c.strip().lower() for c in self.df.columns]

        for required_column in self.REQUIRED_COLUMNS:
            if required_column not in self.df.columns:
                for alt in [
                    required_column,
                    *self.ALTERNATIVE_NAMES.get(required_column, []),
                ]:
                    if alt in normalized_column_names:
                        column_name = self.df.columns[
                            normalized_column_names.index(alt)
                        ]
                        warnings.warn(
                            f"Using column '{column_name}' as '{required_column}' in metadata CSV."
                        )
                        self.df.rename(
                            columns={column_name: required_column}, inplace=True
                        )
                        break
                else:
                    raise MetadataFormatError(
                        f"Metadata must contain column called {required_column}!"
                    )

    def _check_entries_unique(self):
        """
        Check entires of the required columns are unique

        TODO: this will also disallow missing?

        """

        for c in self.UNIQUE_COLUMNS:
            all_entries = self.df[c].tolist()
            observed_entries = []
            for entry in all_entries:
                if entry in observed_entries:
                    raise MetadataFormatError(
                        f"Column {c} must contain only unique entires, but {entry} is duplicated."
                    )
                observed_entries.append(entry)

    def _correct_all_barcodes(self) -> List[str]:
        self.df["barcode"] = [correct_barcode_format(b) for b in self.df["barcode"]]
