import re
from typing import Optional
import warnings
import pandas as pd
from typing import List
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


# --------------------------------------------------------------------------------
# Check validity of various columns in the metadata
#
# --------------------------------------------------------------------------------


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


def check_sample_type_format(sample_type: str, try_to_fix: bool = False) -> str:
    """
    Check that the format of the `sample_type` column is correct, and optionally
    try to fix if it is not

    """

    # Settings
    EXPECTED = ["field", "pos", "neg"]
    KNOWN_POS_CONTROLS = [
        "3D7",
        "Dd2",
        "HB3",
        "GB4",
        "7G8",
        "NF54",
        "FCR3",
    ]  # some smaller ones could be field substrings, e.g. K1, W2

    if pd.isna(sample_type):
        raise MetadataFormatError(
            "Missing information in the 'sample_type' column. Please ensure it is complete."
        )

    if not isinstance(sample_type, str):
        sample_type = str(sample_type)

    sample_type = sample_type.strip()  # safe to do this in all cases
    if sample_type in EXPECTED:
        return sample_type

    if not try_to_fix:
        raise MetadataFormatError(
            f"Found a sample with type '{sample_type}' which is invalid. Please assign one of: {', '.join(EXPECTED)}."
        )

    # Raise a warning and proceed if we are fixing.
    # warnings.warn(
    #     f"Found a sample with type '{sample_type}' which is invalid. Trying to fix..."
    # )

    for e in EXPECTED:
        if sample_type.lower() == e:  # capitalisation issue
            return e
        if sample_type.lower().startswith(e):  # added something to the end
            return e

    for k in KNOWN_POS_CONTROLS:
        if sample_type.lower() == k.lower():
            return "pos"

    for (
        e
    ) in EXPECTED:  # this can be dangerous, only do if other attempts haven't worked
        if sample_type.lower().startswith(e[0]):
            return e

    # Raise if couldn't fix
    raise MetadataFormatError(
        f"Found a sample with type '{sample_type}' which is invalid. Please assign one of: {', '.join(EXPECTED)}."
    )


# --------------------------------------------------------------------------------
# Class(es) to parse metadata table(s) for various use cases, e.g. realtime
# analysis or summarizing
# --------------------------------------------------------------------------------


class MetadataTableParser:
    """
    Parse the `metadata_csv` table, and make sure that it is formatted
    correctly. If not formatted correctly, try and fix it if possible.

    """

    REQUIRED_COLUMNS = ["barcode", "sample_id", "sample_type"]
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
        "sample_type": [
            "sampletype",
            "sample-type",
            "sample type",
            "sampletypes",
            "sample-types",
            "sample types",
        ],
    }

    def __init__(self, metadata_csv: str, include_unclassified: bool = True):
        """
        Load and sanity check the metadata table

        """

        self.csv = metadata_csv
        self.df = pd.read_csv(self.csv, delimiter=get_csv_delimiter(self.csv))

        self._correct_columns()
        self._check_entries_unique()
        self._correct_all_barcodes()

        self.barcodes = self.df["barcode"].tolist()
        self.sample_ids_df = self.df[["barcode", "sample_id"]].set_index("barcode")

        if include_unclassified:
            self.sample_ids_df.loc["unclassified"] = ["unclassified"]
            self.barcodes.append("unclassified")

    def get_sample_id(self, barcode: str) -> Optional[str]:
        if barcode == "unclassified":
            return barcode
        metadata = self.sample_ids_df
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

    def _correct_all_barcodes(self):
        self.df["barcode"] = [correct_barcode_format(b) for b in self.df["barcode"]]


class ExtendedMetadataTableParser(MetadataTableParser):
    """
    Add requirement for sample type, and parse it

    # Should also check we have positive and negative controls in each experiment

    """

    def __init__(self, metadata_csv: str, include_unclassified: bool = True):
        super().__init__(metadata_csv, include_unclassified)

        self._check_sample_type()

    def _check_sample_type(self):
        if "sample_type" not in self.df.columns:
            raise MetadataFormatError(
                "Metadata is missing the column 'sample_type'. "
                "Please create this column and populate it with"
                " 'field' (field sample), 'pos' (positive control) or 'neg' (negative control)"
                " for each sample."
            )
        self.df["sample_type"] = [
            check_sample_type_format(s, try_to_fix=True) for s in self.df["sample_type"]
        ]

        sample_type_counts = self.df.sample_type.value_counts().to_dict()
        if "neg" not in sample_type_counts:
            raise MetadataFormatError("No negative control found for experiment!")
        # if "pos" not in sample_type_counts:
        #     raise MetadataFormatError("No positive control found for experiment!")
