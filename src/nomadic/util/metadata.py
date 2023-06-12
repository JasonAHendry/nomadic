import pandas as pd
from .exceptions import MetadataFormatError


class MetadataTableParser:
    """
    Parse the `metadata_csv` table, and make sure that it is formatted
    correctly

    """

    REQUIRED_COLUMNS = ["barcode", "sample_id"]

    def __init__(self, metadata_csv: str, include_unclassified: bool = True):
        """
        Load and sanity check the metadata table

        """

        self.csv = metadata_csv
        self.df = pd.read_csv(self.csv)

        self._check_for_columns()
        self._check_entries_unique()

        self.barcodes = self.df["barcode"].tolist()
        if include_unclassified:
            self.barcodes.append("unclassified")

    def _check_for_columns(self):
        """
        Check the correct columns are present

        """

        for c in self.REQUIRED_COLUMNS:
            if not c in self.df.columns:
                raise MetadataFormatError(f"Metadata must contain column called {c}!")

    def _check_entries_unique(self):
        """
        Check entires of the required columns are unique

        TODO: this will also disallow missing?

        """

        for c in self.REQUIRED_COLUMNS:
            all_entries = self.df[c].tolist()
            observed_entries = []
            for entry in all_entries:
                if entry in observed_entries:
                    raise MetadataFormatError(
                        f"Column {c} must contain only unique entires, but {entry} is duplicated."
                    )
                observed_entries.append(entry)
