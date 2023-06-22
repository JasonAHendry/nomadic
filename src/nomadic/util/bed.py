import pandas as pd
from dataclasses import dataclass

from nomadic.util.exceptions import BEDFormatError


@dataclass
class BEDRecord:
    """
    Store information from a single line in a BED file

    """

    chrom: str
    start: int
    end: int
    name: str

    # TODO: I am not sure of a simple and general way to handle this
    # for now, don't allow
    # other_fields: List[str|int]

    def __post_init__(self):
        """
        A few simple sanity checks

        """

        if not (self.start > 0) and (self.end > 0):
            raise BEDFormatError(
                "BED file `start` and `end` positions must be all positive."
            )

        if not self.start <= self.end:
            raise BEDFormatError(
                "BED file `start` position must be less or equal to `end`."
            )

    @classmethod
    def from_line(cls, line):
        """
        Create a BED record from a file line

        """

        # Parse string
        fields = line.strip().split("\t")
        if not len(fields) == 4:
            raise BEDFormatError(
                "Region BED file must have four tab-separated columns, "
                f"but finding {len(fields)} columns after parsing."
            )

        # Convert fields to required data types
        chrom = str(fields[0])
        start = int(fields[1])
        end = int(fields[2])
        name = str(fields[3])

        return cls(chrom, start, end, name)


def load_bed_as_dataframe(bed_path: str) -> pd.DataFrame:
    """
    Load a BED file from `bed_path` into a pandas DataFrame

    """

    records = []
    with open(bed_path, "r") as bed:
        for line in bed:
            if line.startswith("#"):
                continue
            if not line.strip():  # handle blank lines, sometimes at end
                continue
            records.append(BEDRecord.from_line(line))
    
    if not records:
        raise BEDFormatError(f"No BED records were found in file {bed_path}.")

    return pd.DataFrame(records)
