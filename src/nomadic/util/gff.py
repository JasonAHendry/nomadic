import gzip
import pandas as pd
from dataclasses import dataclass


def load_gff(gff_path: str) -> pd.DataFrame:
    """Load a gene feature format (.gff) file into a pandas DataFrame"""
    
    if gff_path is None:
        raise FileNotFoundError("No GFF file exists to load.")

    # Define fields
    @dataclass
    class gffEntry:
        seqname: str
        source: str
        feature: str
        start: int
        end: int
        score: str
        strand: str
        frame: str
        attribute: str

    # Prepare to load, handle compression
    if gff_path.endswith(".gz"):
        binary_gff = True
        open_gff = gzip.open
    else:
        binary_gff = False
        open_gff = open

    # Open gff
    entries = []
    with open_gff(gff_path) as gff:

        # Iterate over rows
        for line in gff:
            # Decode as necessary
            if binary_gff:
                line = line.decode()

            # Skip if info
            if line.startswith("#"):
                continue

            # Extract gff fields
            fields = line.strip().split("\t")
            entry = gffEntry(*fields)

            # Store
            entries.append(entry)
    
    # Coerce data types
    gff_df = pd.DataFrame(entries)
    gff_df["start"] = gff_df["start"].astype("int")
    gff_df["end"] = gff_df["end"].astype("int")

    return gff_df