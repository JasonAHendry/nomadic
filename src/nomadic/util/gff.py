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


def parse_attributes(attribute_str: str) -> dict:
    """
    Parse the attributes field of a GFF entry into a dictionary.

    The attributes are expected to be in the format 'key1=value;key2=value2;...'.
    """
    attributes = {}
    if attribute_str:
        for attr in attribute_str.split(";"):
            if "=" in attr:
                key, value = attr.split("=", 1)
                attributes[key.strip()] = value.strip()
            else:
                attributes[attr.strip()] = None  # Handle keys without values
    return attributes


def write_attributes(attributes: dict) -> str:
    """
    Convert a dictionary of attributes back into the GFF attribute string format.

    The output will be in the format 'key1=value;key2=value2;...'.
    """
    return ";".join(f"{k}={v}" for k, v in attributes.items() if v is not None)


def replace_attribute_keys(attributes: dict, replacements: dict) -> dict:
    """
    Replace keys in the attributes dictionary based on a mapping.

    :param attributes: The original attributes dictionary.
    :param replacements: A dictionary mapping old keys to new keys.
    :return: A new dictionary with keys replaced according to the mapping.
    """
    return {replacements.get(k, k): v for k, v in attributes.items()}
