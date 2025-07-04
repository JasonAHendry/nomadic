import os
import subprocess
from typing import Dict


def load_fasta_as_dict(fasta_path: str) -> Dict[str, str]:
    """
    Load a `.fasta` file as a dictionary.

    Args:
        fasta_path (str): Path to the `.fasta` file to be loaded.

    Returns:
        Dict[str, str]: A dictionary with headers as keys and sequences as values.
    """
    dt = {}
    with open(fasta_path, "r") as fasta:
        header = None
        seq = ""
        line = fasta.readline().rstrip()

        while line:
            if line.startswith(">"):
                if header is not None:
                    dt[header] = seq
                header = line
                seq = ""
            else:
                seq += line
            line = fasta.readline().rstrip()

        if header is not None:  # Ensure the last sequence is added to the dictionary
            dt[header] = seq

    return dt


def write_fasta_from_dict(input_dt: Dict[str, str], output_fasta: str) -> None:
    """
    Write a `.fasta` file to `output_fasta` from an input dictionary `input_dt`.

    Args:
        input_dt (Dict[str, str]): Dictionary with headers as keys and sequences as values.
        output_fasta (str): Path where the output `.fasta` file will be written.

    """
    with open(output_fasta, "w") as fasta:
        for header, seq in input_dt.items():
            if not header.startswith(">"):
                header = f">{header}"
            fasta.write(f"{header}\n")
            fasta.write(f"{seq}\n")


def find_lowcomplexity_intervals(fasta_path: str, bed_path: str) -> None:
    """
    Find low-complexity intervals in a FASTA file `fasta_path` using
    the sdust algorithm

    """

    if not os.path.exists(fasta_path):
        raise FileNotFoundError(f"No FASTA file present at {fasta_path}!")

    SUFFIXES = [".fasta", ".fa", ".fna"]
    fasta_suffix = [suffix for suffix in SUFFIXES if fasta_path.endswith(suffix)]
    if not fasta_suffix:
        raise ValueError(
            f"Input file must be FASTA, with one of these suffixes: {', '.join(SUFFIXES)}."
        )

    cmd = f"sdust {fasta_path} > {bed_path}"
    subprocess.run(cmd, check=True, shell=True)
