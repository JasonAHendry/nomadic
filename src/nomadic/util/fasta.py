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
