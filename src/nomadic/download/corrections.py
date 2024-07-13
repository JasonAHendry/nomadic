from typing import List
from dataclasses import dataclass
from nomadic.util.fasta import load_fasta_as_dict, write_fasta_from_dict


@dataclass
class NucleotideChange:
    name: str
    chrom: str
    position: int
    before: str
    after: str


# I want to revert G>C at this position so that 3D7 carries
# an alanine (WT) rather than a glycine (Sulfadoxine R associated)
DHPS = NucleotideChange(
    name="dhps-A437G", chrom="Pf3D7_08_v3", position=549685, before="G", after="C"
)


def update_reference_genome(fasta_path: str, mutations: List[NucleotideChange]) -> None:
    """
    Update a reference genome by reverting mutations; do this IN PLACE

    """
    dt = load_fasta_as_dict(fasta_path)

    # Iterate over mutations
    for mutation in mutations:
        print(f"Correcting: {mutation}...")
        chrom_str = ">" + mutation.chrom

        found_chrom = False
        for key, seq in dt.items():
            if key.startswith(chrom_str):
                found_chrom = True
                break

        if not found_chrom:
            sep = "\n"
            raise ValueError(
                f"Cannot find chromosome {mutation.chrom} in: {sep.join(dt.keys())}!"
            )

        # Replace
        L = len(seq)
        dt[key] = (
            seq[: (mutation.position - 1)] + mutation.after + seq[mutation.position :]
        )
        assert (
            len(dt[key]) == L
        ), f"Something wrong with replacement, sequence length changed {L} != {len(dt[key])}."
        print("Done.")

    print("Overwriting FASTA file...")
    write_fasta_from_dict(dt, fasta_path)
    print("Done.")
    
