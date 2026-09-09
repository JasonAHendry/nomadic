from dataclasses import dataclass, field
from pathlib import Path

metadata_dir = "metadata"
seq_inventory_dir = "seq_inventory"
qc_dir = "quality_control"
variants_dir = "variants"
vcfs_dir = "vcfs"  # inside variants
gene_deletions_dir = "gene_deletions"

sample_inventory_file = "sample_inventory.csv"
metadata_file = "master_metadata.csv"
throughput_file = "samples.by_experiment.csv"

coverage_file = "coverage.csv"
replicates_qc_file = "replicates_qc.csv"
samples_qc_file = "samples_qc.csv"
samples_by_amplicon_qc_file = "samples_amplicons_qc.csv"
experiment_qc_file = "experiments_qc.csv"

aa_changes_file = "aa_changes.csv"
nt_changes_file = "nt_changes.csv"


@dataclass
class DirStructure:
    """
    A class to represent the directory structure of a summary.
    """

    summary_dir: Path

    metadata_dir: Path = field(init=False)
    seq_inventory_dir: Path = field(init=False)
    qc_dir: Path = field(init=False)
    variants_dir: Path = field(init=False)
    vcfs_dir: Path = field(init=False)
    gene_deletions_dir: Path = field(init=False)

    dirs: list[Path] = field(init=False)

    inventory_file: Path = field(init=False)
    metadata_file: Path = field(init=False)
    throughput_file: Path = field(init=False)

    coverage_file: Path = field(init=False)
    replicates_qc_file: Path = field(init=False)
    samples_qc_file: Path = field(init=False)
    samples_by_amplicon_qc_file: Path = field(init=False)
    experiment_qc_file: Path = field(init=False)
    aa_changes_file: Path = field(init=False)
    nt_changes_file: Path = field(init=False)

    def __post_init__(self):
        self.metadata_dir = self.summary_dir / metadata_dir
        self.seq_inventory_dir = self.summary_dir / seq_inventory_dir
        self.qc_dir = self.summary_dir / qc_dir
        self.variants_dir = self.summary_dir / variants_dir
        self.gene_deletions_dir = self.summary_dir / gene_deletions_dir
        self.vcfs_dir = self.summary_dir / variants_dir / vcfs_dir
        self.dirs = [
            self.metadata_dir,
            self.seq_inventory_dir,
            self.qc_dir,
            self.variants_dir,
            self.vcfs_dir,
            # no gene_deletions_dir here because it is optional
        ]

        self.metadata_file = self.metadata_dir / metadata_file

        self.inventory_file = self.seq_inventory_dir / sample_inventory_file
        self.throughput_file = self.seq_inventory_dir / throughput_file

        self.coverage_file = self.qc_dir / coverage_file
        self.replicates_qc_file = self.qc_dir / replicates_qc_file
        self.samples_qc_file = self.qc_dir / samples_qc_file
        self.samples_by_amplicon_qc_file = self.qc_dir / samples_by_amplicon_qc_file
        self.experiment_qc_file = self.qc_dir / experiment_qc_file

        self.aa_changes_file = self.variants_dir / aa_changes_file
        self.nt_changes_file = self.variants_dir / nt_changes_file


def looks_like_summary_dir(path: Path) -> bool:
    """
    Check if the given path looks like a summary directory.
    """
    required_dirs = [
        seq_inventory_dir,
        qc_dir,
    ]
    return all((path / d).exists() for d in required_dirs)
