from dataclasses import dataclass, field
from pathlib import Path

panel_info_dir = "panel_info"
sample_info_dir = "sample_info"
qc_dir = "quality_control"
variants_dir = "variants"
vcfs_dir = "vcfs"
prevalence_dir = "prevalence"
gene_deletions_dir = "gene_deletions"

inventory_file = "inventory.csv"
metadata_file = "metadata.csv"
throughput_file = "throughput.csv"

coverage_file = "summary.coverage.csv"
replicates_qc_file = "summary.replicates_qc.csv"
samples_qc_file = "summary.samples_qc.csv"
samples_by_amplicon_qc_file = "summary.samples_amplicons_qc.csv"
experiment_qc_file = "summary.experiments_qc.csv"

aa_changes_file = "summary.aa_changes.csv"
nt_changes_file = "summary.nt_changes.csv"


@dataclass
class DirStructure:
    """
    A class to represent the directory structure of a summary.

    Attributes:
        summary_dir (Path): The path to the summary directory.
        sample_info_dir (Path): The path to the sample info directory.
        qc_dir (Path): The path to the quality control directory.
        variants_dir (Path): The path to the variants directory.
        prevalence_dir (Path): The path to the prevalence directory.
    """

    summary_dir: Path

    panel_info_dir: Path = field(init=False)
    sample_info_dir: Path = field(init=False)
    qc_dir: Path = field(init=False)
    variants_dir: Path = field(init=False)
    vcfs_dir: Path = field(init=False)
    prevalence_dir: Path = field(init=False)
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
        self.panel_info_dir = self.summary_dir / panel_info_dir
        self.sample_info_dir = self.summary_dir / sample_info_dir
        self.qc_dir = self.summary_dir / qc_dir
        self.variants_dir = self.summary_dir / variants_dir
        self.prevalence_dir = self.summary_dir / prevalence_dir
        self.gene_deletions_dir = self.summary_dir / gene_deletions_dir
        self.vcfs_dir = self.summary_dir / vcfs_dir
        self.dirs = [
            self.panel_info_dir,
            self.sample_info_dir,
            self.qc_dir,
            self.variants_dir,
            self.vcfs_dir,
            self.prevalence_dir,
            # no gene_deletions_dir here because it is optional
        ]

        self.inventory_file = self.sample_info_dir / inventory_file
        self.metadata_file = self.sample_info_dir / metadata_file
        self.throughput_file = self.sample_info_dir / throughput_file

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
        sample_info_dir,
        qc_dir,
        variants_dir,
    ]
    return all((path / d).exists() for d in required_dirs)
