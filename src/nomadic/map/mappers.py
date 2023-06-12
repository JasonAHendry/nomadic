import os
import subprocess
from abc import ABC, abstractmethod


# ================================================================
# Define abstract base class for different mapping
# algorithms
# ================================================================


class MappingAlgorithm(ABC):
    def __init__(self, reference):
        """
        Abstract base class for implementing different mapping algorithms

        Input can be either .fastq files or a .bam file.

        """
        # Set reference
        self.reference = reference

        # Set defaults
        self.map_cmd = None
        self.remap_cmd = ""
        self.input_fastqs = "-"

    def remap_from_bam(self, input_bam):
        """
        Prepare to remap unmapped reads in a .bam file at `input_bam`
        to the `reference`

        """
        self.remap_cmd = f"samtools view -f 0x004 {input_bam}"
        self.remap_cmd += " | samtools fastq | "

    def map_from_fastqs(self, fastq_dir=None, fastq_path=None):
        """
        Prepare to map all .fastq files found in a directory `fastq_dir`

        """
        if fastq_dir is not None:
            fastq_dir = fastq_dir
            fastqs = [
                f"{fastq_dir}/{fastq}"
                for fastq in os.listdir(fastq_dir)
                if fastq.endswith(".fastq") or fastq.endswith(".fastq.gz")
            ]
            self.input_fastqs = " ".join(fastqs)
        elif fastq_path is not None:
            self.input_fastqs = fastq_path
        else:
            raise ValueError("Must set either `fastq_dir` or `fastq_path`.")

    @abstractmethod
    def _define_mapping_command(self, output_bam, flags):
        """
        Define the command for the mapping algorithm
        """
        pass

    def run(self, output_bam, verbose=False):
        """
        Run the mapping algorithm inputs

        """
        # Define the mapping command
        self._define_mapping_command(output_bam)
        if self.map_cmd is None:
            raise ValueError("Must define mapping command before running algorithm.")

        # Combine optional remapping command with mapping command
        cmd = self.remap_cmd + self.map_cmd

        if verbose:
            print(f"Complete command: {cmd}")

        # Run
        subprocess.run(cmd, shell=True, check=True)


# ================================================================
# Concrete implementations of mapping algorithms
#
# ================================================================


class Minimap2(MappingAlgorithm):
    """
    Map long reads with `minimap2`

    """

    def _define_mapping_command(self, output_bam, flags="--eqx --MD"):
        """
        Run minimap2, compress result to .bam file, and sort

        """
        self.map_cmd = "minimap2"
        self.map_cmd += (
            f" -ax map-ont {flags} {self.reference.fasta_path} {self.input_fastqs} |"
        )
        self.map_cmd += " samtools view -S -b - |"
        self.map_cmd += f" samtools sort -o {output_bam}"


class BwaMem(MappingAlgorithm):
    """
    Map short reads with `bwa mem`

    """

    def create_reference_index(self):
        """
        Create an index for bwa

        Produces a set of index files with the same prefix
        as `self.reference.fasta_path`.

        """
        index_cmd = f"bwa index {self.reference.fasta_path}"
        subprocess.run(index_cmd, shell=True, check=True)

    def _define_mapping_command(self, output_bam, flags=""):
        """
        Run bwa, compress result to .bam file, and sort

        """
        self.map_cmd = "bwa mem"
        self.map_cmd += " -R '@RG\\tID:misc\\tSM:pool'"  # ID and SM tags needed for gatk HaplotypeCaller
        self.map_cmd += f" {flags} {self.reference.fasta_path} {self.input_fastqs} |"
        self.map_cmd += " samtools view -S -b - |"
        self.map_cmd += f" samtools sort -o {output_bam}"


# ================================================================
# Create a collection of mapping algorithms
#
# ================================================================


MAPPER_COLLECTION = {"minimap2": Minimap2, "bwa": BwaMem}
