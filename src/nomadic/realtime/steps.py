import os
import json
import uuid
import subprocess
import pandas as pd

from typing import List
from abc import ABC, abstractmethod

from nomadic.util.dirs import produce_dir, ExperimentDirectories
from nomadic.map.mappers import Minimap2
from nomadic.download.references import PlasmodiumFalciparum3D7
from nomadic.util.samtools import samtools_merge, samtools_index, samtools_flagstats


# --------------------------------------------------------------------------------
# Abstract base class
#
# --------------------------------------------------------------------------------


class AnalysisStepRT(ABC):
    """
    Run an analysis step in real-time

    """

    def __init__(self, barcode_name: str, expt_dirs: ExperimentDirectories):
        """
        Definne the barcode name and directories

        """

        self.barcode_name = barcode_name
        self.expt_dirs = expt_dirs
        self.barcode_dir = expt_dirs.get_barcode_dir(barcode_name)

    @abstractmethod
    def run(self, input_files: List[str]):
        """
        Run this analysis step

        """
        pass

    @abstractmethod
    def merge(self):
        """
        Merge the results from a real-time run with the
        already completed results

        """
        pass


# --------------------------------------------------------------------------------
# FASTQ counting step
#
# --------------------------------------------------------------------------------


class FASTQCountRT(AnalysisStepRT):
    """
    Map a set of FASTQ files from a given barcode,
    merge them with the cumulative BAM

    """

    step_name = ""

    def __init__(self, barcode_name: str, expt_dirs: ExperimentDirectories):
        """
        Define initial directories and file paths
        TODO: choose mapping algorithm and reference genome

        """
        super().__init__(barcode_name, expt_dirs)

        self.n_processed_fastq = 0
        self.output_json = (
            f"{self.barcode_dir}/{self.barcode_name}.n_processed_fastq.json"
        )

    def run(self, new_fastqs: List[str]):
        """
        Increase count of processed FASTQs, store in JSON

        """

        self.n_processed_fastq += len(new_fastqs)
        json.dump(
            {"barcode": self.barcode_name, "n_processed_fastq": self.n_processed_fastq},
            open(self.output_json, "w"),
        )

    def merge(self):
        """
        We do not need to merge because we are keeping a running count of the
        total number of FASTQ seen for this barcdoe
        """

        pass


# --------------------------------------------------------------------------------
# Mapping Step
#
# --------------------------------------------------------------------------------


class MappingRT(AnalysisStepRT):
    """
    Map a set of FASTQ files from a given barcode,
    merge them with the cumulative BAM

    """

    step_name = "bams"

    def __init__(self, barcode_name: str, expt_dirs: ExperimentDirectories):
        """
        Define initial directories and file paths
        TODO: choose mapping algorithm and reference genome

        """
        super().__init__(barcode_name, expt_dirs)

        self.step_dir = produce_dir(self.barcode_dir, self.step_name)
        self.inter_dir = produce_dir(self.step_dir, "intermediate")

        self.reference = PlasmodiumFalciparum3D7()
        self.mapper = Minimap2(self.reference)

        self.inter_bams = []
        self.output_bam = (
            f"{self.step_dir}/{self.barcode_name}.{self.reference.name}.final.bam"
        )

    def _get_intermediate_bam_path(self):
        """
        Get path to a new intermediate BAM file

        """

        base_name = f"{self.inter_dir}/{self.barcode_name}.{self.reference.name}"
        n = len(self.inter_bams)
        code = f"n{n:05d}"

        return f"{base_name}.{code}.bam"

    def run(self, new_fastqs: List[str]):
        """
        Map all of the newly observed FASTQ files

        """

        inter_bam = self._get_intermediate_bam_path()
        self.mapper.map_from_fastqs(fastq_path=" ".join(new_fastqs))
        self.mapper.run(output_bam=inter_bam)
        samtools_index(inter_bam)
        self.inter_bams.append(inter_bam)

        return inter_bam

    def merge(self):
        """
        Merge all of the intermediate BAM files, and index
        the final bam

        """

        samtools_merge(bam_files=self.inter_bams, output_bam=self.output_bam)
        samtools_index(self.output_bam)

        return self.output_bam


# --------------------------------------------------------------------------------
# Analysis of FLAG statistics after mapping
#
# --------------------------------------------------------------------------------


class FlagstatsRT(AnalysisStepRT):
    """
    Analyse mapping flag statistics in real-time

    """

    step_name = "qcbams"

    def __init__(self, barcode_name: str, expt_dirs: ExperimentDirectories):
        """
        Define initial directories and file paths
        TODO: choose mapping algorithm and reference genome

        """

        super().__init__(barcode_name, expt_dirs)

        self.step_dir = produce_dir(self.barcode_dir, self.step_name)
        self.inter_dir = produce_dir(self.step_dir, "intermediate")

        self.reference = PlasmodiumFalciparum3D7()

        self.inter_jsons = []
        self.output_json = (
            f"{self.step_dir}/{self.barcode_name}.{self.reference.name}.flagstats.json"
        )

    def _get_intermediate_json_path(self):
        """
        Get path to a new intermediate JSON file

        """

        base_name = f"{self.inter_dir}/{self.barcode_name}.{self.reference.name}"
        n = len(self.inter_jsons)
        code = f"n{n:05d}"

        return f"{base_name}.{code}.json"

    def run(self, inter_bam: str):
        """
        Run `samtools flagstats` for an intermediate BAM file produced
        as part of the real-time analysis

        """

        inter_json = self._get_intermediate_json_path()
        samtools_flagstats(input_bam=inter_bam, output_json=inter_json)
        self.inter_jsons.append(inter_json)

        return inter_json

    def merge(self):
        """
        Merge all the JSON files generated for intermediate BAM files

        """

        dts = [json.load(open(j, "r")) for j in self.inter_jsons]

        keys = dts[0].keys()
        merged = {key: sum([dt[key] for dt in dts]) for key in keys}

        json.dump(merged, open(self.output_json, "w"))

        return self.output_json


# --------------------------------------------------------------------------------
# Analysis of coveragee over BED regions
#
# --------------------------------------------------------------------------------


class BedCovRT(AnalysisStepRT):
    """
    Analyse coverage across a set of regions defined
    by a BED file

    """

    step_name = "bedcov"  # currently unused
    cov_thresh = 100

    def __init__(
        self, barcode_name: str, expt_dirs: ExperimentDirectories, bed_path: str
    ):
        """
        Define initial directories and file paths

        """

        super().__init__(barcode_name, expt_dirs)
        self.reference = (
            PlasmodiumFalciparum3D7()
        )  # TODO: always required for naming, it seems

        self.bed_path = bed_path
        self.output_csv = (
            f"{self.barcode_dir}/{self.barcode_name}.{self.reference.name}.bedcov.csv"
        )

        self.cov_col_name = f"cov_gr{self.cov_thresh}"

    def run(self, input_bam):
        """
        Compute coverage and summary statistics from an `input_bam`

        """

        # Run, writing to temporary BED
        temp_bed = f"{input_bam[:-4]}.temp.{str(uuid.uuid4())[:8]}.bed"
        cmd = f"samtools bedcov -d {self.cov_thresh} -c"
        cmd += f" {self.bed_path} {input_bam} > {temp_bed}"
        subprocess.run(cmd, shell=True, check=True)

        # Load temp BED, compute summaries, write
        df = pd.read_csv(temp_bed, sep="\t", header=None)
        df.columns = [
            "chrom",
            "start",
            "end",
            "name",
            "total_cov",
            self.cov_col_name,
            "n_reads",
        ]
        df.insert(0, "barcode", self.barcode_name)
        df.insert(3, "length", df["end"] - df["start"])
        df.insert(5, "mean_cov", df["total_cov"] / df["length"])
        df.insert(
            7, f"{self.cov_col_name}_per", 100 * df[self.cov_col_name] / df["length"]
        )
        df.to_csv(self.output_csv, index=False)

        # Remove temp BED
        os.remove(temp_bed)

    def merge(self):
        """
        Not needed here, we are not computing from intermediate BAMS

        """
        pass
