import os
import json
import uuid
import glob
import subprocess
import pandas as pd
import shlex

from typing import List
from abc import ABC, abstractmethod

from nomadic.map.mappers import Minimap2
from nomadic.download.references import REFERENCE_COLLECTION
from nomadic.util.dirs import produce_dir
from nomadic.util.experiment import ExperimentDirectories
from nomadic.util.regions import RegionBEDParser
from nomadic.util.samtools import (
    samtools_merge,
    samtools_index,
    samtools_flagstats,
    samtools_depth,
)
from nomadic.util.wrappers import bcftools

import logging

log = logging.getLogger("nomadic")


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


class FASTQProcessedRT(AnalysisStepRT):
    """
    Store the number of processed FASTQ files
    """

    step_name = "fastq"

    def __init__(self, barcode_name: str, expt_dirs: ExperimentDirectories):
        super().__init__(barcode_name, expt_dirs)

        self.step_dir = produce_dir(self.barcode_dir, self.step_name)
        self.incr_dir = produce_dir(self.step_dir, ".incremental")
        self.output_json = f"{self.step_dir}/{self.barcode_name}.n_processed_fastq.json"

    def _get_incremental_json_path(self, incr_id: str):
        """
        Get path to a new incremental json file

        """
        base_name = f"{self.incr_dir}/{self.barcode_name}"
        return f"{base_name}.{incr_id}.n_processed_fastq.json"

    def _get_incremental_jsons(self):
        base_name = f"{self.incr_dir}/{self.barcode_name}"
        return [file for file in glob.glob(f"{base_name}*.n_processed_fastq.json")]

    def run(self, new_fastqs: List[str], incr_id: str):
        """
        Increase count of processed FASTQs, store in JSON

        """
        log.info("Running FASTQ count analysis...")
        json.dump(
            {"barcode": self.barcode_name, "n_processed_fastq": len(new_fastqs)},
            open(self._get_incremental_json_path(incr_id), "w"),
        )

    def merge(self):
        log.debug("Merging FASTQ count files...")
        fastq_data = [json.load(open(j, "r")) for j in self._get_incremental_jsons()]
        n_processed_fastq = sum([data["n_processed_fastq"] for data in fastq_data])
        json.dump(
            {"barcode": self.barcode_name, "n_processed_fastq": n_processed_fastq},
            open(self.output_json, "w"),
        )


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

    def __init__(
        self,
        barcode_name: str,
        expt_dirs: ExperimentDirectories,
        ref_name: str = "Pf3D7",
    ):
        """
        Define initial directories and file paths
        TODO: choose mapping algorithm and reference genome

        """
        super().__init__(barcode_name, expt_dirs)

        self.step_dir = produce_dir(self.barcode_dir, self.step_name)
        self.incr_dir = produce_dir(self.step_dir, ".incremental")

        self.reference = REFERENCE_COLLECTION[ref_name]
        self.mapper = Minimap2(self.reference)

        self.output_bam = (
            f"{self.step_dir}/{self.barcode_name}.{self.reference.name}.final.bam"
        )

    def _get_incremental_bam_path(self, incr_id: str):
        """
        Get path to a new incremental BAM file

        """

        base_name = f"{self.incr_dir}/{self.barcode_name}.{self.reference.name}"
        return f"{base_name}.{incr_id}.bam"

    def _get_incremental_bams(self):
        base_name = f"{self.incr_dir}/{self.barcode_name}.{self.reference.name}"
        return [file for file in glob.glob(f"{base_name}*.bam")]

    def run(self, new_fastqs: List[str], incr_id: str):
        """
        Map all of the newly observed FASTQ files

        """

        log.info("Running real-time mapping...")

        incr_bam = self._get_incremental_bam_path(incr_id)
        self.mapper.map_from_fastqs(fastq_paths=new_fastqs)
        self.mapper.run(output_bam=incr_bam)
        samtools_index(incr_bam)

        return incr_bam

    def merge(self):
        """
        Merge all of the incremental BAM files, and index
        the final bam

        """

        log.debug("Merging incremental BAM files...")

        samtools_merge(
            bam_files=self._get_incremental_bams(), output_bam=self.output_bam
        )
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

    def __init__(
        self,
        barcode_name: str,
        expt_dirs: ExperimentDirectories,
        ref_name: str = "Pf3D7",
    ):
        super().__init__(barcode_name, expt_dirs)

        self.step_dir = produce_dir(self.barcode_dir, self.step_name)
        self.incr_dir = produce_dir(self.step_dir, ".incremental")
        self.ref_name = ref_name

        self.output_json = (
            f"{self.step_dir}/{self.barcode_name}.{self.ref_name}.flagstats.json"
        )

    def _get_incremental_json_path(self, incr_id: str) -> str:
        """
        Get path to a new incremental JSON file

        """

        base_name = f"{self.incr_dir}/{self.barcode_name}.{self.ref_name}"
        return f"{base_name}.{incr_id}.json"

    def _get_incremental_jsons(self) -> list[str]:
        """
        Get all incremental JSON files

        """
        base_name = f"{self.incr_dir}/{self.barcode_name}.{self.ref_name}"
        return [file for file in glob.glob(f"{base_name}*.json")]

    def run(self, incr_bam: str, incr_id: str):
        """
        Run `samtools flagstats` for an incremental BAM file produced
        as part of the real-time analysis

        """

        incr_json = self._get_incremental_json_path(incr_id)
        samtools_flagstats(input_bam=incr_bam, output_json=incr_json)

        return incr_json

    def merge(self):
        """
        Merge all the JSON files generated for incremental BAM files

        """

        log.debug("Merging `samtools flagstat` JSON files...")

        dts = [json.load(open(j, "r")) for j in self._get_incremental_jsons()]

        log.debug(f"Found {len(dts)} to merge.")

        keys = dts[0].keys()
        merged = {key: sum([dt[key] for dt in dts]) for key in keys}

        json.dump(merged, open(self.output_json, "w"))

        return self.output_json


# --------------------------------------------------------------------------------
# Analysis of coverage over BED regions
#
# --------------------------------------------------------------------------------


class RegionCoverage(AnalysisStepRT):
    """
    Analyse coverage across a set of regions defined
    by a BED file

    """

    step_name = "bedcov"  # currently unused
    cov_thresh = 100

    def __init__(
        self,
        barcode_name: str,
        expt_dirs: ExperimentDirectories,
        bed_path: str,
        ref_name: str = "Pf3D7",
    ):
        """
        Define initial directories and file paths

        """

        super().__init__(barcode_name, expt_dirs)
        self.ref_name = ref_name

        self.step_dir = produce_dir(self.barcode_dir, self.step_name)

        self.bed_path = bed_path
        self.output_csv = (
            f"{self.step_dir}/{self.barcode_name}.{self.ref_name}.bedcov.csv"
        )

        self.cov_col_name = f"cov_gr{self.cov_thresh}"

    def run(self, input_bam):
        """
        Compute coverage and summary statistics from an `input_bam`

        """

        # Run, writing to temporary BED
        temp_bed = f"{input_bam[:-4]}.temp.{str(uuid.uuid4())[:8]}.bed"
        cmd = [
            "samtools",
            "bedcov",
            "-d",
            str(self.cov_thresh),
            "-c",
            self.bed_path,
            input_bam,
        ]
        with open(temp_bed, "w") as file:
            subprocess.run(cmd, check=True, stdout=file)

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
        Not needed here, we are not computing from incremental BAMS

        """
        pass


# --------------------------------------------------------------------------------
# Analysis of depth profile across regions
#
# --------------------------------------------------------------------------------


class RegionDepthProfileRT(AnalysisStepRT):
    """
    Analyse coverage across a set of regions defined
    by a BED file

    """

    step_name = "depth"

    def __init__(
        self,
        barcode_name: str,
        expt_dirs: ExperimentDirectories,
        regions: RegionBEDParser,
        ref_name: str,
    ):
        """Initialise output directory and define file names"""

        super().__init__(barcode_name, expt_dirs)

        self.regions = regions
        self.ref_name = ref_name
        self.output_dir = produce_dir(self.barcode_dir, self.step_name)
        self.region_output_dir = produce_dir(self.output_dir, "by_region")
        self.output_csv = (
            f"{self.output_dir}/{self.barcode_name}.{self.ref_name}.depth.csv"
        )

    def run(self, input_bam):
        """
        Compute coverage and summary statistics from an `input_bam`

        """

        for region_name, region_str in self.regions.str_format.items():
            samtools_depth(
                input_bam=input_bam,
                output_path=f"{self.region_output_dir}/{self.barcode_name}.{region_name}.depth",
                region_str=region_str,
            )

    def merge(self):
        """
        Not needed here, we are not computing from incremental BAMS

        TODO: Edge case where there is no coverage for a target?

        """

        dfs = []
        for region_name in self.regions.names:
            df = pd.read_csv(
                f"{self.region_output_dir}/{self.barcode_name}.{region_name}.depth",
                sep="\t",
                header=None,
                names=["chrom", "pos", "depth"],
            )
            df.insert(0, "barcode", self.barcode_name)
            df.insert(0, "name", region_name)
            dfs.append(df)
        merged_df = pd.concat(dfs)
        merged_df.to_csv(self.output_csv, index=False)


# --------------------------------------------------------------------------------
# Variant calling and annnotation
#
# --------------------------------------------------------------------------------
class CallVariantsRTDelve(AnalysisStepRT):
    step_name = "vcfs"

    # Settings
    MAX_DEPTH = 5000
    MIN_DEPTH = 40
    MIN_QUAL = 20
    STRAND_BIAS_ODDS_RATIO = 7
    MODEL_PARAMS = [8, 8, 0.01]
    COMPUTE_BAQ = True

    def __init__(
        self,
        barcode_name: str,
        expt_dirs: ExperimentDirectories,
        bed_path: str,
        ref_name: str = "Pf3D7",
    ):
        """Initialise output directory and define file names"""

        super().__init__(barcode_name, expt_dirs)

        self.bed_path = bed_path
        self.reference = REFERENCE_COLLECTION[ref_name]

        self.output_dir = produce_dir(self.barcode_dir, self.step_name)
        self.output_vcf = (
            f"{self.output_dir}/{self.barcode_name}.{self.reference.name}.vcf.gz"
        )

    def _get_lowcomplexity_filter_command(self) -> str:
        bed_mask_path = self.expt_dirs.regions_bed.replace(
            ".bed", ".lowcomplexity_mask.bed"
        )
        if not os.path.exists(bed_mask_path):
            print("Creating low-complexity mask for amplicons...")
            cmd = [
                "bedtools",
                "intersect",
                "-a",
                self.reference.fasta_mask_path,
                "-b",
                self.expt_dirs.regions_bed,
                "-wa",
            ]
            with open(bed_mask_path, "w") as file:
                subprocess.run(cmd, check=True, stdout=file)
            print("Done.")

        # Masking command
        cmd = (
            "bcftools filter"
            " --mode +"
            " --soft-filter LowComplexity"
            " --set-GTs ."
            f" --mask-file {shlex.quote(bed_mask_path)}"
        )

        return cmd

    def run(self, input_bam: str) -> str:
        """
        Run variant calling with delve

        """
        log.info("Calling variants with delve...")

        filtered_bam_path = input_bam.replace(".bam", ".filtered.bam")

        cmd_filter_input = f"samtools view {shlex.quote(input_bam)} -e '![SA]' -b -o {shlex.quote(filtered_bam_path)}"

        compute_baq = "--compute-baq" if self.COMPUTE_BAQ else ""

        cmd_call = (
            "delve call"
            f" --model-params {','.join([str(v) for v in self.MODEL_PARAMS])}"
            f" --min-cov {self.MIN_DEPTH}"
            f" --min-BQ {self.MIN_QUAL}"
            f" --max-cov {self.MAX_DEPTH}"
            f" {compute_baq}"
            f" -s {shlex.quote(self.barcode_name)}"
            f" -R {shlex.quote(self.bed_path)}"
            f" -f {shlex.quote(self.reference.fasta_path)}"
            " --set-failed-GTs ."
            f" {shlex.quote(filtered_bam_path)}"
        )

        cmd_lowcomplexity_filter = self._get_lowcomplexity_filter_command()
        cmd_sort_index = f"bcftools sort - -Oz -W -o {shlex.quote(self.output_vcf)}"

        cmd = f"{cmd_call} | {cmd_lowcomplexity_filter} | {cmd_sort_index}"

        subprocess.run(cmd_filter_input, check=True, shell=True)
        samtools_index(filtered_bam_path)
        subprocess.run(cmd, check=True, shell=True)

        os.remove(filtered_bam_path)

        return self.output_vcf

    def merge(self):
        pass


class CallVariantsRTBcftools(AnalysisStepRT):
    """
    Call variants in real-time using `bcftools call`

    For now, this class also includes filtering.

    """

    step_name = "vcfs"

    # Settings
    ANNOTATE_MPILEUP = "FORMAT/DP,FORMAT/AD"
    ANNOTATE_CALL = "FORMAT/GQ"
    MAX_DEPTH = 5000
    MIN_DEPTH = 40
    MIN_QUAL = 15

    # Annotation
    AMP_HEADER = (
        "##INFO=<ID=AMP_ID,Number=1,Type=String,Description=Amplicon identifier>"
    )

    def __init__(
        self,
        barcode_name: str,
        expt_dirs: ExperimentDirectories,
        bed_path: str,
        ref_name: str = "Pf3D7",
    ):
        """Initialise output directory and define file names"""

        super().__init__(barcode_name, expt_dirs)

        self.bed_path = bed_path
        self.reference = REFERENCE_COLLECTION[ref_name]

        self.output_dir = produce_dir(self.barcode_dir, self.step_name)
        self.unfiltered_vcf = f"{self.output_dir}/{self.barcode_name}.{self.reference.name}.unfiltered.vcf.gz"
        self.output_vcf = (
            f"{self.output_dir}/{self.barcode_name}.{self.reference.name}.vcf.gz"
        )

    def _get_reheader_command(
        self, input_vcf: str = "-", output_vcf: str = None
    ) -> str:
        """
        Reassign the sample name inside of the VCF to the
        `self.barcode_name`

        Unfortunately, need to load the name from a file
        """

        # Write the file if it doesn't exist
        sample_name_file = f"{self.output_dir}/.sample_name.txt"
        if not os.path.exists(sample_name_file):
            with open(sample_name_file, "w") as file:
                file.write(f"{self.barcode_name}\n")

        cmd = f"bcftools reheader {shlex.quote(input_vcf)} -s {shlex.quote(sample_name_file)}"
        if output_vcf is None:
            return cmd
        return f"{cmd} -o {shlex.quote(output_vcf)}"

    def _get_lowcomplexity_filter_command(self) -> str:
        """
        TODO:
         - Make this optional
         - add input_vcf, output_vcf arguments

        """

        bed_mask_path = self.expt_dirs.regions_bed.replace(
            ".bed", ".lowcomplexity_mask.bed"
        )
        if not os.path.exists(bed_mask_path):
            print("Creating low-complexity mask for amplicons...")
            cmd = [
                "bedtools",
                "intersect",
                "-a",
                self.reference.fasta_mask_path,
                "-b",
                self.expt_dirs.regions_bed,
                "-wa",
            ]
            with open(bed_mask_path, "w") as file:
                subprocess.run(cmd, check=True, stdout=file)
            print("Done.")

        # Masking command
        cmd = (
            "bcftools filter"
            " --mode +"
            " --soft-filter LowComplexity"
            " --set-GTs ."
            f" --mask-file {shlex.quote(bed_mask_path)}"
            " -Ou -"
        )

        return cmd

    def run(self, input_bam: str) -> str:
        """
        Run variant calling with bcftools

        TODO:
        - I probably want to WRITE also to BCF

        NOTES:
        Here I am following recommendations of
        `bcftools` in calling with ONT reads, that is, using
        `-X ont` for `bcftools mpileup`, which sets:

        `-B`  : disable per-base alignment quality
        `-Q5` : Skip bases with base quality < 5  # I probably want this higher
        `--max-BQ 30` : Set the maximum base quality to 30
            - ONT sets homopolymers to 90 for some reason
        `-I`  : Skip indel calling
            - I only report SNPs in dashboard
            - But could be difficult for some variants (K76T)

        Then for `bcftools call`, I use `-P 0.01`,

        `-P` : Prior on mutation rate

        """
        log.info("Calling variants with bcftools...")

        cmd_pileup = (
            "bcftools mpileup -B -I -Q12 --max-BQ 30 -h 100"
            f" --annotate {self.ANNOTATE_MPILEUP}"
            f" --max-depth {self.MAX_DEPTH}"
            f" -f {shlex.quote(self.reference.fasta_path)}"
            f" -Ov {shlex.quote(input_bam)}"
        )

        cmd_call = f"bcftools call -m -P 0.01 -a {self.ANNOTATE_CALL} -Oz - "

        cmd_reheader = self._get_reheader_command(output_vcf=self.unfiltered_vcf)

        cmd = f"{cmd_pileup} | {cmd_call} | {cmd_reheader}"

        subprocess.run(cmd, check=True, shell=True)
        bcftools.index(self.unfiltered_vcf)

        cmd_view = (
            "bcftools view"
            " --max-alleles 2"
            f" -R {shlex.quote(self.bed_path)}"
            f" -Ou {shlex.quote(self.unfiltered_vcf)} "
        )

        cmd_depth_filter = (
            "bcftools filter"
            " --mode +"
            " --soft-filter LowDepth"
            " --set-GTs ."
            f" --exclude 'FORMAT/DP < {self.MIN_DEPTH}'"
            " -Ou - "
        )

        cmd_lowcomplexity_filter = self._get_lowcomplexity_filter_command()

        cmd_qual_filter = (
            "bcftools filter"
            " --mode +"
            " --soft-filter LowQual"
            " --set-GTs ."
            f" --exclude 'QUAL < {self.MIN_QUAL}'"
            f" -Oz -o {shlex.quote(self.output_vcf)} - "
        )

        cmd = f"{cmd_view} | {cmd_depth_filter} | {cmd_lowcomplexity_filter} |{cmd_qual_filter}"

        subprocess.run(cmd, check=True, shell=True)
        bcftools.index(self.output_vcf)
        os.remove(self.unfiltered_vcf)  # don't need unfiltered VCF
        os.remove(self.unfiltered_vcf + ".csi")  # default index type

        return self.output_vcf

    def merge(self):
        pass
