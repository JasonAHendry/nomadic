import os
import json
import uuid
import subprocess
import pandas as pd

from typing import List
from abc import ABC, abstractmethod

from nomadic.map.mappers import Minimap2
from nomadic.download.references import REFERENCE_COLLECTION
from nomadic.util.dirs import produce_dir, ExperimentDirectories
from nomadic.util.regions import RegionBEDParser
from nomadic.util.samtools import (
    samtools_merge,
    samtools_index,
    samtools_flagstats,
    samtools_depth,
)
from nomadic.util.wrappers import bcftools
from nomadic.util.vcf import Consequence

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

        log.info("Running FASTQ count analysis...")
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

    def __init__(self, barcode_name: str, expt_dirs: ExperimentDirectories, ref_name: str="Pf3D7"):
        """
        Define initial directories and file paths
        TODO: choose mapping algorithm and reference genome

        """
        super().__init__(barcode_name, expt_dirs)

        self.step_dir = produce_dir(self.barcode_dir, self.step_name)
        self.inter_dir = produce_dir(self.step_dir, "intermediate")

        self.reference = REFERENCE_COLLECTION[ref_name]
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

        log.info("Running real-time mapping...")

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

        log.debug("Merging intermediate BAM files...")

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

    def __init__(self, barcode_name: str, expt_dirs: ExperimentDirectories, ref_name: str="Pf3D7"):
        """
        Define initial directories and file paths
        TODO: choose mapping algorithm and reference genome

        """

        super().__init__(barcode_name, expt_dirs)

        self.step_dir = produce_dir(self.barcode_dir, self.step_name)
        self.inter_dir = produce_dir(self.step_dir, "intermediate")
        self.ref_name = ref_name

        self.inter_jsons = []
        self.output_json = (
            f"{self.step_dir}/{self.barcode_name}.{self.ref_name}.flagstats.json"
        )

    def _get_intermediate_json_path(self):
        """
        Get path to a new intermediate JSON file

        """
    
        base_name = f"{self.inter_dir}/{self.barcode_name}.{self.ref_name}"
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

        log.debug("Merging `samtools flagstat` JSON files...")

        dts = [json.load(open(j, "r")) for j in self.inter_jsons]

        log.debug(f"Found {len(dts)} to merge.")

        keys = dts[0].keys()
        merged = {key: sum([dt[key] for dt in dts]) for key in keys}

        json.dump(merged, open(self.output_json, "w"))

        return self.output_json


# --------------------------------------------------------------------------------
# Analysis of coverage over BED regions
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
        self, barcode_name: str, expt_dirs: ExperimentDirectories, bed_path: str, ref_name: str="Pf3D7"):
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


# --------------------------------------------------------------------------------
# Analysis of depth profile across regions
#
# --------------------------------------------------------------------------------


class RegionDepthRT(AnalysisStepRT):
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
        ref_name: str="Pf3D7"):
        """Initialise output directory and define file names"""

        super().__init__(barcode_name, expt_dirs)

        self.regions = regions
        self.ref_name = ref_name
        self.output_dir = produce_dir(self.barcode_dir, self.step_name)
        self.region_output_dir = produce_dir(self.output_dir, "by_region")
        self.output_csv = f"{self.output_dir}/{self.barcode_name}.{self.ref_name}.depth.csv"

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
        Not needed here, we are not computing from intermediate BAMS

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


class CallVariantsRT(AnalysisStepRT):
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
        ref_name: str="Pf3D7"
    ):
        """Initialise output directory and define file names"""

        super().__init__(barcode_name, expt_dirs)

        self.bed_path = bed_path
        self.reference = REFERENCE_COLLECTION[ref_name]

        self.output_dir = produce_dir(self.barcode_dir, self.step_name)
        self.unfiltered_vcf = f"{self.output_dir}/{self.barcode_name}.{self.reference.name}.unfiltered.vcf.gz" 
        self.output_vcf = f"{self.output_dir}/{self.barcode_name}.{self.reference.name}.vcf.gz"

    def _get_reheader_command(self, input_vcf: str = "-", output_vcf: str = None) -> str:
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
        
        cmd = f"bcftools reheader {input_vcf} -s {sample_name_file}"
        if output_vcf is None:
            return cmd
        return f"{cmd} -o {output_vcf}"
    
    def _get_lowcomplexity_filter_command(self) -> str:
        """
        TODO: 
         - Make this optional
         - add input_vcf, output_vcf arguments

        """

        bed_mask_path = self.expt_dirs.regions_bed.replace(".bed", ".lowcomplexity_mask.bed")
        if not os.path.exists(bed_mask_path):
            print("Creating low-complexity mask for amplicons...")
            cmd = (
                "bedtools intersect"
                f" -a {self.reference.fasta_mask_path}"
                f" -b {self.expt_dirs.regions_bed}"
                f" -wa > {bed_mask_path}"
            )
            subprocess.run(cmd, shell=True, check=True)
            print("Done.")

        # Masking command
        cmd = (
            "bcftools filter"
            " --mode +"
            " --soft-filter LowComplexity"
            " --set-GTs ."
            f" --mask-file {bed_mask_path}"
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

        cmd_pileup = (
            "bcftools mpileup -B -I -Q12 --max-BQ 30 -h 100"
            f" --annotate {self.ANNOTATE_MPILEUP}"
            f" --max-depth {self.MAX_DEPTH}"
            f" -f {self.reference.fasta_path}"
            f" -Ov {input_bam}"
        )

        cmd_call = (
            "bcftools call -m -P 0.01"
            f" -a {self.ANNOTATE_CALL}"
            " -Oz - "
        )

        cmd_reheader = self._get_reheader_command(output_vcf=self.unfiltered_vcf)

        cmd = f"{cmd_pileup} | {cmd_call} | {cmd_reheader}"

        subprocess.run(cmd, check=True, shell=True)
        bcftools.index(self.unfiltered_vcf)

        cmd_view = (
            "bcftools view"
            " --max-alleles 2"
            f" -R {self.bed_path}"
            f" -Ou {self.unfiltered_vcf} "
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
            f" -Oz -o {self.output_vcf} - "
        )

        cmd = f"{cmd_view} | {cmd_depth_filter} | {cmd_lowcomplexity_filter} |{cmd_qual_filter}"

        subprocess.run(cmd, check=True, shell=True)
        bcftools.index(self.output_vcf)
        os.remove(self.unfiltered_vcf) # don't need unfiltered VCF
        os.remove(self.unfiltered_vcf + ".csi") # default index type

        return self.output_vcf

    def merge(self):
        pass


class AnnotateVariantsRT(AnalysisStepRT):
    """
    Annotate a set of variants in real-time 
    using `bcftools csq`

    TODO:
    - Handle phasing properly for variants
    with multiple consequences

    """

    step_name = "vcfs"
    AMP_HEADER = "##INFO=<ID=AMP_ID,Number=1,Type=String,Description=Amplicon identifier>"

    def __init__(
            self,
            barcode_name: str,
            expt_dirs: ExperimentDirectories,
            bed_path: str,
            ref_name: str="Pf3D7"
        ):
            """Initialise output directory and define file names"""

            super().__init__(barcode_name, expt_dirs)

            self.reference = REFERENCE_COLLECTION[ref_name]

            self.bed_path = bed_path

            self.output_dir = produce_dir(self.barcode_dir, self.step_name)
            self.output_vcf = f"{self.output_dir}/{self.barcode_name}.{self.reference.name}.annotated.vcf.gz"
            self.output_tsv = f"{self.output_dir}/{self.barcode_name}.{self.reference.name}.annotated.tsv"

    
    def _get_annotate_command(self, input_vcf: str="-", output_vcf: str = "") -> str:
        """
        Create a string representing command required to annotate variants with
        their amplicon position
        
        """
        cmd = "bcftools annotate"
        cmd += f" -a {self.bed_path}"
        cmd += " -c CHROM,FROM,TO,AMP_ID"
        cmd += f" -H '{self.AMP_HEADER}'"
        cmd += " -Oz"
        if output_vcf:
            cmd += f" -o {output_vcf}"
        cmd += f" {input_vcf}"

        return cmd
    
    def _get_csq_command(self, input_vcf: str="-", output_vcf: str = "") -> str:
        """
        Create a string representing command required
        to compute variant consequences
        
        """
        cmd = "bcftools csq"
        cmd += f" -f {self.reference.fasta_path}"
        cmd += f" -g {self.reference.gff_standard_path}"
        cmd += " --phase a"
        cmd += " -Oz"
        if output_vcf:
            cmd += f" -o {output_vcf}"
        cmd += f" {input_vcf}"

        return cmd
    

    def run(self, input_vcf: str):
        """
        Annotate variant calls using `bcftools csq`
        
        """

        cmd_annot = self._get_annotate_command(input_vcf)
        cmd_csq = self._get_csq_command("-", output_vcf=self.output_vcf)
        cmd = f"{cmd_annot} | {cmd_csq}"
        subprocess.run(cmd, shell=True, check=True)


    def _convert_to_tsv(self):
        """
        Convert the annotated VCF file to a small TSV file
        in preparation for plotting

        """
        
        fixed = {
            "chrom": "CHROM",
            "pos": "POS",
            "ref": "REF",
            "alt": "ALT",
            "qual": "QUAL",
            "consequence": "BCSQ",
            "amplicon": "AMP_ID"
        }
        called = {
            "gt": "GT",
            "dp": "DP",
            "wsaf": "VAF"
        }

        sep = "\\t"
        cmd_header = f" echo '{sep.join(list(fixed) + list(called))}\n' > {self.output_tsv}"
        sep += "%"
        cmd_query = "bcftools query"
        cmd_query += f" -f '%{sep.join(fixed.values())}\t[%{sep.join(called.values())}]\n'"
        cmd_query += f" {self.output_vcf} >> {self.output_tsv}"
        
        cmd = f"{cmd_header} && {cmd_query}"
        subprocess.run(cmd, shell=True, check=True)

    def _parse_consequences(self):
        """
        Parse the consequenc string in the TSV

        TODO:
        - What happens if there are NO mutations?
        - Then CSQ is an empty list

        """
        
        df = pd.read_csv(self.output_tsv, sep="\t")
        csqs = [
            Consequence.from_string(c)
            for c in df["consequence"]
        ]
        if csqs:
            mut_type, aa_change, strand = zip(*[(c.csq, c.get_concise_aa_change(), c.strand)
                                                for c in csqs])
        else:
            print(f"No mutations passed quality control for {self.barcode_name}.")
            mut_type, aa_change, strand = None, None, None

        df.insert(6, "mut_type", mut_type)
        df.insert(7, "aa_change", aa_change)
        df.insert(8, "strand", strand)
        df.drop("consequence", axis=1, inplace=True)
        df.to_csv(self.output_tsv, sep="\t", index=False)


    def merge(self):
        """
        This is more preparing to merge across barcodes

        """

        self._convert_to_tsv()
        self._parse_consequences()
