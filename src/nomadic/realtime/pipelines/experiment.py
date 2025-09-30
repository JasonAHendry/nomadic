import os
import json
import pandas as pd
import warnings
import subprocess
from abc import ABC, abstractmethod
import shlex

from nomadic.download.references import Reference, PlasmodiumFalciparum3D7
from nomadic.util.experiment import ExperimentDirectories
from nomadic.util.metadata import MetadataTableParser
from nomadic.util.dirs import produce_dir
from nomadic.util.regions import RegionBEDParser
from nomadic.util.wrappers import bcftools
from nomadic.util.vcf import VariantAnnotator


# --------------------------------------------------------------------------------
# Interface for experiment pipelines
#
# TODO:
# - Still some issues with the file names being defined inside methods
# --------------------------------------------------------------------------------


class ExperimentPipelineRT(ABC):
    """
    Interface for complete experiment pipeline, summarising
    information across multiple barcodes

    We implement some concrete _run_<step>() methods that will be the
    same across all subclasses, to reduce code duplication.

    """

    def __init__(
        self,
        metadata: MetadataTableParser,
        expt_dirs: ExperimentDirectories,
        regions: RegionBEDParser,
        reference: Reference = PlasmodiumFalciparum3D7(),
    ):
        """
        Store important metadata as instance attributes

        """

        self.expt_dirs = expt_dirs
        self.metadata = metadata
        self.regions = regions
        self.reference = reference

    @abstractmethod
    def run(self):
        """
        Run the complete pipeline

        """
        pass

    def _run_fastq(self):
        """
        Summarise FASTQ step results across all barcodes

        """
        fastq_dts = []
        for b in self.metadata.barcodes:
            barcode_dir = self.expt_dirs.get_barcode_dir(b)
            try:
                dt = json.load(
                    open(f"{barcode_dir}/fastq/{b}.n_processed_fastq.json", "r")
                )
                fastq_dts.append(dt)
            except FileNotFoundError:
                continue
        df = pd.DataFrame(fastq_dts)
        df = df.join(self.metadata.required_metadata, on="barcode")
        df_path = self.expt_dirs.get_summary_files().fastqs_processed
        df.to_csv(df_path, index=False)

    def _run_qcbams(self):
        """
        Summarise BAM quality control results across all barcodes

        """
        qcbams_dts = []
        for b in self.metadata.barcodes:
            barcode_dir = self.expt_dirs.get_barcode_dir(b)
            try:
                dt = json.load(
                    open(
                        f"{barcode_dir}/qcbams/{b}.{self.reference.name}.flagstats.json"
                    )
                )
                dt["barcode"] = b
                qcbams_dts.append(dt)
            except FileNotFoundError:
                continue
        df = pd.DataFrame(qcbams_dts)
        df = df.join(self.metadata.required_metadata, on="barcode")
        df_path = self.expt_dirs.get_summary_files().read_mapping
        df.to_csv(df_path, index=False)

    def _run_bedcov(self):
        """
        Summarise BED coverage results across all barcodes

        """
        bedcov_dfs = []
        for b in self.metadata.barcodes:
            barcode_dir = self.expt_dirs.get_barcode_dir(b)
            try:
                df = pd.read_csv(
                    f"{barcode_dir}/bedcov/{b}.{self.reference.name}.bedcov.csv"
                )
                bedcov_dfs.append(df)
            except FileNotFoundError:
                continue
        bedcov_df = pd.concat(bedcov_dfs)
        bedcov_df = bedcov_df.join(self.metadata.required_metadata, on="barcode")
        df_path = self.expt_dirs.get_summary_files().region_coverage
        bedcov_df.to_csv(df_path, index=False)

    def _run_depth(self):
        """
        Summarise depth results across all barcodes

        """
        depth_dfs = []
        for b in self.metadata.barcodes:
            barcode_dir = self.expt_dirs.get_barcode_dir(b)
            try:
                df = pd.read_csv(
                    f"{barcode_dir}/depth/{b}.{self.reference.name}.depth.csv"
                )
                depth_dfs.append(df)
            except FileNotFoundError:
                continue
        depth_df = pd.concat(depth_dfs)
        depth_df = depth_df.join(self.metadata.required_metadata, on="barcode")
        df_path = self.expt_dirs.get_summary_files().depth_profiles
        depth_df.to_csv(df_path, index=False)

    def _run_variant(self, caller: str):
        """
        Summarise variant calling results across all barcodes

        """

        # Merge VCF
        vcf_dir = produce_dir(self.expt_dirs.approach_dir, "vcfs")
        vcfs = []
        for b in self.metadata.barcodes:
            if b == "unclassified":
                continue
            barcode_dir = self.expt_dirs.get_barcode_dir(b)
            vcf = f"{barcode_dir}/vcfs/{b}.{self.reference.name}.vcf.gz"
            if not os.path.exists(vcf):
                warnings.warn(f"No VCF file found at: {vcf}")
                continue
            vcfs.append(vcf)

        if len(vcfs) <= 1:
            ### bcftools will complain if we try to merge less than 2 VCFs
            return

        unfiltered_vcf = f"{vcf_dir}/summary.variants.vcf.gz"
        bcftools.merge(vcfs, output_vcf=unfiltered_vcf)

        # Filter variant sites
        filtered_vcf = unfiltered_vcf.replace(".vcf.gz", ".filtered.vcf.gz")
        cmd = (
            "bcftools view"
            " --apply-filters PASS"
            " --types='snps'"
            " --min-alleles 2"
            f" -Oz -o {shlex.quote(filtered_vcf)} {shlex.quote(unfiltered_vcf)}"
        )
        subprocess.run(cmd, check=True, shell=True)

        # Annotate
        annotator = VariantAnnotator(
            input_vcf=filtered_vcf,
            bed_path=self.regions.path,
            reference=self.reference,
            caller=caller,
            output_vcf=filtered_vcf.replace(".vcf.gz", ".annotated.vcf.gz"),
        )
        annotator.run()
        csv_path = self.expt_dirs.get_summary_files().variants
        temp_path = csv_path.replace(".csv", "temp.csv")
        annotator.convert_to_csv(temp_path)

        df = pd.read_csv(temp_path)
        df = df.join(self.metadata.required_metadata, on="barcode")
        df.to_csv(csv_path)

        # Clean-up
        os.remove(temp_path)
        os.remove(filtered_vcf)


# --------------------------------------------------------------------------------
# Concrete experiment pipelines
#
# --------------------------------------------------------------------------------


class BasicExptPipeline(ExperimentPipelineRT):
    """
    Basic experiment pipeline, simply summarises the number
    of FASTQ files processed for each barcode

    """

    def run(self):
        self._run_fastq()


class ExptMappingPipelineRT(ExperimentPipelineRT):
    """
    Experiment pipeline summarise mapping results, coverage
    over amplicons and quality control

    """

    def run(self):
        self._run_fastq()
        self._run_qcbams()
        self._run_bedcov()
        self._run_depth()


class ExptCallingPipelineRT(ExperimentPipelineRT):
    """
    Experiment pipeline summarise mapping results, coverage
    over amplicons and quality control as well as variant calling

    """

    def __init__(
        self,
        metadata: MetadataTableParser,
        expt_dirs: ExperimentDirectories,
        regions: RegionBEDParser,
        caller: str,
        reference: Reference = PlasmodiumFalciparum3D7(),
    ):
        self.caller = caller
        super().__init__(metadata, expt_dirs, regions, reference)

    def run(self):
        self._run_fastq()
        self._run_qcbams()
        self._run_bedcov()
        self._run_depth()
        self._run_variant(self.caller)
