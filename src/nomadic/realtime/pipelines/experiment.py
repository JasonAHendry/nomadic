import json
import pandas as pd

from abc import ABC, abstractmethod

from nomadic.util.metadata import MetadataTableParser
from nomadic.util.dirs import ExperimentDirectories



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
    
    TODO: needs to take reference name

    """

    def __init__(self, metadata: MetadataTableParser, expt_dirs: ExperimentDirectories, ref_name: str="Pf3D7"):
        """
        Store important metadata as instance attributes

        """

        self.expt_dirs = expt_dirs
        self.metadata = metadata
        self.ref_name = ref_name

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
                dt = json.load(open(f"{barcode_dir}/{b}.n_processed_fastq.json", "r"))
                fastq_dts.append(dt)
            except FileNotFoundError:
                continue
        df = pd.DataFrame(fastq_dts)
        df_path = f"{self.expt_dirs.approach_dir}/summary.fastq.csv"
        df.to_csv(df_path, index=False)

    def _run_qcbams(self):
        """
        Summarise BAM quality control results across all barcodes

        """
        qcbams_dts = []
        for b in self.metadata.barcodes:
            barcode_dir = self.expt_dirs.get_barcode_dir(b)
            try:
                dt = json.load(open(f"{barcode_dir}/qcbams/{b}.{self.ref_name}.flagstats.json"))
                dt["barcode"] = b
                qcbams_dts.append(dt)
            except FileNotFoundError:
                continue
        df = pd.DataFrame(qcbams_dts)
        df_path = f"{self.expt_dirs.approach_dir}/summary.bam_flagstats.csv"
        df.to_csv(df_path, index=False)

    def _run_bedcov(self):
        """
        Summarise BED coverage results across all barcodes

        """
        bedcov_dfs = []
        for b in self.metadata.barcodes:
            barcode_dir = self.expt_dirs.get_barcode_dir(b)
            try:
                df = pd.read_csv(f"{barcode_dir}/bedcov/{b}.{self.ref_name}.bedcov.csv")
                bedcov_dfs.append(df)
            except FileNotFoundError:
                continue
        bedcov_df = pd.concat(bedcov_dfs)
        df_path = f"{self.expt_dirs.approach_dir}/summary.bedcov.csv"
        bedcov_df.to_csv(df_path, index=False)

    def _run_depth(self):
        """
        Summarise depth results across all barcodes

        """
        depth_dfs = []
        for b in self.metadata.barcodes:
            barcode_dir = self.expt_dirs.get_barcode_dir(b)
            try:
                df = pd.read_csv(f"{barcode_dir}/depth/{b}.{self.ref_name}.depth.csv")
                depth_dfs.append(df)
            except FileNotFoundError:
                continue
        depth_df = pd.concat(depth_dfs)
        df_path = f"{self.expt_dirs.approach_dir}/summary.depth.csv"
        depth_df.to_csv(df_path, index=False)

    def _run_variant(self):
        """
        Summarise variant calling results across all barcodes

        """
        variant_dfs = []
        for b in self.metadata.barcodes:
            barcode_dir = self.expt_dirs.get_barcode_dir(b)
            try:
                df = pd.read_csv(f"{barcode_dir}/vcfs/{b}.{self.ref_name}.annotated.tsv", sep="\t")
                df.insert(0, "barcode", b)
                variant_dfs.append(df)
            except FileNotFoundError:
                continue
        variant_df = pd.concat(variant_dfs)
        variant_path = f"{self.expt_dirs.approach_dir}/summary.variants.csv"
        variant_df.to_csv(variant_path, index=False)



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
    def run(self):
        self._run_fastq()
        self._run_qcbams()
        self._run_bedcov()
        self._run_depth()
        self._run_variant()

