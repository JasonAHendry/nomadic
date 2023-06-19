import json
import pandas as pd

from abc import ABC, abstractmethod

from nomadic.util.metadata import MetadataTableParser
from nomadic.util.dirs import ExperimentDirectories


class ExperimentPipeline(ABC):
    """
    Interface for complete experiment pipeline, summarising
    information across multiple barcodes

    """

    def __init__(self, metadata: MetadataTableParser, expt_dirs: ExperimentDirectories):
        """
        Store experiment directories and metadata information

        """

        self.expt_dirs = expt_dirs
        self.metadata = metadata

    @abstractmethod
    def run(self):
        pass


class BasicExptPipeline(ExperimentPipeline):
    def run(self):
        """
        Summarise the number of fastqs from all barcodes into a dataframe

        """
        print("Running basic pipeline...")
        fastq_dts = []
        for b in self.metadata.barcodes:
            barcode_dir = self.expt_dirs.get_barcode_dir(b)
            try:
                dt = json.load(open(f"{barcode_dir}/n_processed_fastq.json", "r"))
                fastq_dts.append(dt)
            except FileNotFoundError:
                continue
        df = pd.DataFrame(fastq_dts)
        df_path = f"{self.expt_dirs.approach_dir}/summary.fastq.csv"
        df.to_csv(df_path, index=False)


class MappingInRTExptPipeline(ExperimentPipeline):
    """
    For each analysis step for each barcode, there is ONE key output file,
    the path of which we can get from an abstract method

    Will write several functions describing different types of merging operations, which will
    help here

    Then just get names and merge, write


    """

    def run(self):
        """
        Summarise the number of fastqs from all barcodes into a dataframe

        """

        # MERGE FASTQ stuff
        print("Running basic pipeline...")
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

        # MERGE JSON STUFF
        qcbams_dts = []
        for b in self.metadata.barcodes:
            barcode_dir = self.expt_dirs.get_barcode_dir(b)
            try:
                dt = json.load(open(f"{barcode_dir}/qcbams/{b}.Pf3D7.flagstats.json"))
                dt["barcode"] = b
                qcbams_dts.append(dt)
            except FileNotFoundError:
                continue

        df = pd.DataFrame(qcbams_dts)
        df_path = f"{self.expt_dirs.approach_dir}/summary.bam_flagstats.csv"
        df.to_csv(df_path, index=False)

        # CONCAT DFS
        bedcov_dfs = []
        for b in self.metadata.barcodes:
            barcode_dir = self.expt_dirs.get_barcode_dir(b)
            try:
                df = pd.read_csv(f"{barcode_dir}/{b}.Pf3D7.bedcov.csv")
                bedcov_dfs.append(df)
            except FileNotFoundError:
                continue

        bedcov_df = pd.concat(bedcov_dfs)
        df_path = f"{self.expt_dirs.approach_dir}/summary.bedcov.csv"
        bedcov_df.to_csv(df_path, index=False)


        # CONCAT DEPTH DFS
        depth_dfs = []
        for b in self.metadata.barcodes:
            barcode_dir = self.expt_dirs.get_barcode_dir(b)
            try:
                df = pd.read_csv(f"{barcode_dir}/depth/{b}.depth.csv")
                depth_dfs.append(df)
            except FileNotFoundError:
                continue
        depth_df = pd.concat(depth_dfs)
        df_path = f"{self.expt_dirs.approach_dir}/summary.depth.csv"
        depth_df.to_csv(df_path, index=False)
