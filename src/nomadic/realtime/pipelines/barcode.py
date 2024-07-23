from datetime import datetime
from typing import List
from abc import ABC, abstractmethod
from nomadic.util.dirs import ExperimentDirectories
from nomadic.util.regions import RegionBEDParser
from nomadic.realtime.steps import (
    FASTQCountRT,
    MappingRT,
    FlagstatsRT,
    BedCovRT,
    RegionDepthRT,
    CallVariantsRT,
    AnnotateVariantsRT
)

    
# --------------------------------------------------------------------------------
# Interface for barcode pipelines
#
# TODO:
# - Could implement this as a formal DAG to allow for more complex pipeline
#   step relationships
# - Could have a concrete def _run() that performs steps necessarily shared
#   by all concret pipelines
#
# --------------------------------------------------------------------------------


class BarcodePipelineRT(ABC):
    """
    Interface for per-barcode real-time analysis pipelines

    """

    def __init__(self, barcode_name: str, expt_dirs: ExperimentDirectories, bed_path: str, ref_name: str="Pf3D7"):
        """
        Store important metadata as instance attributes, and define steps

        """

        self.barcode_name = barcode_name
        self.expt_dirs = expt_dirs
        self.barcode_dir = expt_dirs.get_barcode_dir(barcode_name)
        self.ref_name = ref_name

    @abstractmethod
    def _run(self, new_fastq: List[str]) -> None:
        """
        Run analysis steps for the pipeline on new FASTQ files

        """
        pass

    def run(self, new_fastq: List[str]) -> None:
        """
        Wrapper to try and make the stdout easier to read

        TODO: Could log this

        """
        
        t0 = datetime.now().replace(microsecond=0)
        print("")
        print("="*80)
        print(f"Updating: {self.barcode_name}")
        print("-"*80)

        self._run(new_fastq=new_fastq)
        
        t1 = datetime.now().replace(microsecond=0)
        print("-"*80)
        print(f"Done with: {self.barcode_name}")
        print(f"Time elapsed: {t1 - t0}")
        print("="*80)


# --------------------------------------------------------------------------------
# Concrete pipelines
#
# --------------------------------------------------------------------------------


class BarcodeMappingPipelineRT(BarcodePipelineRT):
    """
    Concrete barcode pipeline for mapping reads in real-time, and performing
    basic quality control

    """

    def __init__(self, barcode_name: str, expt_dirs: ExperimentDirectories, bed_path: str, ref_name: str="Pf3D7"):
        """
        Initialise the barcode mapping pipeline

        """
        
        super().__init__(barcode_name, expt_dirs, bed_path)
        
        # Initialise analysis steps
        common = {"barcode_name": barcode_name, "expt_dirs": expt_dirs}
        self.fastq_step = FASTQCountRT(**common)
        self.map_step = MappingRT(**common, ref_name=ref_name)
        self.flagstat_step = FlagstatsRT(**common, ref_name=ref_name)
        self.bedcov_step = BedCovRT(**common, bed_path=bed_path, ref_name=ref_name)
        self.depth_step = RegionDepthRT(**common, regions=RegionBEDParser(bed_path), ref_name=ref_name)

    def _run(self, new_fastq: List[str]) -> None:
        """
        Run mapping and QC from a set of newly generated FASTQ files
        
        """

        inter_bam = self.map_step.run(new_fastq)
        final_bam = self.map_step.merge()

        self.flagstat_step.run(inter_bam)
        self.flagstat_step.merge()

        self.bedcov_step.run(final_bam)
        self.bedcov_step.merge()

        self.depth_step.run(final_bam)
        self.depth_step.merge()

        self.fastq_step.run(new_fastq)


class BarcodeCallingPipelineRT(BarcodePipelineRT):
    """
    Concrete barcode pipeline for mapping reads in real-time, and performing
    basic quality control

    """

    def __init__(self, barcode_name: str, expt_dirs: ExperimentDirectories, bed_path: str, ref_name: str="Pf3D7"):
        """
        Initialise the barcode mapping pipeline

        """
        
        super().__init__(barcode_name, expt_dirs, bed_path)
        
        # Initialise analysis steps
        common = {"barcode_name": barcode_name, "expt_dirs": expt_dirs}
        self.fastq_step = FASTQCountRT(**common)
        self.map_step = MappingRT(**common, ref_name=ref_name)
        self.flagstat_step = FlagstatsRT(**common, ref_name=ref_name)
        self.bedcov_step = BedCovRT(**common, bed_path=bed_path, ref_name=ref_name)
        self.depth_step = RegionDepthRT(**common, regions=RegionBEDParser(bed_path), ref_name=ref_name)
        self.call_step = CallVariantsRT(**common, bed_path=bed_path, ref_name=ref_name)
        #self.annot_step = AnnotateVariantsRT(**common, bed_path=bed_path, ref_name=ref_name)


    def _run(self, new_fastq: List[str]) -> None:
        """
        Run mapping and QC from a set of newly generated FASTQ files
        
        """

        inter_bam = self.map_step.run(new_fastq)
        final_bam = self.map_step.merge()

        self.flagstat_step.run(inter_bam)
        self.flagstat_step.merge()

        self.bedcov_step.run(final_bam)
        self.bedcov_step.merge()

        self.depth_step.run(final_bam)
        self.depth_step.merge()

        final_vcf = self.call_step.run(final_bam)
        self.call_step.merge()

        # self.annot_step.run(final_vcf)
        # self.annot_step.merge()

        self.fastq_step.run(new_fastq)
