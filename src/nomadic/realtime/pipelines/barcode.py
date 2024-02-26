from typing import List

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


class BarcodePipelineRT:
    """
    Run an analysis pipeline for a barcode in real-time

    """

    def __init__(
        self, barcode_name: str, expt_dirs: ExperimentDirectories, bed_path: str
    ):
        """
        Initialise the barcode pipeline

        """

        # Parse user inputs
        self.barcode_name = barcode_name
        self.expt_dirs = expt_dirs
        self.barcode_dir = expt_dirs.get_barcode_dir(barcode_name)

        # Initialise analysis steps
        common = {"barcode_name": barcode_name, "expt_dirs": expt_dirs}
        self.fastq_step = FASTQCountRT(**common)
        self.map_step = MappingRT(**common)
        self.flagstat_step = FlagstatsRT(**common)
        self.bedcov_step = BedCovRT(**common, bed_path)
        self.depth_step = RegionDepthRT(**common, RegionBEDParser(bed_path))
        self.annot_step = AnnotateVariantsRT(**common)

        # HERE
        self.call_step = CallVariantsRT(**common)
        #self.annot_step = None #AnnotateVariantsRT()

    def run(self, new_fastq):
        """
        Run pipeline
        """

        inter_bam = self.map_step.run(new_fastq)
        final_bam = self.map_step.merge()

        self.flagstat_step.run(inter_bam)
        self.flagstat_step.merge()

        self.bedcov_step.run(final_bam)
        self.bedcov_step.merge()

        self.depth_step.run(final_bam)
        self.depth_step.merge()

        final_vcf = self.call_step.run(final_bam) # pileup, call, filter, fill tags
        self.call_step.merge()

        self.annot_step.run(final_vcf)
        self.annot_step.merge()

        self.fastq_step.run(new_fastq)
