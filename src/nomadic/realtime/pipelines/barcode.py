from typing import List

from nomadic.util.dirs import ExperimentDirectories
from nomadic.realtime.steps import MappingRT, FlagstatsRT, BedCovRT


class BarcodePipelineRT:
    """
    Run an analysis pipeline for a barcode in real-time

    """

    def __init__(self, barcode_name: str, 
                 expt_dirs: ExperimentDirectories,
                 bed_path: str
                 ):
        """
        Initialise the barcode pipeline

        """

        # Parse user inputs
        self.barcode_name = barcode_name
        self.expt_dirs = expt_dirs
        self.barcode_dir = expt_dirs.get_barcode_dir(barcode_name)

        # Initialise analysis steps
        self.map_step = MappingRT(barcode_name, expt_dirs)
        self.flagstat_step = FlagstatsRT(barcode_name, expt_dirs)
        self.bedcov_step = BedCovRT(
            barcode_name, 
            expt_dirs, 
            bed_path
        )

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

        # self.depth_step.run(inter_bam)
        # self.depth_step.merge()
