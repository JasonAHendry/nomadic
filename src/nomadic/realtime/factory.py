from typing import List

from nomadic.download.references import REFERENCE_COLLECTION
from nomadic.util.metadata import MetadataTableParser
from nomadic.util.dirs import ExperimentDirectories
from nomadic.util.regions import RegionBEDParser

from nomadic.realtime.watchers import BarcodeWatcher
from nomadic.realtime.pipelines.barcode import (
    BarcodePipelineRT,
    BarcodeMappingPipelineRT,
    BarcodeCallingPipelineRT,
)
from nomadic.realtime.pipelines.experiment import (
    ExperimentPipelineRT,
    ExptMappingPipelineRT,
    ExptCallingPipelineRT,
)
from nomadic.realtime.dashboard.builders import (
    RealtimeDashboardBuilder,
    MappingRTDashboard,
    CallingRTDashboard,
)


class PipelineFactory:
    """
    Given users inputs select the correct versions of the
    barcode pipline, experiment pipeline and dashboard

    """

    def __init__(
        self,
        experiment_name: str,
        metadata: MetadataTableParser,
        regions: RegionBEDParser,
        expt_dirs: ExperimentDirectories,
        fastq_dir: str,
        call: bool = False,
        ref_name: str = "Pf3D7",
    ):
        """
        Store metadata as instance attributes

        """

        self.experiment_name = experiment_name
        self.metadata = metadata
        self.regions = regions
        self.expt_dirs = expt_dirs
        self.fastq_dir = fastq_dir

        if ref_name not in REFERENCE_COLLECTION:
            raise ValueError(
                f"Reference {ref_name} must be in: {','.join(REFERENCE_COLLECTION)}."
            )
        self.ref_name = ref_name
        self.reference = REFERENCE_COLLECTION[ref_name]

        if not isinstance(call, bool):
            raise ValueError("`call` must be a boolean.")

        self.call = call

    def _get_barcode_pipeline(self, barcode_name: str) -> BarcodePipelineRT:
        """
        Get the appropriate barcode pipeline

        """
        kwargs = {
            "barcode_name": barcode_name,
            "expt_dirs": self.expt_dirs,
            "bed_path": self.regions.path,
            "ref_name": self.ref_name,
        }

        if self.call:
            return BarcodeCallingPipelineRT(**kwargs)

        return BarcodeMappingPipelineRT(**kwargs)

    def get_watchers(self) -> List[BarcodeWatcher]:
        """
        Initialise watchers for each barcode, and return them

        """
        return [
            BarcodeWatcher(
                barcode=b,
                fastq_dir=self.fastq_dir,
                barcode_pipeline=self._get_barcode_pipeline(barcode_name=b),
                work_log_path=f"{self.expt_dirs.get_barcode_dir(b)}/.work.log",
            )
            for b in self.metadata.barcodes
        ]

    def get_expt_pipeline(self) -> ExperimentPipelineRT:
        """
        Get the appropriate experiment pipeline

        """
        if self.call:
            return ExptCallingPipelineRT(
                self.metadata, self.expt_dirs, self.regions, self.reference
            )

        return ExptMappingPipelineRT(self.metadata, self.expt_dirs, self.ref_name)

    def get_dashboard(self) -> RealtimeDashboardBuilder:
        """
        Get the appropriate dashboard

        """
        if self.call:
            return CallingRTDashboard(
                expt_name=self.experiment_name,
                regions=self.regions,
                metadata=self.metadata,
                fastq_csv=f"{self.expt_dirs.approach_dir}/summary.fastq.csv",
                flagstats_csv=f"{self.expt_dirs.approach_dir}/summary.bam_flagstats.csv",
                bedcov_csv=f"{self.expt_dirs.approach_dir}/summary.bedcov.csv",
                depth_csv=f"{self.expt_dirs.approach_dir}/summary.depth.csv",
                variant_csv=f"{self.expt_dirs.approach_dir}/summary.variants.csv",
            )

        return MappingRTDashboard(
            expt_name=self.experiment_name,
            regions=self.regions,
            metadata=self.metadata,
            fastq_csv=f"{self.expt_dirs.approach_dir}/summary.fastq.csv",
            flagstats_csv=f"{self.expt_dirs.approach_dir}/summary.bam_flagstats.csv",
            bedcov_csv=f"{self.expt_dirs.approach_dir}/summary.bedcov.csv",
            depth_csv=f"{self.expt_dirs.approach_dir}/summary.depth.csv",
        )
