from typing import List

from nomadic.download.references import REFERENCE_COLLECTION
from nomadic.util.metadata import MetadataTableParser
from nomadic.util.experiment import ExperimentDirectories
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
        caller: str,
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

        self.caller = caller

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

        if self.caller:
            return BarcodeCallingPipelineRT(caller=self.caller, **kwargs)

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
        if self.caller:
            return ExptCallingPipelineRT(
                self.metadata, self.expt_dirs, self.regions, self.caller, self.reference
            )

        return ExptMappingPipelineRT(self.metadata, self.expt_dirs, self.ref_name)

    def get_dashboard(self, *, start_time=None) -> RealtimeDashboardBuilder:
        """
        Get the appropriate dashboard

        """
        summary_files = self.expt_dirs.get_summary_files()
        if self.caller:
            return CallingRTDashboard(
                expt_name=self.experiment_name,
                regions=self.regions,
                metadata=self.metadata,
                fastq_csv=summary_files.fastqs_processed,
                read_mapping_csv=summary_files.read_mapping,
                region_coverage_csv=summary_files.region_coverage,
                depth_profiles_csv=summary_files.depth_profiles,
                variant_csv=summary_files.variants,
                start_time=start_time,
                is_realtime=True,
            )

        return MappingRTDashboard(
            expt_name=self.experiment_name,
            regions=self.regions,
            metadata=self.metadata,
            fastq_csv=summary_files.fastqs_processed,
            read_mapping_csv=summary_files.read_mapping,
            region_coverage_csv=summary_files.region_coverage,
            depth_profiles_csv=summary_files.depth_profiles,
            start_time=start_time,
            is_realtime=True,
        )
