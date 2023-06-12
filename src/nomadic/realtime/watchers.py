import os
from .pipelines.barcode import BarcodePipelineRT


class BarcodeWatcher:
    def __init__(self, barcode_fastq_dir: str, barcode_pipeline: BarcodePipelineRT):
        """
        Watch a barcode directory `barcode_fastq_dir` and run a BarcodePipeline `barcode_pipeline`
        when new FASTQ files are found
        """

        self.barcode_fastq_dir = barcode_fastq_dir
        self.barcode_pipeline = barcode_pipeline

        self.processed_fastqs = set()
        self.unprocessed_fastqs = set()

    @staticmethod
    def _is_fastq(file_name: str) -> bool:
        """
        Is it a FASTQ file?
        """
        return file_name.endswith(".fastq") or file_name.endswith(".fastq.gz")

    def _check_fastqs(self):
        """
        Check if any new FASTQ files have been generated
        """

        # In case the directory hasn't been created yet, i.e. no FASTQs
        if not os.path.exists(self.barcode_fastq_dir):
            return

        # If it has, collect the FASTQs
        observed_fastqs = set(
            [
                f"{self.barcode_fastq_dir}/{file}"
                for file in os.listdir(self.barcode_fastq_dir)
                if self._is_fastq(file)
            ]
        )

        # Are any unprocessed?
        self.unprocessed_fastqs = observed_fastqs.difference(self.processed_fastqs)

    def _all_fastqs_processed(self):
        """
        Set all FASTQ files to processed
        """

        self.processed_fastqs.update(self.unprocessed_fastqs)

    def update(self) -> bool:
        """
        Run the `barcode_pipeline` if any FASTQ files are currently unprocessed
        """

        self._check_fastqs()

        if not self.unprocessed_fastqs:
            return False

        # TODO: What happens if the pipeline *fails*?
        self.barcode_pipeline.run(
            self.unprocessed_fastqs
        )  # note, we take the new FASTQ files
        self._all_fastqs_processed()

        return True
