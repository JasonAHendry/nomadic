import os
import enum
from typing import NamedTuple
import hashlib

from .pipelines.barcode import BarcodePipelineRT


class IncrementInfo(NamedTuple):
    id: str
    files: list[str]


class Action(enum.Enum):
    START = "START"
    END = "END"


class BarcodeWatcher:
    def __init__(
        self,
        barcode_fastq_dir: str,
        barcode_pipeline: BarcodePipelineRT,
        work_log_path: str,
    ):
        """
        Watch a barcode directory `barcode_fastq_dir` and run a BarcodePipeline `barcode_pipeline`
        when new FASTQ files are found
        """

        self.barcode_fastq_dir = barcode_fastq_dir
        self.barcode_pipeline = barcode_pipeline
        # path to work log file, used to recover when application is restarted
        self.work_log_path = work_log_path

        self.processed_fastqs: set[str] = set()

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
        return observed_fastqs.difference(self.processed_fastqs)

    def catch_up_from_work_log(self) -> tuple[int, int]:
        """
        Read the log file and set all FASTQs that have been processed before
        """
        if not os.path.exists(self.work_log_path):
            return 0, 0

        started_and_not_ended = dict()
        with open(self.work_log_path, "r") as log_file:
            number_fully_processed_fastqs = 0
            for line in log_file:
                # TODO verify that files match via hash
                incr_id, action, fastqs = parse_log_line(line.strip())
                # Filename in log is just the basename, so we need to prepend the directory
                fastqs = [
                    f"{self.barcode_fastq_dir}/{fastq}"
                    for fastq in fastqs
                    if self._is_fastq(fastq)
                ]
                if action is Action.START:
                    # Multiple starts for the same increment ID are ok and will be just overriden, e.g. if the application is restarded
                    started_and_not_ended[incr_id] = IncrementInfo(incr_id, fastqs)
                elif action is Action.END:
                    if incr_id not in started_and_not_ended:
                        raise RuntimeError(
                            f"Found END for step {incr_id} without START"
                        )
                    del started_and_not_ended[incr_id]
                    number_fully_processed_fastqs += len(fastqs)
                    self.processed_fastqs.update(fastqs)
                else:
                    raise RuntimeError(f"Unknown action in log: {action}")

        number_reprocessed_fastqs = 0
        for step_info in started_and_not_ended.values():
            # If we have START but no END, processing has somewhere failed in the middle, so reprocess them
            number_reprocessed_fastqs += len(step_info.files)
            self.update_fastq(step_info.files)

        return number_fully_processed_fastqs, number_reprocessed_fastqs

    def update(self) -> bool:
        """
        Run the `barcode_pipeline` if any FASTQ files are currently unprocessed
        """

        unprocessed_fastqs = self._check_fastqs()

        if not unprocessed_fastqs:
            return False

        self.update_fastq(list(unprocessed_fastqs))

        return True

    def update_fastq(self, fastqs: list[str]):
        incr_id = increment_id(fastqs)
        self.write_log_start(fastqs, incr_id)
        self.barcode_pipeline.run(fastqs, incr_id)
        self.write_log_finish(fastqs, incr_id)

        self.processed_fastqs.update(fastqs)

    def write_log_start(self, fastqs: list[str], incr_id: str):
        """
        Write a log entry when the watcher starts processing FASTQs
        :param fastqs: List of FASTQ files being processed
        :param step: Unique step identifier for the log entry
        """

        with open(self.work_log_path, "a") as log_file:
            log_file.write(format_log_line(Action.START, incr_id, fastqs))

    def write_log_finish(self, fastqs: list[str], incr_id: str):
        """
        Write a log entry when the watcher finished processing FASTQs
        """

        with open(self.work_log_path, "a") as log_file:
            log_file.write(format_log_line(Action.END, incr_id, fastqs))


def increment_id(fastqs: list[str]) -> str:
    """
    Generate a unique increment ID based on the sorted list of FASTQ files that are processed for this increment.
    """
    return hashlib.md5(
        "".join(sorted(os.path.basename(fastq) for fastq in fastqs)).encode()
    ).hexdigest()


def format_log_line(action: Action, incr_id: str, fastqs: list[str]) -> str:
    return f"{incr_id}\t{action.value}\t{','.join(os.path.basename(fastq) for fastq in fastqs)}\n"


def parse_log_line(line: str) -> tuple[str, Action, list[str]]:
    """
    Parse a log line into its components
    """
    parts = line.strip().split("\t")
    if len(parts) != 3:
        raise RuntimeError("Invalid log line format")

    incr_id = parts[0]
    action = Action(parts[1])
    fastqs = parts[2].split(",") if parts[2] else []

    return incr_id, action, fastqs
