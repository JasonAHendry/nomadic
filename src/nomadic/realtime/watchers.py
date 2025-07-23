import os
import enum
from typing import NamedTuple, Optional
import hashlib

from nomadic.util import minknow
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
        barcode: str,
        fastq_dir: str,
        barcode_pipeline: BarcodePipelineRT,
        work_log_path: str,
    ):
        """
        Watch a barcode directory and run a BarcodePipeline `barcode_pipeline`
        when new FASTQ files are found
        """

        self.barcode = barcode
        self.fastq_dir = fastq_dir
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

        fastq_barcode_dir = self._fastq_barcode_dir()
        if fastq_barcode_dir is None:
            # Could not find fastq_dir
            return

        # In case the directory hasn't been created yet, i.e. no FASTQs for this barcode
        if not os.path.exists(fastq_barcode_dir):
            return None

        # If it has, collect the FASTQs
        observed_fastqs = set(
            [
                f"{fastq_barcode_dir}/{file}"
                for file in os.listdir(fastq_barcode_dir)
                if self._is_fastq(file)
            ]
        )

        # Are any unprocessed?
        return observed_fastqs.difference(self.processed_fastqs)

    def _fastq_barcode_dir(self) -> Optional[str]:
        if "*" in self.fastq_dir:
            # is the glob to find the fastq_dir
            fastq_dir = minknow.fastq_dir(fastq_dir_glob=self.fastq_dir)
        else:
            fastq_dir = self.fastq_dir

        if fastq_dir is None:
            # Dir is not yet created by minkown or wrong experiment name
            return None

        fastq_barcode_dir = os.path.join(fastq_dir, self.barcode)

        return fastq_barcode_dir

    def catch_up_from_work_log(self) -> tuple[int, int]:
        """
        Read the log file and set all FASTQs that have been processed before
        """
        if not os.path.exists(self.work_log_path):
            return 0, 0

        barcode_fastq_dir = self._fastq_barcode_dir()
        if barcode_fastq_dir is None:
            raise RuntimeError("barcode fastq dir does not exist")

        started_and_not_ended = dict()
        with open(self.work_log_path, "r") as log_file:
            number_fully_processed_fastqs = 0
            for line in log_file:
                # TODO verify that files match via hash
                incr_id, action, fastqs = parse_log_line(line.strip())
                # Filename in log is just the basename, so we need to prepend the directory
                fastqs = [
                    f"{barcode_fastq_dir}/{fastq}"
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
        for incr_info in started_and_not_ended.values():
            # If we have START but no END, processing has somewhere failed in the middle, so reprocess them
            number_reprocessed_fastqs += len(incr_info.files)
            self._process_fastq(incr_info.files)

        return number_fully_processed_fastqs, number_reprocessed_fastqs

    def update(self) -> bool:
        """
        Run the `barcode_pipeline` if any FASTQ files are currently unprocessed
        """

        unprocessed_fastqs = self._check_fastqs()

        if not unprocessed_fastqs:
            return False

        self._process_fastq(list(unprocessed_fastqs))

        return True

    def _process_fastq(self, fastqs: list[str]):
        incr_id = increment_id(fastqs)
        self.write_log_start(fastqs, incr_id)
        self.barcode_pipeline.run(fastqs, incr_id)
        self.write_log_end(fastqs, incr_id)

        self.processed_fastqs.update(fastqs)

    def write_log_start(self, fastqs: list[str], incr_id: str):
        """
        Write a log entry when the watcher starts processing FASTQs
        :param fastqs: List of FASTQ files being processed
        :param step: Unique step identifier for the log entry
        """

        with open(self.work_log_path, "a") as log_file:
            log_file.write(format_log_line(Action.START, incr_id, fastqs))

    def write_log_end(self, fastqs: list[str], incr_id: str):
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
