"""For saving and loading experiment settings"""

import json
import os
from datetime import datetime
from typing import NamedTuple, Optional


class ExperimentSettings(NamedTuple):
    name: str
    start_date: datetime
    call: bool
    fastq_dir: str
    metadata_csv: str
    region_bed: str
    reference_name: str
    n_barcodes: int
    n_regions: int
    minknow_dir: str = None


def parse_date(date_str: str) -> datetime:
    return datetime.fromisoformat(date_str)


def encode_date(date: datetime) -> str:
    return date.isoformat()


def verify_compatible_settings(
    old_settings: ExperimentSettings, new_settings: ExperimentSettings
):
    """Checks if the old and new settings of the experiment are compatible when resuming

    If an incompatible setting is found, raises an exception
    """
    if old_settings.name != new_settings.name:
        raise IncompatibleSettingsError("experiment name")
    if old_settings.fastq_dir != new_settings.fastq_dir:
        raise IncompatibleSettingsError("fastq dir")
    if (old_settings.minknow_dir != new_settings.minknow_dir) and (
        old_settings.minknow_dir is not None
    ):
        raise IncompatibleSettingsError("minknow dir")
    if old_settings.metadata_csv != new_settings.metadata_csv:
        raise IncompatibleSettingsError("metadata csv file")
    if old_settings.region_bed != new_settings.region_bed:
        raise IncompatibleSettingsError("region bed file")
    if old_settings.reference_name != new_settings.reference_name:
        raise IncompatibleSettingsError("reference name")
    if old_settings.n_barcodes != new_settings.n_barcodes:
        raise IncompatibleSettingsError("n barcodes")
    if old_settings.n_regions != new_settings.n_regions:
        raise IncompatibleSettingsError("n regions")
    return None


class IncompatibleSettingsError(Exception):
    def __init__(self, setting: str):
        self.setting = setting

    def __str__(self):
        return f"Incompatible setting: {self.setting}"


def save_settings(path: str, experiment_settings: ExperimentSettings):
    data = experiment_settings._asdict()
    data["start_date"] = encode_date(data["start_date"])
    with open(path, "w") as file:
        json.dump(data, file)


def load_settings(path: str) -> Optional[ExperimentSettings]:
    if not os.path.isfile(path):
        return None
    with open(path, "r") as file:
        data = json.load(file)
        data["start_date"] = parse_date(data["start_date"])
        return ExperimentSettings(**data)
