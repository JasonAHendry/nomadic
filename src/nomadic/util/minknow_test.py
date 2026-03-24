import glob
from pathlib import Path

from nomadic.util.minknow import (
    is_fastq_dir,
    is_minknow_base_dir,
    is_minknow_experiment_dir,
    resolve_minknow_fastq_dirs,
)

test_folder = Path("src/nomadic/util/_test_data")


def test_is_minknow_base_dir_real_minknow_dir():
    assert is_minknow_base_dir(test_folder / "minknow" / "minknow_folder")


def test_is_minknow_base_dir_contains_minknow_experiments():
    assert is_minknow_base_dir(test_folder / "minknow" / "minknow_experiments")


def test_is_minknow_base_dir_empty():
    assert not is_minknow_base_dir(test_folder / "minknow" / "empty")


def test_is_minknow_base_dir_not_exist():
    assert not is_minknow_base_dir(test_folder / "minknow" / "does_not_exist")


def test_is_minknow_experiment_dir_not_exist():
    assert not is_minknow_experiment_dir(
        test_folder / "minknow" / "minknow_experiments" / "does_not_exist"
    )


def test_is_minknow_experiment_dir():
    assert is_minknow_experiment_dir(
        test_folder / "minknow" / "minknow_experiments" / "experiment1"
    )


def test_resolve_minknow_fastq_dirs_from_base_dir():
    minknow_dir, fastq_dir = resolve_minknow_fastq_dirs(
        test_folder / "minknow" / "minknow_experiments", "experiment1"
    )
    assert str(minknow_dir).endswith("/minknow/minknow_experiments/experiment1")
    assert fastq_dir.endswith("/minknow/minknow_experiments/experiment1/*/*/fastq_pass")
    assert len(glob.glob(fastq_dir)) == 1


def test_is_fastq_dir():
    assert is_fastq_dir(
        test_folder
        / "minknow"
        / "minknow_experiments"
        / "experiment1"
        / "Sample"
        / "1234-ABC-123"
        / "fastq_pass"
    )


def test_is_fastq_dir_empty():
    assert not is_fastq_dir(test_folder / "minknow" / "empty")
