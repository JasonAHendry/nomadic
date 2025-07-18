import pytest

from .mappers import encode_input_files


@pytest.mark.parametrize(
    "files,expected",
    [
        ([], "-"),
        (["test"], "test"),
        (["test1", "test2"], "test1 test2"),
        (["test test1", "test test2", "test3"], "'test test1' 'test test2' test3"),
        (["../../test1", "../../test2"], "../../test1 ../../test2"),
    ],
)
def test_encode_input_files(files, expected):
    assert encode_input_files(files) == expected
