import pytest

from nomadic.util.ssh import is_ssh_target, split_ssh_target


@pytest.mark.parametrize(
    ["text", "expected"],
    [
        ["", False],
        ["abc", False],
        ["/root/test", False],
        ["user@", False],
        ["user@server", False],
        ["user@server:/", True],
        ["server:/", True],
        ["user@server:/root/test", True],
        ["server:/root/test", True],
        ["server:/home/user/test", True],
        ["server:/home:/user/test", False],
    ],
)
def test_ssh_target(text, expected):
    assert is_ssh_target(text) == expected


@pytest.mark.parametrize(
    ["text", "expected"],
    [
        # not really a sane input, but let's make sure it splits at the first
        ["server:/test/:/test", ("server", "/test/:/test")],
        ["server:/", ("server", "/")],
        ["user@server:/root/test", ("user@server", "/root/test")],
        ["server:/root/test", ("server", "/root/test")],
        ["server:/home/user/test", ("server", "/home/user/test")],
    ],
)
def test_split_ssh_target(text, expected):
    assert split_ssh_target(text) == expected
