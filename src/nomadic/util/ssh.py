import re
import subprocess
import shlex

import click


def is_ssh_target(s: str) -> bool:
    # matches user@host:/path or host:/path
    return re.match(r"^(?:[^@:\\s]+@)?[^@:\\s]+:/", s) is not None


def split_ssh_target(t: str) -> tuple[str, str]:
    # split into 'user@host' and '/absolute/path'
    idx = t.find(":")
    host = t[:idx]
    path = t[idx + 1 :]
    return host, path


def create_remote_dir(ssh_target: str, verbose: bool = False) -> tuple[bool, str]:
    host, path = split_ssh_target(ssh_target)
    # make sure path is absolute
    if not path.startswith("/"):
        raise ValueError(f"Remote path '{path}' is not absolute")
    path_q = shlex.quote(path.rstrip("/"))
    # create and verify directory on remote
    cmd = f"mkdir -p {path_q} && test -d {path_q}"
    return remote_command(host, cmd, verbose=verbose)


def remote_dir_exists(ssh_target: str, verbose: bool = False) -> tuple[bool, str]:
    host, path = split_ssh_target(ssh_target)
    # make sure path is absolute
    if not path.startswith("/"):
        raise ValueError(f"Remote path '{path}' is not absolute")
    path_q = shlex.quote(path.rstrip("/"))
    # verify directory on remote
    cmd = f"test -d {path_q}"
    return remote_command(host, cmd, verbose=verbose)


def remote_command(host: str, command: str, verbose: bool = False) -> tuple[bool, str]:
    if verbose:
        click.echo(f"Running: ssh {host} {shlex.quote(command)}")
    proc = subprocess.run(["ssh", host, command], text=True, capture_output=True)
    if proc.returncode != 0:
        msg = proc.stderr.strip() or f"ssh exited {proc.returncode}"
        return False, msg
    return True, proc.stdout.strip()
