import logging
import subprocess
from pathlib import Path

import click

from nomadic.util.dirs import produce_dir
from nomadic.util.logging_config import LoggingFascade
from nomadic.util.rsync import selective_rsync

script_dir = Path(__file__).parent.resolve()
log = logging.getLogger(script_dir.stem)


@click.command(short_help="Backup sequence data folder to a local hard disk drive")
@click.option(
    "-b",
    "--backup_dir",
    type=click.Path(exists=True, file_okay=False, dir_okay=True, path_type=Path),
    required=True,
    help="Path to backup folder on external USB drive",
)
@click.option(
    "-m",
    "--minknow_dir",
    type=click.Path(exists=True, file_okay=False, dir_okay=True, path_type=Path),
    default="/var/lib/minknow/data",
    show_default=True,
    help="Path to general output directory of minknow (e.g. `/var/lib/minknow/data`)",
)
def backup(backup_dir: Path, minknow_dir: Path) -> None:
    """
    Backup all sequence data files to a local USB drive.

    """
    # TODO: Unsure how to do logging properly in nomadic
    log = LoggingFascade(logger_name="nomadic")
    # Add in current nomadic workspace
    nomadic_dir = Path.cwd()

    print(f"Backing up nomadic data from {nomadic_dir} to {backup_dir}")
    selective_rsync(
        source_dir=nomadic_dir,
        target_dir=backup_dir,
        recursive=True,
        progressbar=True,
    )
    # TODO: Change print statement to log
    print(
        f"Backing up minknow data from {minknow_dir} and merging into experiment folders in {backup_dir}"
    )
    # Need to identify experimental folders ONLY
    exp_dirs = [f.name for f in (nomadic_dir / "results").iterdir() if f.is_dir()]

    for folder in exp_dirs:
        source_dir = minknow_dir / folder
        target_dir = backup_dir / "results" / folder / "minknow"
        if not target_dir.exists():
            produce_dir(target_dir)
        print(f"   {source_dir} to {target_dir}")
        selective_rsync(
            source_dir=source_dir,
            target_dir=target_dir,
            recursive=True,
            progressbar=True,
        )

    # Ensure that folders are owned by the correct user
    user_id = subprocess.check_output(["id", "-u"]).decode("utf-8").strip()
    group_id = subprocess.check_output(["id", "-g"]).decode("utf-8").strip()
    subprocess.run(["chown", "-R", f"{user_id}:{group_id}", str(backup_dir)])
    print("Backup complete.")
