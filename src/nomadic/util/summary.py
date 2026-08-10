from pathlib import Path


def looks_like_summary_dir(path: Path) -> bool:
    """
    Check if the given path looks like a summary directory.
    """
    required_files = [
        "inventory.csv",
    ]
    return all((path / f).exists() for f in required_files)
