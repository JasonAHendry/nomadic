import subprocess


def test_imports():
    """Ensure that no heavy libraries are imported in cli.py to keep startup time low."""
    output = subprocess.check_output(
        ["python", "-Ximporttime", "-mnomadic.cli"], stderr=subprocess.STDOUT
    )
    output = output.decode("utf-8")

    libs = ["pandas", "numpy", "seaborn"]
    for lib in libs:
        assert lib not in output, (
            f"{lib} should not be imported in cli.py because it slows down the startup time"
        )
