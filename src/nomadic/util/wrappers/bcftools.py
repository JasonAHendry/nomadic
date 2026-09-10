import os
import shlex
import shutil
import subprocess
import uuid


def view(input_vcf, output_vcf, dry_run=False, **kwargs):
    """
    Run `bcftools view` on an `input_vcf`

    params:
        input_vcf: str
            Path to input `.vcf`.
        output_vcf: str
            Path to output `.vcf`.
        dry_run: bool
            Print command instead of running it.
        kwargs: key=value
            Additional arguments will be passed to
            bcftools view as flags; e.g.
            `-<key> <value>`.

    returns
        None

    """

    cmd = f"bcftools view -I {shlex.quote(input_vcf)} "
    cmd += " ".join([f"-{k} {v}" for k, v in kwargs.items()])
    cmd += f" -o {shlex.quote(output_vcf)}"

    if dry_run:
        print(cmd)
        return

    subprocess.run(cmd, shell=True, check=True)


def sort(input_vcf, output_vcf, dry_run=False, **kwargs):
    """
    Run `bcftools sort` on an `input_vcf`

    params:
        input_vcf: str
            Path to input `.vcf`.
        output_vcf: strbcft
            Path to output `.vcf`.
        dry_run: bool
            Print command instead of running it.
        kwargs: key=value
            Additional arguments will be passed to
            bcftools view as flags; e.g.
            `-<key> <value>`.

    returns
        None

    """

    cmd = f"bcftools sort {shlex.quote(input_vcf)} "
    cmd += " ".join([f"-{k} {v}" for k, v in kwargs.items()])
    cmd += f" -o {shlex.quote(output_vcf)}"

    if dry_run:
        print(cmd)
        return

    subprocess.run(cmd, shell=True, check=True)


def index(input_vcf):
    """
    Run `bcftools index`

    params
        input_vcf : str
            VCF file to be indexed.

    returns
        None

    """

    cmd = ["bcftools", "index", "-f", input_vcf]
    subprocess.run(cmd, check=True)


def merge(vcfs: list[str], output_vcf: str):
    """
    Run `bcftools merge` on a list of VCF files

    TODO: Could add arguments to this, but not strictly necessary
    """
    cmd = ["bcftools", "merge", "-F", "x", *vcfs, "-Oz", "-o", output_vcf]
    subprocess.run(cmd, check=True)


def reheader(input_vcf, output_vcf, sample_names):
    """
    Run `bcftools reheader` to change VCF sample name

    params
        input_vcf : str
            Input VCF.
        output_vcf : str
            Output VCF.
        sample_names : list of str
            Name of samples, in a list.

    returns
        None

    """

    # Generate random file name
    sample_file = f"sample_{uuid.uuid4()!s}.txt"

    # Write sample names to file
    with open(sample_file, "w") as fn:
        fn.writelines(f"{sample_name}\n" for sample_name in sample_names)

    # Create temporary file if input and output have same name
    if input_vcf == output_vcf:
        output_vcf = f"{input_vcf}".replace(".vcf", f"{uuid.uuid4()!s}.vcf")
        cleanup = True

    # Construct and run command
    cmd = ["bcftools", "reheader", input_vcf, "-s", sample_file, "-o", output_vcf]
    subprocess.run(cmd, check=True)

    # Remove file
    os.remove(sample_file)

    # Clean up, if necessary
    if cleanup:
        shutil.move(output_vcf, input_vcf)
