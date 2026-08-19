from pathlib import Path
import subprocess

import pytest

from nomadic.download.references import PlasmodiumFalciparum3D7
from nomadic.util.vcf import VariantAnnotator

vcf_test_data_dir = Path("src/nomadic/util/_test_data/vcfs")
snapshot_dir = Path("src/nomadic/util/_snapshots")
bed_test_data_dir = Path("src/nomadic/util/_test_data/beds")


@pytest.fixture
def vcf_annotator():
    reference = PlasmodiumFalciparum3D7()

    bed_path = bed_test_data_dir / "nomadsMVP.amplicons.bed"
    annotator = VariantAnnotator(
        fasta_path=reference.fasta_path,
        gff_path=reference.gff_path,
        bed_path=str(bed_path),
        caller="delve",
    )
    return annotator


def test_annotate_regions(vcf_annotator):
    vcf_path = vcf_test_data_dir / "delve" / "01-a-single-mut-mis.vcf"

    cmd = vcf_annotator._annotate_command(input_vcf=str(vcf_path))

    output = subprocess.check_output(cmd, shell=True).decode("utf-8")

    header_line = (
        "##INFO=<ID=AMP_ID,Number=1,Type=String,Description=Amplicon identifier>"
    )

    has_header = any(line.startswith(header_line) for line in output.splitlines())
    assert has_header, "Annotated VCF does not contain the expected header line."

    variant_lines = [line for line in output.splitlines() if not line.startswith("#")]

    assert len(variant_lines) == 1, "Annotated VCF does not contain 1 variant lines."

    info_field = variant_lines[0].split("\t")[7]

    assert "AMP_ID=dhfr-p1-410" in info_field, (
        "Annotated VCF does not contain the expected AMP_ID in the INFO field."
    )


def test_wsaf_command(vcf_annotator):
    vcf_path = vcf_test_data_dir / "delve" / "01-a-single-mut-mis.vcf"

    cmd = vcf_annotator._fill_wsaf_command(input_vcf=str(vcf_path))

    output = subprocess.check_output(cmd, shell=True).decode("utf-8")

    variant_lines = [line for line in output.splitlines() if not line.startswith("#")]

    assert len(variant_lines) == 1, "Annotated VCF does not contain 1 variant lines."

    format_field = variant_lines[0].split("\t")[8]
    assert "WSAF" in format_field, (
        "Annotated VCF does not contain the expected WSAF in the FORMAT field."
    )

    sample_field = variant_lines[0].split("\t")[9]

    sample_values = sample_field.split(":")

    assert sample_values[format_field.split(":").index("WSAF")] == "0.99", (
        "Annotated VCF does not contain the expected WSAF value in the sample field."
    )


def test_csq_command(vcf_annotator):
    vcf_path = vcf_test_data_dir / "delve" / "01-a-single-mut-mis.vcf"

    cmd = vcf_annotator._csq_command(input_vcf=str(vcf_path))

    output = subprocess.check_output(cmd, shell=True).decode("utf-8")

    variant_lines = [line for line in output.splitlines() if not line.startswith("#")]

    assert len(variant_lines) == 1, "Annotated VCF does not contain 1 variant lines."

    info_field = variant_lines[0].split("\t")[7]
    assert "BCSQ" in info_field, (
        "Annotated VCF does not contain the expected BCSQ in the INFO field."
    )

    bcsq_field = [
        field for field in info_field.split(";") if field.startswith("BCSQ=")
    ][0]
    assert (
        bcsq_field
        == "BCSQ=missense|DHFR-TS|PF3D7_0417200.1|protein_coding|+|51N>51I|748239A>T"
    ), "Annotated VCF does not contain the expected BCSQ value in the INFO field."

    format_field = variant_lines[0].split("\t")[8]
    assert "BCSQ" in format_field, (
        "Annotated VCF does not contain the expected CSQ in the FORMAT field."
    )


@pytest.mark.parametrize(
    "annotated_vcf",
    (vcf_test_data_dir / "annotated").glob("*.vcf"),
    ids=lambda p: p.stem,
)
def test_aa_changes(vcf_annotator: VariantAnnotator, annotated_vcf: Path, snapshot):
    cmd = vcf_annotator._query_aa_changes_command(
        input_vcf=str(annotated_vcf),
    )

    output = subprocess.check_output(cmd, shell=True).decode("utf-8")
    df = vcf_annotator._parse_to_aa_changes(output)

    data = df.to_csv(sep="\t", index=False)
    snapshot.snapshot_dir = snapshot_dir / "aa_changes"
    snapshot.assert_match(data, annotated_vcf.stem + ".tsv")


@pytest.mark.parametrize(
    "annotated_vcf",
    (vcf_test_data_dir / "annotated").glob("*.vcf"),
    ids=lambda p: p.stem,
)
def test_qc(vcf_annotator: VariantAnnotator, annotated_vcf: Path, snapshot):
    cmd = vcf_annotator._query_qc_command(input_vcf=str(annotated_vcf))

    output = subprocess.check_output(cmd, shell=True).decode("utf-8")
    df = vcf_annotator._parse_to_qc(output)

    data = df.to_csv(sep="\t", index=False)
    snapshot.snapshot_dir = snapshot_dir / "qc"
    snapshot.assert_match(data, annotated_vcf.stem + ".tsv")


@pytest.mark.parametrize(
    "annotated_vcf",
    (vcf_test_data_dir / "annotated").glob("*.vcf"),
    ids=lambda p: p.stem,
)
def test_summarize_aa_changes(
    vcf_annotator: VariantAnnotator, annotated_vcf: Path, snapshot
):
    df = vcf_annotator.summarize_aa_changes(input_vcf=str(annotated_vcf))

    data = df.to_csv(sep="\t", index=False)
    snapshot.snapshot_dir = snapshot_dir / "summary"
    snapshot.assert_match(data, annotated_vcf.stem + ".tsv")


@pytest.mark.parametrize(
    "annotated_vcf",
    (vcf_test_data_dir / "annotated").glob("*.vcf"),
    ids=lambda p: p.stem,
)
def test_summarize_nt_changes(
    vcf_annotator: VariantAnnotator, annotated_vcf: Path, snapshot
):
    df = vcf_annotator.summarize_nt_changes(input_vcf=str(annotated_vcf))

    data = df.to_csv(sep="\t", index=False)
    snapshot.snapshot_dir = snapshot_dir / "nt_summary"
    snapshot.assert_match(data, annotated_vcf.stem + ".tsv")
