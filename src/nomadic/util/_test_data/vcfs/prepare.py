import subprocess
from pathlib import Path

from nomadic.download.references import PlasmodiumFalciparum3D7
from nomadic.util.vcf import VariantAnnotator

test_vcf_folders = [Path("src/nomadic/util/_test_data/vcfs/delve")]

output_folder = Path("src/nomadic/util/_test_data/vcfs/annotated")


def main():
    if not output_folder.exists():
        output_folder.mkdir()
    for folder in test_vcf_folders:
        for vcf_file in folder.glob("*.vcf"):
            print(f"Preparing {vcf_file}")
            prepare_vcf(vcf_file)


def prepare_vcf(vcf_file: Path):

    reference = PlasmodiumFalciparum3D7()
    bed_path = Path("src/nomadic/util/_test_data/beds/nomadsMVP.amplicons.bed")
    annotator = VariantAnnotator(
        fasta_path=reference.fasta_path,
        gff_path=reference.gff_path,
        bed_path=str(bed_path),
        caller="delve",
    )

    cmd = f"{annotator._annotate_command(input_vcf=str(vcf_file))} | {annotator._fill_wsaf_command()} | {annotator._csq_command(output_vcf=str(output_folder / vcf_file.name))}"
    print(f"Running command: {cmd}")
    subprocess.run(cmd, shell=True, check=True)


if __name__ == "__main__":
    main()
