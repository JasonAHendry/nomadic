import os
import re
import warnings
import subprocess
import shlex
import pandas as pd
from dataclasses import dataclass
from nomadic.download.references import Reference


# TODO:
# - Setup logging in nomadic
# log = logging.getLogger()


# ================================================================
# Parsing output strings from `bcftools csq`
#
# ================================================================


@dataclass
class Consequence:
    csq: str
    target_gene: str
    transcript: str
    biotype: str

    # Optional, assigned only if within coding sequence
    strand: str = None
    aa_change: str = None
    nt_change: str = None

    # Optional, assigned only if aa_change exists
    concise_aa_change: str = None
    aa_pos: int = -1  # too keep column int in pandas

    def __post_init__(self):
        """
        Clean `aa_change` and `aa_pos`

        """
        if self.aa_change is None:
            return

        AAs = "ARNDCEQGHILKMFPSTWYV"
        stop = r"\*"
        match = re.match(
            f"([0-9]+)([{AAs}|{stop}])(?:>[0-9]+([{AAs}|{stop}]))?", self.aa_change
        )

        if match is None:
            warnings.warn(f"Unable to parse AA change for: {self.aa_change}.")
            return

        aa_pos, from_aa, to_aa = match.groups()
        self.aa_pos = int(aa_pos)

        # Handle synonymous case
        if to_aa is None:
            to_aa = from_aa

        self.concise_aa_change = f"{from_aa}{aa_pos}{to_aa}"

    @classmethod
    def from_string(cls, csq_string: str):
        """
        Parse from the output string

        What to do if "double"?
        Or @
        Or *

        """
        if csq_string == ".":  # intergenic
            return cls(".", ".", ".", ".")

        consequences = csq_string.split(",")

        if any(c.startswith("@") for c in consequences):
            # This means a change for this aa was already recorded elsewhere
            # it might result in the same change, and we can not handle multiple same changes
            return cls(".", ".", ".", ".")
        if len(consequences) > 1:
            warnings.warn(
                f"Found multiple consequences of variant: {csq_string}! Keeping only first."
            )

        fields = consequences[0].split("|")
        assert len(fields) >= 4, f"Failed for {csq_string}"
        assert len(fields) <= 7, f"Failed for {csq_string}"

        return cls(*fields)


# ================================================================
# Annotating a VCF with information about amplicons & effects
#
# ================================================================


class VariantAnnotator:
    AMP_HEADER = (
        "##INFO=<ID=AMP_ID,Number=1,Type=String,Description=Amplicon identifier>"
    )

    def __init__(
        self,
        input_vcf: str,
        bed_path: str,
        reference: Reference,
        caller: str,
        output_vcf: str = None,
    ) -> None:
        # Inputs
        self.input_vcf = input_vcf
        self.bed_path = bed_path
        self.reference = reference  # The GFF path needs to be GFF3 compliant
        self.caller = caller

        # Output
        self.output_vcf = output_vcf
        if self.output_vcf is None:
            self.output_vcf = self.input_vcf.replace(".vcf.gz", ".annotated.vcf.gz")

    def _get_wsaf_command(self, input_vcf: str = "-", output_vcf: str = "") -> str:
        """
        Compute the WSAF for each variant based on allelic depths
        """
        cmd = "bcftools +fill-tags"
        if output_vcf:
            cmd += f" -Oz -o {shlex.quote(output_vcf)}"
        cmd += f" {shlex.quote(input_vcf)}"
        if self.caller == "delve":
            cmd += " -- -t FORMAT/WSAF=FORMAT/MVAF"
        elif self.caller == "bcftools":
            cmd += " -- -t FORMAT/WSAF=1-FORMAT/AD/FORMAT/DP"
        else:
            raise RuntimeError(f"Unknown caller: {self.caller}")

        return cmd

    def _get_annotate_command(self, input_vcf: str = "-", output_vcf: str = "") -> str:
        """
        Create a string representing command required to annotate variants with
        their amplicon position
        """
        cmd = "bcftools annotate"
        cmd += f" -a {shlex.quote(self.bed_path)}"
        cmd += " -c CHROM,FROM,TO,AMP_ID"
        cmd += f" -H '{self.AMP_HEADER}'"
        cmd += " -Oz"
        if output_vcf:
            cmd += f" -o {shlex.quote(output_vcf)}"
        cmd += f" {shlex.quote(input_vcf)}"

        return cmd

    def _get_csq_command(self, input_vcf: str = "-", output_vcf: str = "") -> str:
        """
        Create a string representing command required
        to compute variant consequences
        """
        cmd = "bcftools csq"
        cmd += f" -f {shlex.quote(self.reference.fasta_path)}"
        cmd += f" -g {shlex.quote(self.reference.gff_standard_path)}"
        cmd += " --phase a"
        cmd += " -Oz"
        if output_vcf:
            cmd += f" -o {shlex.quote(output_vcf)}"
        cmd += f" {shlex.quote(input_vcf)}"

        return cmd

    def run(self):
        """
        Annotate variant calls using `bcftools csq`
        """
        cmd_tags = self._get_wsaf_command(
            input_vcf=self.input_vcf,
        )
        cmd_annot = self._get_annotate_command(
            input_vcf="-",
        )
        cmd_csq = self._get_csq_command(input_vcf="-", output_vcf=self.output_vcf)
        cmd = f"{cmd_tags} | {cmd_annot} | {cmd_csq}"
        subprocess.run(cmd, shell=True, check=True)

    def _convert_to_tsv(self, output_tsv: str) -> None:
        """
        Convert from a VCF file into a text file with the indicated
        seperator

        """
        # Define fixed and per-sample fields
        # Note that sample name also gets added.
        fixed = {
            "chrom": "CHROM",
            "pos": "POS",
            "ref": "REF",
            "alt": "ALT",
            "qual": "QUAL",
            "consequence": "BCSQ",
            "amplicon": "AMP_ID",
        }

        if self.caller == "delve":
            called = {"gt": "GT", "dp": "DP", "wsaf": "WSAF"}
        elif self.caller == "bcftools":
            called = {"gt": "GT", "gq": "GQ", "dp": "DP", "wsaf": "WSAF{0}"}
        else:
            raise RuntimeError(f"Unknown caller: {self.caller}")

        # Write header
        sep = "\t"
        cmd_header = f"printf 'barcode{sep}{sep.join(list(fixed) + list(called))}\n' > {shlex.quote(output_tsv)}"
        sepb = f"{sep}%"  # for bcftools, % accesses variable

        # Iterate and query for each sample
        cmd_query = (
            f"for sample in `bcftools query {shlex.quote(self.output_vcf)} -l`; do"
        )
        cmd_query += "  bcftools query -s $sample"
        cmd_query += f' -f "$sample{sepb}{sepb.join(fixed.values())}{sep}[%{sepb.join(called.values())}]\n"'
        cmd_query += f" {shlex.quote(self.output_vcf)} >> {shlex.quote(output_tsv)};"
        cmd_query += " done;"

        cmd = f"{cmd_header} && {cmd_query}"
        subprocess.run(cmd, shell=True, check=True)

    def _parse_consequences(
        self, input_tsv: str, output_file: str, output_sep: str = "\t"
    ):
        """
        Parse the consequenc string in the TSV
        """
        df = pd.read_csv(input_tsv, sep="\t")
        csqs = [Consequence.from_string(c) for c in df["consequence"]]
        if csqs:
            mut_type, aa_change, aa_pos, strand = zip(
                *[(c.csq, c.concise_aa_change, c.aa_pos, c.strand) for c in csqs]
            )
        else:
            warnings.warn(f"No mutations passed quality control for {self.input_vcf}.")
            mut_type, aa_change, aa_pos, strand = None, None, None, None

        df.insert(6, "mut_type", mut_type)
        df.insert(7, "aa_change", aa_change)
        df.insert(8, "aa_pos", aa_pos)
        df.insert(9, "strand", strand)
        df.drop("consequence", axis=1, inplace=True)
        df.sort_values(["barcode", "chrom", "pos"], inplace=True)
        df.to_csv(output_file, sep=output_sep, index=False)

    def convert_to_tsv(self, output_tsv: str = None):
        """
        Convert the annotated VCF file to a TSV; and then
        improve formatting of the consequence field
        """

        if output_tsv is None:
            output_tsv = self.output_vcf.replace(".vcf.gz", ".tsv")
        self._convert_to_tsv(output_tsv)
        self._parse_consequences(input_tsv=output_tsv, output_file=output_tsv)

    def convert_to_csv(self, output_csv: str = None):
        """
        Convert from VCF to a CSV

        NB: We need to go via TSV because `bcftools csq` uses
        a comma as the delimiter when multiple consequences are
        identified for a single variant.

        """

        if output_csv is None:
            output_csv = self.output_vcf.replace(".vcf.gz", ".csv")

        temp_tsv = output_csv.replace(".csv", ".temp.tsv")
        # using tsv file to deal with commas in fields like consequence and multi-allelic alt
        self._convert_to_tsv(temp_tsv)
        self._parse_consequences(
            input_tsv=temp_tsv, output_file=output_csv, output_sep=","
        )
        os.remove(temp_tsv)
