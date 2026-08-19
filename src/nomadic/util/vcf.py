import re
import shlex
import subprocess
from typing import Literal, Optional
from io import StringIO

import numpy as np
import pandas as pd

from line_profiler import profile


_AA_CHANGE_RE = re.compile(r"^(\d+)([A-Za-z])(?:>\1([A-Za-z\*]))?$")
_AA_POS_RE = re.compile(r"\d+")


class VariantAnnotator:
    AMP_HEADER = (
        "##INFO=<ID=AMP_ID,Number=1,Type=String,Description=Amplicon identifier>"
    )

    def __init__(
        self,
        fasta_path: str,
        gff_path: str,
        bed_path: str,
        caller: Literal["delve", "bcftools"],
    ):
        self.bed_path = bed_path
        self.fasta_path = fasta_path
        self.gff_path = gff_path

        if caller == "delve":
            self.wsaf_filltag_str = "FORMAT/WSAF=FORMAT/MVAF"
            self.wsaf_query_tag = "WSAF"
        elif caller == "bcftools":
            self.wsaf_filltag_str = "FORMAT/WSAF=1-FORMAT/AD/FORMAT/DP"
            self.wsaf_query_tag = "WSAF{0}"
        else:
            raise ValueError(f"Unknown caller: {caller}")

    def annotate_variants(self, input_vcf: str, output_vcf: str, *, threads: int = 8):

        cmd_annotate = self._annotate_command(input_vcf=input_vcf, pipe=True)
        cmd_fill_wsaf = self._fill_wsaf_command(pipe=True)
        cmd_csq = self._csq_command(output_vcf=output_vcf, threads=threads)

        cmd = f"{cmd_annotate} | {cmd_fill_wsaf} | {cmd_csq}"

        subprocess.run(cmd, check=True, shell=True)

    @profile
    def summarize_aa_changes(
        self,
        input_vcf: str = "-",
        exclude_amplicons: Optional[list[str]] = None,
        exclude_mutations: Optional[list[str]] = None,
    ) -> pd.DataFrame:
        """
        Query the amino acid changes from a VCF file and return them as a pandas DataFrame
        """
        output = subprocess.check_output(
            self._query_aa_changes_command(input_vcf=input_vcf), shell=True
        ).decode("utf-8")
        df_aa_changes = self._parse_to_aa_changes(
            output, exclude_amplicons=exclude_amplicons
        ).query("csq_type == 'missense'")
        df_aa_changes["aa_pos"] = df_aa_changes["aa_pos"].astype(
            int
        )  # all missense should have an amino acid position, so this should be safe

        output_qc = subprocess.check_output(
            self._query_qc_command(input_vcf=input_vcf), shell=True
        ).decode("utf-8")
        df_qc = self._parse_to_qc(output_qc, exclude_amplicons=exclude_amplicons)

        all_samples = df_qc["barcode"].unique()

        all_mutations = df_aa_changes[
            ["chrom", "amplicon", "gene", "aa_pos", "aa_change"]
        ].drop_duplicates()

        all_mutations["mutation"] = (
            all_mutations["gene"] + "-" + all_mutations["aa_change"]
        )

        if exclude_mutations is not None:
            all_mutations = all_mutations.loc[
                ~all_mutations["mutation"].isin(exclude_mutations)
            ]

        mutation_status_df = pd.merge(
            df_aa_changes[["barcode", "gene", "aa_pos", "aa_change", "nt_change"]],
            df_qc.query("mut_type in ['mixed', 'mutant']")[
                ["barcode", "gene", "aa_pos", "mut_type"]
            ],
            how="inner",
            on=["barcode", "gene", "aa_pos"],
            validate="many_to_one",
        )

        index = pd.DataFrame({"barcode": all_samples}).merge(all_mutations, how="cross")
        result_df = index.merge(
            df_qc.query("mut_type in ['wt', 'failed', 'unphased']")[
                ["barcode", "gene", "aa_pos", "mut_type"]
            ],
            how="left",
            on=["barcode", "gene", "aa_pos"],
            validate="many_to_one",
        )
        result_df["mut_type"] = (
            result_df["mut_type"]
            .astype("category")
            .cat.set_categories(
                ["mutant", "mixed", "absent", "wt", "unphased", "failed"], ordered=True
            )
        )
        result_df = result_df.merge(
            mutation_status_df,
            how="left",
            on=["barcode", "gene", "aa_pos", "aa_change"],
            suffixes=("", "_update"),
            validate="one_to_many",
        )
        result_df["mut_type"] = result_df["mut_type"].combine_first(
            result_df["mut_type_update"]
        )
        result_df.drop(columns=["mut_type_update"], inplace=True)
        result_df["mut_type"] = result_df["mut_type"].fillna("absent")

        # Remove any mutations for which there was never a mutant/mixed call
        result_df = result_df.loc[
            result_df.groupby(["gene", "aa_pos", "aa_change"], dropna=False)[
                "mut_type"
            ].transform(lambda x: (x.isin(["mutant", "mixed"])).any())
        ]

        # Merge in wsaf
        result_df = result_df.merge(
            df_qc[["barcode", "gene", "aa_pos", "mut_wsaf", "abs_wsaf", "depth"]],
            how="left",
            on=["barcode", "gene", "aa_pos"],
            validate="many_to_one",
        )

        result_df["aa_wsaf"] = pd.NA
        mut_mask = result_df["mut_type"].isin(["mutant", "mixed"])
        result_df.loc[mut_mask, "aa_wsaf"] = result_df.loc[mut_mask, "mut_wsaf"]

        abs_mask = result_df["mut_type"].isin(["wt", "absent"])
        result_df.loc[abs_mask, "aa_wsaf"] = result_df.loc[abs_mask, "abs_wsaf"]

        result_df.drop(columns=["mut_wsaf", "abs_wsaf"], inplace=True)

        result_df["n_changes"] = result_df["nt_change"].str.count(r"\+") + 1
        elligible_mask = result_df["mut_type"].eq("mixed") & result_df.query(
            "mut_type == 'mixed'"
        ).groupby(["barcode", "gene", "aa_pos"])["n_changes"].transform("size").gt(1)

        flip_mask = (
            result_df[elligible_mask]
            .groupby(["barcode", "gene", "aa_pos"])["n_changes"]
            .idxmin()
        )

        result_df.loc[flip_mask, "aa_wsaf"] = 1 - result_df.loc[flip_mask, "aa_wsaf"]
        result_df.drop(columns=["n_changes"], inplace=True)

        # Deduplicate samples with same AA change multiple times
        keys = ["barcode", "gene", "aa_change"]
        dup_mask = result_df.duplicated(keys, keep=False)

        result_df.loc[dup_mask, "nt_change"] = (
            result_df.loc[dup_mask].groupby(keys)["nt_change"].transform("|".join)
        )
        result_df.loc[dup_mask, "aa_wsaf"] = (
            result_df.loc[dup_mask].groupby(keys)["aa_wsaf"].transform("sum")
        )

        result_df = result_df.drop_duplicates(keys, keep="first")

        return result_df

    def summarize_nt_changes(
        self,
        input_vcf: str = "-",
        exclude_amplicons: Optional[list[str]] = None,
    ) -> pd.DataFrame:
        output = subprocess.check_output(
            self._query_nt_changes_command(input_vcf=input_vcf), shell=True
        ).decode("utf-8")
        df_nt_changes = self._parse_to_nt_changes(output)
        if exclude_amplicons is not None:
            df_nt_changes = df_nt_changes.loc[
                ~df_nt_changes["amplicon"].isin(exclude_amplicons)
            ]
        df_nt_changes["gt"] = mut_type_from_gt(df_nt_changes["gt"])

        df_nt_changes.loc[df_nt_changes["gt"] == "mutant", "tgt"] = df_nt_changes.loc[
            df_nt_changes["gt"] == "mutant", "tgt"
        ].replace(
            {
                "A/A": "A",
                "C/C": "C",
                "G/G": "G",
                "T/T": "T",
            }
        )

        df_mut = df_nt_changes.loc[
            df_nt_changes["gt"].isin(["mutant", "mixed"]),
            ["barcode", "chrom", "pos", "amplicon", "ref", "tgt", "gt", "wsaf"],
        ]
        df_mut["alt"] = df_mut["tgt"].str.split("/")
        df_mut = df_mut.explode("alt", ignore_index=True)
        df_mut = df_mut.loc[df_mut["ref"] != df_mut["alt"]]
        df_mut = df_mut.drop(columns=["tgt"])

        all_samples = df_nt_changes["barcode"].unique()

        all_mutations = df_mut[
            ["chrom", "pos", "amplicon", "ref", "alt"]
        ].drop_duplicates()

        result_df = pd.merge(
            pd.DataFrame({"barcode": all_samples}).merge(all_mutations, how="cross"),
            df_nt_changes[["barcode", "chrom", "pos", "dp"]],
            how="left",
            on=["barcode", "chrom", "pos"],
            validate="many_to_one",
        )

        result_df = result_df.merge(
            df_nt_changes.query("gt in ['wt', 'failed']")[
                ["barcode", "chrom", "pos", "gt"]
            ],
            how="left",
            on=["barcode", "chrom", "pos"],
            validate="many_to_one",
        )

        result_df = result_df.merge(
            df_nt_changes.query("gt in ['wt']")[["barcode", "chrom", "pos", "wsaf"]],
            how="left",
            on=["barcode", "chrom", "pos"],
            validate="many_to_one",
        )

        result_df["gt"] = (
            result_df["gt"]
            .astype("category")
            .cat.set_categories(
                ["mutant", "mixed", "absent", "wt", "failed"], ordered=True
            )
        )

        result_df = result_df.merge(
            df_mut[["barcode", "chrom", "pos", "ref", "alt", "gt", "wsaf"]],
            how="left",
            on=["barcode", "chrom", "pos", "ref", "alt"],
            suffixes=("", "_update"),
            validate="many_to_one",
        )
        result_df["gt"] = result_df["gt"].combine_first(result_df["gt_update"])
        result_df["wsaf"] = result_df["wsaf"].combine_first(result_df["wsaf_update"])
        result_df.drop(columns=["gt_update"], inplace=True)
        result_df.drop(columns=["wsaf_update"], inplace=True)

        result_df["gt"] = result_df["gt"].fillna("absent")
        result_df.loc[result_df["gt"] == "absent", "wsaf"] = 0

        return result_df

    def _annotate_command(
        self, input_vcf: str = "-", output_vcf: Optional[str] = None, pipe: bool = False
    ) -> str:
        """
        Create a string representing command required to annotate variants with
        their amplicon position
        """
        cmd = "bcftools annotate"
        cmd += f" -a {shlex.quote(self.bed_path)}"
        cmd += " -c CHROM,FROM,TO,AMP_ID"
        cmd += f" -H '{self.AMP_HEADER}'"
        if output_vcf:
            cmd += f" -Oz -o {shlex.quote(output_vcf)}"
        elif pipe:
            cmd += " -Ou"
        cmd += f" {shlex.quote(input_vcf)}"

        return cmd

    def _fill_wsaf_command(
        self, input_vcf: str = "-", output_vcf: Optional[str] = None, pipe: bool = False
    ) -> str:
        """
        Create a string representing command required to fill in the WSAF tag
        in the VCF file, which represents the within-sample allele frequency
        """
        cmd = "bcftools +fill-tags"
        if output_vcf:
            cmd += f" -Oz -o {shlex.quote(output_vcf)}"
        elif pipe:
            cmd += " -Ou"
        cmd += f" {shlex.quote(input_vcf)}"
        cmd += f" -- -t {self.wsaf_filltag_str}"

        return cmd

    def _csq_command(
        self,
        input_vcf: str = "-",
        output_vcf: Optional[str] = None,
        pipe: bool = False,
        threads: int = 8,
    ) -> str:
        """
        Create a string representing command required
        to compute variant consequences
        """
        cmd = "bcftools csq"
        cmd += f" -f {shlex.quote(self.fasta_path)}"
        cmd += f" -g {shlex.quote(self.gff_path)}"
        cmd += " --phase a"
        if output_vcf:
            cmd += " -Oz"
            cmd += f" --threads {threads}"
            cmd += f" -o {shlex.quote(output_vcf)}"
        elif pipe:
            cmd += " -Ou"
        cmd += f" {shlex.quote(input_vcf)}"

        return cmd

    def _query_aa_changes_command(
        self, input_vcf: str = "-", output_tsv: Optional[str] = None
    ) -> str:
        """
        Create a string representing command required to query the consequences
        from a VCF file and output them as a TSV file
        """
        cmd = "bcftools query"
        cmd += " -f '[%SAMPLE\t%CHROM\t%POS\t%AMP_ID\t%TBCSQ{*}\n]'"
        if output_tsv:
            cmd += f" -o {shlex.quote(output_tsv)}"
        cmd += f" {shlex.quote(input_vcf)}"

        return cmd

    def _parse_to_aa_changes(
        self, data: str, exclude_amplicons: Optional[list[str]] = None
    ) -> pd.DataFrame:
        """
        Parse the output of the query command into a pandas DataFrame
        """
        df = pd.read_csv(
            StringIO(data),
            sep="\t",
            names=["barcode", "chrom", "pos", "amplicon", "tbcsq"],
            dtype={
                "barcode": str,
                "chrom": str,
                "pos": int,
                "amplicon": str,
                "tbcsq": str,
            },
            na_values={"tbcsq": "."},
        )

        if exclude_amplicons is not None:
            df = df.loc[~df["amplicon"].isin(exclude_amplicons)]

        rows = []

        for row in df.itertuples(index=False):
            if pd.isna(row.tbcsq):
                continue

            for (
                csq_type,
                gene,
                transcript,
                coding_strand,
                aa_pos,
                aa_change,
                nt_change,
            ) in self._parse_tbcsq(row.tbcsq):
                rows.append(
                    (
                        row.barcode,
                        row.chrom,
                        row.pos,
                        row.amplicon,
                        csq_type,
                        gene,
                        transcript,
                        coding_strand,
                        aa_pos,
                        aa_change,
                        nt_change,
                    )
                )

        result = pd.DataFrame.from_records(
            rows,
            columns=[
                "barcode",
                "chrom",
                "pos",
                "amplicon",
                "csq_type",
                "gene",
                "transcript",
                "coding_strand",
                "aa_pos",
                "aa_change",
                "nt_change",
            ],
        )

        result["aa_pos"] = pd.array(result["aa_pos"], dtype="Int64")

        return result

    def _parse_tbcsq(self, value: str):
        for annotation in value.split(","):
            if annotation.startswith("@"):
                continue

            fields = annotation.split("|", 7)

            # Malformed annotation
            if len(fields) < 4:
                continue

            # If not in gene
            if len(fields) < 7:
                coding_strand = None
                aa_pos = None
                aa_change = None
                nt_change = None
            else:
                coding_strand = fields[4]
                aa_raw = fields[5]

                match = _AA_POS_RE.search(aa_raw)
                aa_pos = int(match.group()) if match else None

                match = _AA_CHANGE_RE.match(aa_raw)
                if match:
                    if match.group(3) is None:
                        aa_change = f"{match.group(2)}{match.group(1)}{match.group(2)}"
                    else:
                        aa_change = f"{match.group(2)}{match.group(1)}{match.group(3)}"
                else:
                    aa_change = None
                nt_change = fields[6]

            yield (
                fields[0],  # csq_type
                fields[1],  # gene
                fields[2],  # transcript
                coding_strand,
                aa_pos,
                aa_change,
                nt_change,
            )

    def _query_qc_command(
        self, input_vcf: str = "-", output_tsv: Optional[str] = None
    ) -> str:
        """
        Create a string representing command required to query the qc
        info from a VCF file and output them as a TSV file
        """
        cmd = "bcftools query"
        cmd += f" -f '[%SAMPLE\t%CHROM\t%POS\t%AMP_ID\t%INFO/BCSQ\t%GT\t%DP\t%{self.wsaf_query_tag}\n]'"
        if output_tsv:
            cmd += f" -o {shlex.quote(output_tsv)}"
        cmd += f" {shlex.quote(input_vcf)}"

        return cmd

    @profile
    def _parse_to_qc(
        self, data: str, exclude_amplicons: Optional[list[str]] = None
    ) -> pd.DataFrame:
        """
        Parse the output of the query command into a pandas DataFrame
        """

        df = pd.read_csv(
            StringIO(data),
            sep="\t",
            names=["barcode", "chrom", "pos", "amplicon", "bcsq", "gt", "dp", "wsaf"],
            na_values={"wsaf": ".", "dp": ".", "bcsq": "."},
            dtype={
                "barcode": str,
                "chrom": str,
                "pos": int,
                "amplicon": str,
                "bcsq": str,
                "gt": str,
                "dp": "Int64",
                "wsaf": "Float64",
            },
        )
        if exclude_amplicons is not None:
            df = df.loc[~df["amplicon"].isin(exclude_amplicons)]

        df["dp"] = df["dp"].fillna(0)

        bcsq_info = extract_bcsq_info(df["bcsq"])
        df[["gene", "aa_pos"]] = bcsq_info[["gene", "aa_pos"]]
        df[["gene", "aa_pos"]] = resolve_bcsq_references(df, bcsq_info["ref_pos"])

        df.dropna(subset=["gene", "aa_pos"], inplace=True)
        df["aa_pos"] = df["aa_pos"].astype(int)
        df["gene"] = df["gene"].astype(str)

        df["mut_type"] = mut_type_from_gt(df["gt"])

        df = df.assign(
            mut_wsaf=df["wsaf"].where(df["mut_type"].isin(["mutant", "mixed"])),
            abs_wsaf=df["wsaf"].where(df["mut_type"].isin(["wt"])),
        )

        return aggregate_qc_by_aa_pos(df)

    def _query_nt_changes_command(
        self, input_vcf: str = "-", output_tsv: Optional[str] = None
    ) -> str:
        """
        Create a string representing command required to query the nucleotide changes
        from a VCF file and output them as a TSV file
        """
        cmd = "bcftools query"
        cmd += " -f '[%SAMPLE\t%CHROM\t%POS\t%REF\t%AMP_ID\t%GT\t%TGT\t%DP\t%WSAF\n]'"
        if output_tsv:
            cmd += f" -o {shlex.quote(output_tsv)}"
        cmd += f" {shlex.quote(input_vcf)}"

        return cmd

    def _parse_to_nt_changes(self, data: str) -> pd.DataFrame:
        """
        Parse the output of the query command into a pandas DataFrame
        """

        df = pd.read_csv(
            StringIO(data),
            sep="\t",
            names=[
                "barcode",
                "chrom",
                "pos",
                "ref",
                "amplicon",
                "gt",
                "tgt",
                "dp",
                "wsaf",
            ],
            na_values={"wsaf": ".", "dp": "."},
            dtype={
                "barcode": str,
                "chrom": str,
                "pos": int,
                "ref": str,
                "amplicon": str,
                "gt": str,
                "tgt": str,
                "dp": "Int64",
                "wsaf": "Float64",
            },
        )

        return df


def extract_aa_pos(bcsq: str) -> Optional[int]:
    """
    Extract the amino acid position from a BCSQ string.
    Returns None if no amino acid change is present.
    """
    if bcsq == ".":
        return None

    bscqs = bcsq.split(",")

    aa_pos = []
    for s in bscqs:
        if s.startswith("@"):
            continue
        items = s.split("|")
        if len(items) < 6:
            return None
        aa_change = items[5]
        if aa_change == ".":
            return None
        m = re.search(r"(\d+)[A-Za-z]>\1+[A-Za-z]", aa_change)
        pos = m.group(1) if m else None
        if pos:
            aa_pos.append(int(pos))

    if not all(pos == aa_pos[0] for pos in aa_pos):
        raise ValueError(f"Inconsistent amino acid positions in BCSQ: {bcsq}")

    if not aa_pos:
        return None

    return aa_pos[0]


def extract_bcsq_info_vec(bcsq: pd.Series) -> pd.DataFrame:
    """
    Extract gene and amino acid position from a BCSQ string.
    Returns a DataFrame with columns: gene, aa_pos.

    Tries to do everything in a vectorized way to be fast.
    """
    annotations = bcsq.str.split(",").explode()
    annotations = annotations[~annotations.str.startswith("@", na=False)]

    parts = annotations.str.split("|", n=6, expand=True)

    parsed = pd.DataFrame(
        {
            "gene": parts[1],
            "aa_pos": parts[5].str.extract(r"(\d+)", expand=False).astype("Int64"),
        },
        index=annotations.index,
    )

    result = parsed.groupby(level=0).agg(
        n_genes=("gene", "nunique"),
        n_positions=("aa_pos", "nunique"),
        gene=("gene", "first"),
        aa_pos=("aa_pos", "first"),
    )

    if result["n_genes"].gt(1).any():
        idx = result["n_genes"].idxmax()
        raise ValueError(f"Inconsistent genes in BCSQ: {parsed.loc[idx]}")

    if result["n_positions"].gt(1).any():
        idx = result["n_positions"].idxmax()
        raise ValueError(
            f"Inconsistent amino acid positions in BCSQ: {parsed.loc[idx]}"
        )

    return result[["gene", "aa_pos"]].reindex(bcsq.index)


@profile
def resolve_bcsq_references(df: pd.DataFrame, ref_positions: pd.Series) -> pd.DataFrame:
    """
    Resolves gene, aa_pos for bcsq @position
    """
    result = df[["gene", "aa_pos"]].copy()

    # Filter to rows where gene and aa_pos are missing, but ref_positions is not empty
    mask = df["aa_pos"].isna() & df["gene"].isna() & ref_positions.notna()

    lookup = df.set_index(["barcode", "chrom", "pos"])[["aa_pos", "gene"]]

    if not lookup.index.is_unique:
        raise ValueError("The combination of chrom and pos is not unique in the VCF")

    refs = pd.DataFrame(
        {
            "barcode": df.loc[mask, "barcode"],
            "chrom": df.loc[mask, "chrom"],
            "ref_pos": ref_positions[mask],
        }
    )

    refs = refs.explode("ref_pos")

    # Create all keys to look up
    keys = pd.MultiIndex.from_arrays(
        [
            refs["barcode"],
            refs["chrom"],
            refs["ref_pos"],
        ],
        names=["barcode", "chrom", "pos"],
    )

    missing_keys = ~keys.isin(lookup.index)
    if missing_keys.any():
        raise ValueError(
            f"Some keys from the VCF are missing in the lookup table: {missing_keys}"
        )

    # Look them all up at once
    resolved = lookup.reindex(keys)
    # Set the index of the resolved DataFrame to match the original DataFrame so we can assign values back
    resolved.index = refs.index

    # TODO checks for consistency of gene and aa_pos for each barcode/chrom/pos combination

    grouped = resolved.groupby(level=0)
    resolved_unique = grouped.agg(
        gene=("gene", "first"),
        aa_pos=("aa_pos", "first"),
    )

    result.loc[resolved_unique.index, ["gene", "aa_pos"]] = resolved_unique[
        ["gene", "aa_pos"]
    ]
    return result


GT_TO_MUT_TYPE = {
    "./.": "failed",
    "0/0": "wt",
    **{
        f"{a}/{b}": "mutant" if a == b else "mixed"
        for a in range(4)
        for b in range(4)
        if (a, b) != (0, 0)
    },
}


def mut_type_from_gt(gts: pd.Series) -> pd.Series:
    result = gts.map(GT_TO_MUT_TYPE)

    unknown = result.isna()

    if unknown.any():
        raise ValueError(f"Unknown mutation type found in VCF: {gts[unknown].unique()}")

    return result.astype("category")


def aggregate_qc_by_aa_pos(df: pd.DataFrame) -> pd.DataFrame:
    # to be able to agg faster
    df = df.assign(
        failed=df["mut_type"].eq("failed"),
        mixed=df["mut_type"].eq("mixed"),
        mutant=df["mut_type"].eq("mutant"),
    )

    result = (
        df.groupby(["barcode", "gene", "aa_pos"])
        .agg(
            failed=("failed", "max"),
            mixed_count=("mixed", "sum"),
            mutant=("mutant", "max"),
            mut_wsaf=("mut_wsaf", "min"),
            abs_wsaf=("abs_wsaf", "max"),
            depth=("dp", "min"),
        )
        .reset_index()
    )

    result["mut_type"] = np.select(
        [
            result["failed"],
            result["mixed_count"] > 1,
            result["mixed_count"] == 1,
            result["mutant"],
        ],
        [
            "failed",
            "unphased",
            "mixed",
            "mutant",
        ],
        default="wt",
    )

    return result[
        [
            "barcode",
            "gene",
            "aa_pos",
            "mut_type",
            "mut_wsaf",
            "abs_wsaf",
            "depth",
        ]
    ]


_AA_POS_RE = re.compile(r"\d+")


@profile
def extract_bcsq_info(bcsq: pd.Series) -> pd.DataFrame:
    """
    Extract gene and amino acid position from BCSQ strings.

    Each BCSQ value may contain multiple comma-separated annotations.
    All non-reference annotations must agree on gene and amino acid position.

    Note: This python function is actually faster than the vectorized verion.
    """

    genes = []
    aa_positions = []
    ref_positions = []

    # Speed optimization, convert once to avoid pandas scalars
    values = bcsq.to_numpy(dtype=object, na_value=None)

    for value in values:
        if value is None:
            genes.append(None)
            aa_positions.append(None)
            ref_positions.append(None)
            continue

        gene = None
        aa_pos = None
        refs = None

        for annotation in value.split(","):
            if annotation.startswith("@"):
                pos = int(annotation[1:])
                if refs is None:
                    refs = [pos]
                else:
                    refs.append(pos)
                continue

            fields = annotation.split("|", 6)

            if len(fields) < 6:
                continue

            current_gene = fields[1]

            match = _AA_POS_RE.search(fields[5])
            current_aa_pos = int(match.group()) if match else None

            if current_gene:
                if gene is None:
                    gene = current_gene
                elif current_gene != gene:
                    raise ValueError(f"Inconsistent genes in BCSQ: {value}")

            if current_aa_pos is not None:
                if aa_pos is None:
                    aa_pos = current_aa_pos
                elif current_aa_pos != aa_pos:
                    raise ValueError(
                        f"Inconsistent amino acid positions in BCSQ: {value}"
                    )

        genes.append(gene)
        aa_positions.append(aa_pos)
        ref_positions.append(refs)

    return pd.DataFrame(
        {
            "gene": genes,
            "aa_pos": pd.array(aa_positions, dtype="Int64"),
            "ref_pos": ref_positions,
        },
        index=bcsq.index,
    )
