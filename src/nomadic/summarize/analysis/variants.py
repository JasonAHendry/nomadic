import shutil
import subprocess
from collections.abc import Iterable
from logging import Logger
from pathlib import Path
from typing import Optional

import pandas as pd
from statsmodels.stats.proportion import proportion_confint

from nomadic.download.references import REFERENCE_COLLECTION
from nomadic.util.timer import Timer
from nomadic.util.vcf import (
    AA_CALL_COL,
    AA_CHANGE_COL,
    AA_POS_COL,
    GENE_COL,
    VariantAnnotator,
)
from nomadic.util.wrappers import bcftools


def load_variants_from_vcfs(
    expt_dirs: list[Path],
    *,
    caller: str,
    temp_dir: Path,
    filtered_vcf: Path,
    annotated_vcf: Path,
    bed_path: Path,
    reference_name: str,
    exclude_amplicons: Optional[list[str]] = None,
    exclude_mutations: Optional[list[str]] = None,
    log: Optional[Logger] = None,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Load variants directly from VCF files, rather than the summary CSVs

    This is more work, but allows us to get all variants that were called, even those that didn't pass filtering and thus don't appear in the summary CSVs.
    """

    timer = Timer("Variant loading from VCFs")
    timer.start()

    if log is not None:
        log.info(f"  Loading variants from VCFs in {len(expt_dirs)} experiments...")

    if any(SAMPLE_SEPERATOR in expt_dir.name for expt_dir in expt_dirs):
        raise ValueError(
            f"Experiment directories can not contain the string '{SAMPLE_SEPERATOR}', as this is used to separate experiment name and barcode when loading from VCFs. Please rename the following directories: {', '.join([d.name for d in expt_dirs if SAMPLE_SEPERATOR in d.name])}."
        )

    if temp_dir.exists():
        # Should be removed by initialization of the output directory, but just in case, we check here
        raise FileExistsError(
            f"Temporary directory {temp_dir} already exists. Please remove it before running this function."
        )
    temp_dir.mkdir()

    timer.time("setup")

    experiment_sample_mapping, temp_vcfs = load_and_reheader_vcfs(
        expt_dirs, output_dir=temp_dir
    )

    timer.time("Loading and reheadering VCF files")
    if log is not None:
        log.info(f"  Merging {len(temp_vcfs)} VCFs...")

    # Now we can merge all the temp VCFs together
    merge_and_filter_vcfs(temp_vcfs, output_path=filtered_vcf)
    timer.time("Merging VCF files")

    if log is not None:
        log.info("  Annotating variants...")

    REFERENCE_COLLECTION[reference_name].confirm_downloaded()
    annotator = VariantAnnotator(
        fasta_path=REFERENCE_COLLECTION[reference_name].fasta_path,
        gff_path=REFERENCE_COLLECTION[reference_name].gff_path,
        bed_path=(str(bed_path)),
        caller=caller,
    )

    annotator.annotate_variants(
        input_vcf=str(filtered_vcf), output_vcf=str(annotated_vcf)
    )

    timer.time("Annotating variants")
    if log is not None:
        log.info("  Summarizing amino acid changes...")

    variant_df = annotator.summarize_aa_changes(
        input_vcf=str(annotated_vcf),
        exclude_amplicons=exclude_amplicons,
        exclude_mutations=exclude_mutations,
    )

    timer.time("Summarizing amino acid changes")

    if log is not None:
        log.info("  Summarizing nt changes...")
    nt_df = annotator.summarize_nt_changes(
        input_vcf=str(annotated_vcf),
        exclude_amplicons=exclude_amplicons,
    )

    timer.time("Summarizing nucleotide changes")

    restore_barcodes(variant_df)
    restore_barcodes(nt_df)

    if not variant_df.empty:
        # Sanity checks that all worked and we have the same samples as before
        # Bcftools has some magic around sample names, so good to check they are consistent with what we expected
        all_samples_after = (
            variant_df.groupby(["expt_name"])["barcode"].unique().to_dict()
        )
        assert set(experiment_sample_mapping.keys()) == set(all_samples_after.keys()), (
            "Mismatch in experiments after loading from VCFs."
        )
        assert all(
            samples == set(all_samples_after[expt_dir])
            for expt_dir, samples in experiment_sample_mapping.items()
        ), "Mismatch in samples after loading from VCFs."

    timer.time("fixing sample names and sanity checking")

    if log is not None:
        log.info("  Done loading variants from VCFs.")
    timer.report()

    shutil.rmtree(temp_dir)

    return variant_df, nt_df


SAMPLE_SEPERATOR = "___"


def encode_barcodes(barcodes: list, *, expt_name: str, vcf_file: str) -> list:
    """Make sample names unique by concating expt name and experiment sample name (the barcode)
        Encode whitespace in sample names to avoid issues with bcftools, which treats whitespace in sample names specially.

    Args:
        barcodes: List of experiment sample names (barcodes).
        expt_name: Name of the experiment.
        vcf_file: Path to the VCF file being processed.

    Returns:
        List of encoded sample names suitable for use with bcftools.
    """
    # Last checks, should be handled already
    assert SAMPLE_SEPERATOR not in expt_name, (
        f"Experiment name {expt_name} cannot contain the string '{SAMPLE_SEPERATOR}'"
    )
    assert all(SAMPLE_SEPERATOR not in s for s in barcodes), (
        f"Samples in VCF file {vcf_file} cannot contain the string '{SAMPLE_SEPERATOR}'"
    )
    encoded_barcodes = [f"{expt_name}{SAMPLE_SEPERATOR}{s}" for s in barcodes]
    encoded_barcodes = [
        s.replace(" ", r"\ ") for s in encoded_barcodes
    ]  # replace space, bcftools treat spaces in samples names special
    return encoded_barcodes


def restore_barcodes(df: pd.DataFrame):
    """Restore original barcodes from encoded sample names.

    Args:
        df: DataFrame containing encoded sample names in the 'barcode' column.

    Returns:
        None. The DataFrame is modified in place to include 'expt_name' and restore the original 'barcode'.
    """
    sample_names = df["barcode"].str.replace(
        r"\ ",
        " ",
        regex=False,
    )  # reverse escaping of spaces

    sample_parts = sample_names.str.split(
        SAMPLE_SEPERATOR,
        n=1,
        expand=True,
        regex=False,
    )

    df.insert(0, "expt_name", sample_parts[0])
    df["barcode"] = sample_parts[1]


def filter_to_analysis_set(
    variant_df: pd.DataFrame,
    *,
    coverage_df: pd.DataFrame,
) -> pd.DataFrame:
    """
    Filter the variant dataframe to only include samples that passed quality control.
    """
    # Merge with the quality control results, then we can subset to the analysis set
    variant_df = pd.merge(
        right=variant_df,
        left=coverage_df.rename({"name": "amplicon"}, axis=1)[
            ["sample_id", "expt_name", "barcode", "chrom", "amplicon", "status"]
        ],
        on=["expt_name", "barcode", "chrom", "amplicon"],
    )

    return variant_df.query("status == 'pass'").drop(columns=["status"])


def filter_variants(vcf: Path, output_path: Path):
    """Filter a VCF file to only include PASS SNPs with at least 2 alleles."""
    subprocess.run(
        [
            "bcftools",
            "view",
            "--apply-filters",
            "PASS",
            "--types",
            "snps",
            "--min-alleles",
            "2",
            "-Oz",
            "-o",
            str(output_path),
            str(vcf),
        ],
        check=True,
    )


def merge_vcfs(vcfs: Iterable[Path], *, threads: int = 8, output_path: Path):
    subprocess.run(
        [
            "bcftools",
            "merge",
            "-Fx",
            "-Oz1",
            "--threads",
            str(threads),
            "--force-single",
            "--write-index",
            "-o",
            str(output_path),
        ]
        + [str(v) for v in vcfs],
        check=True,
    )


def merge_and_filter_vcfs(vcfs: Iterable[Path], *, output_path: Path, threads: int = 8):
    """Merging vcf files and then filtering

    It is much faster to not write the intermediate file if we don't need it.
    """
    merge_cmd = ["bcftools", "merge", "-Fx", "-Ou", "--force-single", *map(str, vcfs)]

    filter_cmd = [
        "bcftools",
        "view",
        "--apply-filters",
        "PASS",
        "--types",
        "snps",
        "--min-alleles",
        "2",
        "-Oz",
        "--threads",
        str(threads),
        "--write-index",
        "-o",
        str(output_path),
        "-",
    ]

    merge = subprocess.Popen(
        merge_cmd,
        stdout=subprocess.PIPE,
    )

    assert merge.stdout is not None

    filtered = subprocess.Popen(
        filter_cmd,
        stdin=merge.stdout,
    )

    merge.stdout.close()

    filter_returncode = filtered.wait()
    merge_returncode = merge.wait()

    if filter_returncode != 0:
        raise subprocess.CalledProcessError(
            filter_returncode,
            filter_cmd,
        )

    if merge_returncode != 0:
        raise subprocess.CalledProcessError(
            merge_returncode,
            merge_cmd,
        )


def load_and_reheader_vcfs(
    expt_dirs: Iterable[Path], *, output_dir: Path
) -> tuple[dict[str, set[str]], list[Path]]:
    """
    Load VCF files from a list of experiment directories, reheader them to make sample names unique, and return a mapping of experiment names to sample names.

    Args:
        expt_dirs: Iterable of paths to experiment directories containing VCF files.
        output_dir: Directory where reheadered VCF files will be written.

    Returns:
        A tuple containing:
            - A dictionary mapping experiment names to sets of sample names.
            - A list of paths to the reheadered VCF files.
    """

    # Record all samples for for sanity check after
    # expt_name -> set of samples in that experiment
    experiment_sample_mapping: dict[str, set[str]] = dict()

    vcfs: list[Path] = []

    for expt_dir in expt_dirs:
        expt_name = expt_dir.name
        vcf_dir = expt_dir / "vcfs"
        vcf_file = vcf_dir / "summary.variants.vcf.gz"
        if not vcf_file.exists():
            raise FileNotFoundError(
                f"VCF file not found at {vcf_file}. Cannot load variants from VCFs."
            )

        experiment_samples = (
            subprocess.check_output(["bcftools", "query", "-l", str(vcf_file)])
            .decode()
            .splitlines()
        )
        experiment_sample_set = set(experiment_samples)
        if len(experiment_sample_set) != len(experiment_samples):
            raise ValueError(
                f"Duplicate sample names found in VCF file {vcf_file}. Sample names must be unique to load from VCFs."
            )
        experiment_sample_mapping[expt_name] = experiment_sample_set
        # Make sample names unique by concating expt name and experiment sample name (the barcode)
        unique_samples = encode_barcodes(
            experiment_samples,
            expt_name=expt_name,
            vcf_file=str(vcf_file),
        )

        ### Reheader vcf file and move
        temp_vcf = output_dir / f"{expt_name}.variants.temp.vcf.gz"
        subprocess.run(
            [
                "bcftools",
                "reheader",
                "-n",
                ",".join(unique_samples),
                "-o",
                str(temp_vcf),
                str(vcf_file),
            ],
            check=True,
        )
        bcftools.index(temp_vcf)
        vcfs.append(temp_vcf)

    return experiment_sample_mapping, vcfs


# Columns used to group amino acid variants for aggregation and filtering.
AA_GROUP_COLUMNS = [
    "chrom",
    GENE_COL,
    AA_POS_COL,
    AA_CHANGE_COL,
]


def remove_false_positives(
    variants_df: pd.DataFrame, min_obs: int = 1, min_aa_wsaf: float = 0.15
):
    """Filter out likely false positive variant calls.

    Remove variants that have only been observed `min_obs` times across all samples in
    the analysis set, and with a WSAF below `min_wsaf`

    This removes likely false-positive mutations, which are usually rare and at low
    WSAF.

    """

    mut = variants_df.loc[variants_df[AA_CALL_COL].isin(["mixed", "mutant"])]
    df = variants_df.merge(
        mut.groupby(AA_GROUP_COLUMNS).agg(
            n_mut=pd.NamedAgg(AA_CALL_COL, len),
            aa_wsaf_max=pd.NamedAgg("aa_wsaf", "max"),
        ),
        on=AA_GROUP_COLUMNS,
        how="left",
    )
    df = df[~(df["n_mut"].le(min_obs) & df["aa_wsaf_max"].lt(min_aa_wsaf))].drop(
        columns=["n_mut", "aa_wsaf_max"]
    )
    return df


def remove_never_observed_variants(variants_df: pd.DataFrame) -> pd.DataFrame:
    """
    Remove variants that have never been observed in any sample in the analysis set.
    """

    mut = variants_df.loc[variants_df[AA_CALL_COL].isin(["mixed", "mutant"])]
    df = variants_df.merge(
        mut.groupby(AA_GROUP_COLUMNS).agg(
            n_mut=pd.NamedAgg(AA_CALL_COL, len),
        ),
        on=AA_GROUP_COLUMNS,
        how="left",
    )
    df = df[df["n_mut"].gt(0)].drop(columns=["n_mut"])
    return df


def compute_variant_prevalence(
    variants_df: pd.DataFrame,
    master_df: Optional[pd.DataFrame] = None,
    additional_groups: Optional[list[str]] = None,
) -> pd.DataFrame:
    """
    Compute the prevalence of each mutation in `variants_df`
    """
    if additional_groups is None:
        additional_groups = []

    if additional_groups:
        assert master_df is not None, (
            "master_df must be provided if additional_groups are used"
        )
        assert all(group in master_df.columns for group in additional_groups), (
            "all additional_groups must be columns in master_df"
        )
        variants_df = variants_df.merge(
            master_df[["sample_id", *additional_groups]],
            on="sample_id",
            how="left",
            validate="m:1",
        )

    passed_types = {"mixed", "mutant", "absent", "wt"}

    # Precompute so we can use fast sum aggregation
    variants_df = variants_df.assign(
        _passed=variants_df[AA_CALL_COL].isin(passed_types),
        _wt=variants_df[AA_CALL_COL].eq("wt"),
        _mixed=variants_df[AA_CALL_COL].eq("mixed"),
        _mut=variants_df[AA_CALL_COL].eq("mutant"),
    )

    prev_df = variants_df.groupby(
        AA_GROUP_COLUMNS + additional_groups, as_index=False
    ).agg(
        n_samples=(AA_CALL_COL, "size"),
        n_passed=("_passed", "sum"),
        n_wt=("_wt", "sum"),
        n_mixed=("_mixed", "sum"),
        n_mut=("_mut", "sum"),
    )
    has_passing_samples = prev_df["n_passed"].ne(0)
    # Compute frequencies
    prev_df.loc[has_passing_samples, "per_wt"] = (
        100
        * prev_df.loc[has_passing_samples, "n_wt"]
        / prev_df.loc[has_passing_samples, "n_passed"]
    )
    prev_df.loc[has_passing_samples, "per_mixed"] = (
        100
        * prev_df.loc[has_passing_samples, "n_mixed"]
        / prev_df.loc[has_passing_samples, "n_passed"]
    )
    prev_df.loc[has_passing_samples, "per_mut"] = (
        100
        * prev_df.loc[has_passing_samples, "n_mut"]
        / prev_df.loc[has_passing_samples, "n_passed"]
    )

    # Compute prevalence
    prev_df.loc[has_passing_samples, "prevalence"] = (
        prev_df.loc[has_passing_samples, "per_mixed"]
        + prev_df.loc[has_passing_samples, "per_mut"]
    )

    # Compute prevalence 95% confidence intervals
    low, high = proportion_confint(
        prev_df.loc[has_passing_samples, "n_mut"]
        + prev_df.loc[has_passing_samples, "n_mixed"],
        prev_df.loc[has_passing_samples, "n_passed"],
        alpha=0.05,
        method="beta",
    )
    prev_df.loc[has_passing_samples, "prevalence_lowci"] = 100 * low
    prev_df.loc[has_passing_samples, "prevalence_highci"] = 100 * high

    return prev_df


def rename_prevalence_by_cols(
    df: pd.DataFrame, old_prefix: str, new_prefix: str
) -> pd.DataFrame:
    """
    Rename prevalence by columns to remove the METADATA_COLUMN_PREFIX

    """
    prevalence_by_cols = [col for col in df.columns if col.startswith(old_prefix)]
    rename_dict = {
        col: col.replace(old_prefix, new_prefix, 1) for col in prevalence_by_cols
    }
    return df.rename(columns=rename_dict)


def remove_never_observed_nt_variants(variants_df: pd.DataFrame) -> pd.DataFrame:
    """
    Remove variants that have never been observed in any sample in the analysis set.
    """

    mut = variants_df.loc[variants_df["gt"].isin(["mixed", "mutant"])]
    df = variants_df.merge(
        mut.groupby(["chrom", "pos", "ref", "alt"]).agg(
            n_mut=pd.NamedAgg("gt", len),
        ),
        on=["chrom", "pos", "ref", "alt"],
        how="left",
    )
    df = df[df["n_mut"].gt(0)].drop(columns=["n_mut"])
    return df
