from pathlib import Path
import subprocess
from typing import Iterable, Optional

import pandas as pd
from statsmodels.stats.proportion import proportion_confint

from nomadic.download.references import REFERENCE_COLLECTION
from nomadic.util.timer import Timer
from nomadic.util.vcf import VariantAnnotator
from nomadic.util.wrappers import bcftools


# These columns are used to define unique variants
VARIANTS_GROUP_COLUMNS = [
    "gene",
    "chrom",
    "aa_pos",
]

# These groups are used to define a unique mutation, e.g. A127E
VARIANTS_MUTATION_COLUMNS = [
    "aa_change",
    "mut_type",
    "mutation",
]


def load_variants_from_vcfs(
    expt_dirs: Iterable[Path],
    *,
    caller: str,
    output_dir: Path,
    bed_path: Path,
    reference_name: str,
) -> pd.DataFrame:
    """
    Load variants directly from VCF files, rather than the summary CSVs

    This is more work, but allows us to get all variants that were called, even those that didn't pass filtering and thus don't appear in the summary CSVs.
    """

    timer = Timer("Variant loading from VCFs")
    timer.start()

    seperator = "___"
    if any(seperator in expt_dir.name for expt_dir in expt_dirs):
        raise ValueError(
            f"Experiment directories can not contain the string '{seperator}', as this is used to separate experiment name and barcode when loading from VCFs. Please rename the following directories: {', '.join([d.name for d in expt_dirs if seperator in d.name])}."
        )

    temp_dir = output_dir / "temp_vcf_processing"
    if temp_dir.exists():
        # Should be removed by initialization of the output directory, but just in case, we check here
        raise FileExistsError(
            f"Temporary directory {temp_dir} already exists. Please remove it before running this function."
        )
    temp_dir.mkdir(exist_ok=True)

    timer.time("setup")

    experiment_sample_mapping, temp_vcfs = load_and_reheader_vcfs(
        expt_dirs, seperator=seperator, output_dir=temp_dir
    )

    timer.time("Loading and reheadering VCF files")

    # Now we can merge all the temp VCFs together
    merged_vcf = output_dir / "summary.variants.vcf.gz"
    merge_vcfs(temp_vcfs, output_path=merged_vcf)

    timer.time("Merging VCF files")

    # Filtering
    filtered_vcf = output_dir / "summary.variants.filtered.vcf.gz"
    filter_variants(merged_vcf, filtered_vcf)

    timer.time("Filtering VCF file")

    REFERENCE_COLLECTION[reference_name].confirm_downloaded()
    annotator = VariantAnnotator(
        input_vcf=str(filtered_vcf),
        bed_path=(str(bed_path)),
        reference=REFERENCE_COLLECTION[reference_name],
        caller=caller,
        output_vcf=str(output_dir / "summary.variants.annotated.vcf.gz"),
    )
    annotator.run()

    timer.time("Annotating variants")

    # TODO saving and loading of this file should be removed
    annotator.convert_to_csv(f"{temp_dir}/summary.variants.merged.csv")

    timer.time("Converting VCF to CSV")

    variant_df = load_variant_summary_csv(f"{temp_dir}/summary.variants.merged.csv")

    timer.time("Loading variants into dataframe")

    # TODO the following code should probably be combined with the anotator
    variant_sample_names = variant_df["barcode"].str.replace(
        r"\ ", " "
    )  # reverse escaping of spaces
    variant_df.insert(0, "expt_name", variant_sample_names.str.split(seperator).str[0])
    variant_df["barcode"] = variant_sample_names.str.split(seperator).str[1]

    # Sanity checks that all worked and we have the same samples as before
    # Bcftools has some magic around sample names, so good to check they are consistent with what we expected
    all_samples_after = variant_df.groupby(["expt_name"])["barcode"].unique().to_dict()
    assert set(experiment_sample_mapping.keys()) == set(all_samples_after.keys()), (
        "Mismatch in experiments after loading from VCFs."
    )
    assert all(
        samples == set(all_samples_after[expt_dir])
        for expt_dir, samples in experiment_sample_mapping.items()
    ), "Mismatch in samples after loading from VCFs."

    timer.time("fixing sample names and sanity checking")

    variant_df = variant_df.query("aa_pos != -1")

    timer.time("Filtering out variants without aa_pos")

    # identifies an aa variant, note this does not include aa_change, as there might be different
    # changes for the same position
    groups = ["expt_name", "barcode", "chrom", "amplicon", "gene", "aa_pos"]
    # wt_df = (
    #     variant_df.groupby(groups)
    #     .agg(
    #         WT=pd.NamedAgg(column="gt_int", aggfunc=lambda x: (x == 0).all()),
    #         ref=("ref", "first"),
    #         alt=("alt", "first"),
    #         pos=("pos", "first"),
    #         qual=("qual", "first"),
    #         dp=("dp", "first"),
    #         wsaf=("wsaf", "first"),
    #         gt=("gt", "first"),
    #         gt_int=("gt_int", "first"),
    #     )
    #     .reset_index()
    #     .query("WT == True")
    # )
    # wt_df.drop(columns=["WT"], inplace=True)
    # wt_df["type"] = "wt"

    # fully vectorized way to get WT variants:
    # Wt is where gt_int is 0 for all variants in the group; if any variant has gt_int != 0, then it's not WT
    is_zero = variant_df["gt_int"].eq(0)
    wt_mask = is_zero.groupby([variant_df[c] for c in groups]).transform("all")
    # TODO, decide how to record gt, wsaf, dp here, currently just taking first
    wt_df = variant_df.loc[wt_mask].groupby(groups, as_index=False).first()
    wt_df["type"] = "wt"

    timer.time("Processing WT variants into final dataframe")

    # filtered_df = (
    #     variant_df.groupby(groups)
    #     .agg(
    #         filtered=pd.NamedAgg(column="gt_int", aggfunc=lambda x: (x == -1).any()),
    #         ref=("ref", "first"),
    #         alt=("alt", "first"),
    #         pos=("pos", "first"),
    #         qual=("qual", "first"),
    #         dp=("dp", "first"),
    #         wsaf=("wsaf", "first"),
    #         gt=("gt", "first"),
    #         gt_int=("gt_int", "first"),
    #     )
    #     .reset_index()
    #     .query("filtered == True")
    # )
    # filtered_df.drop(columns=["filtered"], inplace=True)
    # filtered_df["type"] = "filtered"

    # Fully vectorized way to get filtered variants: filtered is where gt_int is -1 for any variant in the group
    # TODO, decide if we actually want to record this at all
    is_filtered = variant_df["gt_int"].eq(-1)
    filtered_mask = is_filtered.groupby([variant_df[c] for c in groups]).transform(
        "any"
    )
    filtered_df = variant_df.loc[filtered_mask].groupby(groups, as_index=False).first()
    filtered_df["type"] = "filtered"

    timer.time("Processing filtered variants into final dataframe")

    # TODO, decide also here how to record gt, wsaf, dp. Currently just taking the one
    # that is in the row of the vcf where the aa_change was recorded
    mut_df = variant_df.dropna(subset=["aa_change"])
    mut_df["type"] = (
        mut_df.groupby(groups + ["aa_change"])["gt_int"]
        .transform(lambda x: (x == 1).any())
        .map({True: "mixed_mut", False: "mut"})
    )

    timer.time("Processing mutant variants into final dataframe")

    variant_df = pd.concat([wt_df, filtered_df, mut_df]).reset_index(drop=True)

    # drop duplicates if they exist. This should only happen when bcftools record a variant, but we want to filter it
    # drop so that type=filtered stays
    # NOTE this will have to be changed once we allow more than one mut per sample
    pref = variant_df["type"].eq("filtered")
    idx = pref.groupby([variant_df[c] for c in groups]).idxmax()
    variant_df = variant_df.loc[idx].reset_index(drop=True)

    # sanity checks
    assert variant_df.duplicated(subset=groups).sum() == 0, (
        "Duplicate variants found in final dataframe after processing. This should not happen"
    )

    timer.time("Concatenating final dataframe")

    timer.report()

    # shutil.rmtree(temp_dir)

    return variant_df


def load_variant_summary_csv(
    variants_csv: str, define_gene: bool = True
) -> pd.DataFrame:
    """
    Load an clean `summary.variants.csv` data produced by `nomadic`

    """

    # Settings
    NUMERIC_COLUMNS = ["dp", "wsaf"]
    UNPHASED_GT_TO_INT = {
        "./.": -1,
        "0/0": 0,
        "0/1": 1,
        "0/2": 1,
        "0/3": 1,
        "1/1": 2,
        "2/2": 2,
        "3/3": 2,
    }

    # Load
    variants_df = pd.read_csv(variants_csv)

    # Reformat numeric columns to be floats
    for c in NUMERIC_COLUMNS:
        variants_df[c] = [float(v) if v != "." else None for v in variants_df[c]]

    # Reformat unphased genotypes as integers
    variants_df.insert(
        variants_df.columns.tolist().index("gt") + 1,
        "gt_int",
        variants_df["gt"].map(UNPHASED_GT_TO_INT),
    )

    # Optionally reformat amplicon name to gene; assuming like  gene-...
    if define_gene:
        # --- gene ---
        gene_values = variants_df["amplicon"].str.split("-").str[0]

        variants_df.insert(
            variants_df.columns.get_loc("amplicon") + 1,
            "gene",
            gene_values.where(variants_df["amplicon"].notna()),
        )

        # --- mutation ---
        mutation_values = variants_df["gene"] + "-" + variants_df["aa_change"]

        variants_df.insert(
            variants_df.columns.get_loc("gene") + 1,
            "mutation",
            mutation_values.where(
                variants_df["gene"].notna() & variants_df["aa_change"].notna()
            ),
        )

    return variants_df


def load_and_concat_variants(expt_dirs: list[Path]) -> pd.DataFrame:
    """
    Load all of the variant calls for a set of experiment dirs

    Note that because we do note have the unfiltered VCF files, we have to do
    some additional work in order to ensure all mutations are represented across
    all experiments;
    """

    # Load data
    variant_dfs = []
    for expt_dir in expt_dirs:
        variant_csv = f"{expt_dir}/summary.variants.csv"
        variant_df = load_variant_summary_csv(variant_csv)
        variant_df.insert(0, "expt_name", expt_dir.name)
        variant_df.query("barcode != 'unclassified'", inplace=True)
        variant_dfs.append(variant_df)
    variant_df = pd.concat(variant_dfs)

    # Get all unique mutations
    MUT_COLUMNS = [
        "chrom",
        "pos",
        "ref",
        "alt",
        "strand",
        "aa_change",
        "aa_pos",
        "mut_type",
        "mutation",
        "amplicon",
        "gene",
    ]
    uniq_mutation_df = variant_df[MUT_COLUMNS].drop_duplicates()

    # Now we merge these back in for each barcode
    # - By doing a 'right' merge, we make sure all variants are present for each barcode
    # - The 'NaNs' that are present when the variant doesn't exist for that barcode
    #   get filled with zeros, i.e. we assume homozygous reference
    # - In reality, it could have been EITHER 0/0 or ./. (i.e. filtered), but we handle
    #   this afterwards when we merge with QC data;
    full_variant_dfs = []
    for (expt_name, barcode), bdf in variant_df.groupby(["expt_name", "barcode"]):
        mdf = pd.merge(bdf, uniq_mutation_df, on=MUT_COLUMNS, how="right")
        mdf["expt_name"] = expt_name
        mdf["barcode"] = barcode

        # Filled by default with hom reference
        mdf["gt"] = mdf["gt"].fillna("0/0")
        mdf["gt_int"] = mdf["gt_int"].fillna(0.0)
        mdf["wsaf"] = mdf["wsaf"].fillna(0.0)

        full_variant_dfs.append(mdf)

    return pd.concat(full_variant_dfs)


def filter_to_analysis_set(
    variant_df: pd.DataFrame,
    *,
    coverage_df: pd.DataFrame,
    excluded_amplicons: list[str],
    filtered_mutations: list[str],
) -> pd.DataFrame:
    # # Merge with the quality control results, then we can subset to the analysis set
    variant_df = pd.merge(
        left=coverage_df.rename({"name": "amplicon"}, axis=1)[
            ["expt_name", "barcode", "sample_id", "sample_type", "amplicon", "status"]
        ],
        right=variant_df,
        on=["expt_name", "barcode", "amplicon"],
    )

    return (
        variant_df.query("status == 'pass'")
        .query("amplicon not in @excluded_amplicons")
        .query("mutation not in @filtered_mutations")
    )


def filter_variants(vcf: Path, output_path: Path):
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


def merge_vcfs(vcfs: Iterable[Path], *, output_path: Path):
    subprocess.run(
        ["bcftools", "merge", "-Fx", "-Oz", "--force-single", "-o", str(output_path)]
        + [str(v) for v in vcfs],
        check=True,
    )

    bcftools.index(output_path)


def load_and_reheader_vcfs(
    expt_dirs: Iterable[Path], *, seperator: str, output_dir: Path
) -> tuple[dict[str, set[str]], list[Path]]:
    """
    Load and reheader VCFs from multiple experiments to have unique sample names.

    Load VCFs from multiple experiments, reheader them to have unique sample names,
    and return a mapping of experiment names to sample names and a list of paths to the reheadered VCFs."""
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
        # Last checks, should be handled already
        assert seperator not in expt_name, (
            f"Experiment name {expt_name} cannot contain the string '{seperator}'"
        )
        assert all(seperator not in s for s in experiment_samples), (
            f"Samples in VCF file {vcf_file} cannot contain the string '{seperator}'"
        )
        # Make sample names unique by concating expt name and experiment sample name (the barcode)
        unique_samples = [f"{expt_name}{seperator}{s}" for s in experiment_samples]
        unique_samples = [
            s.replace(" ", r"\ ") for s in unique_samples
        ]  # replace space, bcftools treat spaces in samples names special

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


def filter_false_positives(
    variants_df: pd.DataFrame, min_obs: int = 1, min_wsaf: float = 0.15
):
    """Filter out likely false positive variant calls.

    Remove variants that have only been observed `min_obs` times across all samples in
    the analysis set, and with a WSAF below `min_wsaf`

    This removes likely false-positive mutations, which are usually rare and at low
    WSAF.

    """

    mut = variants_df.loc[variants_df["type"].isin(["mixed_mut", "mut"])]
    df = variants_df.merge(
        mut.groupby(VARIANTS_GROUP_COLUMNS + VARIANTS_MUTATION_COLUMNS).agg(
            n_mut=pd.NamedAgg("type", len), wsaf_max=pd.NamedAgg("wsaf", "max")
        ),
        on=VARIANTS_GROUP_COLUMNS + VARIANTS_MUTATION_COLUMNS,
        how="left",
    )
    df = df[~(df["n_mut"].le(min_obs) & df["wsaf_max"].lt(min_wsaf))].drop(
        columns=["n_mut", "wsaf_max"]
    )
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

    agg_aa_change_df = (
        variants_df.loc[variants_df["type"].isin(["mixed_mut", "mut"])]
        .groupby(
            VARIANTS_GROUP_COLUMNS + VARIANTS_MUTATION_COLUMNS + additional_groups,
        )
        .agg(
            n_mixed=pd.NamedAgg("type", lambda x: (x == "mixed_mut").sum()),
            n_mut=pd.NamedAgg("type", lambda x: (x == "mut").sum()),
        )
    )

    groups = (
        variants_df[VARIANTS_GROUP_COLUMNS + additional_groups]
        .drop_duplicates()
        .dropna()
    )
    muts = (
        variants_df[VARIANTS_GROUP_COLUMNS + VARIANTS_MUTATION_COLUMNS]
        .query("mut_type == 'missense'")
        .drop_duplicates()
        .dropna()
    )

    # Build full index so we see also values for groups that have no mutation
    full_index = (
        groups.merge(muts, how="inner", on=VARIANTS_GROUP_COLUMNS)
        .set_index(
            VARIANTS_GROUP_COLUMNS + VARIANTS_MUTATION_COLUMNS + additional_groups
        )
        .index
    )
    # Ensure all n_mut, n_mixed are filled with zeros
    agg_aa_change_df = agg_aa_change_df.reindex(full_index).reset_index().fillna(0)

    agg_aa_pos_df = variants_df.groupby(
        VARIANTS_GROUP_COLUMNS + additional_groups,
        as_index=False,
    ).agg(
        n_samples=pd.NamedAgg("type", "size"),
        n_passed=pd.NamedAgg("type", lambda x: sum(x != "filtered")),
        n_wt=pd.NamedAgg("type", lambda x: sum(x == "wt")),
    )

    prev_df = agg_aa_change_df.merge(
        agg_aa_pos_df,
        on=VARIANTS_GROUP_COLUMNS + additional_groups,
        how="left",
        validate="m:1",
    )

    # Compute frequencies
    prev_df["per_wt"] = 100 * prev_df["n_wt"] / prev_df["n_passed"]
    prev_df["per_mixed"] = 100 * prev_df["n_mixed"] / prev_df["n_passed"]
    prev_df["per_mut"] = 100 * prev_df["n_mut"] / prev_df["n_passed"]

    # Compute prevalence
    prev_df["prevalence"] = prev_df["per_mixed"] + prev_df["per_mut"]

    # Compute prevalence 95% confidence intervals
    low, high = proportion_confint(
        prev_df["n_mut"] + prev_df["n_mixed"],
        prev_df["n_passed"],
        alpha=0.05,
        method="beta",
    )
    prev_df["prevalence_lowci"] = 100 * low
    prev_df["prevalence_highci"] = 100 * high

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
