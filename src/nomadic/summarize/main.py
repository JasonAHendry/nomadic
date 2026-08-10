import glob
import os
import shutil
from typing import Iterable, Optional
from pathlib import Path
import subprocess
import time

import pandas as pd

from nomadic.download.references import REFERENCE_COLLECTION
from nomadic.summarize.analysis.input_data import common_caller
from nomadic.summarize.analysis.variants import (
    rename_prevalence_by_cols,
)
from nomadic.summarize.compute import (
    compute_variant_prevalence,
    filter_false_positives,
    gene_deletion_prevalence_by,
    gene_deletions,
)
from nomadic.summarize.dashboard.builders import BasicSummaryDashboard
from nomadic.summarize.analysis.inventory import (
    add_inventory_status,
    compute_throughput,
    create_inventory_df,
    drop_excluded_samples,
    experiments_in_inventory,
    n_excluded_samples,
    n_field_samples,
)
from nomadic.summarize.analysis.metadata import (
    METADATA_COLUMN_PREFIX,
    get_shared_metadata_columns,
    load_master_metadata,
    master_metadata_from_expts,
    normalize_metadata,
    validate_metadata,
)
from nomadic.summarize.analysis.qc import (
    add_quality_control_columns,
    add_quality_control_status_column,
    amplicons_qc_summary,
    compute_field_coverage_summary,
    create_region_coverage_df,
    experiment_qc_summary,
    replicates_amplicon_qc,
    replicates_qc,
    samples_qc,
)
from nomadic.util.panel import get_panel_settings
from nomadic.util.summary import looks_like_summary_dir
from nomadic.util.vcf import VariantAnnotator
from nomadic.util.workspace import Workspace
from nomadic.util.dirs import produce_dir
from nomadic.util.regions import RegionBEDParser, common_regions
from nomadic.util.experiment import (
    experiment_outputs,
)
from nomadic.util.logging_config import LoggingFascade
from nomadic.util.port import next_free_port
from nomadic.util.summary_settings import (
    Settings,
    load_settings,
)
from nomadic.util.wrappers import bcftools


# --------------------------------------------------------------------------------
# Variant analysis
#
# --------------------------------------------------------------------------------


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
        # variants_df.insert(
        #     variants_df.columns.get_loc("amplicon") + 1,
        #     "gene",
        #     [a.split("-")[0] for a in variants_df["amplicon"]],
        # )
        # variants_df.insert(
        #     variants_df.columns.get_loc("gene") + 1,
        #     "mutation",
        #     [
        #         f"{gene}-{aa_change}"
        #         for gene, aa_change in zip(
        #             variants_df["gene"], variants_df["aa_change"]
        #         )
        #     ],
        # )
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


def load_and_concat_variants(expt_dirs: list[str]) -> pd.DataFrame:
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
        variant_df.insert(0, "expt_name", os.path.basename(expt_dir))
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


def load_variants_from_vcfs(
    expt_dirs: Iterable[str],
    *,
    caller: str,
    output_dir: Path,
    summary_regions: RegionBEDParser,
    reference_name: str,
) -> pd.DataFrame:
    """
    Load variants directly from VCF files, rather than the summary CSVs

    This is more work, but allows us to get all variants that were called, even those that didn't pass filtering and thus don't appear in the summary CSVs.
    """

    timer = Timer()
    timer.start()
    seperator = "___"

    if any(seperator in expt_dir for expt_dir in expt_dirs):
        raise ValueError(
            f"Experiment directories can not contain the string '{seperator}', as this is used to separate experiment name and barcode when loading from VCFs. Please rename the following directories: {', '.join([d for d in expt_dirs if seperator in d])}."
        )

    temp_dir = output_dir / "temp_vcf_processing"
    temp_dir.mkdir(exist_ok=True)

    # Record all samples for for sanity check after
    # expt_name -> set of samples in that experiment
    experiment_sample_mapping: dict[str, set[str]] = dict()

    timer.time("setup")

    temp_vcfs = []

    for expt_dir in expt_dirs:
        expt_name = os.path.basename(expt_dir)
        vcf_dir = Path(expt_dir) / "vcfs"
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
        unique_samples = [f"{expt_name}{seperator}{s}" for s in experiment_samples]
        unique_samples = [
            s.replace(" ", r"\ ") for s in unique_samples
        ]  # replace space, bcftools treat spaces in samples names special

        ### Reheader vcf file and move
        temp_vcf = temp_dir / f"{Path(expt_dir).name}.variants.temp.vcf.gz"
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
        temp_vcfs.append(temp_vcf)

    timer.time("Loading and reheadering VCF files")

    # Now we can merge all the temp VCFs together
    merged_vcf = output_dir / "summary.variants.vcf.gz"

    subprocess.run(
        ["bcftools", "merge", "-Fx", "-Oz", "--force-single", "-o", str(merged_vcf)]
        + [str(v) for v in temp_vcfs],
        check=True,
    )

    bcftools.index(merged_vcf)
    timer.time("Merging VCF files")

    # Filtering
    filtered_vcf = output_dir / "summary.variants.filtered.vcf.gz"
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
            str(filtered_vcf),
            str(merged_vcf),
        ],
        check=True,
    )

    timer.time("Filtering VCF file")

    REFERENCE_COLLECTION[reference_name].confirm_downloaded()
    annotator = VariantAnnotator(
        input_vcf=str(filtered_vcf),
        bed_path=summary_regions.path,
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

    shutil.rmtree(temp_dir)

    return variant_df


# --------------------------------------------------------------------------------
# Main
#
# --------------------------------------------------------------------------------
def main(
    *,
    workspace: Workspace,
    output_dir: Path,
    expt_dirs: list[str],
    summary_name: str,
    metadata_path: Optional[Path],
    settings_file_path: Path,
    maps: list[str],
    show_dashboard: bool = True,
    prevalence_by: list[str],
    no_master_metadata: bool = False,
    qc_min_coverage: int,
    qc_max_contam: float,
    qc_replicate_passing_threshold: float,
    host: str = "127.0.0.1",
    port: Optional[int] = None,
) -> None:
    """
    Define the main function for the summary analysis
    """
    assert (metadata_path is not None) or no_master_metadata

    ###############
    # Log setup
    ##############
    log = LoggingFascade(logger_name="nomadic")
    log.info("Input parameters:")
    log.info(f"  Summary Name: {summary_name}")
    if no_master_metadata:
        log.info("  No master metadata will be used.")
    else:
        log.info(f"  Master metadata: {metadata_path}")
    log.info(f"  Setting file: {settings_file_path}")
    log.info(f"  Found {len(expt_dirs)} experiment directories.")
    log.info(f"  Output directory: {output_dir}")

    if output_dir.exists():
        if looks_like_summary_dir(output_dir):
            shutil.rmtree(output_dir)
        else:
            raise ValueError(
                f"Output directory {output_dir} already exists and does not look like a summary directory. Please remove it or choose a different output directory."
            )
    produce_dir(str(output_dir))

    # Check experiments are complete
    if not expt_dirs:
        raise ValueError("No experiment directories found to summarize.")

    log.info("Data status:")
    expts = [
        experiment_outputs(expt_dir, allow_missing_files=["depth", "fastq"])
        for expt_dir in expt_dirs
    ]
    log.info(f"  All {len(expts)} experiments are complete.")

    # Check experiments are consistent
    regions = common_regions([expt.regions for expt in expts])
    log.info("  All experiments use the same regions.")
    caller = common_caller([expt.caller for expt in expts])
    log.info(f"  All experiments use same variant caller: {caller}")
    reference_name = common_reference_name([expt.reference_name for expt in expts])
    log.info(f"  All experiments use same reference genome: {reference_name}")

    if expts:
        panel_name = regions.name
    else:
        panel_name = "Unknown"
    log.info(f"  Panel used: {panel_name}")
    panel_settings = get_panel_settings(panel_name)
    log.info(f"Loaded panel settings for panel '{panel_settings.name}'.")

    settings: Settings = Settings()
    if settings_file_path.exists():
        settings = load_settings(settings_file_path)
        log.info(f"  Loaded summary settings from {settings_file_path}.")

    ####################
    # Load metadata and check for consistency
    ####################
    full_inventory_df = create_inventory_df(expts)
    if metadata_path is not None and not no_master_metadata:
        master_metadata_df = load_master_metadata(metadata_path, settings=settings)
    else:
        # create metadata from experiment meta data files
        shared_columns = get_shared_metadata_columns([expt.metadata for expt in expts])
        log.info(
            f"  Found {len(shared_columns)} non-essential shared columns across all experiment metadata files: {', '.join(shared_columns)}"
        )
        master_metadata_df = master_metadata_from_expts(
            expts, shared_columns=shared_columns
        )
        log.info(
            f"  Created master metadata from experiment metadata files with {len(master_metadata_df)} samples."
        )

    master_metadata_df = master_metadata_df.pipe(normalize_metadata).pipe(
        validate_metadata
    )
    full_inventory_df = add_inventory_status(full_inventory_df, master_metadata_df)

    n_excluded = n_excluded_samples(full_inventory_df)
    full_inventory_df.to_csv(f"{output_dir}/inventory.csv", index=False)

    inventory_df = drop_excluded_samples(full_inventory_df)

    # Filter experiment directories to only those that are in the master metadata, i.e. that have at least one included sample
    expt_dirs = experiments_in_inventory(inventory_df, expt_dirs)
    if n_field_samples(inventory_df) == 0:
        log.info("No included field samples, exiting...")
        return

    # Throughput data
    log.info("Overall sequencing throughput:")
    throughput, throughput_df = compute_throughput(inventory_df)
    log.info(f"  Positive controls: {throughput.n_pos}")
    log.info(f"  Negative controls: {throughput.n_neg}")
    log.info(f"  Fields samples sequenced (total): {throughput.n_field_total}")
    log.info(f"  Field samples (unique): {throughput.n_field_unique}")
    log.info(f"  Excluded samples: {n_excluded}")
    throughput_df.to_csv(f"{output_dir}/summary.throughput.csv", index=True)

    ############################
    # Quality control
    ############################

    # Now let's evaluate coverage
    coverage_df = (
        create_region_coverage_df(expt_dirs, inventory_df)
        .pipe(
            add_quality_control_columns,
            min_coverage=qc_min_coverage,
            max_contam=qc_max_contam,
        )
        .pipe(add_quality_control_status_column)
    )

    log.info("Amplicon-Sample QC Statistics:")
    field_coverage_summary = compute_field_coverage_summary(coverage_df)
    low_cov_perc = 100 * field_coverage_summary.n_lowcov / field_coverage_summary.n
    contam_perc = 100 * field_coverage_summary.n_contam / field_coverage_summary.n
    pass_perc = 100 * field_coverage_summary.n_pass / field_coverage_summary.n
    log.info(
        f"  Coverage below <{qc_min_coverage}x: {field_coverage_summary.n_lowcov} ({low_cov_perc:.2f}%)"
    )
    log.info(
        f"  Contamination >{qc_max_contam}: {field_coverage_summary.n_contam} ({contam_perc}%)"
    )
    log.info(f"  Passing QC: {field_coverage_summary.n_pass} ({pass_perc}%)")
    coverage_df.to_csv(f"{output_dir}/summary.coverage.csv", index=False)

    replicates_qc_df = replicates_qc(coverage_df, qc_replicate_passing_threshold)
    replicates_qc_df.to_csv(f"{output_dir}/summary.replicates_qc.csv", index=False)

    samples_summary_df = samples_qc(master_metadata_df, replicates_qc_df)
    samples_summary_df.to_csv(f"{output_dir}/summary.samples_qc.csv", index=False)

    replicates_amplicon_qc_df = replicates_amplicon_qc(coverage_df)
    samples_by_amplicon_summary_df = amplicons_qc_summary(
        master_metadata_df, replicates_amplicon_qc_df
    )
    samples_by_amplicon_summary_df.to_csv(
        f"{output_dir}/summary.samples_amplicons_qc.csv", index=False
    )

    experiment_qc_summary_df = experiment_qc_summary(replicates_amplicon_qc_df)
    experiment_qc_summary_df.to_csv(
        f"{output_dir}/summary.experiments_qc.csv", index=False
    )
    log.info("Quality control complete.")

    # --------------------------------------------------------------------------------
    # Let's move onto to variant calling results
    #
    # --------------------------------------------------------------------------------

    log.info("Loading variants...")
    # variant_df = load_and_concat_variants(expt_dirs)

    # if "sample_id" in variant_df.columns:
    #     variant_df.drop(columns=["sample_id"], inplace=True)

    timer = Timer()
    timer.start()
    variant_df = load_variants_from_vcfs(
        expt_dirs,
        caller=caller,
        output_dir=output_dir,
        summary_regions=regions,
        reference_name="Pf3D7",
    )
    timer.time("Loading and annotating variants from VCFs")

    variant_df.to_csv(f"{output_dir}/intermediate.summary.variants.csv", index=False)

    # # Merge with the quality control results, then we can subset to the analysis set
    variant_df = pd.merge(
        left=coverage_df.rename({"name": "amplicon"}, axis=1)[
            ["expt_name", "barcode", "sample_id", "sample_type", "amplicon", "status"]
        ],
        right=variant_df,
        on=["expt_name", "barcode", "amplicon"],
    )
    timer.time("Merging variants with coverage data")

    log.info("Filtering to analysis set...")
    remove_amplicons = panel_settings.excluded_amplicons  # noqa: F841 later used in query
    remove_mutations = panel_settings.filtered_mutations  # noqa: F841 later used in query
    analysis_df = (
        variant_df.query("status == 'pass'")
        .query("amplicon not in @remove_amplicons")
        .query("mutation not in @remove_mutations")
    )
    timer.time("Filtering to analysis set")

    # Filter out false positives
    analysis_df = filter_false_positives(analysis_df, min_obs=1, min_wsaf=0.33)
    analysis_df.to_csv(f"{output_dir}/summary.variants.analysis_set.csv", index=False)
    timer.time("Filtering false positives")

    # Then we will compute prevalence
    prev_df = compute_variant_prevalence(analysis_df)
    prev_df.to_csv(f"{output_dir}/summary.variants.prevalence.csv", index=False)
    timer.time("Computing variant prevalence")

    for col in prevalence_by:
        prev_by_col_df = compute_variant_prevalence(
            analysis_df, master_metadata_df, [METADATA_COLUMN_PREFIX + col]
        ).pipe(rename_prevalence_by_cols, METADATA_COLUMN_PREFIX, "by_")
        prev_by_col_df.to_csv(
            f"{output_dir}/summary.variants.prevalence-{col}.csv", index=False
        )

    # --------------------------------------------------------------------------------
    # Gene deletion analysis
    #
    # --------------------------------------------------------------------------------

    if panel_settings.deletion_genes:
        log.info("Calculate gene deletions...")
        gene_deletion_df = gene_deletions(coverage_df, panel_settings.deletion_genes)
        gene_deletion_df.to_csv(f"{output_dir}/summary.gene_deletions.csv", index=False)

        prev_gen_deletions_df = gene_deletion_prevalence_by(
            gene_deletion_df, master_metadata_df, []
        )
        prev_gen_deletions_df.to_csv(
            f"{output_dir}/summary.gene-deletions.prevalence.csv", index=False
        )

        for col in prevalence_by:
            prev_gen_deletion_by_col_df = gene_deletion_prevalence_by(
                gene_deletion_df, master_metadata_df, [METADATA_COLUMN_PREFIX + col]
            ).pipe(rename_prevalence_by_cols, METADATA_COLUMN_PREFIX, "by_")
            prev_gen_deletion_by_col_df.to_csv(
                f"{output_dir}/summary.gene-deletions.prevalence-{col}.csv", index=False
            )

    master_metadata_df.to_csv(f"{output_dir}/metadata.csv", index=False)

    log.info("Copy relevant files to summary output directory...")

    # Copy relevant files
    if expts:
        shutil.copy(
            expts[0].regions.path,
            os.path.join(output_dir, os.path.basename(expts[0].regions.path)),
        )
    for map_name in maps:
        file = Path(workspace.path) / "maps" / f"{map_name}.geojson"
        if file.exists():
            shutil.copy(file, f"{output_dir}/{map_name.split('-')[-1]}.geojson")
    coords_file = f"{workspace.get_metadata_dir()}/{summary_name}.coords.csv"
    if os.path.isfile(coords_file):
        shutil.copy(
            coords_file,
            os.path.join(output_dir, "coords.csv"),
        )
    if os.path.isfile(settings_file_path):
        shutil.copy(
            settings_file_path,
            os.path.join(output_dir, "settings.yaml"),
        )

    log.info("Summary analysis complete.")

    timer.report()

    if show_dashboard:
        view(output_dir, summary_name, host=host, port=port)


def view(input_dir: Path, summary_name: str, host: str, port: Optional[int]) -> None:
    """
    View the summary dashboard for a given summary
    """
    print(f'View summary dashboard for "{summary_name}".')
    settings: Settings = Settings()
    settings_file = Path(f"{input_dir}/settings.yaml")
    if settings_file.exists():
        print(f"Loading settings from {settings_file}...")
        settings = load_settings(settings_file)

    bed_files = glob.glob(f"{input_dir}/*.bed")
    if bed_files:
        panel_name = os.path.basename(bed_files[0]).split(".")[0]
        print(f"Use panel name from regions BED file: {panel_name}")
        amplicons = RegionBEDParser(bed_files[0]).names
    else:
        raise ValueError("No regions BED file found in summary directory.")

    panel_settings = get_panel_settings(panel_name)
    amplicon_sets = panel_settings.amplicon_sets
    deletion_genes = panel_settings.deletion_genes

    print(f"Load data from {input_dir}...")

    dashboard = BasicSummaryDashboard(
        summary_name,
        throughput_csv=f"{input_dir}/summary.throughput.csv",
        samples_csv=f"{input_dir}/summary.samples_qc.csv",
        samples_amplicons_csv=f"{input_dir}/summary.samples_amplicons_qc.csv",
        coverage_csv=f"{input_dir}/summary.experiments_qc.csv",
        analysis_csv=f"{input_dir}/summary.variants.analysis_set.csv",
        gene_deletions_csv=f"{input_dir}/summary.gene_deletions.csv",
        master_csv=f"{input_dir}/metadata.csv",
        geojson_glob=f"{input_dir}/*.geojson",
        location_coords_csv=f"{input_dir}/coords.csv",
        settings=settings,
        amplicons=amplicons,
        amplicon_sets=amplicon_sets,
        deletion_genes=deletion_genes,
    )

    print("")
    print("Launching dashboard (press CNTRL+C to exit):")
    print("")
    if port is None:
        port = next_free_port(8050)
    debug = bool(os.getenv("NOMADIC_DEBUG"))
    dashboard.run(debug=debug, auto_open=not debug, host=host, port=port)


class Timer:
    """
    Simple timer class for measuring execution time of code blocks
    """

    def __init__(self):
        self.start_time = None
        self.timings = {}

    def start(self):
        self.start_time = time.time()

    def time(self, name: str):
        if self.start_time is None:
            raise RuntimeError("time can only be called after start")
        end_time = time.time()
        elapsed_time = end_time - self.start_time
        self.timings[name] = elapsed_time
        self.start_time = end_time

    def report(self):
        print("Execution time report:")
        for name, timing in self.timings.items():
            print(f"  {name}: {timing:.2f} seconds")
