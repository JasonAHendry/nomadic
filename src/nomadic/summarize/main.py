import glob
import os
import shutil
from pathlib import Path
from typing import Optional

from nomadic.summarize.analysis.deletions import (
    gene_deletion_prevalence_by,
    gene_deletions,
)
from nomadic.summarize.analysis.input_data import common_caller, common_reference_name
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
from nomadic.summarize.analysis.variants import (
    compute_variant_prevalence,
    filter_to_analysis_set,
    load_variants_from_vcfs,
    remove_false_positives,
    remove_never_observed_nt_variants,
    remove_never_observed_variants,
    rename_prevalence_by_cols,
)
from nomadic.summarize.dashboard.builders import BasicSummaryDashboard
from nomadic.summarize.dir_structure import DirStructure, looks_like_summary_dir
from nomadic.util.dirs import produce_dir
from nomadic.util.experiment import (
    experiment_outputs,
)
from nomadic.util.logging_config import LoggingFascade
from nomadic.util.panel import get_panel_settings
from nomadic.util.port import next_free_port
from nomadic.util.regions import RegionBEDParser, common_regions
from nomadic.util.summary_settings import (
    Settings,
    load_settings,
)
from nomadic.util.timer import Timer
from nomadic.util.workspace import Workspace

# --------------------------------------------------------------------------------
# Variant analysis
#
# --------------------------------------------------------------------------------


# --------------------------------------------------------------------------------
# Main
#
# --------------------------------------------------------------------------------
def main(
    *,
    workspace: Workspace,
    output_dir: Path,
    expt_dirs: list[Path],
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
        experiment_outputs(
            expt_dir, allow_missing_files=["depth", "fastq", "nt_changes"]
        )
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

    summary_dir_structure = DirStructure(summary_dir=output_dir)

    for dir in summary_dir_structure.dirs:
        produce_dir(dir)

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
    full_inventory_df.to_csv(summary_dir_structure.inventory_file, index=False)

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
    throughput_df.to_csv(summary_dir_structure.throughput_file, index=True)

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
    coverage_df.to_csv(summary_dir_structure.coverage_file, index=False)

    replicates_qc_df = replicates_qc(coverage_df, qc_replicate_passing_threshold)
    replicates_qc_df.to_csv(summary_dir_structure.replicates_qc_file, index=False)

    samples_summary_df = samples_qc(master_metadata_df, replicates_qc_df)
    samples_summary_df.to_csv(summary_dir_structure.samples_qc_file, index=False)

    replicates_amplicon_qc_df = replicates_amplicon_qc(coverage_df)
    samples_by_amplicon_summary_df = amplicons_qc_summary(
        master_metadata_df, replicates_amplicon_qc_df
    )
    samples_by_amplicon_summary_df.to_csv(
        summary_dir_structure.samples_by_amplicon_qc_file, index=False
    )

    experiment_qc_summary_df = experiment_qc_summary(replicates_amplicon_qc_df)
    experiment_qc_summary_df.to_csv(
        summary_dir_structure.experiment_qc_file, index=False
    )
    log.info("Quality control complete.")

    # --------------------------------------------------------------------------------
    # Let's move onto variant calling results
    #
    # --------------------------------------------------------------------------------
    timer = Timer("Variant analysis")
    timer.start()

    log.info("Loading variants...")
    aa_changes_df, nt_changes_df = load_variants_from_vcfs(
        expt_dirs,
        caller=caller,
        temp_dir=Path(output_dir) / "temp_vcf_processing",
        filtered_vcf=summary_dir_structure.vcfs_dir
        / "summary.variants.filtered.vcf.gz",
        annotated_vcf=summary_dir_structure.vcfs_dir
        / "summary.variants.annotated,vcf.gz",
        bed_path=Path(regions.path),
        reference_name=reference_name,
        exclude_amplicons=panel_settings.excluded_amplicons,
        exclude_mutations=panel_settings.filtered_mutations,
        log=log,
    )
    timer.time("Loading and annotating variants from VCFs")

    # # Merge with the quality control results, then we can subset to the analysis set
    aa_changes_df = filter_to_analysis_set(
        aa_changes_df,
        coverage_df=coverage_df,
    )
    timer.time("Filtering to analysis set")
    nt_changes_df = filter_to_analysis_set(
        nt_changes_df,
        coverage_df=coverage_df,
    )
    timer.time("Filtering nt set")

    # Filter out false positives
    aa_changes_df = remove_false_positives(aa_changes_df, min_obs=1, min_aa_wsaf=0.33)
    timer.time("Filtering false positives")
    aa_changes_df = remove_never_observed_variants(aa_changes_df)
    nt_changes_df = remove_never_observed_nt_variants(nt_changes_df)
    timer.time("Removing never observed variants")

    if panel_settings.amplicon_sets:
        for set_name, amplicons in panel_settings.amplicon_sets.items():
            set_df = aa_changes_df[aa_changes_df["amplicon"].isin(amplicons)]
            set_df.to_csv(
                summary_dir_structure.variants_dir
                / f"summary.aa_changes.{set_name.lower()}.csv",
                index=False,
            )
    else:
        aa_changes_df.to_csv(
            summary_dir_structure.aa_changes_file,
            index=False,
        )
    timer.time("Writing aa_changes to CSV")

    nt_changes_df.to_csv(
        summary_dir_structure.nt_changes_file,
        index=False,
    )
    timer.time("Writing nt changes to CSV")

    # Then we will compute prevalence
    prev_df = compute_variant_prevalence(aa_changes_df)
    timer.time("Computing variant prevalence")
    prev_df.to_csv(
        summary_dir_structure.prevalence_dir / "summary.aa_changes.prevalence.csv",
        index=False,
    )
    timer.time("Writing aa changes prevalence to CSV")

    for col in prevalence_by:
        prev_by_col_df = compute_variant_prevalence(
            aa_changes_df, master_metadata_df, [METADATA_COLUMN_PREFIX + col]
        ).pipe(rename_prevalence_by_cols, METADATA_COLUMN_PREFIX, "by_")
        prev_by_col_df.to_csv(
            summary_dir_structure.prevalence_dir
            / f"summary.aa_changes.prevalence-{col}.csv",
            index=False,
        )

    # --------------------------------------------------------------------------------
    # Gene deletion analysis
    #
    # --------------------------------------------------------------------------------

    if panel_settings.deletion_genes:
        log.info("Calculate gene deletions...")
        produce_dir(summary_dir_structure.gene_deletions_dir)
        gene_deletion_df = gene_deletions(coverage_df, panel_settings.deletion_genes)
        gene_deletion_df.to_csv(
            summary_dir_structure.gene_deletions_dir / "summary.gene_deletions.csv",
            index=False,
        )

        prev_gen_deletions_df = gene_deletion_prevalence_by(
            gene_deletion_df, master_metadata_df, []
        )
        prev_gen_deletions_df.to_csv(
            summary_dir_structure.gene_deletions_dir
            / "summary.gene-deletions.prevalence.csv",
            index=False,
        )

        for col in prevalence_by:
            prev_gen_deletion_by_col_df = gene_deletion_prevalence_by(
                gene_deletion_df, master_metadata_df, [METADATA_COLUMN_PREFIX + col]
            ).pipe(rename_prevalence_by_cols, METADATA_COLUMN_PREFIX, "by_")
            prev_gen_deletion_by_col_df.to_csv(
                summary_dir_structure.gene_deletions_dir
                / f"summary.gene-deletions.prevalence-{col}.csv",
                index=False,
            )

    master_metadata_df.pipe(
        rename_prevalence_by_cols, METADATA_COLUMN_PREFIX, ""
    ).to_csv(summary_dir_structure.metadata_file, index=False)

    log.info("Copy relevant files to summary output directory...")

    # Copy relevant files
    if expts:
        shutil.copy(
            expts[0].regions.path,
            summary_dir_structure.panel_info_dir
            / os.path.basename(expts[0].regions.path),
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

    dir = DirStructure(summary_dir=input_dir)

    settings: Settings = Settings()
    settings_file = Path(f"{input_dir}/settings.yaml")
    if settings_file.exists():
        print(f"Loading settings from {settings_file}...")
        settings = load_settings(settings_file)

    bed_file = next(dir.panel_info_dir.glob("*.bed"), None)
    if bed_file is not None:
        panel_name = bed_file.stem.removesuffix(".amplicons")
        print(f"Use panel name from regions BED file: {panel_name}")
        amplicons = RegionBEDParser(str(bed_file)).names
    else:
        raise ValueError("No regions BED file found in summary directory.")

    panel_settings = get_panel_settings(panel_name)
    amplicon_sets = panel_settings.amplicon_sets
    deletion_genes = panel_settings.deletion_genes

    print(f"Load data from {input_dir}...")

    dashboard = BasicSummaryDashboard(
        summary_name,
        throughput_csv=str(dir.throughput_file),
        samples_csv=str(dir.samples_qc_file),
        samples_amplicons_csv=str(dir.samples_by_amplicon_qc_file),
        experiment_qc_csv=str(dir.experiment_qc_file),
        aa_changes_csvs=list(glob.glob(f"{dir.variants_dir}/summary.aa_changes*.csv")),
        gene_deletions_csv=str(dir.gene_deletions_dir / "summary.gene_deletions.csv"),
        master_csv=str(dir.metadata_file),
        geojson_glob=f"{input_dir}/*.geojson",
        location_coords_csv=f"{input_dir}/coords.csv",
        settings=settings,
        amplicons=amplicons,
        amplicon_sets=amplicon_sets,
        deletion_genes=deletion_genes,
    )

    print()
    print("Launching dashboard (press CNTRL+C to exit):")
    print()
    if port is None:
        port = next_free_port(8050)
    debug = bool(os.getenv("NOMADIC_DEBUG"))
    dashboard.run(debug=debug, auto_open=not debug, host=host, port=port)
