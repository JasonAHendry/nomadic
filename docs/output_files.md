
*Nomadic* produces several files that provide information about the quality of your sequencing run, and the variants that were detected.


*Nomadic* produces output files in two contexts:

- **Realtime experiment outputs**, written to `results/<expt_name>/` by
	`nomadic realtime`.
- **Summarize outputs**, written to a summary directory by `nomadic summarize`
	when combining multiple experiments. See [Summarizing Experiments](summary.md)
	for the workflow and dashboard.

## Realtime experiment output files

### `summary.read_mapping.csv`

This file contains read-mapping statistics for each barcode. The values come
from [`samtools flagstat`](https://www.htslib.org/doc/samtools-flagstat.html)
and summarize how well reads mapped to the reference genome.

| Column | Description |
| --- | --- |
| `barcode` | Sample barcode. |
| `sample_id` | Sample ID. |
| `n_total` | Total number of read alignments. |
| `n_mapped` | Number of mapped read alignments. |
| `n_unmapped` | Number of unmapped reads. |
| `n_primary` | Number of reads mapping uniquely. |
| `n_secondary` | Number of reads mapping to more than one location. |
| `n_supplementary` | Number of reads mapping as chimeras. |

For more information about read mapping, see [Understanding the Dashboard](understand.md#read-mapping-statistics).

### `summary.region_coverage.csv`

This file contains coverage over each amplicon in each sample. The values
come from [`samtools bedcov`](https://www.htslib.org/doc/samtools-bedcov.html).

| Column | Description |
| --- | --- |
| `barcode` | Sample barcode. |
| `sample_id` | Sample ID. |
| `chrom` | Chromosome of the amplicon. |
| `start` | Start position of the amplicon. |
| `end` | End position of the amplicon. |
| `length` | Length of the amplicon. |
| `name` | Amplicon name from the fourth column of the BED file. |
| `n_reads` | Number of reads mapping to the amplicon. |
| `mean_cov` | Mean coverage over the amplicon. |
| `cov_gr100` | Number of positions within the amplicon having greater than 100x coverage. |
| `per_cov_gr100` | Percentage of positions within the amplicon having greater than 100x coverage. |
| `total_cov` | Total coverage over the amplicon, equal to mean coverage times length. |

For more information about region coverage, see [Understanding the Dashboard](understand.md#region-coverage-statistics).

### `summary.aa_changes.csv` {#realtime-aa-changes}

This file reports missense amino-acid changes identified in each sample. Every
observed mutation is represented for every sample, with calls of `mutant`,
`mixed`, `absent`, `wt`, `unphased`, or `failed` as appropriate.

| Column | Description |
| --- | --- |
| `barcode` | Sample barcode. |
| `chrom` | Chromosome containing the mutation. |
| `amplicon` | Amplicon containing the mutation. |
| `gene` | Gene affected by the mutation. |
| `aa_pos` | Amino-acid position. |
| `aa_change` | Amino-acid change. |
| `aa_call` | Amino-acid call for the sample. See [Amino-acid calls (`aa_call`)](#amino-acid-calls-aa_call). |
| `aa_dp` | Sequencing depth at the amino-acid position, calculated as the minimum depth among the nucleotide positions making up the codon. |
| `aa_wsaf` | Within-sample allele frequency for the amino-acid call. |
| `nt_change` | Nucleotide change or changes underlying the amino-acid change, separated by `+` signs if multiple. |

#### Amino-acid calls (`aa_call`)

| Value | Meaning |
| --- | --- |
| `mutant` | The mutation is present in a homozygous or monoclonal form. |
| `mixed` | The mutation is present in only a subset of reads. |
| `absent` | This amino-acid change is not present, but another amino-acid change occurs at the same position. |
| `wt` | The sample matches the reference at that amino-acid position. |
| `unphased` | More than one mixed nucleotide call exists at the same amino-acid position, so the correct amino-acid change cannot be determined. |
| `failed` | There is not enough reliable information to classify the site. |

### `summary.nt_changes.csv` {#realtime-nt-changes}

This file reports nucleotide changes identified in each sample. Every observed
nucleotide mutation is represented for every sample.

| Column | Description |
| --- | --- |
| `barcode` | Sample barcode. |
| `chrom` | Chromosome containing the mutation. |
| `pos` | Reference position of the mutation. |
| `amplicon` | Amplicon containing the mutation. |
| `ref` | Reference nucleotide. |
| `alt` | Alternative nucleotide. |
| `gene` | Gene containing the position, when annotated. |
| `aa_pos` | Amino-acid position, when annotated. |
| `dp` | Sequencing depth at the position. |
| `gt` | Nucleotide call: `mutant`, `mixed`, `absent`, `wt`, or `failed`. See [Nucleotide calls (`gt`)](#nucleotide-calls-gt). |
| `wsaf` | Within-sample allele frequency. |

#### Nucleotide calls (`gt`)

| Value | Meaning |
| --- | --- |
| `mutant` | The genotype is homozygous or monoclonal alternative at this position. |
| `mixed` | The genotype is heterozygous or otherwise mixed at this position. |
| `absent` | The variant is not present, but another variant occurs at the same position. |
| `wt` | The position matches the reference genotype. |
| `failed` | The genotype could not be determined reliably. |

These call categories are also used by the dashboard heatmaps. For more
information about variant calling, see [Understanding the Dashboard](understand.md#variant-calling).

## Summarize output files

The files below are written by `nomadic summarize` when it combines completed
experiments. Only the CSV outputs are described here. The `<set>` placeholder
means an amplicon set configured for the panel, and `<column>` means a metadata
column supplied with `--prevalence-by`. Files containing either placeholder are
created only when the corresponding option or panel configuration is used.

### `seq_inventory/sample_inventory.csv`

This file lists the samples found across the input experiments and records
whether each sample is included, excluded, or a control.

| Column | Description |
| --- | --- |
| `expt_name` | Experiment name. |
| `barcode` | Sample barcode. |
| `sample_id` | Sample identifier. |
| `sample_type` | Sample type, such as `field`, `pos`, or `neg`. |
| `status` | Inventory status: `included`, `excluded`, or `control`. |
| Other metadata columns | Metadata columns carried into the inventory when available. |

### `seq_inventory/samples.by_experiment.csv`

This file is a throughput cross-tabulation. Its index identifies an experiment
or an aggregate row such as `All` or `included_expts`. The columns are
`All`, `pos`, `neg`, `field_included`, and `field_unique`.

### `quality_control/coverage.csv`

This file contains one row per amplicon and barcode after merging experiment
coverage with the sample inventory. It is the detailed QC input for the other
QC summaries.

| Column | Description |
| --- | --- |
| `expt_name` | Experiment name. |
| `barcode` | Sample barcode. |
| `sample_id` | Sample identifier. |
| `sample_type` | Sample type. |
| `chrom`, `start`, `end`, `length`, `name` | Amplicon coordinates and name. |
| `n_reads`, `mean_cov`, `cov_gr100`, `per_cov_gr100`, `total_cov` | Amplicon coverage measurements from the experiment output. |
| `mean_cov_neg` | Mean coverage of the negative control for the experiment and amplicon. |
| `fail_lowcov` | Whether mean coverage is below the configured minimum. |
| `fail_contam_rel` | Whether relative contamination exceeds the configured maximum. |
| `fail_contam_abs` | Whether negative-control coverage reaches the configured absolute threshold. |
| `fail_contam` | Whether the amplicon fails either contamination check. |
| `passing` | Whether the amplicon passes coverage and contamination QC. |
| `status` | Detailed status, such as `pass`, `lowcov`, `contam`, `duplicate`, or `control`. |

### `quality_control/replicates_qc.csv`

This file summarizes QC for each field-sample replicate, identified by
experiment, barcode, and sample ID.

| Column | Description |
| --- | --- |
| `expt_name`, `barcode`, `sample_id` | Replicate identifiers. |
| `n_amplicons` | Number of amplicons evaluated. |
| `n_passing` | Number of passing amplicons. |
| `n_fail_contam` | Number of amplicons failing contamination QC. |
| `n_fail_lowcov` | Number of amplicons failing low-coverage QC. |
| `passing` | Whether the replicate meets the configured passing-amplicon threshold. |

### `quality_control/samples_qc.csv`

This file summarizes all samples in the master metadata, including samples
that were not sequenced.

| Column | Description |
| --- | --- |
| `sample_id` | Sample identifier. |
| `n_replicates` | Number of sequenced replicates. |
| `n_passing` | Number of passing replicates. |
| `status` | `passing`, `failing`, or `not_sequenced`. |

### `quality_control/samples_amplicons_qc.csv`

This file summarizes replicate QC for each sample and amplicon. It contains
`sample_id`, `name`, `n_replicates`, `n_passing`, and `status`, where `status`
uses the same values as `samples_qc.csv`.

### `quality_control/experiments_qc.csv`

This file summarizes field-sample QC for each experiment and amplicon. It
contains `expt_name`, `name`, `mean_cov_field`, `mean_cov_neg`, `n_field`,
`n_field_contam`, `n_field_lowcov`, `n_field_passing`,
`per_field_contam`, `per_field_lowcov`, and `per_field_passing`.

### `variants/aa_changes.csv` and `variants/aa_changes.<set>.csv`

These files contain amino-acid calls combined across experiments and filtered
to the QC analysis set. They use the same columns and call categories as
[the realtime amino-acid changes file](#realtime-aa-changes). When the panel defines
amplicon sets, one file is written per set instead of the unsuffixed file.

### `variants/nt_changes.csv`

This file contains nucleotide calls combined across experiments and filtered to
the QC analysis set. It uses the same columns and call categories as
[the realtime nucleotide changes file](#realtime-nt-changes).

### `variants/prevalence.aa_changes.csv` and related files

`prevalence.aa_changes.csv` reports amino-acid change prevalence across the
included samples. The corresponding `<set>` file is written for each configured
amplicon set. When `--prevalence-by <column>` is used, the output is named
`prevalence.aa_changes.by-<column>.csv`, or
`prevalence.aa_changes.<set>.by-<column>.csv` for an amplicon set.

| Column | Description |
| --- | --- |
| `chrom`, `amplicon`, `gene`, `aa_pos`, `aa_change` | Amino-acid change identifiers. |
| `n_samples` | Number of sample-level calls included in the group. |
| `n_passed` | Number of calls with an interpretable result. |
| `n_wt` | Number of wild-type calls. |
| `n_mixed` | Number of mixed calls. |
| `n_mut` | Number of mutant calls. |
| `per_wt`, `per_mixed`, `per_mut` | Percentages among passed calls. |
| `prevalence` | Percentage of passed calls that are mixed or mutant. |
| `prevalence_lowci`, `prevalence_highci` | Lower and upper bounds of the 95% beta confidence interval for prevalence. |
| `by_<column>` | Group value when prevalence is stratified by metadata. |

The summary variant tables remove variants never observed in the analysis set
and apply panel-specific filtering before prevalence is calculated.
