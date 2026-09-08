
## Output files

*Nomadic* produces several files that provide information about the quality of your sequencing run, and the variants that were detected. They can be found in the results directory (`results/<expt_name>`), and are described briefly below:

### `summary.read_mapping.csv`
The `summary.read_mapping.csv` contains information about read mapping for each sample.

Each row corresponds to a sample, which can be identified by its barcode (e.g. `barcode01`). The values provided come from [`samtools flagstats`](https://www.htslib.org/doc/samtools-flagstat.html) and are counts of read mapping flags inside of the sample's BAM file. Together they summarise how well your sequencing reads mapped to your reference genome.

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

For more information about read mapping, please see [Understanding the Dashboard](understand.md#read-mapping-statistics).

### `summary.region_coverage.csv`
The `summary.region_coverage.csv` file contains information about sequencing coverage over each amplicon in each sample.


Each row contains information about coverage over a specific amplicon in a specific sample. The sample is indicated by its barcode (e.g. `barcode01`) and the amplicon is indicated by its `name` (as well as position and length). The values provided come from [`samtools bedcov`](https://www.htslib.org/doc/samtools-bedcov.html).

| Column | Description |
| --- | --- |
| `barcode` | Sample barcode. |
| `sample_id` | Sample ID. |
| `chrom` | Chromosome of amplicon. |
| `start` | Start position of amplicon. |
| `end` | End position of amplicon. |
| `length` | Length of amplicon. |
| `name` | Name of amplicon. This comes from the fourth column of BED file used when running `nomadic realtime` (e.g. `-b` flag). |
| `n_reads` | Number of reads mapping to the amplicon. |
| `mean_cov` | Mean coverage over the amplicon. |
| `cov_gr100` | Number of positions within the amplicon having >100x coverage. |
| `per_cov_gr100` | Percentage of positions within the amplicon having >100x coverage. |
| `total_cov` | Total coverage over amplicon, i.e. mean coverage times length |

For more information about region coverage, please see [Understanding the Dashboard](understand.md#region-coverage-statistics).


### `summary.aa_changes.csv`
The `summary.aa_changes.csv` file reports missense amino-acid changes identified in each sample. Every observed mutation is represented for every sample, with calls of `mutant`, `mixed`, `absent`, `wt`, `unphased`, or `failed` as appropriate.

| Column | Description |
| --- | --- |
| `barcode` | Sample barcode. |
| `chrom` | Chromosome containing the mutation. |
| `amplicon` | Amplicon containing the mutation. |
| `gene` | Gene affected by the mutation. |
| `aa_pos` | Amino-acid position. |
| `aa_change` | Amino-acid change. |
| `aa_call` | Amino-acid call for the sample. See [Amino-acid calls (`aa_call`)](#amino-acid-calls-aa_call). |
| `aa_dp` | Sequencing depth at the amino-acid position. This is calculated as the minimum depth among all nucleotide positions that make up the codon for the amino-acid. |
| `aa_wsaf` | Within-sample allele frequency for the amino-acid call. This is calculated from the wsaf of the nucleotide positions that make up the codon for the amino-acid. |
| `nt_change` | Nucleotide change or changes underlying the amino-acid change, separated by + signs if multiple. |

#### Amino-acid calls (`aa_call`)

The `aa_call` column in `summary.aa_changes.csv` is a sample-level summary for a particular amino-acid change.

| Value | Meaning |
| --- | --- |
| `mutant` | The mutation is present in the sample in a homozygous/monoclonal form. All informative reads support the alternative amino acid. |
| `mixed` | The mutation is present but only in a subset of reads. This indicates a mixed or heterozygous call within the sample. |
| `absent` | This amino-acid change was considered but is not present in the sample, but there is another amino-acid change at the same position. |
| `wt` | The sample matches the reference at that amino-acid position (wild type) homozygously/monoclonally. |
| `unphased` | More than one mixed nucleotide call exists at the same amino-acid position. As we don't phase yet, we can not determine the correct amino-acid change. |
| `failed` | There was not enough reliable information to classify the site (e.g., low coverage, poor-quality reads, or a strand bias); the call failed or the site was not callable. |

### `summary.nt_changes.csv`
The `summary.nt_changes.csv` file reports nucleotide changes identified in each sample. Every observed nucleotide mutation is represented for every sample.

| Column | Description |
| --- | --- |
| `barcode` | Sample barcode. |
| `chrom` | Chromosome containing the mutation. |
| `pos` | Reference position of the mutation. |
| `amplicon` | Amplicon containing the mutation. |
| `ref` | Reference nucleotide. |
| `alt` | Alternative nucleotide. |
| `dp` | Sequencing depth at the position. This is the number of reads included in that variant call after filtering for quality and mapping criteria. |
| `gt` | Nucleotide call: `mutant`, `mixed`, `absent`, `wt`, or `failed`. (See [Nucleotide calls](#nucleotide-calls-gt)) |
| `wsaf` | Within-sample allele frequency. (see [Understanding the Dashboard](understand.md#variant-calling) for more details) |

#### Nucleotide calls (`gt`)

The `gt` column in `summary.nt_changes.csv` records the call for a specific nucleotide position.

| Value | Meaning |
| --- | --- |
| `mutant` | The genotype is homozygous/monoclonal alternative at this position. |
| `mixed` | The genotype is heterozygous or otherwise mixed at this position.|
| `absent` | The variant is not present in this sample at this position, but there is another variant at the same position. |
| `wt` | The position is wild type and matches the reference genotype. |
| `failed` | The genotype could not be determined reliably from the sequencing data (e.g., low coverage, poor-quality reads, or a strand bias). |

These values are the same categories used by the dashboard heatmaps and the variant-calling summaries, and they provide a compact way to interpret whether a sample carries a mutation, is mixed, or lacks enough information for a call.

For more information about variant calling, please see [Understanding the Dashboard](understand.md#variant-calling).
