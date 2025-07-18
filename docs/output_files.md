
## Output files

*Nomadic* produces several files that provide information about the quality of your sequencing run, and the variants that were detected. They can be found in the results directory (`results/<expt_name>`), and are described briefly below:

### `summary.bam_flagstats.csv`
The `summary.bam_flagstats.csv` contains information about read mapping for each sample.

Each row corresponds to a sample, which can be identified by its barcode (e.g. `barcode01`). The values provided come from [`samtools flagstats`](https://www.htslib.org/doc/samtools-flagstat.html) and are counts of read mapping flags inside of the sample's BAM file. Together they summarise how well your sequencing reads mapped to your reference genome.

| Column | Description |
| --- | --- |
| `barcode` | Sample barcode. |
| `n_total` | Total number of read alignments. |
| `n_mapped` | Number of mapped read alignments. |
| `n_unmapped` | Number of unmapped reads. |
| `n_primary` | Number of reads mapping uniquely. |
| `n_secondary` | Number of reads mapping to more than one location. |
| `n_chimera` | Number of reads mapping as chimeras. |


### `summary.bedcov.csv`
The `summary.bedcov.csv` file contains information about sequencing coverage over each amplicon in each sample.


Each row contains information about coverage over a specific amplicon in a specific sample. The sample is indicated by its barcode (e.g. `barcode01`) and the amplicon is indicated by its `name` (as well as position and length). The values provided come from [`samtools bedcov`](https://www.htslib.org/doc/samtools-bedcov.html).

| Column | Description |
| --- | --- |
| `barcode` | Sample barcode. |
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


### `summary.variants.csv`
The `summary.variants.csv`contains preliminary information about the variants identified in each sample.


Each row contains information about the genotype, depth, quality, and within-sample allele frequency (WSAF) of a specific single-nucleotide polymorphism (SNP) in a specific sample. For all samples the same set of SNPs are described. The set of SNPs described includes all SNPs where at *at least one* sample carried the alternative allele. Note this file is only generated when the `nomadic realtime ... --call` flag is used.

| Column | Description |
| --- | --- |
| `barcode` | Sample barcode. |
| `chrom` | Chromosome of SNP. |
| `pos` | Position of SNP. |
| `ref` | Reference nucleotide for SNP. |
| `alt` | Alternative nucleotide for SNP. |
| `qual` | Variant quality score of SNP. |
| `mut_type` | Type of mutation caused by SNP, e.g. synonymous or non-synonymous. |
| `aa_change` | Amino acid change caused by SNP. For synonymous mutations, we still report (e.g. `V380V`). |
| `aa_pos` | Amino acid number containing the SNP. |
| `strand` | Strand of gene containing the SNP. |
| `amplicon` | Name of amplicon containing SNP.  This comes from the fourth column of BED file used when running `nomadic realtime` (e.g. `-b` flag). |
| `gt` | Called SNP genotype for the sample. Can be reference (`0/0`), heterozygous (`0/1`), homozygous (`1/1`) or failed QC (`./.`). **Note: these are from `bcftools call` and assume a diploid genome.** |
| `gq` | SNP genotype quality for the sample. Note this is different than variant quality (`qual`) as it refers to the quality of the genotype call, rather than whether or not the site is variable. |
| `dp` | Sequencing depth over the SNP. Equivalent to coverage. |
| `wsaf` | Within-sample alternative allele frequency (`ad_alt / (ad_alt + ad_ref)`), where `ad_ref` and `ad_alt` are the depths of the reference and alternative allele. |


## Sharing data
The results folder for an experiment (`results/<expt_name>`) contains all the outputs from *Nomadic*. 

Only the CSV files starting `summary` and the `metadata` folder are required to relaunch the dashboard with the `nomadic dashboard` command (see [Basic Usage](basic.md)).