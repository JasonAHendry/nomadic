<p align="center"><img src="misc/nomadic_logo-01.png" width="500"></p>

## Overview
*Nomadic* supports real-time mapping and analysis of amplicon-based nanopore sequencing, rendering the data to a browser-based dashboard.

## Features
- [x] Real-time read mapping with [*Minimap2*](https://github.com/lh3/minimap2)
- [x] Real-time sample quality control and amplicon coverage evaluation
- [x] Real-time variant calling with [*bcftools*](https://github.com/samtools/bcftools). These calls are preleiminary; treat with caution.
- [x] Different reference genomes
- [x] Different amplicon panels 

## Install

#### Requirements
<details>
  
To install `nomadic`, you will need:
- The version control software [git](https://github.com/git-guides/install-git)
- The package manager [conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html) or [mamba](https://mamba.readthedocs.io/en/latest/installation.html) 
  - Mamba is faster and is recommended
    
</details>

#### Steps

**1.  Clone the repository:**
```
git clone https://github.com/JasonAHendry/nomadic.git
cd nomadic
```
\
**2.  Install the depedendencies with conda:**
```
conda env create -f environments/run.yml
```
or equivalently, with mamba:
```
mamba env create -f environments/run.yml
```
\
**3. Install `nomadic` and remaining dependencies:**
```
conda activate nomadic
pip install -e .
```
\
**4. Test your installation.**
In the terminal, you should see available commands by typing:
```
nomadic --help
```


## Basic usage

**A. Download your reference genome** 

*Nomadic* performs real-time mapping to a reference genome. Start by downloading the reference genome for your target organism, e.g.:
```
nomadic download -r Pf3D7
```
For the 3D7 reference genome of *P. falciparum*.

\
**B. Run `nomadic realtime`**

Once sequencing has started, you can perform real-time analysis using the `nomadic realtime` command as follows:

```
nomadic realtime \
 -e <your_experiment_name> \
 -f <path/to/fastq_pass> \
 -m <path/to/metadata.csv> \
 -b <path/to/your_regions.bed> \
--call
```

Once you run this command, you should get a dashboard link in your terminal, something like:

```
Dash is running on http://127.0.0.1:8050/
```

Copy and paste the link `http://127.0.0.1:8050/` into your web browser to view the dashboard. 

\
Flag information:
| Flag | Description | Required / Optional |
| ---    | --- | --- |
| ` -e ` | Your experiment name. For example, `2023-06-12_nomads8`. This is used to create the output directory. | Required |
| ` -f ` | Path to the `fastq_pass` directory that will be created by `MinKNOW` or `Guppy`. If your experiment has multiple samples, this folder will typically contain folders for each sample, inside of which there are `.fastq` or `.fastq.gz` files. | Required |
| ` -m ` | Path to a metadata CSV file. This file *must* have a `barcode` and `sample_id` column; and both must contain only unique entries. See [here](example_data/metadata/sample_info.csv) for an example. | Required |
| ` -b ` | Path to a BED file defining your amplicons of interest. This BED file *must* contain a forth column for region name information, and the names must be unique. They are used to generate plots and compute summary statistics. See [here](example_data/beds/nomads8.amplicons.bed) for an example. | Required |
| ` -c ` | Or `--call`. If provided, variant calling will be performed. | Optional |

For a full running example look at `scripts/run_realtime.sh`

## Output files

*Nomadic* produces several files that provide information about the quality of your sequencing run, and the variants that were detected. They can be found in the results directory (`results/<expt_name>`), and are described briefly below:

**`summary.bam_flagstats.csv`** contains information about read mapping for each sample.
<details>

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

</details>

**`summary.bedcov.csv`** contains information about sequencing coverage over each amplicon in each sample.
<details>

Each row contains information about coverage over a specific amplicon in a specific sample. The sample is indicated by its barcode (e.g. `barcode01`) and the amplicon is indicated by its `name` (as well as position and length).

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
| `cov_gr100` | Number of positions within the amplicon having $>100\times$ coverage. |
| `per_cov_gr100` | Percentage of positions within the amplicon having $>100\times$ coverage. |
| `total_cov` | Total coverage over amplicon, i.e. mean coverage times length |

</details>

**`summary.variants.csv`** contains preliminary information about the variants identified in each sample.
<details>

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

</details>

## Testing
I have created a script to simulate a small nanopore sequencing run, that allows you to test `nomadic realtime` without having an actual sequencing experiment running. To try this, first run: 

```
./scripts/run_realtime.sh
```

Then, in a second terminal window, and run:
```
python scripts/simulate_sequencing.py
```

## Acknowledgements
This work was funded by the Bill and Melinda Gates Foundation (INV-003660, INV-048316).




