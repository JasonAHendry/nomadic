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

## Testing
I have created a script to simulate a small nanopore sequencing run, that allows you to test `nomadic3 realtime` without having an actual sequencing experiment running. To try this, first run: 

```
./scripts/run_realtime.sh
```

Then, in a second terminal window, and run:
```
python scripts/simulate_sequencing.py
```

## Acknowledgements
This work was funded by the Bill and Melinda Gates Foundation (INV-003660, INV-048316).




