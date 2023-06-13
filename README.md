# nomadic3
## Overview
`nomadic3` supports real-time mapping and analysis of amplicon-based nanopore sequencing, rendering the data to a browser-based dashboard.

## Install

#### Requirements

To install `nomadic3`, you will need:
- The version control software [git](https://github.com/git-guides/install-git)
- The package manager [conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html) or [mamba](https://mamba.readthedocs.io/en/latest/installation.html) 
  - Mamba is faster and is recommended

#### Steps

**1.  Clone the repository:**
```
git clone https://github.com/JasonAHendry/nomadic3.git
cd nomadic3
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
**3. Install `nomadic3` and remaining dependencies:**
```
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

`nomadic3` performs real-time mapping to a reference genome. Start by downloading the reference genome for your target organism, e.g.:
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
 -b <path/to/your_regions.bed>
```

Once you run this command, you should get a dashboard link in your terminal, something like:

```
Dash is running on http://127.0.0.1:8050/
```

Copy and paste the link `http://127.0.0.1:8050/` into your web browser to view the dashboard. 

\
Flag information:
- `-e`: Your experiment name. For example, '2023-06-12_seq-nomads8'. This is used to create the output directory.
- `-f`: Path to the `fastq_pass` directory that will be created by `MinKNOW` or `Guppy`. If your experiment has multiple samples, this folder will typically contain folders for each sample, inside of which there are `.fastq` or `.fastq.gz` files.
- `-m`: Path to a metadata CSV file. 
  - This file *must* have a `barcode` and `sample_id` column.
  - Both `barcode` and `sample_id` columns must have only unique entries.
  - See `example_data/metadata/sample_info.csv` for an example.
- `-b`: Path to a BED file defining your regions of interest.
  - Note the BED file specification [here](https://en.wikipedia.org/wiki/BED_(file_format)).
  - This BED file *must* contain a forth column for region name information, and the names must be unique. They are used to generate plots and compute summary statistics.
  - See `resources/beds/nomads8.amplicons.bed` for an example.

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




