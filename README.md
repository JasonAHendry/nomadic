[![build](https://github.com/JasonAHendry/nomadic/actions/workflows/build.yml/badge.svg)](https://github.com/JasonAHendry/nomadic/actions/workflows/build.yml)
[![View - Docs](https://img.shields.io/badge/View-Docs-blue?logo=materialformkdocs&logoColor=blue)](https://jasonahendry.github.io/nomadic/)
[![Conda Downloads](https://img.shields.io/conda/dn/bioconda/nomadic?color=green&link=https%3A%2F%2Fanaconda.org%2Fbioconda%2Fnomadic)](https://anaconda.org/bioconda/nomadic)
![OS - Linux | OSX](https://img.shields.io/badge/OS-Linux_|_OSX-informational)
<p align="center"><img src="docs/img/home/nomadic_logo.png" width="500"></p>

## Overview
*Nomadic* is a real-time bioinformatics pipeline and dashboard for nanopore sequencing data. While sequencing is still ongoing, it performs read mapping and sample quality control, as well as variant calling and annotation. This information is displayed in real-time to a graphical dashboard that has interactive features.

## Features
 - Real-time read mapping with Minimap2.
 - Real-time sample quality control and amplicon coverage evaluation.
 - Real-time variant calling with bcftools. These calls are preleiminary; treat with caution.
 - Support for different reference genomes or amplicons panels.

 ## Installation
 Nomadic can be installed from [bioconda](https://anaconda.org/bioconda/nomadic)
 ```
 conda install bioconda::nomadic
 ```

 ## Quickstart
Navigate to a directory where your nomadic files should live and setup a workspace with
```
nomadic start pfalciparum
```

Afterwards, navigate to the newly created workspace, create a [metadata file](https://jasonahendry.github.io/nomadic/basic/#using-nomadic-for-real-time-analysis) and start your experiment.

```
cd nomadic
nomadic realtime <expt_name>
```

For more detailed information, see our [Documentation](https://jasonahendry.github.io/nomadic)

## Development
For development, clone the repo, create a conda environment and activate it

```
conda env create -f environments/dev.yml
conda activate nomadic-dev
```

Install the package locally in development mode

```
pip install -e .
```

## Testing
I have created a script to simulate a small nanopore sequencing run, that allows you to test `nomadic realtime` without having an actual sequencing experiment running. To try this, first run: 

```
./scripts/run_realtime.sh
```

Then, in a second terminal window, and run:
```
python scripts/simulate_sequencing.py
```

Additionally you can run our automated test from the root with

```
pytest
```

## Building the documentation
The documentation was created using [MkDocs](https://www.mkdocs.org/). You can serve it locally with:

```
conda env create -f environments/dev.yml
mkdocs serve
```

## Acknowledgements
This work was funded by the Bill and Melinda Gates Foundation (INV-003660, INV-048316).




