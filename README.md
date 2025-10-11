[![build](https://github.com/JasonAHendry/nomadic/actions/workflows/build.yml/badge.svg)](https://github.com/JasonAHendry/nomadic/actions/workflows/build.yml)
[![View - Docs](https://img.shields.io/badge/View-Docs-blue?logo=materialformkdocs&logoColor=blue)](https://jasonahendry.github.io/nomadic/)
[![Conda Downloads](https://img.shields.io/conda/dn/bioconda/nomadic?color=green&link=https%3A%2F%2Fanaconda.org%2Fbioconda%2Fnomadic)](https://anaconda.org/bioconda/nomadic)
![OS - Linux | OSX](https://img.shields.io/badge/OS-Linux_|_OSX-informational)
[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io/recipes/nomadic/README.html)

<p align="center"><img src="docs/img/home/nomadic_logo.png" width="500"></p>

## Overview
*Nomadic* is a real-time bioinformatics pipeline and dashboard for nanopore sequencing data. While sequencing is still ongoing, it performs read mapping and sample quality control, as well as variant calling and annotation. This information is displayed in real-time to a graphical dashboard that has interactive features.

Please visit our [documentation](https://jasonahendry.github.io/nomadic) to learn more.

## Features
 - [x] Real-time read mapping with Minimap2.
 - [x] Real-time sample quality control and amplicon coverage evaluation.
 - [x] Real-time variant calling with bcftools. These calls are preleiminary; treat with caution.
 - [x] Support for different reference genomes or amplicons panels.

 ## Installation
 *Nomadic* can be installed from [bioconda](https://anaconda.org/bioconda/nomadic)
 ```
 conda install bioconda::nomadic
 ```

 ## Quickstart
Navigate to a directory where your nomadic files should live and setup a workspace with

```
nomadic start pfalciparum
```

Afterwards, navigate to the newly created workspace, create a [metadata file](https://jasonahendry.github.io/nomadic/basic/#using-nomadic-for-real-time-analysis) and start your experiment:

```
cd nomadic
nomadic realtime <expt_name>
```

For more detailed information, see our [documentation](https://jasonahendry.github.io/nomadic).

## Development
If you would like to develop *Nomadic* you can install it from source. First, clone the github repository

```
git clone https://github.com/JasonAHendry/nomadic.git
```
and then create the development conda environment, and activated it:

```
conda env create -f environments/dev.yml
conda activate nomadic-dev
```

Finally, install the package locally in development mode:

```
pip install -e .
```

Please note that if new dependencies are added, you will have to update your conda environment. For this, run:

```
conda env update -f environments/dev.yml
```

## Testing
To test *Nomadic*, we have written a small python script that simulates nanopore sequencing:

```
python scripts/simulate_sequencing.py
```

This will move FASTQ files into the directory `example_data/minknow/fastq_pass`, which you can process with `nomadic realtime`.

Additionally, if you followed our development [instructions](#development) you can test with:

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




