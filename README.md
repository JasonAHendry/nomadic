[![build](https://github.com/JasonAHendry/nomadic/actions/workflows/build.yml/badge.svg)](https://github.com/JasonAHendry/nomadic/actions/workflows/build.yml)
[![View - Docs](https://img.shields.io/badge/View-Docs-blue?logo=materialformkdocs&logoColor=blue)](https://jasonahendry.github.io/nomadic/)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/nomadic/badges/downloads.svg)](https://anaconda.org/bioconda/nomadic)
<p align="center"><img src="docs/img/home/nomadic_logo.png" width="500"></p>

## Overview
*Nomadic* is a real-time bioinformatics pipeline and dashboard for nanopore sequencing data. While sequencing is still ongoing, it performs read mapping and sample quality control, as well as variant calling and annotation. This information is displayed in real-time to a graphical dashboard that has interactive features.

## Testing
I have created a script to simulate a small nanopore sequencing run, that allows you to test `nomadic realtime` without having an actual sequencing experiment running. To try this, first run: 

```
./scripts/run_realtime.sh
```

Then, in a second terminal window, and run:
```
python scripts/simulate_sequencing.py
```

## Documentation
The documentation was created using [MkDocs](https://www.mkdocs.org/). You can serve it locally with:

```
conda env create -f environments/dev.yml
mkdocs serve
```

## Acknowledgements
This work was funded by the Bill and Melinda Gates Foundation (INV-003660, INV-048316).




