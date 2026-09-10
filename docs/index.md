#

<p align="center">
<img src="img/home/nomadic_logo.png" alt="nomadic" width="80%">
</p>

---

## Overview
*Nomadic* is a real-time bioinformatics pipeline and dashboard for nanopore sequencing data. While sequencing is still ongoing, it performs read mapping and sample quality control, as well as variant calling and annotation. This information is displayed in real-time to a graphical dashboard that has interactive features.


<p align="center">
<img src="img/home/dashboard-workflow.png" alt="dashboard" width="80%">
</p>


It was designed to work with amplicon sequencing data from the NOMADS-MVP protocol, which targets a panel of genes important for the control of *Plasmodium falciparum* malaria (see [Basic Usage](basic.md)). However, it was coded flexibly and works with other organisms or amplicon panels (see [Advanced Usage](advanced.md)).

<br>

<p align="center">
<img src="img/home/nomadic_in_kenya.jpg" alt="example" class="bordered-img" width="75%">
</p>
<!-- *Nomadic* being used to process *P. falciparum* data in Kisian, Kenya. For more cool pictures, see our gallery. -->
<br>

## Features

- [x] Real-time read mapping with [_Minimap2_](https://github.com/lh3/minimap2).
- [x] Real-time sample quality control and amplicon coverage evaluation.
- [x] Real-time variant calling with [_delve_](https://github.com/berndbohmeier/delve) or [_bcftools_](https://github.com/samtools/bcftools).
- [x] Summarizing results with qc and prevalence calculations.
- [x] Support for different reference genomes or amplicons panels.

## Resources

- The NOMADS-MVP protocol is available in [English](https://www.protocols.io/view/nomads-mvp-rapid-genomic-surveillance-of-malaria-w-kxygxy284l8j/v1) and [French](https://www.protocols.io/view/surveillance-g-nomique-du-paludisme-par-la-m-thod-q26g75b5qlwz/v1).
- Read our paper [Mwenda, M., Mosler, K., Bohmeier, B. et al. Continental-scale genomic surveillance of Plasmodium falciparum malaria across sub-Saharan Africa with rapid nanopore sequencing. Nat Commun 17, 4218 (2026)](https://doi.org/10.1038/s41467-026-72358-z).


## Acknowledgements

This work was funded by the Bill and Melinda Gates Foundation (INV-003660, INV-048316).
