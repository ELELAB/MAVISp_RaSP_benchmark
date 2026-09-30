Cancer Systems Biology, Technical University of Denmark, 2800, Lyngby, Denmark
Cancer Structural Biology, Danish Cancer Institute, 2100, Copenhagen, Denmark
Repository associated to the publication:

> MAVISp: A Modular Structure-Based Framework for Genomic Variant Interpretation
> Matteo Arnaudi, Ludovica Beltrame, Kristine Degn, Mattia Utichi, Simone Scrima,
> Pablo Sanchez Izquierdo, Karolina Krzesinska, Francesca Maselli, Terezia Dorcakova,
> Jordan Safer, Alberte Heering Estad, Katrine Meldgard, Philipp Becker, Julie Bruun Brockhoff,
> Amalie Drud Nielsen, Valentina Sora, Alberto Pettenella, Jeremy Vinhas,
> Peter Wad Sackett, Claudia Cava, Anna Rohlin, Mef Nilbert, Sumaiya Iqbal, Matteo Lambrughi,
> Matteo Tiberti, Elena Papaleo. 
> bioRxiv https://doi.org/10.1101/2022.10.22.513328


## MAVISp_RaSP_benchmark

This repository contains our benchmarking of changes of folding free energy
upon mutation predicted by RaSP against those predicted by Rosetta using
data from MAVISp. This was done to assess whether we could confidently
replace costly Rosetta predictions with RaSP ones in the context of the
MAVISp framework.

These results are discussed in Supplementary Text S1 and its figures of the
MAVISp paper.

## Reproducing

### Requirements

Reproducing the results requires the `R` interpreter with some open source
R packages installed. We performed our calculations with the following setup:

- Ubuntu Linux 22.04, installed on a x86-64 server
- R 4.3.1
- the following R packages, that are explicitely loaded in the code:
  - tidyverse 2.0.0
  - tidymodels 1.1.1
  - corrr 0.4.4
  - forcats 1.0.0
  - ggnewscale 0.4.10

we also include a Conda environment definition file in the repository, which
defines the full environment in which we ran this code.

### Steps to reproduce

The steps that follow are designed to be performed on Linux or macOS,
using the terminal, with the `git` software installed

1. Create a local copy of this repository and enter its directory:
```
git clone https://github.com/ELELAB/MAVISp_RaSP_benchmark.git
cd MAVISp_RaSP_benchmark
```

2. (Optional) you can build a Conda environment using our provided definition,
to reproduce as closely as possible the environment we have used during
development. In order to do so, you need to have `anaconda` or `miniconda`
installed in your system and have the `conda` command available. To reproduce
our environment on a Linux system, run:

```
conda env create -p ./env -f environment.yaml
```

you will then need to activate the environment:

```
conda activate ./env
```

3. Run our R script to reproduce the analysis:

```
Rscript rasp_ros_bench.r rasp_benchmark_20122023/
```

