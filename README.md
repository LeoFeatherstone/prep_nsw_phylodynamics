# NSW PrEP Phylodynamics

This repository contains files for BEAST analyses and plotting in Obeng et al. (2026).

## Directory Structure
- `sequences/` - FASTA sequence files
- `xml/` - BEAST XML configuration files
- `scripts/` - R scripts for data processing and plotting
- `results/` - Output from analyses
  - `tables/` - Summary statistics
  - `figures/` - Generated plots

## Steps
1. Prepare sequence files using `scripts/sequences.R`.
2. Run BEAST analyses using the XML files in `xml/`.
    + Run with `./run_beast.sh` or individually: `beast2 -overwrite -seed 4321 xml/<xml_file>`
3. Plots & tables using `scripts/results.R`.

## Notes
- `xml/re_averaged.xml`: Birth-death skyline analyses modelling $R_e$ over all clusters
- `xml/dpp.xml`: Pooled clusters with Dirichlet process prior on $R_e$

- R packages required:
    + `tidyverse`, `latex2exp`, `beastio`, `readxl`, `patchwork`, `knitr`