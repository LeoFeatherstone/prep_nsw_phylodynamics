# NSW PrEP Phylodynamics

This repositroy contains files for BEAST analyses and plotting in Obeng et al. (2026).

## Steps
1. Prepare sequence files using `seq.R`.
2. Run BEAST analyses using the `.xml` files.
    + Each run with `beast2 -overwrite -seed 4321 <xml_file>`.
3. Plot results using `plots.R`.

## Notes
- `.xml` files for BEAST analyses.
    + `re_average.xml`: Birth-death skyline analyses modelling $R_e$ over all clusters before PrEP, during the rollout, during lockdown, and post lockdown
    + `cluster_{37,161,268}.xml`: Birth-Death skyline for top-3 largest clusters
    + `dpp.xml`: Again pooling clusters but with Dirchlet process prior on $R_e$

- `R` packages required: 
    + `tidyverse`, `latex2exp`, `beastio`, `readxl`, `patchwork`, `knitr`

## About analyses
Each analysis can be run with, or with `run_beast.sh`:
```bash
arr=(re_averaged.xml cluster_161.xml cluster_168.xml cluster_37.xml dpp.xml)
for xml in "${arr[@]}"; do
    beast2 -overwrite -seed 4321 -beagle $xml
done
```