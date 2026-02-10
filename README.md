# README

This repositroy contains files for BEAST analyses and plotting in Obeng et al. (2025).

## Steps
1. Prepare sequence files using `seq.R`.
2. Run BEAST analyses using the `.xml` files.
    + Each run with `beast2 -overwrite -seed 4321 <xml_file>`.
3. Plot results using `plots.R`.

## Notes
- `.xml` files for BEAST analyses.
    + `BD_averaged.{xml,log}`: Are v2.6 versions of the averaged $R_e$ analyses. TODO: update to v2.7.7 in `re_average.xml`

- `R` packages required: 
    + `tidyverse`, `latex2exp`, `beastio`, `readxl`, `patchwork`

## About analyses
Each analysis can be run with:
```bash
arr=(re_averaged.xml cluster_161.xml cluster_168.xml cluster_37.xml dpp.xml)
for xml in "${arr[@]}"; do
    beast2 -overwrite -seed 4321 -beagle $xml
done
```

- `BD_averaged.xml` / `re_average.xml`: Birth-death skyline analyses modelling $R_e$ over all clusters before PrEP, during the rollout, during lockdown, and post lockdown.

- `cluster_161.xml`, `cluster_168.xml`, `cluster_37.xml`: Birth-death skyline analyses with the same intervals as abbove for the top thre largest clusters.

- `dpp.xml` : Fits dirichlet process to all datasets similarly to the averaged $R_e$ analyes. TODO: fit Dirichlet process to each interval? (current barrier is .xml)