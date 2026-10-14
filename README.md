# Effective (*N~e~*) and census (*N~T~*) population sizes responses to fishing in 15 North Sea fish species, simulated with OSMOSE.

#                Analysis and figure scripts.

### About

R scripts to analyse effective (*N~e~*) and census (*N~T~*) population sizes and the *N~e~/N~T~* ratio under varying fishing pressure for 15 North Sea fish species, using outputs of the multispecies eco-genetic model OSMOSE.

### Background

This repository contains the R scripts used to process OSMOSE simulation outputs and produce the analyses and figures of the manuscript **Genetic drift can be faster than the decline in abundance in harvested fish communities**. We examine how effective population size (*N~e~*), census population size (*N~T~*), the *N~e~/N~T~* ratio and the components of *N~e~* respond to a gradient of fishing pressure (*θ~F~*) across 15 North Sea fish species. *N~e~* is estimated with the Hill–Felsenstein framework, and its drivers are decomposed following Waples et al. (2011): generation length (*T*), variance in reproductive success (*V~k.~*) and longevity (*L~T~*).

### This folder is separated in two sections:

- The first section is about OSMOSE simulations and the version of the configurations used per simulation. This section also includes the extraction of the simultions outputs from the supercomputer we used, to produce light RData files, easily transferable.

  - This section includes the following files and folders:

    - `list_of_simu_F.txt` (list of the simulations launched, all the fishing scenarios tested).

    - `new_config_Fmax_Fstart` (folder with all the configuration files called by `list_of_simu_F.txt`).

    - `RESTART` (needed for OSMOSE simulations).

    - `R-datarmor` (folder with all the scripts used on the calculation machine).

    - `output/DATASETS` (folder with all the datasets extracted from the supercomputer - datasets were compressed to fit Github size requirements).

- The second section is about the analysis of the simulations

  - This section includes the following files and folders:

    - `DATA` (complementary information, not needed to have the scripts running).

    - `R` (scripts used to produce plots - if used, plots will be in output/PLOTS)

### Bibliography

Waples, R. S., Do, C., & Chopelet, J. (2011). Calculating *N* ~e~ and *N* ~e~ / *N* in age-structured populations: A hybrid Felsenstein-Hill approach. *Ecology*, *92*(7), 1513–1522. <https://doi.org/10.1890/10-1796.1>
