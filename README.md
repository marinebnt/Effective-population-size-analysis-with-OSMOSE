# Effective and census population sizes responses to fishing in 15 North Sea fish species, simulated with OSMOSE.

#                Simulations and Analysis

### About

R scripts to analyse effective (*N~e~*) and census (*N~T~*) population sizes and the *N~e~/N~T~* ratio under varying fishing pressure for 15 North Sea fish species, using outputs of the multispecies eco-genetic model OSMOSE.

### Background

This repository contains the R scripts used to launch OSMOSE simulation. It also includes the scripts used to process the simulations' outputs of the manuscript **Genetic drift can increase faster than abundance declines in harvested fish: insights from a multispecies ecosystem model**. We examine how effective population size (*N~e~*), census population size (*N~T~*), the *N~e~/N~T~* ratio and the components of *N~e~* respond to a gradient of fishing pressure (*θ~F~*) across 15 North Sea fish species. *N~e~* is estimated with the Hill–Felsenstein framework, and its drivers are decomposed following Waples et al. (2011): generation length (*T*), variance in reproductive success (*V~k.~*) and longevity (*L~T~*).

### This folder is separated in two sections:

- The first section is about OSMOSE simulations and the configuration versions used for each simulation. For each species, we simulated the change in its fishing coefficient *θ~F~*, and replicated the simulation 20 times. This section also includes the extraction of the simultions outputs from the supercomputer we used, to produce light RData files, easily transferable:
 
    - `R-datarmor/prepare_scenarios_for_simulations` (from the initial configuration: writes the configuration including different fishing scenarios, stored in `new_config_Fmax_Fstart`).

    - `new_config_Fmax_Fstart` (folder with all the configuration files called by `list_of_simu_F.txt`).

    - `list_of_simu_F.txt` (list of the simulations launched, for all the fishing scenarios tested).

    - `RESTART` (needed for OSMOSE simulations, simulations’ spin-up lasting 80 years, and used as the initial condition for five 120-year-long fishing scenarios (see Morell et al. (2026) for details)).

    - `R-datarmor` (folder with all the scripts used on the calculation machine, which also calls functions from `R/functions`).

    - `output/DATASETS` (folder with all the datasets extracted from the supercomputer - datasets were compressed to fit Github size requirements).

- The second section is about the analysis of the simulations:

    - `R` (scripts used to produce figures - if used, plots will be stored in output/PLOTS)

    - `DATA` (complementary information, not needed to have the scripts running).

### Bibliography
Morell, A., Shin, Y.-J., Barrier, N., Travers-Trolet, M., Halouani, G., & Ernande, B. (2023). Bioen-OSMOSE: A bioenergetic marine ecosystem model with physiological response to temperature and oxygen. Progress in Oceanography, 216, 103064. https://doi.org/10.1016/j.pocean.2023.103064 
Waples, R. S., Do, C., & Chopelet, J. (2011). Calculating *N* ~e~ and *N* ~e~ / *N* in age-structured populations: A hybrid Felsenstein-Hill approach. *Ecology*, *92*(7), 1513–1522. <https://doi.org/10.1890/10-1796.1>
