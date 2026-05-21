
# osmose-ns-v4.4.1_clim/R/functions/dependencies.R
required_packages <- c("osmose", "stringr", "reshape2", "cowplot", "patchwork",
                       "minpack.lm", "ggplot2", "dplyr", "tidyr", "ncdf4", 
                       "gridExtra", "tidyverse", "matrixStats")

invisible(lapply(required_packages, library, character.only=T))


library(osmose)
library(stringr)
library(reshape2)
library(cowplot)
library(patchwork)
library(minpack.lm)
library(ggplot2)
library(dplyr)
library(tidyr)
library(ncdf4)
library(gridExtra)   # pour arranger plusieurs plots
# library(GGally)
library(tidyverse)
library(matrixStats)
library(ggpubr)

source(paste0("R/functions/functions_for_Ne.R"))
source(file.path("R/functions/functions_Ne.R"))
#source(file.path("osmose-ns-v4.4.1_clim/R/functions/save_plots_ne.R"))
source(file.path("R/functions/functions_temporalNe.R"))
#source(file.path("osmose-ns-v4.4.1_clim/R/functions/save_plots_ne_tempo.R"))
source(file.path("R/functions/Data_abund_mat_by_age.R"))
source(file.path("R/functions/function_INIT1_launch_simu_if_missing.R"))
source(file.path("R/functions/function_INIT2_produceNe.R"))
source(file.path("R/functions/ways_to_estimate_Nc.R"))
source(file.path("R/functions/functions_Nb2.R"))
source(file = "R/functions/auxiliary_for_lists.R")


# 0. OSMOSE configuration -------------------------------------------------

version    = "4.4.1"
javaPath   = "C:/Users/mbeneat/Documents/osmose/osmose-private/inst/java"
jarFile    = file.path(javaPath, sprintf("osmose-%s-jar-with-dependencies.jar", version))

# 1. Model Configuration --------------------------------------------------
configFile = file.path(main) # path to main configuration file
conf = read_osmose(input=configFile)
seeding_y <- get_par(conf, "population.seeding.year.max")

# 2. Load config parameters  --------------------------------------------------

cat("extract configuration parameters\n")
dt       <- get_par(conf, "simulation.time.ndtPerYear")/get_par(conf, "output.recordfrequency.ndt")
ny       <- get_par(conf, "simulation.time.nyear")-get_par(conf, "output.start.year")
time_vec <- seq(0,ny, by=1/dt)[-1] # if year 0 is removed
years    <- time_vec[length(time_vec)]
nb_rep   <- get_par(conf, "simulation.nsimulation")
nb_sp    <- get_par(conf, "simulation.nspecies") - 1 # remove shrimp
species  <- unlist(as.matrix(get_par(conf, "species.name.")[1:nb_sp]))
sp_names <- unlist(as.matrix(get_par(conf, "species.name.")[1:nb_sp]))
species_i <- rownames(as.matrix(get_par(conf, "species.name.")[1:nb_sp]))
sex_ratio <- get_par(conf, "sexratio")
ncohort   <- years
max_age   <- get_par(conf, "age.max")
dim_vec  <- c("abundanceByAge", "eggNumber", "mortalityByAge", "schoolNumber", "meanSizeByAge")

