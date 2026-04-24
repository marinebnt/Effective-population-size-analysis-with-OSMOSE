library(osmose)
library(stringr)
library(reshape2)
library(cowplot)
library(patchwork)
library(dplyr)

source(paste0("osmose-ns-v4.4.1_clim/R/functions/functions_for_Ne.R"))
source(file.path("osmose-ns-v4.4.1_clim/R/functions/functions_Ne.R"))
source(file.path("osmose-ns-v4.4.1_clim/R/functions/save_plots_ne.R"))
source(file.path("osmose-ns-v4.4.1_clim/R/functions/functions_temporalNe.R"))
source(file.path("osmose-ns-v4.4.1_clim/R/functions/save_plots_ne_tempo.R"))
source(file.path("osmose-ns-v4.4.1_clim/R/functions/Data_abund_mat_by_age.R"))
source(file.path("osmose-ns-v4.4.1_clim/R/functions/function_INIT1_launch_simu_if_missing.R"))
source(file.path("osmose-ns-v4.4.1_clim/R/functions/function_INIT2_produceNe.R"))
source(file.path("osmose-ns-v4.4.1_clim/R/functions/ways_to_estimate_Nc.R"))
source(file.path("osmose-ns-v4.4.1_clim/R/functions/functions_Nb.R"))

# 0. OSMOSE configuration -------------------------------------------------

version    = "4.4.1"
javaPath   = "../osmose-private/inst/java"
jarFile    = file.path(javaPath, sprintf("osmose-%s-jar-with-dependencies.jar", version))

# 1. Model Configuration --------------------------------------------------
configDir  = "osmose-ns-v4.4.1_clim"
#main       = "New_param_Eco_Bioen_neutral.R"#_neutral.R"
configFile = file.path(configDir, main) # path to main configuration file
conf       = read_osmose(input=configFile)
seeding_y  = get_par(conf, "population.seeding.year.max")

# 2. Load config parameters  --------------------------------------------------

cat("extract configuration parameters\n")
dt       <- get_par(conf, "simulation.time.ndtPerYear")/get_par(conf, "output.recordfrequency.ndt")
ny       <- get_par(conf, "simulation.time.nyear")-get_par(conf, "output.start.year")
time_vec <- seq(0,ny, by=1/dt)[-1] # if year 0 is removed
years    <- time_vec[length(time_vec)]
nb_rep   <- get_par(conf, "simulation.nsimulation")
nb_sp    <- get_par(conf, "simulation.nspecies")
species  <- unlist(as.matrix(get_par(conf, "species.name.")[1:nb_sp]))
sp_names <- unlist(as.matrix(get_par(conf, "species.name.")[1:nb_sp]))
species_i <- rownames(as.matrix(get_par(conf, "species.name.")[1:nb_sp]))
sex_ratio <- get_par(conf, "sexratio")
ncohort   <- years
max_age   <- get_par(conf, "age.max")
dim_vec  <- c("abundanceByAge", "eggNumber", "mortalityByAge", "schoolNumber", "meanSizeByAge")
