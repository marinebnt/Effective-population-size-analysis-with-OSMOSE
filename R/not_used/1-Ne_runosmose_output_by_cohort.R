#* marine beneat
#* 15/05/2025
#* Following Waples 2011
#* Idea : produce an output of the model giving the effective population size over time
#* Suggestion : run it in background and save it on the go, it is heavyier as you increase the species number

library(osmose)
library(stringr)
library(reshape2)

# 0. OSMOSE configuration -------------------------------------------------

version    = "4.4.1"
# javaPath = "C:/Users/cleaa/OSMOSE/osmose-private/inst/java"
javaPath   = "C:/Users/mbeneat/Documents/osmose/osmose-private/inst/java"
jarFile    = file.path(javaPath, sprintf("osmose-%s-jar-with-dependencies.jar", version))
# jarFile    = file.path("calibration_med05_test/.osmose.jar")

# 1. Model Configuration --------------------------------------------------

configDir  = "osmose-ns-v4.4.1_clim"
main       = "New_param_Eco_Bioen_neutral.R"#_neutral.R"
simulation = "output_etaequalsone_bioen_neutral"

configFile = file.path(configDir, main) # path to main configuration file
outputDir  = file.path(configDir, "output", simulation)

conf = read_osmose(input=configFile)

# 2. collect the model outputs configuration --------------------------------------------------

#update_osmose(input = configFile, osmose = jarFile, version = "4.4.1", verbose = TRUE)

options = "-Xmx3g -Xms1g"
# run_osmose(input = configFile, output = outputDir, osmose = jarFile, version=version, options=options)
# 
cat("==> The following steps will work, only if you set output.recordfrequency.ndt = 1\n")
cat("Need to activate the following outputs : meanSizeByAge, mortalityByAge, abundanceByAge, eggNumber")
ns1 = read_osmose(path = outputDir)

# 3. Transform the outputs into a list of the necessary traits --------------------------------------------------  
#[trait, sp, time, cohort, rep]
dt       <- get_par(conf, "simulation.time.ndtPerYear")/get_par(conf, "output.recordfrequency.ndt")
ny       <- get_par(conf, "simulation.time.nyear")-get_par(conf, "output.start.year")
time_vec <- seq(0,ny, by=1/dt)[-1] # if year 0 is removed
nb_rep   <- get_par(conf, "simulation.nsimulation")
nb_sp    <- get_par(conf, "simulation.nspecies")
sp_names <- unlist(as.matrix(get_par(conf, "species.name.")[1:nb_sp]))

# vector of the variables needed per cohort, per species, per replicate, to estimate the effective population size
dim_vec  <- c("abundanceByAge", "eggNumber", "mortalityByAge", "schoolNumber")

# creates a file stored in the outputDir_Ne folder, reorganizing the outputs of osmose with each column corresponding to a **cohort** (all individuals born at same time)
Ne_version = "restarting_everything"
source(paste0(configDir, "/R/functions/functions_for_Ne.R"))
estimate_for_Ne(outputDir_Simu=outputDir, outputDir_Ne=file.path(configDir, "output", "Ne", Ne_version), dim_vec = dim_vec, species=sp_names, nb_rep = nb_rep)


