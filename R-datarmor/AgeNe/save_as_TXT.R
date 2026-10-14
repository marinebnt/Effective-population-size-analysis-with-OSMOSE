# create a .txt file per cohort, as input for AgeNe
# 05/09/25
# MB




library(ggplot2)
library(osmose)
library(stringr)
library(reshape2)
library(cowplot)
library(patchwork)

# 0. OSMOSE configuration -------------------------------------------------

version    = "4.4.1"
# javaPath = "C:/Users/cleaa/OSMOSE/osmose-private/inst/java"
javaPath   = "C:/Users/mbeneat/Documents/osmose/osmose-private/inst/java"
jarFile    = file.path(javaPath, sprintf("osmose-%s-jar-with-dependencies.jar", version))
# jarFile    = file.path("calibration_med05_test/.osmose.jar")

# 1. Model Configuration --------------------------------------------------

configDir  = "osmose-ns-v4.4.1_clim"
main       = "New_param_Eco_Bioen_neutral.R"#_neutral.R"
configFile = file.path(configDir, main) # path to main configuration file
conf = read_osmose(input=configFile)
seeding_y <- get_par(conf, "population.seeding.year.max")

source(paste0(configDir, "/R/functions/functions_for_Ne.R"))
source(file.path(configDir, "R/functions/functions_Ne.R"))
source(file.path(configDir, "R/functions/save_plots_ne.R"))
source(file.path(configDir, "R/functions/functions_temporalNe.R"))
source(file.path(configDir, "R/functions/save_plots_ne_tempo.R"))
source(file.path(configDir, "R/functions/Data_abund_mat_by_age.R"))


# 2. Loop over the multipliers  --------------------------------------------------

# multipliers = c("output_etaequalsone_bioen_neutral", 2, 3, 4, 5, 6)
multipliers = c(1,2,3,4,5,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9)
multipliers = c(1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9)


for (m in multipliers){
  simulation = paste0("output_m", m)
  outputDir  = file.path(configDir, "output", simulation)
  if(!dir.exists(outputDir)){
    cat(m, "\n")
    parameters = paste0("-Psimulation.nschool.multiplier=", m)
    options = paste0("-Xmx3g -Xms1g")
    # run_osmose(input = configFile, output = outputDir, osmose = jarFile, version=version, options=options, parameters = parameters, force=T)
  }
  else{next}
}



# 3. Launch loop to have trait variables per cohort + Save as txt  --------------------------------------------------
#* Rules : 
# LINE 1
# The first line is free, so include comments or useful information here. This information will appear in the output.
#*
# LINE 2
# Haploid model (see file „haploid.txt‟ for example)
# Requires 2 positive integers: MaxAge (age in years of oldest individuals), and N1 (number of newborns each year).
#*
# LINE 3 and + 
# each line is data for one age, up to MaxAge years. So, the total number of input lines for each population is 2+MaxAge.
# Haploid model
# The first column for each of these lines is the age, followed by age-specific survival rate (sx), birth rate (bx), and Poisson factor.
# A sample line (for age 3) looks like this: 3 0.5 2 1
#* 




rep=1
incr_multiplier <- 0
multipliers=1
m=1
var_to_store <- c("mean_Ne", "mean_Ne_tempo", "mean_htz", "quantile95_Ne", "quantile5_Ne", 
                 "quantile95_Ne_tempo", "quantile5_Ne_tempo", "quantile95_htz","quantile5_htz")

for (m in multipliers){
  
  # initialisation
  simulation = paste0("output_m", m)
  Ne_version = paste0("School_multiplier_", m)
  if (!(m %in% c(1,2,3,4,5,6,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9))){simulation = paste0(m)}
  cat("\n\t---> folder ", simulation, "\n")
  incr_multiplier <- incr_multiplier + 1
  ns1 = NULL
  outputDir  = file.path(configDir, "output", simulation)
  outputDir_Ne =file.path(configDir, "output", "Ne","by_multiplier", Ne_version)
  dir.create(outputDir_Ne)
  dim_vec  <- c("abundanceByAge", "eggNumber", "mortalityByAge", "schoolNumber")
  ns1 = read_osmose(path = outputDir)
  
  cat("extract configuration parameters\n")
  dt       <- get_par(conf, "simulation.time.ndtPerYear")/get_par(conf, "output.recordfrequency.ndt")
  ny       <- get_par(conf, "simulation.time.nyear")-get_par(conf, "output.start.year")
  time_vec <- seq(0,ny, by=1/dt)[-1] # if year 0 is removed
  years    <- time_vec[length(time_vec)]
  nb_rep   <- get_par(conf, "simulation.nsimulation")
  nb_sp    <- get_par(conf, "simulation.nspecies")
  species <- unlist(as.matrix(get_par(conf, "species.name.")[1:nb_sp]))
  species_i <- rownames(as.matrix(get_par(conf, "species.name.")[1:nb_sp]))
  sex_ratio <- get_par(conf, "sexratio")
  ncohort   <- years
  
  # effective population sizes
  cat("prepare osmose data to fit the by-cohort requisites\n")
  # if (!file.exists(paste0(outputDir_Ne, "/Data_for_Ne.RData"))){
  cat("producing by cohort raw data for Ne and temporal Ne estimations\n")
  estimate_for_Ne(OUTPUT = ns1, outputDir_Simu=outputDir, outputDir_Ne = outputDir_Ne, dim_vec = dim_vec, species=species, nb_rep = nb_rep)
  # }  else {cat("=> loading already existing by-cohort info\n")}
  data_Ne_raw <- get(load(paste0(outputDir_Ne, "/Data_for_Ne.RData")))
  # if (!file.exists(paste0(outputDir_Ne, "/Ne.RData"))){
  cat("producing Ne data\n")
  data_Ne <-  SimDesign::quiet(estimate_Ne(input_dir = outputDir, outputDir_Ne = outputDir_Ne, data_Ne = data_Ne_raw, species = species, nb_rep = nb_rep))
  
  for (sp in species){
    for (coho in 1:ncohort){
      #LINE 1
      title <- c(sp, ", cohort ",coho)
      
      #LINE 2 
      age_max <- as.numeric(which(data_Ne$abundanceByAge[[sp]][,coho]==0)[1])-1
      N1 <- round(data_Ne$abundanceByAge[[sp]][1,coho]*10^-3)
        
      if(age_max == 0){next}
      
      #LINE 3 : age, followed by age-specific survival rate (sx), birth rate (bx), and Poisson factor.
      age <- 1:age_max
      sx  <- as.numeric(data_Ne$survival_rate[[sp]]$rep_1[,coho][-as.numeric(which(data_Ne$abundanceByAge[[sp]][,coho]==0))])
      bx  <- data_Ne$b_rate[[sp]]$rep_1[,coho][-as.numeric(which(data_Ne$abundanceByAge[[sp]][,coho]==0))]
      P   <- rep(1, length(age))
      
      Tab <- data.frame(age, sx, bx, P) # line 3 and  +
      Tab <- rbind(c(title, ""), c(age_max, N1, "", ""), Tab) # line 2
      colnames(Tab) <- NULL
      
      # write .txt file in a dedicated repo
      dir_out <- paste0(configDir, "/output/TXT_life_table/m_", m, "/", sp)
      dir.create(dir_out, recursive = T)
      write.table(x = Tab,
                  file = paste0(dir_out, "/", coho, ".txt"),  # Save the file as pirates.txt
                  sep = "\t",
                  row.names = F,
                  quote=FALSE) 
    }
  }
}





