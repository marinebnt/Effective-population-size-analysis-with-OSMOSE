# Test multiple school size multipliers and then explore the effect on Ne & else
# 
# marine beneat
# 02/06/25
# 






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
configFile = file.path(configDir, main) # path to main configuration file
conf = read_osmose(input=configFile)
seeding_y <- get_par(conf, "population.seeding.year.max")
# 2. Loop over the multipliers  --------------------------------------------------

# multipliers = c("output_etaequalsone_bioen_neutral", 2, 3, 4, 5, 6)
multipliers = c(1,2,3,4,5,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9)
multipliers = c(1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9)

repo="Simu_MULTI"

for (m in multipliers){
  simulation = paste0(repo, "/output_m", m)
  outputDir  = file.path(configDir, "output", simulation)
  if(!dir.exists(outputDir)){
    cat(m, "\n")
    parameters = paste0("-Psimulation.nschool.multiplier=", m)
    options = paste0("-Xmx3g -Xms1g")
    run_osmose(input = configFile, output = outputDir, osmose = jarFile, version=version, options=options, parameters = parameters, force=T)
  }
  else{next}
}

# 3. Plot  --------------------------------------------------

source(paste0(configDir, "/R/functions/functions_for_Ne.R"))
source(file.path(configDir, "R/functions/functions_Ne.R"))
source(file.path(configDir, "R/functions/save_plots_ne.R"))
source(file.path(configDir, "R/functions/functions_temporalNe.R"))
source(file.path(configDir, "R/functions/save_plots_ne_tempo.R"))
source(file.path(configDir, "R/functions/Data_abund_mat_by_age.R"))

rep=1
multipliers = c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2,3,4,5)
incr_multiplier <- 0

var_to_plot <- c("mean_Ne", "mean_Ne_tempo", "mean_htz", "95quantile_Ne", "5quantile_Ne", 
                 "95quantile_Ne_tempo", "5quantile_Ne_tempo", "95quantile_htz","5quantile_htz")

var_list <- array(dim = c(length(var_to_plot), length(multipliers) , 16), 
                  dimnames = list(var_to_plot,
                                  multipliers,
                                  unlist(get_par(conf, "species.name"))[1:16]))



for (m in multipliers){
  k=0
  plot <- list()
  m_vec <- c()
  
  # initialisation
  simulation = paste0(repo, "/output_m", m)
  if (!(m %in% c(1,2,3,4,5,6,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9))){simulation = paste0(m)}
  cat("\n\t---> folder ", simulation, "\n")
  incr_multiplier <- incr_multiplier + 1
  ns1 = NULL
  Ne_version = paste0("School_multiplier_", m)
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
  sp_names <- unlist(as.matrix(get_par(conf, "species.name.")[1:nb_sp]))
  species_i <- rownames(as.matrix(get_par(conf, "species.name.")[1:nb_sp]))
  sex_ratio <- get_par(conf, "sexratio")
  ncohort   <- years
  max_age   <- get_par(conf, "age.max")
  
  # effective population sizes
  cat("\nprepare osmose data to fit the by-cohort requisites\n")
  data_Ne_raw <- estimate_for_Ne(outputDir_Simu=outputDir, outputDir_Ne = outputDir_Ne, dim_vec = dim_vec, species=species, nb_rep = nb_rep, OUTPUT_data=ns1)
  cat("producing Ne data")
  data_Ne <- estimate_Ne(outputDir_Ne=outputDir_Ne, data_Ne_raw, species=species, nb_rep, input_dir=outputDir)
  cat("producing temporal Ne data\n")
  tempoNe <- estimate_temporalNe(outputDir_Simu=outputDir, data_Ne, species, conf, seeding_y)
  
  for(sp in species){
    k=k+1
    cohort_time_corr <- data_Ne[["cohort_time_step_corresponding"]]
    co_to_include <- which(colSums(data_Ne[["longevity"]][[sp]][[paste0("rep_", rep)]] )>0)
    
    # Nc
    Nc <- as.data.frame(data_Ne[["abundanceByAge"]][[sp]])
    # Add time as a column
    Nc$age <- as.numeric(rownames(Nc))
    # Reshape data to long format
    df_nc_long_i <- pivot_longer(Nc, cols = -age, names_to = "Cohort", values_to = "Nc")
    
    # Ne
    Ne <- as.data.frame(data_Ne[["Ne"]][[sp]][[paste0("rep_", nb_rep)]])
    # Add time as a column
    Ne$age <- as.numeric(rownames(Ne))
    # Reshape data to long format
    df_ne_long_i <- pivot_longer(Ne, cols = -age, names_to = "Cohort", values_to = "Ne")
    # match the cohort age with the time step to which it corresponds
    for (co in 0:c(ncohort-1)){
      ID_co <- which(cohort_time_corr==co, arr.ind = T)[,1]
      df_ne_long_i[grep(paste0("\\bcohort_", co, "\\b"), df_ne_long_i$Cohort)[1:length(ID_co)], "Time"] <- ID_co
    }
    
    # merge Nc and Ne dataset
    data_tot <- full_join(df_nc_long_i, df_ne_long_i)
    for (t in min(data_tot$Time, na.rm = T):max(data_tot$Time, na.rm = T)){
      data <- data_tot[which(data_tot$Time == t),]
      data_tot[which(data_tot$Time == t),"Ne_w_mean"] <- weighted.mean(x=data$Ne, w=data$Nc, na.rm = T)
    }
    data_tot <- data_tot[which(data_tot$age == min(data_tot$age)),]
    
    # Ne tempo data
    df      <- data.frame(tempoNe[sp])
    df$Time <- as.numeric(rownames(df))
    colnames(df) <- c("Ne_tempo", "Time")
    
    # merge Ne, Nc, weighted Ne, tempo Ne
    data_tot <-full_join(data_tot, df)
    data_tot <- data_tot[-which(is.na(data_tot$Ne_tempo)),]
    
    multip <-  (max(data_tot$Ne_tempo, na.rm = T)/max(data_tot$Ne_w_mean, na.rm = T))
   
    if(is.nan(multip)){cat("next")
      plot[[k]] <- NULL
      m_vec <- c(m_vec, NA)
    next}
    
    if(multip == 0){cat("next")
      plot[[k]] <- NULL
      m_vec <- c(m_vec, NA)
    next}
    
    multip <-  (max(data_tot$Ne_tempo, na.rm = T)/max(data_tot$Ne_w_mean, na.rm = T))
    
    # Plot using ggplot2
    a <- ggplot(data_tot, aes(x = Time)) +
      geom_line(aes(y = Ne_w_mean)) +
      geom_line(aes(y = Ne_tempo/multip), col="orange", lty="longdash", size=0.75)+
      labs(title = paste0(sp, ": Comparing Functional traits Ne weighted by Nc\nand Temporal Ne"), x = "Time Step") +
      theme_minimal()+ 
      theme(legend.position="none")+
      scale_y_continuous(
        # Features of the first axis
        name = "weighted Functional traits Ne (black)",
        # Add a second axis and specify its features
        sec.axis = sec_axis(~.*multip, name="Temporal Ne (orange)")
      ) 
    par(mfrow = c(3,3))
    
    plot[[k]] <- a
    m_vec <- c(m_vec, multip)
  }
  
  # STORE THE PLOTS
  dir.create(paste0(configDir, "/output/Ne/PLOT_compareNe_and_Ne_temp/"))
  pdf(file = paste0(configDir, "/output/Ne/PLOT_compareNe_and_Ne_temp/BYSP", m, ".pdf"), width = 10)
  # pdf(file = paste0(paste0(pathNCDF, "_FINAL_plot/", "tot", ".pdf")), width=25, height=15)
  par(mfrow=c(1,1))
  for (k in 1:length(plot)) {
    multip <- m_vec[k]
    print(plot[[k]])
  }
  dev.off()
  
  
}

