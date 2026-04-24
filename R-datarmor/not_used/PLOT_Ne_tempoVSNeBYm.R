


source(file.path("osmose-ns-v4.4.1_clim/R/functions/source.R"))
repo        = "Simu_MULTI"
multip_name = "m"
multipliers = na.omit(as.numeric(stringr::str_remove(list.files("osmose-ns-v4.4.1_clim/output/Simu_MULTI"), "output_m")))


# function to prepare environment : launch simulation if repo missing
function_INIT1(MULTIP=multipliers, REPO=repo)


# prepare environment for this particular plot
rep=1
incr_multiplier <- 0
k <- 0
dt_Ne_tempo <- c(20)


multi_vec <-c()
Nc_Ne <- list()
for (sp in species){
  k=0
  for (m in multipliers){
    k=k+1
    
    # initialisation
    simulation = paste0("Simu_MULTI/output_m", m)
    Ne_version = paste0("School_multiplier_", m)
    if (!(m %in% c(1,2,3,4,5,6,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9))){simulation = paste0(m)}
    cat("\n\t---> folder ", simulation, "\n")
    incr_multiplier <- incr_multiplier + 1
    ns1 = NULL
    outputDir_Simu  = file.path(configDir, "output", simulation)
    outputDir_Ne    = file.path(configDir, "output", "Ne","by_multiplier", Ne_version)
    dir.create(outputDir_Ne)
    dim_vec  <- c("abundanceByAge", "eggNumber", "mortalityByAge", "schoolNumber")
    ns1 = read_osmose(path = outputDir_Simu)
    
    # function to prepare environment : 
    INIT2 <- function_INIT2(outputDir_Simu, outputDir_Ne, dim_vec, species, nb_rep, OUTPUT_data=ns1 , conf, seeding_y, REPO=repo)
    data_Ne_raw  <- INIT2[[1]]
    data_Ne      <- INIT2[[2]]
    tempoNe      <- INIT2[[3]]
  
    cat(sp, "\n")
    cohort_time_corr <- data_Ne[["cohort_time_step_corresponding"]]
    co_to_include    <- which(colSums(data_Ne[["longevity"]][[sp]][[paste0("rep_", rep)]] )>0)

    # Ne data
    df <- as.data.frame(data_Ne[["Ne"]][[sp]][[paste0("rep_", rep)]])
    # Add time as a column
    df$age <- as.numeric(rownames(df))
    # Reshape data to long format
    df_long_Ne <- pivot_longer(df, cols = -age, names_to = "Cohort", values_to = "Ne")

    # match the cohort age with the time step to which it corresponds
    for (co in 0:c(ncohort-3)){ # removing the last 3 cohorts that are incomplete
      ID_co <- which(cohort_time_corr==co, arr.ind = T)[,1]
      df_long_Ne[grep(paste0("\\bcohort_", co, "\\b"), df_long_Ne$Cohort)[1:length(ID_co)], "Time"] <- ID_co
    }

    # clean df_long_Ne
    df_long_Ne <- df_long_Ne[which(df_long_Ne$age==1),]

    # Ne tempo dataset
    df      <- as.data.frame(tempoNe[sp])
    df$Time <- as.numeric(rownames(df))
    colnames(df) <- c("Ne_tempo", "Time")

    # bind datasets
    df_long_tot <- dplyr::left_join(df_long_Ne, df)
    df_long_tot <- df_long_tot[-which(is.na(df_long_tot$Ne_tempo)),]
    if(any(which(df_long_tot$Ne == 0))){
      df_long_tot[-which(df_long_tot$Ne == 0),] -> df_long_tot
    }
    if(any(is.na(df_long_tot$Ne))){
      df_long_tot[-which(is.na(df_long_tot$Ne)), ] -> df_long_tot
    }
    if(!length(df_long_tot$Time > seeding_y)>0){
      cat("no data in the required time and age categories\n")
      next
    }
    multip <- (max(df_long_tot$Ne_tempo, na.rm = T)/max(df_long_tot$Ne, na.rm = T))

    a <- ggplot(df_long_tot, aes(x = Time)) +
      geom_point(aes(y = Ne), size=3) +
      geom_smooth(method = loess, aes(y = Ne), lwd=0.7)+
      geom_point(aes(y = Ne_tempo/(multip)), size=2, col="red") +
      geom_smooth(method = loess, aes(y=Ne_tempo/(multip)), lwd=0.7, col="red")+
      labs(title = paste0(sp, ": Comparing Ne and Ne_tempo,\n multiplier = ",m), x = "Time Step", y = "Ne") +
      theme_minimal()+
      theme(legend.position="none")+
      scale_y_continuous(
        # Features of the first axis
        name = "Ne (blue/black)",
        # Add a second axis and specify its features
        sec.axis = sec_axis(~.*multip, name="Ne_tempo (red)")
      )
    multi_vec <- c(multi_vec, multip)
    Nc_Ne[[k]] <- a
  }
  # STORE THE PLOTS
  dir.create(paste0(configDir, "/output/Ne/PLOT_compareNe_and_Ne_temp/"))
  pdf(file = paste0(configDir, "/output/Ne/PLOT_compareNe_and_Ne_temp/BYm", sp,".pdf"), width = 10)
  # pdf(file = paste0(paste0(pathNCDF, "_FINAL_plot/", "tot", ".pdf")), width=25, height=15)
  par(mfrow=c(1,1))
  for (j in 1:length(Nc_Ne)) {
    multip <- multi_vec[j]
    print(Nc_Ne[j])
  }
  dev.off()
}

# 
# 
# library(osmose)
# library(stringr)
# library(reshape2)
# 
# # 0. OSMOSE configuration -------------------------------------------------
# 
# version    = "4.4.1"
# # javaPath = "C:/Users/cleaa/OSMOSE/osmose-private/inst/java"
# javaPath   = "C:/Users/mbeneat/Documents/osmose/osmose-private/inst/java"
# jarFile    = file.path(javaPath, sprintf("osmose-%s-jar-with-dependencies.jar", version))
# # jarFile    = file.path("calibration_med05_test/.osmose.jar")
# 
# # 1. Model Configuration --------------------------------------------------
# 
# configDir  = "osmose-ns-v4.4.1_clim"
# main       = "New_param_Eco_Bioen_neutral.R"#_neutral.R"
# configFile = file.path(configDir, main) # path to main configuration file
# conf = read_osmose(input=configFile)
# seeding_y <- get_par(conf, "population.seeding.year.max")
# # 2. Loop over the multipliers  --------------------------------------------------
# 
# # multipliers = c("output_etaequalsone_bioen_neutral", 2, 3, 4, 5, 6)
# multipliers = c(1,2,3,4,5,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9)
# multipliers = c(1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9)
# 
# for (m in multipliers){
#   simulation = paste0("output_m", m)
#   outputDir  = file.path(configDir, "output", simulation)
#   if(!dir.exists(outputDir)){
#     cat(m, "\n")
#     parameters = paste0("-Psimulation.nschool.multiplier=", m)
#     options = paste0("-Xmx3g -Xms1g")
#     # run_osmose(input = configFile, output = outputDir, osmose = jarFile, version=version, options=options, parameters = parameters, force=T)
#   }
#   else{next}
# }
# 
# # 3. Plot  --------------------------------------------------
# 
# source(paste0(configDir, "/R/functions/functions_for_Ne.R"))
# source(file.path(configDir, "R/functions/functions_Ne.R"))
# source(file.path(configDir, "R/functions/save_plots_ne.R"))
# source(file.path(configDir, "R/functions/functions_temporalNe.R"))
# source(file.path(configDir, "R/functions/save_plots_ne_tempo.R"))
# 
# rep=1
# multipliers = c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,2,3,4,5)
# multipliers = c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.1,2,3,4)
# incr_multiplier <- 0
# var_list <- array(dim = c(3, length(multipliers), 16), dimnames = list(c("var_Ne_max", "var_Ne_mean", "var_htz"),
#                                                                        multipliers,
#                                                                        unlist(get_par(conf, "species.name"))[1:16]))
# 
# multi_vec <-c()
# Nc_Ne <- list()
# species <- c("AtlanticHerring",  "AtlanticMackerel", "Sandeel",          "EuropeanSprat",    "NorwayPout",       
#              "EuropeanPlaice"  ,"CommonSole" ,      "Saithe"     ,      "AtlanticCod"  ,    "Haddock"  ,
#              "HorseMackerel" ,   "Whiting"  ,"CommonDab" ,       "GreyGurnard"   ,   "Hake" ,            "Shrimp" )         
# for (sp in species[1]){
#   k=0
#   for (m in multipliers){
#     k=k+1
#     # initialisation
#     simulation = paste0("output_m", m)
#     if (!(m %in% c(1,2,3,4,5,6,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9))){simulation = paste0(m)}
#     cat("---> folder ", simulation, "\n")
#     incr_multiplier <- incr_multiplier + 1
#     ns1 = NULL
#     Ne_version = paste0("School_multiplier_", m)
#     outputDir  = file.path(configDir, "output", simulation)
#     outputDir_Ne =file.path(configDir, "output", "Ne","by_multiplier", Ne_version)
#     dir.create(outputDir_Ne)
#     dim_vec  <- c("abundanceByAge", "eggNumber", "mortalityByAge", "schoolNumber")
#     ns1 = read_osmose(path = outputDir)
#     cat("extract configuration parameters\n")
#     dt       <- get_par(conf, "simulation.time.ndtPerYear")/get_par(conf, "output.recordfrequency.ndt")
#     ny       <- get_par(conf, "simulation.time.nyear")-get_par(conf, "output.start.year")
#     time_vec <- seq(0,ny, by=1/dt)[-1] # if year 0 is removed
#     years    <- time_vec[length(time_vec)]
#     nb_rep   <- get_par(conf, "simulation.nsimulation")
#     nb_sp    <- get_par(conf, "simulation.nspecies")
#     species <- unlist(as.matrix(get_par(conf, "species.name.")[1:nb_sp]))
#     species_i <- rownames(as.matrix(get_par(conf, "species.name.")[1:nb_sp]))
#     sex_ratio <- get_par(conf, "sexratio")
#     ncohort   <- years
#     
#     # effective population sizes
#     cat("prepare osmose data to fit the by-cohort requisites\n")
#     # if (!file.exists(paste0(outputDir_Ne, "/Data_for_Ne.RData"))){
#       cat("producing raw by cohort data for Ne and temporal Ne estimation\n")
#       estimate_for_Ne(outputDir_Simu=outputDir, outputDir_Ne = outputDir_Ne, dim_vec = dim_vec, species=species, nb_rep = nb_rep)
#     # }  else {cat("=> loading already existing by-cohort info\n")}
#     data_Ne_raw <- get(load(paste0(outputDir_Ne, "/Data_for_Ne.RData")))
#     # if (!file.exists(paste0(outputDir_Ne, "/Ne.RData"))){
#       cat("producing Ne data\n")
#       estimate_Ne(data_Ne=data_Ne_raw, species, nb_rep)
#     # }  else {cat("=> loading already existing Ne info\n")}
#     data_Ne <- get(load(paste0(outputDir_Ne, "/Ne.RData")))
#     # if (!file.exists(paste0(outputDir_Ne, "/temporalNe.RData"))){
#       cat("producing temporal Ne data\n")
#       estimate_temporalNe(outputDir_Simu=outputDir, data_Ne, species, conf, seeding_y)
#     # }  else {cat("=> loading already existing temporal Ne info\n")}
#     tempoNe <- get(load(paste0(outputDir_Ne, "/temporalNe.RData")))
#     
#     cat(sp, "\n")
#     cohort_time_corr <- data_Ne[["cohort_time_step_corresponding"]]
#     co_to_include <- which(colSums(data_Ne[["longevity"]][[sp]][[paste0("rep_", rep)]] )>0)
#     
#     # Ne data
#     df <- as.data.frame(data_Ne[["Ne"]][[sp]][[paste0("rep_", rep)]])
#     # Add time as a column
#     df$age <- as.numeric(rownames(df))
#     # Reshape data to long format
#     df_long_Ne <- pivot_longer(df, cols = -age, names_to = "Cohort", values_to = "Ne")
#     
#     # match the cohort age with the time step to which it corresponds
#     for (co in 0:c(ncohort-3)){ # removing the last 3 cohorts that are incomplete
#       ID_co <- which(cohort_time_corr==co, arr.ind = T)[,1]
#       df_long_Ne[grep(paste0("\\bcohort_", co, "\\b"), df_long_Ne$Cohort)[1:length(ID_co)], "Time"] <- ID_co
#     }
#     
#     # clean df_long_Ne
#     df_long_Ne <- df_long_Ne[which(df_long_Ne$age==1),]
#     
#     # Ne tempo dataset
#     df      <- as.data.frame(tempoNe[sp])
#     df$Time <- as.numeric(rownames(df))
#     colnames(df) <- c("Ne_tempo", "Time")
#     
#     # bind datasets
#     df_long_tot <- dplyr::left_join(df_long_Ne, df)
#     df_long_tot <- df_long_tot[-which(is.na(df_long_tot$Ne_tempo)),]
#     if(any(which(df_long_tot$Ne == 0))){
#       df_long_tot[-which(df_long_tot$Ne == 0),] -> df_long_tot
#     }
#     if(any(is.na(df_long_tot$Ne))){
#       df_long_tot[-which(is.na(df_long_tot$Ne)), ] -> df_long_tot
#     }
#     if(!length(df_long_tot$Time > seeding_y)>0){
#       cat("no data in the required time and age categories\n")
#       next
#     }
#     multip <- (max(df_long_tot$Ne_tempo, na.rm = T)/max(df_long_tot$Ne, na.rm = T))
#     
#     a <- ggplot(df_long_tot, aes(x = Time)) +
#       geom_point(aes(y = Ne), size=3) +
#       geom_smooth(method = loess, aes(y = Ne), lwd=0.7)+
#       geom_point(aes(y = Ne_tempo/(multip)), size=2, col="red") +
#       geom_smooth(method = loess, aes(y=Ne_tempo/(multip)), lwd=0.7, col="red")+
#       labs(title = paste0(sp, ": Comparing Ne and Ne_tempo,\n multiplier = ",m), x = "Time Step", y = "Ne") +
#       theme_minimal()+ 
#       theme(legend.position="none")+
#       scale_y_continuous(
#         # Features of the first axis
#         name = "Ne (blue/black)",
#         # Add a second axis and specify its features
#         sec.axis = sec_axis(~.*multip, name="Ne_tempo (red)")
#       ) 
#     multi_vec <- c(multi_vec, multip)
#     Nc_Ne[[k]] <- a
#   }
#   # STORE THE PLOTS
#   dir.create(paste0(configDir, "/output/Ne/PLOT_compareNe_and_Ne_temp/"))
#   pdf(file = paste0(configDir, "/output/Ne/PLOT_compareNe_and_Ne_temp/BYm", sp,".pdf"), width = 10)
#   # pdf(file = paste0(paste0(pathNCDF, "_FINAL_plot/", "tot", ".pdf")), width=25, height=15)
#   par(mfrow=c(1,1))
#   for (j in 1:length(Nc_Ne)) {
#     multip <- multi_vec[j]
#     print(Nc_Ne[j])
#   }
#   dev.off()
# }
# 
# 
# 
# 
