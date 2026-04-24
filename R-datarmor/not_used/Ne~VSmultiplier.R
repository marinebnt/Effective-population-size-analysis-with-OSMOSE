# Test Simu_MULTIple school size multipliers and then explore the effect on Ne & else
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
# multipliers = c(1,2,3,4,5,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9)
multipliers = c(1,0.1,2,3,4,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9)

multipliers = c(0.5,1,1.5,2,3,4,5)
for (m in multipliers[1]){
  simulation = paste0("Simu_MULTI/TESToutput_m", m)
  outputDir  = file.path(configDir, "output", simulation)
  if(!dir.exists(outputDir)){
    cat(m, "\n")
    parameters = paste0("-Psimulation.nschool.Simu_MULTIplier=", m)
    options = paste0("-Xmx3g -Xms1g")
    
    
    
    
    
    # In your R script before launching OSMOSE
    cat("Checking configuration parameters:\n")
    cat("output.mature.abundance.byage.enabled =", 
        get_par(conf, "output.mature.abundance.byage.enabled"), "\n")
    cat("output.start.year =", get_par(conf, "output.start.year"), "\n")
    cat("Age distribution min/max:", 
        get_par(conf, "output.distrib.byAge.min"), "/",
        get_par(conf, "output.distrib.byAge.max"), "\n")
    
    
    run_osmose(input = configFile, output = outputDir, osmose = jarFile, version=version, options=options, parameters = parameters, force=T)
  }  else{next}
}

# 3. Plot  --------------------------------------------------

source(paste0(configDir, "/R/functions/functions_for_Ne.R"))
source(file.path(configDir, "R/functions/functions_Ne.R"))
source(file.path(configDir, "R/functions/save_plots_ne.R"))
source(file.path(configDir, "R/functions/functions_temporalNe.R"))
source(file.path(configDir, "R/functions/save_plots_ne_tempo.R"))
source(file.path(configDir, "R/functions/Data_abund_mat_by_age.R"))

rep=1
# multipliers = c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2,3,4,5)
incr_Simu_MULTIplier <- 0
var_list <- array(dim = c(3, length(multipliers), 16), dimnames = list(c("var_Ne_max", "var_Ne_mean", "var_htz"),
                                                                       multipliers,
                                                                       unlist(get_par(conf, "species.name"))[1:16]))




# multipliers = "output_etaequalsone_bioen_neutral"
for (m in multipliers[2]){

  # initialisation
  simulation = paste0("Simu_MULTI/output_m", m)
  if (!(m %in% c(1,2,3,4,5,6,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9))){simulation = paste0(m)}
  cat("\n\t---> folder ", simulation, "\n")
  incr_Simu_MULTIplier <- incr_Simu_MULTIplier + 1
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
  sp_names <- unlist(as.matrix(get_par(conf, "species.name.")[1:nb_sp]))
  species  <- unlist(as.matrix(get_par(conf, "species.name.")[1:nb_sp]))
  species_i <- rownames(as.matrix(get_par(conf, "species.name.")[1:nb_sp]))
  sex_ratio <- get_par(conf, "sexratio")
  ncohort   <- years
  max_age   <- get_par(conf, "age.max")

  # effective population sizes
  # cat("prepare osmose data to fit the by-cohort requisites\n")
  # if (!file.exists(paste0(outputDir_Ne, "/Data_for_Ne.RData"))){
    cat("producing raw by cohort data for Ne and temporal Ne estimation\n")
    estimate_for_Ne(outputDir_Simu=outputDir, outputDir_Ne = outputDir_Ne, dim_vec = dim_vec, species=species, nb_rep = nb_rep, OUTPUT=ns1)
  # }  else {cat("=> loading already existing by-cohort info\n")}
  data_Ne_raw <- get(load(paste0(outputDir_Ne, "/Data_for_Ne.RData")))
  
  # if (!file.exists(paste0(outputDir_Ne, "/Ne.RData"))){
    cat("producing Ne data\n")
    data_Ne <-estimate_Ne(data_Ne_shape = data_Ne_raw, outputDir_Ne = outputDir_Ne, input_dir=outputDir , species = sp, nb_rep = nb_rep)
  # }  else {cat("=> loading already existing Ne info\n")}
  data_Ne <- get(load(paste0(outputDir_Ne, "/Ne.RData")))
  # if (!file.exists(paste0(outputDir_Ne, "/temporalNe.RData"))){
    cat("producing temporal Ne data\n")
    estimate_temporalNe(outputDir_Simu=outputDir, data_Ne, species, conf, seeding_y)
  # }  else {cat("=> loading already existing temporal Ne info\n")}
  tempoNe <- get(load(paste0(outputDir_Ne, "/temporalNe.RData")))
  
  # #plots
  # cat("save plots\n")
  # SimDesign::quiet(save_plots_ne(data_Ne, config=ns1, outputDir_Simu=paste0(outputDir), 
  #                                outputDir_Ne, seeding_y))
  # 
  # SimDesign::quiet(save_plots_ne_tempo(data_Ne_t = tempoNe, data_Ne = data_Ne, config=ns1, outputDir_Simu=paste0(outputDir), 
  #                                paste0(outputDir_Ne), seeding_y))
  # 
  # cat("save variability Ne per species per Simu_MULTIplier, to plot it after (only for rep1)\n")
  # rep_fix=1
  # for(sp in species){
  #   # VAR NE
  #   Ne <- data_Ne[["Ne"]][[sp]][[paste0("rep_", rep_fix)]]
  #   N  <- data_Ne[["abund_scaled"]][[sp]][[paste0("rep_", rep_fix)]]
  #   variability_over_time <- apply(Ne, 1, sd)
  #   variability_over_time <- apply(Ne, 1, sd)
  #   mx <- max(variability_over_time)
  #   mn <- mean(variability_over_time, na.rm = T)
  #   
  #   # VAR HTZ
  #   if (!file.exists(paste0(outputDir_Ne, "/htz", sp, ".RData"))){
  #     cat("saving htz data\n")
  #     data <- ncdf4::nc_open(paste0(outputDir, "/Genetic/ns_observedHtz-", sp, "Simu0.nc"))
  #     data_nc <- ncdf4::ncvar_get(data, "observed_heterozygosity")
  #     trait1_data <- data_nc[, 1, ]
  #     ncdf4::nc_close(data)
  #     mean_htz <- colMeans(trait1_data)
  #     time <- 0:c(dim(data_nc)[3]-1)
  #     data_htz <- data.frame(Time=time, HTZ=mean_htz)
  #     save(data_htz, file=paste0(outputDir_Ne, "/htz", sp, ".RData"))
  #   }  else {
  #     data_htz <- get(load(paste0(outputDir_Ne, "/htz", sp, ".RData")))
  #     cat("=> loading already existing htz info\n")}
  #   var_htz <- sd(data_htz$HTZ[-which(is.nan(data_htz$HTZ))])
  #   
  #   var_list["var_Ne_max", incr_Simu_MULTIplier, sp] <- c(mx)
  #   var_list["var_Ne_mean", incr_Simu_MULTIplier, sp] <- c(mn)
  #   var_list["var_htz", incr_Simu_MULTIplier, sp] <- c(mn)
  # }
}


# 
# 
# 
# plot <- list()
# coef_vec <- c()
# for(sp in species){
#   cat(sp)
#   par(mfrow=c(1,2))
#   df_var <- as.data.frame(var_list[,,sp])
#   colnames(df_var) <- factor(colnames(df_var), levels=c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2,3,4,5))
#   df_var <- t(rbind("School_number_Simu_MULTIplier"=multipliers, df_var))
#   df_var <- df_var[,order(colnames(df_var))]
#   coef <- (max(df_var[, "var_htz"], na.rm = T)/max(df_var[,"var_Ne_max"], na.rm = T)) #which(df_long_tot$age==1 & 
#   coef_vec <- c(coef_vec, coef)
#   
#   if(any(is.nan(df_var))){next}
#   
#   plot[[sp]] <- ggplot(df_var, aes(x=School_number_Simu_MULTIplier))+
#     geom_point(aes(y=var_Ne_max))+
#     geom_point(aes(y=var_Ne_mean))+
#     geom_point(aes(y=var_htz/coef), col="red")+
#     scale_y_continuous(
#       # Features of the first axis
#       name = "Max and min std deviation of Ne",
#       # Add a second axis and specify its features
#       sec.axis = sec_axis(~.*coef, name=paste0("HTZ std deviation")))
#     labs(title = paste0("Mean and Max std deviation of Ne and Htz \nover ",  sp," cohorts, according to the \nnumber of school Simu_MULTIplier"))
# 
#   # plot(mx_var_list[[sp]], ylim = c(min(mn_var_list[[sp]], na.rm = T),
#   #                                  max(mx_var_list[[sp]], na.rm = T)))
#   # points(mn_var_list[[sp]], col="red")
#   # # readline(prompt="Press [enter] to continue")
# }
# 
# 
# 
# # STORE THE PLOTS
# pdf(file = paste0(configDir, "/output/Ne/plot_var_Simu_MULTIplier_longer.pdf"), width = 10)
# # pdf(file = paste0(paste0(pathNCDF, "_FINAL_plot/", "tot", ".pdf")), width=25, height=15)
# par(mfrow=c(1,1))
# for (k in 1:length(plot)) {
#   print(plot[k])
# }
# dev.off()




