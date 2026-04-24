#* marine beneat
#* 15/05/2025
#* Following Waples 2011
#* Idea : produce an output of the model giving the effective population size over time
#* Suggestion : run it in background and save it on the go, it is heavyier as you increase the species number
# Ne is estimated from a combination of parameters estimated as outputs of output_Ne.R

library(osmose)
library(stringr)
library(reshape2)
library(ggplot2)
library(tidyr)

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

Ne_version = "restarting_everything"

# 2. Read inputs ------------------------

data_Ne_raw <- get(load(paste0(configDir, "/output/Ne/", Ne_version, "/Data_for_Ne.RData")))
ns1 = read_osmose(path = outputDir)

dt       <- get_par(conf, "simulation.time.ndtPerYear")/get_par(conf, "output.recordfrequency.ndt")
ny       <- get_par(conf, "simulation.time.nyear")-get_par(conf, "output.start.year")
time_vec <- seq(0,ny, by=1/dt)[-1] # if year 0 is removed
nb_rep   <- get_par(conf, "simulation.nsimulation")
nb_replicates = nb_rep
nb_sp    <- get_par(conf, "simulation.nspecies")
species <- unlist(as.matrix(get_par(conf, "species.name.")[1:nb_sp]))
seeding_y <- get_par(conf, "population.seeding.year.max")

# 3. Ne estimate

source(file.path(configDir, "R/functions/functions_Ne.R"))
data_Ne <- estimate_Ne(data_Ne=data_Ne_raw, species, nb_replicates)


# 3. Ne TEMPORAL estimate

source(file.path(configDir, "R/functions/functions_temporalNe.R"))
data_Ne_t <- estimate_temporalNe(species = species, seeding_y = seeding_y, conf=conf, data_Ne = data_Ne_raw, 
                               outputDir_Simu = file.path(configDir, "output", simulation))









# old plots scripts

# for(sp in species){
#   rep=1
#   cat(sp, "\n")
#   Ne <- as.matrix(data_Ne[["Ne"]][[sp]][[paste0("rep_", rep)]])
#   if(all(is.na(Ne))){next}
#   Ne[which(is.nan(Ne), arr.ind = T)] <- NA
#   variability_over_time <- apply(Ne, 1, sd)
#   cat(var(na.omit(variability_over_time)), "\n") #homoscédasticité
#   print(plot(variability_over_time))
#   readline(prompt="Press [enter] to continue")
# }
# 
# 
# # ou plot en 3D
# 
# # Plotting the cohorts Ne
# 
# # Assuming your data is in a data frame called `df`
# # Rows: time steps, Columns: cohort abundances
# # Example: df <- data.frame(Cohort1 = ..., Cohort2 = ..., ...)
# 
# par(mfrow = c(3,3))
# k=0
# Nc_Ne <- list()
# multi_vec <- c()
# for(sp in species){
#   cat(sp, "\n")
#   k=k+1
#   
#   # Nc data 
#   data <- as.data.frame(ns1[["abundance"]])
#   time <- rownames(data)
#   data <- data[,grep(sp, colnames(data))]
#   data_Nc_frame <- data.frame(Time=as.numeric(time), Nc=as.numeric(data))
#   
#   # Ne data
#   df <- as.data.frame(data_Ne[["Ne"]][[sp]][[paste0("rep_", rep)]])
#   if(all(is.na(df))){next}
#   # Add time as a column
#   df$Time <- as.numeric(rownames(df))
#   # Reshape data to long format
#   df_long <- pivot_longer(df, cols = -Time, names_to = "Cohort", values_to = "Ne")
#   
#   # bind datasets
#   df_long_tot <- dplyr::left_join(df_long, data_Nc_frame, by = "Time")
#   multip <- (max(df_long_tot$Nc, na.rm = T)/max(df_long_tot$Ne, na.rm = T))
#   
#   if(max(df_long_tot$Ne, na.rm = T)=="Inf"){
#     cat("max Ne is Infinite \n")
#     next
#   }
#   
#   # Plot using ggplot2
# 
#   a <- ggplot(df_long_tot, aes(x = Time)) +
#     geom_line(aes(y = Ne, colour = Cohort)) +
#     geom_line(aes(y = Nc/multip), lwd = 0.75)+
#     labs(title = paste0(sp, ": Comparing Ne and Nc"), x = "Time Step", y = "Ne") +
#     theme_minimal()+ 
#     theme(legend.position="none")+
#     scale_y_continuous(
#       # Features of the first axis
#       name = "Ne (color by cohort)",
#       # Add a second axis and specify its features
#       sec.axis = sec_axis(~.*multip, name="Total species Nc (black)")
#     ) 
#   par(mfrow = c(3,3))
#   print(a)
#   readline(prompt="Press [enter] to continue")
#  
#   multi_vec <- c(multi_vec, multip)
#   Nc_Ne[[k]] <- a
# }
# 
# 
# par(mfrow = c(3,3))
# k=0
# HTZ_Ne <- list()
# coef_vec <- c()
# nm_vec <- c()
# for(sp in species){
#   cat(sp, "\n")
#   k=k+1
#   
#   # Nc data 
#   data <- ncdf4::nc_open(paste0(outputDir, "/Genetic/ns_observedHtz-", sp, "Simu0.nc"))
#   data_nc <- ncdf4::ncvar_get(data, "observed_heterozygosity")
#   trait1_data <- data_nc[, 1, ]
#   mean_htz <- colMeans(trait1_data)
#   time <- 0:c(dim(data_nc)[3]-1)
#   data_htz <- data.frame(Time=time, HTZ=mean_htz)
#   
#   # Ne data
#   df <- as.data.frame(data_Ne[["Ne"]][[sp]][[paste0("rep_", rep)]])
#   if(all(is.na(df))){next}
#   # Add time as a column
#   df$Time <- as.numeric(rownames(df))
#   # Reshape data to long format
#   df_long <- pivot_longer(df, cols = -Time, names_to = "Cohort", values_to = "Ne")
#   
#   # bind datasets
#   df_long_tot <- dplyr::left_join(df_long, data_htz, by = "Time")
#   coef <- (max(df_long_tot$HTZ, na.rm = T)/max(df_long_tot$Ne, na.rm = T))
#   lm_htz <- lm(data_htz$HTZ~data_htz$Time)
#   slope <- lm_htz$coefficients[2]
#   
#   if(max(df_long_tot$Ne, na.rm = T)=="Inf"){
#     cat("max Ne is Infinite \n")
#     next
#   }
#   
#   # Plot using ggplot2
#   a1 <- ggplot(df_long_tot, aes(x = Time)) +
#     geom_line(aes(y = Ne, colour = Cohort)) +
#     geom_line(aes(y = HTZ/coef), lwd = 0.75)+
#     geom_smooth(method = lm, aes(y = HTZ/coef), lwd=0.7)+
#     labs(title = "Comparing observed HTZ and Ne" , x = "Time Step", y = "Ne") +
#     theme_minimal()+ 
#     theme(legend.position="none")+
#     scale_y_continuous(
#       # Features of the first axis
#       name = "Ne (color by cohort)",
#       # Add a second axis and specify its features
#       sec.axis = sec_axis(~.*coef, name=paste0("Mean species HTZ (black), slope lm model (blue)= ", slope))
#     ) 
#   par(mfrow = c(3,3))
#   print(a1)
#   # readline(prompt="Press [enter] to continue")
#   
#   readline(prompt="Press [enter] to continue")
#   coef_vec <- c(coef_vec, coef)
#   HTZ_Ne[[k]] <- a1
#   nm_vec[k] <- sp
#   
# }
# 
# 
# # Save as .pdf
# pdf(file = paste0(configDir, "/output/Ne/", Ne_version, "/plot_Nc_Ne.pdf"), width = 10)
# # pdf(file = paste0(paste0(pathNCDF, "_FINAL_plot/", "tot", ".pdf")), width=25, height=15)
# par(mfrow=c(1,1))
# for (k in 1:length(Nc_Ne)) {
#   if(length(Nc_Ne[[k]])==0){
#     next
#   }
#   multip <- multi_vec[k]
#   coef <- coef_vec[k]
#   nm <- nm_vec[k]
#   print(Nc_Ne[[k]])
# }
# dev.off()
# 
# # Save as .pdf
# pdf(file = paste0(configDir, "/output/Ne/", Ne_version, "/plot_HTZ_Ne.pdf"), width = 10)
# # pdf(file = paste0(paste0(pathNCDF, "_FINAL_plot/", "tot", ".pdf")), width=25, height=15)
# par(mfrow=c(1,1))
# for (k in 1:length(Nc_Ne)) {
#   if(length(Nc_Ne[[k]])==0){
#     next
#   }
#   multip <- multi_vec[k]
#   coef <- coef_vec[k]
#   print(HTZ_Ne[[k]])
# }
# dev.off()
# 
# # Save as .pdf
# pdf(file = paste0(configDir, "/output/Ne/", Ne_version, "/plot_Ne.pdf"), width = 10)
# # pdf(file = paste0(paste0(pathNCDF, "_FINAL_plot/", "tot", ".pdf")), width=25, height=15)
# par(mfrow=c(1,1))
# for (k in 1:length(Nc_Ne)) {
#   if(length(Nc_Ne[[k]])==0){
#     next
#   }
#   multip <- multi_vec[k]
#   coef <- coef_vec[k]
#   print(cowplot::plot_grid(Nc_Ne[[k]], HTZ_Ne[[k]], ncol=2, nrow=1))
# }
# dev.off()
# 
