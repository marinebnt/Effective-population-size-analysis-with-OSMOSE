
main = "New_gpu_conf.R"
source("osmose-ns-v4.4.1_clim/R/functions/source.R")

base_dir <- "osmose-ns-v4.4.1_clim/output/PLOTS/All_clean_genetic_drift_indicators/INIT/"
folders  <- list.dirs(base_dir, recursive = F)
INIT2 <- get(load(list.files(folders[1], full.names = T)))

folders_traits_y <- "osmose-ns-v4.4.1_clim/output/PLOTS/Trait_variables_by_year/RData/"
folders_traits_y <- list.dirs(folders_traits_y, recursive = F)

folders_traits_c <- "osmose-ns-v4.4.1_clim/output/PLOTS/Trait_variables_by_cohort/RData/"
folders_traits_c  <- list.dirs(folders_traits_c, recursive = F)

# Auxiliar functions ------------------------------------------------------

# load
source("osmose-ns-v4.4.1_clim/R/functions/auxiliary_for_lists.R")
source("osmose-ns-v4.4.1_clim/R/functions/functions_output_SimuRData_into_df_for_plots.R")



# Loop for Fishing Pressure -----------------------------------------------


pathout <- paste0(configDir, "/output/PLOTS/Htz_vs_Ne_variables")

# from RData output of the loops producing the plots of Trait by time and Genetic drift indicators, to df
# produce a dataframe with the right format : Ne variables and htz variables together

species_df=species
list_data_traits <- data_trait_and_genetic_into_df(Time1=150, Time2=160, pathoutput=paste0(configDir, "/output/PLOTS/RDATA"), 
                                                   FORCE = F, species = species_df)

# produce a dataframe with the right format : Ne variables and htz variables together
species_df=species
list_data_traits_ratio <- data_traitRATIO_and_genetic_into_df(Time1=150, Time2=160, pathoutput=paste0(configDir, "/output/PLOTS/RDATA"), 
                                                              FORCE = F, species = species_df)


# list of outputs 
list_outputs <- list("traits"=list_data_traits, "trait_ratio"=list_data_traits_ratio)








pathoutputRDATA = paste0(configDir, "/output/PLOTS/RDATA")

for (f in names(list_outputs)){
  
  genet_drift_data_full       <- list_outputs[[f]][["genet_drift_data"]]
  trait_data_full             <- list_outputs[[f]][["trait_data"]]

  # preparing a dataframe for the regression
  genet_drift_data_t <- bind_rows(genet_drift_data_full, .id = "variable")
  genet_drift_data_t <- rbind(genet_drift_data_t, trait_data_full$htz_list[, -which(colnames(trait_data_full$htz_list) == "SP_ref")])
  genet_drift_data_t$F <- paste0("Fmax=", (genet_drift_data_t$Fmax), " Fstart=", (genet_drift_data_t$Fstart))
  genet_drift_data_t <- genet_drift_data_t[, -which(colnames(genet_drift_data_t) %in% c("Fmax", "Fstart"))]

  # convert to wide format
  genet_drift_data_t_summary_wide <- pivot_wider(genet_drift_data_t, names_from = variable, values_from = value)

  # standardize the variables
  genet_drift_data_t_summary_wide <- genet_drift_data_t_summary_wide %>% 
    mutate(
      htz_list = scale(htz_slope),
      Ne_N_1 = scale(Ne_N_1),
      Ne_N_alpha = scale(Ne_N_alpha),
      Ne_temporal = scale(Ne_temporal)
    )
  genet_drift_data_t_summary_wide$Species <- as.factor(genet_drift_data_t_summary_wide$Species)
  genet_drift_data_t_summary_wide$F <- as.factor(genet_drift_data_t_summary_wide$F)

  # correlation plot
  for (spp_i in unique(genet_drift_data_t_summary_wide$Species)){
    pdf(paste0(pathout, "/correlation_plot/Htz_Ne_correlation_plot_", spp_i, "_", f, ".pdf"))
    for (f_i in unique(genet_drift_data_t_summary_wide$F)){  
      genet_drift_data_t_summary_wide_i <- genet_drift_data_t_summary_wide[which(genet_drift_data_t_summary_wide$F == f_i),]
      genet_drift_data_t_summary_wide_i <- genet_drift_data_t_summary_wide_i[which(genet_drift_data_t_summary_wide_i$Species == spp_i),]
      my_data <- genet_drift_data_t_summary_wide_i[, c("htz_list", "Ne_N_1", "Ne_N_alpha", "Ne_temporal", "Nb")]
      if (any(apply(my_data, 2, function(x) sum(!is.na(x))<3))) my_data <- my_data[, -which(apply(my_data, 2, function(x) sum(!is.na(x))<3))]
      if (nrow(na.omit(my_data)) < 2 | ncol(my_data)==0) next
      chart.Correlation(my_data, histogram=TRUE, pch=19, main=paste0("Correlation plot for ", spp_i, "\nat ", f_i))
    }
    dev.off()
  }

  # regressions
  # lm(htz_list ~ Ne_N_1 + Ne_N_alpha + Ne_temporal, data = genet_drift_data_t_summary_wide) -> model_slopes
  # lmer(htz_list ~ Ne_N_1 + Ne_N_alpha + Ne_temporal + (1|Species) + (1|F), data = genet_drift_data_t_summary_wide) -> model_slopes

}
