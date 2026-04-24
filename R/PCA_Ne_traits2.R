
main = "New_gpu_conf.R"
source("osmose-ns-v4.4.1_clim/R/functions/source.R")
base_dir <- "osmose-ns-v4.4.1_clim/output/PLOTS/All_clean_genetic_drift_indicators/INIT/"
folders  <- list.dirs(base_dir, recursive = F)
#INIT2 <- get(load(list.files(folders[1], full.names = T)))

folders_traits_y <- "osmose-ns-v4.4.1_clim/output/PLOTS/Trait_variables_by_year2/RData/"
folders_traits_y <- list.dirs(folders_traits_y, recursive = F)

folders_traits_c <- "osmose-ns-v4.4.1_clim/output/PLOTS/Trait_variables_by_cohort2/RData/"
folders_traits_c  <- list.dirs(folders_traits_c, recursive = F)


folders_traits_y_a <- "osmose-ns-v4.4.1_clim/output/PLOTS/Trait_variables_by_year_by_age2/RData/"
folders_traits_y_a <- list.dirs(folders_traits_y_a, recursive = F)

folders_traits_c_a <- "osmose-ns-v4.4.1_clim/output/PLOTS/Trait_variables_by_cohort_by_age2/RData/"
folders_traits_c_a  <- list.dirs(folders_traits_c_a, recursive = F)



# Auxiliar functions ------------------------------------------------------

# load
source("osmose-ns-v4.4.1_clim/R/functions/auxiliary_for_lists.R")
source("osmose-ns-v4.4.1_clim/R/functions/functions_output_SimuRData_into_df_for_plots.R")






# Loop for Fishing Pressure -----------------------------------------------


pathout <- paste0(configDir, "/output/PLOTS/Htz_vs_Ne_variables")

# produce a dataframe with the right format : Ne variables and htz variables together
species_df=species
# from RData output of the loops producing the plots of Trait by time and Genetic drift indicators, to df
# produce a dataframe with the right format : Ne variables and htz variables together
Time1=150
Time2=160
species_df=species
list_data_traits <- data_trait_and_genetic_into_df(Time1=150, Time2=160, pathoutput=paste0(configDir, "/output/PLOTS/RDATA"), 
                                                   FORCE = F, species = species_df)

# produce a dataframe with the right format : Ne variables and htz variables together
list_data_traits_ratio <- data_trait_and_genetic_into_df(Time1=150, Time2=160, pathoutput=paste0(configDir, "/output/PLOTS/RDATA"), 
                                                         FORCE = F, species = species_df, RATIO=T)
# list of outputs 
list_outputs <- list("traits"=list_data_traits, "trait_ratio"=list_data_traits_ratio)

# # Bind together the trait data coming from similar dimensions : cohort or year data together (by genetic drift dimension)

pathoutputRDATA = paste0(configDir, "/output/PLOTS/RDATA")
# data_full_merged <- dataset_attach_all_traits_genetic_drift(pathoutputRDATA, list_outputs=list("traits"=list_data_traits), 
#                                                             trait_cv = list("folders_traits_c"=folders_traits_c, "folders_traits_y"=folders_traits_y))
dataRATIO_full_merged <- dataset_attach_all_traits_genetic_drift(pathoutputRDATA, list_outputs=list("ratio"=list_data_traits_ratio), 
                                                                 trait_cv = list("folders_traits_c"=folders_traits_c, "folders_traits_y"=folders_traits_y), FORCE=F)

# 
# data_full_merged_age <- dataset_attach_all_traits_genetic_drift(pathoutputRDATA, list_outputs = list("traits"=list_data_traits), s_pecies = species,
#                                                                 trait_cv = list("folders_traits_c_a"=folders_traits_c_a, "folders_traits_y_a"=folders_traits_c_a),
#                                                                 keep_cohort_or_age_dimension = T, FORCE=F)

# 
# # Dataset ready for modeling --------------------------------------------
# # scaling as a standardized loop
# data_mean_sc = scale_dataset_list(pathoutputRDATA, data_mean=data_full_merged, FORCE=FALSE)
dataRATIO_mean_sc = scale_dataset_list(pathoutputRDATA, data_mean=dataRATIO_full_merged, FORCE=T)
# 
# # new dataset of the difference between Fishing pressure 1,1 and Else--------------------
# diff_data = relative_dataset(data_mean_sc, pathoutputRDATA, FORCE=T)
diffRATIO_data = relative_dataset(dataRATIO_mean_sc, pathoutputRDATA, FORCE=T, RATIO=T)



############ PLOT ************************

plot_PCA <- function (data, DIFF = FALSE, Traits_raw = FALSE){
  for (j in names(data)){
    P <- list()
    dat_merged_w_sc <- data[[j]]
    if ("variable" %in% colnames(dat_merged_w_sc)) dat_merged_w_sc <- pivot_wider(dat_merged_w_sc, names_from = "variable", values_from = "value")
    # ACP
    dat_merged_pca <- dat_merged_w_sc
    dat_merged_pca$Fmax <- as.numeric(dat_merged_pca$Fmax)
    dat_merged_pca$Fstart <- as.numeric(dat_merged_pca$Fstart)
    dat_merged_pca$replicate <- as.numeric(dat_merged_pca$replicate)
    
    if (Traits_raw == TRUE) {
      ID_keep_dat <- c(grep(colnames(dat_merged_pca), pattern="cohort"), 
                       grep(colnames(dat_merged_pca), pattern="year"), 
                       which(colnames(dat_merged_pca) %in% c("replicate", "Fstart", "Fmax", "Species")), 
                       which(colnames(dat_merged_pca) %in% c("mean_Ne_N_alpha", "mean_ratio_Ne_N_alpha", "mean_Ne_N_1", "mean_ratio_Ne_N_1",
                                                             "mean_Nb", "mean_ratio_Nb", "htz_slope", "mean_Ne_temporal", "mean_ratio_Ne_temporal",
                                                             "diff_mean_Ne_N_alpha", "diff_mean_ratio_Ne_N_alpha", "diff_mean_Ne_N_1", "diff_mean_ratio_Ne_N_1",
                                                             "diff_mean_Nb", "diff_mean_ratio_Nb", "diff_htz_slope", "diff_mean_Ne_temporal", "diff_mean_ratio_Ne_temporal"))
      )
      dat_merged_pca <- dat_merged_pca[, ID_keep_dat]
    }
    
    
    colnm_div_genet <- which(colnames(dat_merged_pca) %in% c("mean_Ne_N_alpha", "mean_ratio_Ne_N_alpha", "mean_Ne_N_1", "mean_ratio_Ne_N_1",
                                                             "mean_Nb", "mean_ratio_Nb", "htz_slope", "mean_Ne_temporal", "mean_ratio_Ne_temporal",
                                                             "diff_mean_Ne_N_alpha", "diff_mean_ratio_Ne_N_alpha", "diff_mean_Ne_N_1", "diff_mean_ratio_Ne_N_1",
                                                             "diff_mean_Nb", "diff_mean_ratio_Nb", "diff_htz_slope", "diff_mean_Ne_temporal", "diff_mean_ratio_Ne_temporal"))
    for (div_genet in colnm_div_genet){
      
      
      data_PCA <- dat_merged_pca
      data_PCA <- data_PCA[,-colnm_div_genet]
      assign(colnames(dat_merged_pca)[div_genet], dat_merged_pca[div_genet])
      
      pca1 <- PCA(data_PCA |> dplyr::select(where(is.numeric)),
                  scale.unit = TRUE,
                  graph = FALSE)
      
      for (spi in unique(data_PCA$Species)){
        ind <- as.data.frame(pca1$ind$coord)
        ind$Species <- dat_merged_pca$Species
        ind[colnames(dat_merged_pca)[div_genet]] <- get(colnames(dat_merged_pca)[div_genet])
        ind[colnames(dat_merged_pca)[div_genet]][which(is.na(ind[colnames(dat_merged_pca)[div_genet]])),] <- 0
        ind <- ind[-which(ind[colnames(dat_merged_pca)[div_genet]]==0),]
        ind <- ind[which(ind$Species %in% spi),]
        
        Var <- as.data.frame(pca1$var$coord)
        scale <- mean(unlist(as.vector(Var[,1:2])), na.rm=T)/quantile(unlist(as.vector(ind[,1:2])), na.rm=T, 0.70)
        if (scale>1) scale = 1/(2*scale)
        if (abs(scale) < 0.05) scale = 3*scale
        Var <- Var/scale
        Var$var <- rownames(Var)
        
        P[[spi]] <- ggplot() +  
          stat_summary_hex(fun = mean, bins = 25, data=ind, 
                           aes(Dim.1, Dim.2, z=.data[[colnames(dat_merged_pca)[div_genet]]])) +
          geom_point(data=ind, 
                     aes(Dim.1, Dim.2, color = .data[[colnames(dat_merged_pca)[div_genet]]]), alpha=0.5) +
          
          geom_segment(data = Var, aes(x=0, y=0, xend = Dim.1, yend = Dim.2), arrow = arrow(length = unit(0.1, "inches"))
          )+  
          geom_text(data = Var, aes(x = Dim.1, y = Dim.2, label = var), size=3)+
          scale_color_viridis_c() +
          theme_minimal() +
          geom_abline(slope = 0, intercept = 0, linetype = "dashed", color = "grey") +
          geom_hline(yintercept = 0, linetype = "dashed", color = "grey") +
          labs(title = paste0("PCA of traits vs ", colnames(dat_merged_pca)[div_genet]),
               subtitle = spi,
               x = paste0("PC1, ", round(pca1$eig[1,2], 2), "%"),
               y = paste0("PC2, ", round(pca1$eig[2,2], 2), "%"))
        
      }
      pathoutputplots = paste0(configDir, "/output/PLOTS/Htz_vs_Ne_variables/plots")
      if(!dir.exists(pathoutputplots)) dir.create(pathoutputplots, recursive = T)
      name_pdf = paste0(pathoutputplots, "/", 
                        ifelse(grepl(colnames(dat_merged_pca)[div_genet], pattern="htz"), paste0(j, colnames(dat_merged_pca)[div_genet]), colnames(dat_merged_pca)[div_genet]),
                        ifelse(Traits_raw==TRUE, "_raw_traits", ""),ifelse(DIFF==TRUE, "_DIFF", ""), 
                        "_one_genet_drift_PCA_htz_vs_traits.pdf")
      pdf(name_pdf, width = 10, height = 8)
      print(P)
      dev.off()
    }
  }
}
plot_PCA(data=dataRATIO_full_merged)
# plot_PCA(data=diffRATIO_data, DIFF=TRUE)
plot_PCA(data=dataRATIO_full_merged, Traits_raw = T)
# plot_PCA(data=diffRATIO_data, DIFF=TRUE, Traits_raw = T)

plot_PCA_sep_species <- function (data, DIFF = FALSE, Traits_raw = FALSE){
  for (j in names(data)){
    P <- list()
    dat_merged_w_sc <- data[[j]]
    if ("variable" %in% colnames(dat_merged_w_sc)) dat_merged_w_sc <- pivot_wider(dat_merged_w_sc, names_from = "variable", values_from = "value")
    # ACP
    dat_merged_pca <- dat_merged_w_sc
    dat_merged_pca$Fmax <- as.numeric(dat_merged_pca$Fmax)
    dat_merged_pca$Fstart <- as.numeric(dat_merged_pca$Fstart)
    dat_merged_pca$replicate <- as.numeric(dat_merged_pca$replicate)
    
    if (Traits_raw == TRUE) {
      ID_keep_dat <- c(grep(colnames(dat_merged_pca), pattern="cohort"), 
                       grep(colnames(dat_merged_pca), pattern="year"), 
                       which(colnames(dat_merged_pca) %in% c("replicate", "Fstart", "Fmax", "Species")), 
                       which(colnames(dat_merged_pca) %in% c("mean_Ne_N_alpha", "mean_ratio_Ne_N_alpha", "mean_Ne_N_1", "mean_ratio_Ne_N_1",
                                                             "mean_Nb", "mean_ratio_Nb", "htz_slope", "mean_Ne_temporal", "mean_ratio_Ne_temporal",
                                                             "diff_mean_Ne_N_alpha", "diff_mean_ratio_Ne_N_alpha", "diff_mean_Ne_N_1", "diff_mean_ratio_Ne_N_1",
                                                             "diff_mean_Nb", "diff_mean_ratio_Nb", "diff_htz_slope", "diff_mean_Ne_temporal", "diff_mean_ratio_Ne_temporal"))
      )
      dat_merged_pca <- dat_merged_pca[, ID_keep_dat]
    }
    
    
    colnm_div_genet <- which(colnames(dat_merged_pca) %in% c("mean_Ne_N_alpha", "mean_ratio_Ne_N_alpha", "mean_Ne_N_1", "mean_ratio_Ne_N_1",
                                                             "mean_Nb", "mean_ratio_Nb", "htz_slope", "mean_Ne_temporal", "mean_ratio_Ne_temporal",
                                                             "diff_mean_Ne_N_alpha", "diff_mean_ratio_Ne_N_alpha", "diff_mean_Ne_N_1", "diff_mean_ratio_Ne_N_1",
                                                             "diff_mean_Nb", "diff_mean_ratio_Nb", "diff_htz_slope", "diff_mean_Ne_temporal", "diff_mean_ratio_Ne_temporal"))
    for (div_genet in colnm_div_genet){
      
      
      data_PCA <- dat_merged_pca
      data_PCA <- data_PCA[,-colnm_div_genet]
      assign(colnames(dat_merged_pca)[div_genet], dat_merged_pca[div_genet])
      
      
      for (spi in unique(data_PCA$Species)){
        
        data_PCA_sp <- data_PCA[which(data_PCA$Species == spi), ]
        pca1 <- PCA(data_PCA_sp |> dplyr::select(where(is.numeric)),
                    scale.unit = TRUE,
                    graph = FALSE)
        
        ind <- as.data.frame(pca1$ind$coord)
        ind$Species <- data_PCA_sp$Species
        ind[colnames(dat_merged_pca)[div_genet]] <- get(colnames(dat_merged_pca)[div_genet])[which(dat_merged_pca$Species == spi),1]
        if (any((is.na(ind[colnames(dat_merged_pca)[div_genet]])))) ind[colnames(dat_merged_pca)[div_genet]][which(is.na(ind[colnames(dat_merged_pca)[div_genet]])),] <- 0
        if (any((ind[colnames(dat_merged_pca)[div_genet]]==0)))     ind <- ind[-which(ind[colnames(dat_merged_pca)[div_genet]]==0),]
        
        Var <- as.data.frame(pca1$var$coord)
        scale <- mean(unlist(as.vector(Var[,1:2])), na.rm=T)/quantile(unlist(as.vector(ind[,1:2])), na.rm=T, 0.70)
        if (scale>1) scale = 1/(2*scale)
        if (abs(scale) < 0.05) scale = 3*scale
        Var <- Var/scale
        Var$var <- rownames(Var)
        
        P[[spi]] <- ggplot() +  
          stat_summary_hex(fun = mean, bins = 25, data=ind, 
                           aes(Dim.1, Dim.2, z=.data[[colnames(dat_merged_pca)[div_genet]]])) +
          geom_point(data=ind, 
                     aes(Dim.1, Dim.2, color = .data[[colnames(dat_merged_pca)[div_genet]]]), alpha=0.5) + #
          
          geom_segment(data = Var, aes(x=0, y=0, xend = Dim.1, yend = Dim.2), arrow = arrow(length = unit(0.1, "inches"))
          )+  
          geom_text(data = Var, aes(x = Dim.1, y = Dim.2, label = var), size=3)+
          scale_color_viridis_c() +
          theme_minimal() +
          geom_abline(slope = 0, intercept = 0, linetype = "dashed", color = "grey") +
          geom_hline(yintercept = 0, linetype = "dashed", color = "grey") +
          labs(title = paste0("PCA of traits vs ", colnames(dat_merged_pca)[div_genet]),
               subtitle = spi,
               x = paste0("PC1, ", round(pca1$eig[1,2], 2), "%"),
               y = paste0("PC2, ", round(pca1$eig[2,2], 2), "%"))
        
      }
      pathoutputplots = paste0(configDir, "/output/PLOTS/Htz_vs_Ne_variables/plots")
      if(!dir.exists(pathoutputplots)) dir.create(pathoutputplots, recursive = T)
      name_pdf = paste0(pathoutputplots, "/SEP_SPECIES_", 
                        ifelse(grepl(colnames(dat_merged_pca)[div_genet], pattern="htz"), paste0(j, colnames(dat_merged_pca)[div_genet]), colnames(dat_merged_pca)[div_genet]),
                        ifelse(Traits_raw==TRUE, "_raw_traits", ""),ifelse(DIFF==TRUE, "_DIFF", ""), 
                        "_one_genet_drift_PCA_htz_vs_traits.pdf")
      pdf(name_pdf, width = 10, height = 8)
      print(P)
      dev.off()
    }
  }
}
plot_PCA_sep_species(data=dataRATIO_full_merged)
# plot_PCA_sep_species(data=diffRATIO_data, DIFF = T)
plot_PCA_sep_species(data=dataRATIO_full_merged, Traits_raw = T)
# plot_PCA_sep_species(data=diffRATIO_data, DIFF = T, Traits_raw = T)



# ALL TRAITS, all dimensions, all species
plot_PCA_sep_species <- function (data, DIFF = FALSE, Traits_raw = FALSE){
  for (j in names(data)){
    P <- list()
    dat_merged_w_sc <- data[[j]]
    if ("variable" %in% colnames(dat_merged_w_sc)) dat_merged_w_sc <- pivot_wider(dat_merged_w_sc, names_from = "variable", values_from = "value")
    # ACP
    dat_merged_pca <- dat_merged_w_sc
    dat_merged_pca$Fmax <- as.numeric(dat_merged_pca$Fmax)
    dat_merged_pca$Fstart <- as.numeric(dat_merged_pca$Fstart)
    dat_merged_pca$replicate <- as.numeric(dat_merged_pca$replicate)

    if (Traits_raw == TRUE) {
      ID_keep_dat <- c(grep(colnames(dat_merged_pca), pattern="cohort"),
                       grep(colnames(dat_merged_pca), pattern="year"),
                       which(colnames(dat_merged_pca) %in% c("replicate", "Fstart", "Fmax", "Species")),
                       which(colnames(dat_merged_pca) %in% c("mean_Ne_N_alpha", "mean_ratio_Ne_N_alpha", "mean_Ne_N_1", "mean_ratio_Ne_N_1",
                                                             "mean_Nb", "mean_ratio_Nb", "htz_slope", "mean_Ne_temporal", "mean_ratio_Ne_temporal",
                                                             "diff_mean_Ne_N_alpha", "diff_mean_ratio_Ne_N_alpha", "diff_mean_Ne_N_1", "diff_mean_ratio_Ne_N_1",
                                                             "diff_mean_Nb", "diff_mean_ratio_Nb", "diff_htz_slope", "diff_mean_Ne_temporal", "diff_mean_ratio_Ne_temporal"))
      )
      dat_merged_pca <- dat_merged_pca[, ID_keep_dat]
    }


    colnm_div_genet <- which(colnames(dat_merged_pca) %in% c("mean_Ne_N_alpha", "mean_ratio_Ne_N_alpha", "mean_Ne_N_1", "mean_ratio_Ne_N_1",
                                                             "mean_Nb", "mean_ratio_Nb", "htz_slope", "mean_Ne_temporal", "mean_ratio_Ne_temporal",
                                                             "diff_mean_Ne_N_alpha", "diff_mean_ratio_Ne_N_alpha", "diff_mean_Ne_N_1", "diff_mean_ratio_Ne_N_1",
                                                             "diff_mean_Nb", "diff_mean_ratio_Nb", "diff_htz_slope", "diff_mean_Ne_temporal", "diff_mean_ratio_Ne_temporal"))
    for (div_genet in colnm_div_genet){


      data_PCA <- dat_merged_pca
      data_PCA <- data_PCA[,-colnm_div_genet]
      assign(colnames(dat_merged_pca)[div_genet], dat_merged_pca[div_genet])


      # for (spi in unique(data_PCA$Species)){

      data_PCA_sp <- data_PCA#[which(data_PCA$Species == spi), ]
      pca1 <- PCA(data_PCA_sp |> dplyr::select(where(is.numeric)),
                  scale.unit = TRUE,
                  graph = FALSE)

      ind <- as.data.frame(pca1$ind$coord)
      ind$Species <- data_PCA_sp$Species
      ind[colnames(dat_merged_pca)[div_genet]] <- get(colnames(dat_merged_pca)[div_genet])#[which(dat_merged_pca$Species == spi),1]
      if (any((is.na(ind[colnames(dat_merged_pca)[div_genet]])))) ind[colnames(dat_merged_pca)[div_genet]][which(is.na(ind[colnames(dat_merged_pca)[div_genet]])),] <- 0
      if (any((ind[colnames(dat_merged_pca)[div_genet]]==0)))     ind <- ind[-which(ind[colnames(dat_merged_pca)[div_genet]]==0),]

      Var <- as.data.frame(pca1$var$coord)
      scale <- mean(unlist(as.vector(Var[,1:2])), na.rm=T)/quantile(unlist(as.vector(ind[,1:2])), na.rm=T, 0.70)
      if (scale>1) scale = 1/(2*scale)
      Var <- Var/scale
      Var$var <- rownames(Var)

      #[[spi]]
      P<- ggplot() +
        stat_summary_hex(fun = mean, bins = 25, data=ind,
                         aes(Dim.1, Dim.2, z=.data[[colnames(dat_merged_pca)[div_genet]]])) +
        geom_point(data=ind,
                   aes(Dim.1, Dim.2, color = Species), alpha=0.5) + #.data[[colnames(dat_merged_pca)[div_genet]]]

        geom_segment(data = Var, aes(x=0, y=0, xend = Dim.1, yend = Dim.2), arrow = arrow(length = unit(0.1, "inches"))
        )+
        geom_text(data = Var, aes(x = Dim.1, y = Dim.2, label = var), size=3)+
        scale_color_viridis_d() +
        theme_minimal() +
        geom_abline(slope = 0, intercept = 0, linetype = "dashed", color = "grey") +
        geom_hline(yintercept = 0, linetype = "dashed", color = "grey") +
        labs(title = paste0("PCA of traits vs ", colnames(dat_merged_pca)[div_genet]),
             #subtitle = spi,
             x = paste0("PC1, ", round(pca1$eig[1,2], 2), "%"),
             y = paste0("PC2, ", round(pca1$eig[2,2], 2), "%"))

      #}
      pathoutputplots = paste0(configDir, "/output/PLOTS/Htz_vs_Ne_variables/plots")
      if(!dir.exists(pathoutputplots)) dir.create(pathoutputplots, recursive = T)
      name_pdf = paste0(pathoutputplots, "/SEP_SPECIES_",
                        ifelse(grepl(colnames(dat_merged_pca)[div_genet], pattern="htz"), paste0(j, colnames(dat_merged_pca)[div_genet]), colnames(dat_merged_pca)[div_genet]),
                        ifelse(Traits_raw==TRUE, "_raw_traits", ""),ifelse(DIFF==TRUE, "_DIFF", ""),
                        # "_one_genet_drift_PCA_htz_vs_traits.pdf"
                        "ALL.pdf")
      pdf(name_pdf, width = 10, height = 8)
      print(P)
      dev.off()
    }
  }
}