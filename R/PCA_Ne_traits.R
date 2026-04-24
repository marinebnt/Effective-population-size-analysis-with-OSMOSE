
main = "New_gpu_conf.R"
source("osmose-ns-v4.4.1_clim/R/functions/source.R")
base_dir <- "osmose-ns-v4.4.1_clim/output/PLOTS/All_clean_genetic_drift_indicators/INIT/"
folders  <- list.dirs(base_dir, recursive = F)
INIT2 <- get(load(list.files(folders[1], full.names = T)))

folders_traits_y <- "osmose-ns-v4.4.1_clim/output/PLOTS/Trait_variables_by_year/RData/"
folders_traits_y <- list.dirs(folders_traits_y, recursive = F)

folders_traits_c <- "osmose-ns-v4.4.1_clim/output/PLOTS/Trait_variables_by_cohort/RData/"
folders_traits_c  <- list.dirs(folders_traits_c, recursive = F)


folders_traits_y_a <- "osmose-ns-v4.4.1_clim/output/PLOTS/Trait_variables_by_year_by_age/RData/"
folders_traits_y_a <- list.dirs(folders_traits_y_a, recursive = F)

folders_traits_c_a <- "osmose-ns-v4.4.1_clim/output/PLOTS/Trait_variables_by_cohort_by_age/RData/"
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
data_full_merged <- dataset_attach_all_traits_genetic_drift(pathoutputRDATA, list_outputs=list("traits"=list_data_traits), 
                                                            trait_cv = list("folders_traits_c"=folders_traits_c, "folders_traits_y"=folders_traits_y))
dataRATIO_full_merged <- dataset_attach_all_traits_genetic_drift(pathoutputRDATA, list_outputs=list("ratio"=list_data_traits_ratio), 
                                                                 trait_cv = list("folders_traits_c"=folders_traits_c, "folders_traits_y"=folders_traits_y))


data_full_merged_age <- dataset_attach_all_traits_genetic_drift(pathoutputRDATA, list_outputs = list("traits"=list_data_traits), s_pecies = species,
                                                                trait_cv = list("folders_traits_c_a"=folders_traits_c_a, "folders_traits_y_a"=folders_traits_c_a),
                                                                keep_cohort_or_age_dimension = T, FORCE=F)

# 
# # Dataset ready for modeling --------------------------------------------
# # scaling as a standardized loop
data_mean_sc = scale_dataset_list(pathoutputRDATA, data_mean=data_full_merged[[1]], FORCE=FALSE)
dataRATIO_mean_sc = scale_dataset_list(pathoutputRDATA, data_mean=dataRATIO_full_merged[[1]], FORCE=FALSE)
# 
# # new dataset of the difference between Fishing pressure 1,1 and Else--------------------
diff_data = relative_dataset(data_mean_sc, pathoutputRDATA, FORCE=F)
diffRATIO_data = relative_dataset(dataRATIO_mean_sc, pathoutputRDATA, FORCE=FALSE, RATIO=T)


# ANALYSIS --------------------


# Modeling and PCA --------------------------------------------------------

p <- list()
for (j in names(data_mean_sc)){

    dat_merged_w_sc <- data_mean_sc[[j]]

    # model1 <- NULL
    # modelmer <- NULL
    # model1 <- lme4::lmer(formula = as.formula(paste0("htz_slope", " ~", 
    #         paste(colnames(dat_merged_w_sc)[-which(colnames(dat_merged_w_sc) %in% 
    #         c("Species", "replicate", "Fmax", "Fstart", "htz_slope"))], collapse="+"), 
    #         "+ (1|Fmax) + (1|Fstart) + (1|replicate) + (1|Species)")), data=dat_merged_w_sc)
    # modelmer <- buildmer(as.formula(paste0("htz_slope", " ~", 
    #         paste(colnames(dat_merged_w_sc)[-which(colnames(dat_merged_w_sc) %in% 
    #         c("Species", "replicate", "Fmax", "Fstart", "htz_slope"))], collapse="+"), 
    #         "+ (1|Fmax) + (1|Fstart) + (1|replicate) + (1|Species)")), data=dat_merged_w_sc,)

    for (spi in unique(dat_merged_w_sc$Species)){
      dat_merged_w_sc$Species <- as.character(dat_merged_w_sc$Species)
      dat_merged_w_sc_sp <- dat_merged_w_sc[which(dat_merged_w_sc$Species == spi), ]

      # ACP
      dat_merged_pca <- dat_merged_w_sc_sp
      dat_merged_pca$Fmax <- as.numeric(dat_merged_pca$Fmax)
      dat_merged_pca$Fstart <- as.numeric(dat_merged_pca$Fstart)
      dat_merged_pca$replicate <- as.numeric(dat_merged_pca$replicate)
      dat_merged_pca <- dat_merged_pca[, -which(colnames(dat_merged_pca) %in% c("htz_slope"))]

      pca1 <- PCA(dat_merged_pca |> dplyr::select(where(is.numeric)),
              scale.unit = TRUE,
              graph = FALSE)

      # # Variance expliquée
      # fviz_eig(pca1)

      # # Individus
      # fviz_pca_ind(pca1,
      #             geom = "point",
      #             pointsize = 2,
      #             labels = NULL,)

      # # Variables
      # fviz_pca_var(pca1)

      # Biplot
      var_to_plot <- grep("Mpred", colnames(dat_merged_w_sc_sp))
      
      p[[spi]] <- fviz_pca_biplot(pca1, label = "var", psize = 2,
                  col.ind = dat_merged_w_sc_sp[[var_to_plot]],
                      repel = TRUE, 
                      gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"))+
          labs(title = paste0("PCA of traits vs Htz_slope for ", spi, " species - dataset: ", j),
               color = "Mean Mpred per age")
    
      # #with factors
      # pca2 <- PCA(dat_merged_pca |> dplyr::select(where(is.numeric)),
      #         scale.unit = TRUE,
      #         graph = FALSE)

      # fviz_pca_ind(pca2,
      #             addEllipses = TRUE,
      #             ellipse.level = 0.95)

      # # by group
      # dat_merged_pca |> 
      # group_split(Species) |> 
      # map(~ PCA(select(.x, where(is.numeric)),
      #             scale.unit = TRUE,
      #             graph = FALSE))
  }
  pathoutputplots = paste0(configDir, "/output/PLOTS/Htz_vs_Ne_variables/plots")
  pdf(paste0(pathoutputplots, "/MPRED_PCA_htz_vs_traits_", j, ".pdf"), width = 10, height = 8)
  print(p)
  dev.off()
}






a <- data_mean_sc
a_long <- NULL
for (a_n in names(a)){
  if (is.null(a_long)){
    a_long <- a[[a_n]]
  } else {
    datatoadd1 <- a[[a_n]][,-which(colnames(a[[a_n]]) %in% colnames(a_long))]
    datatoadd  <- cbind(a[[a_n]][,which(colnames(a[[a_n]]) %in% c("Species", "replicate", "Fmax", "Fstart"))], datatoadd1)
    if (!is.null(datatoadd1)) a_long <- full_join(a_long, datatoadd)
  }
}

# boxplot(a_long$htz_slope~a_long$Species)
# boxplot(a_long$Ne_N_alpha~a_long$Species)
# boxplot(a_long$Ne_N_1~a_long$Species)
# boxplot(a_long$Ne_temporal~a_long$Species)
# boxplot(a_long$Nb~a_long$Species)

q <- list()
for (spi in unique(a_long$Species)){

  # a_long_t <- a_long[, which(colnames(a_long) %in% c("Species", "Fmax", "Fstart", "replicate",
  #                      "Ne_N_1", "Ne_N_alpha", "Nb", "tempoNe", "htz_slope"))]
  a_long_t <- a_long

  a_long_t$Species <- as.character(a_long_t$Species)
  # a_long_sp <- a_long_t[which(a_long_t$Species == spi), ]
  a_long_sp <- a_long_t

  # ACP
  dat_merged_pca <- a_long_sp
  dat_merged_pca$Fmax <- as.numeric(dat_merged_pca$Fmax)
  dat_merged_pca$Fstart <- as.numeric(dat_merged_pca$Fstart)
  dat_merged_pca$replicate <- as.numeric(dat_merged_pca$replicate)
  # dat_merged_pca <- dat_merged_pca[, -which(colnames(dat_merged_pca) %in% c("htz_slope"))]
  

  pca1 <- PCA(dat_merged_pca |> dplyr::select(where(is.numeric)),
          scale.unit = TRUE,
          graph = FALSE)

  # # Variance expliquée
  # fviz_eig(pca1)

  # # Individus
  # fviz_pca_ind(pca1,
  #             geom = "point",
  #             pointsize = 2,
  #             labels = NULL,)

  # # Variables
  # fviz_pca_var(pca1)

  # Biplot
  q[[spi]] <- fviz_pca_biplot(pca1, label = "var", psize = 2,
              col.ind = a_long_sp$Species)+
              labs(title = paste0("PCA of traits vs Htz_slope for ", spi, " species - dataset: ", j),
             color = "HTZ_slope")
              # ,
      #             repel = TRUE, 
      #             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"))+
      # labs(title = paste0("PCA of traits vs Htz_slope for ", spi, " species - dataset: ", j),
      #       color = "HTZ_slope")

Var <- as.data.frame(pca1$var$coord)*5
Var$var <- rownames(Var)
ind <- as.data.frame(pca1$ind$coord)
ind$Species <- a_long_sp$Species
ind$htz <- a_long_sp$htz_slope
ind <- ind[which(ind$Species %in% c("Whiting", "Hake", "Saithe", "AtlanticCod", "AtlanticMackerel", "GreyGurnard")),]
ind <- ind[-which(ind$Species %in% c("EuropeanSprat","Sandeel")),]


r <- ggplot() +  
  geom_point(data=ind, aes(Dim.1, Dim.2, color = Species), alpha=0.6) +
  geom_segment(data = Var, aes(x=0, y=0, xend = Dim.1, yend = Dim.2), arrow = arrow(length = unit(0.1, "inches"))
               )+  
  geom_text(data = Var, aes(x = Dim.1, y = Dim.2, label = var), size=3)+

  scale_color_viridis_d() +
  theme_minimal() +
  geom_abline(slope = 0, intercept = 0, linetype = "dashed", color = "grey") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey") +
  labs(title = paste0("PCA of traits vs Htz_slope - dataset: ", j),
       x = "PC1",
       y = "PC2")


  # #with factors
  # pca2 <- PCA(dat_merged_pca |> dplyr::select(where(is.numeric)),
  #         scale.unit = TRUE,
  #         graph = FALSE)

  # fviz_pca_ind(pca2,
  #             addEllipses = TRUE,
  #             ellipse.level = 0.95)

  # # by group
  # dat_merged_pca |> 
  # group_split(Species) |> 
  # map(~ PCA(select(.x, where(is.numeric)),
  #             scale.unit = TRUE,
  #             graph = FALSE))
}
pathoutputplots = paste0(configDir, "/output/PLOTS/Htz_vs_Ne_variables/plots")
pdf(paste0(pathoutputplots, "/_PCA_htz_vs_traits_", j, ".pdf"), width = 10, height = 8)
print(q)
dev.off()

pdf(paste0(pathoutputplots, "/SpeciesA_colored_PCA_htz_vs_traits_", j, ".pdf"), width = 10, height = 8)
print(r)
dev.off()

a_long
model1 <- NULL
modelmer <- NULL
model1 <- lme4::lmer(formula = as.formula(paste0("htz_slope", " ~",
        paste(colnames(a_long)[-which(colnames(dat_merged_w_sc) %in%
        c("Species", "replicate", "Fmax", "Fstart", "htz_slope"))], collapse="+"),
        "+ (1|Fmax) + (1|Fstart) + (1|replicate) + (1|Species)")), data=a_long)
modelmer <- buildmer(as.formula(paste0("htz_slope", " ~",
        paste(colnames(a_long)[-which(colnames(dat_merged_w_sc) %in%
        c("Species", "replicate", "Fmax", "Fstart", "htz_slope"))], collapse="+"),
        "+ (1|Fmax) + (1|Fstart) + (1|replicate) + (1|Species)")), data=a_long)




## PLOT PCA ---------------------------------- 


for (j in names(data_mean_sc)){
  P <- list()
  dat_merged_w_sc <- data_mean_sc[[j]]
  
  # ACP
  dat_merged_pca <- dat_merged_w_sc
  dat_merged_pca$Fmax <- as.numeric(dat_merged_pca$Fmax)
  dat_merged_pca$Fstart <- as.numeric(dat_merged_pca$Fstart)
  dat_merged_pca$replicate <- as.numeric(dat_merged_pca$replicate)
  dat_merged_pca <- dat_merged_pca[, -which(colnames(dat_merged_pca) %in% c("htz_slope", paste0("mean_", names(data_mean_sc))))]
  
  pca1 <- PCA(dat_merged_pca |> dplyr::select(where(is.numeric)),
              scale.unit = TRUE,
              graph = FALSE)
    
  Var <- as.data.frame(pca1$var$coord)*5
  Var$var <- rownames(Var)
  
  for (spi in unique(dat_merged_w_sc$Species)){
    ind <- as.data.frame(pca1$ind$coord)
    ind$Species <- dat_merged_pca$Species
    ind$htz <- dat_merged_w_sc$htz_slope
    ind[, paste0("mean", j)] <- dat_merged_w_sc[,paste0("mean_", j)]
    ind <- ind[which(ind$Species %in% spi),]
    
    P[[spi]] <- ggplot() +  
      stat_summary_hex(fun = mean, bins = 25, data=ind, 
                       aes(Dim.1, Dim.2, z=.data[[paste0("mean", j)]])) +
      geom_point(data=ind, 
                 aes(Dim.1, Dim.2, color = .data[[paste0("mean", j)]]), alpha=0.5) +
      
      geom_segment(data = Var, aes(x=0, y=0, xend = Dim.1, yend = Dim.2), arrow = arrow(length = unit(0.1, "inches"))
      )+  
      geom_text(data = Var, aes(x = Dim.1, y = Dim.2, label = var), size=3)+
      
      scale_color_viridis_c() +
      theme_minimal() +
      geom_abline(slope = 0, intercept = 0, linetype = "dashed", color = "grey") +
      geom_hline(yintercept = 0, linetype = "dashed", color = "grey") +
      labs(title = paste0("PCA of traits vs", j),
           subtitle = spi,
           x = paste0("PC1, ", round(pca1$eig[1,2], 2), "%"),
           y = paste0("PC2, ", round(pca1$eig[2,2], 2), "%"))
    
  }
  pathoutputplots = paste0(configDir, "/output/PLOTS/Htz_vs_Ne_variables/plots")
  pdf(paste0(pathoutputplots, "/one_genet_drift_PCA_htz_vs_traits_", j, ".pdf"), width = 10, height = 8)
  print(P)
  dev.off()
}



## PLOT PCA FRACTION ---------------------------------- 


for (j in names(data_mean_sc)){
  P <- list()
  dat_merged_w_sc <- data_mean_sc[[j]]
  # ACP
  dat_merged_pca <- dat_merged_w_sc
  dat_merged_pca$Fmax <- as.numeric(dat_merged_pca$Fmax)
  dat_merged_pca$Fstart <- as.numeric(dat_merged_pca$Fstart)
  dat_merged_pca$replicate <- as.numeric(dat_merged_pca$replicate)
  dat_merged_pca <- dat_merged_pca[, -which(colnames(dat_merged_pca) %in% c("htz_slope", paste0("mean_", names(data_mean_sc))))]
  
  colNe <- c(grep(colnames(dat_merged_w_sc), pattern="Ne"))
  if (length(colNe)==0){colNe <- grep(colnames(dat_merged_w_sc), pattern="Nb", ignore.case = "egg")}
  colNtot <- c(grep(colnames(dat_merged_w_sc), pattern="NT"))
  if (any(grepl("alpha", colnames(dat_merged_w_sc)))) colNtot <- c(grep(colnames(dat_merged_w_sc), pattern="Numerator_N_alpha"))
  if (any(grepl("1", colnames(dat_merged_w_sc))))     colNtot <- c(grep(colnames(dat_merged_w_sc), pattern="Numerator_N_1"))
  frac <- as.numeric(unlist(as.data.frame(dat_merged_w_sc[,colNe])/as.data.frame(dat_merged_w_sc[,colNtot])))
  
  pca1 <- PCA(dat_merged_pca |> dplyr::select(where(is.numeric)),
              scale.unit = TRUE,
              graph = FALSE)
  
  Var <- as.data.frame(pca1$var$coord)*5
  Var$var <- rownames(Var)
  
  for (spi in unique(dat_merged_w_sc$Species)){
    ind <- as.data.frame(pca1$ind$coord)
    ind$Species <- dat_merged_pca$Species
    ind$htz <- dat_merged_w_sc$htz_slope
    ind[, paste0("ratio", j)] <- frac
    ind <- ind[which(ind$Species %in% spi),]
    
    P[[spi]] <- ggplot() +  
      stat_summary_hex(fun = mean, bins = 25, data=ind, 
                       aes(Dim.1, Dim.2, z=.data[[paste0("ratio", j)]])) +
      geom_point(data=ind, 
                 aes(Dim.1, Dim.2, color = .data[[paste0("ratio", j)]]), alpha=0.5) +
      
      geom_segment(data = Var, aes(x=0, y=0, xend = Dim.1, yend = Dim.2), arrow = arrow(length = unit(0.1, "inches"))
      )+  
      geom_text(data = Var, aes(x = Dim.1, y = Dim.2, label = var), size=3)+
      
      scale_color_viridis_c() +
      theme_minimal() +
      geom_abline(slope = 0, intercept = 0, linetype = "dashed", color = "grey") +
      geom_hline(yintercept = 0, linetype = "dashed", color = "grey") +
      labs(title = paste0("PCA of traits vs ratio ", j),
           subtitle = spi,
           x = paste0("PC1, ", round(pca1$eig[1,2], 2), "%"),
           y = paste0("PC2, ", round(pca1$eig[2,2], 2), "%"))
    
  }
  pathoutputplots = paste0(configDir, "/output/PLOTS/Htz_vs_Ne_variables/plots")
  pdf(paste0(pathoutputplots, "/FRACone_genet_drift_PCA_htz_vs_traits_", j, "RATIO.pdf"), width = 10, height = 8)
  print(P)
  dev.off()
}


## PLOT DIFF ---------------------------------- 




for (j in names(data_mean_sc)){
  P <- list()
  dat_merged_w_sc <- diff_data[[j]]
  # ACP
  dat_merged_pca <- dat_merged_w_sc
  dat_merged_pca$Fmax <- as.numeric(dat_merged_pca$Fmax)
  dat_merged_pca$Fstart <- as.numeric(dat_merged_pca$Fstart)
  dat_merged_pca$replicate <- as.numeric(dat_merged_pca$replicate)
  dat_merged_pca_wide <- pivot_wider(dat_merged_pca, names_from = "variable", values_from = "value")
  dat_merged_pca <- dat_merged_pca_wide[, -which(colnames(dat_merged_pca_wide) %in% c("diff_htz_slope", paste0("diff_mean_", names(data_mean_sc))))]
  # colNe <- c(grep(colnames(dat_merged_pca), pattern="Ne"))
  # if (length(colNe)==0){colNe <- grep(colnames(dat_merged_pca), pattern="Nb", ignore.case = "egg")}
  # colNtot <- c(grep(colnames(dat_merged_pca), pattern="NT"))
  # if (any(grepl("alpha", colnames(dat_merged_pca)))) colNtot <- c(grep(colnames(dat_merged_pca), pattern="Numerator_N_alpha"))
  # if (any(grepl("1", colnames(dat_merged_pca))))     colNtot <- c(grep(colnames(dat_merged_pca), pattern="Numerator_N_1"))
  # dat_merged_pca$frac <- as.numeric(unlist(as.data.frame(dat_merged_pca[,colNe])/as.data.frame(dat_merged_pca[,colNtot])))
  # dat_merged_pca <- dat_merged_pca[, -c(colNe)]

  pca1 <- PCA(dat_merged_pca |> dplyr::select(where(is.numeric)),
              scale.unit = TRUE,
              graph = FALSE)

  Var <- as.data.frame(pca1$var$coord)*5
  Var$var <- rownames(Var)

  for (spi in unique(dat_merged_w_sc$Species)){
    ind <- as.data.frame(pca1$ind$coord)
    ind$Species <- dat_merged_pca$Species
    ind$htz <- dat_merged_pca_wide$diff_htz_slope
    ind[, paste0("mean", j)] <- dat_merged_pca_wide[,paste0("diff_mean_", j)]
    ind <- ind[which(ind$Species %in% spi),]

    P[[spi]] <- ggplot() +  
      stat_summary_hex(fun = mean, bins = 25, data=ind, 
                       aes(Dim.1, Dim.2, z=.data[[paste0("mean", j)]])) +
      geom_point(data=ind, 
                 aes(Dim.1, Dim.2, color = .data[[paste0("mean", j)]]), alpha=0.5) +
      
      geom_segment(data = Var, aes(x=0, y=0, xend = Dim.1, yend = Dim.2), arrow = arrow(length = unit(0.1, "inches"))
      )+  
      geom_text(data = Var, aes(x = Dim.1, y = Dim.2, label = var), size=3)+
      
      scale_color_viridis_c() +
      theme_minimal() +
      geom_abline(slope = 0, intercept = 0, linetype = "dashed", color = "grey") +
      geom_hline(yintercept = 0, linetype = "dashed", color = "grey") +
      labs(title = paste0("PCA of traits vs", j),
           subtitle = spi,
           x = paste0("PC1, ", round(pca1$eig[1,2], 2), "%"),
           y = paste0("PC2, ", round(pca1$eig[2,2], 2), "%"))
    
  }
  pathoutputplots = paste0(configDir, "/output/PLOTS/Htz_vs_Ne_variables/plots")
  pdf(paste0(pathoutputplots, "/DIFFone_genet_drift_PCA_htz_vs_traits_", j, ".pdf"), width = 10, height = 8)
  print(P)
  dev.off()
}






### PCA when dataset is relative + ratio


for (j in names(diffRATIO_data)){
  P <- list()
  dat_merged_w_sc <- diff_data[[j]]
  dat_merged_w_sc <- pivot_wider(dat_merged_w_sc, names_from = "variable", values_from = "value")
  # ACP
  dat_merged_pca <- dat_merged_w_sc
  dat_merged_pca$Fmax <- as.numeric(dat_merged_pca$Fmax)
  dat_merged_pca$Fstart <- as.numeric(dat_merged_pca$Fstart)
  dat_merged_pca$replicate <- as.numeric(dat_merged_pca$replicate)
  dat_merged_pca <- dat_merged_pca[, -which(colnames(dat_merged_pca) %in% 
                      c("diff_htz_slope", "diff_frac",paste0("diff_mean_", names(data_mean_sc))))]
  
  colNe <- c(grep(colnames(dat_merged_w_sc), pattern="Ne"))
  if (length(colNe)==0){colNe <- grep(colnames(dat_merged_w_sc), pattern="Nb", ignore.case = "egg")}
  colNtot <- c(grep(colnames(dat_merged_w_sc), pattern="NT"))
  if (any(grepl("alpha", colnames(dat_merged_w_sc)))) colNtot <- c(grep(colnames(dat_merged_w_sc), pattern="Numerator_N_alpha"))
  if (any(grepl("1", colnames(dat_merged_w_sc))))     colNtot <- c(grep(colnames(dat_merged_w_sc), pattern="Numerator_N_1"))
  frac <- as.numeric(unlist(as.data.frame(dat_merged_w_sc[,colNe])/as.data.frame(dat_merged_w_sc[,colNtot])))

  pca1 <- PCA(dat_merged_pca |> dplyr::select(where(is.numeric)),
              scale.unit = TRUE,
              graph = FALSE)

  Var <- as.data.frame(pca1$var$coord)*5
  Var$var <- rownames(Var)

  for (spi in unique(dat_merged_w_sc$Species)){
    ind <- as.data.frame(pca1$ind$coord)
    ind$Species <- dat_merged_pca$Species
    ind$htz <- dat_merged_w_sc$diff_htz_slope
    ind[, paste0("ratio_diff_", j)] <- frac
    ind <- ind[which(ind$Species %in% spi),]

    P[[spi]] <- ggplot() +  
      stat_summary_hex(fun = mean, bins = 25, data=ind, 
                       aes(Dim.1, Dim.2, z=.data[[paste0("ratio_diff_", j)]])) +
      geom_point(data=ind, 
                 aes(Dim.1, Dim.2, color = .data[[paste0("ratio_diff_", j)]]), alpha=0.5) +
      
      geom_segment(data = Var, aes(x=0, y=0, xend = Dim.1, yend = Dim.2), arrow = arrow(length = unit(0.1, "inches"))
      )+  
      geom_text(data = Var, aes(x = Dim.1, y = Dim.2, label = var), size=3)+
      
      scale_color_viridis_c() +
      theme_minimal() +
      geom_abline(slope = 0, intercept = 0, linetype = "dashed", color = "grey") +
      geom_hline(yintercept = 0, linetype = "dashed", color = "grey") +
      labs(title = paste0("PCA of traits vs ratio_diff_", j),
           subtitle = spi,
           x = paste0("PC1, ", round(pca1$eig[1,2], 2), "%"),
           y = paste0("PC2, ", round(pca1$eig[2,2], 2), "%"))

  }
  pathoutputplots = paste0(configDir, "/output/PLOTS/Htz_vs_Ne_variables/plots")
  pdf(paste0(pathoutputplots, "/DIFF_RATIOone_genet_drift_PCA_htz_vs_traits_", j, "RATIO.pdf"), width = 10, height = 8)
  print(P)
  dev.off()
}
