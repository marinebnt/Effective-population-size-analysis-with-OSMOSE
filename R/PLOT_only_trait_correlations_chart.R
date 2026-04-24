

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







pathoutputRDATA = paste0("osmose-ns-v4.4.1_clim/output/PLOTS/Htz_vs_Ne_variables/RDATA")



# Loop for Fishing Pressure -----------------------------------------------


pathout <- paste0(configDir, "/output/PLOTS/Htz_vs_Ne_variables")

# from RData output of the loops producing the plots of Trait by time and Genetic drift indicators, to df
# produce a dataframe with the right format : Ne variables and htz variables together
Time1=150
Time2=160
species_df=species
list_data_traits <- data_trait_and_genetic_into_df(Time1=150, Time2=160, pathoutput=paste0(configDir, "/output/PLOTS/RDATA"), 
                                                   FORCE = F, species = species_df)

# produce a dataframe with the right format : Ne variables and htz variables together
species_df=species
list_data_traits_ratio <- data_traitRATIO_and_genetic_into_df(Time1=150, Time2=160, pathoutput=paste0(configDir, "/output/PLOTS/RDATA"), 
                                                              FORCE = F, species = species_df)


# list of outputs 
list_outputs <- list("traits"=list_data_traits, "trait_ratio"=list_data_traits_ratio)


# Bind together the trait data coming from similar dimensions : cohort or year data together (by genetic drift dimension)
data_full_merged = NULL
pathoutputRDATA = paste0(configDir, "/output/PLOTS/RDATA")
if (!file.exists(paste0(pathoutputRDATA, "/PLOT_only_trait_correlations_chart.RData"))){
  for (f in names(list_outputs)[1]){
    
    genet_drift_data_full       <- list_outputs[[f]][["genet_drift_data"]]
    trait_data_full             <- list_outputs[[f]][["trait_data"]]
  
    for (SP in species){
      cat(SP)
      for (species_plotted in SP){
        cat(species_plotted, "\n")
        for (i in names(genet_drift_data_full)){
          cat("   ", i, "\n")
          # loading trait data      
          coho_data = TRUE
          ID_fold_y <- grep(folders_traits_y, pattern=species_plotted)
          ID_fold_c <- grep(folders_traits_c, pattern=species_plotted)
          if (length(ID_fold_c) ==0) { 
            cat("No trait file for cohort found for ", SP, "\n")
            coho_data=FALSE
          }
          # choosing trait data suited to compare genetic drift data with 
          if (any(grepl(i, pattern="Ne_N_1"),grepl(i, pattern="Ne_N_alpha"))) {
            if (!coho_data){
              #cat("No coho data, skipping ",genetic_drift_i, "\n")
              next}
            trait_file_c <- list.files(folders_traits_c[ID_fold_c], full.names = T) #names(t_c)[which(t_c==3)]
            trait_var_c <- as.tibble(get(load(trait_file_c)))
            trait_var_c <- trait_var_c[, -which(colnames(trait_var_c)=="var_F_for_sp")]
            trait_df <- trait_var_c
            if ("TRAIT" %in% colnames(trait_df)) trait_df$TRAIT <- paste0("cohort_", trait_df$TRAIT)
            if ("variable" %in% colnames(trait_df)) trait_df$variable <- paste0("cohort_", trait_df$TRAIT)
          } else {
            trait_file_y <- list.files(folders_traits_y[ID_fold_y], full.names = T) #names(t_c)[which(t_c==3)]
            trait_var_y <- get(load(trait_file_y)) 
            trait_df <- trait_var_y
            if ("TRAIT" %in% colnames(trait_df)) trait_df$TRAIT <- paste0("year_", trait_df$TRAIT)
            if ("variable" %in% colnames(trait_df)) trait_df$variable <- paste0("year_", trait_df$TRAIT)
          }
          if (any(trait_df$time==0)) trait_df <- trait_df[-which(trait_df$time==0),]
          if (length(trait_df$time)>0) {
            trait_df$cohort <- trait_df$time
            trait_df$time <- NULL
          }
          if (length(trait_df$rep)>0) {
            trait_df$replicate <- trait_df$rep
            trait_df$rep <- NULL
          }
          if (nrow(trait_df) ==0) {
            #at("No trait data for ", species_plotted, " with Fmax=", Fmax, " and Fstart=", Fstart, "\n")
            next
          }
          #cat("===========> here I am artificially patching an issue \n")
          trait_df$Fmax <-  unlist(stringr::str_split(unlist(stringr::str_split(trait_df$Fish_press, pattern="Fmax="))[seq(by=2, 2, 2*nrow(trait_df))], pattern = " "))[seq(by=2, 1, 2*nrow(trait_df))]
          trait_df$Fstart <-  unlist(stringr::str_split(trait_df$Fish_press, pattern="Fstart="))[seq(by=2, 2, 2*nrow(trait_df))]
          colnames(trait_df)[which(colnames(trait_df) == "TRAIT")] <- "variable"
          # when the species plotted is the same as the reference one for varying Fmax and Fstart, it means that the 
          trait_df <- trait_df[, -which(colnames(trait_df) %in% c("Fish_press"))]
          trait_df <- trait_df[which(trait_df$Species == species_plotted), ]
          trait_df <- trait_df[(as.numeric(trait_df$cohort) %in% Time1:Time2),]
          trait_df <- trait_df %>%
            group_by(replicate, variable, Fmax, Fstart) %>%
            summarise(across(where(is.numeric), mean, na.rm = TRUE)) %>%
            ungroup()
          trait_df$Species <- species_plotted
          
          # adding traits to the dataset : coming from Ne estimation process
          i_names_geneti_n <- ifelse(i %in% c("Ne_N_1", "Ne_N_alpha"), 1, ifelse(i %in% c("tempoNe"), 3,2))
          if (length(trait_data_full[[i_names_geneti_n]]$SP_ref)>1) trait_data_full[[i_names_geneti_n]]$SP_ref <- NULL
          trait_df <- rbind(trait_df, trait_data_full[[i_names_geneti_n]])
          if (!any(grepl(unique(trait_df$variable), pattern="htz"))) trait_df <- rbind(trait_df, trait_data_full$htz_list[, which(colnames(trait_data_full$htz_list) %in% colnames(trait_df))])
          if (i_names_geneti_n == "tempoNe") trait_df <- rbind(trait_df, trait_data_full[["data_Nb"]][, which(colnames(trait_data_full$data_Nb) == "NT")]) # because Nb and Netempo same reference : yearly
          # keeping only Ne data of interest
          genet_drift_di = genet_drift_data_full[[i]]
          genet_drift_di <- genet_drift_di[which(genet_drift_di$Species == species_plotted), ]
          if(nrow(genet_drift_di) == 0) next
          
          # attach together trait and genetic drift data
          trait_df <- trait_df[, which(colnames(trait_df) %in% colnames(genet_drift_di))]
          dat_merged <- rbind(genet_drift_di, trait_df)
          dat_merged_w <- pivot_wider(dat_merged, names_from = variable, values_from = value)
          
          # are there replicates ? 
          qq <- trait_df %>%
            dplyr::group_by(replicate, Fmax, Fstart, Species, variable) %>%
            dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
            dplyr::filter(n > 1L) 
          if (!is.null(qq)) stop_for_problems("duplicates")
          
          # remove missing columns 
          if (nrow(dat_merged_w) ==0) {
            cat("No merged data for ", species_plotted, " with Fmax=", Fmax, " and Fstart=", Fstart, "\n")
            next}
          
          if (is.null(data_full_merged[[i]])){
            data_full_merged[[i]] <- dat_merged_w
          } else {
            data_full_merged[[i]] <- rbind(data_full_merged[[i]], dat_merged_w)
          }

        } 
      }# sp plotted
    }#species
  } 
  save(data_full_merged, file=paste0(pathoutputRDATA, "/PLOT_only_trait_correlations_chart.RData"))
}else {
data_full_merged = get(load(paste0(pathoutputRDATA, "/PLOT_only_trait_correlations_chart.RData")))
}


# Dataset ready for modeling --------------------------------------------
# if (!file.exists(paste0(pathoutputRDATA, "/PLOT_only_trait_correlations_chart_SECONDPART.RData"))){

  data_full_merged -> data_mean

  # scaling as a standardized loop
  data_mean_sc <- list()
  for (j in names(data_mean)){
    data_mean_sc[[j]] <- data_mean[[j]]
    data_mean_sc[[j]][, -which(colnames(data_mean_sc[[j]]) %in%
                                 c("Species", "replicate", "Fmax", "Fstart"))] <- apply(data_mean_sc[[j]][, -which(colnames(data_mean_sc[[j]]) %in%
                                                                                                                     c("Species", "replicate", "Fmax", "Fstart"))], 2, function(x) scale(x))
    data_mean_sc[[j]][,as.numeric(names(table(which(is.na(data_mean_sc[[j]]), arr.ind=TRUE)[,2])))] <-
      lapply(data_mean_sc[[j]][,as.numeric(names(table(which(is.na(data_mean_sc[[j]]), arr.ind=TRUE)[,2])))],
             function(x) if(all(is.na(x))) NULL else x)
    data_mean_sc[[j]] <- na.omit(data_mean_sc[[j]])
  }

  save(data_mean_sc, file=paste0(pathoutputRDATA, "/PLOT_only_trait_correlations_chart_SECONDPART.RData"))
# }else {
#   data_mean_sc = get(load(paste0(pathoutputRDATA, "/PLOT_only_trait_correlations_chart_SECONDPART.RData")))
# }





#-------------------
# PLOT COrrelation plot ----------------------------------------


p <- list()
for (j in names(data_mean_sc)){

  dat_merged_w_sc <- data_mean_sc[[j]]
  dat_merged_w_sc$F_p <- paste0("Fmax=", dat_merged_w_sc$Fmax, "_Fstart=", dat_merged_w_sc$Fstart)


  for (spi in unique(dat_merged_w_sc$Species)){
    dat_merged_w_sc$Species <- as.character(dat_merged_w_sc$Species)
    dat_merged_w_sc_sp <- dat_merged_w_sc[which(dat_merged_w_sc$Species == spi), ]

    pathoutputplots = paste0(configDir, "/output/PLOTS/Htz_vs_Ne_variables/plots")
    pdf(paste0(pathoutputplots, "/CORR_htz_vs_traits_", spi, "_",j, ".pdf"), width = 10, height = 8)

    for (F_i in unique(dat_merged_w_sc_sp$F_p)){
      dat_merged_w_sc_sp_F <- dat_merged_w_sc_sp[which(dat_merged_w_sc_sp$F_p == F_i), ]
      my_data <- dat_merged_w_sc_sp_F[, -which(colnames(dat_merged_w_sc_sp_F) %in% c("Species", "Fmax", "Fstart", "replicate", "F_p"))]
      if (any(apply(my_data, 2, function(x) sum(!is.na(x))<3))) my_data <- my_data[, -which(apply(my_data, 2, function(x) sum(!is.na(x))<3))]
      if (ncol(my_data) < 2) next
      htz_col <- c(grep(colnames(my_data), pattern="htz"), grep(colnames(my_data), pattern="Ne"))
      chart.Correlation(R=as.data.frame(my_data[,union(htz_col, seq(1,floor(ncol(my_data)/2)))]), histogram=TRUE, se=T, pch=19, main=paste0("Correlation plot for ", spi, "\nat ", F_i))
    }
    for (F_i in unique(dat_merged_w_sc$F_p)){
      dat_merged_w_sc_sp_F <- dat_merged_w_sc[which(dat_merged_w_sc_sp$F_p == F_i), ]
      my_data <- dat_merged_w_sc_sp_F[, -which(colnames(dat_merged_w_sc_sp_F) %in% c("Species", "Fmax", "Fstart", "replicate", "F_p"))]
      if (any(apply(my_data, 2, function(x) sum(!is.na(x))<3))) my_data <- my_data[, -which(apply(my_data, 2, function(x) sum(!is.na(x))<3))]
      if (ncol(my_data) < 2) next
      htz_col <- c(grep(colnames(my_data), pattern="htz"), grep(colnames(my_data), pattern="Ne"))
      chart.Correlation(R=as.data.frame(my_data[,union(htz_col, seq(floor(ncol(my_data)/2), ncol(my_data)))]), histogram=TRUE, se=T, pch=19, main=paste0("Correlation plot for ", spi, "\nat ", F_i))
    }

    dev.off()
  }
}

# plist_correlation plot with ratios

p <- list()
for (j in names(data_mean_sc)){

  dat_merged_w_sc <- data_mean_sc[[j]]
  if (j == "Ne_temporal") dat_merged_w_sc <- left_join(dat_merged_w_sc, data_mean_sc[["Nb"]][, which(colnames(data_mean_sc[["Nb"]]) %in%
                                                                      c("Fmax", "Fstart", "replicate", "Species", "mean_NT"))])
  dat_merged_w_sc$F_p <- paste0("Fmax=", dat_merged_w_sc$Fmax, "_Fstart=", dat_merged_w_sc$Fstart)

  for (spi in unique(dat_merged_w_sc$Species)){
    dat_merged_w_sc$Species <- as.character(dat_merged_w_sc$Species)
    dat_merged_w_sc_sp <- dat_merged_w_sc[which(dat_merged_w_sc$Species == spi), ]
    colNe   <- c(grep(colnames(dat_merged_w_sc), pattern="Ne", ignore.case = "egg"), grep(colnames(dat_merged_w_sc), pattern="Nb"))
    colNtot <- c(grep(colnames(dat_merged_w_sc), pattern="NT"))
    if (any(grepl("alpha", colnames(dat_merged_w_sc_sp)))) colNtot <- c(grep(colnames(dat_merged_w_sc_sp), pattern="Numerator_N_alpha"))
    if (any(grepl("1", colnames(dat_merged_w_sc_sp))))     colNtot <- c(grep(colnames(dat_merged_w_sc_sp), pattern="Numerator_N_1"))
    dat_merged_w_sc_sp$frac <- as.numeric(unlist(as.data.frame(dat_merged_w_sc_sp[,colNe])/as.data.frame(dat_merged_w_sc_sp[,colNtot])))
    dat_merged_w_sc_sp <- dat_merged_w_sc_sp[, -c(colNe)]

    pathoutputplots = paste0(configDir, "/output/PLOTS/Htz_vs_Ne_variables/plots")
    pdf(paste0(pathoutputplots, "/CORR_FRAC_htz_vs_traits_", spi, "_",j, ".pdf"), width = 10, height = 8)

    for (F_i in unique(dat_merged_w_sc$F_p)){
      dat_merged_w_sc_sp_F <- dat_merged_w_sc[which(dat_merged_w_sc_sp$F_p == F_i), ]
      my_data <- dat_merged_w_sc_sp_F[, -which(colnames(dat_merged_w_sc_sp_F) %in% c("Species", "Fmax", "Fstart", "replicate", "F_p"))]
      if (any(apply(my_data, 2, function(x) sum(!is.na(x))<3))) my_data <- my_data[, -which(apply(my_data, 2, function(x) sum(!is.na(x))<3))]
      if (ncol(my_data) < 2) next
      htz_col <- c(grep(colnames(my_data), pattern="htz"), grep(colnames(my_data), pattern="Ne"))
      chart.Correlation(R=as.data.frame(my_data[,union(htz_col, seq(1,floor(ncol(my_data)/2)))]), histogram=TRUE, se=T, pch=19, main=paste0("Correlation plot for ", spi, "\nat ", F_i))
    }
    for (F_i in unique(dat_merged_w_sc$F_p)){
      dat_merged_w_sc_sp_F <- dat_merged_w_sc[which(dat_merged_w_sc_sp$F_p == F_i), ]
      my_data <- dat_merged_w_sc_sp_F[, -which(colnames(dat_merged_w_sc_sp_F) %in% c("Species", "Fmax", "Fstart", "replicate", "F_p"))]
      if (any(apply(my_data, 2, function(x) sum(!is.na(x))<3))) my_data <- my_data[, -which(apply(my_data, 2, function(x) sum(!is.na(x))<3))]
      if (ncol(my_data) < 2) next
      htz_col <- c(grep(colnames(my_data), pattern="htz"), grep(colnames(my_data), pattern="Ne"))
      chart.Correlation(R=as.data.frame(my_data[,union(htz_col, seq(floor(ncol(my_data)/2), ncol(my_data)))]), histogram=TRUE, se=T, pch=19, main=paste0("Correlation plot for ", spi, "\nat ", F_i))
    }

    dev.off()
  }
}


# plot_correlation plot with relative values (relative to F1-1)
# 
p <- list()
for (j in names(data_mean_sc)){

  dat_merged_w_sc <- data_mean_sc[[j]]
  if (j == "Ne_temporal") dat_merged_w_sc <- left_join(dat_merged_w_sc, data_mean_sc[["Nb"]][, which(colnames(data_mean_sc[["Nb"]]) %in%
                                                                                                       c("Fmax", "Fstart", "replicate", "Species", "mean_NT"))])
  dat_merged_w_sc$F_p <- paste0("Fmax=", dat_merged_w_sc$Fmax, "_Fstart=", dat_merged_w_sc$Fstart)


  for (spi in unique(dat_merged_w_sc$Species)){

    dat_merged_w_sc$Species <- as.character(dat_merged_w_sc$Species)
    dat_merged_w_sc_sp <- dat_merged_w_sc[which(dat_merged_w_sc$Species == spi), ]
    colNe   <- c(grep(colnames(dat_merged_w_sc), pattern="Ne", ignore.case = "egg"), grep(colnames(dat_merged_w_sc), pattern="Nb"))
    colNtot <- c(grep(colnames(dat_merged_w_sc), pattern="NT"))
    if (any(grepl("alpha", colnames(dat_merged_w_sc_sp)))) colNtot <- c(grep(colnames(dat_merged_w_sc_sp), pattern="Numerator_N_alpha"))
    if (any(grepl("1", colnames(dat_merged_w_sc_sp))))     colNtot <- c(grep(colnames(dat_merged_w_sc_sp), pattern="Numerator_N_1"))
    if (length(dat_merged_w_sc_sp[,colNtot])<1) {
      cat("dat_merged_w_sc_sp[,colNtot] is empty \n")
      next
    }
    dat_merged_w_sc_sp$frac <- as.numeric(unlist(as.data.frame(dat_merged_w_sc_sp[,colNe])/as.data.frame(dat_merged_w_sc_sp[,colNtot])))
    dat_merged_w_sc_sp <- dat_merged_w_sc_sp[, -c(colNe)]

    ref <- dat_merged_w_sc[which(dat_merged_w_sc_sp$F_p == paste0("Fmax=", 1, "_Fstart=", 1)), ]
    ref <-  ref[,-which(colnames(ref) %in% c('Fmax', "Fstart", "F_p"))]
    colnames(ref)[-which(colnames(ref) %in% c("replicate", 'Fmax', "Fstart", "Species", "F_p"))] <- paste0("ref_", colnames(ref[, -which(colnames(ref) %in% c("replicate", 'Fmax', "Fstart", "Species", "F_p"))]))
    pathoutputplots = paste0(configDir, "/output/PLOTS/Htz_vs_Ne_variables/plots")
    pdf(paste0(pathoutputplots, "/CORR_RELATIVE_htz_vs_traits_", spi, "_",j, ".pdf"), width = 10, height = 8)

    for (F_i in unique(dat_merged_w_sc$F_p)){

      # difference between reference data and other data ()
      dat_merged_w_sc_sp_F <- dat_merged_w_sc[which(dat_merged_w_sc_sp$F_p == F_i), ]
      dat_merged_w_sc_sp_F <- left_join(dat_merged_w_sc_sp_F, ref, by=c("Species", "replicate"))
      for (nm_ref in colnames(dat_merged_w_sc_sp_F)[-grep(colnames(dat_merged_w_sc_sp_F), pattern="ref_")]) {
        if (nm_ref %in% c("replicate", 'Fmax', "Fstart", "Species", "F_p")) next
        dat_merged_w_sc_sp_F[, paste0("diff", nm_ref)] <- dat_merged_w_sc_sp_F[, paste0("ref_", nm_ref)] - dat_merged_w_sc_sp_F[, nm_ref]
      }
      my_data <- dat_merged_w_sc_sp_F[, grep("diff", x = colnames(dat_merged_w_sc_sp_F))]
      if (any(apply(my_data, 2, function(x) sum(!is.na(x))<3))) my_data <- my_data[, -which(apply(my_data, 2, function(x) sum(!is.na(x))<3))]
      if (ncol(my_data) < 2) next
      htz_col <- c(grep(colnames(my_data), pattern="htz"), grep(colnames(my_data), pattern="Ne"))
      chart.Correlation(R=as.data.frame(my_data[,union(htz_col, seq(1,floor(ncol(my_data)/2)))]), histogram=TRUE, se=T, pch=19, main=paste0("Correlation plot for ", spi, "\nat ", F_i))
    }
    for (F_i in unique(dat_merged_w_sc$F_p)){
      dat_merged_w_sc_sp_F <- dat_merged_w_sc[which(dat_merged_w_sc_sp$F_p == F_i), ]
      my_data <- dat_merged_w_sc_sp_F[, -which(colnames(dat_merged_w_sc_sp_F) %in% c("Species", "Fmax", "Fstart", "replicate", "F_p"))]
      if (any(apply(my_data, 2, function(x) sum(!is.na(x))<3))) my_data <- my_data[, -which(apply(my_data, 2, function(x) sum(!is.na(x))<3))]
      if (ncol(my_data) < 2) next
      htz_col <- c(grep(colnames(my_data), pattern="htz"), grep(colnames(my_data), pattern="Ne"))
      chart.Correlation(R=as.data.frame(my_data[,union(htz_col, seq(floor(ncol(my_data)/2), ncol(my_data)))]), histogram=TRUE, se=T, pch=19, main=paste0("Correlation plot for ", spi, "\nat ", F_i))
    }

    dev.off()
  }
}


