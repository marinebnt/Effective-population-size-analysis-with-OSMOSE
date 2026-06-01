library(dplyr)
library(tidyr)
library(ggpubr)

main = "New_gpu_conf.R" # config file
source("osmose-ns-v4.4.1_clim/R/functions/source.R") # load config and key parameters and functions
source("C:/Users/mbeneat/Documents/osmose/Chap2/osmose-ns/osmose-ns-v4.4.1_clim/R/functions/auxiliary_for_lists.R")
source("C:/Users/mbeneat/Documents/osmose/Chap2/osmose-ns/osmose-ns-v4.4.1_clim/R/functions/functions_key_traits.R")


path <- c("C:/Users/mbeneat/Documents/osmose/Chap2/osmose-ns/osmose-ns-v4.4.1_clim/output/PLOTS/FINAL_Ne_Htz2")
folders = list.dirs(paste0(path, "/INIT/"), full.names = T, recursive = F)




## PLOT 0- Biomass vs Fmax & keep it on the side -------------------
pathoutputRDATA = paste0(configDir, "/output/PLOTS/KEY_TRAITS/RDATA")
load(file=paste0(pathoutputRDATA, "/Biomass_plots.RData"))
rel_biomassplot=list_bio$rel_biomassplot
biomassplot=list_bio$biomassplot


### produce dataset reading multiple .RData to make 1 dataframe, stored as a .RData file ----------
sp_i = "Hake"
if (!file.exists(paste0(path, "/Complete_INIT2", sp_i, ".RData"))){
  full_dataset = list()
  full_dataset_age = list()
  for (f_i in folders){
    cat(f_i, "\n")
    
    Fmax   <- stringr::str_split(stringr::str_split(f_i, pattern="Fmax")[[1]][2], pattern="_Fstart")[[1]][1]
    Fstart <- stringr::str_split(stringr::str_split(f_i, pattern="Fstart")[[1]][2], pattern="_sp")[[1]][1]
    
    if (Fstart != 1) next
    
    load(list.files(f_i, full.names = T))
    dataNe <- INIT2$data_Ne
    
    for (sp in names(dataNe)[1]){
      
      cat(sp, "\n")
      
      f_i_sp <- dataNe[[sp]]
      
      for (dim in names(f_i_sp)){
        
        f_i_sp_d <- f_i_sp[[dim]]
        
        if (is.data.frame(f_i_sp_d) || is.matrix(f_i_sp_d)){
          
          if (is.matrix(f_i_sp_d)) f_i_sp_d <- as.data.frame(f_i_sp_d)
          f_i_sp_d$cohort <- rownames(f_i_sp_d)
          f_i_sp_d <- pivot_longer(f_i_sp_d, names_to = "replicate", values_to = "value", cols=-cohort)
          f_i_sp_d$Species <- sp
          f_i_sp_d$variable <- dim
          f_i_sp_d$Fmax <- Fmax
          f_i_sp_d$Fstart <- Fstart
          
          if (is.null(full_dataset[[sp]])){
            full_dataset[[sp]] <- f_i_sp_d
          } else {
            full_dataset[[sp]] <- rbind(full_dataset[[sp]], f_i_sp_d)
          }
          
        } # if matrix or df
        if (!is.data.frame(f_i_sp_d) & !is.matrix(f_i_sp_d)){
          
          f_i_sp_d = lapply(names(f_i_sp_d), FUN = function(x) {
            f_i_sp_d[[x]] <- as.data.frame(f_i_sp_d[[x]])
            f_i_sp_d[[x]][, "replicate"] = stringr::str_remove(x, pattern="rep_")
            f_i_sp_d[[x]][, "age"] = rownames(f_i_sp_d[[x]])
            return(f_i_sp_d[[x]])
          })
          f_i_sp_d = do.call(rbind, f_i_sp_d)
          
          f_i_sp_d <- pivot_longer(f_i_sp_d, names_to = "cohort", values_to = "value", cols=-c(replicate, age))
          f_i_sp_d[, "cohort"] = stringr::str_remove(as.matrix(f_i_sp_d[, "cohort"]), pattern="cohort_")
          f_i_sp_d$Species <- sp
          f_i_sp_d$variable <- dim
          f_i_sp_d$Fmax <- stringr::str_split(stringr::str_split(f_i, pattern="Fmax")[[1]][2], pattern="_Fstart")[[1]][1]
          f_i_sp_d$Fstart <- stringr::str_split(stringr::str_split(f_i, pattern="Fstart")[[1]][2], pattern="_sp")[[1]][1]
          
          if (is.null(full_dataset_age[[sp]])){
            full_dataset_age[[sp]] <- f_i_sp_d
          } else {
            full_dataset_age[[sp]] <- rbind(full_dataset_age[[sp]], f_i_sp_d)
          }
          
        } # if list
      }# variable dimension 
    }# species dimension
  }# folder dimension (Fmax*Fstart)
  
  for (sp_i in names(full_dataset_age)){
     # MERGE raw variables without age dimension AND variables with age dimension --------------
    
   # produce a dataset without age dimension : weighted mean
    weight_df <- full_dataset_age[[sp_i]] %>%
      dplyr::filter(variable == "abundanceByAge") %>%
      dplyr::select(Species, replicate, Fstart, Fmax, cohort, age, weight_value = value)
    data_with_w <- full_dataset_age[[sp_i]] %>%
      dplyr::left_join(weight_df, by = c("Species", "replicate", "cohort", "age", "Fmax", "Fstart"), relationship = "many-to-many")
    data_with_w <- data_with_w[-which(data_with_w$weight_value == 0), ]
    data_with_w$age <- as.numeric(data_with_w$age)
    data_out <- data_with_w %>%
      dplyr::group_by(Species, replicate, Fstart, Fmax, cohort, variable) %>%
      dplyr::summarise(
        value = stats::weighted.mean(value, w = weight_value, na.rm = TRUE),
        .groups = "drop"
      )
    final_dataset <- rbind(full_dataset[[sp_i]], data_out)
    
   # produce a dataste with the age dimension
    full_dataset_add_age <- lapply(names(table(full_dataset_age[[sp_i]]$age)), FUN = function(x){
      out = as.data.frame(full_dataset[[sp_i]])
      out$age = x
      return(out)
    })
    full_dataset_add_age <- do.call(rbind, full_dataset_add_age)
    
    final_dataset_age <- rbind(full_dataset_add_age, full_dataset_age[[sp_i]])
   
  # save all 
    complete = list(final_dataset, final_dataset_age)
    
    save(complete, file=paste0(path, "/Complete_INIT2", sp_i, ".RData"))
  }
} else {
  cat("you have to load the dataset species by species")
}




# PLOTs : understand Vx variability --------------------
Time1=150
Time2=160
PATHOUT  = "osmose-ns-v4.4.1_clim/output/PLOTS/KEY_TRAITS/"
for (relative_to_ref in c(TRUE, FALSE)){
  # for the plot
  plot_biom = ifelse(relative_to_ref, "rel_biomassplot", "biomassplot")
  plot_biom = ifelse(relative_to_ref, "biomassplot", "biomassplot")
  
  for (sp in species){
    p = list()
    cat(sp)
    
    load(file=paste0(path, "/Complete_INIT2", sp, ".RData"))
    final_dataset_age = complete[[2]]
    final_dataset_age$cohort <- stringr::str_remove(final_dataset_age$cohort, "cohort_")
    final_dataset_age$replicate <- stringr::str_remove(final_dataset_age$replicate, "rep_")
    final_dataset_age <- final_dataset_age[which(final_dataset_age$cohort %in% Time1:Time2),]
    rm(complete)
    
    df_trait_sp <- final_dataset_age %>% filter(Species == sp)
    df_trait_sp <- df_trait_sp %>% filter(Fstart == 1)
    
    var_i    <- unique(df_trait_sp$variable)
    # compare all traits to their value /N1 to understand where is the increase coming from 
    for (v in var_i){
      if (v == "Numerator_N_1") {
        df_trait_t <- df_trait_sp[which(df_trait_sp$variable == "Numerator_N_1"), ]
        df_trait_t$value <- 0
        df_trait_t$variable <- paste0(v, "\nDivNumerator_N_1")
      } else{ 
        df_trait_t    <- estimateDiv(trait_data_full_i = df_trait_sp, var=paste0(v, "\nDivNumerator_N_1"))
      }
      df_trait_sp   <- rbind(df_trait_sp, df_trait_t)
    }
    var_i    <- unique(df_trait_sp$variable)
    
    ref_data <- df_trait_sp %>%
      filter(Fmax == 1 & Fstart == 1) %>%
      dplyr::select(variable, value, replicate, age, cohort) %>%
      dplyr::rename(ref_value = value)
    
    for (Fmax_i in c("0.5", "0.75", "1", "1.25", "1.5")){
      cat(Fmax_i, " ")
      df_trait_sp_a <- df_trait_sp %>% filter(Fmax == Fmax_i)
      
      if (sum(is.nan(df_trait_sp_a$value))>0) df_trait_sp_a <- df_trait_sp_a[-which(is.nan(df_trait_sp_a$value)), ]
      if (sum(is.infinite(df_trait_sp_a$value))>0) df_trait_sp_a <- df_trait_sp_a[-which(is.infinite(df_trait_sp_a$value)), ]
      
      # conditions to stop
      if (nrow(df_trait_sp_a)==0) next
      if (all(is.na(df_trait_sp_a$value) | df_trait_sp_a$value == 0)) {
        next
      }
      
      # Use the function
      if (length(var_i)>3){
        k=1
        id_to_plot = NULL
        for (i in 1:c((length(var_i)+1)/2)){
          
          # produce plot
          id_to_plot = c(k, ifelse((k+1)>length(var_i), 1, k+1))
          data_split <- df_trait_sp_a[which(df_trait_sp_a$variable %in% var_i[id_to_plot] ),]
          ref_data_split <- ref_data[which(ref_data$variable %in% var_i[id_to_plot]),]
          if (all(is.na(data_split$value))) next
          result <- plot_trait_vs_fmax(df_trait_sp = data_split, trait_name = var_i[id_to_plot],
                                       Time1, Time2, relative_to_ref = relative_to_ref, spec=sp, 
                                       age_dim=T, reference_input=ref_data, FORCElm=T)
          p[[sp]][[paste(var_i[id_to_plot], collapse="")]][[Fmax_i]]         <- result$plot
          k=k+2
        }
      }
    }#age
    
    # save plot
    suffix_relative <- paste0(ifelse(relative_to_ref, "_RELATIVE", ""))
    if (!dir.exists(paste0(PATHOUT, "/byage", sp))) dir.create(paste0(PATHOUT, "/byage", sp))
    if (is.null(p[[sp]])) next
    pdf(file = paste0(PATHOUT, "/byage", sp, "/VxRAWplot_", "data_Ne", sp, "vsFmax", Fmax_i, suffix_relative, ".pdf"), 
        width = 15, height = 8)
    for (trait_plotted in names(p[[sp]])){
      plot_list <- c(p[[sp]][[trait_plotted]], list(get(plot_biom)[[sp]]))
      plot_list <- lapply(plot_list, function(x) x + theme_minimal(base_size = 11) + theme(legend.position = "bottom"))
      combined_plot <- ggarrange(plotlist = plot_list, ncol = 3,
                                 nrow = ceiling((length(plot_list))/3),
                                 common.legend = TRUE, legend = "bottom")
      print(combined_plot)
    } 
    dev.off()
   

  }#sp
}


Time1=150
Time2=160
traits_sp_to_plot = c("Deltasquare\nDivNumerator_N_1", "Vk_point", 
                      "k_point", "kx_bar", 
                      "k_point\nDivNumerator_N_1", "kx_bar\nDivNumerator_N_1", 
                      "DeadSum\nDivNumerator_N_1", "Vk_point",
                      "colSumSSDix", "colSumSSDgx", 
                      "colSumSSDix\nDivNumerator_N_1", "colSumSSDgx\nDivNumerator_N_1")
PATHOUT  = "osmose-ns-v4.4.1_clim/output/PLOTS/KEY_TRAITS/"
for (relative_to_ref in c(TRUE, FALSE)){
  # for the plot
  plot_biom = ifelse(relative_to_ref, "rel_biomassplot", "biomassplot")
  plot_biom = ifelse(relative_to_ref, "biomassplot", "biomassplot")
  
  for (sp in species){
    p = list()
    cat(sp)
    
    load(file=paste0(path, "/Complete_INIT2", sp, ".RData"))
    final_dataset_age = complete[[2]]
    final_dataset_age$cohort <- stringr::str_remove(final_dataset_age$cohort, "cohort_")
    final_dataset_age$replicate <- stringr::str_remove(final_dataset_age$replicate, "rep_")
    final_dataset_age <- final_dataset_age[which(final_dataset_age$cohort %in% Time1:Time2),]
    rm(complete)
    
    df_trait_sp <- final_dataset_age %>% filter(Species == sp)
    df_trait_sp <- df_trait_sp %>% filter(Fstart == 1)
    
    #add delta square
    df_trait_t    <- estimateSquare(trait_data_full_i = df_trait_sp, var=paste0("Deltasquare"))
    df_trait_sp   <- rbind(df_trait_sp, df_trait_t)
    
    #add sum by age of SSDgx and SSDix
    df_trait_t    <- estimatecolSum(trait_data_full_i = df_trait_sp, var=paste0("colSumSSDix"))
    df_trait_sp   <- rbind(df_trait_sp, df_trait_t)
    df_trait_t    <- estimatecolSum(trait_data_full_i = df_trait_sp, var=paste0("colSumSSDgx"))
    df_trait_sp   <- rbind(df_trait_sp, df_trait_t)
    
    var_i    <- unique(df_trait_sp$variable)
    # compare all traits to their value /N1 to understand where is the increase coming from 
    for (v in var_i){
      if (v == "Numerator_N_1") {
        df_trait_t <- df_trait_sp[which(df_trait_sp$variable == "Numerator_N_1"), ]
        df_trait_t$value <- 0
        df_trait_t$variable <- paste0(v, "\nDivNumerator_N_1")
      } else{ 
        df_trait_t    <- estimateDiv(trait_data_full_i = df_trait_sp, var=paste0(v, "\nDivNumerator_N_1"))
      }
      df_trait_sp   <- rbind(df_trait_sp, df_trait_t)
    }
    var_i    <- unique(df_trait_sp$variable)
    
    ref_data <- df_trait_sp %>%
      filter(Fmax == 1 & Fstart == 1) %>%
      dplyr::select(variable, value, replicate, age, cohort) %>%
      dplyr::rename(ref_value = value)
    
    for (Fmax_i in c("0.5", "0.75", "1", "1.25", "1.5")){
      cat(Fmax_i, " ")
      df_trait_sp_a <- df_trait_sp %>% filter(Fmax == Fmax_i)
      
      if (sum(is.nan(df_trait_sp_a$value))>0) df_trait_sp_a <- df_trait_sp_a[-which(is.nan(df_trait_sp_a$value)), ]
      if (sum(is.infinite(df_trait_sp_a$value))>0) df_trait_sp_a <- df_trait_sp_a[-which(is.infinite(df_trait_sp_a$value)), ]
      
      # conditions to stop
      if (nrow(df_trait_sp_a)==0) next
      if (all(is.na(df_trait_sp_a$value) | df_trait_sp_a$value == 0)) {
        next
      }
      
      # Use the function
      var_i = traits_sp_to_plot 
      if (length(var_i)>3){
        k=1
        id_to_plot = NULL
        for (i in 1:c((length(var_i)+1)/2)){
          
          # produce plot
          id_to_plot = c(k, ifelse((k+1)>length(var_i), 1, k+1))
          data_split <- df_trait_sp_a[which(df_trait_sp_a$variable %in% var_i[id_to_plot] ),]
          ref_data_split <- ref_data[which(ref_data$variable %in% var_i[id_to_plot]),]
          if (all(is.na(data_split$value))) next
          result <- plot_trait_vs_fmax(df_trait_sp = data_split, trait_name = var_i[id_to_plot],
                                       Time1, Time2, relative_to_ref = relative_to_ref, spec=sp, 
                                       age_dim=T, reference_input=ref_data, FORCElm=T)
          p[[sp]][[paste(var_i[id_to_plot], collapse="")]][[Fmax_i]]         <- result$plot
          k=k+2
        }
      }
    }#age
    
    # save plot
    suffix_relative <- paste0(ifelse(relative_to_ref, "_RELATIVE", ""))
    if (!dir.exists(paste0(PATHOUT, "/byage", sp))) dir.create(paste0(PATHOUT, "/byage", sp))
    if (is.null(p[[sp]])) next
    pdf(file = paste0(PATHOUT, "/byage", sp, "/VxRAW2plot_", "data_Ne", sp, "vsFmax", Fmax_i, suffix_relative, ".pdf"), 
        width = 15, height = 8)
    for (trait_plotted in names(p[[sp]])){
      plot_list <- c(p[[sp]][[trait_plotted]], list(get(plot_biom)[[sp]]))
      plot_list <- lapply(plot_list, function(x) x + theme_minimal(base_size = 11) + theme(legend.position = "bottom"))
      combined_plot <- ggarrange(plotlist = plot_list, ncol = 3,
                                 nrow = ceiling((length(plot_list))/3),
                                 common.legend = TRUE, legend = "bottom")
      print(combined_plot)
    } 
    dev.off()
    
    
  }#sp
}




# age delta²/N and Death by age/N
Time1=150
Time2=160
traits_sp_to_plot = c("Deltasquare\nDivNumerator_N_1", 
                      "Dead\nDivNumerator_N_1")
PATHOUT  = "osmose-ns-v4.4.1_clim/output/PLOTS/KEY_TRAITS/"
for (relative_to_ref in c(TRUE, FALSE)){
  # for the plot
  plot_biom = ifelse(relative_to_ref, "rel_biomassplot", "biomassplot")
  plot_biom = ifelse(relative_to_ref, "biomassplot", "biomassplot")
  
  for (sp in species){
    p = list()
    cat(sp)
    
    load(file=paste0(path, "/Complete_INIT2", sp, ".RData"))
    final_dataset_age = complete[[2]]
    final_dataset_age$cohort <- stringr::str_remove(final_dataset_age$cohort, "cohort_")
    final_dataset_age$replicate <- stringr::str_remove(final_dataset_age$replicate, "rep_")
    final_dataset_age <- final_dataset_age[which(final_dataset_age$cohort %in% Time1:Time2),]
    # rm(complete)
    
    df_trait_sp <- final_dataset_age %>% filter(Species == sp)
    df_trait_sp <- df_trait_sp %>% filter(Fstart == 1)
    
    #add delta square
    df_trait_t    <- estimateSquare(trait_data_full_i = df_trait_sp, var=paste0("Deltasquare"))
    df_trait_sp   <- rbind(df_trait_sp, df_trait_t)
    
    #add sum by age of SSDgx and SSDix
    df_trait_t    <- estimatecolSum(trait_data_full_i = df_trait_sp, var=paste0("colSumSSDix"))
    df_trait_sp   <- rbind(df_trait_sp, df_trait_t)
    df_trait_t    <- estimatecolSum(trait_data_full_i = df_trait_sp, var=paste0("colSumSSDgx"))
    df_trait_sp   <- rbind(df_trait_sp, df_trait_t)
    
    var_i    <- unique(df_trait_sp$variable)
    # compare all traits to their value /N1 to understand where is the increase coming from 
    for (v in var_i){
      if (v == "Numerator_N_1") {
        df_trait_t <- df_trait_sp[which(df_trait_sp$variable == "Numerator_N_1"), ]
        df_trait_t$value <- 0
        df_trait_t$variable <- paste0(v, "\nDivNumerator_N_1")
      } else{ 
        df_trait_t    <- estimateDiv(trait_data_full_i = df_trait_sp, var=paste0(v, "\nDivNumerator_N_1"))
      }
      df_trait_sp   <- rbind(df_trait_sp, df_trait_t)
    }
    var_i    <- unique(df_trait_sp$variable)
    
    ref_data <- df_trait_sp %>%
      filter(Fstart == 1 & Fmax==1) %>%
      dplyr::select(variable, value, replicate, age, cohort) %>%
      dplyr::rename(ref_value = value)
    df_trait_sp_a <- df_trait_sp
      
    if (sum(is.nan(df_trait_sp_a$value))>0) df_trait_sp_a <- df_trait_sp_a[-which(is.nan(df_trait_sp_a$value)), ]
    if (sum(is.infinite(df_trait_sp_a$value))>0) df_trait_sp_a <- df_trait_sp_a[-which(is.infinite(df_trait_sp_a$value)), ]
    
    # conditions to stop
    if (nrow(df_trait_sp_a)==0) next
    if (all(is.na(df_trait_sp_a$value) | df_trait_sp_a$value == 0)) {
      next
    }

    
    # produce plot
    var_i = traits_sp_to_plot 
    df_trait_sp_trait <- df_trait_sp_a[which(df_trait_sp_a$variable %in% var_i),]
    ref_data_split <- ref_data[which(ref_data$variable %in% var_i),]
    if (all(is.na(data_split$value))) next
    for (t in var_i){
      result <- plot_trait_facet_fmax(df_trait_sp = df_trait_sp_trait[which(df_trait_sp_trait$variable == t),], trait_name = t,
                                      Time1, Time2, relative_to_ref = relative_to_ref, spec=sp, 
                                      age_dim=T, reference_input=ref_data[which(ref_data$variable == t),], FORCElm=T, title = "simplified")
      p[[sp]][[t]][[Fmax_i]]         <- result$plot
    }
    
    
    
    # save plot
    suffix_relative <- paste0(ifelse(relative_to_ref, "_RELATIVE", ""))
    if (!dir.exists(paste0(PATHOUT, "/byage", sp))) dir.create(paste0(PATHOUT, "/byage", sp))
    # pdf(file = paste0(PATHOUT, "/byage", sp, "/VxRAW2plot_", "data_Ne", sp, "vsFmax", Fmax_i, suffix_relative, ".pdf"), 
        # width = 15, height = 8)
    for (trait_plotted in names(p[[sp]])){
      plot_list <- c(p[[sp]][[trait_plotted]])
      plot_list <- lapply(plot_list, function(x) x + theme_minimal(base_size = 11) + theme(legend.position = "bottom"))
      combined_plot[[trait_plotted]] <- ggarrange(plotlist = plot_list, ncol = 3,
                                 nrow = ceiling((length(plot_list))/3),
                                 common.legend = TRUE, legend = "bottom")
      print(combined_plot[[trait_plotted]])
    } 
    # dev.off()
  }#sp
}
save(combined_plot, file=paste0(pathoutputRDATA, "/deaths_deltasquare_byage.RData"))


### PLOT : mortality decomposition --------------------------

source("osmose-ns-v4.4.1_clim/R/functions/function_data_traits_summary.R")
pathout <- paste0("osmose-ns-v4.4.1_clim/output/PLOTS/Trait_variables_by_cohort_by_age2/otherplots")
if (!dir.exists(pathout)) dir.create(pathout, recursive = T)

# load dataset
data_complete_sp = NULL
data_complete_c = NULL
base_path = "osmose-ns-v4.4.1_clim/output/PLOTS/Trait_variables_by_cohort_by_age2"
folder_name <- basename(base_path)
suffix <- gsub("Trait_variables_by_", "", folder_name)
suffix <- gsub("_", "_", suffix)  # keep underscores
save_dir <- file.path(dirname(dirname(base_path)), "RDATA")
save_file <- file.path(save_dir, paste0("summary_trait_", suffix, ".RData"))
if (!dir.exists(save_dir)){dir.create(save_dir, recursive = T)}
if (!file.exists(save_file)){
  data_complete_c <- data_traits_summary(base_path)
} else {
  data_complete_c = get(load(file=save_file))
}

# dimensions dataset required : [time, Sp, value, rep, Fish_press, TRAIT]
head(data_complete_c)
data_complete_c <- data_complete_c[, which(colnames(data_complete_c) %in% c("age", "time", "Species", "value", "TRAIT", "replicate", "Fish_press"))]
data_complete_c[grep(data_complete_c$TRAIT, pattern="fec"), "TRAIT"] <- "mean_of_eggNumber_per_age"
unique(data_complete_c$TRAIT)

# adding natural mortality dimension
data_complete_c[which((data_complete_c$TRAIT) == "mean_of_Mnat_per_age"), "TRAIT"] <- "mean_of_MaddbyY"
data_complete_c <- pivot_wider(data = data_complete_c, names_from = "TRAIT", values_from = "value")
data_complete_c <- data_complete_c[,-grep(colnames(data_complete_c), pattern="mean_of_Dis_per_age")]
data_complete_c <- data_complete_c[,-grep(colnames(data_complete_c), pattern="mean_of_Z_per_age")]
data_complete_c$mean_of_MnatbyY <- data_complete_c$mean_of_MaddbyY + data_complete_c$mean_of_Mpred_per_age + data_complete_c$mean_of_F_per_age + data_complete_c$mean_of_Mstar_per_age
data_complete_c <- pivot_longer(data = data_complete_c, names_to = "TRAIT", values_to = "value", cols=-colnames(data_complete_c)[1:4])

# changing colnames to match expectations
colnames(data_complete_c)[which(colnames(data_complete_c) == "TRAIT")] <- "variable"
data_complete_c$Fmax   <- as.numeric(unlist(stringr::str_remove(unlist(stringr::str_split(data_complete_c$Fish_press, " Fstart"))[seq(1, 2*length(data_complete_c$Fish_press), by=2)], "Fmax=")))
data_complete_c$Fstart <- as.numeric(unlist(stringr::str_split(data_complete_c$Fish_press, " Fstart="))[seq(2, 2*length(data_complete_c$Fish_press), by=2)])
data_complete_c <- (data_complete_c[,-which(colnames(data_complete_c) == "Fish_press")])



PATHOUT  = "osmose-ns-v4.4.1_clim/output/PLOTS/KEY_TRAITS/"
for (relative_to_ref in c(TRUE, FALSE)){
  # for the plot
  plot_biom = ifelse(relative_to_ref, "rel_biomassplot", "biomassplot")
  plot_biom = ifelse(relative_to_ref, "biomassplot", "biomassplot")
  
  for (sp in species){
    p = list()
    cat(sp)
    
    df_trait_sp <- data_complete_c %>% filter(Species == sp)
    df_trait_sp <- df_trait_sp %>% filter(Fstart == 1)
    
    df_trait_t    <- estimateSum(trait_data_full_i = df_trait_sp, var=paste0("mean_of_F_per_age", "\nsum", "mean_of_MnatbyY"))
    df_trait_t$variable    <- "sum_mortality_rates"
    df_trait_sp   <- rbind(df_trait_sp, df_trait_t)
    
    var_i    <- unique(df_trait_sp$variable)
    # compare all mortality rates to their value *abundance to have the number of individuals to understand where is the increase coming from 
    for (v in var_i){
      if (any(grepl(pattern="mean_of_M", v) | grepl(pattern="mean_of_F", v) | grepl(pattern="sum_mortality_rates", v))) {
        df_trait_t    <- estimateX(trait_data_full_i = df_trait_sp, var=paste0(v, "\nx", "sum_of_abundanceByAge_per_age"))
        df_trait_sp   <- rbind(df_trait_sp, df_trait_t)
      }
    }
    var_i    <- unique(df_trait_sp$variable)
    
    
    ref_data <- df_trait_sp %>%
      filter(Fmax == 1 & Fstart == 1) %>%
      dplyr::select(variable, value, replicate, age) %>%
      dplyr::rename(ref_value = value)
    
    for (Fmax_i in c("0.5", "0.75", "1", "1.25", "1.5")){
      cat(Fmax_i, " ")
      df_trait_sp_a <- df_trait_sp %>% filter(Fmax == Fmax_i)
      
      if (sum(is.nan(df_trait_sp_a$value))>0) df_trait_sp_a <- df_trait_sp_a[-which(is.nan(df_trait_sp_a$value)), ]
      if (sum(is.infinite(df_trait_sp_a$value))>0) df_trait_sp_a <- df_trait_sp_a[-which(is.infinite(df_trait_sp_a$value)), ]
      
      # conditions to stop
      if (nrow(df_trait_sp_a)==0) next
      if (all(is.na(df_trait_sp_a$value) | df_trait_sp_a$value == 0)) {
        next
      }
      
      # Use the function
      if (length(var_i)>3){
        k=1
        id_to_plot = NULL
        for (i in 1:c((length(var_i)+1)/2)){
          
          # produce plot
          id_to_plot = c(k, ifelse((k+1)>length(var_i), 1, k+1))
          data_split <- df_trait_sp_a[which(df_trait_sp_a$variable %in% var_i[id_to_plot] ),]
          ref_data_split <- ref_data[which(ref_data$variable %in% var_i[id_to_plot]),]
          if (all(is.na(data_split$value))) next
          result <- plot_trait_vs_fmax(df_trait_sp = data_split, trait_name = var_i[id_to_plot],
                                       Time1=150, Time2=160, relative_to_ref = relative_to_ref, spec=sp, 
                                       age_dim=T, reference_input=ref_data)
          p[[sp]][[paste(var_i[id_to_plot], collapse="")]][[Fmax_i]]         <- result$plot
          k=k+2
        }
      }
    }#age
    
    # save plot
    suffix_relative <- paste0(ifelse(relative_to_ref, "_RELATIVE", ""))
    if (!dir.exists(paste0(PATHOUT, "/byage", sp))) dir.create(paste0(PATHOUT, "/byage", sp))
    if (is.null(p[[sp]])) next
    pdf(file = paste0(PATHOUT, "/byage", sp, "/VxMortal_", "data_Ne", sp, "vsFmax", Fmax_i, suffix_relative, ".pdf"), 
        width = 15, height = 8)
    for (trait_plotted in names(p[[sp]])){
      plot_list <- c(p[[sp]][[trait_plotted]], list(get(plot_biom)[[sp]]))
      plot_list <- lapply(plot_list, function(x) x + theme_minimal(base_size = 11) + theme(legend.position = "bottom"))
      combined_plot <- ggarrange(plotlist = plot_list, ncol = 3,
                                 nrow = ceiling((length(plot_list))/3),
                                 common.legend = TRUE, legend = "bottom")
      print(combined_plot)
    } 
    dev.off()
    
    
  }#sp
}


