library(ggplot2)
library(gridExtra)   # pour arranger plusieurs plots
library(patchwork)   # alternative moderne (recommandée)
library(stringr)
library(dplyr)
library(GGally)
library(tidyverse)
library(lme4)
library("PerformanceAnalytics")
library(legendry)
library(ggpubr)

format_slope_subtitle <- function(df, x_col = "Fmax", y_col = "mean", group_col = "variable") {
  df <- df[!is.na(df[[x_col]]) & !is.na(df[[y_col]]), ]
  if (nrow(df) < 2) return("Slope (lm mean ~ Fmax): NA")
  groups <- unique(df[[group_col]])
  slope_texts <- character(0)
  for (g in groups) {
    dfg <- df[df[[group_col]] == g, ]
    if (nrow(dfg) < 2 || length(unique(dfg[[x_col]])) < 2) {
      slope <- NA
      p-val <- NA
    } else {
      slope <- coef(lm(dfg[[y_col]] ~ dfg[[x_col]]))[2]
      p_val <- summary(lm(dfg[[y_col]] ~ dfg[[x_col]]))$coefficients[2,4]
    }
    slope_texts <- c(slope_texts, paste0(g, "=", formatC(slope, digits = 3, format = "f"), "\n(p-val=", formatC(p_val, digits = 3, format = "f"), ")"))
  }
  paste0("\nSlope (lm mean ~ Fmax): ", paste(slope_texts, collapse = "; "))
}

# Function to plot traits vs Fmax with optional relative-to-reference option
plot_trait_vs_fmax <- function(df_trait_sp, pathout, trait_name, Time1, Time2, relative_to_ref = FALSE, spec) {
  # Extract reference (Fmax=1, Fstart=1) if relative mode is on
  df_trait_sp_F        <- df_trait_sp[which(df_trait_sp$Species %in% spec),]
  df_trait_sp_F$Fstart <- as.numeric(df_trait_sp_F$Fstart)
  df_trait_sp_F$Fmax   <- as.numeric(df_trait_sp_F$Fmax)
  if (relative_to_ref) {
    ref_data <- df_trait_sp_F %>%
      filter(Fmax == 1 & Fstart == 1) %>%
      dplyr::select(variable, value, replicate) %>%
      dplyr::rename(ref_value = value)
    ref_data <- na.omit(ref_data)
    if (nrow(ref_data)>1)  {
      # Subtract reference from all means
      df_trait_sp_F <- df_trait_sp_F %>%
        left_join(ref_data, by = c("replicate", "variable"), relationship = "many-to-one") %>%
        mutate(value =(value - ref_value)/ref_value) %>%
        dplyr::select(-ref_value)
    }
  }
  df_trait_sp_F <- df_trait_sp_F %>%
    group_by(Fmax, Fstart, variable) %>%
    summarise(mean = mean(value, na.rm = TRUE),
              sd   = sd(value, na.rm = TRUE))
  
  scale = 1
  var_scale_regular = NA
  
  mean_var <- aggregate(df_trait_sp_F[, "mean"], list(df_trait_sp_F$variable), mean)
  var_scale_regular = mean_var[which(mean_var[, 2] > (0.8 * (max(mean_var[, 2])))), 1]
  if (length(var_scale_regular) > 0) scale = mean(mean_var[-which.max(mean_var[, 2]), "mean"], na.rm = T) / quantile(mean(mean_var[which.max(mean_var[, 2]), "mean"], na.rm = T), 0.7)
  if (is.na(scale)) scale = 1
  if (relative_to_ref & abs(scale) < 0.01) scale = sign(scale)*0.05
  if (length(unique(df_trait_sp_F$variable)) > 1) {
      df_trait_sp_F[-which(df_trait_sp_F$variable %in% var_scale_regular), "mean"] <-
      df_trait_sp_F[-which(df_trait_sp_F$variable %in% var_scale_regular), "mean"] / scale
  }
  
  df_trait_sp_F$Fmax     = as.numeric(df_trait_sp_F$Fmax)
  df_trait_sp_F$variable = as.factor(df_trait_sp_F$variable)
  df_trait_sp_F$Fstart   = as.factor(df_trait_sp_F$Fstart)
  df_trait_sp_F$mean     = as.numeric(df_trait_sp_F$mean)

  slope_subtitle <- format_slope_subtitle(df = df_trait_sp_F, x_col = "Fmax", y_col = "mean", group_col = "variable")
  
  Fstart_names <- as_labeller(c(
    `0.5` = "Fstart 0.5",
    `0.75` = "Fstart 0.75",
    `1` = "Fstart 1",
    `1.25` = "Fstart 1.25",
    `1.5` = "Fstart 1.5"
  ))
  z_axis <- mean_var[ifelse(any(mean_var[, 2] >(0.8 * (max(mean_var[, 2])))), which((mean_var[, 2] <=(0.8 * (max(mean_var[, 2]))))), 0), 1]
  if (length(z_axis)>0) {name_y <- mean_var[,1][-which(mean_var[,1] == z_axis)]
  }else  name_y <- paste(mean_var[,1] , collapse=" & ")
  p <- ggplot(df_trait_sp_F, aes(x = Fmax, y = mean,
    colour = variable, shape = Fstart, group = interaction(Fstart, variable)
  )) +
    geom_smooth(method = "loess", se = TRUE,
      aes(group = interaction(variable, Fstart))
    ) +
    geom_point(size = 2) +
    theme_bw() +
    facet_wrap(~Fstart, labeller = Fstart_names) +
    labs(
      title = paste0(ifelse(relative_to_ref, "Relative ", ""), paste(trait_name, collapse =" & "), " over ", Time1, " and ", Time2),
      x = "Fmax",
      y =  name_y, 
      subtitle = paste(spec, slope_subtitle, sep = " | ")
    ) +
    theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
    scale_y_continuous(
      sec.axis = sec_axis(~. * scale, 
        name = z_axis
      )
    )
  
  return(list(plot = p, scale = scale))
}

# Function to plot traits with secondary axis scaling
plot_trait_vs_fmax_simple <- function(df_trait_sp, pathout, sp_name, var_names, Time1, Time2, relative_to_ref = FALSE, scale_factor = 1) {
  # Extract reference (Fmax=1, Fstart=1) if relative mode is on
  if (relative_to_ref) {
    ref_data <- df_trait_sp %>%
      filter(Fmax == 1 & Fstart == "1") %>%
      dplyr::select(variable, mean) %>%
      dplyr::rename(ref_mean = mean)
    
    # Subtract reference from all means
    df_trait_sp <- df_trait_sp %>%
      left_join(ref_data[, -which(colnames(ref_data) %in% c("Fmax", "Fstart"))], by = "variable", relationship = "many-to-one") %>%
      mutate(mean =(mean - ref_mean)/ref_mean) %>%
      dplyr::select(-ref_mean)
  }

  
  p <- ggplot(df_trait_sp, aes(x = Fmax, group = variable, colour=variable)) +
    geom_ribbon(aes(ymin = mean - sd, ymax = mean + sd),
      alpha = 0.5, fill = "grey70"
    ) +
    geom_smooth(method = "loess", se = TRUE,
                aes(group = interaction(variable, Fstart), y=mean)
    ) +
    geom_point(aes(pch = as.factor(variable), y = mean)) +
    theme_bw() +
    labs(
      title = paste0(paste(var_names, collapse = " and "), " mean\nover ", Time1, " and ", Time2),
      x = "Fmax",
      y = ifelse(relative_to_ref, paste(var_names, "- RELATIVE", collapse = "\n"), paste(var_names, collapse = "\n")),
      subtitle = sp_name
    ) +
    theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
    scale_y_continuous(
      sec.axis = sec_axis(~. * scale_factor, name = "T_Vk")
    )
  
  return(p)
}

main = "New_gpu_conf.R"
source("osmose-ns-v4.4.1_clim/R/functions/source.R")
base_dir <- "osmose-ns-v4.4.1_clim/output/PLOTS/FINAL_Ne_Htz/INIT"
folders  <- list.dirs(base_dir, recursive = F)
INIT2 <- get(load(list.files(folders[1], full.names = T)))

folders_traits_y <- "osmose-ns-v4.4.1_clim/output/PLOTS/Trait_variables_by_year2/RData/"
folders_traits_y <- list.dirs(folders_traits_y, recursive = F)

folders_traits_c <- "osmose-ns-v4.4.1_clim/output/PLOTS/Trait_variables_by_cohort2/RData/"
folders_traits_c  <- list.dirs(folders_traits_c, recursive = F)

folders_traits_y_a <- "osmose-ns-v4.4.1_clim/output/PLOTS/Trait_variables_by_year_by_age2/RData/"
folders_traits_y_a <- list.dirs(folders_traits_y_a, recursive = F)

folders_traits_c_a <- "osmose-ns-v4.4.1_clim/output/PLOTS/Trait_variables_by_cohort_by_age2/RData/"
folders_traits_c_a  <- list.dirs(folders_traits_c_a, recursive = F)

 # pathout <- paste0(configDir, "/output/PLOTS/Htz_vs_Ne_variables")


# load Ne dataset ------------------------------------------------------
# load functions
source("osmose-ns-v4.4.1_clim/R/functions/auxiliary_for_lists.R")
source("osmose-ns-v4.4.1_clim/R/functions/functions_output_SimuRData_into_df_for_plots.R")





pathoutputRDATA <-paste0(configDir, "/outpout/RData")
pathout <- paste0(configDir, "/output/PLOTS/Ne_vs_F")
if (!dir.exists(pathout)) dir.create(pathout, recursive = T)
# from RData output of the loops producing the plots of Trait by time and Genetic drift indicators, to df
# produce a dataframe with the right format : Ne variables and htz variables together
Time1=150
Time2=160
species_df=species
# list_data_traits <- data_trait_and_genetic_into_df(Time1=150, Time2=160, pathoutput=paste0(configDir, "/output/PLOTS/RDATA"), 
#                                                    FORCE =F, species = species_df)

# produce a dataframe with the right format : Ne variables and htz variables together
list_data_traits_ratio <- data_trait_and_genetic_into_df(Time1=150, Time2=160, pathoutput=paste0(configDir, "/output/PLOTS/RDATA"), 
                                                         FORCE = F, species = species_df, RATIO=T)
# list of outputs 
# list_outputs <- list("traits"=list_data_traits, "trait_ratio"=list_data_traits_ratio)

# # Bind together the trait data coming from similar dimensions : cohort or year data together (by genetic drift dimension)

pathoutputRDATA = paste0(configDir, "/output/PLOTS/RDATA")
# data_full_merged <- dataset_attach_all_traits_genetic_drift(pathoutputRDATA, list_outputs=list("traits"=list_data_traits), 
#                                                             trait_cv = list("folders_traits_c"=folders_traits_c, "folders_traits_y"=folders_traits_y), FORCE=F)
dataRATIO_full_merged <- dataset_attach_all_traits_genetic_drift(pathoutputRDATA, list_outputs=list("ratio"=list_data_traits_ratio),
                                                                 trait_cv = list("folders_traits_c"=folders_traits_c, "folders_traits_y"=folders_traits_y), FORCE=F)
# 
# 
# data_full_merged_age <- dataset_attach_all_traits_genetic_drift(pathoutputRDATA, list_outputs = list("traits"=list_data_traits), s_pecies = species,
#                                                                 trait_cv = list("folders_traits_c_a"=folders_traits_c_a, "folders_traits_y_a"=folders_traits_c_a),
#                                                                 keep_cohort_or_age_dimension = T, FORCE=F)
# 
# # 
# # # Dataset ready for modeling --------------------------------------------
# # # scaling as a standardized loop
# data_mean_sc = scale_dataset_list(pathoutputRDATA = pathoutputRDATA, data_mean=data_full_merged, FORCE=T)
dataRATIO_mean_sc = scale_dataset_list(pathoutputRDATA, data_mean=dataRATIO_full_merged, FORCE=F)
# 
# # new dataset of the difference between Fishing pressure 1,1 and Else--------------------
# diff_data = relative_dataset(data_mean_sc, pathoutputRDATA, FORCE=T)
diffRATIO_data = relative_dataset(dataRATIO_mean_sc, pathoutputRDATA, FORCE=F, RATIO=T)








# Plot the biomass dynamics : to insert in other plots -------------------
biomassplot = NULL
rel_biomassplot = NULL
for (nm in names(dataRATIO_full_merged)){
  d <- dataRATIO_full_merged[[nm]]
  d <- d[d$Fstart==1, ]
  id_biom <- grep(pattern='biomass', colnames(d))
  name_id <- colnames(d)[id_biom]
  for (sp in species){
    d_i_i <- d[which(d$Species == sp), ]
    if (any(is.na(d_i_i[, name_id]))) d_i_i <- d_i_i[-which(is.na(d_i_i[, name_id])), ]
    d_i_i$Iszero <- ifelse(d_i_i[,name_id] == 0, "Absolute value is 0", "Different from 0")
    
    
    biomassplot[[sp]] <- ggplot(d_i_i, aes(y=.data[[name_id]], x=Fmax)) +
      geom_point(aes(colour = Iszero))+
      labs(title = "Biomass per fishing pressure", subtitle = sp)
    
    ref_data <- d_i_i %>%
      filter(Fmax == 1 & Fstart == 1) %>%
      dplyr::select(.data[[name_id]], replicate) %>%
      dplyr::rename(ref_value = .data[[name_id]])
    # Subtract reference from all means
    rel_d_i_i <- d_i_i %>%
      left_join(ref_data, by = c("replicate"), relationship = "many-to-one") %>%
      mutate(relative_biomass =(.data[[name_id]] - ref_value)/ref_value) %>%
      dplyr::select(-ref_value)
    rel_biomassplot[[sp]] <- ggplot(rel_d_i_i, aes(y=relative_biomass, x=Fmax)) +
      geom_point(aes(colour = Iszero))+
      labs(title = "Relatve biomass per fishing pressure", subtitle = sp)
  }
}

  














# PLOT NE BY FISHING PRESSURE ---------------------

#summary <- get(load(file=paste0(pathoutputRDATA, "/summary.RData")))
#summary$Fmax <- paste0("Fmax=", as.numeric(summary$Fmax))
#summary$Fstart <- paste0("Fstart=", as.numeric(summary$Fstart))

# list_data_traits_ratio

# Period range : 
# Time1=150
# Time2=160

# to plot only Ne and not Ne/N, change ↓ list_data_traits_ratio to list_data_traits

genet_drift_data_full       <- list_data_traits_ratio[["genet_drift_data"]]
trait_data_full             <- list_data_traits_ratio[["trait_data"]]

relative_to_ref = FALSE  # Set to TRUE to plot relative to Fmax==1 and Fstart==1
for (relative_to_ref in c(FALSE)){
scale_v = list()
  for (i in (1:length(genet_drift_data_full))[-3]){
    summary = genet_drift_data_full[[i]]
    names_substr <- substr(summary$variable, nchar(summary$variable) - 3 + 1, nchar(summary$variable))
    for (TRAIT in unique(names_substr)){
        df_trait <- summary[grep(summary$variable, pattern = TRAIT),]
        var_i <- unique(df_trait$variable)
        
        p = list()
        for (sp in unique(df_trait$Species)){
          df_trait_sp <- df_trait %>% filter(Species == sp)
          df_trait_sp <- df_trait_sp %>% filter(Fstart == 1)
          if (sum(is.nan(df_trait_sp$value))>0) df_trait_sp <- df_trait_sp[-which(is.nan(df_trait_sp$value)), ]
          if (sum(is.infinite(df_trait_sp$value))>0) df_trait_sp <- df_trait_sp[-which(is.infinite(df_trait_sp$value)), ]
          
          # Use the function
          result <- plot_trait_vs_fmax(df_trait_sp = df_trait_sp, pathout = pathout, trait_name = var_i, 
                                       Time1, Time2, relative_to_ref = relative_to_ref, spec=sp)
          p[[sp]] <- result$plot
          scale_v[[sp]] <- result$scale
        }#sp
        trait_i = unique(var_i)
        if (length(unique(var_i))) trait_i = paste0(unique(var_i)[1], "_WITH_")
        suffix_relative <- ifelse(relative_to_ref, "_RELATIVE", "")
        pdf(file = paste0(pathout, "/plot_F1_", trait_i, "vsFmax", suffix_relative, ".pdf"), width = 10, height = 8)
        par(mfrow=c(1,1))
        # ggarrange(
        for (sp in 1:length(scale_v)){
          scale=scale_v[[sp]] 
          print(p[[sp]])
        }#, ncol=6, nrow=3)
        dev.off()
    }
  }
}





# PLOT ALL genet_drift_data_full VARIABLES ON SAME PAGE PER SPECIES ---------------------

genet_drift_data_full       <- list_data_traits_ratio[["genet_drift_data"]]
genet_drift_data_full       <- c(genet_drift_data_full, "htz_slope" = list(list_data_traits_ratio[["trait_data"]]$htz_list))

for (relative_to_ref_i in c(TRUE, FALSE)){
  p_list <- list()
  for (sp in unique(genet_drift_data_full[[1]]$Species)) {
    # p_list <- list()
    
    for (i in (1:length(genet_drift_data_full))) {
      summary <- genet_drift_data_full[[i]]
      df_trait_sp <- summary %>% filter(Species == sp)
      df_trait_sp <- df_trait_sp %>% filter(Fstart == 1)
      
      if (nrow(df_trait_sp) == 0) next
      
      if (sum(is.nan(df_trait_sp$value)) > 0) df_trait_sp <- df_trait_sp[-which(is.nan(df_trait_sp$value)), ]
      if (sum(is.infinite(df_trait_sp$value)) > 0) df_trait_sp <- df_trait_sp[-which(is.infinite(df_trait_sp$value)), ]
      
      if (nrow(df_trait_sp) == 0) next
      
      var_i <- unique(df_trait_sp$variable)
      
      if (length(var_i)>2) {
        for (i in 1:c(length(var_i)/2)) {
          var_i_subset <- var_i[((i-1)*2+1):(i*2)]
          result <- plot_trait_vs_fmax(df_trait_sp = df_trait_sp[which(df_trait_sp$variable %in% var_i_subset),], 
                                       pathout =pathout, trait_name = var_i_subset, 
                                     Time1, Time2, relative_to_ref = relative_to_ref_i, spec = sp)
          p_list[[sp]][[var_i_subset[1]]] <- result$plot
        }
      } else {
        result <- plot_trait_vs_fmax(df_trait_sp = df_trait_sp, pathout = pathout, trait_name = var_i, 
                                     Time1, Time2, relative_to_ref = relative_to_ref_i, spec = sp)
        p_list[[sp]][[names(genet_drift_data_full)[i]]] <- result$plot
      }
    }
  }
  if (length(p_list) > 0) {
    plot_biom = ifelse(relative_to_ref_i, "rel_biomassplot", "biomassplot")
    suffix_relative <- paste0("_Combined", ifelse(relative_to_ref_i, "relative_", ""))
    pdf(file = paste0(pathout, "/plot_all_genet_drift_variables_", suffix_relative, ".pdf"), width = 14, height = 12)
    for (sp in species){
      combined_plot <- ggarrange(plotlist = c(p_list[[sp]], list(get(plot_biom)[[sp]])), ncol = 2, nrow = 3)
      print(combined_plot)
    }
    dev.off()
  }
}



# PLOT genet_drift_data_full VARIABLES WITH T_Vk ---------------------

genet_drift_data_full       <- list_data_traits_ratio[["genet_drift_data"]]
trait_data_full             <- list_data_traits_ratio[["trait_data"]]
traits_list_to_plot <- list(c("generation_length", "V"))

for (relative_to_ref in c(FALSE)){
  scale_v = list()
  for (i in (1:length(genet_drift_data_full))){
    summary = genet_drift_data_full[[i]]
    names_substr <- substr(summary$variable, nchar(summary$variable) - 3 + 1, nchar(summary$variable))
    for (TRAIT in unique(names_substr)){
      df_trait <- summary[grep(summary$variable, pattern = TRAIT),]
      var_i <- unique(df_trait$variable)
      
      # Calculate T_Vk
      vec_i <- c(grep(trait_data_full[[names(genet_drift_data_full)[i]]]$variable, pattern=traits_list_to_plot[[1]][1]), 
                 grep(trait_data_full[[names(genet_drift_data_full)[i]]]$variable, pattern=traits_list_to_plot[[1]][2]))
      df_trait_t <- trait_data_full[[names(genet_drift_data_full)[i]]][vec_i,]
      df_trait_t <- pivot_wider(df_trait_t, names_from = "variable", values_from = "value") 
      if (ncol(df_trait_t) != 6) next
      df_trait_t[, "T_Vk"] <- df_trait_t[, grep(pattern = traits_list_to_plot[[1]][1], colnames(df_trait_t))]/
        df_trait_t[, grep(pattern = traits_list_to_plot[[1]][2], colnames(df_trait_t))]
      df_trait_t <- pivot_longer(df_trait_t, cols = colnames(df_trait_t)[-c(1:4)], values_to = "value", names_to = "variable")
      
      df_trait <- rbind(df_trait, df_trait_t[which(df_trait_t$variable == "T_Vk"),])
      
      p = list()
      for (sp in unique(df_trait$Species)){
        df_trait_sp <- df_trait %>% filter(Species == sp)
        df_trait_sp <- df_trait_sp %>% filter(Fstart == 1)
        if (sum(is.nan(df_trait_sp$value))>0) df_trait_sp <- df_trait_sp[-which(is.nan(df_trait_sp$value)), ]
        if (sum(is.infinite(df_trait_sp$value))>0) df_trait_sp <- df_trait_sp[-which(is.infinite(df_trait_sp$value)), ]
        
        if (nrow(df_trait_sp)==0) next
    
        # Use the function
        for (tri in unique(df_trait$variable)){
          result <- plot_trait_vs_fmax(df_trait_sp = df_trait_sp[which(df_trait_sp$variable %in% c(tri, "T_Vk")), ], pathout = pathout, trait_name = c(tri, "T_Vk"), 
                                       Time1, Time2, relative_to_ref = relative_to_ref, spec=sp)
          p[[sp]][[tri]] <- result$plot
          scale_v[[sp]] <- result$scale
        }
      }#sp
      trait_i = unique(var_i)
      if (length(unique(var_i))) trait_i = paste0(unique(var_i)[1], "_WITH_")
      suffix_relative <- ifelse(relative_to_ref, "_RELATIVE", "")
      pdf(file = paste0(pathout, "/plot_", names(genet_drift_data_full)[i], "_T_Vk_", trait_i, "vsFmax", suffix_relative, ".pdf"), width = 10, height = 8)
      par(mfrow=c(1,1))
      for (sp in species){
        scale=scale_v[[sp]]
        combined_plot <- ggarrange(plotlist = p[[sp]], ncol = 2, nrow = ceiling(length(p[[sp]])/2))
        print(combined_plot)
      }
      dev.off()
    }
  }
}


# PLOT list(c("generation_length", "V")) BY FISHING PRESSURE ---------------------

#summary <- get(load(file=paste0(pathoutputRDATA, "/summary.RData")))
#summary$Fmax <- paste0("Fmax=", as.numeric(summary$Fmax))
#summary$Fstart <- paste0("Fstart=", as.numeric(summary$Fstart))

# list_data_traits_ratio

# Period range : 
# Time1=150
# Time2=160

genet_drift_data_full       <- list_data_traits_ratio[["genet_drift_data"]]
trait_data_full             <- list_data_traits_ratio[["trait_data"]]
trait_data_full$data_Nb     <- rbind(trait_data_full$data_Nb, trait_data_full$tempoNe)
# genet_drift_data_full[["htz_slope"]] <- trait_data_full$htz_slope

traits_list_to_plot <- list(c("generation_length", "V"))
traits_list_to_plot <- list(c("generation_length", "V"))

for (i in 1:length(trait_data_full)){
  
  summary = trait_data_full[[i]]
  for (TRAIT in 1:length(traits_list_to_plot)){
    
    df_trait=NULL
    vec_i = c()
    i_i = c()
    for (j in 1:length(traits_list_to_plot[[TRAIT]])){
      vec_i <- c(vec_i, grep(pattern = traits_list_to_plot[[TRAIT]][j], summary$variable))
      i_i <- c(i_i, any(grepl(pattern = traits_list_to_plot[[TRAIT]][j], summary$variable)))
    }
    if (length(vec_i)==0) next
    
    df_trait <- summary[vec_i,]
    if (!(length(traits_list_to_plot[[TRAIT]]) == length(unique(df_trait$variable)))) next
    df_trait_t <- pivot_wider(df_trait, names_from = "variable", values_from = "value") 
    df_trait_t[, "T_Vk"] <- df_trait_t[, grep(pattern = traits_list_to_plot[[TRAIT]][1], colnames(df_trait_t))]/
                               df_trait_t[, grep(pattern = traits_list_to_plot[[TRAIT]][2], colnames(df_trait_t))]
    df_trait <- pivot_longer(df_trait_t, cols = colnames(df_trait_t)[-c(1:4)], values_to = "value", names_to = "variable")
    
    p = list()
    q = list()
    scale_v = c()
    for (relative_to_ref in c(TRUE, FALSE)){
      for (sp in unique(df_trait$Species)){
        df_trait_sp <- df_trait %>% filter(Species == sp)
        df_trait_sp <- df_trait_sp %>% filter(Fstart == 1)
        df_trait_sp <- df_trait_sp %>%
          group_by(Fmax, Fstart, variable) %>%
          summarise(mean = mean(value, na.rm = TRUE),
                    sd   = sd(value, na.rm = TRUE))   
        scale = mean(df_trait_sp$mean[grep(pattern =  "T_Vk", df_trait_sp$variable)], na.rm=T)/
                mean(df_trait_sp$mean[grep(pattern = traits_list_to_plot[[TRAIT]][2], df_trait_sp$variable)], na.rm=T)
        scale_v[sp] = scale
        df_trait_sp$mean[grep(pattern =  "T_Vk", df_trait_sp$variable)] <- df_trait_sp$mean[grep(pattern =  "T_Vk", df_trait_sp$variable)]/scale
        df_trait_sp$sd[grep(pattern =  "T_Vk", df_trait_sp$variable)] <- df_trait_sp$sd[grep(pattern =  "T_Vk", df_trait_sp$variable)]/scale
        
        # Use the function
        p[[sp]] <- plot_trait_vs_fmax_simple(df_trait_sp = df_trait_sp, pathout =  pathout, sp_name = sp, 
                                             var_names = unique(df_trait_sp$variable), Time1 = Time1, Time2 = Time2, 
                                              relative_to_ref = relative_to_ref, scale_factor = scale)
      }
      suffix_relative <- ifelse(relative_to_ref, "_RELATIVE", "")
      pdf(file = paste0(pathout, "/plot_", names(genet_drift_data_full)[i], "_", paste(unique(df_trait_sp$variable), collapse ="_"), "_vsFmax", suffix_relative, ".pdf"), width = 10, height = 8)
      par(mfrow=c(1,1))
      for (sp in names(scale_v)){
        scale=scale_v[[sp]]
        print(p[[sp]])
      }
      dev.off()
    }
  }
}



# PLOT list(c("mortality"), c("year_mortality"), "cohort_sum_of_mortalityByAge_per_age") BY FISHING PRESSURE ---------------------

#summary <- get(load(file=paste0(pathoutputRDATA, "/summary.RData")))
#summary$Fmax <- paste0("Fmax=", as.numeric(summary$Fmax))
#summary$Fstart <- paste0("Fstart=", as.numeric(summary$Fstart))

# list_data_traits_ratio

# Period range : 
# Time1=150
# Time2=160

data_full_merged

# genet_drift_data_full[["htz_slope"]] <- trait_data_full$htz_slope

traits_list_to_plot <- list(c("mortality"), c("year_mortality"), c("cohort_sum_of_mortalityByAge_per_age")) 
traits_list_to_plot <- list(c("mortality"), c("year_mortality"), c("cohort_sum_of_mortalityByAge_per_age"))
# traits_list_to_plot <- list(c("generation_length", "V"))

for (i in 1:length(dataRATIO_full_merged)){
  
  summary = dataRATIO_full_merged[[i]]
  for (TRAIT in 1:length(traits_list_to_plot)){
    
    if (ncol(summary)>6) summary <- pivot_longer(data=summary, values_to = "value", cols = colnames(summary)[-c(1:4)], names_to = "variable")
    df_trait=NULL
    p = NULL
    vec_i = c()
    i_i = c()
    for (j in 1:length(traits_list_to_plot[[TRAIT]])){
      vec_i <- c(vec_i, which(summary$variable %in% traits_list_to_plot[[TRAIT]][j]))
      i_i <- c(i_i, any(grepl(pattern = traits_list_to_plot[[TRAIT]][j], summary$variable)))
    }
    if (length(vec_i)==0) next
    
    df_trait <- summary[vec_i,]
    # if (!(length(traits_list_to_plot[[TRAIT]]) == length(unique(df_trait$variable)))) next
    
    p = list()
    q = list()
    for (relative_to_ref in c(TRUE, FALSE)){
      for (sp in unique(df_trait$Species)){
        df_trait_sp <- df_trait %>% filter(Species == sp)
        df_trait_sp <- df_trait_sp %>% filter(Fstart == 1)
        df_trait_sp <- df_trait_sp %>%
          group_by(Fmax, Fstart, variable) %>%
          summarise(mean = mean(value, na.rm = TRUE),
                    sd   = sd(value, na.rm = TRUE))   
        
        # Use the function
        p[[sp]] <- plot_trait_vs_fmax_simple(df_trait_sp, pathout, sp, unique(df_trait_sp$variable), Time1, Time2, 
                                              relative_to_ref = relative_to_ref, scale_factor = 1)
      }
      suffix_relative <- ifelse(relative_to_ref, "_RELATIVE", "")
      pdf(file = paste0(pathout, "/plot_", names(dataRATIO_full_merged)[i], (traits_list_to_plot)[[TRAIT]][1], "_", paste(unique(df_trait_sp$variable), collapse ="_"), "_vsFmax", suffix_relative, ".pdf"), width = 10, height = 8)
      par(mfrow=c(1,1))
      for (sp in names(p)){
        scale=scale_v[[sp]]
        print(p[[sp]])
      }
      dev.off()
    }
  }
}
  





# PLOT T & Ne BY FISHING PRESSURE ---------------------



genet_drift_data_full       <- list_data_traits_ratio[["genet_drift_data"]]
trait_data_full             <- list_data_traits_ratio[["trait_data"]]
trait_data_full$data_Nb     <- rbind(trait_data_full$data_Nb, trait_data_full$tempoNe)
traits_list_to_plot <- list(c("generation_length", "V"))
relative_to_ref = FALSE  # Set to TRUE to plot relative to Fmax==1 and Fstart==1

# with htz
for (relative_to_ref in c(FALSE)){
  scale_v = list()
  for (i in (1:length(genet_drift_data_full))){
    summary = rbind(genet_drift_data_full[[i]], trait_data_full[[names(genet_drift_data_full)[i]]])
    names_substr <- substr(genet_drift_data_full[[i]]$variable, nchar(summary$variable) - 3 + 1, nchar(summary$variable))
    for (TRAIT in unique(names_substr)){
      df_trait <- summary[grep(summary$variable, pattern = TRAIT),]
      var_i <- unique(df_trait$variable)
      df_trait <- df_trait[grep(df_trait$variable, pattern="ratio"),]
      df_trait <- rbind(trait_data_full$htz_list, df_trait)
      
      vec_i <- c(grep(summary$variable, pattern=traits_list_to_plot[[1]][1]), grep(summary$variable, pattern=traits_list_to_plot[[1]][2]))
      df_trait_t <- summary[vec_i,]
      if (!(length(traits_list_to_plot[[1]]) == length(unique(df_trait_t$variable)))) next
      df_trait_t <- pivot_wider(df_trait_t, names_from = "variable", values_from = "value") 
      df_trait_t[, "T_Vk"] <- df_trait_t[, grep(pattern = traits_list_to_plot[[1]][1], colnames(df_trait_t))]/
        df_trait_t[, grep(pattern = traits_list_to_plot[[1]][2], colnames(df_trait_t))]
      df_trait_t <- pivot_longer(df_trait_t, cols = colnames(df_trait_t)[-c(1:4)], values_to = "value", names_to = "variable")
      
      df_trait <- rbind(df_trait, df_trait_t[which(df_trait_t$variable == "T_Vk"),])
      
      p = list()
      for (sp in unique(df_trait$Species)){
        df_trait_sp <- df_trait %>% filter(Species == sp)
        df_trait_sp <- df_trait_sp %>% filter(Fstart == 1)
        if (sum(is.nan(df_trait_sp$value))>0) df_trait_sp <- df_trait_sp[-which(is.nan(df_trait_sp$value)), ]
        if (sum(is.infinite(df_trait_sp$value))>0) df_trait_sp <- df_trait_sp[-which(is.infinite(df_trait_sp$value)), ]

        # Use the function
        result <- plot_trait_vs_fmax(df_trait_sp = df_trait_sp, pathout = pathout, trait_name = var_i, 
                                     Time1, Time2, relative_to_ref = relative_to_ref, spec=sp)
        p[[sp]] <- result$plot
        scale_v[[sp]] <- result$scale
      }#sp
      trait_i = unique(var_i)
      if (length(unique(var_i))) trait_i = paste0(unique(var_i)[1], "_WITH_")
      suffix_relative <- ifelse(relative_to_ref, "_RELATIVE", "")
      pdf(file = paste0(pathout, "/plot_T_Vk_Ne_loess_", trait_i, "vsFmax", suffix_relative, ".pdf"), width = 10, height = 8)
      par(mfrow=c(1,1))
      for (sp in 1:length(scale_v)){
        scale=scale_v[[sp]]
        print(p[[sp]])
      }
      dev.off()
    }
  }
}
#without htz
for (relative_to_ref in c(FALSE)){
  scale_v = list()
  for (i in (1:length(genet_drift_data_full))){
    summary = rbind(genet_drift_data_full[[i]], trait_data_full[[names(genet_drift_data_full)[i]]])
    names_substr <- substr(genet_drift_data_full[[i]]$variable, nchar(summary$variable) - 3 + 1, nchar(summary$variable))
    for (TRAIT in unique(names_substr)){
      df_trait <- summary[grep(summary$variable, pattern = TRAIT),]
      var_i <- unique(df_trait$variable)
      df_trait <- df_trait[grep(df_trait$variable, pattern="ratio"),]
      
      vec_i <- c(grep(summary$variable, pattern=traits_list_to_plot[[1]][1]), grep(summary$variable, pattern=traits_list_to_plot[[1]][2]))
      df_trait_t <- summary[vec_i,]
      if (!(length(traits_list_to_plot[[1]]) == length(unique(df_trait_t$variable)))) next
      df_trait_t <- pivot_wider(df_trait_t, names_from = "variable", values_from = "value") 
      df_trait_t[, "T_Vk"] <- df_trait_t[, grep(pattern = traits_list_to_plot[[1]][1], colnames(df_trait_t))]/
        df_trait_t[, grep(pattern = traits_list_to_plot[[1]][2], colnames(df_trait_t))]
      df_trait_t <- pivot_longer(df_trait_t, cols = colnames(df_trait_t)[-c(1:4)], values_to = "value", names_to = "variable")
      
      df_trait <- rbind(df_trait, df_trait_t[which(df_trait_t$variable == "T_Vk"),])
      
      p = list()
      for (sp in unique(df_trait$Species)){
        df_trait_sp <- df_trait %>% filter(Species == sp)
        df_trait_sp <- df_trait_sp %>% filter(Fstart == 1)
        if (sum(is.nan(df_trait_sp$value))>0) df_trait_sp <- df_trait_sp[-which(is.nan(df_trait_sp$value)), ]
        if (sum(is.infinite(df_trait_sp$value))>0) df_trait_sp <- df_trait_sp[-which(is.infinite(df_trait_sp$value)), ]
        
        # Use the function
        result <- plot_trait_vs_fmax(df_trait_sp = df_trait_sp, pathout = pathout, trait_name = var_i, 
                                     Time1, Time2, relative_to_ref = relative_to_ref, spec=sp)
        p[[sp]] <- result$plot
        scale_v[[sp]] <- result$scale
      }#sp
      trait_i = unique(var_i)
      if (length(unique(var_i))) trait_i = paste0(unique(var_i)[1], "_WITH_")
      suffix_relative <- ifelse(relative_to_ref, "_RELATIVE", "")
      pdf(file = paste0(pathout, "/plot_T_Vk_Ne_loess_", trait_i, "vsFmax", suffix_relative, ".pdf"), width = 10, height = 8)
      par(mfrow=c(1,1))
      for (sp in 1:length(scale_v)){
        scale=scale_v[[sp]]
        print(p[[sp]])
      }
      dev.off()
    }
  }
}




# PLOT TRAITS without T_Vk


a_full_merged

# genet_drift_data_full[["htz_slope"]] <- trait_data_full$htz_slope

traits_list_to_plot <- list(c("mortality"), c("year_mortality"), c("cohort_sum_of_mortalityByAge_per_age")) 
traits_list_to_plot <- list(c("mortality"), c("year_mortality"), c("cohort_sum_of_mortalityByAge_per_age"))
traits_list_to_plot <- list(c("mortality"), c("year_mortality", "year_abundanceByAge"), c("cohort_sum_of_mortalityByAge_per_age", "cohort_sum_of_abundanceByAge_per_y"))
traits_list_to_plot <- list(c("mortality"), c("year_eggNumberByAge", "year_abundanceByAge"), c("cohort_sum_of_eggNumberByAge_per_y", "cohort_sum_of_abundanceByAge_per_y"))
traits_list_to_plot <- list(c("mortality"), c("year_eggNumberByAge", "year_abundanceByAge"), c("cohort_sum_of_eggNumberByAge_per_y", "cohort_sum_of_abundanceByAge_per_y"))
traits_list_to_plot <- list(c("mortality"), c("year_MnatbyY", "year_abundanceByAge"), c("cohort_mean_of_Mnat_per_y", "cohort_sum_of_abundanceByAge_per_y"))
traits_list_to_plot <- list(c("mortality"), c("year_MpredbyY", "year_abundanceByAge"), c("cohort_mean_of_Mpred_per_y", "cohort_sum_of_abundanceByAge_per_y"))
traits_list_to_plot <- list(c("mortality"), c("year_MstarbyY", "year_abundanceByAge"), c("cohort_mean_of_Mstar_per_y", "cohort_sum_of_abundanceByAge_per_y"))



# traits_list_to_plot <- list(c("generation_length", "V"))
scale_v = c()
for (i in 1:length(dataRATIO_full_merged)){
  
  summary = dataRATIO_full_merged[[i]]
  for (TRAIT in 1:length(traits_list_to_plot)){
    
    if (ncol(summary)>6) summary <- pivot_longer(data=summary, values_to = "value", cols = colnames(summary)[-c(1:4)], names_to = "variable")
    df_trait=NULL
    p = NULL
    vec_i = c()
    i_i = c()
    for (j in 1:length(traits_list_to_plot[[TRAIT]])){
      vec_i <- c(vec_i, which(summary$variable %in% traits_list_to_plot[[TRAIT]][j]))
      i_i <- c(i_i, any(grepl(pattern = traits_list_to_plot[[TRAIT]][j], summary$variable)))
    }
    if (length(vec_i)==0) next
    
    df_trait <- summary[vec_i,]
    # if (!(length(traits_list_to_plot[[TRAIT]]) == length(unique(df_trait$variable)))) next
    
    p = list()
    q = list()
    for (relative_to_ref in c(TRUE, FALSE)){
      for (sp in unique(df_trait$Species)){
        df_trait_sp <- df_trait %>% filter(Species == sp)
        df_trait_sp <- df_trait_sp %>% filter(Fstart == 1)
        # df_trait_sp <- df_trait_sp %>%
        #   group_by(Fmax, Fstart, variable) %>%
        #   summarise(mean = mean(value, na.rm = TRUE),
        #             sd   = sd(value, na.rm = TRUE))   
        
        # Use the function
        result <- plot_trait_vs_fmax(df_trait_sp = df_trait_sp, pathout = pathout, spec = sp, trait_name = unique(df_trait_sp$variable), 
                                             Time1, Time2, 
                                             relative_to_ref = relative_to_ref)
        p[[sp]] = result[[1]]
        scale_v[[sp]] = result[[2]]
      }
      suffix_relative <- ifelse(relative_to_ref, "_RELATIVE", "")
      pdf(file = paste0(pathout, "/plot_", names(dataRATIO_full_merged)[i], (traits_list_to_plot)[[TRAIT]][1], "_", paste(unique(df_trait_sp$variable), collapse ="_"), "_vsFmax", suffix_relative, ".pdf"), width = 10, height = 8)
      par(mfrow=c(1,1))
      for (spi in unique(df_trait$Species)){
        scale=scale_v[[spi]]
        print(p[[spi]])
      }
      dev.off()
    }
  }
}



