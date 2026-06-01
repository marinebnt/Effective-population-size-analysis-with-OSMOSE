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


# functions
main = "New_gpu_conf.R" # config file
source("osmose-ns-v4.4.1_clim/R/functions/source.R") # load config and key parameters and functions
source("osmose-ns-v4.4.1_clim/R/functions/functions_key_traits.R") # plots against Fmax
source("osmose-ns-v4.4.1_clim/R/functions/auxiliary_for_lists.R") # from list to df 
source("osmose-ns-v4.4.1_clim/R/functions/functions_output_SimuRData_into_df_for_plots.R") # to transform dataset from RData to long df


# paths
PATHOUT  = "osmose-ns-v4.4.1_clim/output/PLOTS/KEY_TRAITS/"
if (!dir.exists(PATHOUT)) dir.create(PATHOUT, recursive = TRUE)
base_dir <- "osmose-ns-v4.4.1_clim/output/PLOTS/FINAL_Ne_Htz/INIT"
folders  <- list.dirs(base_dir, recursive = F) # files INIT2 with raw Ne datasets

folders_traits_y <- "osmose-ns-v4.4.1_clim/output/PLOTS/Trait_variables_by_year_BYSP//RData/"
folders_traits_y <- list.dirs(folders_traits_y, recursive = F)
folders_traits_c <- "osmose-ns-v4.4.1_clim/output/PLOTS/Trait_variables_by_cohort/RData/"
folders_traits_c  <- list.dirs(folders_traits_c, recursive = F)
folders_traits_y_a <- "osmose-ns-v4.4.1_clim/output/PLOTS/Trait_variables_by_year_by_age2/RData/"
folders_traits_y_a <- list.dirs(folders_traits_y_a, recursive = F)
folders_traits_c_a <- "osmose-ns-v4.4.1_clim/output/PLOTS/Trait_variables_by_cohort_by_age2/RData/"
folders_traits_c_a  <- list.dirs(folders_traits_c_a, recursive = F)

# load or format dataset
pathoutputRDATA = paste0(configDir, "/output/PLOTS/KEY_TRAITS/RDATA")
if(!dir.exists(pathoutputRDATA)) dir.create(pathoutputRDATA, recursive = T)
Time1=150
Time2=160
species_df=species
    # df by loading the trait datasets by Fish_press (by year, by cohort, by year by age, by cohort by age) and binding them together
list_data_traits_ratio <- data_trait_and_genetic_into_df(folders_in = folders, Time1=150, Time2=160, pathoutput=pathoutputRDATA, 
                                                         FORCE = F, species = species_df, RATIO=T)
    # df of only the variables used to estimate Ne, by age (mean over cohorts 150-160 is weighted by abundance)
only_Ne_traits <- extract_INIT2_data_Ne_weighted_by_age(pathoutputRDATA, folders_in=folders, Time1=150, Time2=160, species = species, FORCE=F)
#     # df by binding all genetic drift indic + trait datasets (by year, by cohort, by year by age, by cohort by age)
dataRATIO_full_merged <- dataset_attach_all_traits_genetic_drift(pathoutputRDATA, list_outputs=list("ratio"=list_data_traits_ratio),
                                                                  trait_cv = list("folders_traits_c"=folders_traits_c, "folders_traits_y"=folders_traits_y), 
                                                                  FORCE=F)
#     # df by binding all genetic drift indic + trait datasets (by year, by cohort, by year by age, by cohort by age) WITH AGE dimension                                                                
dataRATIO_full_merged_age <- dataset_attach_all_traits_genetic_drift(pathoutputRDATA, list_outputs = list("ratio"=list_data_traits_ratio), s_pecies = species,
                                                                trait_cv = list("folders_traits_c_a"=folders_traits_c_a, "folders_traits_y_a"=folders_traits_y_a),
                                                                keep_cohort_or_age_dimension = T, FORCE=F)
#     # scale the dataset so trait values can be plotted together with a common secondary axis, and compute the relative change vs reference for the ratio dataset.
dataRATIO_mean_sc = scale_dataset_list(pathoutputRDATA, data_mean=dataRATIO_full_merged, FORCE=F)
    # extract only the relative change vs reference for the ratio dataset, to be plotted with traits.
diffRATIO_data = relative_dataset(dataRATIO_mean_sc, pathoutputRDATA, FORCE=F, RATIO=T)
# 
# 


# BIPLOTS -------------------------------------------------------

## PLOT 0- Biomass vs Fmax & keep it on the side -------------------
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
list_bio <- list("rel_biomassplot"=rel_biomassplot, "biomassplot"=biomassplot)
save(list_bio, file=paste0(pathoutputRDATA, "/Biomass_plots.RData"))
# load(file=paste0(paste0(configDir, "/output/PLOTS/RDATA"), "/Biomass_plots.RData"))
rel_biomassplot=list_bio$rel_biomassplot
biomassplot=list_bio$biomassplot

## PLOT 1- NE/N VS FMAX + Biomass VS FMAX ------------------

genet_drift_data_full       <- list_data_traits_ratio[["genet_drift_data"]]
genet_drift_data_full       <- c(genet_drift_data_full, "htz_slope" = list(list_data_traits_ratio[["trait_data"]]$htz_list))
dataRATIO_full_merged_age
for (relative_to_ref_i in c(TRUE, FALSE)[1]){
  p_list <- list()
  for (sp in unique(dataRATIO_full_merged_age[[1]]$Species)[1]) {

    for (i in (1:length(dataRATIO_full_merged_age$data_full_merged))[1]) {
      summary <- dataRATIO_full_merged_age$data_full_merged[[i]]
      df_trait_sp <- summary %>% filter(Species == sp)
      df_trait_sp <- df_trait_sp %>% filter(Fstart == 1)

      if (nrow(df_trait_sp) == 0) next

      if (sum(is.nan(df_trait_sp$value)) > 0) df_trait_sp <- df_trait_sp[-which(is.nan(df_trait_sp$value)), ]
      if (sum(is.infinite(df_trait_sp$value)) > 0) df_trait_sp <- df_trait_sp[-which(is.infinite(df_trait_sp$value)), ]

      if (nrow(df_trait_sp) == 0) next

      var_i <- c("Ne_N_1", "Ne_N_1/Numerator_N_1", "mean_Vk_point")

      if (length(var_i)>2) {
        for (i in 1:c(length(var_i)/2)) {
          var_i_subset <- var_i[((i-1)*2+1):(i*2)]
          result <- plot_trait_vs_fmax(df_trait_sp = df_trait_sp[which(df_trait_sp$variable %in% var_i_subset),],
                                       trait_name = var_i_subset,
                                     Time1, Time2, relative_to_ref = relative_to_ref_i, spec = sp, age_dim = F, 
                                     proportional = F, sp_dim = F, weight=build_weight_df(df_trait_sp, key_cols = c("replicate", "Fstart", "Fmax", "Species")))
          p_list[[sp]][[var_i_subset[1]]] <- result$plot
        }
      } else {
        result <- plot_trait_vs_fmax(df_trait_sp = df_trait_sp,
                                     trait_name = var_i,
                                     Time1, Time2, relative_to_ref = relative_to_ref_i, spec = sp)
        p_list[[sp]][[names(genet_drift_data_full)[i]]] <- result$plot
      }
    }
  }
  if (length(p_list) > 0) {
    plot_biom = ifelse(relative_to_ref_i, "rel_biomassplot", "biomassplot")
    plot_biom = ifelse(relative_to_ref_i, "biomassplot", "biomassplot")
    suffix_relative <- paste0("_Combined", ifelse(relative_to_ref_i, "relative_", ""))
    pdf(file = paste0(PATHOUT, "/plot_NE_FMAX_", suffix_relative, ".pdf"), width = 14, height = 12)
    for (sp in species){
      combined_plot <- ggarrange(plotlist = c(p_list[[sp]], list(get(plot_biom)[[sp]])), ncol = 2, nrow = 3)
      print(combined_plot)
    }
    dev.off()
  }
}




## PLOT 2- 1-T & Vk	              2- N & T/Vk	      3- Ne & Ne/N ou Nb & Nb/N ------------------

traits_to_plot <- list("T&Vk" = c("generation", "mean_V"), "N&T/Vk" =c("N", "T/Vk"), "Ne&Ne/N"=c("Ne", "Nb", "/"))

data = dataRATIO_full_merged
data =  list_data_traits_ratio
for (relative_to_ref in c(TRUE)){
  if (!any(grepl(names(data), pattern="Ne"))) {
    tempo <- c(names(data[[1]])[1], names(data[[1]])[2], names(data[[1]])[3])
    data <- list(do.call(rbind, list(data[[1]][[names(data[[1]])[1]]], data[[2]][[names(data[[1]])[1]]])),
                 do.call(rbind, list(data[[1]][[names(data[[1]])[2]]], data[[2]][[names(data[[1]])[2]]])),
                 do.call(rbind, list(data[[1]][[names(data[[1]])[3]]], data[[2]][[names(data[[1]])[3]]])))
    names(data) <- tempo

  }
  scale_v = list()
  for (genetic_drift_data in (1:length(data))){
    p = list()
    for (traits_selected in names(traits_to_plot)) {
        # names_substr <- substr(summary$variable, nchar(summary$variable) - 3 + 1, nchar(summary$variable))

        # select dataset
        summary = data[[genetic_drift_data]]
        if (names(data)[genetic_drift_data] == "tempoNe") {
          if (ncol(data[["data_Nb"]])>6)  summary <- full_join(summary, data[["data_Nb"]][,c(1:4,grep(colnames(data[["data_Nb"]]), pattern="mean_V"))])
          if (ncol(data[["data_Nb"]])==6) summary <- rbind(summary, data[["data_Nb"]][grep(data[["data_Nb"]]$variable, pattern="mean_V"),])
        }
        if (names(data)[genetic_drift_data] == "data_Nb") {
          if (ncol(data[["data_Nb"]])>6)  summary <- full_join(summary, data[["tempoNe"]][,c(1:4,grep(colnames(data[["tempoNe"]]), pattern="generation"))])
          if (ncol(data[["data_Nb"]])==6) summary <- rbind(summary, data[["tempoNe"]][grep(data[["tempoNe"]]$variable, pattern="generation"),])
        }
        if (ncol(summary)>6) summary <- pivot_longer(summary, cols = colnames(summary)[-c(1:4)], names_to = "variable", values_to = "value")
        if (any(grepl(traits_selected, pattern="T/Vk"))) {
            df_trait_t <- estimateT_Vk(trait_data_full_i = summary)
            summary <- rbind(summary, df_trait_t[which(df_trait_t$variable == "T/Vk"),])
        }

        df_trait <- summary[as.vector(unlist(sapply(traits_to_plot[[traits_selected]], FUN = function(x) grep(summary$variable, pattern = x)))),]
        if (nrow(df_trait)<1) next
        var_i    <- unique(df_trait$variable)
        if (any(grepl(var_i, pattern=c("of")))) var_i <- var_i[-grep(var_i, pattern=c("of"))]
        if (any(grepl(var_i, pattern=c("egg")))) var_i <- var_i[-grep(var_i, pattern=c("egg"))]
        id_mean <- (grep(var_i, pattern="mean"))
        which_to_remove <- NULL
        if (any(grepl(var_i[id_mean], pattern=c("Numerator"))) & traits_selected != "N&T/Vk") which_to_remove <- grep(var_i[id_mean], pattern=c("Numerator"))
        if (any(grepl(var_i[id_mean], pattern=c("NT"))) & traits_selected != "N&T/Vk") which_to_remove <- grep(var_i[id_mean], pattern=c("NT"))
        if (length(which_to_remove)>0) var_i <- var_i[-id_mean[which_to_remove]]
        if (traits_selected == "N&T/Vk" & length(var_i)>2) var_i <- c(var_i[-c(grep(pattern="Nb", var_i), grep(pattern="Ne", var_i))])
        for (sp in unique(df_trait$Species)){
            df_trait_sp <- df_trait %>% filter(Species == sp)
            df_trait_sp <- df_trait_sp %>% filter(Fstart == 1)
            df_trait_sp <- df_trait_sp %>% filter(variable %in% var_i)
            if (sum(is.nan(df_trait_sp$value))>0) df_trait_sp <- df_trait_sp[-which(is.nan(df_trait_sp$value)), ]
            if (sum(is.infinite(df_trait_sp$value))>0) df_trait_sp <- df_trait_sp[-which(is.infinite(df_trait_sp$value)), ]

            # conditions to stop
            if (nrow(df_trait_sp)==0) next
            if (all(is.na(df_trait_sp$value))) {
              p[[sp]][[traits_selected]] <- NULL
              next}

            # Use the function
            if (length(var_i)>2){
              k=1
              id_to_plot = NULL
              for (i in 1:c((length(var_i)+1)/2)){
                id_to_plot = k:(k+1)
                if (traits_selected == "N&T/Vk") id_to_plot=c(which(var_i == "T/Vk"),
                                    ifelse(i==1, which(var_i != "T/Vk")[1], which(var_i != "T/Vk")[2]))
                data_split <- df_trait_sp[which(df_trait_sp$variable %in% var_i[id_to_plot] ),]
                if (all(is.na(data_split$value))) next
                result <- plot_trait_vs_fmax(df_trait_sp = data_split, trait_name = var_i[id_to_plot],
                                             Time1, Time2, relative_to_ref = relative_to_ref, spec=sp)
                p[[sp]][[paste(var_i[id_to_plot], collapse="")]]       <- result$plot
                scale_v[[sp]][[paste(var_i[id_to_plot], collapse="")]] <- result$scale
                k=k+2
              }
            } else {
              result <- plot_trait_vs_fmax(df_trait_sp = df_trait_sp, trait_name = var_i,
                                           Time1, Time2, relative_to_ref = relative_to_ref, spec=sp)
              p[[sp]][[traits_selected]]       <- result$plot
              scale_v[[sp]][[traits_selected]] <- result$scale
            }
        }#sp
    } # trait selected in the list
        trait_i = unique(var_i)
    if (length(unique(var_i))) trait_i = paste0(unique(var_i)[1], "_WITH_")
    plot_biom = ifelse(relative_to_ref, "rel_biomassplot", "biomassplot")
    plot_biom = ifelse(relative_to_ref, "biomassplot", "biomassplot")
    suffix_relative <- ifelse(relative_to_ref, "_RELATIVE", "")
    pdf(file = paste0(PATHOUT, "/plot_", names(data)[genetic_drift_data], "_T_Vk_N_", "vsFmax", suffix_relative, ".pdf"), width = 15, height = 10)
    par(mfrow=c(1,1))
    for (sp in species){
        combined_plot <- ggarrange(plotlist = c(p[[sp]], list(get(plot_biom)[[sp]])), ncol = 2, nrow=ceiling((length(p[[sp]])+1)/2))
        print(combined_plot)
    }
    dev.off()
  } # genetic drift dim
}# relative




## PLOT 3- 1-Z(a) & F(a)	        2- l(a) & Fec(a)	------------------
traits_to_plot <- list("surv_fec" = c("surv_curve", "fecundity"), "FecxSurv"=c("FecxSurv"), "Z_F" = c("Z", "F"), "Mnat_F" = c("Mnat", "F"))

for (relative_to_ref in c(TRUE)){
  scale_v = list()
  for (genetic_drift_data in (1:length(dataRATIO_full_merged_age[[2]]))[1]){
    p = list()
    for (traits_selected in names(traits_to_plot)) {
      cat(traits_selected, "\n")
      # select dataset
      summary = dataRATIO_full_merged_age[[2]][[genetic_drift_data]]
      summary$var_F_for_sp <- NULL
      if ("cohort_mean_of_Mnat_per_age" %in% summary$variable) summary[which(summary$variable == "cohort_mean_of_Mnat_per_age"),"variable"] <- "cohort_mean_of_Madd_per_age"
      if ("year_mean_of_Mnat_per_age" %in% summary$variable) summary[which(summary$variable == "year_mean_of_Mnat_per_age"),"variable"] <- "year_mean_of_Madd_per_age"
      if ("cohort_mean_of_Z_per_age" %in% summary$variable) summary <- summary[-which(summary$variable == "cohort_mean_of_Z_per_age"),]
      if ("year_mean_of_Z_per_age" %in% summary$variable) summary <- summary[-which(summary$variable == "year_mean_of_Z_per_age"),]
      if (ncol(summary)>8) summary <- pivot_longer(summary, cols = colnames(summary)[-c(1:4)], names_to = "variable", values_to = "value")
      if (any(grepl(traits_selected, pattern="T/Vk"))) {
        df_trait_t <- estimateT_Vk(trait_data_full_i = summary)
        summary <- rbind(summary, df_trait_t[which(df_trait_t$variable == "T/Vk"),])
      }
      if (any(grepl(traits_selected, pattern="Z"))) {
        df_trait_t <- estimateZ(trait_data_full_i = summary)
        summary <- rbind(summary, df_trait_t[which(df_trait_t$variable == "Z"),])
      }
      if (any(grepl(traits_selected, pattern="Mnat"))) {
        df_trait_t <- estimateMnat(trait_data_full_i = summary)
        summary <- rbind(summary, df_trait_t[which(df_trait_t$variable == "Mnat"),])
      }      
      if (any(grepl(traits_selected, pattern="Fec/Surv"))) {
        df_trait_t <- estimateFonS(trait_data_full_i = summary)
        summary <- rbind(summary, df_trait_t[which(df_trait_t$variable == "Fec/Surv"),])
      }
      if (any(grepl(traits_selected, pattern="FecxSurv"))) {
        df_trait_t <- estimateFS(trait_data_full_i = summary)
        summary <- rbind(summary, df_trait_t[which(df_trait_t$variable == "FecxSurv"),])
      }
      if (any(grepl(traits_selected, pattern="surv"))) {
        df_trait_t <- estimateSurv(trait_data_full_i = summary)
        summary <- rbind(summary, df_trait_t[which(df_trait_t$variable == "surv_curve"),])
        df_trait_t[which(df_trait_t$variable == "surv_curve"),"value"] <- log10(df_trait_t[which(df_trait_t$variable == "surv_curve"),"value"])
        df_trait_t[which(df_trait_t$variable == "surv_curve"),"variable"] <- "surv_curve (log10)"
      }

      df_trait <- summary[as.vector(unlist(sapply(traits_to_plot[[traits_selected]], FUN = function(x) grep(summary$variable, pattern = x)))),]
      if (nrow(df_trait)<1) next
      var_i    <- unique(df_trait$variable)

      # if (any(grepl(var_i, pattern=c("of")))) var_i <- var_i[-grep(var_i, pattern=c("of"))]
      # if (any(grepl(var_i, pattern=c("egg")))) var_i <- var_i[-grep(var_i, pattern=c("egg"))]
      # id_mean <- (grep(var_i, pattern="mean"))
      # which_to_remove <- NULL
      # if (any(grepl(var_i[id_mean], pattern=c("Numerator"))) & traits_selected != "N&T/Vk") which_to_remove <- grep(var_i[id_mean], pattern=c("Numerator"))
      # if (any(grepl(var_i[id_mean], pattern=c("NT"))) & traits_selected != "N&T/Vk") which_to_remove <- grep(var_i[id_mean], pattern=c("NT"))
      # if (length(which_to_remove)>0) var_i <- var_i[-id_mean[which_to_remove]]
      # if (traits_selected == "N&T/Vk" & length(var_i)>2) var_i <- c(var_i[-c(grep(pattern="Nb", var_i), grep(pattern="Ne", var_i))])
      
      for (sp in unique(df_trait$Species)){
        cat(sp)
        df_trait_sp <- df_trait %>% filter(Species == sp)
        df_trait_sp <- df_trait_sp %>% filter(Fstart == 1)
        df_trait_sp <- df_trait_sp %>% filter(variable %in% var_i)
        
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
            p[[sp]][[traits_selected]][[Fmax_i]] <- NULL
            next
          }

          # Use the function
          # if (length(var_i)>2){
          #   k=1
          #   id_to_plot = NULL
          #   for (i in 1:c((length(var_i)+1)/2)){
          #     id_to_plot = k:(k+1)
          #     if (traits_selected == "N&T/Vk") id_to_plot=c(which(var_i == "T/Vk"),
          #                                                   ifelse(i==1, which(var_i != "T/Vk")[1], which(var_i != "T/Vk")[2]))
          #     data_split <- df_trait_sp_a[which(df_trait_sp_a$variable %in% var_i[id_to_plot] ),]
          #     if (all(is.na(data_split$value))) next
          #     result <- plot_trait_vs_fmax(df_trait_sp = data_split, trait_name = var_i[id_to_plot],
          #                                  Time1, Time2, relative_to_ref = relative_to_ref, spec=sp, age_dim=T, reference_input=ref_data)
          #     p[[sp]][[paste(var_i[id_to_plot], collapse="")]][[Fmax_i]]       <- result$plot
          #     scale_v[[sp]][[paste(var_i[id_to_plot], collapse="")]][[Fmax_i]] <- result$scale
          #     k=k+2
          
          #   }
          # } else {
            result <- plot_trait_vs_fmax(df_trait_sp = df_trait_sp_a, trait_name = var_i,
                                         Time1, Time2, relative_to_ref = relative_to_ref,
                                         spec=sp, age_dim=T, reference_input=ref_data)
            p[[sp]][[traits_selected]][[Fmax_i]]       <- result$plot
            scale_v[[sp]][[traits_selected]][[Fmax_i]] <- result$scale
          # }
        }#age
      }#sp
    } # trait selected in the list
    plot_biom = ifelse(relative_to_ref, "rel_biomassplot", "biomassplot")
    plot_biom = ifelse(relative_to_ref, "biomassplot", "biomassplot")
    suffix_relative <- ifelse(relative_to_ref, "_RELATIVE", "")
    for (sp in species){
      if (!dir.exists(paste0(PATHOUT, "/byage", sp))) dir.create(paste0(PATHOUT, "/byage", sp))
      pdf(file = paste0(PATHOUT, "/byage", sp, "/ZFECSplot_", names(dataRATIO_full_merged_age[[2]])[genetic_drift_data], sp, "vsFmax", suffix_relative, ".pdf"), 
          width = 15, height = 8)
      if (is.null(p[[sp]])) next
      for (traits_selected in names(p[[sp]])) {
        plot_list <- c(p[[sp]][[traits_selected]], list(get(plot_biom)[[sp]]))
        plot_list <- lapply(plot_list, function(x) x + theme_minimal(base_size = 11) + theme(legend.position = "bottom"))
        combined_plot <- ggarrange(plotlist = plot_list, ncol = 3,
                                   nrow = ceiling((length(plot_list))/3),
                                   common.legend = TRUE, legend = "bottom")
        print(combined_plot)
      }
      dev.off()
    }
  } # genetic drift dim
}# relative


## PLOT 4- 1-Abund(a) & Egg(a)	------------------

traits_to_plot <- list("ab_Egg" = c("abundanceBy", "egg"), "abun_Biom" = c("abundanceBy", "biomassBy"), 
                       "abund_size" = c("abundanceBy", "meanSizeByAge"), "fec_weig"=c("fecundity", "weight", "abundanceBy"),
                       "fec_ab"=c("fecundity", "abundanceBy"))

for (relative_to_ref in c(TRUE)){
  scale_v = list()
  for (genetic_drift_data in (1:length(dataRATIO_full_merged_age[[2]]))[1]){
    p = list()
    for (traits_selected in names(traits_to_plot)) {
      cat(traits_selected, "\n")
      # select dataset
      summary = dataRATIO_full_merged_age[[2]][[genetic_drift_data]]
      summary$var_F_for_sp <- NULL
      if ("cohort_mean_of_Mnat_per_age" %in% summary$variable) summary[which(summary$variable == "cohort_mean_of_Mnat_per_age"),"variable"] <- "cohort_mean_of_Madd_per_age"
      if ("year_mean_of_Mnat_per_age" %in% summary$variable) summary[which(summary$variable == "year_mean_of_Mnat_per_age"),"variable"] <- "year_mean_of_Madd_per_age"
      if ("cohort_mean_of_Z_per_age" %in% summary$variable) summary <- summary[-which(summary$variable == "cohort_mean_of_Z_per_age"),]
      if ("year_mean_of_Z_per_age" %in% summary$variable) summary <- summary[-which(summary$variable == "year_mean_of_Z_per_age"),]
      if (ncol(summary)>8) summary <- pivot_longer(summary, cols = colnames(summary)[-c(1:4)], names_to = "variable", values_to = "value")
      if (any(grepl(traits_selected, pattern="T/Vk"))) {
        df_trait_t <- estimateT_Vk(trait_data_full_i = summary)
        summary <- rbind(summary, df_trait_t[which(df_trait_t$variable == "T/Vk"),])
      }
      if (any(grepl(traits_selected, pattern="Z"))) {
        df_trait_t <- estimateZ(trait_data_full_i = summary)
        summary <- rbind(summary, df_trait_t[which(df_trait_t$variable == "Z"),])
      }
      if (any(grepl(traits_selected, pattern="Mnat"))) {
        df_trait_t <- estimateMnat(trait_data_full_i = summary)
        summary <- rbind(summary, df_trait_t[which(df_trait_t$variable == "Mnat"),])
      }      
      if (any(grepl(traits_selected, pattern="surv"))) {
        df_trait_t <- estimateSurv(trait_data_full_i = summary)
        summary <- rbind(summary, df_trait_t[which(df_trait_t$variable == "surv_curve"),])
      }
      if (any(grepl(traits_selected, pattern="wei"))) {
        df_trait_t <- estimateWei(trait_data_full_i = summary)
        summary <- rbind(summary, df_trait_t[which(df_trait_t$variable == "weight"),])
      }
      
      df_trait <- summary[as.vector(unlist(sapply(traits_to_plot[[traits_selected]], FUN = function(x) grep(summary$variable, pattern = x)))),]
      if (nrow(df_trait)<1) next
      var_i    <- unique(df_trait$variable)
      
      # if (any(grepl(var_i, pattern=c("of")))) var_i <- var_i[-grep(var_i, pattern=c("of"))]
      # if (any(grepl(var_i, pattern=c("egg")))) var_i <- var_i[-grep(var_i, pattern=c("egg"))]
      # id_mean <- (grep(var_i, pattern="mean"))
      # which_to_remove <- NULL
      # if (any(grepl(var_i[id_mean], pattern=c("Numerator"))) & traits_selected != "N&T/Vk") which_to_remove <- grep(var_i[id_mean], pattern=c("Numerator"))
      # if (any(grepl(var_i[id_mean], pattern=c("NT"))) & traits_selected != "N&T/Vk") which_to_remove <- grep(var_i[id_mean], pattern=c("NT"))
      # if (length(which_to_remove)>0) var_i <- var_i[-id_mean[which_to_remove]]
      # if (traits_selected == "N&T/Vk" & length(var_i)>2) var_i <- c(var_i[-c(grep(pattern="Nb", var_i), grep(pattern="Ne", var_i))])
      
      for (sp in unique(df_trait$Species)){
        cat(sp)
        df_trait_sp <- df_trait %>% filter(Species == sp)
        df_trait_sp <- df_trait_sp %>% filter(Fstart == 1)
        df_trait_sp <- df_trait_sp %>% filter(variable %in% var_i)
        
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
            p[[sp]][[traits_selected]][[Fmax_i]] <- NULL
            next
          }
          
        # Use the function
          k=1
          id_to_plot = NULL
          result <- plot_trait_vs_fmax(df_trait_sp = df_trait_sp_a, trait_name = var_i,
                                       Time1, Time2, relative_to_ref = relative_to_ref, spec=sp, age_dim=T, reference_input=ref_data)
          p[[sp]][[traits_selected]][[Fmax_i]]       <- result$plot
          scale_v[[sp]][[traits_selected]][[Fmax_i]] <- result$scale
        }#age
      }#sp
    } # trait selected in the list
    plot_biom = ifelse(relative_to_ref, "rel_biomassplot", "biomassplot")
    plot_biom = ifelse(relative_to_ref, "biomassplot", "biomassplot")
    suffix_relative <- ifelse(relative_to_ref, "_RELATIVE", "")
    for (sp in species){
      if (!dir.exists(paste0(PATHOUT, "/byage", sp))) dir.create(paste0(PATHOUT, "/byage", sp))
      pdf(file = paste0(PATHOUT, "/byage", sp, "/ABUNDplot_", names(dataRATIO_full_merged_age[[2]])[genetic_drift_data], sp, "vsFmax", suffix_relative, ".pdf"), 
          width = 15, height = 8)
      if (is.null(p[[sp]])) next
      for (traits_selected in names(p[[sp]])) {
        plot_list <- c(p[[sp]][[traits_selected]], list(get(plot_biom)[[sp]]))
        plot_list <- lapply(plot_list, function(x) x + theme_minimal(base_size = 11) + theme(legend.position = "bottom"))
        combined_plot <- ggarrange(plotlist = plot_list, ncol = 3,
                                   nrow = ceiling((length(plot_list))/3),
                                   common.legend = TRUE, legend = "bottom")
        print(combined_plot)
      }
      dev.off()
    }
  } # genetic drift dim
}# relative


## PLOT 5- 1-Mpred(a) & F(a)	    2- Mfor(a) & F(a)	3- Mstarv(a) & F(a)	    4- Madd(a) & F(a)------------------

traits_to_plot <- list("Mpred&F"=c("Mpred", "F"), "Mfor&F"=c("Mfor", "F"), "Mstar&F"=c("Mstar", "F"), "Madd&F"=c("Madd", "F"))

for (relative_to_ref in c(TRUE)){
  scale_v = list()
  for (genetic_drift_data in (1:length(dataRATIO_full_merged_age[[2]]))[1]){
    p = list()
    for (traits_selected in names(traits_to_plot)) {
      cat(traits_selected, "\n")
      # select dataset
      summary = dataRATIO_full_merged_age[[2]][[genetic_drift_data]]
      summary$var_F_for_sp <- NULL
      if ("cohort_mean_of_Mnat_per_age" %in% summary$variable) summary[which(summary$variable == "cohort_mean_of_Mnat_per_age"),"variable"] <- "cohort_mean_of_Madd_per_age"
      if ("year_mean_of_Mnat_per_age" %in% summary$variable) summary[which(summary$variable == "year_mean_of_Mnat_per_age"),"variable"] <- "year_mean_of_Madd_per_age"
      if ("cohort_mean_of_Z_per_age" %in% summary$variable) summary <- summary[-which(summary$variable == "cohort_mean_of_Z_per_age"),]
      if ("year_mean_of_Z_per_age" %in% summary$variable) summary <- summary[-which(summary$variable == "year_mean_of_Z_per_age"),]
      if (ncol(summary)>8) summary <- pivot_longer(summary, cols = colnames(summary)[-c(1:4)], names_to = "variable", values_to = "value")
      if (any(grepl(traits_selected, pattern="T/Vk"))) {
        df_trait_t <- estimateT_Vk(trait_data_full_i = summary)
        summary <- rbind(summary, df_trait_t[which(df_trait_t$variable == "T/Vk"),])
      }
      if (any(grepl(traits_selected, pattern="Z"))) {
        df_trait_t <- estimateZ(trait_data_full_i = summary)
        summary <- rbind(summary, df_trait_t[which(df_trait_t$variable == "Z"),])
      }

      df_trait <- summary[as.vector(unlist(sapply(traits_to_plot[[traits_selected]], FUN = function(x) grep(summary$variable, pattern = x)))),]
      if (nrow(df_trait)<1) next
      var_i    <- unique(df_trait$variable)

      # if (any(grepl(var_i, pattern=c("of")))) var_i <- var_i[-grep(var_i, pattern=c("of"))]
      # if (any(grepl(var_i, pattern=c("egg")))) var_i <- var_i[-grep(var_i, pattern=c("egg"))]
      # id_mean <- (grep(var_i, pattern="mean"))
      # which_to_remove <- NULL
      # if (any(grepl(var_i[id_mean], pattern=c("Numerator"))) & traits_selected != "N&T/Vk") which_to_remove <- grep(var_i[id_mean], pattern=c("Numerator"))
      # if (any(grepl(var_i[id_mean], pattern=c("NT"))) & traits_selected != "N&T/Vk") which_to_remove <- grep(var_i[id_mean], pattern=c("NT"))
      # if (length(which_to_remove)>0) var_i <- var_i[-id_mean[which_to_remove]]
      # if (traits_selected == "N&T/Vk" & length(var_i)>2) var_i <- c(var_i[-c(grep(pattern="Nb", var_i), grep(pattern="Ne", var_i))])


      for (sp in unique(df_trait$Species)){
        cat(sp)
        df_trait_sp <- df_trait %>% filter(Species == sp)
        df_trait_sp <- df_trait_sp %>% filter(Fstart == 1)
        df_trait_sp <- df_trait_sp %>% filter(variable %in% var_i)
        
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
            p[[sp]][[traits_selected]][[Fmax_i]] <- NULL
            next
          }
          
          # Use the function
          if (length(var_i)>2){
            k=1
            id_to_plot = NULL
            for (i in 1:c((length(var_i)+1)/2)){
              id_to_plot = k:(k+1)
              if (traits_selected == "N&T/Vk") id_to_plot=c(which(var_i == "T/Vk"),
                                                            ifelse(i==1, which(var_i != "T/Vk")[1], which(var_i != "T/Vk")[2]))
              data_split <- df_trait_sp_a[which(df_trait_sp_a$variable %in% var_i[id_to_plot] ),]
              if (all(is.na(data_split$value))) next
              result <- plot_trait_vs_fmax(df_trait_sp = data_split, trait_name = var_i[id_to_plot],
                                           Time1, Time2, relative_to_ref = relative_to_ref, spec=sp, age_dim=T, reference_input=ref_data)
              p[[sp]][[paste(var_i[id_to_plot], collapse="")]][[Fmax_i]]       <- result$plot
              scale_v[[sp]][[paste(var_i[id_to_plot], collapse="")]][[Fmax_i]] <- result$scale
              k=k+2
            }
          } else {
            result <- plot_trait_vs_fmax(df_trait_sp = df_trait_sp_a, trait_name = var_i,
                                         Time1, Time2, relative_to_ref = relative_to_ref, spec=sp, age_dim=T, reference_input=ref_data)
            p[[sp]][[traits_selected]][[Fmax_i]]       <- result$plot
            scale_v[[sp]][[traits_selected]][[Fmax_i]] <- result$scale
          }
        }#age
      }#sp
      
      # for (sp in unique(df_trait$Species)){
      #   cat(sp)
      #   df_trait_sp <- df_trait %>% filter(Species == sp)
      #   df_trait_sp <- df_trait_sp %>% filter(Fstart == 1)
      #   df_trait_sp <- df_trait_sp %>% filter(variable %in% var_i)
      #   for (age_i in unique(df_trait_sp$age)){
      #     cat(age_i, " ")
      #     df_trait_sp_a <- df_trait_sp %>% filter(age == age_i)
      # 
      #     if (sum(is.nan(df_trait_sp_a$value))>0) df_trait_sp_a <- df_trait_sp_a[-which(is.nan(df_trait_sp_a$value)), ]
      #     if (sum(is.infinite(df_trait_sp_a$value))>0) df_trait_sp_a <- df_trait_sp_a[-which(is.infinite(df_trait_sp_a$value)), ]
      # 
      #     # conditions to stop
      #     if (nrow(df_trait_sp_a)==0) next
      #     if (all(is.na(df_trait_sp_a$value) | df_trait_sp_a$value == 0)) {
      #       p[[sp]][[traits_selected]][[age_i]] <- NULL
      #       next
      #     }
      # 
      #     # Use the function
      #     if (length(var_i)>2){
      #       k=1
      #       id_to_plot = NULL
      #       for (i in 1:c((length(var_i)+1)/2)){
      #         id_to_plot = k:(k+1)
      #         if (traits_selected == "N&T/Vk") id_to_plot=c(which(var_i == "T/Vk"),
      #                                                       ifelse(i==1, which(var_i != "T/Vk")[1], which(var_i != "T/Vk")[2]))
      #         data_split <- df_trait_sp_a[which(df_trait_sp_a$variable %in% var_i[id_to_plot] ),]
      #         if (all(is.na(data_split$value))) next
      #         result <- plot_trait_vs_fmax(df_trait_sp = data_split, trait_name = var_i[id_to_plot],
      #                                      Time1, Time2, relative_to_ref = relative_to_ref, spec=sp)
      #         p[[sp]][[paste(var_i[id_to_plot], collapse="")]][[age_i]]       <- result$plot
      #         scale_v[[sp]][[paste(var_i[id_to_plot], collapse="")]][[age_i]] <- result$scale
      #         k=k+2
      #       }
      #     } else {
      #       result <- plot_trait_vs_fmax(df_trait_sp = df_trait_sp_a, trait_name = var_i,
      #                                    Time1, Time2, relative_to_ref = relative_to_ref, spec=sp)
      #       p[[sp]][[traits_selected]][[age_i]]       <- result$plot
      #       scale_v[[sp]][[traits_selected]][[age_i]] <- result$scale
      #     }
      #   }#age
      # }#sp
    } # trait selected in the list
    plot_biom = ifelse(relative_to_ref, "rel_biomassplot", "biomassplot")
    plot_biom = ifelse(relative_to_ref, "biomassplot", "biomassplot")
    suffix_relative <- ifelse(relative_to_ref, "_RELATIVE", "")
    for (sp in species){
      if (!dir.exists(paste0(PATHOUT, "/byage", sp))) dir.create(paste0(PATHOUT, "/byage", sp))
      pdf(file = paste0(PATHOUT, "/byage", sp, "/MORTplot_", names(dataRATIO_full_merged_age[[2]])[genetic_drift_data], sp, "vsFmax", suffix_relative, ".pdf"),
          width = 15, height = 8)
      if (is.null(p[[sp]])) next
      for (traits_selected in names(p[[sp]])) {
        plot_list <- c(p[[sp]][[traits_selected]], list(get(plot_biom)[[sp]]))
        plot_list <- lapply(plot_list, function(x) x + theme_minimal(base_size = 11) + theme(legend.position = "bottom"))
        combined_plot <- ggarrange(plotlist = plot_list, ncol = 3,
                                   nrow = ceiling((length(plot_list))/3),
                                   common.legend = TRUE, legend = "bottom")
        print(combined_plot)
      }
      dev.off()
    }
  } # genetic drift dim
}# relative

traits_to_plot <- list("Mpred&Mnat"=c("Mpred", "Mnat"), "Mfor&Mnat"=c("Mfor", "Mnat"), "Mstar&Mnat"=c("Mstar", "Mnat"), "Madd&Mnat"=c("Madd", "Mnat"))

for (relative_to_ref in c(TRUE)){
  scale_v = list()
  for (genetic_drift_data in (1:length(dataRATIO_full_merged_age[[2]]))[1]){
    p = list()
    for (traits_selected in names(traits_to_plot)) {
      cat(traits_selected, "\n")
      # select dataset
      summary = dataRATIO_full_merged_age[[2]][[genetic_drift_data]]
      summary$var_F_for_sp <- NULL
      if ("cohort_mean_of_Mnat_per_age" %in% summary$variable) summary[which(summary$variable == "cohort_mean_of_Mnat_per_age"),"variable"] <- "cohort_mean_of_Madd_per_age"
      if ("year_MnatbyY" %in% summary$variable) summary[which(summary$variable == "year_MnatbyY"),"variable"] <- "year_MaddbyY"
      if ("cohort_mean_of_Z_per_age" %in% summary$variable) summary <- summary[-which(summary$variable == "cohort_mean_of_Z_per_age"),]
      if ("year_mean_of_Z_per_age" %in% summary$variable) summary <- summary[-which(summary$variable == "year_mean_of_Z_per_age"),]
      if (ncol(summary)>8) summary <- pivot_longer(summary, cols = colnames(summary)[-c(1:4)], names_to = "variable", values_to = "value")
      if (any(grepl(traits_selected, pattern="T/Vk"))) {
        df_trait_t <- estimateT_Vk(trait_data_full_i = summary)
        summary <- rbind(summary, df_trait_t[which(df_trait_t$variable == "T/Vk"),])
      }
      if (any(grepl(traits_selected, pattern="Z"))) {
        df_trait_t <- estimateZ(trait_data_full_i = summary)
        summary <- rbind(summary, df_trait_t[which(df_trait_t$variable == "Z"),])
      }
      if (any(grepl(traits_selected, pattern="Mnat"))) {
        df_trait_t <- estimateMnat(trait_data_full_i = summary)
        summary <- rbind(summary, df_trait_t[which(df_trait_t$variable == "Mnat"),])
      }
      df_trait <- summary[as.vector(unlist(sapply(traits_to_plot[[traits_selected]], FUN = function(x) grep(summary$variable, pattern = x)))),]
      if (nrow(df_trait)<1) next
      var_i    <- unique(df_trait$variable)
      
      # if (any(grepl(var_i, pattern=c("of")))) var_i <- var_i[-grep(var_i, pattern=c("of"))]
      # if (any(grepl(var_i, pattern=c("egg")))) var_i <- var_i[-grep(var_i, pattern=c("egg"))]
      # id_mean <- (grep(var_i, pattern="mean"))
      # which_to_remove <- NULL
      # if (any(grepl(var_i[id_mean], pattern=c("Numerator"))) & traits_selected != "N&T/Vk") which_to_remove <- grep(var_i[id_mean], pattern=c("Numerator"))
      # if (any(grepl(var_i[id_mean], pattern=c("NT"))) & traits_selected != "N&T/Vk") which_to_remove <- grep(var_i[id_mean], pattern=c("NT"))
      # if (length(which_to_remove)>0) var_i <- var_i[-id_mean[which_to_remove]]
      # if (traits_selected == "N&T/Vk" & length(var_i)>2) var_i <- c(var_i[-c(grep(pattern="Nb", var_i), grep(pattern="Ne", var_i))])
      
      for (sp in unique(df_trait$Species)){
        cat(sp)
        df_trait_sp <- df_trait %>% filter(Species == sp)
        df_trait_sp <- df_trait_sp %>% filter(Fstart == 1)
        df_trait_sp <- df_trait_sp %>% filter(variable %in% var_i)
        
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
            p[[sp]][[traits_selected]][[Fmax_i]] <- NULL
            next
          }
          
          # Use the function
          if (length(var_i)>2){
            k=1
            id_to_plot = NULL
            for (i in 1:c((length(var_i)+1)/2)){
              id_to_plot = k:(k+1)
              if (traits_selected == "N&T/Vk") id_to_plot=c(which(var_i == "T/Vk"),
                                                            ifelse(i==1, which(var_i != "T/Vk")[1], which(var_i != "T/Vk")[2]))
              data_split <- df_trait_sp_a[which(df_trait_sp_a$variable %in% var_i[id_to_plot] ),]
              if (all(is.na(data_split$value))) next
              result <- plot_trait_vs_fmax(df_trait_sp = data_split, trait_name = var_i[id_to_plot],
                                           Time1, Time2, relative_to_ref = relative_to_ref, spec=sp, age_dim=T, reference_input=ref_data)
              p[[sp]][[paste(var_i[id_to_plot], collapse="")]][[Fmax_i]]       <- result$plot
              scale_v[[sp]][[paste(var_i[id_to_plot], collapse="")]][[Fmax_i]] <- result$scale
              k=k+2
            }
          } else {
            result <- plot_trait_vs_fmax(df_trait_sp = df_trait_sp_a, trait_name = var_i,
                                         Time1, Time2, relative_to_ref = relative_to_ref, spec=sp, age_dim=T, reference_input=ref_data)
            p[[sp]][[traits_selected]][[Fmax_i]]       <- result$plot
            scale_v[[sp]][[traits_selected]][[Fmax_i]] <- result$scale
          }
        }#age
      }#sp
      
      # for (sp in unique(df_trait$Species)){
      #   cat(sp)
      #   df_trait_sp <- df_trait %>% filter(Species == sp)
      #   df_trait_sp <- df_trait_sp %>% filter(Fstart == 1)
      #   df_trait_sp <- df_trait_sp %>% filter(variable %in% var_i)
      #   for (age_i in unique(df_trait_sp$age)){
      #     cat(age_i, " ")
      #     df_trait_sp_a <- df_trait_sp %>% filter(age == age_i)
      #     
      #     if (sum(is.nan(df_trait_sp_a$value))>0) df_trait_sp_a <- df_trait_sp_a[-which(is.nan(df_trait_sp_a$value)), ]
      #     if (sum(is.infinite(df_trait_sp_a$value))>0) df_trait_sp_a <- df_trait_sp_a[-which(is.infinite(df_trait_sp_a$value)), ]
      #     
      #     # conditions to stop
      #     if (nrow(df_trait_sp_a)==0) next
      #     if (all(is.na(df_trait_sp_a$value) | df_trait_sp_a$value == 0)) {
      #       p[[sp]][[traits_selected]][[age_i]] <- NULL
      #       next
      #     }
      #     
      #     # Use the function
      #     if (length(var_i)>2){
      #       k=1
      #       id_to_plot = NULL
      #       for (i in 1:c((length(var_i)+1)/2)){
      #         id_to_plot = k:(k+1)
      #         if (traits_selected == "N&T/Vk") id_to_plot=c(which(var_i == "T/Vk"),
      #                                                       ifelse(i==1, which(var_i != "T/Vk")[1], which(var_i != "T/Vk")[2]))
      #         data_split <- df_trait_sp_a[which(df_trait_sp_a$variable %in% var_i[id_to_plot] ),]
      #         if (all(is.na(data_split$value))) next
      #         result <- plot_trait_vs_fmax(df_trait_sp = data_split, trait_name = var_i[id_to_plot],
      #                                      Time1, Time2, relative_to_ref = relative_to_ref, spec=sp)
      #         p[[sp]][[paste(var_i[id_to_plot], collapse="")]][[age_i]]       <- result$plot
      #         scale_v[[sp]][[paste(var_i[id_to_plot], collapse="")]][[age_i]] <- result$scale
      #         k=k+2
      #       }
      #     } else {
      #       result <- plot_trait_vs_fmax(df_trait_sp = df_trait_sp_a, trait_name = var_i,
      #                                    Time1, Time2, relative_to_ref = relative_to_ref, spec=sp)
      #       p[[sp]][[traits_selected]][[age_i]]       <- result$plot
      #       scale_v[[sp]][[traits_selected]][[age_i]] <- result$scale
      #     }
      #   }#age
      # }#sp
    } # trait selected in the list
    plot_biom = ifelse(relative_to_ref, "rel_biomassplot", "biomassplot")
    plot_biom = ifelse(relative_to_ref, "biomassplot", "biomassplot")
    suffix_relative <- ifelse(relative_to_ref, "_RELATIVE", "")
    for (sp in species){
      if (!dir.exists(paste0(PATHOUT, "/byage", sp))) dir.create(paste0(PATHOUT, "/byage", sp))
      pdf(file = paste0(PATHOUT, "/byage", sp, "/MORTMNATplot_", names(dataRATIO_full_merged_age[[2]])[genetic_drift_data], sp, "vsFmax", suffix_relative, ".pdf"),
          width = 15, height = 8)
      if (is.null(p[[sp]])) next
      for (traits_selected in names(p[[sp]])) {
        plot_list <- c(p[[sp]][[traits_selected]], list(get(plot_biom)[[sp]]))
        plot_list <- lapply(plot_list, function(x) x + theme_minimal(base_size = 11) + theme(legend.position = "bottom"))
        combined_plot <- ggarrange(plotlist = plot_list, ncol = 3,
                                   nrow = ceiling((length(plot_list))/3),
                                   common.legend = TRUE, legend = "bottom")
        print(combined_plot)
      }
      dev.off()
    }
  } # genetic drift dim
}# relative

## PLOT 6- DECOUPLING Vk elements -----------------------

traits_to_plot <- list("Dead_cumFecx" = c("Dead", "cumfecundity", "Deadxcumfecundity"), 
                       "Dead_cumFecx^2" = c("Dead", "cumfecunditysquare", "Deadxcumfecunditysquare"), 
                       "sumDead_cumFecx^2" = c("Dead", "deltasquare", "Deadxcumfecundity", "Deadxcumfecunditysquare", "Deadxcumfecundity\nsumDeadxcumfecunditysquare"))

traits_to_plot <- list("Dead_cumFecx" = c("Dead", "cumfecundity", "Deadxcumfecundity"), 
                       "Dead_k_pointdeltax^2" = c("Dead", "k_pointdeltasquare", "Deadxk_pointdeltasquare"), 
                       "sumDead_k_pointdeltax^2Div" = c("Dead", "cumfecunditysquare", "k_pointdeltasquare", "Deadxk_pointdeltasquare", "Deadxcumfecunditysquare", 
                                                     "Deadxcumfecundity\nsumDeadxk_pointdeltasquare", "Deadxcumfecundity\nsumDeadxk_pointdeltasquare\nDivage1abundanceBy"))

for (relative_to_ref in c(TRUE)){
  scale_v = list()
  for (genetic_drift_data in (1:length(dataRATIO_full_merged_age[[2]]))[1]){
    p = list()
    for (traits_selected in names(traits_to_plot)) {
      cat(traits_selected, "\n")
      # select dataset
      summary = dataRATIO_full_merged_age[[2]][[genetic_drift_data]]
      summary$age <- as.numeric(as.matrix(summary$age))
      summary$var_F_for_sp <- NULL
      if ("cohort_mean_of_Mnat_per_age" %in% summary$variable) summary[which(summary$variable == "cohort_mean_of_Mnat_per_age"),"variable"] <- "cohort_mean_of_Madd_per_age"
      if ("year_mean_of_Mnat_per_age" %in% summary$variable) summary[which(summary$variable == "year_mean_of_Mnat_per_age"),"variable"] <- "year_mean_of_Madd_per_age"
      if ("cohort_mean_of_Z_per_age" %in% summary$variable) summary <- summary[-which(summary$variable == "cohort_mean_of_Z_per_age"),]
      if ("year_mean_of_Z_per_age" %in% summary$variable) summary <- summary[-which(summary$variable == "year_mean_of_Z_per_age"),]
      if (ncol(summary)>8) summary <- pivot_longer(summary, cols = colnames(summary)[-c(1:4)], names_to = "variable", values_to = "value")
      if (any(grepl(traits_selected, pattern="T/Vk"))) {
        df_trait_t <- estimateT_Vk(trait_data_full_i = summary)
        summary <- rbind(summary, df_trait_t[which(df_trait_t$variable == "T/Vk"),])
      }
      if (any(grepl(traits_selected, pattern="Dead"))) {
        df_trait_t <- estimateDead(trait_data_full_i = summary)
        df_trait_t[which(is.na(df_trait_t$value)),"value"] <- 0
        summary    <- rbind(summary, df_trait_t[which(df_trait_t$variable == "Dead"),])
      }
      if (any(grepl(traits_selected, pattern="k_point"))) {
        # kx_bar
        var="cumfecundity"
        df_trait_t <- estimateCumSum(trait_data_full_i = summary, var=var)
        summary    <- rbind(summary, df_trait_t)
        # Death*kx_bar (=cumfecundity)
        var="Deadxcumfecundity"
        for (v in var[ifelse(length(var)>1, -length(var), 1)]){
          df_trait_t <- estimateX(trait_data_full_i = summary, var=v)
          summary    <- rbind(summary, df_trait_t)
        }
        # colSums Death*kx_bar
        df_trait_t <- estimatecolSum(trait_data_full_i = summary, var="colSumDeadxcumfecundity")
        summary    <- rbind(summary, df_trait_t)
        # Age 1 abund
        df_trait_t <- estimateage1(trait_data_full_i = summary, var="age1abundanceBy")
        summary    <- rbind(summary, df_trait_t)
        # div colSums Death*kx_bar by Age 1 abund
        df_trait_t <- estimateDiv(trait_data_full_i = summary, var="colSumDeadxcumfecundity\nDivage1abundanceBy")
        summary    <- rbind(summary, df_trait_t)
        # store
        summary[which(summary$variable == "colSumDeadxcumfecundity\nDivage1abundanceBy"),"variable"] <- "k_point"
      }
      if (any(grepl(traits_selected, pattern="delta"))) {
        var = "cumfecundity\nsumminusk_point"
        minusk_point <- summary[which(summary$variable == "k_point"), ]
        minusk_point$value <- -minusk_point$value
        minusk_point$variable <- "minusk_point"
        summary <- rbind(summary, minusk_point)
        df_trait_t <- estimateSum(trait_data_full_i = summary, var=var)
        summary <- rbind(summary, df_trait_t)
        summary[which(summary$variable == "cumfecundity\nsumminusk_point"),"variable"] <- "k_pointdelta"
      }
      if (any(grepl(traits_selected, pattern="Z"))) {
        df_trait_t <- estimateZ(trait_data_full_i = summary)
        summary <- rbind(summary, df_trait_t[which(df_trait_t$variable == "Z"),])
      }
      if (any(grepl(traits_selected, pattern="Mnat"))) {
        df_trait_t <- estimateMnat(trait_data_full_i = summary)
        summary <- rbind(summary, df_trait_t[which(df_trait_t$variable == "Mnat"),])
      }      
      if (any(grepl(traits_selected, pattern="surv"))) {
        df_trait_t <- estimateSurv(trait_data_full_i = summary)
        summary    <- rbind(summary, df_trait_t[which(df_trait_t$variable == "surv_curve"),])
        df_trait_t[which(df_trait_t$variable == "surv_curve"),"value"] <- log10(df_trait_t[which(df_trait_t$variable == "surv_curve"),"value"])
        df_trait_t[which(df_trait_t$variable == "surv_curve"),"variable"] <- "surv_curve (log10)"
      }
      if (any(grepl(traits_selected, pattern="cum"))) {
        var=traits_to_plot[[traits_selected]][grep(traits_to_plot[[traits_selected]], pattern="cum")][1]
        if(!(var %in% summary$variable)) {
          df_trait_t <- estimateCumSum(trait_data_full_i = summary, var=var)
          summary    <- rbind(summary, df_trait_t)
        }
      }
      if (any(grepl(traits_selected, pattern="\\^2"))) {
        var=traits_to_plot[[traits_selected]][grep(traits_to_plot[[traits_selected]], pattern="square")]
        for (v in var[ifelse(length(var)>1, -length(var), 1)]){
          if(!(v %in% summary$variable) & length(str_subset(pattern = c("[x]"), string =  v, negate=T))>0) {
            df_trait_t <- estimateSquare(trait_data_full_i = summary, var=v)
            summary    <- rbind(summary, df_trait_t)
          }
        }
      }
      if (any(grepl(traits_selected, pattern="x"))) {
        var=traits_to_plot[[traits_selected]][grep(traits_to_plot[[traits_selected]], pattern="x")]
        for (v in var[ifelse(length(var)>1, -length(var), 1)]){
          if(!(v %in% summary$variable) & length(str_subset(pattern = c("\n"), string =  v, negate=T))>0) {
            df_trait_t <- estimateX(trait_data_full_i = summary, var=v)
            summary    <- rbind(summary, df_trait_t)
          }
        }
      }
      if (any(grepl(traits_selected, pattern="sum"))) {
        var=traits_to_plot[[traits_selected]][grep(traits_to_plot[[traits_selected]], pattern="sum")]
        for (v in var[ifelse(length(var)>1, -length(var), 1)]){
          if(!(v %in% summary$variable)) {
            df_trait_t <- estimateSum(trait_data_full_i = summary, var=v)
            summary    <- rbind(summary, df_trait_t)
          }
        }
      }
      if (any(grepl(traits_selected, pattern="Div"))) {
        var=traits_to_plot[[traits_selected]][grep(traits_to_plot[[traits_selected]], pattern="Div")]
        df_trait_t <- estimateDiv(trait_data_full_i = summary, var=var)
        summary    <- rbind(summary, df_trait_t)
      }
      
      df_trait <- summary[as.vector(unlist(sapply(traits_to_plot[[traits_selected]], FUN = function(x) grep(summary$variable, pattern = x)))),]
      if (nrow(df_trait)<1) next
      var_i    <- unique(df_trait$variable)
      
      for (sp in unique(df_trait$Species)){
        cat(sp)
        df_trait_sp <- df_trait %>% filter(Species == sp)
        df_trait_sp <- df_trait_sp %>% filter(Fstart == 1)
        df_trait_sp <- df_trait_sp %>% filter(variable %in% var_i)
        
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
            p[[sp]][[traits_selected]][[Fmax_i]] <- NULL
            next
          }
          if (length(var_i)>3){
            k=1
            id_to_plot = NULL
            for (i in 1:c((length(var_i)+1)/2)){
              id_to_plot = c(k, ifelse((k+1)>length(var_i), 1, k+1))
              if (traits_selected == "N&T/Vk") id_to_plot=c(which(var_i == "T/Vk"),
                                                            ifelse(i==1, which(var_i != "T/Vk")[1], which(var_i != "T/Vk")[2]))
              data_split <- df_trait_sp_a[which(df_trait_sp_a$variable %in% var_i[id_to_plot] ),]
              ref_data_split <- ref_data[which(ref_data$variable %in% var_i[id_to_plot]),]
              if (all(is.na(data_split$value))) next
              result <- plot_trait_vs_fmax(df_trait_sp = data_split, trait_name = var_i[id_to_plot],
                                           Time1, Time2, relative_to_ref = relative_to_ref, spec=sp, 
                                           age_dim=T, reference_input=ref_data_split)
              p[[sp]][[paste(var_i[id_to_plot], collapse="")]][[Fmax_i]]         <- result$plot
              scale_v[[sp]][[paste(var_i[id_to_plot], collapse="")]][[Fmax_i]]   <- result$scale
              k=k+2
            }
          } else {
            result <- plot_trait_vs_fmax(df_trait_sp = df_trait_sp_a, trait_name = var_i,
                                         Time1, Time2, relative_to_ref = relative_to_ref,
                                         spec=sp, age_dim=T, reference_input=ref_data)
            p[[sp]][[traits_selected]][[Fmax_i]]         <- result$plot
            scale_v[[sp]][[traits_selected]][[Fmax_i]]   <- result$scale
          }
          # result <- plot_trait_vs_fmax(df_trait_sp = df_trait_sp_a, trait_name = var_i,
          #                              Time1, Time2, relative_to_ref = relative_to_ref,
          #                              spec=sp, age_dim=T, reference_input=ref_data)
          # p[[sp]][[traits_selected]][[Fmax_i]]       <- result$plot
          # scale_v[[sp]][[traits_selected]][[Fmax_i]] <- result$scale
          # }
          
          if ((traits_selected) == "sumDead_k_pointdeltax^2Div") {
            id_to_plot = which(var_i=="Deadxcumfecundity\nsumDeadxk_pointdeltasquare\nDivage1abundanceBy")
            data_split <- df_trait_sp_a[which(df_trait_sp_a$variable %in% var_i[id_to_plot] ),]
            ref_data_split <- ref_data[which(ref_data$variable %in% var_i[id_to_plot]),]
            if (all(is.na(data_split$value))) next
            result <- plot_trait_vs_fmax(df_trait_sp = data_split[which(data_split$variable == var_i[id_to_plot]), ], trait_name = var_i[id_to_plot],
                                         Time1, Time2, relative_to_ref = relative_to_ref, spec=sp, 
                                         age_dim=T, reference_input=ref_data_split)
            vk_point <- data_split[which(data_split$variable == var_i[id_to_plot]), ] %>%
              group_by(replicate) %>%
              summarize(vk_point = sum(value))
            vk_point <- mean(vk_point$vk_point, na.rm=T)
            result$plot <- result$plot + 
              annotate(paste0="Vk_point = ", vk_point)
              
            p[[sp]][[paste(var_i[id_to_plot], collapse="")]][[Fmax_i]]         <- result$plot
            scale_v[[sp]][[paste(var_i[id_to_plot], collapse="")]][[Fmax_i]]   <- result$scale
            k=k+2
          }
        }#age
      }#sp
    } # trait selected in the list
    plot_biom = ifelse(relative_to_ref, "rel_biomassplot", "biomassplot")
    plot_biom = ifelse(relative_to_ref, "biomassplot", "biomassplot")
    suffix_relative <- ifelse(relative_to_ref, "_RELATIVE", "")
    for (sp in species){
      if (!dir.exists(paste0(PATHOUT, "/byage", sp))) dir.create(paste0(PATHOUT, "/byage", sp))
      pdf(file = paste0(PATHOUT, "/byage", sp, "/Vk_par_plot_", names(dataRATIO_full_merged_age[[2]])[genetic_drift_data], 
                        sp, "vsFmax", suffix_relative, ".pdf"), 
          width = 15, height = 8)
      if (is.null(p[[sp]])) next
      for (traits_selected in names(p[[sp]])) {
        plot_list <- c(p[[sp]][[traits_selected]], list(get(plot_biom)[[sp]]))
        # plot_list <- lapply(plot_list, function(x) x + theme_minimal(base_size = 11) + theme(legend.position = "bottom"))
        combined_plot <- ggarrange(plotlist = plot_list, ncol = 3,
                                   nrow = ceiling((length(plot_list))/3),
                                   common.legend = TRUE, legend = "bottom")
        print(combined_plot)
      }
      dev.off()
    }
  } # genetic drift dim
}# relative



 # FAIL
## PLOT 7- DECOUPLING Vk elements ------------------

only_Ne_traits


traits_to_plot <- list("Dead_cumFecx" = c("Dead", "cumfecundity", "Deadxcumfecundity"), 
                       "Dead_k_pointdeltax^2" = c("Dead", "k_pointdeltasquare", "Deadxk_pointdeltasquare"), 
                       "sumDead_k_pointdeltax^2Div" = c("Dead", "cumfecunditysquare", "k_pointdeltasquare", "Deadxk_pointdeltasquare", "Deadxcumfecunditysquare", 
                                                        "Deadxcumfecundity\nsumDeadxk_pointdeltasquare", "Deadxcumfecundity\nsumDeadxk_pointdeltasquare\nDivage1abundanceBy"))

# "Dead_eggNumber" = c("Dead", "eggNumber"), "ab_Div_egg" = c("abundanceBy", "eggNumber\nDivabundanceByAge"), 
# "surv" = c("survivorship_curve", "survival_rate"), 
                       #"V^2" = c("Vsquare")
traits_to_plot <- list("sumV_k_barxDeath^2age1DivCS" = c("Vsquare", "VsquarexDead", "VxDead", "age1abundanceBy", "VsquarexDead\nsumVxDead", 
                                                         "VsquarexDead\nsumVxDead\nDivage1abundanceBy", "colSumVsquarexDead\nsumVxDead\nDivage1abundanceBy"), 
                       # "VxDead" = c("V", "VxDead"),
                       "Dead_eggNumber" = c("Dead", "eggNumber"), 
                       "ab_Div_egg" = c("abundanceBy", "eggNumber\nDivabundanceByAge"), 
                       "surv" = c("survivorship_curve", "survival_rate")#, 
)

for (relative_to_ref in c(FALSE)){
  scale_v = list()
  p = list()
  for (traits_selected in names(traits_to_plot)) {
    traits_to_plot2 <- traits_to_plot
    cat(traits_selected, "\n")
    # select dataset
    summary = only_Ne_traits
    summary$var_F_for_sp <- NULL
    if ("cohort_mean_of_Mnat_per_age" %in% summary$variable) summary[which(summary$variable == "cohort_mean_of_Mnat_per_age"),"variable"] <- "cohort_mean_of_Madd_per_age"
    if ("year_mean_of_Mnat_per_age" %in% summary$variable) summary[which(summary$variable == "year_mean_of_Mnat_per_age"),"variable"] <- "year_mean_of_Madd_per_age"
    if ("cohort_mean_of_Z_per_age" %in% summary$variable) summary <- summary[-which(summary$variable == "cohort_mean_of_Z_per_age"),]
    if ("year_mean_of_Z_per_age" %in% summary$variable) summary <- summary[-which(summary$variable == "year_mean_of_Z_per_age"),]
    if (ncol(summary)>8) summary <- pivot_longer(summary, cols = colnames(summary)[-c(1:4)], names_to = "variable", values_to = "value")
    if (any(grepl(traits_selected, pattern="\\^2"))) {
      var=traits_to_plot2[[traits_selected]][grep(traits_to_plot2[[traits_selected]], pattern="square")]
      for (v in var[1]){
        if(!(v %in% summary$variable) & length(str_subset(pattern = c("[x]"), string =  v, negate=T))>0) {
          if (v == "Vsquare"){ # it is not Vx that has to be ² to form Delta², but (kx-k_bar)² normally equal to (Vx-1)²
            tempo      <- summary[which(summary$variable == "Vx"),]
            tempo$variable <- paste0("minus1V")
            tempo$value    <- tempo$value - 1 
            summary    <- rbind(summary, tempo)
            df_trait_t <- estimateSquare(trait_data_full_i = summary, var=paste0("minus1", v))
            df_trait_t$variable <- v
            summary    <- rbind(summary, df_trait_t)
            summary    <- summary[-grep(summary$variable, pattern= "minus"),]
          } else {
            df_trait_t <- estimateSquare(trait_data_full_i = summary, var=v)
            summary    <- rbind(summary, df_trait_t)
          }
        }
      }
    }      
    if (any(grepl(traits_selected, pattern="age1"))) {
      var=traits_to_plot2[[traits_selected]][grep(traits_to_plot2[[traits_selected]], pattern="age1")]
      for (v in var[ifelse(length(var)>1, -length(var), 1)]){
        if(!(v %in% summary$variable) & length(str_subset(pattern = c("[x]"), string =  v, negate=T))>0) {
          df_trait_t <- estimateage1(trait_data_full_i = summary, var=v)
          summary    <- rbind(summary, df_trait_t)
        }
      }
    }
    if (any(grepl(traits_selected, pattern="x"))) {
      var=traits_to_plot2[[traits_selected]][grep(traits_to_plot2[[traits_selected]], pattern="x")]
      for (v in var[ifelse(length(var)>1,-c(length(var)), 1)]){
        if(!(v %in% summary$variable) & length(str_subset(pattern = c("\n"), string =  v, negate=T))>0) {
          df_trait_t <- estimateX(trait_data_full_i = summary, var=v)
          summary    <- rbind(summary, df_trait_t)
        }
      }
    }
    if (any(grepl(traits_selected, pattern="sum"))) {
      var=traits_to_plot2[[traits_selected]][grep(traits_to_plot2[[traits_selected]], pattern="sum")]
      for (v in var[ifelse(length(var)>1, -length(var), 1)]){
        if(!(v %in% summary$variable) & (length(str_subset(pattern = c("Div"), string =  v, negate=T))>0)) {
          df_trait_t <- estimateSum(trait_data_full_i = summary, var=v)
          summary    <- rbind(summary, df_trait_t)
        }
      }
    }
    if (any(grepl(traits_selected, pattern="Div"))) {
      var=traits_to_plot2[[traits_selected]][grep(traits_to_plot2[[traits_selected]], pattern="Div")]
      for (v in var[ifelse(length(var)>1, -length(var), 1)]){
        if(!(v %in% summary$variable) & (length(str_subset(pattern = c("colSum"), string =  v, negate=T))>0)) {
          df_trait_t <- estimateDiv(trait_data_full_i = summary, var=v)
          summary    <- rbind(summary, df_trait_t)
        }
      }
    }
    if (any(grepl(traits_selected, pattern="CS"))) {
      var=traits_to_plot[[traits_selected]][grep(traits_to_plot2[[traits_selected]], pattern="colSum")]
      df_trait_t <- estimatecolSum(trait_data_full_i = summary, var=var)
      summary    <- rbind(summary, df_trait_t)
      if (var == "colSumVsquarexDead\nsumVxDead\nDivage1abundanceBy"){
        df_trait_t <- estimateDiv(trait_data_full_i = summary, var=paste0(var, "\nDivage1abundanceBy"))
        df_trait_t$variable <- "colzero"
        df_trait_t$value    <- 0
        summary    <- rbind(summary, df_trait_t)
        traits_to_plot2[[traits_selected]] <- c(traits_to_plot2[[traits_selected]], "colzero")
      }
    }
    df_trait <- summary[as.vector(unlist(sapply(traits_to_plot2[[traits_selected]], FUN = function(x) which(summary$variable == x)))),]
    if (nrow(df_trait)<1) next
    var_i    <- unique(df_trait$variable)
    
    # compare all traits to their value /N1 to understand where is the increase coming from 
    for (v in var_i[-c(grep(var_i, pattern="age1"), grep(var_i, pattern="colzero"))]){
      df_trait_t <- estimateDiv(trait_data_full_i = df_trait, var=paste0(v, "\nDivage1abundanceBy"))
      df_trait   <- rbind(df_trait, df_trait_t)
    }
    var_i    <- unique(df_trait$variable)
    # if (any(grepl(var_i, pattern=c("of")))) var_i <- var_i[-grep(var_i, pattern=c("of"))]
    # if (any(grepl(var_i, pattern=c("egg")))) var_i <- var_i[-grep(var_i, pattern=c("egg"))]
    # id_mean <- (grep(var_i, pattern="mean"))
    # which_to_remove <- NULL
    # if (any(grepl(var_i[id_mean], pattern=c("Numerator"))) & traits_selected != "N&T/Vk") which_to_remove <- grep(var_i[id_mean], pattern=c("Numerator"))
    # if (any(grepl(var_i[id_mean], pattern=c("NT"))) & traits_selected != "N&T/Vk") which_to_remove <- grep(var_i[id_mean], pattern=c("NT"))
    # if (length(which_to_remove)>0) var_i <- var_i[-id_mean[which_to_remove]]
    # if (traits_selected == "N&T/Vk" & length(var_i)>2) var_i <- c(var_i[-c(grep(pattern="Nb", var_i), grep(pattern="Ne", var_i))])
    
    for (sp in unique(df_trait$Species)[1:7]){
      cat(sp)
      df_trait_sp <- df_trait %>% filter(Species == sp)
      df_trait_sp <- df_trait_sp %>% filter(Fstart == 1)
      df_trait_sp <- df_trait_sp %>% filter(variable %in% var_i)
      
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
          p[[sp]][[traits_selected]][[Fmax_i]] <- NULL
          next
        }
        
        # Use the function
        if (length(var_i)>3){
          k=1
          id_to_plot = NULL
          for (i in 1:c((length(var_i)+1)/2)){
            id_to_plot = c(k, ifelse((k+1)>length(var_i), 1, k+1))
            if (traits_selected == "N&T/Vk") id_to_plot=c(which(var_i == "T/Vk"),
                                                          ifelse(i==1, which(var_i != "T/Vk")[1], which(var_i != "T/Vk")[2]))
            data_split <- df_trait_sp_a[which(df_trait_sp_a$variable %in% var_i[id_to_plot] ),]
            ref_data_split <- ref_data[which(ref_data$variable %in% var_i[id_to_plot]),]
            if (all(is.na(data_split$value))) next
            result <- plot_trait_vs_fmax(df_trait_sp = data_split, trait_name = var_i[id_to_plot],
                                         Time1, Time2, relative_to_ref = relative_to_ref, spec=sp, 
                                         age_dim=T, reference_input=ref_data)
            p[[sp]][[paste(var_i[id_to_plot], collapse="")]][[Fmax_i]]         <- result$plot
            scale_v[[sp]][[paste(var_i[id_to_plot], collapse="")]][[Fmax_i]]   <- result$scale
            k=k+2
          }
        } else {
          result <- plot_trait_vs_fmax(df_trait_sp = df_trait_sp_a, trait_name = var_i,
                                       Time1, Time2, relative_to_ref = relative_to_ref,
                                       spec=sp, age_dim=T, reference_input=ref_data)
          p[[sp]][[traits_selected]][[Fmax_i]]         <- result$plot
          scale_v[[sp]][[traits_selected]][[Fmax_i]]   <- result$scale
        }
      }#age
    }#sp
  } # trait selected in the list
  plot_biom = ifelse(relative_to_ref, "rel_biomassplot", "biomassplot")
  plot_biom = ifelse(relative_to_ref, "biomassplot", "biomassplot")
  suffix_relative <- ifelse(relative_to_ref, "_RELATIVE", "")
  for (sp in names(p)){
    if (!dir.exists(paste0(PATHOUT, "/byage", sp))) dir.create(paste0(PATHOUT, "/byage", sp))
    pdf(file = paste0(PATHOUT, "/byage", sp, "/VxEMERGEDplot_", "data_Ne", sp, "vsFmax", suffix_relative, ".pdf"), 
        width = 15, height = 8)
    if (is.null(p[[sp]])) next
    for (traits_selected in names(p[[sp]])) {
      plot_list <- c(p[[sp]][[traits_selected]], list(get(plot_biom)[[sp]]))
      plot_list <- lapply(plot_list, function(x) x + theme_minimal(base_size = 11) + theme(legend.position = "bottom"))
      combined_plot <- ggarrange(plotlist = plot_list, ncol = 3,
                                 nrow = ceiling((length(plot_list))/3),
                                 common.legend = TRUE, legend = "bottom")
      print(combined_plot)
    }
    dev.off()
  }
}# relative


a <- only_Ne_traits[which(only_Ne_traits$Species == "Hake"),]
a <- a[which(a$Fmax == 1),]
a <- a[which(a$Fstart == 1),]
a <- a[-which(a$age == 1),]

boxplot(a$value~a$age)





# PCA --- -------------------------------------------------------

## PLOT raw traits -------------------
plot_PCA(data=dataRATIO_mean_sc, Traits_raw = T, pathoutputplots=paste0(PATHOUT, "/PCA/raw_traits"))
plot_PCA_sep_species(data=dataRATIO_mean_sc, Traits_raw = T, pathoutputplots=paste0(PATHOUT, "/PCA/raw_traits"))
plot_PCA_ALL(data=dataRATIO_mean_sc, Traits_raw = T, pathoutputplots=paste0(PATHOUT, "/PCA/raw_traits"))

## PLOT calculated traits ------------------- OK
plot_PCA(data=dataRATIO_mean_sc, Traits_Ne = T, pathoutputplots=paste0(PATHOUT, "/PCA/Ne_traits"))
plot_PCA_sep_species(data=dataRATIO_mean_sc, Traits_Ne = T, pathoutputplots=paste0(PATHOUT, "/PCA/Ne_traits"))
plot_PCA_ALL(data=dataRATIO_mean_sc, Traits_Ne = T, pathoutputplots=paste0(PATHOUT, "/PCA/Ne_traits"))

## PLOT raw + calculated  -------------------
plot_PCA(data=dataRATIO_mean_sc, pathoutputplots=paste0(PATHOUT, "/PCA/alltraits"))
plot_PCA_sep_species(data=dataRATIO_mean_sc, pathoutputplots=paste0(PATHOUT, "/PCA/alltraits"))
plot_PCA_ALL(data=dataRATIO_mean_sc, pathoutputplots=paste0(PATHOUT, "/PCA/alltraits"))



# MORTALITY comparison --------------------------------------------

data_for_mort_indic <- prepare_only_simutraits_dataset(path_simutraits="osmose-ns-v4.4.1_clim/output/PLOTS/Trait_variables_by_year_by_age/", natural_mort =T)
age_mat <- prepare_only_simutraits_dataset(path_simutraits="osmose-ns-v4.4.1_clim/output/PLOTS/Trait_variables_by_year_BYSP/", natural_mort =T)
age_mat <- age_mat[which(age_mat$TRAIT=="age_Mature"),]






