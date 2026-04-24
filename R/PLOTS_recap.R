# the idea is to store separately all the key parameters per species into plots
# and to produce one pdf sheet per species with all the key variables to have it all summarized at one glance




#### Initialisation --------------
library(legendry)
library(ggpubr)

# functions and environment 
main = "New_gpu_conf.R" # config file
source("osmose-ns-v4.4.1_clim/R/functions/source.R") # load config and key parameters and functions
source("osmose-ns-v4.4.1_clim/R/functions/functions_key_traits.R") # plots against Fmax
source("osmose-ns-v4.4.1_clim/R/functions/auxiliary_for_lists.R") # from list to df 
source("osmose-ns-v4.4.1_clim/R/functions/functions_output_SimuRData_into_df_for_plots.R") # to transform dataset from RData to long df
source("osmose-ns-v4.4.1_clim/R/functions/function_data_traits_summary.R")

# paths
PATHOUT  = "osmose-ns-v4.4.1_clim/output/PLOTS/KEY_TRAITS/"
if (!dir.exists(PATHOUT)) dir.create(PATHOUT, recursive = TRUE)
base_dir <- "osmose-ns-v4.4.1_clim/output/PLOTS/FINAL_Ne_Htz2/INIT"
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
dataRATIO_full_merged <- dataset_attach_all_traits_genetic_drift(pathoutput=pathoutputRDATA, list_outputs=list("ratio"=list_data_traits_ratio),
                                                                 trait_cv = list("folders_traits_c_a"=folders_traits_c_a, "folders_traits_y_a"=folders_traits_y_a),
                                                                 FORCE=F, keep_cohort_or_age_dimension = T)






#### PLOT 0- Biomass vs Fmax & keep it on the side -------------------
biomassplot = NULL
rel_biomassplot = NULL
# for (nm in names(dataRATIO_full_merged)[1]){
#   d <- dataRATIO_full_merged[[nm]]
#   d <- d[d$Fstart==1, ]
#   id_biom <- grep(pattern='biomass', colnames(d))
#   name_id <- colnames(d)[id_biom]
#   for (sp in species){
#     d_i_i <- d[which(d$Species == sp), ]
#     if (any(is.na(d_i_i[, name_id]))) d_i_i <- d_i_i[-which(is.na(d_i_i[, name_id])), ]
#     d_i_i$Iszero <- ifelse(d_i_i[,name_id] < 100, "Abundance <100", "Abundance >100")
# 
# 
#     biomassplot[[sp]] <- ggplot(d_i_i, aes(y=.data[[name_id]], x=Fmax)) +
#       geom_point(aes(colour = Iszero))+
#       labs(title = "Biomass per fishing pressure", subtitle = sp)
# 
#     ref_data <- d_i_i %>%
#       filter(Fmax == 1 & Fstart == 1) %>%
#       dplyr::select(.data[[name_id]], replicate) %>%
#       dplyr::rename(ref_value = .data[[name_id]])
#     # Subtract reference from all means
#     rel_d_i_i <- d_i_i %>%
#       left_join(ref_data, by = c("replicate"), relationship = "many-to-one") %>%
#       mutate(relative_biomass =(.data[[name_id]] - ref_value)/ref_value) %>%
#       dplyr::select(-ref_value)
#     rel_d_i_i$Iszero_relative <- ifelse(rel_d_i_i$relative_biomass < 0.001, "Relative value is 0", "Different from 0")
#     rel_biomassplot[[sp]] <- ggplot(rel_d_i_i, aes(y=relative_biomass, x=Fmax)) +
#       geom_point(aes(colour = Iszero_relative))+
#       labs(title = "Relative biomass per fishing pressure", subtitle = sp)
#   }
# }
# list_bio <- list("rel_biomassplot"=rel_biomassplot, "biomassplot"=biomassplot)
# save(list_bio, file=paste0(pathoutputRDATA, "/Biomass_plots.RData"))
load(file=paste0(pathoutputRDATA, "/Biomass_plots.RData"))
rel_biomassplot=list_bio$rel_biomassplot
biomassplot=list_bio$biomassplot



#### PLOT 1- NE/N VS FMAX + Biomass VS FMAX ------------------


traits_to_plot <- list("T&Vk$T/Vk." = c("generation", "mean_V", "T/Vk."), "Ne&Ne/N"=c("Ne", "Numerator", "/"))

data =  dataRATIO_full_merged[[1]]

Ne_plot_list = list()
# for (relative_to_ref in c(TRUE, FALSE, "proportional")){
#   proportional = FALSE
#   if (relative_to_ref == "proportional") {
#     relative_to_ref = FALSE
#     proportional = TRUE
#   }
#   relative_to_ref = as.logical(relative_to_ref)
#   if (is.data.frame(data)) data <- list("data_Ne"=data)
#   # for the plot
#   if (!any(grepl(names(data), pattern="Ne"))) {
#     tempo <- c(names(data[[1]])[1], names(data[[1]])[2], names(data[[1]])[3])
#     if (length(data)==3) data <- list(do.call(rbind, list(data[[1]][[names(data[[1]])[1]]], data[[2]][[names(data[[1]])[1]]])),
#                                       do.call(rbind, list(data[[1]][[names(data[[1]])[2]]], data[[2]][[names(data[[1]])[2]]])),
#                                       do.call(rbind, list(data[[1]][[names(data[[1]])[3]]], data[[2]][[names(data[[1]])[3]]])))
#     if (length(data)==1)data <- list(do.call(rbind, list(data[[1]][[names(data[[1]])[1]]], data[[2]][[names(data[[1]])[1]]])))
#     names(data) <- tempo[1:length(data)]
#   }
#   scale_v = list()
#   for (genetic_drift_data in (1:length(data))[1]){
#     for (traits_selected in names(traits_to_plot)) {
#       # names_substr <- substr(summary$variable, nchar(summary$variable) - 3 + 1, nchar(summary$variable))
# 
#       # select dataset
#       summary = data[[genetic_drift_data]]
#       if (names(data)[genetic_drift_data] == "tempoNe") {
#         if (ncol(data[["data_Nb"]])>6)  summary <- full_join(summary, data[["data_Nb"]][,c(1:4,grep(colnames(data[["data_Nb"]]), pattern="mean_V"))])
#         if (ncol(data[["data_Nb"]])==6) summary <- rbind(summary, data[["data_Nb"]][grep(data[["data_Nb"]]$variable, pattern="mean_V"),])
#       }
#       if (names(data)[genetic_drift_data] == "data_Nb") {
#         if (ncol(data[["data_Nb"]])>6)  summary <- full_join(summary, data[["tempoNe"]][,c(1:4,grep(colnames(data[["tempoNe"]]), pattern="generation"))])
#         if (ncol(data[["data_Nb"]])==6) summary <- rbind(summary, data[["tempoNe"]][grep(data[["tempoNe"]]$variable, pattern="generation"),])
#       }
#       summary <- as_long_if_needed(summary, columns_keep_pivot = colnames(summary)[-which(colnames(summary) %in% c("replicate", "Fmax", "Fstart", "Species", "age", "cohort"))])
#       if (any(grepl(traits_selected, pattern="T/Vk."))) {
#         df_trait_t <- estimateT_Vk(trait_data_full_i = summary)
#         summary <- rbind(summary, df_trait_t[which(df_trait_t$variable == "T/Vk."),])
#       }
# 
#       WEIGHT = summary[which(summary$variable == "mean_Numerator_N_1"),]
#       colnames(WEIGHT)[ncol(WEIGHT)] <- "weight"
# 
#       df_trait <- summary[as.vector(unlist(sapply(traits_to_plot[[traits_selected]], FUN = function(x) grep(summary$variable, pattern = x)))),]
#       if (nrow(df_trait)<1) next
#       if (length(table(df_trait$variable)) != length(traits_to_plot[[traits_selected]])) {
#         if (traits_selected == "Ne&Ne/N") df_trait <- df_trait[which(df_trait$variable %in% c("mean_Numerator_N_1", "Ne_N_1", "Ne_N_1/Numerator_N_1")),]
#       }
#       var_i    <- unique(df_trait$variable)
# 
#       # if (any(grepl(var_i, pattern=c("of")))) var_i <- var_i[-grep(var_i, pattern=c("of"))]
#       # if (any(grepl(var_i, pattern=c("egg")))) var_i <- var_i[-grep(var_i, pattern=c("egg"))]
#       # id_mean <- (grep(var_i, pattern="mean"))
#       # which_to_remove <- NULL
#       # if (any(grepl(var_i[id_mean], pattern=c("Numerator"))) & traits_selected != "N&T/Vk") which_to_remove <- grep(var_i[id_mean], pattern=c("Numerator"))
#       # if (any(grepl(var_i[id_mean], pattern=c("NT"))) & traits_selected != "N&T/Vk") which_to_remove <- grep(var_i[id_mean], pattern=c("NT"))
#       # if (length(which_to_remove)>0) var_i <- var_i[-id_mean[which_to_remove]]
#       # if (traits_selected == "N&T/Vk" & length(var_i)>2) var_i <- c(var_i[-c(grep(pattern="Nb", var_i), grep(pattern="Ne", var_i))])
# 
# 
# 
#       for (sp in unique(df_trait$Species)){
#         df_trait_sp <- df_trait %>% filter(Species == sp, Fstart == 1, variable %in% var_i)
#         if (sum(is.nan(df_trait_sp$value))>0) df_trait_sp <- df_trait_sp[-which(is.nan(df_trait_sp$value)), ]
#         if (sum(is.infinite(df_trait_sp$value))>0) df_trait_sp <- df_trait_sp[-which(is.infinite(df_trait_sp$value)), ]
# 
#         # conditions to stop
#         if (nrow(df_trait_sp)==0) next
#         if (all(is.na(df_trait_sp$value))) {
#           Ne_plot_list[[paste0(as.character(relative_to_ref), as.character(proportional))]][[sp]][[traits_selected]] <- NULL
#           next}
# 
#         # Use the function
#         if (length(var_i)>2){
#           k=1
#           id_to_plot = NULL
#           for (i in 1:c((length(var_i)+1)/3)){
#             id_to_plot = k:(k+2)
#             # if (traits_selected == "N&T/Vk") id_to_plot=c(which(var_i == "T/Vk"),
#             #                                               ifelse(i==1, which(var_i != "T/Vk")[1], which(var_i != "T/Vk")[2]))
#             data_split <- df_trait_sp[which(df_trait_sp$variable %in% var_i[id_to_plot] ),]
#             data_split$variable <- improve_vriable_names(data_split$variable)
#             data_split$variable <- improve_vriable_names(data_split$variable)
#             t <- improve_vriable_names(var_i[id_to_plot])
#             if (all(is.na(data_split$value))) next
#             result <- plot_trait_vs_fmax(df_trait_sp = data_split, trait_name = t, age_dim = F,weight = WEIGHT,
#                                          Time1, Time2, relative_to_ref = relative_to_ref, spec=sp, proportional=proportional,
#                                          position_legend = "right")
#             Ne_plot_list[[paste0(as.character(relative_to_ref), as.character(proportional))]][[sp]][[paste(t, collapse="")]]<- result$plot
#             scale_v[[sp]][[paste(t, collapse="")]] <- result$scale
#             k=k+2
#           }
#         } else {
#           df_trait_sp_i <- df_trait_sp
#           df_trait_sp_i$variable <- improve_vriable_names(df_trait_sp_i$variable)
#           t <- improve_vriable_names(var_i)
#           result <- plot_trait_vs_fmax(df_trait_sp = df_trait_sp, trait_name = t, age_dim = F, weight = WEIGHT,
#                                        Time1, Time2, relative_to_ref = relative_to_ref, spec=sp, proportional = proportional,
#                                        position_legend = "right")
#           Ne_plot_list[[paste0(as.character(relative_to_ref), as.character(proportional))]][[sp]][[t]]       <- result$plot
#           scale_v[[sp]][[t]] <- result$scale
#         }
#       }#sp
#     } # trait selected in the list
#     # trait_i = unique(var_i)
#     # if (length(unique(var_i))) trait_i = paste0(unique(var_i)[1], "_WITH_")
#     # plot_biom = ifelse(relative_to_ref, "rel_biomassplot", "biomassplot")
#     # plot_biom = ifelse(relative_to_ref, "biomassplot", "biomassplot")
#     # suffix_relative <- ifelse(relative_to_ref, "_RELATIVE", "")
#     # pdf(file = paste0(PATHOUT, "/plot_", names(data)[genetic_drift_data], "_T_Vk_N_", "vsFmax", suffix_relative, ".pdf"), width = 15, height = 10)
#     # par(mfrow=c(1,1))
#     # for (sp in species){
#     #   combined_plot <- ggarrange(plotlist = c(Ne_plot_list[[as.character(relative_to_ref)]][[sp]], list(get(plot_biom)[[sp]])), ncol = 2, nrow=ceiling((length(p[[sp]])+1)/2))
#     #   print(combined_plot)
#     # }
#     # dev.off()
#   } # genetic drift dim
# }# relative
# save(Ne_plot_list, file=paste0(pathoutputRDATA, "/Ne_plots.RData"))
load(file=paste0(pathoutputRDATA, "/Ne_plots.RData"))




#### PLOT 2- SSDg/N1 SSDi/N1 VS FMAX  ------------------


Time1=150
Time2=160

PATHOUT  = "osmose-ns-v4.4.1_clim/output/PLOTS/KEY_TRAITS/"
path <- c("C:/Users/mbeneat/Documents/osmose/Chap2/osmose-ns/osmose-ns-v4.4.1_clim/output/PLOTS/FINAL_Ne_Htz2")


### produce dataset reading multiple .RData to make 1 dataframe, stored as a .RData file ----------
sp_i = "Hake"
FORCE=F
if (!file.exists(paste0(path, "/Complete_INIT2", sp_i, ".RData")) || FORCE ==T){
  full_dataset = list()
  full_dataset_age = list()
  for (f_i in folders){
    cat(f_i, "\n")

    Fmax   <- stringr::str_split(stringr::str_split(f_i, pattern="Fmax")[[1]][2], pattern="_Fstart")[[1]][1]
    Fstart <- stringr::str_split(stringr::str_split(f_i, pattern="Fstart")[[1]][2], pattern="_sp")[[1]][1]

    if (Fstart != 1) next

    load(list.files(f_i, full.names = T))
    dataNe <- INIT2$data_Ne

    for (sp in names(dataNe)){

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

traits_sp_to_plot = c("colSumSSDix\nDivNumerator_N_1", "colSumSSDgx\nDivNumerator_N_1", "Vk_point",
                      "colSumDeltasquare\nDivNumerator_N_1", "DeadSum", "colSumDeltasquare\nDivNumerator_N_1\nxDead")
SSDplot = list()
# for (sp in species){
#   cat(sp)
# 
#   load(file=paste0(path, "/Complete_INIT2", sp, ".RData"))
#   final_dataset_age = complete[[2]]
#   final_dataset_age <- final_dataset_age %>%
#     dplyr::mutate(dplyr::across(c(cohort, replicate), ~ stringr::str_remove(.x, "^(cohort_|rep_)")))
#   final_dataset_age <- final_dataset_age[which(final_dataset_age$cohort %in% Time1:Time2),]
#   rm(complete)
# 
#   df_trait_sp <- final_dataset_age %>% filter(Species == sp, Fstart == 1)
# 
#   #add delta square
#   df_trait_t    <- estimateSquare(trait_data_full_i = df_trait_sp, var=paste0("Deltasquare"))
#   df_trait_sp   <- rbind(df_trait_sp, df_trait_t)
# 
#   #add sum by age of SSDgx and SSDix
#   df_trait_t    <- estimatecolSum(trait_data_full_i = df_trait_sp, var=paste0("colSumSSDix"))
#   df_trait_sp   <- rbind(df_trait_sp, df_trait_t)
#   df_trait_t    <- estimatecolSum(trait_data_full_i = df_trait_sp, var=paste0("colSumSSDgx"))
#   df_trait_sp   <- rbind(df_trait_sp, df_trait_t)
#   df_trait_t    <- estimatecolSum(trait_data_full_i = df_trait_sp, var=paste0("colSumDeltasquare"))
#   df_trait_sp   <- rbind(df_trait_sp, df_trait_t)
# 
#   var_i <- unique(df_trait_sp$variable)
#   # compare all traits to their value /N1 to understand where is the increase coming from
#   for (v in var_i){
#     if (v != "Numerator_N_1") {
#       df_trait_t    <- estimateDiv(trait_data_full_i = df_trait_sp, var=paste0(v, "\nDivNumerator_N_1"))
#       df_trait_sp   <- rbind(df_trait_sp, df_trait_t)
#     }
#   }
# 
#   #add the SSDgx estimation by hand
#   df_trait_t    <- estimateX(trait_data_full_i = df_trait_sp, var=paste0("Deltasquare\nDivNumerator_N_1\nxDead"))
#   df_trait_sp   <- rbind(df_trait_sp, df_trait_t)
#   df_trait_t    <- estimatecolSum(trait_data_full_i = df_trait_sp, var=paste0("colSumDeltasquare\nDivNumerator_N_1\nxDead"))
#   df_trait_sp   <- rbind(df_trait_sp, df_trait_t)
# 
#   # reduce dataset dimensions
#   var_i = traits_sp_to_plot
#   df_trait_sp  <- df_trait_sp[which(df_trait_sp$variable %in% c(var_i, "abundanceByAge")), ]
# 
#   # reduce dimension number by averaging (weighted mean)
#   weight_df <- build_weight_df(df = df_trait_sp, weight_var = "abundanceByAge",
#                                key_cols = c("replicate", "Species", "Fmax", "Fstart"))
#   df_trait_sp <- weighted_reduce(
#     df = df_trait_sp,
#     weight_df = weight_df,
#     join_cols = c("Species", "replicate", "Fmax", "Fstart"),
#     group_cols = c("Species", "replicate", "Fstart", "Fmax", "variable"),
#     drop_age_zero = TRUE
#   )
# 
#   for (relative_to_ref in c(FALSE, TRUE, "proportional")){
#     proportional = FALSE
#     if (relative_to_ref == "proportional") {
#       relative_to_ref = FALSE
#       proportional = TRUE
#     }
#     relative_to_ref = as.logical(relative_to_ref)
#     # for the plot
#     plot_biom = ifelse(relative_to_ref, "rel_biomassplot", "biomassplot")
#     plot_biom = ifelse(relative_to_ref, "biomassplot", "biomassplot")
# 
#     ref_data <- df_trait_sp %>%
#       filter(Fmax == 1 & Fstart == 1) %>%
#       dplyr::select(variable, value, replicate) %>%
#       dplyr::rename(ref_value = value)
# 
#     df_trait_sp <- clean_values(df_trait_sp)
# 
#     # conditions to stop
#     if (nrow(df_trait_sp)==0) next
#     if (all(is.na(df_trait_sp$value) | df_trait_sp$value == 0)) {
#       next
#     }
# 
#     # Use the function
#     k=1
#     id_to_plot = NULL
#     for (i in seq_len(2)){
#       if (i==1 ) id_to_plot = 1:3
#       if (i==2 ) id_to_plot = 4:5
# 
#       var_ii <- improve_vriable_names(var_i[id_to_plot])
#       # produce plot
#       data_split <- df_trait_sp[which(df_trait_sp$variable %in% var_i[id_to_plot] ),]
#       data_split$variable <- improve_vriable_names(data_split$variable)
#       # data_split <- data_split[, -which(colnames(data_split) == "Fstart")]
#       ref_data_split <- ref_data[which(ref_data$variable %in% var_i[id_to_plot]),]
#       ref_data_split$variable <- improve_vriable_names(ref_data_split$variable)
#       # ref_data_split <- ref_data_split[, -which(colnames(ref_data_split) == "Fstart")]
#       # p_l <- list(c(0.35, 0.35), c(0.35, 0.70))[[i]]
#       if (all(is.na(data_split$value))) next
#       result <- plot_trait_vs_fmax(df_trait_sp = data_split, trait_name = var_ii,
#                                    Time1, Time2, relative_to_ref = relative_to_ref, spec=sp,
#                                    age_dim=F, reference_input=ref_data_split, FORCElm=F,
#                                    proportional = proportional, weight = weight_df, position_legend="auto")
# 
#       SSDplot[[paste0(as.character(relative_to_ref), as.character(proportional))]][[sp]][[paste(var_ii, collapse="")]]   <- result$plot
#       k=k+2
#     }
#   } # relative to ref
#   # # save plot
#   # suffix_relative <- paste0(ifelse(relative_to_ref, "_RELATIVE", ""))
#   # if (!dir.exists(paste0(PATHOUT, "/byage", sp))) dir.create(paste0(PATHOUT, "/byage", sp))
#   # if (is.null(p[[sp]])) next
#   # pdf(file = paste0(PATHOUT, "/byage", sp, "/VxRAW2plot_", "data_Ne", sp, "vsFmax", Fmax_i, suffix_relative, ".pdf"),
#   #     width = 15, height = 8)
#   # for (trait_plotted in names(p[[sp]])){
#   #   plot_list <- c(SSDplot[[as.character(relative_to_ref)]][[sp]][[trait_plotted]], list(get(plot_biom)[[sp]]))
#   #   plot_list <- lapply(plot_list, function(x) x + theme_minimal(base_size = 11) + theme(legend.position = "bottom"))
#   #   combined_plot <- ggarrange(plotlist = plot_list, ncol = 3,
#   #                              nrow = ceiling((length(plot_list))/3),
#   #                              common.legend = TRUE, legend = "bottom")
#   #   print(combined_plot)
#   # }
#   # dev.off()
# }#sp
# save(SSDplot, file=paste0(pathoutputRDATA, "/SSDplots.RData"))
load(file=paste0(pathoutputRDATA, "/SSDplots.RData"))

#### PLOT 2bis: delta² and deaths by age ---------------
# 
# # age delta²/N and Death by age
combined_plot =list()
Time1=150
Time2=160
traits_sp_to_plot = list(c("Deltasquare","Dead"),c("Deltasquare\nDivNumerator_N_1","Dead\nDivNumerator_N_1"),
                      c("DeltasquarexDead", "DeltasquarexDead\nDivNumerator_N_1"), 
                      c("survival_rate"), c("Delta"))
PATHOUT  = "osmose-ns-v4.4.1_clim/output/PLOTS/KEY_TRAITS/"
path <- c("C:/Users/mbeneat/Documents/osmose/Chap2/osmose-ns/osmose-ns-v4.4.1_clim/output/PLOTS/FINAL_Ne_Htz2")
plots_by_age = list()
# for (sp in species){
#   cat(sp)
# 
#   load(file=paste0(path, "/Complete_INIT2", sp, ".RData"))
#   final_dataset_age = complete[[2]]
#   final_dataset_age <- final_dataset_age %>%
#     dplyr::mutate(dplyr::across(c(cohort, replicate), ~ stringr::str_remove(.x, "^(cohort_|rep_)")))
#   final_dataset_age <- final_dataset_age[which(final_dataset_age$cohort %in% Time1:Time2),]
#   rm(complete)
# 
#   # remove ages that are not included
#   final_dataset_age$age <- as.numeric(final_dataset_age$age)
#   ages_to_remove <- table(final_dataset_age[which(final_dataset_age$variable == "Dead" & final_dataset_age$value ==0),"age"])
#   ages_to_remove <- ages_to_remove[which(ages_to_remove == max(ages_to_remove))]
#   final_dataset_age <- final_dataset_age[-which(final_dataset_age$age %in% as.numeric(names(ages_to_remove))), ]
# 
#   df_trait_sp <- final_dataset_age %>% filter(Species == sp, Fstart == 1)
# 
#   #add delta square
#   df_trait_t    <- estimateSquare(trait_data_full_i = df_trait_sp, var=paste0("Deltasquare"))
#   df_trait_sp   <- rbind(df_trait_sp, df_trait_t)
# 
#   #add delta square
#   df_trait_t    <- estimateX(trait_data_full_i = df_trait_sp, var=paste0("DeltasquarexDead"))
#   df_trait_sp   <- rbind(df_trait_sp, df_trait_t)
# 
#   #add sum by age of SSDgx and SSDix
#   df_trait_t    <- estimatecolSum(trait_data_full_i = df_trait_sp, var=paste0("colSumSSDix"))
#   df_trait_sp   <- rbind(df_trait_sp, df_trait_t)
#   df_trait_t    <- estimatecolSum(trait_data_full_i = df_trait_sp, var=paste0("colSumSSDgx"))
#   df_trait_sp   <- rbind(df_trait_sp, df_trait_t)
#   df_trait_t    <- estimatecolSum(trait_data_full_i = df_trait_sp, var=paste0("colSumDeltasquare"))
#   df_trait_sp   <- rbind(df_trait_sp, df_trait_t)
# 
#   var_i <- unique(df_trait_sp$variable)
# 
#   # keep abundance by age variable on the side
#   WEIGHT = df_trait_sp[grep(df_trait_sp$variable, pattern="abundanceByAge"),]
# 
#   # compare all traits to their value /N1 to understand where is the increase coming from
#   df_trait_t    <- estimatecolSum(trait_data_full_i = df_trait_sp, var=paste0("colSumDeltasquarexDead"))
#   df_trait_sp   <- rbind(df_trait_sp, df_trait_t)
# 
#   for (v in var_i){
#     if (v != "Numerator_N_1") {
#       df_trait_t    <- estimateDiv(trait_data_full_i = df_trait_sp, var=paste0(v, "\nDivNumerator_N_1")) # div by population size
#       df_trait_sp   <- rbind(df_trait_sp, df_trait_t)
#     }
# 
#   }
# 
#   ssdg <-  add_SSDg(df=df_trait_sp, weight=WEIGHT)
#   df_trait_sp <- left_join(df_trait_sp, ssdg, by=c("Fmax", "Fstart", "variable", "Species"), relationship = "many-to-many")
#   
# 
#   var_i <- unique(df_trait_sp$variable)
# 
#   for (relative_to_ref in c(TRUE, FALSE, "proportional")){
#     proportional = FALSE
#     if (relative_to_ref == "proportional") {
#       relative_to_ref = FALSE
#       proportional = TRUE
#     }
#     relative_to_ref = as.logical(relative_to_ref)
#     # for the plot
#     plot_biom = ifelse(relative_to_ref, "rel_biomassplot", "biomassplot")
#     plot_biom = ifelse(relative_to_ref, "biomassplot", "biomassplot")
# 
#     ref_data <- df_trait_sp %>%
#       filter(Fstart == 1 & Fmax==1) %>%
#       dplyr::select(variable, value, replicate, age) %>%
#       dplyr::rename(ref_value = value)
#     df_trait_sp_a <- df_trait_sp
# 
#     df_trait_sp_a <- clean_values(df_trait_sp_a)
# 
#     # conditions to stop
#     if (nrow(df_trait_sp_a)==0) next
#     if (all(is.na(df_trait_sp_a$value) | df_trait_sp_a$value == 0)) {
#       next
#     }
# 
# 
#     # produce plot
#     var_i <- unlist(traits_sp_to_plot)
#     df_trait_sp_trait <- df_trait_sp_a[which(df_trait_sp_a$variable %in% var_i),]
#     ref_data_split <- ref_data[which(ref_data$variable %in% var_i),]
#     if (all(is.na(df_trait_sp_trait$value))) next
#     for (t in traits_sp_to_plot){
#       df_trait_sp_plot <- df_trait_sp_trait[which(df_trait_sp_trait$variable == t),]
#       df_trait_sp_plot$variable <- improve_vriable_names(df_trait_sp_plot$variable)
#       t <- improve_vriable_names(t)
#       df_trait_sp_plot$variable <- as.factor(df_trait_sp_plot$variable)
#       factor_order_forced = factor(df_trait_sp_plot$variable,levels=c(levels(df_trait_sp_plot$variable), c('', ' ', '  ')[1:((2*2)-length(levels(df_trait_sp_plot$variable)))]))
#       df_trait_sp_plot$variable <- factor_order_forced
#       result <- plot_trait_vs_fmax(df_trait_sp=df_trait_sp_plot, spec = sp, trait_name = t, age_dim = T,
#                                          Time1=150, Time2=160, relative_to_ref = relative_to_ref,
#                                          #scale_factor = 1, #
#                                          ncol=2, nrow=2, 
#                                          weight=WEIGHT, facet=F, 
#                                          proportional = proportional, position_legend="auto")
#       plots_by_age[[paste0(as.character(relative_to_ref), as.character(proportional))]][[sp]][[paste(t, collapse="_")]]<- result
#     }
#   }#relative/proportional/none
# }#sp
# save(plots_by_age, file=paste0(pathoutputRDATA, "/deaths_deltasquare_byage.RData"))
deaths_deltasquare_byage <- get(load(file=paste0(pathoutputRDATA, "/deaths_deltasquare_byage.RData")))



#### PLOT 3: mortality decomposition --------------------------

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
FORCE=F
if (!file.exists(save_file) || FORCE==T){
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
data_complete_c[which((data_complete_c$TRAIT) == "mean_of_Mnat_per_age"), "TRAIT"] <- "mean_of_Madd_per_age"
data_complete_c <- pivot_wider(data = data_complete_c, names_from = "TRAIT", values_from = "value")
data_complete_c <- data_complete_c[,-grep(colnames(data_complete_c), pattern="mean_of_Dis_per_age")]
data_complete_c <- data_complete_c[,-grep(colnames(data_complete_c), pattern="mean_of_Z_per_age")]
data_complete_c$mean_of_Mnat_per_age <- data_complete_c$mean_of_Madd_per_age + data_complete_c$mean_of_Mpred_per_age + data_complete_c$mean_of_Mstar_per_age
data_complete_c <- pivot_longer(data = data_complete_c, names_to = "TRAIT", values_to = "value", cols=-colnames(data_complete_c)[1:4])

# changing colnames to match expectations
colnames(data_complete_c)[which(colnames(data_complete_c) == "TRAIT")] <- "variable"
data_complete_c$Fmax   <- as.numeric(unlist(stringr::str_remove(unlist(stringr::str_split(data_complete_c$Fish_press, " Fstart"))[seq(1, 2*length(data_complete_c$Fish_press), by=2)], "Fmax=")))
data_complete_c$Fstart <- as.numeric(unlist(stringr::str_split(data_complete_c$Fish_press, " Fstart="))[seq(2, 2*length(data_complete_c$Fish_press), by=2)])
data_complete_c <- (data_complete_c[,-which(colnames(data_complete_c) == "Fish_press")])



PATHOUT  = "osmose-ns-v4.4.1_clim/output/PLOTS/KEY_TRAITS/"
mortality_plot = list()
# for (sp in species){
#   cat(sp)
# 
#   df_trait_sp <- data_complete_c %>% filter(Species == sp, Fstart == 1)
# 
#   df_trait_t    <- estimateSum(trait_data_full_i = df_trait_sp, var=paste0("mean_of_F_per_age", "\nsum", "mean_of_Mnat_per_age"))
#   df_trait_t$variable    <- "sum_mortality_rates"
#   df_trait_sp   <- rbind(df_trait_sp, df_trait_t)
# 
#   var_i <- unique(df_trait_sp$variable)
#   # compare all mortality rates to their value *abundance to have the number of individuals to understand where is the increase coming from
#   for (v in var_i){
#     if (any(grepl(pattern="mean_of_M", v) | grepl(pattern="mean_of_F", v) | grepl(pattern="sum_mortality_rates", v))) {
#       df_trait_t    <- estimateMnumber(trait_data_full_i = df_trait_sp, var=paste0(v))
#       df_trait_sp   <- rbind(df_trait_sp, df_trait_t)
#     }
#   }
#   var_i <- unique(df_trait_sp$variable)
#   for (v in var_i){
#     if (any(grepl(pattern="mean_of_M", v) | grepl(pattern="mean_of_F", v))) {
#       df_trait_t    <- estimateDiv(trait_data_full_i = df_trait_sp, var=paste0(v, "\nDiv", "sum_mortality_rates"))
#       df_trait_sp   <- rbind(df_trait_sp, df_trait_t)
#     }
#   }
#   df_trait_t    <- estimatecolSum(trait_data_full_i = df_trait_sp, var=paste0("colSumall_deaths"))
#   df_trait_sp   <- rbind(df_trait_sp, df_trait_t)
#   var_i <- unique(df_trait_sp$variable)
#   for (v in var_i){
#     if (any(grepl(pattern="deaths", v) & c(grepl(pattern="M", v)) | (grepl(pattern="F_deaths", v)))) {
#       df_trait_t    <- estimateDiv(trait_data_full_i = df_trait_sp, var=paste0(v, "\nDiv", "colSumall_deaths"))
#       df_trait_t2   <- estimatecolSum(trait_data_full_i = df_trait_t, var=paste0("colSum", (df_trait_t$variable)[1]))
#       df_trait_sp   <- rbind(df_trait_sp, df_trait_t)
#       df_trait_sp   <- rbind(df_trait_sp, df_trait_t2)
#     }
#   }
#   for (v in var_i){
#     if (v == "all_deaths") {
#       df_trait_t2   <- estimatecolSum(trait_data_full_i = df_trait_t, var=paste0("colSum", v))
#       df_trait_sp   <- rbind(df_trait_sp, df_trait_t2)
#     }
#   }
#   var_i <- unique(df_trait_sp$variable)
#   df_trait_sp_a <- df_trait_sp
#   # reduce dimension number
#   weight_df <- build_weight_df(df = df_trait_sp_a, weight_var = "mean_of_abundanceByAge_per_age",
#                                key_cols = c("replicate", "Species", "Fmax", "Fstart"))
#   df_trait_sp_a <- weighted_reduce(
#     df = df_trait_sp_a,
#     weight_df = weight_df,
#     join_cols = c("Species", "replicate", "Fmax", "Fstart"),
#     group_cols = c("Species", "replicate", "Fstart", "Fmax", "variable"),
#     drop_age_zero = TRUE
#   )
# 
#   # reference dataset
#   ref_data <- df_trait_sp_a %>%
#     filter(Fmax == 1 & Fstart == 1) %>%
#     dplyr::select(variable, value, replicate, Species) %>%
#     dplyr::rename(ref_value = value)
# 
# 
#   df_trait_sp_a <- clean_values(df_trait_sp_a)
# 
#   # conditions to stop
#   if (nrow(df_trait_sp_a)==0) next
#   if (all(is.na(df_trait_sp_a$value) | df_trait_sp_a$value == 0)) {
#     next
#   }
#   for (relative_to_ref in c(TRUE, FALSE, "proportional")){
#     proportional = FALSE
#     if (relative_to_ref == "proportional") {
#       relative_to_ref = FALSE
#       proportional = TRUE
#     }
#     relative_to_ref <- as.logical(relative_to_ref)
# 
#     # for the plot
#     plot_biom = ifelse(relative_to_ref, "rel_biomassplot", "biomassplot")
#     plot_biom = ifelse(relative_to_ref, "biomassplot", "biomassplot")
# 
#     result =NULL
#     var_i = c("sum_mortality_rates", "colSumall_deaths")
#     data_split <- df_trait_sp_a[which(df_trait_sp_a$variable %in% var_i ),]
#     ref_data_split <- ref_data[which(ref_data$variable %in% var_i),]
#     data_split$variable <- improve_vriable_names(data_split$variable)
#     var_i <- improve_vriable_names(var_i)
#     ref_data_split$variable <- improve_vriable_names(ref_data_split$variable)
#     if (all(is.na(data_split$value))) next
#     result <- plot_trait_vs_fmax(df_trait_sp = data_split, trait_name = var_i,
#                                  Time1=150, Time2=160, relative_to_ref = relative_to_ref, spec=sp,
#                                  age_dim=F, reference_input=ref_data_split, proportional = proportional,
#                                  weight=weight_df,
#                                  chosen_title=paste0(ifelse(proportional, "Proportional ", ""),
#                                                      ifelse(relative_to_ref, "Relative ", ""),
#                                                      "Total Mortality rates and Deaths numbers\naverage years ", Time1, " to ",Time2))
#     mortality_plot[[paste0(as.character(relative_to_ref), as.character(proportional))]][[sp]][[stringr::str_remove(paste(var_i, collapse=""), pattern="\n")]]    <- result$plot
# 
#     result =NULL
#     var_i = c("mean_of_Mstar_per_age\nDivsum_mortality_rates", "mean_of_Mpred_per_age\nDivsum_mortality_rates",
#               "mean_of_F_per_age\nDivsum_mortality_rates", "mean_of_Madd_per_age\nDivsum_mortality_rates")
#     data_split <- df_trait_sp_a[which(df_trait_sp_a$variable %in% var_i ),]
#     ref_data_split <- ref_data[which(ref_data$variable %in% var_i),]
#     data_split$variable <- improve_vriable_names(data_split$variable)
#     var_i <- improve_vriable_names(var_i)
#     ref_data_split$variable <- improve_vriable_names(ref_data_split$variable)
#     if (all(is.na(data_split$value))) next
#     result <- plot_trait_vs_fmax(df_trait_sp = data_split, trait_name = var_i,
#                                  Time1=150, Time2=160, relative_to_ref = relative_to_ref, spec=sp,
#                                  age_dim=F, reference_input=ref_data_split, proportional = proportional,
#                                  weight=weight_df,
#                                  chosen_title=paste0(ifelse(proportional, "Proportional ", ""),
#                                                      ifelse(relative_to_ref, "Relative ", ""),
#                                                     "Mortality rates\naverage years ", Time1, " to ",Time2))
#     mortality_plot[[paste0(as.character(relative_to_ref), as.character(proportional))]][[sp]][[stringr::str_remove(paste(var_i, collapse=""), pattern="\n")]]  <- result$plot
# 
#     result =NULL
#     var_i = c("colSumMstar_deaths\nDivcolSumall_deaths", "colSumF_deaths\nDivcolSumall_deaths",
#               "colSumMpred_deaths\nDivcolSumall_deaths", "colSumMadd_deaths\nDivcolSumall_deaths")
#     data_split <- df_trait_sp_a[which(df_trait_sp_a$variable %in% var_i ),]
#     ref_data_split <- ref_data[which(ref_data$variable %in% var_i),]
#     data_split$variable <- improve_vriable_names(data_split$variable)
#     var_i <- improve_vriable_names(var_i)
#     ref_data_split$variable <- improve_vriable_names(ref_data_split$variable)
#     if (all(is.na(data_split$value))) next
#     result <- plot_trait_vs_fmax(df_trait_sp = data_split, trait_name = var_i,
#                                  Time1=150, Time2=160, relative_to_ref = relative_to_ref, spec=sp,
#                                  age_dim=F, reference_input=ref_data_split, proportional = proportional,
#                                  weight=weight_df,
#                                  chosen_title=paste0(ifelse(proportional, "Proportional ", ""),
#                                                      ifelse(relative_to_ref, "Relative ", ""),
#                                                     "Death number per Mortality rate\naverage years ", Time1, " to ",Time2))
#     mortality_plot[[paste0(as.character(relative_to_ref), as.character(proportional))]][[sp]][[stringr::str_remove(paste(var_i, collapse=""), pattern="\n")]]  <- result$plot
# 
#     # # save plot
#     # suffix_relative <- paste0(ifelse(relative_to_ref, "_RELATIVE", ""))
#     # if (!dir.exists(paste0(PATHOUT, "/byage", sp))) dir.create(paste0(PATHOUT, "/byage", sp))
#     # if (is.null(p[[sp]])) next
#     # pdf(file = paste0(PATHOUT, "/byage", sp, "/VxMortal_", "data_Ne", sp, "vsFmax", Fmax_i, suffix_relative, ".pdf"),
#     #     width = 15, height = 8)
#     # for (trait_plotted in names(p[[sp]])){
#     #   plot_list <- c(p[[sp]][[trait_plotted]], list(get(plot_biom)[[sp]]))
#     #   plot_list <- lapply(plot_list, function(x) x + theme_minimal(base_size = 11) + theme(legend.position = "bottom"))
#     #   combined_plot <- ggarrange(plotlist = plot_list, ncol = 3,
#     #                              nrow = ceiling((length(plot_list))/3),
#     #                              common.legend = TRUE, legend = "bottom")
#     #   print(combined_plot)
#     # }
#     # dev.off()
#   }#sp
# }
# save(mortality_plot, file=paste0(pathoutputRDATA, "/mortalityplots.RData"))
load(file=paste0(pathoutputRDATA, "/mortalityplots.RData"))
 


# mortality & age 
PATHOUT  = "osmose-ns-v4.4.1_clim/output/PLOTS/KEY_TRAITS/"
mortality_plot_age = list()
# for (sp in species){
#   cat(sp)
# 
#   df_trait_sp <- data_complete_c %>% filter(Species == sp, Fstart == 1)
# 
#   df_trait_t    <- estimateSum(trait_data_full_i = df_trait_sp, var=paste0("mean_of_F_per_age", "\nsum", "mean_of_Mnat_per_age"))
#   df_trait_t$variable    <- "sum_mortality_rates"
#   df_trait_sp   <- rbind(df_trait_sp, df_trait_t)
# 
#   var_i <- unique(df_trait_sp$variable)
#   # compare all mortality rates to their value *abundance to have the number of individuals to understand where is the increase coming from
#   for (v in var_i){
#     if (any(grepl(pattern="mean_of_M", v) | grepl(pattern="mean_of_F", v) | grepl(pattern="sum_mortality_rates", v))) {
#       df_trait_t    <- estimateMnumber(trait_data_full_i = df_trait_sp, var=paste0(v))
#       df_trait_sp   <- rbind(df_trait_sp, df_trait_t)
#     }
#   }
#   var_i <- unique(df_trait_sp$variable)
#   for (v in var_i){
#     if (any(grepl(pattern="mean_of_M", v) | grepl(pattern="mean_of_F", v))) {
#       # cat(v)
#       df_trait_t    <- estimateDiv(trait_data_full_i = df_trait_sp, var=paste0(v, "\nDiv", "sum_mortality_rates"))
#       df_trait_sp   <- rbind(df_trait_sp, df_trait_t)
#     }
#   }
#   df_trait_t    <- estimatecolSum(trait_data_full_i = df_trait_sp, var=paste0("colSumall_deaths"))
#   df_trait_sp   <- rbind(df_trait_sp, df_trait_t)
#   var_i <- unique(df_trait_sp$variable)
#   for (v in var_i){
#     if (any(grepl(pattern="deaths", v) & c(grepl(pattern="M", v)) | (grepl(pattern="F_deaths", v)) | (v=="all_deaths"))) {
#       df_trait_t    <- estimateDiv(trait_data_full_i = df_trait_sp, var=paste0(v, "\nDiv", "colSumall_deaths"))
#       df_trait_t2   <- estimatecolSum(trait_data_full_i = df_trait_t, var=paste0("colSum", (df_trait_t$variable)[1]))
#       df_trait_sp   <- rbind(df_trait_sp, df_trait_t)
#       df_trait_sp   <- rbind(df_trait_sp, df_trait_t2)
#     }
#   }
#   df_trait_t    <- estimateage1(trait_data_full_i = df_trait_sp, var=paste0("age1mean_of_abundanceByAge_per_age"))
#   df_trait_sp   <- rbind(df_trait_sp, df_trait_t)
#   df_trait_t    <- estimateDiv(trait_data_full_i = df_trait_sp, var=paste0("all_deaths\nDiv", "age1mean_of_abundanceByAge_per_age"))
#   df_trait_sp   <- rbind(df_trait_sp, df_trait_t)
#   for (v in var_i){
#     if (v == "all_deaths") {
#       df_trait_t2   <- estimatecolSum(trait_data_full_i = df_trait_t, var=paste0("colSum", v))
#       df_trait_sp   <- rbind(df_trait_sp, df_trait_t2)
#     }
#   }
#   var_i <- unique(df_trait_sp$variable)
#   df_trait_sp_a <- df_trait_sp
#   # reduce dimension number
#   weight_df <- build_weight_df(df = df_trait_sp, weight_var = "mean_of_abundanceByAge_per_age",
#                                key_cols = c("replicate", "Species", "Fmax", "Fstart", "age"))
#   df_trait_sp <- weighted_reduce(
#     df = df_trait_sp,
#     weight_df = weight_df,
#     join_cols = c("Species", "replicate", "Fmax", "Fstart", "age"),
#     group_cols = c("Species", "replicate", "Fstart", "Fmax", "variable", "age"),
#     drop_age_zero = TRUE
#   )
# 
#   # reference dataset
#   ref_data <- df_trait_sp_a %>%
#     filter(Fmax == 1 & Fstart == 1) %>%
#     dplyr::select(variable, value, replicate, age, Species) %>%
#     dplyr::rename(ref_value = value)
# 
# 
#   df_trait_sp_a <- clean_values(df_trait_sp_a)
# 
#   # conditions to stop
#   if (nrow(df_trait_sp_a)==0) next
#   if (all(is.na(df_trait_sp_a$value) | df_trait_sp_a$value == 0)) {
#     next
#   }
# 
# 
#   for (relative_to_ref in c(TRUE, FALSE, "proportional")){
#     proportional = FALSE
#     if (relative_to_ref == "proportional") {
#       relative_to_ref = FALSE
#       proportional = TRUE
#     }
#     relative_to_ref <- as.logical(relative_to_ref)
# 
#     # for the plot
#     plot_biom = ifelse(relative_to_ref, "rel_biomassplot", "biomassplot")
#     plot_biom = ifelse(relative_to_ref, "biomassplot", "biomassplot")
# 
# 
#     result =NULL
#     var_i = c("all_deaths\nDivcolSumall_deaths", "sum_mortality_rates")
#     data_split <- df_trait_sp_a[which(df_trait_sp_a$variable %in% var_i ),]
#     ref_data_split <- ref_data[which(ref_data$variable %in% var_i),]
#     data_split$variable <- improve_vriable_names(data_split$variable)
#     var_i <- improve_vriable_names(var_i)
#     ref_data_split$variable <- improve_vriable_names(ref_data_split$variable)
#     if (all(is.na(data_split$value))) next
#     result <- plot_trait_vs_fmax(df_trait_sp = data_split, trait_name = var_i,
#                                  Time1=150, Time2=160, relative_to_ref = relative_to_ref, spec=sp,
#                                  age_dim=T, reference_input=ref_data_split, proportional = proportional,
#                                  weight=weight_df, 
#                                  chosen_title=paste0(ifelse(proportional, "Proportional ", ""),
#                                                      ifelse(relative_to_ref, "Relative ", ""),
#                                                      "Total Mortality rates and Deaths numbers\naverage years ", Time1, " to ",Time2))
#     mortality_plot_age[[paste0(as.character(relative_to_ref), as.character(proportional))]][[sp]][[stringr::str_remove(paste(var_i, collapse=""), pattern="\n")]]    <- result$plot
# 
#     # 
#     # result =NULL
#     # var_i = c("sum_mortality_rates", "all_deaths")
#     # data_split <- df_trait_sp_a[which(df_trait_sp_a$variable %in% var_i ),]
#     # ref_data_split <- ref_data[which(ref_data$variable %in% var_i),]
#     # data_split$variable <- improve_vriable_names(data_split$variable)
#     # var_i <- improve_vriable_names(var_i)
#     # ref_data_split$variable <- improve_vriable_names(ref_data_split$variable)
#     # if (all(is.na(data_split$value))) next
#     # result <- plot_trait_vs_fmax(df_trait_sp = data_split, trait_name = var_i,
#     #                              Time1=150, Time2=160, relative_to_ref = relative_to_ref, spec=sp,
#     #                              age_dim=T, reference_input=ref_data_split, proportional = proportional,
#     #                              weight=weight_df,
#     #                              chosen_title=paste0(ifelse(proportional, "Proportional ", ""),
#     #                                                  ifelse(relative_to_ref, "Relative ", ""),
#     #                                                  "Total Mortality rates and Deaths numbers\naverage years ", Time1, " to ",Time2))
#     # mortality_plot_age[[paste0(as.character(relative_to_ref), as.character(proportional))]][[sp]][[stringr::str_remove(paste(var_i, collapse=""), pattern="\n")]]    <- result$plot
#     # 
#     # 
#     result =NULL
#     var_i = c("mean_of_Mstar_per_age\nDivsum_mortality_rates", "mean_of_Mpred_per_age\nDivsum_mortality_rates",
#               "mean_of_F_per_age\nDivsum_mortality_rates", "mean_of_Madd_per_age\nDivsum_mortality_rates")
#     data_split <- df_trait_sp_a[which(df_trait_sp_a$variable %in% var_i ),]
#     ref_data_split <- ref_data[which(ref_data$variable %in% var_i),]
#     data_split$variable <- improve_vriable_names(data_split$variable)
#     var_i <- improve_vriable_names(var_i)
#     ref_data_split$variable <- improve_vriable_names(ref_data_split$variable)
#     if (all(is.na(data_split$value))) next
#     result <- plot_trait_vs_fmax(df_trait_sp = data_split, trait_name = var_i,
#                                  Time1=150, Time2=160, relative_to_ref = relative_to_ref, spec=sp,
#                                  age_dim=T, reference_input=ref_data_split, proportional = proportional,
#                                  weight=weight_df,
#                                  chosen_title=paste0(ifelse(proportional, "Proportional ", ""),
#                                                      ifelse(relative_to_ref, "Relative ", ""),
#                                                      "Mortality rates\naverage years ", Time1, " to ",Time2))
#     mortality_plot_age[[paste0(as.character(relative_to_ref), as.character(proportional))]][[sp]][[stringr::str_remove(paste(var_i, collapse=""), pattern="\n")]]  <- result$plot
# 
#     result =NULL
#     var_i = c("Mstar_deaths\nDivcolSumall_deaths", "F_deaths\nDivcolSumall_deaths",
#               "Mpred_deaths\nDivcolSumall_deaths", "Madd_deaths\nDivcolSumall_deaths")
#     data_split <- df_trait_sp_a[which(df_trait_sp_a$variable %in% var_i ),]
#     ref_data_split <- ref_data[which(ref_data$variable %in% var_i),]
#     data_split$variable <- improve_vriable_names(data_split$variable)
#     var_i <- improve_vriable_names(var_i)
#     ref_data_split$variable <- improve_vriable_names(ref_data_split$variable)
#     if (all(is.na(data_split$value))) next
#     result <- plot_trait_vs_fmax(df_trait_sp = data_split, trait_name = var_i,
#                                  Time1=150, Time2=160, relative_to_ref = relative_to_ref, spec=sp,
#                                  age_dim=T, reference_input=ref_data_split, proportional = proportional,
#                                  weight=weight_df,
#                                  chosen_title=paste0(ifelse(proportional, "Proportional ", ""),
#                                                      ifelse(relative_to_ref, "Relative ", ""),
#                                                      "Death number per Mortality rate\naverage years ", Time1, " to ",Time2))
#     mortality_plot_age[[paste0(as.character(relative_to_ref), as.character(proportional))]][[sp]][[stringr::str_remove(paste(var_i, collapse=""), pattern="\n")]]  <- result$plot
# 
#     # # save plot
#     # suffix_relative <- paste0(ifelse(relative_to_ref, "_RELATIVE", ""))
#     # if (!dir.exists(paste0(PATHOUT, "/byage", sp))) dir.create(paste0(PATHOUT, "/byage", sp))
#     # if (is.null(p[[sp]])) next
#     # pdf(file = paste0(PATHOUT, "/byage", sp, "/VxMortal_", "data_Ne", sp, "vsFmax", Fmax_i, suffix_relative, ".pdf"),
#     #     width = 15, height = 8)
#     # for (trait_plotted in names(p[[sp]])){
#     #   plot_list <- c(p[[sp]][[trait_plotted]], list(get(plot_biom)[[sp]]))
#     #   plot_list <- lapply(plot_list, function(x) x + theme_minimal(base_size = 11) + theme(legend.position = "bottom"))
#     #   combined_plot <- ggarrange(plotlist = plot_list, ncol = 3,
#     #                              nrow = ceiling((length(plot_list))/3),
#     #                              common.legend = TRUE, legend = "bottom")
#     #   print(combined_plot)
#     # }
#     # dev.off()
#   }#relative
# }#sp
# save(mortality_plot_age, file=paste0(pathoutputRDATA, "/mortalityplots_age.RData"))
load(file=paste0(pathoutputRDATA, "/mortalityplots_age.RData"))


#### PLOT : combine all and SAVE ----------------
for (relative_to_ref in c(TRUE, FALSE, "proportional")){
  proportional = FALSE
  if (relative_to_ref == "proportional") {
    relative_to_ref = FALSE
    proportional = TRUE
  }
  relative_to_ref <- as.logical(relative_to_ref)
  
  plot_name_ch=c()
  plot_name <- paste0(as.character(relative_to_ref), as.character(proportional))
  if (plot_name=="TRUEFALSE") plot_name_ch <- c("relative")
  if (plot_name=="FALSEFALSE") plot_name_ch <- c("")
  if (plot_name=="FALSETRUE") plot_name_ch <- c("proportional")
  
  plot_biom = ifelse(relative_to_ref, "rel_biomassplot", "biomassplot")
  plot_biom = ifelse(relative_to_ref, "biomassplot", "biomassplot")
  
  
  
  
  windowsFonts(myCustomWindowsFontName=windowsFont("DejaVu Sans"))
  cairo_pdf(file = paste0(PATHOUT, "/Clean_for_paper/Ne_Vk_T_", plot_name_ch, ".pdf"),
      width = 747/72, height = 2*450/72, onefile=T)
  for (sp in species){
    plot_list <- c(Ne_plot_list[[plot_name]][[sp]][2:1])
    # plot_list <- lapply(plot_list, function(x) x + theme_minimal(base_size = 11) + theme(legend.position = "bottom"))
    combined_plot <- ggarrange(plotlist = plot_list,
                               ncol = 1,
                               nrow = 2)
    combined_plot <- annotate_figure(combined_plot, top = text_grob(paste0("Explaining ", sp, " ", plot_name_ch, " variation of Ne/N"),
                                                                    color = "grey60", face = "bold", size = 14))
    print(combined_plot)
  }
  dev.off()
  

  windowsFonts(myCustomWindowsFontName=windowsFont("DejaVu Sans"))
  cairo_pdf(file = paste0(PATHOUT, "/Clean_for_paper/Ne_final", plot_name_ch, ".pdf"),
      width = 747/72, height = 747/72, onefile=T)
  for (sp in species){
    plot_list <- c(Ne_plot_list[[plot_name]][["Whiting"]][2], 
                   Ne_plot_list[[plot_name]][["Saithe"]][2],
                   Ne_plot_list[[plot_name]][["EuropeanSprat"]][2],
                   Ne_plot_list[[plot_name]][["Haddock"]][2])
    # plot_list <- lapply(plot_list, function(x) x + theme_minimal(base_size = 11) + theme(legend.position = "bottom"))
    combined_plot <- ggarrange(plotlist = plot_list,
                               ncol = 2,
                               nrow = 2, common.legend = T)
    combined_plot <- annotate_figure(combined_plot)
    print(combined_plot)
  }
  dev.off()

  
  windowsFonts(myCustomWindowsFontName=windowsFont("DejaVu Sans"))
  cairo_pdf(file = paste0(PATHOUT, "/Clean_for_paper/Vk_final", plot_name_ch, ".pdf"),
            width = 747/72, height = 747/72, onefile=T)
  for (sp in species){
    plot_list <- c(Ne_plot_list[[plot_name]][["Whiting"]][1], 
                   Ne_plot_list[[plot_name]][["Saithe"]][1],
                   Ne_plot_list[[plot_name]][["EuropeanSprat"]][1],
                   Ne_plot_list[[plot_name]][["Haddock"]][1])
    # plot_list <- lapply(plot_list, function(x) x + theme_minimal(base_size = 11) + theme(legend.position = "bottom"))
    combined_plot <- ggarrange(plotlist = plot_list,
                               ncol = 2,
                               nrow = 2, common.legend = T)
    combined_plot <- annotate_figure(combined_plot)
    print(combined_plot)
  }
  dev.off()
  
  
  

  pdf(file = paste0(PATHOUT, "/Clean_for_paper/Deltasquare_Death", plot_name_ch, ".pdf"),
      width = 10, height = 10)
  for (sp in species){
    plot_list <- c(SSDplot[[plot_name]][[sp]][1],
                  deaths_deltasquare_byage[[plot_name]][[sp]][1:2])
    # plot_list <- lapply(plot_list, function(x) x + theme_minimal(base_size = 11) + theme(legend.position = "bottom"))
    combined_plot <- ggarrange(plotlist = plot_list,
                               ncol = 2,
                               nrow = 2)
    combined_plot <- annotate_figure(combined_plot, top = text_grob(paste0("Explaining ", sp, " ", plot_name_ch, " variation of Ne/N"),
                                                                    color = "grey60", face = "bold", size = 14))
    print(combined_plot)
  }
  dev.off()

  
  pdf(file = paste0(PATHOUT, "/Clean_for_paper/Deltasquare_Death", plot_name_ch, ".pdf"),
      width = 10, height = 10)
  for (sp in species){
    plot_list <- c(SSDplot[[plot_name]][[sp]][1],
                   deaths_deltasquare_byage[[plot_name]][[sp]][1:3])
    # plot_list <- lapply(plot_list, function(x) x + theme_minimal(base_size = 11) + theme(legend.position = "bottom"))
    combined_plot <- ggarrange(plotlist = plot_list,
                               ncol = 2,
                               nrow = 2)
    combined_plot <- annotate_figure(combined_plot, top = text_grob(paste0("Explaining ", sp, " ", plot_name_ch, " variation of Ne/N"),
                                                                    color = "grey60", face = "bold", size = 14))
    print(combined_plot)
  }
  dev.off()
  

  # pdf(file = paste0(PATHOUT, "/Clean_for_paper/Mortality", plot_name_ch, ".pdf"),
  #     width = 8, height = 10)
  # for (sp in species){
  #   plot_list <- c(mortality_plot[[plot_name]][[sp]])
  #   # plot_list <- lapply(plot_list, function(x) x + theme_minimal(base_size = 11) + theme(legend.position = "bottom"))
  #   combined_plot <- ggarrange(plotlist = plot_list,
  #                              ncol = 2,
  #                              nrow = 2)
  #   combined_plot <- annotate_figure(combined_plot, top = text_grob(paste0("Explaining ", sp, " ", plot_name_ch, " variation of Ne/N"),
  #                                                                   color = "grey60", face = "bold", size = 14))
  #   print(combined_plot)
  # }
  # dev.off()  
  # 
  # pdf(file = paste0(PATHOUT, "/Clean_for_paper/Mortality_rates_", plot_name_ch, ".pdf"),
  #     width = 8, height = 10)
  # for (sp in species){
  #   plot_list <- c(mortality_plot_age[[plot_name]][[sp]])
  #   # plot_list <- lapply(plot_list, function(x) x + theme_minimal(base_size = 11) + theme(legend.position = "bottom"))
  #   combined_plot <- ggarrange(plotlist = plot_list,
  #                              ncol = 2,
  #                              nrow = 2)
  #   combined_plot <- annotate_figure(combined_plot, top = text_grob(paste0("Explaining ", sp, " ", plot_name_ch, " variation of Ne/N"),
  #                                                                   color = "grey60", face = "bold", size = 14))
  #   print(combined_plot)
  # }
  # dev.off()

  # pdf(file = paste0(PATHOUT, "/Clean_for_paper/Survival_Biomass", plot_name_ch, ".pdf"),
  #     width = 10, height = 10)
  # for (sp in species){
  #   plot_list <- c(deaths_deltasquare_byage[[plot_name]][[sp]][3],
  #                  list(get(plot_biom)[[sp]]))
  #   # plot_list <- lapply(plot_list, function(x) x + theme_minimal(base_size = 11) + theme(legend.position = "bottom"))
  #   combined_plot <- ggarrange(plotlist = plot_list,
  #                              ncol = 1,
  #                              nrow = 2)
  #   combined_plot <- annotate_figure(combined_plot, top = text_grob(paste0("Explaining ", sp," ", plot_name_ch, " variation of Ne/N"),
  #                                                                   color = "grey60", face = "bold", size = 14))
  #   print(combined_plot)
  # }
  # dev.off()
}

  
# 
# # 
# # 
# # 
# # ##### PLOT 4 : separately : Survival and mortality by age -------
# 
# pathout <- paste0("osmose-ns-v4.4.1_clim/output/PLOTS/Trait_variables_by_cohort_by_age2/otherplots")
# if (!dir.exists(pathout)) dir.create(pathout, recursive = T)
# 
# # load dataset
# data_complete_sp = NULL
# data_complete_c = NULL
# base_path = "osmose-ns-v4.4.1_clim/output/PLOTS/Trait_variables_by_cohort_by_age2"
# folder_name <- basename(base_path)
# suffix <- gsub("Trait_variables_by_", "", folder_name)
# suffix <- gsub("_", "_", suffix)  # keep underscores
# save_dir <- file.path(dirname(dirname(base_path)), "RDATA")
# save_file <- file.path(save_dir, paste0("summary_trait_", suffix, ".RData"))
# if (!dir.exists(save_dir)){dir.create(save_dir, recursive = T)}
# FORCE=F
# if (!file.exists(save_file) || FORCE==T){
#   data_complete_c <- data_traits_summary(base_path)
# } else {
#   data_complete_c = get(load(file=save_file))
# }
# 
# # dimensions dataset required : [time, Sp, value, rep, Fish_press, TRAIT]
# head(data_complete_c)
# data_complete_c <- data_complete_c[, which(colnames(data_complete_c) %in% c("age", "time", "Species", "value", "TRAIT", "replicate", "Fish_press"))]
# data_complete_c[grep(data_complete_c$TRAIT, pattern="fec"), "TRAIT"] <- "mean_of_eggNumber_per_age"
# unique(data_complete_c$TRAIT)
# 
# # adding natural mortality dimension
# data_complete_c[which((data_complete_c$TRAIT) == "mean_of_Mnat_per_age"), "TRAIT"] <- "mean_of_Madd_per_age"
# data_complete_c <- pivot_wider(data = data_complete_c, names_from = "TRAIT", values_from = "value")
# data_complete_c <- data_complete_c[,-grep(colnames(data_complete_c), pattern="mean_of_Dis_per_age")]
# data_complete_c <- data_complete_c[,-grep(colnames(data_complete_c), pattern="mean_of_Z_per_age")]
# data_complete_c$mean_of_Mnat_per_age <- data_complete_c$mean_of_Madd_per_age + data_complete_c$mean_of_Mpred_per_age + data_complete_c$mean_of_Mstar_per_age
# data_complete_c <- pivot_longer(data = data_complete_c, names_to = "TRAIT", values_to = "value", cols=-colnames(data_complete_c)[1:4])
# 
# # changing colnames to match expectations
# colnames(data_complete_c)[which(colnames(data_complete_c) == "TRAIT")] <- "variable"
# data_complete_c$Fmax   <- as.numeric(unlist(stringr::str_remove(unlist(stringr::str_split(data_complete_c$Fish_press, " Fstart"))[seq(1, 2*length(data_complete_c$Fish_press), by=2)], "Fmax=")))
# data_complete_c$Fstart <- as.numeric(unlist(stringr::str_split(data_complete_c$Fish_press, " Fstart="))[seq(2, 2*length(data_complete_c$Fish_press), by=2)])
# data_complete_c <- (data_complete_c[,-which(colnames(data_complete_c) == "Fish_press")])
# 
# 
# 
# PATHOUT  = "osmose-ns-v4.4.1_clim/output/PLOTS/KEY_TRAITS/"
# surv_Mpred = list()
# data_full <- NULL
# for (relative_to_ref in c(TRUE, FALSE)){
# 
#   # for the plot
#   plot_biom = ifelse(relative_to_ref, "rel_biomassplot", "biomassplot")
#   plot_biom = ifelse(relative_to_ref, "biomassplot", "biomassplot")
# 
#   for (sp in species){
#     cat(sp)
# 
#     df_trait_sp <- data_complete_c %>% filter(Species == sp)
#     df_trait_sp <- df_trait_sp %>% filter(Fstart == 1)
# 
#     df_trait_t    <- estimateSum(trait_data_full_i = df_trait_sp, var=paste0("mean_of_F_per_age", "\nsum", "mean_of_Mnat_per_age"))
#     df_trait_t$variable    <- "sum_mortality_rates"
#     df_trait_sp   <- rbind(df_trait_sp, df_trait_t)
# 
#     df_trait_t    <- estimateSurv(trait_data_full_i = df_trait_sp)
#     df_trait_sp   <- rbind(df_trait_sp, df_trait_t)
# 
#     var_i    <- unique(df_trait_sp$variable)
#     # compare all mortality rates to their value *abundance to have the number of individuals to understand where is the increase coming from
#     for (v in var_i){
#       if (any(grepl(pattern="mean_of_M", v) | grepl(pattern="mean_of_F", v) | grepl(pattern="sum_mortality_rates", v))) {
#         df_trait_t    <- estimateMnumber(trait_data_full_i = df_trait_sp, var=paste0(v))
#         df_trait_sp   <- rbind(df_trait_sp, df_trait_t)
#       }
#     }
#     var_i    <- unique(df_trait_sp$variable)
#     for (v in var_i){
#       if (any(grepl(pattern="deaths", v))) {
#         df_trait_t    <- estimatecolSum(trait_data_full_i = df_trait_sp, var=paste0("colSum", v))
#         df_trait_sp   <- rbind(df_trait_sp, df_trait_t)
#       }
#     }
#     var_i    <- unique(df_trait_sp$variable)
# 
#     # reduce dimension number
#     weight_df <- df_trait_sp %>%
#       dplyr::filter(variable == "mean_of_abundanceByAge_per_age") %>%
#       dplyr::select(Species, replicate, Fstart, Fmax, age, weight_value = value)
#     data_with_w <- df_trait_sp %>%
#       dplyr::left_join(weight_df, by = c("Species", "replicate", "age", "Fmax", "Fstart"), relationship = "many-to-many")
#     data_with_w <- data_with_w[-which(data_with_w$weight_value == 0), ]
#     data_with_w$age <- as.numeric(data_with_w$age)
#     if(any(data_with_w$age == 0)) data_with_w <- data_with_w[-which(data_with_w$age == 0), ]
#     df_trait_sp_a <- data_with_w %>%
#       dplyr::group_by(Species, replicate, Fstart, Fmax, variable, age) %>%
#       dplyr::summarise(
#         value = stats::weighted.mean(value, w = weight_value, na.rm = TRUE),
#         .groups = "drop")
# 
#     ref_data <- df_trait_sp_a %>%
#       filter(Fmax == 1 & Fstart == 1) %>%
#       dplyr::select(variable, value, replicate, age) %>%
#       dplyr::rename(ref_value = value)
# 
# 
#     if (sum(is.nan(df_trait_sp_a$value))>0) df_trait_sp_a <- df_trait_sp_a[-which(is.nan(df_trait_sp_a$value)), ]
#     if (sum(is.infinite(df_trait_sp_a$value))>0) df_trait_sp_a <- df_trait_sp_a[-which(is.infinite(df_trait_sp_a$value)), ]
# 
#     # conditions to stop
#     if (nrow(df_trait_sp_a)==0) next
#     if (all(is.na(df_trait_sp_a$value) | df_trait_sp_a$value == 0)) {
#       next
#     }
#     df_all <- df_trait_sp_a
#     for (Fmax_i in c("0.5", "0.75", "1", "1.25", "1.5")){
#       df_trait_sp_a <- df_all[which(df_all$Fmax == Fmax_i),]
# 
#       var_i = c("mean_of_Mpred_per_age", "Mpred_deaths", "surv_curve")
#       data_split <- df_trait_sp_a[which(df_trait_sp_a$variable %in% var_i ),]
#       ref_data_split <- ref_data[which(ref_data$variable %in% var_i),]
#       if (all(is.na(data_split$value))) next
#       result <- plot_trait_vs_fmax(df_trait_sp = data_split, trait_name = var_i,
#                                    Time1=150, Time2=160, relative_to_ref = relative_to_ref, spec=sp,
#                                    age_dim=T, reference_input=ref_data_split)
#       surv_Mpred[[sp]][[paste(var_i, collapse="")]][[Fmax_i]]    <- result$plot
# 
#       data_split$relative_to_ref <- relative_to_ref
#       if (is.null(data_full)){
#         data_full <- data_split
#       } else {
#         data_full <- rbind(data_full,  data_split)
#       }
#     }
# 
#     # save plot
#     suffix_relative <- paste0(ifelse(relative_to_ref, "_RELATIVE", ""))
#     if (!dir.exists(paste0(PATHOUT, "/byage", sp))) dir.create(paste0(PATHOUT, "/byage", sp))
#     pdf(file = paste0(PATHOUT, "/byage", sp, "/surv_Mpred_", "data_Ne", sp, "vsFmax", Fmax_i, suffix_relative, ".pdf"),
#         width = 15, height = 8)
#     for (trait_plotted in names(surv_Mpred[[sp]])){
#       plot_list <- c(surv_Mpred[[sp]][[trait_plotted]], list(get(plot_biom)[[sp]]))
#       plot_list <- lapply(plot_list, function(x) x + theme_minimal(base_size = 11) + theme(legend.position = "bottom"))
#       combined_plot <- ggarrange(plotlist = plot_list, ncol = 3,
#                                  nrow = ceiling((length(plot_list))/3),
#                                  common.legend = TRUE, legend = "bottom")
#       print(combined_plot)
#     }
#     dev.off()
#   }#sp
# }
# full_mort_surv_plot = list()
# for (relative_to_ref in c(TRUE, FALSE)){
#   # for the plot
#   plot_biom = ifelse(relative_to_ref, "rel_biomassplot", "biomassplot")
#   plot_biom = ifelse(relative_to_ref, "biomassplot", "biomassplot")
#   suffix_relative <- paste0(ifelse(relative_to_ref, "_RELATIVE", ""))
#   data_full_ref <- data_full[which(data_full$relative_to_ref == relative_to_ref), ]
#   for (trait_plotted in names(table(data_full$variable))){
#     data_full_i <- data_full_ref[which(data_full_ref$variable == trait_plotted), ]
#     ref_data <- data_full_i %>%
#       filter(Fmax == 1 & Fstart == 1) %>%
#       dplyr::select(variable, value, replicate, age, Species) %>%
#       dplyr::rename(ref_value = value)
#     for (Fmax_i in c("0.5","0.75","1","1.25","1.5")){
#       data_full_t <- data_full_i[which(data_full_i$Fmax == Fmax_i), ]
#       result <- plot_trait_vs_fmax(df_trait_sp = data_full_t, trait_name = trait_plotted,
#                                    Time1=150, Time2=160, relative_to_ref = F, spec=species,
#                                    age_dim=T, sp_dim = T, reference_input=ref_data)
#       full_mort_surv_plot[[as.character(relative_to_ref)]][[trait_plotted]][[Fmax_i]] <- result$plot
#     }
#   }#trait
#   pdf(file = paste0(PATHOUT, "/comparing_sp_traits_", "data_Ne", sp, "vsFmax", Fmax_i, suffix_relative, ".pdf"),
#       width = 15, height = 12)
#   for (trait_plotted in names( full_mort_surv_plot[[as.character(relative_to_ref)]])){
#     combined_plot <- ggarrange(plotlist = full_mort_surv_plot[[as.character(relative_to_ref)]][[trait_plotted]], ncol = 3,
#                                nrow = ceiling((length(full_mort_surv_plot[[as.character(relative_to_ref)]][[trait_plotted]]))/3),
#                                common.legend = TRUE, legend = "bottom")
#     print(combined_plot)
#   }
#   dev.off()
# }# relaitve



abund <- data_complete_c[which(data_complete_c$variable == "mean_of_abundanceByAge_per_age"),]
abund$Fmax <- as.factor(abund$Fmax)
ggplot(abund, aes(x=age, y=value, group=interaction(Fmax, Species), color=Fmax)) + 
  geom_line() +
  theme_bw() +
  facet_wrap(~Species)

max_age_sp <- abund %>%
  group_by(Species, replicate, Fmax) %>%
  summarize(age_min = abund[which(abund$value == 0)[1], "age"])
logabund <- abund
logabund <- logabund[-which(logabund$value == 0),]
logabund$value <- log(logabund$value)
logabund$replicate <- as.factor(logabund$replicate)
logabund$value[which(is.infinite(logabund$value))] <- 0
# logabund <- logabund %>%
#   group_by(Species, Fmax, age) %>%
#   summarize(mean=mean(value))
logabund <- logabund[which(logabund$Fstart == 1), ]
p <- list()
for (rep in levels(logabund$replicate)){
  p[[rep]] <- ggplot(logabund[which(logabund$replicate==rep),], aes(x=age, y=value, group=interaction(Fmax, Species), color=Fmax)) + 
    geom_line() +
    theme_bw() +
    facet_wrap(~Species, scales = "free") + 
    labs(title = "mean Log abundance by age per species and Fmax")
}

pdf(file = paste0(PATHOUT, "/Clean_for_paper/abundance_by_sp.pdf"),
    width=12, height = 9)
print(p)
dev.off()
