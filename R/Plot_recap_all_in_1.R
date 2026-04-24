# the idea is to store separately all the key parameters per species into plots
# and to produce one pdf sheet per species with all the key variables to have it all summarized at one glance




#### Initialisation --------------
library(legendry)
library(ggpubr)

# functions and environment
main = "New_gpu_conf.R" # config file
source("R/functions/source.R") # load config and key parameters and functions
source("R/functions/functions_key_traits.R") # plots against Fmax
source("R/functions/auxiliary_for_lists.R") # from list to df
# source("R/functions/functions_output_SimuRData_into_df_for_plots.R") # to transform dataset from RData to long df
source("R/functions/function_data_traits_summary.R")

# paths
simu_repo <- "C:/Users/mbeneat/Documents/osmose/Chap2/osmose-ns/osmose-ns-v4.4.1_clim/"
path_long <- paste0(simu_repo, "/output/PLOTS/FINAL_Ne_Htz_mature/")
path <- paste0(simu_repo, "/output/PLOTS/FINAL_Ne_Htz_mature_short/")
PATHOUT  = "output/PLOTS/KEY_TRAITS_MATURE/"
if (!dir.exists(PATHOUT)) dir.create(PATHOUT, recursive = TRUE)
pathoutputRDATA = paste0(PATHOUT, "/RDATA_MATURE")
if(!dir.exists(pathoutputRDATA)) dir.create(pathoutputRDATA, recursive = T)
base_dir <- paste0(path, "INIT")
folders  <- list.dirs(base_dir, recursive = F) # files INIT2 with raw Ne datasets

folders_traits_y <- paste0(simu_repo, "output/PLOTS/Trait_variables_by_year_BYSP/RData/")
folders_traits_y <- list.dirs(folders_traits_y, recursive = F)
folders_traits_c <- paste0(simu_repo, "output/PLOTS/Trait_variables_by_cohort/RData/")
folders_traits_c  <- list.dirs(folders_traits_c, recursive = F)
folders_traits_y_a <- paste0(simu_repo, "output/PLOTS/Trait_variables_by_year_by_age2/RData/")
folders_traits_y_a <- list.dirs(folders_traits_y_a, recursive = F)
folders_traits_c_a <- paste0(simu_repo, "output/PLOTS/Trait_variables_by_cohort_by_age2/RData/")
folders_traits_c_a  <- list.dirs(folders_traits_c_a, recursive = F)

#### DATASETS load or format dataset ------------------------
Time1=150
Time2=160
species_df=species
# DATA NE and related traits, dimensions:, replicate, Fstart, Fmax, Species
### df by loading the trait datasets by Fish_press (by year, by cohort, by year by age, by cohort by age) and binding them together
# list_data_traits_ratio <- data_trait_and_genetic_into_df(folders_in = folders, Time1=150, Time2=160, pathoutput=pathoutputRDATA,
#                                                          FORCE = F, species = species_df, RATIO=T)
# 
# # DATA NE + TRAITS dimensions : (age), replicate, Fstart, Fmax, Species
# dataRATIO_full_merged_age <- dataset_attach_all_traits_genetic_drift(pathoutput=pathoutputRDATA, list_outputs=list("ratio"=list_data_traits_ratio),
#                                                                  trait_cv = list("folders_traits_c_a"=folders_traits_c_a, "folders_traits_y_a"=folders_traits_y_a),
#                                                                  FORCE=F, keep_cohort_or_age_dimension = T)
# dataRATIO_full_merged <- dataset_attach_all_traits_genetic_drift(pathoutput=pathoutputRDATA, list_outputs=list("ratio"=list_data_traits_ratio),
#                                                                  trait_cv = list("folders_traits_c"=folders_traits_c, "folders_traits_y"=folders_traits_y),
#                                                                  FORCE=F, keep_cohort_or_age_dimension = F)
# RAW Ne & associated variables, dimensions : cohort, age, replicate, Fstart, Fmax, Species
# ### produce dataset reading multiple .RData to make 1 dataframe, stored as a .RData file ----------
# sp_i = "Hake"
# FORCE=F
# full_dataset = list()
# full_dataset_age = list()
# for (f_i in folders){
#   cat(f_i, "\n")
# 
#   Fmax   <- as.numeric(stringr::str_split(stringr::str_split(f_i, pattern="Fmax")[[1]][2], pattern="_Fstart")[[1]][1])
#   Fstart <- as.numeric(stringr::str_split(stringr::str_split(f_i, pattern="Fstart")[[1]][2], pattern="_sp")[[1]][1])
# 
#   if (Fstart != 1) next
# 
#   load(list.files(f_i, full.names = T))
#   dataNe <- INIT2$data_Ne
# 
#   for (sp in names(dataNe)[1]){
#     if (file.exists(paste0(path, "/Complete_INIT2", sp_i, ".RData")) && FORCE !=T){
#       next
#     }else{
# 
#       cat(sp, "\n")
# 
#       f_i_sp <- dataNe[[sp]]
# 
#       for (dim in names(f_i_sp)){
# 
#         f_i_sp_d <- f_i_sp[[dim]]
# 
#         if (is.data.frame(f_i_sp_d) || is.matrix(f_i_sp_d)){
# 
#           if (is.matrix(f_i_sp_d)) f_i_sp_d <- as.data.frame(f_i_sp_d)
#           f_i_sp_d$cohort <- rownames(f_i_sp_d)
#           f_i_sp_d <- pivot_longer(f_i_sp_d, names_to = "replicate", values_to = "value", cols=-cohort)
#           f_i_sp_d$Species <- sp
#           f_i_sp_d$variable <- dim
#           f_i_sp_d$Fmax <- Fmax
#           f_i_sp_d$Fstart <- Fstart
# 
#           if (is.null(full_dataset[[sp]])){
#             full_dataset[[sp]] <- f_i_sp_d
#           } else {
#             full_dataset[[sp]] <- rbind(full_dataset[[sp]], f_i_sp_d)
#           }
# 
#         } # if matrix or df
#         if (!is.data.frame(f_i_sp_d) & !is.matrix(f_i_sp_d)){
# 
#           f_i_sp_d = lapply(names(f_i_sp_d), FUN = function(x) {
#             f_i_sp_d[[x]] <- as.data.frame(f_i_sp_d[[x]])
#             f_i_sp_d[[x]][, "replicate"] = stringr::str_remove(x, pattern="rep_")
#             f_i_sp_d[[x]][, "age"] = rownames(f_i_sp_d[[x]])
#             return(f_i_sp_d[[x]])
#           })
#           f_i_sp_d = do.call(rbind, f_i_sp_d)
# 
#           f_i_sp_d <- pivot_longer(f_i_sp_d, names_to = "cohort", values_to = "value", cols=-c(replicate, age))
#           f_i_sp_d[, "cohort"] = stringr::str_remove(as.matrix(f_i_sp_d[, "cohort"]), pattern="cohort_")
#           f_i_sp_d$Species <- sp
#           f_i_sp_d$variable <- dim
#           f_i_sp_d$Fmax <- stringr::str_split(stringr::str_split(f_i, pattern="Fmax")[[1]][2], pattern="_Fstart")[[1]][1]
#           f_i_sp_d$Fstart <- stringr::str_split(stringr::str_split(f_i, pattern="Fstart")[[1]][2], pattern="_sp")[[1]][1]
# 
#           if (is.null(full_dataset_age[[sp]])){
#             full_dataset_age[[sp]] <- f_i_sp_d
#           } else {
#             full_dataset_age[[sp]] <- rbind(full_dataset_age[[sp]], f_i_sp_d)
#           }
# 
#         } # if list
#       }# variable dimension
#     }# only if file does not exist yet
#   }# species dimension
# }# folder dimension (Fmax*Fstart)
# #----- MERGE raw variables without age dimension AND variables with age dimension ---#
# for (sp_i in names(full_dataset_age)){
#   if (file.exists(paste0(path, "/Complete_INIT2", sp_i, ".RData")) && FORCE !=T){
#     next
#   }else{
# 
#   # produce a dataset without age dimension : weighted mean
#   weight_df <- full_dataset_age[[sp_i]] %>%
#     dplyr::filter(variable == "abundanceByAge") %>%
#     dplyr::select(Species, replicate, Fstart, Fmax, cohort, age, weight_value = value)
#   data_with_w <- full_dataset_age[[sp_i]] %>%
#     dplyr::left_join(weight_df, by = c("Species", "replicate", "cohort", "age", "Fmax", "Fstart"), relationship = "many-to-many")
#   data_with_w <- data_with_w[-which(data_with_w$weight_value == 0), ]
#   data_with_w$age <- as.numeric(data_with_w$age)
#   data_out <- data_with_w %>%
#     dplyr::group_by(Species, replicate, Fstart, Fmax, cohort, variable) %>%
#     dplyr::summarise(
#       value = stats::weighted.mean(value, w = weight_value, na.rm = TRUE),
#       .groups = "drop"
#     )
#   final_dataset <- rbind(full_dataset[[sp_i]], data_out)
# 
#   # produce a dataste with the age dimension
#   full_dataset_add_age <- lapply(names(table(full_dataset_age[[sp_i]]$age)), FUN = function(x){
#     out = as.data.frame(full_dataset[[sp_i]])
#     out$age = x
#     return(out)
#   })
#   full_dataset_add_age <- do.call(rbind, full_dataset_add_age)
# 
#   final_dataset_age <- rbind(full_dataset_add_age, full_dataset_age[[sp_i]])
# 
#   # save all
#   complete = list(final_dataset, final_dataset_age)
# 
#   save(complete, file=paste0(path, "/Complete_INIT2", sp_i, ".RData"))
#   }
# }



# RAW TRAITS, dimensions :  age, replicate, Fish_press, Species, TRAIT
base_path = paste0(simu_repo, "output/PLOTS/Trait_variables_by_cohort_by_age2")
pathout   = paste0(base_path, "/otherplots")
if (!dir.exists(pathout)) dir.create(pathout, recursive = T)

# load dataset
data_complete_sp = NULL
data_complete_c = NULL
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

# avoid the estimation bias of very populations of very low size
# if (any(((data_complete_c[which(data_complete_c$TRAIT == "mean_of_abundanceByAge_per_age"), "value"])<10))){
#   data_complete_c <- pivot_wider(data_complete_c, names_from = "TRAIT", values_from = "value")
#   data_complete_c <- data_complete_c[-which((data_complete_c$mean_of_abundanceByAge_per_age)<10), ]
#   data_complete_c <- pivot_longer(data_complete_c, cols = -c("age", "replicate", "Fish_press", 'Species'), 
#                                   names_to = "TRAIT", values_to = "value")
# }

# removing age 0
data_complete_c <- data_complete_c[-which(data_complete_c$age == 0), ]
# correction
data_complete_c[grep(data_complete_c$TRAIT, pattern="fec"), "TRAIT"] <- "mean_of_eggNumber_per_age"

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




# _  ------------------
#### PLOT 0- Biomass vs Fmax & keep it on the side -------------------
biomassplot = NULL
rel_biomassplot = NULL
list_relative_biomass <- list()
for (nm in "data_Ne"){
  d <- data_complete_c
  d <- d[d$Fstart==1, ]
  id_biom <- c(grep(pattern='biomass', d$variable), grep(pattern='abundance', d$variable))
  for (sp in species){
    d_i_i <- d[id_biom, ]
    d_i_i <- d_i_i[which(d_i_i$Species == sp), ]
    d_i_i <- pivot_wider(d_i_i, values_from = "value", names_from = "variable")
    if (any(is.na(d_i_i$mean_of_biomassByAge_per_age))) d_i_i <- d_i_i[-which(is.na(d_i_i[, "mean_of_biomassByAge_per_age"])), ]
    d_i_i$Iszero <- ifelse(d_i_i$mean_of_abundanceByAge_per_age < 100, "Abundance <100", "Abundance >100")

    biomassplot[[sp]] <- ggplot(d_i_i, aes(y=mean_of_biomassByAge_per_age, x=Fmax)) +
      geom_point(aes(colour = Iszero))+
      labs(title = "Biomass per fishing pressure", subtitle = sp)

    ref_data <- d_i_i %>%
      filter(Fmax == 1 & Fstart == 1) %>%
      dplyr::select(mean_of_biomassByAge_per_age, replicate, age) %>%
      dplyr::rename(ref_value = mean_of_biomassByAge_per_age)
    ref_data_percent <- d_i_i %>%
      filter(Fmax == 0.5 & Fstart == 1) %>%
      dplyr::select(mean_of_biomassByAge_per_age, replicate, age, mean_of_abundanceByAge_per_age) %>%
      dplyr::rename(ref_value_percent = mean_of_biomassByAge_per_age) %>%
      dplyr::rename(ref_value_percent_abund = mean_of_abundanceByAge_per_age) 
    # Subtract reference from all means
    rel_d_i_i <- d_i_i %>%
      left_join(ref_data, by = c("replicate", "age"), relationship = "many-to-one") %>%
      mutate(relative_biomass = (mean_of_biomassByAge_per_age - ref_value)/ref_value) %>%
      dplyr::select(-ref_value)    
    rel_d_i_i_perc <- d_i_i %>%
      left_join(ref_data_percent, by = c("replicate", "age"), relationship = "many-to-one") %>%
      mutate(relative_biomass_percent = 100*mean_of_biomassByAge_per_age/ref_value_percent-100) %>%
      mutate(relative_abundance_percent = 100*mean_of_abundanceByAge_per_age/ref_value_percent_abund-100) %>%
      dplyr::select(-ref_value_percent)
    list_relative_biomass[[sp]] <- rel_d_i_i_perc
    rel_d_i_i$Iszero_relative <- ifelse(abs(rel_d_i_i$relative_biomass) < 0.001, "Relative value is 0", "Different from 0")
    rel_biomassplot[[sp]] <- ggplot(rel_d_i_i, aes(y=relative_biomass, x=Fmax)) +
      geom_point(aes(colour = Iszero_relative))+
      labs(title = "Relative biomass per fishing pressure", subtitle = sp)
    
  }
}
list_bio <- list("rel_biomassplot"=rel_biomassplot, "biomassplot"=biomassplot)
save(list_bio, file=paste0(pathoutputRDATA, "/Biomass_plots.RData"))
load(file=paste0(pathoutputRDATA, "/Biomass_plots.RData"))
rel_biomassplot=list_bio$rel_biomassplot
biomassplot=list_bio$biomassplot

## Estimate the % in decline (in biomass and abundance) from Fmax=1.5 to Fmax = 0.5
# df_relative_biomass <- do.call(rbind, list_relative_biomass)
# summary_df_relative_biomass <- df_relative_biomass %>%
#   group_by(Fmax, Species)%>%
#   # summarize(weighted.mean = weighted.mean(w=mean_of_abundanceByAge_per_age, x=relative_biomass_percent, na.rm=T))
#   summarize(percent_difference_withFmax0.5 = mean(relative_biomass_percent, na.rm=T))
# summary_df_relative_biomass <- summary_df_relative_biomass %>%
#   dplyr::mutate(percent_difference_withFmax0.5=round(percent_difference_withFmax0.5, 2))
# write.csv(summary_df_relative_biomass, file=paste0(PATHOUT, "/Clean_for_paper/Biomass_declinetrends.csv"))
# 
# summary_df_relative_abunance <- df_relative_biomass %>%
#   group_by(Fmax, Species)%>%
#   summarize(percent_difference_withFmax0.5 = mean(relative_abundance_percent, na.rm=T))
# summary_df_relative_abunance <- summary_df_relative_abunance %>%
#   dplyr::mutate(percent_difference_withFmax0.5=round(percent_difference_withFmax0.5, 2))
# write.csv(summary_df_relative_abunance, file=paste0(PATHOUT, "/Clean_for_paper/Abund_declinetrends.csv"))

## Estimate the % in decline (in Ne) from Fmax=1.5 to Fmax = 0.5
# list_relative_Ne <- list()
# for (nm in "data_Ne"){
#   for (sp in species){
#     cat(sp)
#     ## load and filter data
#     load(file=paste0(path, "/Complete_INIT2", sp, ".RData"))
#     final_dataset_age = complete[[2]]
#     final_dataset_age <- final_dataset_age %>%
#       dplyr::mutate(dplyr::across(c(cohort, replicate), ~ stringr::str_remove(.x, "^(cohort_|rep_)")))
#     final_dataset_age <- final_dataset_age[which(final_dataset_age$cohort %in% Time1:Time2),]
#     rm(complete)
#     
#     # remove ages that are not included
#     final_dataset_age$age <- as.numeric(final_dataset_age$age)
#     ages_to_remove <- table(final_dataset_age[which(final_dataset_age$variable == "Dead" & final_dataset_age$value ==0),"age"])
#     ages_to_remove <- ages_to_remove[which(ages_to_remove == max(ages_to_remove))]
#     final_dataset_age <- final_dataset_age[-which(final_dataset_age$age %in% as.numeric(names(ages_to_remove))), ]
#     final_dataset_age = pivot_wider(final_dataset_age, values_from = "value", names_from = "variable")
#     
#     ref_data_percent <- final_dataset_age %>%
#       filter(Fmax == 0.5 & Fstart == 1) %>%
#       dplyr::select(Ne_N_alpha, replicate, cohort, age) %>%
#       dplyr::rename(ref_value_percent = Ne_N_alpha) 
#     # Subtract reference from all means
#     rel_d_i_i_perc <- final_dataset_age %>%
#       left_join(ref_data_percent, by = c("replicate", "age", "cohort"), relationship = "many-to-one") %>%
#       mutate(relative_Ne_percent = 100*Ne_N_alpha/ref_value_percent-100) %>%
#       dplyr::select(-ref_value_percent)
#     list_relative_Ne[[sp]] <- rel_d_i_i_perc
#   }
# }
# df_relative_Ne <- do.call(rbind, list_relative_Ne)
# summary_df_relative_abunance <- df_relative_Ne %>%
#   group_by(Fmax, Species)%>%
#   summarize(percent_difference_withFmax0.5 = mean(relative_Ne_percent, na.rm=T))
# summary_df_relative_abunance <- summary_df_relative_abunance %>%
#   dplyr::mutate(percent_difference_withFmax0.5=round(percent_difference_withFmax0.5, 2))
# write.csv(summary_df_relative_abunance, file=paste0(PATHOUT, "/Clean_for_paper/Ne_declinetrends.csv"), )




#### PLOT 1- NE/N VS FMAX + Biomass VS FMAX ------------------

traits_to_plot <- list(c("mean_generation_length", "mean_Vk_point", "T/Vk."), c("Ne_N_alpha", "mean_Numerator_N_alpha", "Ne_N_alpha\nDivmean_Numerator_N_alpha"), 
                       c("Ne_N_alpha", "mean_Numerator_N_alpha"), c("mean_Numerator_N_alpha", "T/Vk."), 
                       c("Ne_N_alpha\nDivmean_Numerator_N_alpha", "Ne_N_alpha\nDivmean_NT"))


Ne_plot_list = list()
for (sp in species){
  cat(sp)
  
  ## load and filter data
  load(file=paste0(path, "/Complete_INIT2", sp, ".RData"))
  final_dataset_age = complete[[2]]
  final_dataset_age <- final_dataset_age %>%
    dplyr::mutate(dplyr::across(c(cohort, replicate), ~ stringr::str_remove(.x, "^(cohort_|rep_)")))
  final_dataset_age <- final_dataset_age[which(final_dataset_age$cohort %in% Time1:Time2),]
  rm(complete)
  
  # remove ages that are not included
  final_dataset_age$age <- as.numeric(final_dataset_age$age)
  ages_to_remove <- table(final_dataset_age[which(final_dataset_age$variable == "Dead" & final_dataset_age$value ==0),"age"])
  ages_to_remove <- ages_to_remove[which(ages_to_remove == max(ages_to_remove))]
  final_dataset_age <- final_dataset_age[-which(final_dataset_age$age %in% as.numeric(names(ages_to_remove))), ]
  
  # # avoid the estimation bias of very populations of very low size
  # if (any(((final_dataset_age[which(final_dataset_age$variable == "abundanceByAge"), "value"])<10))){
  #   final_dataset_age <- pivot_wider(final_dataset_age, names_from = "variable", values_from = "value")
  #   final_dataset_age <- final_dataset_age[-which((final_dataset_age$abundanceByAge)<10), ]
  #   final_dataset_age <- pivot_longer(final_dataset_age, cols = -c("age", "replicate", "Fmax", "Fstart", 'Species', "cohort"), names_to = "variable", values_to = "value")
  # }
  
  ## dataset
  data <-  final_dataset_age %>% filter(Species == sp, Fstart == 1)
  
  if (is.data.frame(data)) data <- list("data_Ne"=data)
  # for the plot
  if (!any(grepl(names(data), pattern="Ne"))) {
    tempo <- c(names(data[[1]])[1], names(data[[1]])[2], names(data[[1]])[3])
    if (length(data)==3) data <- list(do.call(rbind, list(data[[1]][[names(data[[1]])[1]]], data[[2]][[names(data[[1]])[1]]])),
                                      do.call(rbind, list(data[[1]][[names(data[[1]])[2]]], data[[2]][[names(data[[1]])[2]]])),
                                      do.call(rbind, list(data[[1]][[names(data[[1]])[3]]], data[[2]][[names(data[[1]])[3]]])))
    if (length(data)==1)data <- list(do.call(rbind, list(data[[1]][[names(data[[1]])[1]]], data[[2]][[names(data[[1]])[1]]])))
    names(data) <- tempo[1:length(data)]
  }
  scale_v = list()
  for (genetic_drift_data in (1:length(data))[1]){

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

    summary <- as_long_if_needed(summary, columns_keep_pivot = colnames(summary)[-which(colnames(summary) %in% c("replicate", "Fmax", "Fstart", "Species", "age", "cohort"))])
    df_trait_t <- estimateT_Vk(trait_data_full_i = summary)
    summary <- rbind(summary, df_trait_t)

    # add NT
    df_trait_t    <- estimateX(trait_data_full_i = summary, var=paste0("abundanceByAgexsurvivorship_curve"))
    df_trait_t$variable <- "NT"
    summary       <- rbind(summary, df_trait_t)      
    
    # add the Ne/NT
    df_trait_t    <- estimateDiv(trait_data_full_i = summary, var=paste0("Ne_N_alpha\nDivNT"))
    summary       <- rbind(summary, df_trait_t)

    # add the Ne/Nalpha
    df_trait_t    <- estimateDiv(trait_data_full_i = summary, var=paste0("Ne_N_alpha\nDivNumerator_N_alpha"))
    summary       <- rbind(summary, df_trait_t)

    WEIGHT = summary[which(summary$variable == "abundanceByAge"),]
    colnames(WEIGHT)[which(colnames(WEIGHT) == "value")] <- "weight"


    if (nrow(summary)<1) next
    df_trait_sp <- summary %>% filter(Species == sp, Fstart == 1, variable %in% unlist(traits_to_plot))
    for (relative_to_ref in c(TRUE, FALSE, "proportional")){  
      
      proportional = FALSE
      if (relative_to_ref == "proportional") {
        relative_to_ref = FALSE
        proportional = TRUE
      }
      relative_to_ref = as.logical(relative_to_ref)

      for (traits_selected in traits_to_plot) {

        df_trait <- df_trait_sp[which(df_trait_sp$variable %in% traits_selected),]

        if (sum(is.nan(df_trait$value))>0) df_trait <- df_trait[-which(is.nan(df_trait$value)), ]
        if (sum(is.infinite(df_trait$value))>0) df_trait <- df_trait[-which(is.infinite(df_trait$value)), ]

        # conditions to stop
        if (nrow(df_trait)==0) next
        if (all(is.na(df_trait$value))) {
          Ne_plot_list[[paste0(as.character(relative_to_ref), as.character(proportional))]][[sp]][[traits_selected]] <- NULL
          next}

        # Use the function
        # if (length(var_i)>2){
        k=1
        id_to_plot = NULL
        # for (i in 1:c((length(var_i)+1)/3)){
        # id_to_plot = k:(k+2)
        # if (traits_selected == "N&T/Vk") id_to_plot=c(which(var_i == "T/Vk"),
        #                                               ifelse(i==1, which(var_i != "T/Vk")[1], which(var_i != "T/Vk")[2]))
        df_trait$variable <- improve_vriable_names(df_trait$variable)
        t <- improve_vriable_names(traits_selected)
        if (all(is.na(df_trait$value))) next
        result <- plot_trait_vs_fmax(df_trait_sp = df_trait, trait_name = t, age_dim = F,weight = WEIGHT,
                                     Time1, Time2, relative_to_ref = relative_to_ref, spec=sp, proportional=proportional,
                                     position_legend = "right")
        Ne_plot_list[[paste0(as.character(relative_to_ref), as.character(proportional))]][[sp]][[paste(traits_selected, collapse="")]]<- result$plot
        # scale_v[[sp]][[paste(t, collapse="")]] <- result$scale
        k=k+2

        # } else {
        #   df_trait_sp_i <- df_trait_sp
        #   df_trait_sp_i$variable <- improve_vriable_names(df_trait_sp_i$variable)
        #   t <- improve_vriable_names(var_i)
        #   result <- plot_trait_vs_fmax(df_trait_sp = df_trait_sp, trait_name = t, age_dim = F, weight = WEIGHT,
        #                                Time1, Time2, relative_to_ref = relative_to_ref, spec=sp, proportional = proportional)
        #   Ne_plot_list[[paste0(as.character(relative_to_ref), as.character(proportional))]][[sp]][[t]]       <- result$plot
        #   scale_v[[sp]][[t]] <- result$scale
        # }
      } # trait selected in the list
    }#sp
  } # genetic drift dim
}# relative
save(Ne_plot_list, file=paste0(pathoutputRDATA, "/Ne_plots.RData"))
load(file=paste0(pathoutputRDATA, "/Ne_plots.RData"))


#### PLOT 2- SSDg and SSDg/N1 and by mortality : only the area under the curve ------------


combined_plot =list()
Time1=150
Time2=160
plot_Deathsdelta_age = list()
plot_Abund_age = list()
plots_by_age_delta= list()
inset.tot.F.Nat = list()
plot_Deathsdelta_age2 = list()
plot_ssdg_age <- list()
for (sp in species){
  cat(sp)

  load(file=paste0(path_long, "/Complete_INIT2", sp, ".RData"))
  final_dataset_age = complete[[2]]
  final_dataset_age <- final_dataset_age %>%
    dplyr::mutate(dplyr::across(c(cohort, replicate), ~ stringr::str_remove(.x, "^(cohort_|rep_)")))
  final_dataset_age <- final_dataset_age[which(final_dataset_age$cohort %in% Time1:Time2),]
  rm(complete)
  
  # remove ages that are not included
  final_dataset_age$age <- as.numeric(final_dataset_age$age)
  ages_to_remove <- table(final_dataset_age[which(final_dataset_age$variable == "Dead" & final_dataset_age$value ==0),"age"])
  ages_to_remove <- ages_to_remove[which(ages_to_remove == max(ages_to_remove))]
  final_dataset_age <- final_dataset_age[-which(final_dataset_age$age %in% as.numeric(names(ages_to_remove))), ]
  
  # avoid the estimation bias of very populations of very low size
  # if (any(((final_dataset_age[which(final_dataset_age$variable == "abundanceByAge"), "value"])<10))){
  #   final_dataset_age <- pivot_wider(final_dataset_age, names_from = "variable", values_from = "value")
  #   final_dataset_age <- final_dataset_age[-which((final_dataset_age$abundanceByAge)<10), ]
  #   final_dataset_age <- pivot_longer(final_dataset_age, cols = -c("age", "replicate", "Fmax", "Fstart", 'Species', "cohort"), names_to = "variable", values_to = "value")
  # }
  # 
  
  ########### SSDG ----
  df_trait_sp <- final_dataset_age %>% filter(Species == sp, Fstart == 1)

  # keep abundance by age variable on the side
  WEIGHT = df_trait_sp[which(df_trait_sp$variable=="abundanceByAge"),]

  # add the SSDgx/Nalpha
  # df_trait_t    <- estimateDiv(trait_data_full_i = df_trait_sp, var=paste0("SSDgx\nDiv", "Numerator_N_alpha"))
  # df_trait_sp   <- rbind(df_trait_sp, df_trait_t)
  
  # all \nDiv Numerator_N_alpha
  var_i <- unique(df_trait_sp$variable)
  for (v in var_i){
    if (v != "Numerator_N_alpha") {
      df_trait_t    <- estimateDiv(trait_data_full_i = df_trait_sp, var=paste0(v, "\nDivNumerator_N_alpha")) # div by population size
      df_trait_sp   <- rbind(df_trait_sp, df_trait_t)
    }
  }

  # produce and store plots
  for (relative_to_ref in c(FALSE)){
    proportional = FALSE
    if (relative_to_ref == "proportional") {
      relative_to_ref = FALSE
      proportional = TRUE
    }
    relative_to_ref = as.logical(relative_to_ref)
    plot_name_ch=c()
    plot_name <- paste0(as.character(relative_to_ref), as.character(proportional))
    if (plot_name=="TRUEFALSE") plot_name_ch <- c("Relative")
    if (plot_name=="FALSEFALSE") plot_name_ch <- c("Absolute")
    if (plot_name=="FALSETRUE") plot_name_ch <- c("Proportional")
    # for the plot
    plot_biom = ifelse(relative_to_ref, "rel_biomassplot", "biomassplot")
    plot_biom = ifelse(relative_to_ref, "biomassplot", "biomassplot")

    v= c("SSDgx\nDivNumerator_N_alpha")
    df_trait_sp_plot <- df_trait_sp[which(df_trait_sp$variable %in% v),]
    ssdg <-  add_SSDg(df=df_trait_sp_plot, weight=WEIGHT)
    ssdg$variable <- improve_vriable_names(ssdg$variable)
    df_trait_sp_plot$variable  <- improve_vriable_names(df_trait_sp_plot$variable)
    t <- improve_vriable_names(v)
    
    inset.tot.F.Nat[[plot_name]][[sp]][[paste(t, collapse="_")]] <- ggplot(ssdg, group = (variable)) +
      geom_boxplot(aes(y=SSDg, x=Fmax)) +
      labs(title=sp, y="")+ #, y="SSD\u0262"
      theme(strip.background = element_rect(fill="white"),
            panel.background = element_rect(fill = "white", colour = "grey50"),
            plot.title = element_text(hjust = 0.5),
            axis.text.x = element_text(angle = 60, hjust = 1),
            axis.text = element_text(colour = 1, size = 9),
            panel.border = element_rect(colour = "black", fill=NA),
            strip.text = element_text(size = 10))
    
    v= c("SSDgx", "SSDgx\nDivNumerator_N_alpha")
    df_trait_sp_plot <- df_trait_sp[which(df_trait_sp$variable %in% v),]
    ssdg <-  add_SSDg(df=df_trait_sp_plot, weight=WEIGHT)
    ssdg$variable <- improve_vriable_names(ssdg$variable)
    df_trait_sp_plot$variable  <- improve_vriable_names(df_trait_sp_plot$variable)
    t <- improve_vriable_names(v)

    inset.tot.F.Nat[[plot_name]][[sp]][[paste(t, collapse="_")]] <- ggplot(ssdg, group = (variable)) +
      geom_boxplot(aes(y=SSDg, x=Fmax)) +
      facet_wrap(~variable, scales="free", ncol=2) +
      labs(title=sp, y="")+ #, y="SSD\u0262"
      theme(strip.background = element_rect(fill="white"),
            panel.background = element_rect(fill = "white", colour = "grey50"),
            plot.title = element_text(hjust = 0.5),
            axis.text.x = element_text(angle = 60, hjust = 1),
            axis.text = element_text(colour = 1, size = 9),
            panel.border = element_rect(colour = "black", fill=NA),
            strip.text = element_text(size = 10))
    
    result <- plot_trait_vs_fmax(df_trait_sp_plot, trait_name = t, relative_to_ref = relative_to_ref, 
                                 Time1 = Time1, Time2=Time2, weight=WEIGHT, 
                                 spec = sp, age_dim = T, sp_dim = F, reference_input = NULL, 
                                 FORCElm = FALSE, nrow = 1, facet = T, position_legend = "right")
    plot_ssdg_age[[plot_name]][[sp]][[paste(t, collapse="_")]] <- result
    
    v= c("SSDgx", "SSDgx_F", "SSDgx_Mpred", "SSDgx_Mnat")
    df_trait_sp_plot <- df_trait_sp[which(df_trait_sp$variable %in% v),]
    df_trait_sp_plot$variable <- factor(df_trait_sp_plot$variable, levels = v)
    ssdg <-  add_SSDg(df=df_trait_sp_plot, weight=WEIGHT)
    ssdg$variable <- improve_vriable_names(ssdg$variable)
    t <- improve_vriable_names(v)
    ssdg$variable <- factor(ssdg$variable, levels = t)
    
    inset.tot.F.Nat[[plot_name]][[sp]][[paste(t, collapse="_")]] <- ggplot(ssdg, group = (variable)) +
      geom_boxplot(aes(y=SSDg, x=Fmax)) +
      facet_wrap(~variable, scales="free", ncol=2) +
      labs(title=sp, y="")+ #, y="SSD\u0262"
      theme(strip.background = element_rect(fill="white"),
            panel.background = element_rect(fill = "white", colour = "grey50"),
            plot.title = element_text(hjust = 0.5),
            axis.text.x = element_text(angle = 60, hjust = 1),
            axis.text = element_text(colour = 1, size = 9),
            panel.border = element_rect(colour = "black", fill=NA),
            strip.text = element_text(size = 10))    
    
    v= c("SSDgx\nDivNumerator_N_alpha", "SSDgx_F\nDivNumerator_N_alpha", 
         "SSDgx_Mpred\nDivNumerator_N_alpha", "SSDgx_Mnat\nDivNumerator_N_alpha")
    df_trait_sp_plot <- df_trait_sp[which(df_trait_sp$variable %in% v),]
    ssdg <-  add_SSDg(df=df_trait_sp_plot, weight=WEIGHT)
    ssdg$variable <- improve_vriable_names(ssdg$variable)
    t <- improve_vriable_names(v)
    ssdg$variable <- factor(ssdg$variable, levels = t)
    
    inset.tot.F.Nat[[plot_name]][[sp]][[paste(t, collapse="_")]] <- ggplot(ssdg, group = (variable)) +
      geom_boxplot(aes(y=SSDg, x=Fmax)) +
      facet_wrap(~variable, scales="free", ncol=2) +
      labs(title=sp, y="")+ #, y="SSD\u0262"
      theme(strip.background = element_rect(fill="white"),
            panel.background = element_rect(fill = "white", colour = "grey50"),
            plot.title = element_text(hjust = 0.5),
            axis.text.x = element_text(angle = 60, hjust = 1),
            axis.text = element_text(colour = 1, size = 9),
            panel.border = element_rect(colour = "black", fill=NA),
            strip.text = element_text(size = 10))   
    
    v= c("Deltasquare", "SSDgx_Mfor", "SSDgx_Mnat")
    df_trait_sp_plot <- df_trait_sp[which(df_trait_sp$variable %in% v),]
    ssdg <-  add_SSDg(df=df_trait_sp_plot, weight=WEIGHT)
    ssdg$variable <- improve_vriable_names(ssdg$variable)
    t <- improve_vriable_names(v)
    ssdg$variable <- factor(ssdg$variable, levels = t)
    
    inset.tot.F.Nat[[plot_name]][[sp]][[paste(t, collapse="_")]] <- ggplot(ssdg, group = (variable)) +
      geom_boxplot(aes(y=SSDg, x=Fmax)) +
      facet_wrap(~variable, scales="free", ncol=2) +
      labs(title=sp, y="")+ #, y="SSD\u0262"
      theme(strip.background = element_rect(fill="white"),
            panel.background = element_rect(fill = "white", colour = "grey50"),
            plot.title = element_text(hjust = 0.5),
            axis.text.x = element_text(angle = 60, hjust = 1),
            axis.text = element_text(colour = 1, size = 9),
            panel.border = element_rect(colour = "black", fill=NA),
            strip.text = element_text(size = 10))
    
  }

  
  # ######### Delta²xMortalitiesDead area under the curves -----------
  # 
  # df_trait_sp <- final_dataset_age %>% filter(Species == sp, Fstart == 1)
  # # remove cohort dimension
  # df_trait_sp <- df_trait_sp %>%
  #   mutate(age = as.numeric(age))%>%
  #   mutate(replicate = as.numeric(replicate)-1)%>% # to make the datasets match : unify the replicate range
  #   mutate(Fmax = as.numeric(Fmax))%>%
  #   mutate(Fstart = as.numeric(Fstart))
  # weight_df <- build_weight_df(df_trait_sp, weight_var = "abundanceByAge", key_cols = c("age", "Fmax", "Fstart", "Species", "replicate"))
  # df_trait_sp <- weighted_reduce(
  #   df = df_trait_sp,
  #   weight_df = weight_df,
  #   join_cols = c("Species", "replicate", "Fmax", "Fstart", "age"),
  #   group_cols = c("Species", "replicate", "Fstart", "Fmax", "variable", "age"),
  #   drop_age_zero = TRUE
  # )
  # 
  # # keep abundance by age variable on the side
  # WEIGHT = df_trait_sp[which(df_trait_sp$variable=="abundanceByAge"),]
  # 
  # # combine all traits: Mortality + Ne traits : both need to be pivot_wider to be combined
  # df_trait_sp <- pivot_wider(df_trait_sp, names_from = "variable", values_from = "value")
  # data_complete_c_i <- data_complete_c %>%
  #   filter(Species == sp)
  # data_complete_c_i <- pivot_wider(data_complete_c_i, names_from = "variable", values_from = "value")
  # df_trait_sp <- left_join(data_complete_c_i,df_trait_sp, by=c("age", "Fmax", "Fstart", "Species", "replicate"))
  # if (any(df_trait_sp$age == 0, na.rm=T)) df_trait_sp <- df_trait_sp[-which(df_trait_sp$age == 0), ]
  # df_trait_sp <- pivot_longer(data = df_trait_sp, cols = -which(colnames(df_trait_sp) %in% c("age", "Fmax", "Fstart", "Species", "replicate")),
  #                             names_to = "variable", values_to = "value")
  # df_trait_sp <- na.omit(df_trait_sp)
  # 
  # # total mortality rate 
  # df_trait_t    <- estimateSum(trait_data_full_i = df_trait_sp, var=("mean_of_F_per_age\nsummean_of_Mnat_per_age"))
  # df_trait_t$variable <- "total_M"
  # df_trait_sp   <- rbind(df_trait_sp, df_trait_t)
  # 
  # # add deaths estimates per mortaltiy source
  # var_i <- unique(df_trait_sp$variable)
  # for (v in var_i){
  #   if (any(grepl(pattern="mean_of_M", v) | grepl(pattern="mean_of_F", v) | grepl(pattern="total_M", v))) {
  #     df_trait_t    <- estimateMnumber(trait_data_full_i = df_trait_sp, var=paste0(v))
  #     df_trait_sp   <- rbind(df_trait_sp, df_trait_t)
  #   }
  # }
  # 
  # df_trait_t <- estimateSum(trait_data_full_i = df_trait_sp, var = paste0("Mnat_deaths\nsumF_deaths"))
  # df_trait_t$variable <- "all_deaths"
  # df_trait_sp   <- rbind(df_trait_sp, df_trait_t)
  # 
  # df_trait_t    <- estimateX(trait_data_full_i = df_trait_sp, var=paste0("survival_ratexabundanceByAge"))
  # #add delta square
  # df_trait_t    <- estimateSquare(trait_data_full_i = df_trait_sp, var=paste0("Deltasquare"))
  # df_trait_sp   <- rbind(df_trait_sp, df_trait_t)
  # 
  # #add delta square times dead number
  # var_i <- unique(df_trait_sp$variable)
  # for (v in var_i){
  #   if (any(grepl(pattern="deaths", v) | (v=="Dead"))) {
  #     cat(v)
  #     df_trait_t    <- estimateX(trait_data_full_i = df_trait_sp, var=paste0("Deltasquarex", v))
  #     df_trait_sp   <- rbind(df_trait_sp, df_trait_t)
  #   }
  # }
  # # divide all by Nalpha
  # var_i <- unique(df_trait_sp$variable)
  # for (v in var_i){
  #   if (v != "Numerator_N_alpha") {
  #     cat(v)
  #     df_trait_t    <- estimateDiv(trait_data_full_i = df_trait_sp, var=paste0(v, "\nDivNumerator_N_alpha")) # div by population size
  #     df_trait_sp   <- rbind(df_trait_sp, df_trait_t)
  #   }
  # }
  # 
  # # produce and store plots
  # for (relative_to_ref in c(FALSE)){
  #   proportional = FALSE
  #   if (relative_to_ref == "proportional") {
  #     relative_to_ref = FALSE
  #     proportional = TRUE
  #   }
  #   relative_to_ref = as.logical(relative_to_ref)
  #   plot_name_ch=c()
  #   plot_name <- paste0(as.character(relative_to_ref), as.character(proportional))
  #   if (plot_name=="TRUEFALSE") plot_name_ch <- c("Relative")
  #   if (plot_name=="FALSEFALSE") plot_name_ch <- c("Absolute")
  #   if (plot_name=="FALSETRUE") plot_name_ch <- c("Proportional")
  #   # for the plot
  #   plot_biom = ifelse(relative_to_ref, "rel_biomassplot", "biomassplot")
  #   plot_biom = ifelse(relative_to_ref, "biomassplot", "biomassplot")
  # 
  #   v= names(table(df_trait_sp$variable[which((grepl(df_trait_sp$variable, pattern="death") | grepl(df_trait_sp$variable, pattern="Dead")) &
  #                                 !grepl(df_trait_sp$variable, pattern="\nDiv") &
  #                                   grepl(df_trait_sp$variable, pattern="Deltasquare"))]))
  #   v <- v[c(grep("Dead", v), grep("nat", v), grep("F", v), grep("pred", v))]
  #   df_trait_sp_plot <- df_trait_sp[which(df_trait_sp$variable %in% c(v)),]#, "SSDgx", "SSDgx\nDivNumerator_N_alpha")),]
  #   df_trait_sp_plot$Fmax <- as.factor(df_trait_sp_plot$Fmax)
  #   df_trait_sp_plot$Fstart <- as.factor(df_trait_sp_plot$Fstart)
  #   WEIGHT$Fmax <- as.factor(WEIGHT$Fmax)
  #   WEIGHT$Fstart <- as.factor(WEIGHT$Fstart)
  #   WEIGHT$replicate <- as.numeric(WEIGHT$replicate)
  #   ssdg <-  add_SSDg(df=df_trait_sp_plot, weight=WEIGHT)
  #   ssdg$variable <- improve_vriable_names(ssdg$variable)
  #   t <- improve_vriable_names(v)
  # 
  #   inset.tot.F.Nat[[plot_name]][[sp]][[paste(t, collapse="_")]] <- ggplot(ssdg, group = (variable)) +
  #     geom_boxplot(aes(y=SSDg, x=Fmax)) +
  #     facet_wrap(~variable, scales="free", ncol=2, nrow=2) +
  #     labs(title=sp, y="")+ #, y="SSD\u0262"
  #     theme(strip.background = element_rect(fill="white"),
  #           panel.background = element_rect(fill = "white", colour = "grey50"),
  #           plot.title = element_text(hjust = 0.5),
  #           axis.text.x = element_text(angle = 60, hjust = 1),
  #           axis.text = element_text(colour = 1, size = 9),
  #           panel.border = element_rect(colour = "black", fill=NA),
  #           strip.text = element_text(size = 10))
  # 
  # 
  # 
  #   v= names(table(df_trait_sp$variable[which((grepl(df_trait_sp$variable, pattern="death") | grepl(df_trait_sp$variable, pattern="Dead")) &
  #                                               grepl(df_trait_sp$variable, pattern="\nDiv") &
  #                                               grepl(df_trait_sp$variable, pattern="Deltasquare"))]))
  #   v <- v[c(grep("Dead", v), grep("nat", v), grep("F", v), grep("pred", v))]
  #   df_trait_sp_plot <- df_trait_sp[which(df_trait_sp$variable %in% v),]
  #   df_trait_sp_plot$Fmax <- as.factor(df_trait_sp_plot$Fmax)
  #   df_trait_sp_plot$Fstart <- as.factor(df_trait_sp_plot$Fstart)
  #   WEIGHT$Fmax <- as.factor(WEIGHT$Fmax)
  #   WEIGHT$Fstart <- as.factor(WEIGHT$Fstart)
  #   WEIGHT$replicate <- as.numeric(WEIGHT$replicate)
  #   ssdg <-  add_SSDg(df=df_trait_sp_plot, weight=WEIGHT)
  #   ssdg$variable <- improve_vriable_names(ssdg$variable)
  #   t <- improve_vriable_names(v)
  # 
  #   inset.tot.F.Nat[[plot_name]][[sp]][[paste(t, collapse="_")]] <- ggplot(ssdg, group = (variable)) +
  #     geom_boxplot(aes(y=SSDg, x=Fmax)) +
  #     facet_wrap(~variable, scales="free", ncol=2, nrow=2) +
  #     labs(title=sp, y="")+ #, y="SSD\u0262"
  #     theme(strip.background = element_rect(fill="white"),
  #           panel.background = element_rect(fill = "white", colour = "grey50"),
  #           plot.title = element_text(hjust = 0.5),
  #           axis.text.x = element_text(angle = 60, hjust = 1),
  #           axis.text = element_text(colour = 1, size = 9),
  #           panel.border = element_rect(colour = "black", fill=NA),
  #           strip.text = element_text(size = 10))
  # }
}#sp
save(inset.tot.F.Nat, file=paste0(pathoutputRDATA, "/inset.tot.F.Nat.RData"))
load(file=paste0(pathoutputRDATA, "/inset.tot.F.Nat.RData"))
save(plot_ssdg_age, file=paste0(pathoutputRDATA, "/plot_ssdg_age.RData"))
load(file=paste0(pathoutputRDATA, "/plot_ssdg_age.RData"))



#### PLOT 3- Species functional traits ? Delta² Dead and Ne VS FMAX and AGE ------------------

Traitplot = list()
Traitplot_age = list()
for (sp in species){
  cat(sp)

  load(file=paste0(path, "/Complete_INIT2", sp, ".RData"))
  final_dataset_age = complete[[2]]
  final_dataset_age <- final_dataset_age %>%
    dplyr::mutate(dplyr::across(c(cohort, replicate), ~ stringr::str_remove(.x, "^(cohort_|rep_)")))
  final_dataset_age <- final_dataset_age[which(final_dataset_age$cohort %in% Time1:Time2),]
  rm(complete)

  # remove ages that are not included
  final_dataset_age$age <- as.numeric(final_dataset_age$age)
  ages_to_remove <- table(final_dataset_age[which(final_dataset_age$variable == "Dead" & final_dataset_age$value ==0),"age"])
  ages_to_remove <- ages_to_remove[which(ages_to_remove == max(ages_to_remove))]
  final_dataset_age <- final_dataset_age[-which(final_dataset_age$age %in% as.numeric(names(ages_to_remove))), ]

  # # avoid the estimation bias of very populations of very low size
  # if (any(((final_dataset_age[which(final_dataset_age$variable == "abundanceByAge"), "value"])<10))){
  #   final_dataset_age <- pivot_wider(final_dataset_age, names_from = "variable", values_from = "value")
  #   final_dataset_age <- final_dataset_age[-which((final_dataset_age$abundanceByAge)<10), ]
  #   final_dataset_age <- pivot_longer(final_dataset_age, cols = -c("age", "replicate", "Fmax", "Fstart", 'Species', "cohort"), names_to = "variable", values_to = "value")
  # }
  
  #### traits to plot ####
  traits_sp_to_plot = list(c("eggNumber_sc", "eggNumber_sc\nDivNumerator_N_alpha", "b_rate_sc", 
                             "kx_bar", "Delta", "Deltasquare", 
                             "Dead", "Dead\nDivNumerator_N_alpha", "Numerator_N_alpha"))
  df_trait_sp <- final_dataset_age %>% filter(Species == sp, Fstart == 1)
  rm(final_dataset_age)

  #add delta square
  df_trait_t    <- estimateSquare(trait_data_full_i = df_trait_sp, var=paste0("Deltasquare"))
  df_trait_sp   <- rbind(df_trait_sp, df_trait_t)
  
  #Devide by Nalpha
  for (v in c("eggNumber_sc\nDivNumerator_N_alpha",  "Dead\nDivNumerator_N_alpha")){
    df_trait_t    <- estimateDiv(trait_data_full_i = df_trait_sp, var=v)
    df_trait_sp   <- rbind(df_trait_sp, df_trait_t)
  }

  # weight
  weight_df <- build_weight_df(df = df_trait_sp, weight_var = "abundanceByAge",
                               key_cols = c("replicate", "Species", "Fmax", "Fstart", "age", "cohort"))

  
  
  
  for (trait_i in traits_sp_to_plot){

    cat(trait_i)

    # reduce dataset dimensions
    df_trait_sp_plot  <- df_trait_sp[which(df_trait_sp$variable %in% c(trait_i)), ]

    for (relative_to_ref in c(FALSE, TRUE, "proportional")){

      cat(relative_to_ref)

      proportional = FALSE
      if (relative_to_ref == "proportional") {
        relative_to_ref = FALSE
        proportional = TRUE
      }
      relative_to_ref = as.logical(relative_to_ref)
      # for the plot
      plot_biom = ifelse(relative_to_ref, "rel_biomassplot", "biomassplot")
      # plot_biom = ifelse(relative_to_ref, "biomassplot", "biomassplot")
      #
      # ref_data <- df_trait_sp_plot %>%
      #   filter(Fmax == 1 & Fstart == 1) %>%
      #   dplyr::select(variable, value, replicate, cohort, age) %>%
      #   dplyr::rename(ref_value = value)

      df_trait_sp_plot <- clean_values(df_trait_sp_plot)

      # conditions to stop
      if (nrow(df_trait_sp)==0) next
      if (all(is.na(df_trait_sp_plot$value) | df_trait_sp_plot$value == 0)) {
        next
      }

      var_ii                        <- improve_vriable_names(trait_i)
      df_trait_sp_plot$variable     <- improve_vriable_names(df_trait_sp_plot$variable)
      # ref_data$variable             <- improve_vriable_names(ref_data$variable)

      if (all(is.na(df_trait_sp_plot$value))) next
      result <- plot_trait_vs_fmax(df_trait_sp = df_trait_sp_plot, trait_name = var_ii,
                                   Time1, Time2, relative_to_ref = relative_to_ref, spec=sp,
                                   ncol=3, nrow=3,
                                   facet=T, 
                                   age_dim=T, FORCElm=F, #reference_input=ref_data,
                                   proportional = proportional, weight = weight_df, position_legend="right")
      Traitplot_age[[paste0(as.character(relative_to_ref), as.character(proportional))]][[sp]][[paste(var_ii, collapse="")]] <- result$plot

      result <- plot_trait_vs_fmax(df_trait_sp = df_trait_sp_plot, trait_name = var_ii,
                                   Time1, Time2, relative_to_ref = relative_to_ref, spec=sp,
                                   ncol=3, nrow=3,
                                   facet=T, 
                                   age_dim=F,  FORCElm=F, #reference_input=ref_data,
                                   proportional = proportional, weight = weight_df, position_legend="right")

      Traitplot[[paste0(as.character(relative_to_ref), as.character(proportional))]][[sp]][[paste(var_ii, collapse="")]]   <- result$plot
    }# relative to ref
  }# trait
}#sp


save(Traitplot, file=paste0(pathoutputRDATA, "/Traitplot.RData"))
load(file=paste0(pathoutputRDATA, "/Traitplot.RData"))
save(Traitplot_age, file=paste0(pathoutputRDATA, "/Traitplot_age.RData"))
load(file=paste0(pathoutputRDATA, "/Traitplot_age.RData"))



#### PLOT 4 : Mortality decomposition-----------------

combined_plot =list()
Time1=150
Time2=160
plot_Mortality_decomposition <- list()
for (sp in species){
  cat(sp)
  
  load(file=paste0(path, "/Complete_INIT2", sp, ".RData"))
  final_dataset_age = complete[[2]]
  final_dataset_age <- final_dataset_age %>%
    dplyr::mutate(dplyr::across(c(cohort, replicate), ~ stringr::str_remove(.x, "^(cohort_|rep_)")))
  final_dataset_age <- final_dataset_age[which(final_dataset_age$cohort %in% Time1:Time2),]
  rm(complete)
  
  # remove ages that are not included
  final_dataset_age$age <- as.numeric(final_dataset_age$age)
  ages_to_remove <- table(final_dataset_age[which(final_dataset_age$variable == "Dead" & final_dataset_age$value ==0),"age"])
  ages_to_remove <- ages_to_remove[which(ages_to_remove == max(ages_to_remove))]
  final_dataset_age <- final_dataset_age[-which(final_dataset_age$age %in% as.numeric(names(ages_to_remove))), ]

  # # avoid the estimation bias of very populations of very low size
  # if (any(((final_dataset_age[which(final_dataset_age$variable == "abundanceByAge"), "value"])<10))){
  #   final_dataset_age <- pivot_wider(final_dataset_age, names_from = "variable", values_from = "value")
  #   final_dataset_age <- final_dataset_age[-which((final_dataset_age$abundanceByAge)<10), ]
  #   final_dataset_age <- pivot_longer(final_dataset_age, cols = -c("age", "replicate", "Fmax", "Fstart", 'Species', "cohort"), names_to = "variable", values_to = "value")
  # }
  
  ### dataset
  df_trait_sp <- final_dataset_age %>% filter(Species == sp, Fstart == 1)
  
  # remove cohort dimension
  df_trait_sp <- df_trait_sp %>%
    mutate(age = as.numeric(age))%>%
    mutate(replicate = as.numeric(replicate)-1)%>% # to make the datasets match : unify the replicate range
    mutate(Fmax = as.numeric(Fmax))%>%
    mutate(Fstart = as.numeric(Fstart))
  weight_df <- build_weight_df(df_trait_sp, weight_var = "abundanceByAge", key_cols = c("age", "Fmax", "Fstart", "Species", "replicate"))
  df_trait_sp <- weighted_reduce(
    df = df_trait_sp,
    weight_df = weight_df,
    join_cols = c("Species", "replicate", "Fmax", "Fstart", "age"),
    group_cols = c("Species", "replicate", "Fstart", "Fmax", "variable", "age"),
    drop_age_zero = TRUE
  )
  
  # combine all traits: Mortality + Ne traits : both need to be pivot_wider to be combined
  df_trait_sp <- pivot_wider(df_trait_sp, names_from = "variable", values_from = "value")
  data_complete_c_i <- data_complete_c %>%
    filter(Species == sp)
  data_complete_c_i <- pivot_wider(data_complete_c_i, names_from = "variable", values_from = "value")
  df_trait_sp <- full_join(data_complete_c_i,df_trait_sp, by=c("age", "Fmax", "Fstart", "Species", "replicate"))
  if (any(df_trait_sp$age == 0)) df_trait_sp <- df_trait_sp[-which(df_trait_sp$age == 0), ]
  df_trait_sp <- pivot_longer(data = df_trait_sp, cols = -which(colnames(df_trait_sp) %in% c("age", "Fmax", "Fstart", "Species", "replicate")),
                              names_to = "variable", values_to = "value")
  df_trait_sp <- na.omit(df_trait_sp)
  
  # add deaths estimates per mortaltiy source
  df_trait_t    <- estimateSum(trait_data_full_i = df_trait_sp, var=paste0("mean_of_F_per_age", "\nsum", "mean_of_Mnat_per_age"))
  df_trait_t$variable    <- "total_M"
  df_trait_sp   <- rbind(df_trait_sp, df_trait_t)
  var_i <- unique(df_trait_sp$variable)
  for (v in var_i){
    if (any(grepl(pattern="mean_of_M", v) | grepl(pattern="mean_of_F", v) | grepl(pattern="sum_mortality_rates", v))) {
      df_trait_t    <- estimateMnumber(trait_data_full_i = df_trait_sp, var=paste0(v))
      df_trait_sp   <- rbind(df_trait_sp, df_trait_t)
    }
  }  
  # for (v in var_i){
  #   if (any(grepl(pattern="mean_of_M", v) | grepl(pattern="mean_of_F", v) | grepl(pattern="sum_mortality_rates", v))) {
  #     df_trait_t    <- estimateMnumber(trait_data_full_i = df_trait_sp, var=paste0(v), var2="survivorship_curve")
  #     df_trait_sp   <- rbind(df_trait_sp, df_trait_t)
  #   }
  # }
  df_trait_t    <- estimatecolSum(trait_data_full_i = df_trait_sp, var=paste0("colSumall_deaths"))
  df_trait_sp   <- rbind(df_trait_sp, df_trait_t)
  var_i <- unique(df_trait_sp$variable)
  for (v in var_i){
    if (any(grepl(pattern="death", v) & c(grepl(pattern="M", v)) | (grepl(pattern="F_death", v)) | (v=="all_deaths")| (v=="all_deaths_rate"))) {
      df_trait_t    <- estimateDiv(trait_data_full_i = df_trait_sp, var=paste0(v, "\nDiv", "colSumall_deaths"))
      df_trait_t2   <- estimatecolSum(trait_data_full_i = df_trait_t, var=paste0("colSum", (df_trait_t$variable)[1]))
      df_trait_sp   <- rbind(df_trait_sp, df_trait_t)
      df_trait_sp   <- rbind(df_trait_sp, df_trait_t2)
    }
  }
  
  #add delta square
  df_trait_t    <- estimateSquare(trait_data_full_i = df_trait_sp, var=paste0("Deltasquare"))
  df_trait_sp   <- rbind(df_trait_sp, df_trait_t)
  
  #add delta square
  for (v in var_i){
    if (any(grepl(pattern="deaths", v) | (v=="Dead"))) {
      df_trait_t    <- estimateX(trait_data_full_i = df_trait_sp, var=paste0("Deltasquarex", v))
      df_trait_sp   <- rbind(df_trait_sp, df_trait_t)
    }
  }
  
  #add sum by age of SSDgx and SSDix
  df_trait_t    <- estimatecolSum(trait_data_full_i = df_trait_sp, var=paste0("colSumSSDix"))
  df_trait_sp   <- rbind(df_trait_sp, df_trait_t)
  df_trait_t    <- estimatecolSum(trait_data_full_i = df_trait_sp, var=paste0("colSumSSDgx"))
  df_trait_sp   <- rbind(df_trait_sp, df_trait_t)
  df_trait_t    <- estimatecolSum(trait_data_full_i = df_trait_sp, var=paste0("colSumDeltasquare"))
  df_trait_sp   <- rbind(df_trait_sp, df_trait_t)
  
  var_i <- unique(df_trait_sp$variable)
  
  # keep abundance by age variable on the side
  WEIGHT = df_trait_sp[which(df_trait_sp$variable=="abundanceByAge"),]
  
  # compare all traits to their value /N1 to understand where is the increase coming from
  df_trait_t    <- estimatecolSum(trait_data_full_i = df_trait_sp, var=paste0("colSumDeltasquarexDead"))
  df_trait_sp   <- rbind(df_trait_sp, df_trait_t)
  
  for (v in var_i){
    if (v != "Numerator_N_alpha") {
      df_trait_t    <- estimateDiv(trait_data_full_i = df_trait_sp, var=paste0(v, "\nDivNumerator_N_alpha")) # div by population size
      df_trait_sp   <- rbind(df_trait_sp, df_trait_t)
    }
  }
  
  var_i <- unique(df_trait_sp$variable)
  
  for (relative_to_ref in c(FALSE)){
    proportional = FALSE
    if (relative_to_ref == "proportional") {
      relative_to_ref = FALSE
      proportional = TRUE
    }
    relative_to_ref = as.logical(relative_to_ref)
    plot_name_ch=c()
    plot_name <- paste0(as.character(relative_to_ref), as.character(proportional))
    if (plot_name=="TRUEFALSE") plot_name_ch <- c("Relative")
    if (plot_name=="FALSEFALSE") plot_name_ch <- c("Absolute")
    if (plot_name=="FALSETRUE") plot_name_ch <- c("Proportional")
    # for the plot
    plot_biom = ifelse(relative_to_ref, "rel_biomassplot", "biomassplot")
    plot_biom = ifelse(relative_to_ref, "biomassplot", "biomassplot")
    
    ref_data <- df_trait_sp %>%
      filter(Fstart == 1 & Fmax==1) %>%
      dplyr::select(variable, value, replicate, age) %>%
      dplyr::rename(ref_value = value)
    df_trait_sp_a <- df_trait_sp
    
    df_trait_sp_a <- clean_values(df_trait_sp_a)
    
    # conditions to stop
    if (nrow(df_trait_sp_a)==0) next
    if (all(is.na(df_trait_sp_a$value) | df_trait_sp_a$value == 0)) next
    
    # produce plot
    v <- names(table(df_trait_sp_a$variable))[c(which(grepl(names(table(df_trait_sp_a$variable)), pattern="Deltasquare") &
                                                        c(grepl(names(table(df_trait_sp_a$variable)), pattern="M") | grepl(names(table(df_trait_sp$variable)), pattern="F")) 
                                                      & !grepl(names(table(df_trait_sp_a$variable)), pattern="Div") & !grepl(names(table(df_trait_sp_a$variable)), pattern="deaths_rate") &
                                                        grepl(names(table(df_trait_sp_a$variable)), pattern="death")
                                                      & !grepl(names(table(df_trait_sp_a$variable)), pattern="colSum")),
                                                which(names(table(df_trait_sp_a$variable)) == "Deltasquarexall_deaths"))]
    v2 <- names(table(df_trait_sp_a$variable))[c(which(grepl(names(table(df_trait_sp_a$variable)), pattern="Deltasquare") & 
                                                         c(grepl(names(table(df_trait_sp_a$variable)), pattern="M") | grepl(names(table(df_trait_sp$variable)), pattern="F"))
                                                       & grepl(names(table(df_trait_sp_a$variable)), pattern="Div") & !grepl(names(table(df_trait_sp_a$variable)), pattern="deaths_rate")
                                                       & !grepl(names(table(df_trait_sp_a$variable)), pattern="colSum")),
                                                 which(names(table(df_trait_sp_a$variable)) == "Deltasquarexall_deaths\nDivNumerator_N_alpha"))]
    
    df_trait_sp_plot <- df_trait_sp_a[which(df_trait_sp_a$variable %in% v),]
    ssdg <-  add_SSDg(df=df_trait_sp_plot, weight=WEIGHT)
    df_trait_sp_plot <- left_join(df_trait_sp_plot, ssdg, by=c("Fmax", "Fstart", "variable", "Species"), relationship = "many-to-many")
    id_FNatTot <- which(df_trait_sp_plot$variable %in% c("Deltasquarexall_deaths", 
                                                         "DeltasquarexF_deaths", 
                                                         "DeltasquarexMnat_deaths"))
    df_trait_sp_plot$variable <- improve_vriable_names(df_trait_sp_plot$variable)
    t <- improve_vriable_names(v)
    result <- plot_trait_vs_age_simple(df_trait_sp=df_trait_sp_plot, sp_name=sp,
                                       var_names=t, Time1=150, Time2=160, relative_to_ref = relative_to_ref,
                                       scale_factor = 1, weight=WEIGHT, nrow=3, ncol=3,
                                       proportional = proportional, position_legend="auto", facet=T,
                                       chosen_title = paste0(plot_name_ch, " total variance of reproductive success per mortality rates sources per age,\n",
                                                             "averaged over years ",  Time1, " and ", Time2))
    if ("SSDg"%in%colnames(df_trait_sp_plot)){
      df_trait_sp_plot$Fmax <- as.factor(df_trait_sp_plot$Fmax)
      inset.plot <- ggplot(df_trait_sp_plot, group = (variable)) +
        geom_boxplot(aes(y=SSDg, x=Fmax)) +
        facet_wrap(~variable, scales="free") +
        theme(strip.background = element_rect(fill="white"), 
              panel.background = element_rect(fill = "white", colour = "grey50"),
              plot.title = element_text(hjust = 0.5), 
              axis.text.x = element_text(angle = 60, hjust = 1),
              axis.text = element_text(colour = 1, size = 9),
              panel.border = element_rect(colour = "black", fill=NA))
      
      inset.tot.F.Nat[[plot_name]][[sp]][[paste(t, collapse="_")]] <- ggplot(df_trait_sp_plot[id_FNatTot,], group = (variable)) +
        geom_boxplot(aes(y=SSDg, x=Fmax)) +
        facet_wrap(~variable, scales="free", ncol=1) +
        labs(title=sp)+ #, y="SSD\u0262"
        theme(strip.background = element_rect(fill="white"), 
              panel.background = element_rect(fill = "white", colour = "grey50"),
              plot.title = element_text(hjust = 0.5), 
              axis.text.x = element_text(angle = 60, hjust = 1),
              axis.text = element_text(colour = 1, size = 9),
              panel.border = element_rect(colour = "black", fill=NA))
      
      result <- result + theme(legend.position = c(0.53, 0))
      result <- ggdraw() +
        draw_plot(result) +
        draw_plot(inset.plot, x = 0.65, y = 0.0, width = 0.3, height = 0.32)
    }
    plot_Mortality_decomposition[[plot_name]][[sp]][[paste(t, collapse="_")]]<- result
    
    
    df_trait_sp_plot <- df_trait_sp_a[which(df_trait_sp_a$variable %in% v2),]
    ssdg <-  add_SSDg(df=df_trait_sp_plot, weight=WEIGHT)
    df_trait_sp_plot <- left_join(df_trait_sp_plot, ssdg, by=c("Fmax", "Fstart", "variable", "Species"), relationship = "many-to-many")
    id_FNatTot <- which(df_trait_sp_plot$variable %in% c("Deltasquarexall_deaths\nDivNumerator_N_alpha", 
                                                         "DeltasquarexF_deaths\nDivNumerator_N_alpha", 
                                                         "DeltasquarexMnat_deaths\nDivNumerator_N_alpha"))
    df_trait_sp_plot$variable <- improve_vriable_names(df_trait_sp_plot$variable)
    t <- improve_vriable_names(v2)
    result <- plot_trait_vs_age_simple(df_trait_sp=df_trait_sp_plot, sp_name=sp,
                                       var_names=t, Time1=150, Time2=160, relative_to_ref = relative_to_ref,
                                       scale_factor = 1, weight=WEIGHT, nrow=3, ncol=3,
                                       proportional = proportional, position_legend="auto", facet=T,
                                       chosen_title = paste0(plot_name_ch, " individual variance of reproductive success per mortality sources per age,\n",
                                                             "averaged over years ",  Time1, " and ", Time2))
    if ("SSDg"%in%colnames(df_trait_sp_plot)){
      df_trait_sp_plot$Fmax <- as.factor(df_trait_sp_plot$Fmax)
      inset.plot <- ggplot(df_trait_sp_plot, group = (variable)) +
        geom_boxplot(aes(y=SSDg, x=Fmax)) +
        facet_wrap(~variable, scales="free") +
        theme(strip.background = element_rect(fill="white"), 
              panel.background = element_rect(fill = "white", colour = "grey50"),
              plot.title = element_text(hjust = 0.5), 
              axis.text.x = element_text(angle = 60, hjust = 1),
              axis.text = element_text(colour = 1, size = 9),
              panel.border = element_rect(colour = "black", fill=NA))
      
      result <- result + theme(legend.position = c(0.53, 0))
      result <- ggdraw() +
        draw_plot(result) +
        draw_plot(inset.plot, x = 0.65, y = 0.0, width = 0.3, height = 0.32)
    }
    plot_Mortality_decomposition[[plot_name]][[sp]][[paste(t, collapse="_")]]<- result
  }#relative/proportional/none
}#sp

save(plot_Mortality_decomposition, file=paste0(pathoutputRDATA, "/plot_Mortality_decomposition.RData"))
plot_Mortality_decomposition <- get(load(file=paste0(pathoutputRDATA, "/plot_Mortality_decomposition.RData")))


#### PLOT : combine all and SAVE ----------------

## mortality decomposition -----------------
for (relative_to_ref in c(FALSE)){ #, FALSE, "proportional"
  p=list()
  proportional = FALSE
  if (relative_to_ref == "proportional") {
    relative_to_ref = FALSE
    proportional = TRUE
  }
  relative_to_ref = as.logical(relative_to_ref)
  plot_name_ch=c()
  plot_name <- paste0(as.character(relative_to_ref), as.character(proportional))
  if (plot_name=="TRUEFALSE") plot_name_ch <- c("relative")
  if (plot_name=="FALSEFALSE") plot_name_ch <- c("absolute")
  if (plot_name=="FALSETRUE") plot_name_ch <- c("proportional")
  
  windowsFonts(myCustomWindowsFontName=windowsFont("DejaVu Sans"))
  cairo_pdf(onefile=T, file = paste0(PATHOUT, "/Clean_for_paper/plotMortality_decomposition_", plot_name_ch, ".pdf"), width = 15, height = 10)
  for (sp in species){ p[[sp]] <-(plot_Mortality_decomposition[[plot_name]][[sp]][1])}
  print(p)
  dev.off()
  
  windowsFonts(myCustomWindowsFontName=windowsFont("DejaVu Sans"))
  cairo_pdf(onefile=T, file = paste0(PATHOUT, "/Clean_for_paper/plotMortality_decomposition_onN_", plot_name_ch, ".pdf"), width = 15, height = 10)
  for (sp in species){p[[sp]] <-(plot_Mortality_decomposition[[plot_name]][[sp]][2])}
  print(p)
  dev.off()
}


## explain why Ne ↑ --------------
for (relative_to_ref in c(TRUE, FALSE, "proportional")){
  if (!dir.exists(paste0(PATHOUT, "/Clean_for_paper"))) dir.create(paste0(PATHOUT, "/Clean_for_paper"))
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
  suffix_relative <- plot_name_ch
  
  windowsFonts(myCustomWindowsFontName=windowsFont("DejaVu Sans"))
  cairo_pdf(file = paste0(PATHOUT, "/Clean_for_paper/Ne_Vk_T_whyNeincrease", plot_name_ch, ".pdf"),
            width =12, height = 8, onefile=T)
  for (sp in species){
    plot_list <- c(Ne_plot_list[[plot_name]][[sp]][2:5])
    # plot_list <- lapply(plot_list, function(x) x + theme_minimal(base_size = 11) + theme(legend.position = "bottom"))
    combined_plot <- ggarrange(plotlist = plot_list,
                               ncol = 2,
                               nrow = 2)
    combined_plot <- annotate_figure(combined_plot, top = text_grob(paste0("Explaining ", sp, " ", plot_name_ch, " variation of Ne/N"),
                                                                    color = "grey60", face = "bold", size = 14))
    print(combined_plot)
  }
  dev.off()  
}


## Functional trait groups -----------
for (relative_to_ref in c(TRUE, FALSE, "proportional")){
  if (!dir.exists(paste0(PATHOUT, "/Clean_for_paper"))) dir.create(paste0(PATHOUT, "/Clean_for_paper"))
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
  suffix_relative <- plot_name_ch
  
  windowsFonts(myCustomWindowsFontName=windowsFont("DejaVu Sans"))
  cairo_pdf(file = paste0(PATHOUT, "/Clean_for_paper/Traitplot", plot_name_ch, ".pdf"),
            width =10, height = 8, onefile=T)
  for (sp in species){
    plot_list <- c(Traitplot[[plot_name]][[sp]])                          
    print(plot_list)
  }
  dev.off()  
  
  
  windowsFonts(myCustomWindowsFontName=windowsFont("DejaVu Sans"))
  cairo_pdf(file = paste0(PATHOUT, "/Clean_for_paper/Traitplot_age", plot_name_ch, ".pdf"),
            width =10, height = 8, onefile=T)
  for (sp in species){
    plot_list <- c(Traitplot_age[[plot_name]][[sp]])
    print(plot_list)
  }
  dev.off()  
  
}

## Ne and T/Vk------------------
for (relative_to_ref in c(TRUE, FALSE, "proportional")){
  if (!dir.exists(paste0(PATHOUT, "/Clean_for_paper"))) dir.create(paste0(PATHOUT, "/Clean_for_paper"))
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
  suffix_relative <- plot_name_ch
  
  windowsFonts(myCustomWindowsFontName=windowsFont("DejaVu Sans"))
  cairo_pdf(file = paste0(PATHOUT, "/Clean_for_paper/Ne_Vk_T_", plot_name_ch, ".pdf"),
            width = 747/72, height = 2*450/72, onefile=T)
  for (sp in species){
    plot_list <- c(Ne_plot_list[[plot_name]][[sp]][1:2])
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
  plot_list <- c(Ne_plot_list[[plot_name]][["Whiting"]][2], 
                 Ne_plot_list[[plot_name]][["AtlanticMackerel"]][2],
                 Ne_plot_list[[plot_name]][["Sandeel"]][2],
                 Ne_plot_list[[plot_name]][["HorseMackerel"]][2])
  # plot_list <- lapply(plot_list, function(x) x + theme_minimal(base_size = 11) + theme(legend.position = "bottom"))
  combined_plot <- ggarrange(plotlist = plot_list,
                             ncol = 2,
                             nrow = 2, common.legend = T)
  combined_plot <- annotate_figure(combined_plot)
  print(combined_plot)
  dev.off()
  
  
  windowsFonts(myCustomWindowsFontName=windowsFont("DejaVu Sans"))
  cairo_pdf(file = paste0(PATHOUT, "/Clean_for_paper/Vk_final", plot_name_ch, ".pdf"),
            width = 747/72, height = 747/72, onefile=T)
  plot_list <- c(Ne_plot_list[[plot_name]][["Whiting"]][1], 
                 Ne_plot_list[[plot_name]][["AtlanticMackerel"]][1],
                 Ne_plot_list[[plot_name]][["Sandeel"]][1],
                 Ne_plot_list[[plot_name]][["HorseMackerel"]][1])
  # plot_list <- lapply(plot_list, function(x) x + theme_minimal(base_size = 11) + theme(legend.position = "bottom"))
  combined_plot <- ggarrange(plotlist = plot_list,
                             ncol = 2,
                             nrow = 2, common.legend = T)
  combined_plot <- annotate_figure(combined_plot)
  print(combined_plot)
  dev.off() 
  
    # pdf(file = paste0(PATHOUT, "/Clean_for_paper/Deltasquare_Death", suffix_relative, ".pdf"),
  #     width = 10, height = 10)
  # for (sp in species){
  #   plot_list <- c(SSDplot[[plot_name]][[sp]][1],
  #                  deaths_deltasquare_byage[[plot_name]][[sp]][c(1:2)])
  #   # plot_list <- lapply(plot_list, function(x) x + theme_minimal(base_size = 11) + theme(legend.position = "bottom"))
  #   combined_plot <- ggarrange(plotlist = plot_list,
  #                              ncol = 2,
  #                              nrow = 2)
  #   combined_plot <- annotate_figure(combined_plot, top = text_grob(paste0("Explaining ", sp, " ", suffix_relative, " variation of Ne/N"),
  #                                                                   color = "grey60", face = "bold", size = 14))
  #   print(combined_plot)
  # }
  # dev.off()

  
  # pdf(file = paste0(PATHOUT, "/Clean_for_paper/Mortality", suffix_relative, ".pdf"),
  #     width = 8, height = 10)
  # for (sp in species){
  #   plot_list <- c(mortality_plot[[plot_name]][[sp]])
  #   # plot_list <- lapply(plot_list, function(x) x + theme_minimal(base_size = 11) + theme(legend.position = "bottom"))
  #   combined_plot <- ggarrange(plotlist = plot_list,
  #                              ncol = 2,
  #                              nrow = 2)
  #   combined_plot <- annotate_figure(combined_plot, top = text_grob(paste0("Explaining ", sp, " ", suffix_relative, " variation of Ne/N"),
  #                                                                   color = "grey60", face = "bold", size = 14))
  #   print(combined_plot)
  # }
  # dev.off()
  # 
  # pdf(file = paste0(PATHOUT, "/Clean_for_paper/Survival_Biomass", suffix_relative, ".pdf"),
  #     width = 10, height = 10)
  # for (sp in species){
  #   plot_list <- c(deaths_deltasquare_byage[[plot_name]][[sp]][3],
  #                  list(get(plot_biom)[[sp]]))
  #   # plot_list <- lapply(plot_list, function(x) x + theme_minimal(base_size = 11) + theme(legend.position = "bottom"))
  #   combined_plot <- ggarrange(plotlist = plot_list,
  #                              ncol = 1,
  #                              nrow = 2)
  #   combined_plot <- annotate_figure(combined_plot, top = text_grob(paste0("Explaining ", sp," ", suffix_relative, " variation of Ne/N"),
  #                                                                   color = "grey60", face = "bold", size = 14))
  #   print(combined_plot)
  # }
  # dev.off()
}


### SSD ----
for (relative_to_ref in c(FALSE)){
  if (!dir.exists(paste0(PATHOUT, "/Clean_for_paper"))) dir.create(paste0(PATHOUT, "/Clean_for_paper"))
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
  suffix_relative <- plot_name_ch

  
  windowsFonts(myCustomWindowsFontName=windowsFont("DejaVu Sans"))
  cairo_pdf(file = paste0(PATHOUT, "/Clean_for_paper/SSDg_tot_", plot_name_ch, ".pdf"),
            width = 6, height =6, onefile=T)
  for (sp in species){
    plot_list <- c(inset.tot.F.Nat[[plot_name]][[sp]])
    # plot_list <- lapply(plot_list, function(x) x + theme_minimal(base_size = 11) + theme(legend.position = "bottom"))
    print(plot_list)
  }
  dev.off()
  
  
  windowsFonts(myCustomWindowsFontName=windowsFont("DejaVu Sans"))
  cairo_pdf(file = paste0(PATHOUT, "/Clean_for_paper/SSDgonN_final_", plot_name_ch, ".pdf"),
            width = 10, height =6, onefile=T)
  plot_list <- c(inset.tot.F.Nat[[plot_name]][["Whiting"]][1], 
                 inset.tot.F.Nat[[plot_name]][["AtlanticMackerel"]][1],
                 inset.tot.F.Nat[[plot_name]][["Sandeel"]][1],
                 inset.tot.F.Nat[[plot_name]][["EuropeanPlaice"]][1])
  # plot_list <- lapply(plot_list, function(x) x + theme_minimal(base_size = 11) + theme(legend.position = "bottom"))
  combined_plot <- ggarrange(plotlist = plot_list,
                             ncol = 2,
                             nrow = 2, common.legend = T)
  combined_plot <- annotate_figure(combined_plot)
  print(combined_plot)
  dev.off()  
  
  windowsFonts(myCustomWindowsFontName=windowsFont("DejaVu Sans"))
  cairo_pdf(file = paste0(PATHOUT, "/Clean_for_paper/SSDgSSDgonN_final_", plot_name_ch, ".pdf"),
            width = 10, height =6, onefile=T)
  plot_list <- c(inset.tot.F.Nat[[plot_name]][["Whiting"]][2], 
                 inset.tot.F.Nat[[plot_name]][["AtlanticMackerel"]][2],
                 inset.tot.F.Nat[[plot_name]][["Sandeel"]][2],
                 inset.tot.F.Nat[[plot_name]][["EuropeanPlaice"]][2])
  # plot_list <- lapply(plot_list, function(x) x + theme_minimal(base_size = 11) + theme(legend.position = "bottom"))
  combined_plot <- ggarrange(plotlist = plot_list,
                             ncol = 2,
                             nrow = 2, common.legend = T)
  combined_plot <- annotate_figure(combined_plot)
  print(combined_plot)
  dev.off()  
  
  windowsFonts(myCustomWindowsFontName=windowsFont("DejaVu Sans"))
  cairo_pdf(file = paste0(PATHOUT, "/Clean_for_paper/DeltasquareDead_final_", plot_name_ch, ".pdf"),
            width = 14, height = 10, onefile=T)
  plot_list <- c(inset.tot.F.Nat[[plot_name]][["Whiting"]][3], 
                 inset.tot.F.Nat[[plot_name]][["AtlanticMackerel"]][3],
                 inset.tot.F.Nat[[plot_name]][["Sandeel"]][3],
                 inset.tot.F.Nat[[plot_name]][["EuropeanPlaice"]][3])
  # plot_list <- lapply(plot_list, function(x) x + theme_minimal(base_size = 11) + theme(legend.position = "bottom"))
  combined_plot <- ggarrange(plotlist = plot_list,
                             ncol = 2,
                             nrow = 2, common.legend = T)
  combined_plot <- annotate_figure(combined_plot)
  print(combined_plot)
  dev.off()
  
  windowsFonts(myCustomWindowsFontName=windowsFont("DejaVu Sans"))
  cairo_pdf(file = paste0(PATHOUT, "/Clean_for_paper/DeltasquareDead_finalonN_", plot_name_ch, ".pdf"),
            width = 14, height = 10, onefile=T)
  plot_list <- c(inset.tot.F.Nat[[plot_name]][["Whiting"]][4], 
                 inset.tot.F.Nat[[plot_name]][["AtlanticMackerel"]][4],
                 inset.tot.F.Nat[[plot_name]][["Sandeel"]][4],
                 inset.tot.F.Nat[[plot_name]][["EuropeanPlaice"]][4])
  # plot_list <- lapply(plot_list, function(x) x + theme_minimal(base_size = 11) + theme(legend.position = "bottom"))
  combined_plot <- ggarrange(plotlist = plot_list,
                             ncol = 2,
                             nrow = 2, common.legend = T)
  combined_plot <- annotate_figure(combined_plot)
  print(combined_plot)
  dev.off()  
  
  windowsFonts(myCustomWindowsFontName=windowsFont("DejaVu Sans"))
  cairo_pdf(file = paste0(PATHOUT, "/Clean_for_paper/DeltasquareDead_final_increasing2_", plot_name_ch, ".pdf"),
            width = 14, height = 10, onefile=T)
  plot_list <- c(inset.tot.F.Nat[[plot_name]][["AtlanticMackerel"]][3], 
                 inset.tot.F.Nat[[plot_name]][["GreyGurnard"]][3],
                 inset.tot.F.Nat[[plot_name]][["CommonDab"]][3],
                 inset.tot.F.Nat[[plot_name]][["EuropeanPlaice"]][3])
  # plot_list <- lapply(plot_list, function(x) x + theme_minimal(base_size = 11) + theme(legend.position = "bottom"))
  combined_plot <- ggarrange(plotlist = plot_list,
                             ncol = 2,
                             nrow = 2, common.legend = T)
  combined_plot <- annotate_figure(combined_plot)
  print(combined_plot)
  dev.off()
  
  windowsFonts(myCustomWindowsFontName=windowsFont("DejaVu Sans"))
  cairo_pdf(file = paste0(PATHOUT, "/Clean_for_paper/DeltasquareDead_finalonN_increasing2_", plot_name_ch, ".pdf"),
            width = 14, height = 10, onefile=T)
  plot_list <- c(inset.tot.F.Nat[[plot_name]][["AtlanticMackerel"]][4], 
                 inset.tot.F.Nat[[plot_name]][["GreyGurnard"]][4],
                 inset.tot.F.Nat[[plot_name]][["CommonDab"]][4],
                 inset.tot.F.Nat[[plot_name]][["EuropeanPlaice"]][4])
  # plot_list <- lapply(plot_list, function(x) x + theme_minimal(base_size = 11) + theme(legend.position = "bottom"))
  combined_plot <- ggarrange(plotlist = plot_list,
                             ncol = 2,
                             nrow = 2, common.legend = T)
  combined_plot <- annotate_figure(combined_plot)
  print(combined_plot)
  dev.off()
  
  
  windowsFonts(myCustomWindowsFontName=windowsFont("DejaVu Sans"))
  par(mfrow=c(2,2))
  cairo_pdf(file = paste0(PATHOUT, "/Clean_for_paper/plot_ssdg_age", plot_name_ch, ".pdf"),
            width = 14, height = 7, onefile=T)
  for (sp in c("Whiting", "AtlanticMackerel", "Sandeel", "HorseMackerel")){
    plot_list <- c(plot_ssdg_age[[plot_name]][[sp]][1])
    print(plot_list)
  }
  dev.off()
}



# _  ------------------
# NOT USED ↓ ---------------


#### PLOT +4: mortalities and delta² and survival ----------------
# 
# combined_plot =list()
# Time1=150
# Time2=160
# plot_Deathsdelta_age = list()
# plot_Abund_age = list()
# plots_by_age_delta= list()
# inset.tot.F.Nat = list()
# plot_Deathsdelta_age2 = list()
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
#   
#   
#   ## PLOT DELTA²xDEAD_i/N and compare the mortality source of variation in the results -----------------
#   traits_sp_to_plot = list(c("Deltasquare","Dead"),
#                            c("DeltasquarexDead", "DeltasquarexDead\nDivNumerator_N_alpha"), 
#                            c("survival_rate"))
#   traits_sp_to_plot <- list(c("Deltasquare","Dead"),c("Deltasquare\nDivNumerator_N_alpha","Dead\nDivNumerator_N_alpha"),
#                            c("DeltasquarexDead", "DeltasquarexDead\nDivNumerator_N_alpha"), 
#                            c("survival_rate"), c("kx_point", "kx_bar"))
#   df_trait_sp <- final_dataset_age %>% filter(Species == sp, Fstart == 1)
#   
#   # avoid the estimation bias of very populations of very low size
#   if (any(((df_trait_sp[which(df_trait_sp$variable == "abundanceByAge"), "value"])<10))){
#     df_trait_sp <- pivot_wider(df_trait_sp, names_from = "variable", values_from = "value")
#     df_trait_sp <- df_trait_sp[-which((df_trait_sp$abundanceByAge)<10), ]
#     df_trait_sp <- pivot_longer(df_trait_sp, cols = -c("age", "replicate", "Fmax", "Fstart", 'Species', "cohort"), names_to = "variable", values_to = "value")
#   }
#   
#   # remove cohort dimension
#   df_trait_sp <- df_trait_sp %>%
#     mutate(age = as.numeric(age))%>%
#     mutate(replicate = as.numeric(replicate)-1)%>% # to make the datasets match : unify the replicate range
#     mutate(Fmax = as.numeric(Fmax))%>%
#     mutate(Fstart = as.numeric(Fstart))
#   weight_df <- build_weight_df(df_trait_sp, weight_var = "abundanceByAge", key_cols = c("age", "Fmax", "Fstart", "Species", "replicate"))
#   df_trait_sp <- weighted_reduce(
#     df = df_trait_sp,
#     weight_df = weight_df,
#     join_cols = c("Species", "replicate", "Fmax", "Fstart", "age"),
#     group_cols = c("Species", "replicate", "Fstart", "Fmax", "variable", "age"),
#     drop_age_zero = TRUE
#   )
#   
#   # combine all traits: Mortality + Ne traits : both need to be pivot_wider to be combined
#   df_trait_sp <- pivot_wider(df_trait_sp, names_from = "variable", values_from = "value")
#   data_complete_c_i <- data_complete_c %>%
#     filter(Species == sp)
#   data_complete_c_i <- pivot_wider(data_complete_c_i, names_from = "variable", values_from = "value")
#   df_trait_sp <- full_join(data_complete_c_i,df_trait_sp, by=c("age", "Fmax", "Fstart", "Species", "replicate"))
#   if (any(df_trait_sp$age == 0)) df_trait_sp <- df_trait_sp[-which(df_trait_sp$age == 0), ]
#   df_trait_sp <- pivot_longer(data = df_trait_sp, cols = -which(colnames(df_trait_sp) %in% c("age", "Fmax", "Fstart", "Species", "replicate")),
#                               names_to = "variable", values_to = "value")
#   df_trait_sp <- na.omit(df_trait_sp)
#   
#   # add deaths estimates per mortaltiy source
#   df_trait_t    <- estimateSum(trait_data_full_i = df_trait_sp, var=paste0("mean_of_F_per_age", "\nsum", "mean_of_Mnat_per_age"))
#   df_trait_t$variable    <- "sum_mortality_rates"
#   df_trait_sp   <- rbind(df_trait_sp, df_trait_t)
#   var_i <- unique(df_trait_sp$variable)
#   for (v in var_i){
#     if (any(grepl(pattern="mean_of_M", v) | grepl(pattern="mean_of_F", v) | grepl(pattern="sum_mortality_rates", v))) {
#       df_trait_t    <- estimateMnumber(trait_data_full_i = df_trait_sp, var=paste0(v))
#       df_trait_sp   <- rbind(df_trait_sp, df_trait_t)
#     }
#   }  
#   # for (v in var_i){
#   #   if (any(grepl(pattern="mean_of_M", v) | grepl(pattern="mean_of_F", v) | grepl(pattern="sum_mortality_rates", v))) {
#   #     df_trait_t    <- estimateMnumber(trait_data_full_i = df_trait_sp, var=paste0(v), var2="survivorship_curve")
#   #     df_trait_sp   <- rbind(df_trait_sp, df_trait_t)
#   #   }
#   # }
#   df_trait_t    <- estimatecolSum(trait_data_full_i = df_trait_sp, var=paste0("colSumall_deaths"))
#   df_trait_sp   <- rbind(df_trait_sp, df_trait_t)
#   var_i <- unique(df_trait_sp$variable)
#   for (v in var_i){
#     if (any(grepl(pattern="death", v) & c(grepl(pattern="M", v)) | (grepl(pattern="F_death", v)) | (v=="all_deaths")| (v=="all_deaths_rate"))) {
#       df_trait_t    <- estimateDiv(trait_data_full_i = df_trait_sp, var=paste0(v, "\nDiv", "colSumall_deaths"))
#       df_trait_t2   <- estimatecolSum(trait_data_full_i = df_trait_t, var=paste0("colSum", (df_trait_t$variable)[1]))
#       df_trait_sp   <- rbind(df_trait_sp, df_trait_t)
#       df_trait_sp   <- rbind(df_trait_sp, df_trait_t2)
#     }
#   }
#   
#   #add delta square
#   df_trait_t    <- estimateSquare(trait_data_full_i = df_trait_sp, var=paste0("Deltasquare"))
#   df_trait_sp   <- rbind(df_trait_sp, df_trait_t)
#   
#   #add delta square
#   for (v in var_i){
#     if (any(grepl(pattern="deaths", v) | (v=="Dead"))) {
#       df_trait_t    <- estimateX(trait_data_full_i = df_trait_sp, var=paste0("Deltasquarex", v))
#       df_trait_sp   <- rbind(df_trait_sp, df_trait_t)
#     }
#   }
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
#   WEIGHT = df_trait_sp[which(df_trait_sp$variable=="abundanceByAge"),]
#   
#   # compare all traits to their value /N1 to understand where is the increase coming from
#   df_trait_t    <- estimatecolSum(trait_data_full_i = df_trait_sp, var=paste0("colSumDeltasquarexDead"))
#   df_trait_sp   <- rbind(df_trait_sp, df_trait_t)
#   
#   for (v in var_i){
#     if (v != "Numerator_N_alpha") {
#       df_trait_t    <- estimateDiv(trait_data_full_i = df_trait_sp, var=paste0(v, "\nDivNumerator_N_alpha")) # div by population size
#       df_trait_sp   <- rbind(df_trait_sp, df_trait_t)
#     }
#   }
#   
#   var_i <- unique(df_trait_sp$variable)
#   
#   for (relative_to_ref in c(FALSE)){
#     proportional = FALSE
#     if (relative_to_ref == "proportional") {
#       relative_to_ref = FALSE
#       proportional = TRUE
#     }
#     relative_to_ref = as.logical(relative_to_ref)
#     plot_name_ch=c()
#     plot_name <- paste0(as.character(relative_to_ref), as.character(proportional))
#     if (plot_name=="TRUEFALSE") plot_name_ch <- c("Relative")
#     if (plot_name=="FALSEFALSE") plot_name_ch <- c("Absolute")
#     if (plot_name=="FALSETRUE") plot_name_ch <- c("Proportional")
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
#     if (all(is.na(df_trait_sp_a$value) | df_trait_sp_a$value == 0)) next
#     
#   # produce plot
#     v <- names(table(df_trait_sp_a$variable))[c(which(grepl(names(table(df_trait_sp_a$variable)), pattern="Deltasquare") &
#                                                         c(grepl(names(table(df_trait_sp_a$variable)), pattern="M") | grepl(names(table(df_trait_sp$variable)), pattern="F")) 
#                                                       & !grepl(names(table(df_trait_sp_a$variable)), pattern="Div") & !grepl(names(table(df_trait_sp_a$variable)), pattern="deaths_rate") &
#                                                         grepl(names(table(df_trait_sp_a$variable)), pattern="death")
#                                                       & !grepl(names(table(df_trait_sp_a$variable)), pattern="colSum")),
#                                                 which(names(table(df_trait_sp_a$variable)) == "Deltasquarexall_deaths"))]
#     v2 <- names(table(df_trait_sp_a$variable))[c(which(grepl(names(table(df_trait_sp_a$variable)), pattern="Deltasquare") & 
#                                                          c(grepl(names(table(df_trait_sp_a$variable)), pattern="M") | grepl(names(table(df_trait_sp$variable)), pattern="F"))
#                                                        & grepl(names(table(df_trait_sp_a$variable)), pattern="Div") & !grepl(names(table(df_trait_sp_a$variable)), pattern="deaths_rate")
#                                                        & !grepl(names(table(df_trait_sp_a$variable)), pattern="colSum")),
#                                                  which(names(table(df_trait_sp_a$variable)) == "Deltasquarexall_deaths\nDivNumerator_N_alpha"))]
#     
#     df_trait_sp_plot <- df_trait_sp_a[which(df_trait_sp_a$variable %in% v),]
#     ssdg <-  add_SSDg(df=df_trait_sp_plot, weight=WEIGHT)
#     df_trait_sp_plot <- left_join(df_trait_sp_plot, ssdg, by=c("Fmax", "Fstart", "variable", "Species"), relationship = "many-to-many")
#     id_FNatTot <- which(df_trait_sp_plot$variable %in% c("Deltasquarexall_deaths", 
#                                                          "DeltasquarexF_deaths", 
#                                                          "DeltasquarexMnat_deaths"))
#     df_trait_sp_plot$variable <- improve_vriable_names(df_trait_sp_plot$variable)
#     t <- improve_vriable_names(v)
#     result <- plot_trait_vs_age_simple(df_trait_sp=df_trait_sp_plot, sp_name=sp,
#                                        var_names=t, Time1=150, Time2=160, relative_to_ref = relative_to_ref,
#                                        scale_factor = 1, weight=WEIGHT, nrow=3, ncol=3,
#                                        proportional = proportional, position_legend="auto", facet=T,
#                                        chosen_title = paste0(plot_name_ch, " total variance of reproductive success per mortality rates sources per age,\n",
#                                                              "averaged over years ",  Time1, " and ", Time2))
#     if ("SSDg"%in%colnames(df_trait_sp_plot)){
#       df_trait_sp_plot$Fmax <- as.factor(df_trait_sp_plot$Fmax)
#       inset.plot <- ggplot(df_trait_sp_plot, group = (variable)) +
#         geom_boxplot(aes(y=SSDg, x=Fmax)) +
#         facet_wrap(~variable, scales="free") +
#         theme(strip.background = element_rect(fill="white"), 
#               panel.background = element_rect(fill = "white", colour = "grey50"),
#               plot.title = element_text(hjust = 0.5), 
#               axis.text.x = element_text(angle = 60, hjust = 1),
#               axis.text = element_text(colour = 1, size = 9),
#               panel.border = element_rect(colour = "black", fill=NA))
#       
#       inset.tot.F.Nat[[plot_name]][[sp]][[paste(t, collapse="_")]] <- ggplot(df_trait_sp_plot[id_FNatTot,], group = (variable)) +
#         geom_boxplot(aes(y=SSDg, x=Fmax)) +
#         facet_wrap(~variable, scales="free", ncol=1) +
#         labs(title=sp)+ #, y="SSD\u0262"
#         theme(strip.background = element_rect(fill="white"), 
#               panel.background = element_rect(fill = "white", colour = "grey50"),
#               plot.title = element_text(hjust = 0.5), 
#               axis.text.x = element_text(angle = 60, hjust = 1),
#               axis.text = element_text(colour = 1, size = 9),
#               panel.border = element_rect(colour = "black", fill=NA))
#       
#       result <- result + theme(legend.position = c(0.53, 0))
#       result <- ggdraw() +
#         draw_plot(result) +
#         draw_plot(inset.plot, x = 0.65, y = 0.0, width = 0.3, height = 0.32)
#     }
#     plot_Deathsdelta_age[[plot_name]][[sp]][[paste(t, collapse="_")]]<- result
#     
#     
#     df_trait_sp_plot <- df_trait_sp_a[which(df_trait_sp_a$variable %in% v2),]
#     ssdg <-  add_SSDg(df=df_trait_sp_plot, weight=WEIGHT)
#     df_trait_sp_plot <- left_join(df_trait_sp_plot, ssdg, by=c("Fmax", "Fstart", "variable", "Species"), relationship = "many-to-many")
#     id_FNatTot <- which(df_trait_sp_plot$variable %in% c("Deltasquarexall_deaths\nDivNumerator_N_alpha", 
#                                                          "DeltasquarexF_deaths\nDivNumerator_N_alpha", 
#                                                          "DeltasquarexMnat_deaths\nDivNumerator_N_alpha"))
#     df_trait_sp_plot$variable <- improve_vriable_names(df_trait_sp_plot$variable)
#     t <- improve_vriable_names(v2)
#     result <- plot_trait_vs_age_simple(df_trait_sp=df_trait_sp_plot, sp_name=sp,
#                                        var_names=t, Time1=150, Time2=160, relative_to_ref = relative_to_ref,
#                                        scale_factor = 1, weight=WEIGHT, nrow=3, ncol=3,
#                                        proportional = proportional, position_legend="auto", facet=T,
#                                        chosen_title = paste0(plot_name_ch, " individual variance of reproductive success per mortality sources per age,\n",
#                                                              "averaged over years ",  Time1, " and ", Time2))
#     if ("SSDg"%in%colnames(df_trait_sp_plot)){
#       df_trait_sp_plot$Fmax <- as.factor(df_trait_sp_plot$Fmax)
#       inset.plot <- ggplot(df_trait_sp_plot, group = (variable)) +
#         geom_boxplot(aes(y=SSDg, x=Fmax)) +
#         facet_wrap(~variable, scales="free") +
#         theme(strip.background = element_rect(fill="white"), 
#               panel.background = element_rect(fill = "white", colour = "grey50"),
#               plot.title = element_text(hjust = 0.5), 
#               axis.text.x = element_text(angle = 60, hjust = 1),
#               axis.text = element_text(colour = 1, size = 9),
#               panel.border = element_rect(colour = "black", fill=NA))
#       
#       inset.tot.F.Nat[[plot_name]][[sp]][[paste(t, collapse="_")]] <- ggplot(df_trait_sp_plot[id_FNatTot,], group = (variable)) +
#         geom_boxplot(aes(y=SSDg, x=Fmax)) +
#         facet_wrap(~variable, scales="free", ncol=1) +
#         labs(title=sp)+ #, y="SSD\u0262"
#         theme(strip.background = element_rect(fill="white"), 
#               panel.background = element_rect(fill = "white", colour = "grey50"),
#               plot.title = element_text(hjust = 0.5), 
#               axis.text.x = element_text(angle = 60, hjust = 1),
#               axis.text = element_text(colour = 1, size = 9),
#               panel.border = element_rect(colour = "black", fill=NA))
#       
#       result <- result + theme(legend.position = c(0.53, 0))
#       result <- ggdraw() +
#         draw_plot(result) +
#         draw_plot(inset.plot, x = 0.65, y = 0.0, width = 0.3, height = 0.32)
#     }
#     plot_Deathsdelta_age2[[plot_name]][[sp]][[paste(t, collapse="_")]]<- result
#   }
#   
#   
#   
#   #### PLOT SURVIVAL ----------------------------
#   df_trait_sp <- final_dataset_age %>% filter(Species == sp, Fstart == 1)
#   
#   # avoid the estimation bias of very populations of very low size
#   if (any(((df_trait_sp[which(df_trait_sp$variable == "abundanceByAge"), "value"])<10))){
#     df_trait_sp <- pivot_wider(df_trait_sp, names_from = "variable", values_from = "value")
#     df_trait_sp <- df_trait_sp[-which((df_trait_sp$abundanceByAge)<10), ]
#     df_trait_sp <- pivot_longer(df_trait_sp, cols = -c("age", "replicate", "Fmax", "Fstart", 'Species', "cohort"), names_to = "variable", values_to = "value")
#   }
#   
#   # remove cohort dimension
#   df_trait_sp <- df_trait_sp %>%
#     mutate(age = as.numeric(age))%>%
#     mutate(replicate = as.numeric(replicate)-1)%>% # to make the datasets match : unify the replicate range
#     mutate(Fmax = as.numeric(Fmax))%>%
#     mutate(Fstart = as.numeric(Fstart))
#   weight_df <- build_weight_df(df_trait_sp, weight_var = "abundanceByAge", key_cols = c("age", "Fmax", "Fstart", "Species", "replicate"))
#   df_trait_sp <- weighted_reduce(
#     df = df_trait_sp,
#     weight_df = weight_df,
#     join_cols = c("Species", "replicate", "Fmax", "Fstart", "age"),
#     group_cols = c("Species", "replicate", "Fstart", "Fmax", "variable", "age"),
#     drop_age_zero = TRUE
#   )
#   
#   # combine all traits: Mortality + Ne traits : both need to be pivot_wider to be combined
#   df_trait_sp <- pivot_wider(df_trait_sp, names_from = "variable", values_from = "value")
#   data_complete_c_i <- data_complete_c %>%
#     filter(Species == sp)
#   data_complete_c_i <- pivot_wider(data_complete_c_i, names_from = "variable", values_from = "value")
#   df_trait_sp <- full_join(data_complete_c_i,df_trait_sp, by=c("age", "Fmax", "Fstart", "Species", "replicate"))
#   if (any(df_trait_sp$age == 0)) df_trait_sp <- df_trait_sp[-which(df_trait_sp$age == 0), ]
#   df_trait_sp <- na.omit(df_trait_sp)
#   df_trait_sp <- pivot_longer(data = df_trait_sp, cols = -which(colnames(df_trait_sp) %in% c("age", "Fmax", "Fstart", "Species", "replicate")), 
#                               names_to = "variable", values_to = "value")
#   
#   #add lxbx
#   df_trait_t    <- estimateX(trait_data_full_i=df_trait_sp, var=paste0("survivorship_curvexb_rate"))
#   df_trait_sp   <- rbind(df_trait_sp, df_trait_t)
#   df_trait_t    <- estimateX(trait_data_full_i=df_trait_sp, var=paste0("survivorship_curvexb_rate_sc"))
#   df_trait_sp   <- rbind(df_trait_sp, df_trait_t)
#   
#   var_i <- unique(df_trait_sp$variable)
#   
#   # colSumDeltasquarexDead
#   df_trait_t    <- estimatecolSum(trait_data_full_i = df_trait_sp, var=paste0("colSumDeltasquarexDead"))
#   df_trait_sp   <- rbind(df_trait_sp, df_trait_t)  
#   
#   # Bx = bxNx
#   df_trait_t    <- estimateX(trait_data_full_i = df_trait_sp, var=paste0("b_ratexabundanceByAge"))
#   df_trait_t$variable <- "B_age"
#   df_trait_sp   <- rbind(df_trait_sp, df_trait_t)  
#   
#   # Bx = agebxNx
#   df_trait_t    <- df_trait_sp[which(df_trait_sp$variable == "B_age"), ]
#   df_trait_t$value    <- df_trait_t$value[which(df_trait_t$variable == "B_age")]*df_trait_sp$age[which(df_trait_t$variable == "B_age")]
#   df_trait_t$variable <- "ageB_age (generation length)"
#   df_trait_sp   <- rbind(df_trait_sp, df_trait_t)  
#   
#   # Generation length
#   df_trait_t    <- estimateDiv(trait_data_full_i = df_trait_sp, var=paste0("ageB_age (generation length)\nDivNumerator_N_alpha"))
#   df_trait_sp   <- rbind(df_trait_sp, df_trait_t) 
#   df_trait_t    <- estimatecolSum(trait_data_full_i = df_trait_sp, var=paste0("colSumxBx\nDivNumerator_N_alpha"))
#   df_trait_t$variable <- "Generation length"
#   df_trait_sp   <- rbind(df_trait_sp, df_trait_t)  
#   
#   # Nx = lxNx
#   df_trait_t    <- estimateX(trait_data_full_i = df_trait_sp, var=paste0("survivorship_curvexNumerator_N_alpha"))
#   df_trait_sp   <- rbind(df_trait_sp, df_trait_t)
#   
#   
#   # delta square
#   df_trait_t    <- estimateSquare(trait_data_full_i = df_trait_sp, var=paste0("Deltasquare"))
#   df_trait_sp   <- rbind(df_trait_sp, df_trait_t)
#   
#   var_i <- unique(df_trait_sp$variable)
#   # second kind of traits vector construction
#   for (v in var_i){
#     if (v %in% c("k_point", "kx_bar", "Dead", "Deltasquare")){
#       df_trait_t    <- estimateDiv(trait_data_full_i = df_trait_sp, var=paste0(v, "\nDivNumerator_N_alpha")) # div by population size
#       df_trait_sp   <- rbind(df_trait_sp, df_trait_t)
#     }
#   }
#   
#   
#   # keep abundance by age variable on the side
#   WEIGHT = df_trait_sp[which(df_trait_sp$variable=="abundanceByAge"),]
#   
#   var_i <- unique(df_trait_sp$variable)
#   
#   for (relative_to_ref in c(FALSE)){
#     proportional = FALSE
#     if (relative_to_ref == "proportional") {
#       relative_to_ref = FALSE
#       proportional = TRUE
#     }
#     relative_to_ref = as.logical(relative_to_ref)
#     plot_name_ch=c()
#     plot_name <- paste0(as.character(relative_to_ref), as.character(proportional))
#     if (plot_name=="TRUEFALSE") plot_name_ch <- c("Relative")
#     if (plot_name=="FALSEFALSE") plot_name_ch <- c("Absolute")
#     if (plot_name=="FALSETRUE") plot_name_ch <- c("Proportional")
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
#     v <- list( 
#       unique(names(table(df_trait_sp_a$variable)[c(
#         which(names(table(df_trait_sp_a$variable))=="b_rate"), 
#         which(names(table(df_trait_sp_a$variable))=="B_age"), 
#         # which(names(table(df_trait_sp_a$variable))=="ageB_age (generation length)"), 
#         grep(names(table(df_trait_sp_a$variable)), pattern="survivorship"))])),
#       
#       unique(names(table(df_trait_sp_a$variable)[c(which(names(table(df_trait_sp_a$variable)) == "k_point"), 
#              which(names(table(df_trait_sp_a$variable)) == "kx_bar"), 
#              which(names(table(df_trait_sp_a$variable)) == "k_point\nDivNumerator_N_alpha"), 
#              which(names(table(df_trait_sp_a$variable)) == "kx_bar\nDivNumerator_N_alpha"))])),    
#     
#       unique(names(table(df_trait_sp_a$variable)[c(which(names(table(df_trait_sp_a$variable)) == "Deltasquare"), 
#              which(names(table(df_trait_sp_a$variable)) == "Deltasquare\nDivNumerator_N_alpha"), 
#              which(names(table(df_trait_sp_a$variable)) == "Dead"), 
#              which(names(table(df_trait_sp_a$variable)) == "Dead\nDivNumerator_N_alpha"))]))
#     )
#     # v <- v[-which(v == "survivorship_curvexNumerator_N_alpha")]
#     
#     
#     if (all(is.na(df_trait_sp_a$value))) next
#     for (v_ii in 1:length(v)){
#       v_i = v[[v_ii]]
#       df_trait_sp_plot <- df_trait_sp_a[which(df_trait_sp_a$variable %in% v_i),]
#       t <- improve_vriable_names(v_i)
#       df_trait_sp_plot$variable <- improve_vriable_names(df_trait_sp_plot$variable)
#       df_trait_sp_plot$variable <- factor(as.factor(df_trait_sp_plot$variable),levels=c(t))
#       df_trait_sp_plot <- df_trait_sp_plot %>%
#         group_by(variable, Fmax, Species, Fstart) %>% 
#         mutate(q99 = quantile(value, 0.99)) %>%
#         filter(value <= q99) %>%
#         select(-c(q99))
#       ncol=ifelse(length(t)>3, 3, 2)
#       result <- plot_trait_vs_fmax(df_trait_sp=df_trait_sp_plot, spec=sp, age_dim = T, sp_dim = F,  
#                                    trait_name=t, Time1=150, Time2=160, relative_to_ref = relative_to_ref,
#                                    # scale_factor = 1, 
#                                    weight=WEIGHT, nrow=ceiling(length(t)/ncol), ncol=ncol, 
#                                    proportional = proportional, position_legend="auto", facet=T, scale_z=F,
#                                    chosen_title = paste0(plot_name_ch, " abundance and related variables per age, per Fmax,\n", 
#                                                          "averaged over years ",  Time1, " and ", Time2))
#       plot_Abund_age[[plot_name]][[sp]][[paste(t, collapse="_")]]<- result
#     }#plot
#   }#relative/proportional/none
#   
# }#sp
# 
# save(plot_Deathsdelta_age, file=paste0(pathoutputRDATA, "/plot_Deathsdelta_age.RData"))
# plot_Deathsdelta_age <- get(load(file=paste0(pathoutputRDATA, "/plot_Deathsdelta_age.RData")))
# for (relative_to_ref in c(FALSE)){ #, FALSE, "proportional"
#   proportional = FALSE
#   if (relative_to_ref == "proportional") {
#     relative_to_ref = FALSE
#     proportional = TRUE
#   }
#   relative_to_ref = as.logical(relative_to_ref)
#   plot_name_ch=c()
#   plot_name <- paste0(as.character(relative_to_ref), as.character(proportional))
#   if (plot_name=="TRUEFALSE") plot_name_ch <- c("relative")
#   if (plot_name=="FALSEFALSE") plot_name_ch <- c("")
#   if (plot_name=="FALSETRUE") plot_name_ch <- c("proportional")
#   
#   windowsFonts(myCustomWindowsFontName=windowsFont("DejaVu Sans"))
#   cairo_pdf(onefile=T, file = paste0(PATHOUT, "/Clean_for_paper/DeathsDeltasquareonN_2_", plot_name_ch, ".pdf"), width = 15, height = 10)
#   p <-(plot_Deathsdelta_age[[plot_name]])
#   print(p)
#   dev.off()
# }
# 
# 
# save(plot_Deathsdelta_age2, file=paste0(pathoutputRDATA, "/plot_Deathsdelta_age2.RData"))
# plot_Deathsdelta_age2 <- get(load(file=paste0(pathoutputRDATA, "/plot_Deathsdelta_age2.RData")))
# for (relative_to_ref in c(FALSE)){ #, FALSE, "proportional"
#   proportional = FALSE
#   if (relative_to_ref == "proportional") {
#     relative_to_ref = FALSE
#     proportional = TRUE
#   }
#   relative_to_ref = as.logical(relative_to_ref)
#   plot_name_ch=c()
#   plot_name <- paste0(as.character(relative_to_ref), as.character(proportional))
#   if (plot_name=="TRUEFALSE") plot_name_ch <- c("relative")
#   if (plot_name=="FALSEFALSE") plot_name_ch <- c("absolute")
#   if (plot_name=="FALSETRUE") plot_name_ch <- c("proportional")
#   
#   windowsFonts(myCustomWindowsFontName=windowsFont("DejaVu Sans"))
#   cairo_pdf(onefile=T, file = paste0(PATHOUT, "/Clean_for_paper/DeathsDeltasquareRAWonN_2_", plot_name_ch, ".pdf"), width = 15, height = 10)
#   p <-(plot_Deathsdelta_age2[[plot_name]])
#   print(p)
#   dev.off()
# }
# 
# save(plot_Abund_age, file=paste0(pathoutputRDATA, "/plot_Abund_age.RData"))
# plot_Abund_age <- get(load(file=paste0(pathoutputRDATA, "/plot_Abund_age.RData")))
# for (relative_to_ref in c(FALSE)){ #, FALSE, "proportional"
#   proportional = FALSE
#   if (relative_to_ref == "proportional") {
#     relative_to_ref = FALSE
#     proportional = TRUE
#   }
#   relative_to_ref = as.logical(relative_to_ref)
#   plot_name_ch=c()
#   plot_name <- paste0(as.character(relative_to_ref), as.character(proportional))
#   if (plot_name=="TRUEFALSE") plot_name_ch <- c("relative")
#   if (plot_name=="FALSEFALSE") plot_name_ch <- c("")
#   if (plot_name=="FALSETRUE") plot_name_ch <- c("proportional")
#   for (trai in names(plot_Abund_age[[plot_name]][["Hake"]])){
#     windowsFonts(myCustomWindowsFontName=windowsFont("DejaVu Sans"))
#     cairo_pdf(file = paste0(PATHOUT, "/Clean_for_paper/Abund_age_lxbx_", plot_name_ch, ".pdf"), width = 15, height = 10, onefile = T)
#     for (sp in species)  p <-(plot_Abund_age[[plot_name]][[sp]][[trai]])
#     print(p)
#     dev.off()
#   }
# }
# 
# 
# 
# 
# 
# 
# #### PLOT 2- SSDg/N1 SSDi/N1 VS FMAX  ------------------
# 
# SSDplot = list()
# plots_by_age= list()
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
#   #### 1- SSD and Vk and DeadSum
#   traits_sp_to_plot = c("colSumSSDix\nDivNumerator_N_alpha", "colSumSSDgx\nDivNumerator_N_alpha", "Vk_point",
#                         "colSumDeltasquare\nDivNumerator_N_alpha", "DeadSum", "colSumDeltasquare\nDivNumerator_N_alpha\nxDead")
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
#       df_trait_t    <- estimateDiv(trait_data_full_i = df_trait_sp, var=paste0(v, "\nDivNumerator_N_alpha"))
#       df_trait_sp   <- rbind(df_trait_sp, df_trait_t)
#     }
#   }
#   
#   #add the SSDgx estimation by hand
#   df_trait_t    <- estimateX(trait_data_full_i = df_trait_sp, var=paste0("Deltasquare\nDivNumerator_N_alpha\nxDead"))
#   df_trait_sp   <- rbind(df_trait_sp, df_trait_t)
#   df_trait_t    <- estimatecolSum(trait_data_full_i = df_trait_sp, var=paste0("colSumDeltasquare\nDivNumerator_N_alpha\nxDead"))
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
#   
#   
#   
#   #### 2- DELTA, DEAD, DEALTA², survival rate
#   traits_sp_to_plot = list(c("Deltasquare","Dead"),
#                            c("DeltasquarexDead", "DeltasquarexDead\nDivNumerator_N_alpha"),
#                            c("survival_rate"))
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
#       df_trait_t    <- estimateDiv(trait_data_full_i = df_trait_sp, var=paste0(v, "\nDivNumerator_N_alpha")) # div by population size
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
#                                    Time1=150, Time2=160, relative_to_ref = relative_to_ref,
#                                    #scale_factor = 1, #
#                                    ncol=2, nrow=2,
#                                    weight=WEIGHT, facet=F,
#                                    proportional = proportional, position_legend="auto")
#       plots_by_age[[paste0(as.character(relative_to_ref), as.character(proportional))]][[sp]][[paste(t, collapse="_")]]<- result$plot
#     }
#   }#relative/proportional/none
# }#sp
# 
# save(SSDplot, file=paste0(pathoutputRDATA, "/SSDplots.RData"))
# load(file=paste0(pathoutputRDATA, "/SSDplots.RData"))
# 
# save(plots_by_age, file=paste0(pathoutputRDATA, "/deaths_deltasquare_byage.RData"))
# deaths_deltasquare_byage <- get(load(file=paste0(pathoutputRDATA, "/deaths_deltasquare_byage.RData")))
# 
# 
# # 
# # #### PLOT 3: mortality decomposition --------------------------
# # mortality_plot = list()
# # for (sp in species){
# #   cat(sp)
# #   
# #   df_trait_sp <- data_complete_c %>% filter(Species == sp, Fstart == 1)
# #   
# #   df_trait_t    <- estimateSum(trait_data_full_i = df_trait_sp, var=paste0("mean_of_F_per_age", "\nsum", "mean_of_Mnat_per_age"))
# #   df_trait_t$variable    <- "sum_mortality_rates"
# #   df_trait_sp   <- rbind(df_trait_sp, df_trait_t)
# #   
# #   var_i <- unique(df_trait_sp$variable)
# #   # compare all mortality rates to their value *abundance to have the number of individuals to understand where is the increase coming from
# #   for (v in var_i){
# #     if (any(grepl(pattern="mean_of_M", v) | grepl(pattern="mean_of_F", v) | grepl(pattern="sum_mortality_rates", v))) {
# #       df_trait_t    <- estimateMnumber(trait_data_full_i = df_trait_sp, var=paste0(v))
# #       df_trait_sp   <- rbind(df_trait_sp, df_trait_t)
# #     }
# #   }
# #   var_i <- unique(df_trait_sp$variable)
# #   for (v in var_i){
# #     if (any(grepl(pattern="mean_of_M", v) | grepl(pattern="mean_of_F", v))) {
# #       df_trait_t    <- estimateDiv(trait_data_full_i = df_trait_sp, var=paste0(v, "\nDiv", "sum_mortality_rates"))
# #       df_trait_sp   <- rbind(df_trait_sp, df_trait_t)
# #     }
# #   }
# #   df_trait_t    <- estimatecolSum(trait_data_full_i = df_trait_sp, var=paste0("colSumall_deaths"))
# #   df_trait_sp   <- rbind(df_trait_sp, df_trait_t)
# #   var_i <- unique(df_trait_sp$variable)
# #   for (v in var_i){
# #     if (any(grepl(pattern="deaths", v) & c(grepl(pattern="M", v)) | (grepl(pattern="F_deaths", v)))) {
# #       df_trait_t    <- estimateDiv(trait_data_full_i = df_trait_sp, var=paste0(v, "\nDiv", "colSumall_deaths"))
# #       df_trait_t2   <- estimatecolSum(trait_data_full_i = df_trait_t, var=paste0("colSum", (df_trait_t$variable)[1]))
# #       df_trait_sp   <- rbind(df_trait_sp, df_trait_t)
# #       df_trait_sp   <- rbind(df_trait_sp, df_trait_t2)
# #     }
# #   }
# #   for (v in var_i){
# #     if (v == "all_deaths") {
# #       df_trait_t2   <- estimatecolSum(trait_data_full_i = df_trait_t, var=paste0("colSum", v))
# #       df_trait_sp   <- rbind(df_trait_sp, df_trait_t2)
# #     }
# #   }
# #   var_i <- unique(df_trait_sp$variable)
# #   df_trait_sp_a <- df_trait_sp
# #   # reduce dimension number
# #   weight_df <- build_weight_df(df = df_trait_sp_a, weight_var = "mean_of_abundanceByAge_per_age",
# #                                key_cols = c("replicate", "Species", "Fmax", "Fstart"))
# #   df_trait_sp_a <- weighted_reduce(
# #     df = df_trait_sp_a,
# #     weight_df = weight_df,
# #     join_cols = c("Species", "replicate", "Fmax", "Fstart"),
# #     group_cols = c("Species", "replicate", "Fstart", "Fmax", "variable"),
# #     drop_age_zero = TRUE
# #   )
# #   
# #   # reference dataset
# #   ref_data <- df_trait_sp_a %>%
# #     filter(Fmax == 1 & Fstart == 1) %>%
# #     dplyr::select(variable, value, replicate, Species) %>%
# #     dplyr::rename(ref_value = value)
# #   
# #   
# #   df_trait_sp_a <- clean_values(df_trait_sp_a)
# #   
# #   # conditions to stop
# #   if (nrow(df_trait_sp_a)==0) next
# #   if (all(is.na(df_trait_sp_a$value) | df_trait_sp_a$value == 0)) {
# #     next
# #   }
# #   for (relative_to_ref in c(TRUE, FALSE, "proportional")){
# #     proportional = FALSE
# #     if (relative_to_ref == "proportional") {
# #       relative_to_ref = FALSE
# #       proportional = TRUE
# #     }
# #     relative_to_ref <- as.logical(relative_to_ref)
# #     
# #     # for the plot
# #     plot_biom = ifelse(relative_to_ref, "rel_biomassplot", "biomassplot")
# #     plot_biom = ifelse(relative_to_ref, "biomassplot", "biomassplot")
# #     
# #     result =NULL
# #     var_i = c("sum_mortality_rates", "colSumall_deaths")
# #     data_split <- df_trait_sp_a[which(df_trait_sp_a$variable %in% var_i ),]
# #     ref_data_split <- ref_data[which(ref_data$variable %in% var_i),]
# #     data_split$variable <- improve_vriable_names(data_split$variable)
# #     var_i <- improve_vriable_names(var_i)
# #     ref_data_split$variable <- improve_vriable_names(ref_data_split$variable)
# #     if (all(is.na(data_split$value))) next
# #     result <- plot_trait_vs_fmax(df_trait_sp = data_split, trait_name = var_i,
# #                                  Time1=150, Time2=160, relative_to_ref = relative_to_ref, spec=sp,
# #                                  age_dim=F, reference_input=ref_data_split, proportional = proportional,
# #                                  weight=weight_df,
# #                                  chosen_title=paste0(ifelse(proportional, "Proportional ", ""),
# #                                                      ifelse(relative_to_ref, "Relative ", ""),
# #                                                      "Total Mortality rates and Deaths numbers\naverage years ", Time1, " to ",Time2))
# #     mortality_plot[[paste0(as.character(relative_to_ref), as.character(proportional))]][[sp]][[stringr::str_remove(paste(var_i, collapse=""), pattern="\n")]]    <- result$plot
# #     
# #     result =NULL
# #     var_i = c("mean_of_Mstar_per_age\nDivsum_mortality_rates", "mean_of_Mpred_per_age\nDivsum_mortality_rates",
# #               "mean_of_F_per_age\nDivsum_mortality_rates", "mean_of_Madd_per_age\nDivsum_mortality_rates")
# #     data_split <- df_trait_sp_a[which(df_trait_sp_a$variable %in% var_i ),]
# #     ref_data_split <- ref_data[which(ref_data$variable %in% var_i),]
# #     data_split$variable <- improve_vriable_names(data_split$variable)
# #     var_i <- improve_vriable_names(var_i)
# #     ref_data_split$variable <- improve_vriable_names(ref_data_split$variable)
# #     if (all(is.na(data_split$value))) next
# #     result <- plot_trait_vs_fmax(df_trait_sp = data_split, trait_name = var_i,
# #                                  Time1=150, Time2=160, relative_to_ref = relative_to_ref, spec=sp,
# #                                  age_dim=F, reference_input=ref_data_split, proportional = proportional,
# #                                  weight=weight_df,
# #                                  chosen_title=paste0(ifelse(proportional, "Proportional ", ""),
# #                                                      ifelse(relative_to_ref, "Relative ", ""),
# #                                                      "Mortality rates\naverage years ", Time1, " to ",Time2))
# #     mortality_plot[[paste0(as.character(relative_to_ref), as.character(proportional))]][[sp]][[stringr::str_remove(paste(var_i, collapse=""), pattern="\n")]]  <- result$plot
# #     
# #     result =NULL
# #     var_i = c("colSumMstar_deaths\nDivcolSumall_deaths", "colSumF_deaths\nDivcolSumall_deaths",
# #               "colSumMpred_deaths\nDivcolSumall_deaths", "colSumMadd_deaths\nDivcolSumall_deaths")
# #     data_split <- df_trait_sp_a[which(df_trait_sp_a$variable %in% var_i ),]
# #     ref_data_split <- ref_data[which(ref_data$variable %in% var_i),]
# #     data_split$variable <- improve_vriable_names(data_split$variable)
# #     var_i <- improve_vriable_names(var_i)
# #     ref_data_split$variable <- improve_vriable_names(ref_data_split$variable)
# #     if (all(is.na(data_split$value))) next
# #     result <- plot_trait_vs_fmax(df_trait_sp = data_split, trait_name = var_i,
# #                                  Time1=150, Time2=160, relative_to_ref = relative_to_ref, spec=sp,
# #                                  age_dim=F, reference_input=ref_data_split, proportional = proportional,
# #                                  weight=weight_df,
# #                                  chosen_title=paste0(ifelse(proportional, "Proportional ", ""),
# #                                                      ifelse(relative_to_ref, "Relative ", ""),
# #                                                      "Death number per Mortality rate\naverage years ", Time1, " to ",Time2))
# #     mortality_plot[[paste0(as.character(relative_to_ref), as.character(proportional))]][[sp]][[stringr::str_remove(paste(var_i, collapse=""), pattern="\n")]]  <- result$plot
# #   }#sp
# # }
# # save(mortality_plot, file=paste0(pathoutputRDATA, "/mortalityplots.RData"))
# # load(file=paste0(pathoutputRDATA, "/mortalityplots.RData"))
# 
# 
# # 
# # # mortality & age 
# # PATHOUT  = "output/PLOTS/KEY_TRAITS/"
# # mortality_plot_age = list()
# # for (sp in species){
# #   cat(sp)
# #   
# #   df_trait_sp <- data_complete_c %>% filter(Species == sp, Fstart == 1)
# #   
# #   df_trait_t    <- estimateSum(trait_data_full_i = df_trait_sp, var=paste0("mean_of_F_per_age", "\nsum", "mean_of_Mnat_per_age"))
# #   df_trait_t$variable    <- "sum_mortality_rates"
# #   df_trait_sp   <- rbind(df_trait_sp, df_trait_t)
# #   
# #   var_i <- unique(df_trait_sp$variable)
# #   # compare all mortality rates to their value *abundance to have the number of individuals to understand where is the increase coming from
# #   for (v in var_i){
# #     if (any(grepl(pattern="mean_of_M", v) | grepl(pattern="mean_of_F", v) | grepl(pattern="sum_mortality_rates", v))) {
# #       df_trait_t    <- estimateMnumber(trait_data_full_i = df_trait_sp, var=paste0(v))
# #       df_trait_sp   <- rbind(df_trait_sp, df_trait_t)
# #     }
# #   }
# #   var_i <- unique(df_trait_sp$variable)
# #   for (v in var_i){
# #     if (any(grepl(pattern="mean_of_M", v) | grepl(pattern="mean_of_F", v)| (v=="sum_mortality_rates"))) {
# #       df_trait_t    <- estimateDiv(trait_data_full_i = df_trait_sp, var=paste0(v, "\nDiv", "sum_mortality_rates"))
# #       df_trait_sp   <- rbind(df_trait_sp, df_trait_t)
# #     }
# #   }
# #   df_trait_t    <- estimatecolSum(trait_data_full_i = df_trait_sp, var=paste0("colSumall_deaths"))
# #   df_trait_sp   <- rbind(df_trait_sp, df_trait_t)
# #   var_i <- unique(df_trait_sp$variable)
# #   for (v in var_i){
# #     if (any(grepl(pattern="deaths", v) & c(grepl(pattern="M", v)) | (grepl(pattern="F_deaths", v)) | (v=="all_deaths"))) {
# #       df_trait_t    <- estimateDiv(trait_data_full_i = df_trait_sp, var=paste0(v, "\nDiv", "colSumall_deaths"))
# #       df_trait_t2   <- estimatecolSum(trait_data_full_i = df_trait_t, var=paste0("colSum", (df_trait_t$variable)[1]))
# #       df_trait_sp   <- rbind(df_trait_sp, df_trait_t)
# #       df_trait_sp   <- rbind(df_trait_sp, df_trait_t2)
# #     }
# #   }
# #   df_trait_t    <- estimateage1(trait_data_full_i = df_trait_sp, var=paste0("age1mean_of_abundanceByAge_per_age"))
# #   df_trait_sp   <- rbind(df_trait_sp, df_trait_t)
# #   df_trait_t    <- estimateDiv(trait_data_full_i = df_trait_sp, var=paste0("all_deaths\nDiv", "age1mean_of_abundanceByAge_per_age"))
# #   df_trait_sp   <- rbind(df_trait_sp, df_trait_t)
# #   for (v in var_i){
# #     if (v == "all_deaths") {
# #       df_trait_t2   <- estimatecolSum(trait_data_full_i = df_trait_t, var=paste0("colSum", v))
# #       df_trait_sp   <- rbind(df_trait_sp, df_trait_t2)
# #     }
# #   }
# #   var_i <- unique(df_trait_sp$variable)
# #   df_trait_sp_a <- df_trait_sp
# #   # reduce dimension number
# #   weight_df <- build_weight_df(df = df_trait_sp, weight_var = "mean_of_abundanceByAge_per_age",
# #                                key_cols = c("replicate", "Species", "Fmax", "Fstart", "age"))
# #   df_trait_sp <- weighted_reduce(
# #     df = df_trait_sp,
# #     weight_df = weight_df,
# #     join_cols = c("Species", "replicate", "Fmax", "Fstart", "age"),
# #     group_cols = c("Species", "replicate", "Fstart", "Fmax", "variable", "age"),
# #     drop_age_zero = TRUE
# #   )
# #   
# #   # reference dataset
# #   ref_data <- df_trait_sp_a %>%
# #     filter(Fmax == 1 & Fstart == 1) %>%
# #     dplyr::select(variable, value, replicate, age, Species) %>%
# #     dplyr::rename(ref_value = value)
# #   
# #   
# #   df_trait_sp_a <- clean_values(df_trait_sp_a)
# #   
# #   # conditions to stop
# #   if (nrow(df_trait_sp_a)==0) next
# #   if (all(is.na(df_trait_sp_a$value) | df_trait_sp_a$value == 0)) {
# #     next
# #   }
# #   
# #   
# #   for (relative_to_ref in c(TRUE, FALSE, "proportional")){
# #     proportional = FALSE
# #     if (relative_to_ref == "proportional") {
# #       relative_to_ref = FALSE
# #       proportional = TRUE
# #     }
# #     relative_to_ref <- as.logical(relative_to_ref)
# #     
# #     # for the plot
# #     plot_biom = ifelse(relative_to_ref, "rel_biomassplot", "biomassplot")
# #     plot_biom = ifelse(relative_to_ref, "biomassplot", "biomassplot")
# #     
# #     
# #     result =NULL
# #     var_i = c("all_deaths\nDivcolSumall_deaths", )
# #     data_split <- df_trait_sp_a[which(df_trait_sp_a$variable %in% var_i ),]
# #     ref_data_split <- ref_data[which(ref_data$variable %in% var_i),]
# #     data_split$variable <- improve_vriable_names(data_split$variable)
# #     var_i <- improve_vriable_names(var_i)
# #     ref_data_split$variable <- improve_vriable_names(ref_data_split$variable)
# #     if (all(is.na(data_split$value))) next
# #     result <- plot_trait_vs_fmax(df_trait_sp = data_split, trait_name = var_i,
# #                                  Time1=150, Time2=160, relative_to_ref = relative_to_ref, spec=sp,
# #                                  age_dim=T, reference_input=ref_data_split, proportional = proportional,
# #                                  weight=weight_df,
# #                                  chosen_title=paste0(ifelse(proportional, "Proportional ", ""),
# #                                                      ifelse(relative_to_ref, "Relative ", ""),
# #                                                      "Total Mortality rates and Deaths numbers\naverage years ", Time1, " to ",Time2))
# #     mortality_plot_age[[paste0(as.character(relative_to_ref), as.character(proportional))]][[sp]][[stringr::str_remove(paste(var_i, collapse=""), pattern="\n")]]    <- result$plot
# #     
# #     result =NULL
# #     var_i = c("Mstar_deaths\nDivcolSumall_deaths", "F_deaths\nDivcolSumall_deaths",
# #               "Mpred_deaths\nDivcolSumall_deaths", "Madd_deaths\nDivcolSumall_deaths")
# #     data_split <- df_trait_sp_a[which(df_trait_sp_a$variable %in% var_i ),]
# #     ref_data_split <- ref_data[which(ref_data$variable %in% var_i),]
# #     data_split$variable <- improve_vriable_names(data_split$variable)
# #     var_i <- improve_vriable_names(var_i)
# #     ref_data_split$variable <- improve_vriable_names(ref_data_split$variable)
# #     if (all(is.na(data_split$value))) next
# #     result <- plot_trait_vs_fmax(df_trait_sp = data_split, trait_name = var_i,
# #                                  Time1=150, Time2=160, relative_to_ref = relative_to_ref, spec=sp,
# #                                  age_dim=T, reference_input=ref_data_split, proportional = proportional,
# #                                  weight=weight_df,
# #                                  chosen_title=paste0(ifelse(proportional, "Proportional ", ""),
# #                                                      ifelse(relative_to_ref, "Relative ", ""),
# #                                                      "Death number per Mortality rate\naverage years ", Time1, " to ",Time2))
# #     mortality_plot_age[[paste0(as.character(relative_to_ref), as.character(proportional))]][[sp]][[stringr::str_remove(paste(var_i, collapse=""), pattern="\n")]]  <- result$plot
# #   }#relative
# # }#sp
# # save(mortality_plot_age, file=paste0(pathoutputRDATA, "/mortalityplots_age.RData"))
# # load(file=paste0(pathoutputRDATA, "/mortalityplots_age.RData"))
# 
# 
# 
