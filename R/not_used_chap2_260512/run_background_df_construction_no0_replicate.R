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
path <- paste0(simu_repo, "/output/PLOTS/FINAL_Ne_Htz_mature/")
# path <- paste0(simu_repo, "/output/PLOTS/FINAL_Ne_Htz_mature_short/")
PATHOUT  = "output/PLOTS/KEY_TRAITS_MATURE/"
if (!dir.exists(PATHOUT)) dir.create(PATHOUT, recursive = TRUE)
pathoutputRDATA = paste0(PATHOUT, "/RDATA_MATURE")
if(!dir.exists(pathoutputRDATA)) dir.create(pathoutputRDATA, recursive = T)
base_dir <- paste0(path, "INIT")
folders  <- list.dirs(base_dir, recursive = F) # files INIT2 with raw Ne datasets
# folders_long <- list.dirs(paste0(path_long, "INIT"), recursive = F)

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


## RAW Ne & associated variables, dimensions : cohort, age, replicate, Fstart, Fmax, Species
### produce dataset reading multiple .RData to make 1 dataframe, stored as a .RData file ----------
sp_i = "Hake"
FORCE=T
full_dataset = list()
full_dataset_age = list()
dataset="data_Ne"
name_out = ifelse(dataset=="data_Ne", "Complete_INIT2", "Complete_INIT2_raw")
folders <- list(folders) #, folders_long)
path = list(path) #, path_long)
path = path[[ifelse(dataset=="data_Ne", 1, 2)]]
for (f_i in folders[[ifelse(dataset=="data_Ne", 1, 2)]]){
  cat(f_i, "\n")
  
  Fmax   <- as.numeric(stringr::str_split(stringr::str_split(f_i, pattern="Fmax")[[1]][2], pattern="_Fstart")[[1]][1])
  Fstart <- as.numeric(stringr::str_split(stringr::str_split(f_i, pattern="Fstart")[[1]][2], pattern="_sp")[[1]][1])
  
  if (Fstart != 1) next
  
  load(list.files(f_i, full.names = T))
  dataNe <- INIT2[[dataset]]
  rm(INIT2)
  
  for (sp in names(dataNe)[1]){
    if (file.exists(paste0(path, "/", name_out, sp_i, ".RData")) && FORCE !=T){
      next
    }else{
      
      cat(sp, "\n")
      
      if(dataset == "data_Ne") {
        f_i_sp <- dataNe[[sp]]
        f_i_sp[["NT"]] <- do.call(cbind, (lapply(f_i_sp$abundanceByAge, FUN=function(x) colSums(x))))
      }
      if(dataset != "data_Ne") f_i_sp <- dataNe
      
      
      for (dim in names(f_i_sp)){
        if (dim == "cohort_time_step_corresponding"){next}
        f_i_sp_d <- f_i_sp[[dim]]
        if(dataset != "data_Ne") {
          sp = names(f_i_sp[[dim]])[1]
          f_i_sp_d <- list()
          for (rep in 1:dim(f_i_sp[[dim]][[sp]])[3]){ # all the dimensions variables are considered from age 1, thus Ne is estimated by variables from age 1, except for Numerator_N_alpha
            f_i_sp_d[[paste0("rep_", rep)]] <- f_i_sp[[dim]][[sp]][which(as.numeric(rownames(f_i_sp[[dim]][[sp]]))>=1), ,paste0("rep_", rep)]
          }
        }
        
        if (is.data.frame(f_i_sp_d) || is.matrix(f_i_sp_d)){
          
          if (is.matrix(f_i_sp_d)) f_i_sp_d <- as.data.frame(f_i_sp_d)
          f_i_sp_d$cohort <- rownames(f_i_sp_d)
          f_i_sp_d$cohort = stringr::str_remove(as.matrix(f_i_sp_d[, "cohort"]), pattern="cohort_")
          f_i_sp_d <- pivot_longer(f_i_sp_d, names_to = "replicate", values_to = "value", cols=-cohort)
          f_i_sp_d$replicate = stringr::str_remove(as.matrix(f_i_sp_d[, "replicate"]), pattern="rep_")
          f_i_sp_d <- f_i_sp_d[which(f_i_sp_d$cohort %in% c(Time1:Time2)), ]
          f_i_sp_d$Species <- sp
          f_i_sp_d$variable <- dim
          f_i_sp_d$Fmax <- Fmax
          f_i_sp_d$Fstart <- Fstart
          
          # harmonic mean
          data_Ne_Nalpha = NULL
          vector_harmonic_mean = c("Ne_N_alpha", "Numerator_N_alpha", "NT")
          vector_harmonic_mean <- c(vector_harmonic_mean, improve_vriable_names(vector_harmonic_mean))
          if (any(f_i_sp_d$variable %in% vector_harmonic_mean & "cohort" %in% colnames(f_i_sp_d))){
            data_Ne_Nalpha <- f_i_sp_d [which(f_i_sp_d$variable %in% vector_harmonic_mean), ]
            f_i_sp_d   <- f_i_sp_d [-which(f_i_sp_d$variable %in% vector_harmonic_mean), ]
            data_Ne_Nalpha <- data_Ne_Nalpha %>%
              group_by(Fmax, Fstart, variable, replicate, Species) %>%
              summarise(value = psych::harmonic.mean(value, na.rm = TRUE))
          }
          # arithmetic mean mean
          f_i_sp_d <- f_i_sp_d %>%
            group_by(Fmax, Fstart, variable, replicate, Species) %>%
            summarise(value = mean(value, na.rm = TRUE))
          # bind arithmetic and harmonic mean dataset
          if (!is.null(data_Ne_Nalpha)) {
            f_i_sp_d = rbind(data_Ne_Nalpha,  f_i_sp_d)
          }
          
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
          f_i_sp_d <- f_i_sp_d[which(f_i_sp_d$cohort %in% c(Time1:Time2)), ]
          f_i_sp_d$Species <- sp
          f_i_sp_d$variable <- dim
          f_i_sp_d$Fmax <- stringr::str_split(stringr::str_split(f_i, pattern="Fmax")[[1]][2], pattern="_Fstart")[[1]][1]
          f_i_sp_d$Fstart <- stringr::str_split(stringr::str_split(f_i, pattern="Fstart")[[1]][2], pattern="_sp")[[1]][1]
          
          
          # mean of the cohorts : arithmetic or harmonic
          # harmonic mean
          data_Ne_Nalpha = NULL
          vector_harmonic_mean = c("Ne_N_alpha", "Numerator_N_alpha", "NT",
                                   "Ne_N_alpha/Numerator_N_alpha")
          vector_harmonic_mean <- c(vector_harmonic_mean, improve_vriable_names(vector_harmonic_mean))
          if (any(f_i_sp_d$variable %in% vector_harmonic_mean & "cohort" %in% colnames(f_i_sp_d))){
            data_Ne_Nalpha <- f_i_sp_d [which(f_i_sp_d$variable %in% vector_harmonic_mean), ]
            f_i_sp_d   <- f_i_sp_d [-which(f_i_sp_d$variable %in% vector_harmonic_mean), ]
            data_Ne_Nalpha <- data_Ne_Nalpha %>%
              group_by(age, Fmax, Fstart, variable, replicate, Species) %>%
              summarise(value = psych::harmonic.mean(value, na.rm = TRUE))
          }
          # bind arithmetic and harmonic mean dataset
          if (!is.null(data_Ne_Nalpha)) {
            data_Ne_Nalpha$cohort <- 1
            df_trait_sp = rbind(data_Ne_Nalpha, df_trait_sp)
          }
          # arithmetic mean mean
          f_i_sp_d <- f_i_sp_d %>%
            group_by(age, Fmax, Fstart, variable, replicate, Species) %>%
            summarise(value = mean(value, na.rm = TRUE))
          
          
          if (is.null(full_dataset_age[[sp]])){
            full_dataset_age[[sp]] <- f_i_sp_d
          } else {
            full_dataset_age[[sp]] <- rbind(full_dataset_age[[sp]], f_i_sp_d)
          }
          
        } # if list
      }# variable dimension
    }# only if file does not exist yet
  }# species dimension
}# folder dimension (Fmax*Fstart)
#----- MERGE raw variables without age dimension AND variables with age dimension ---#
for (sp_i in names(full_dataset_age)){
  if (file.exists(paste0(path, "/", name_out, sp_i, ".RData")) && FORCE !=T){
    next
  }else{
    
    # produce a dataset without age dimension : weighted mean
    # weight_df <- full_dataset_age[[sp_i]] %>%
    #   dplyr::filter(variable == "abundanceByAge") %>%
    #   dplyr::select(Species, replicate, Fstart, Fmax, cohort, age, weight_value = value)
    # data_with_w <- full_dataset_age[[sp_i]] %>%
    #   dplyr::left_join(weight_df, by = c("Species", "replicate", "cohort", "age", "Fmax", "Fstart"), relationship = "many-to-many")
    # data_with_w <- data_with_w[-which(data_with_w$weight_value == 0), ]
    # data_with_w$age <- as.numeric(data_with_w$age)
    data_out <- full_dataset_age[[sp_i]] %>%
      dplyr::group_by(Species, replicate, Fstart, Fmax, variable) %>%
      dplyr::summarise(
        value = mean(value, na.rm=T), #stats::weighted.mean(value, w = weight_value, na.rm = TRUE),
        .groups = "drop"
      )
    
    data_out$Fmax <- as.numeric(data_out$Fmax)
    data_out$Fstart <- as.numeric(data_out$Fstart)
    final_dataset <- rbind(full_dataset[[sp_i]], data_out)
    
    # produce a dataste with the age dimension
    # full_dataset_add_age <- lapply(names(table(full_dataset_age[[sp_i]]$age)), FUN = function(x){
    #   out = as.data.frame(full_dataset_age[[sp_i]])
    #   out$age = x
    #   return(out)
    # })
    # full_dataset_add_age <- do.call(rbind, full_dataset_age)
    
    final_dataset_age <- (full_dataset_age[[sp_i]])
    
    
    # save all
    complete = list("complete_without_age"=final_dataset, "complete_with_age"=final_dataset_age)
    
    save(complete, file=paste0(path, "/", name_out, sp_i, ".RData"))
  }
}




# 
# 
# #### Initialisation --------------
# library(legendry)
# library(ggpubr)
# 
# # functions and environment
# main = "New_gpu_conf.R" # config file
# source("R/functions/source.R") # load config and key parameters and functions
# source("R/functions/functions_key_traits.R") # plots against Fmax
# source("R/functions/auxiliary_for_lists.R") # from list to df
# # source("R/functions/functions_output_SimuRData_into_df_for_plots.R") # to transform dataset from RData to long df
# source("R/functions/function_data_traits_summary.R")
# 
# # paths
# simu_repo <- "C:/Users/mbeneat/Documents/osmose/Chap2/osmose-ns/osmose-ns-v4.4.1_clim/"
# path <- paste0(simu_repo, "/output/PLOTS/FINAL_Ne_Htz_mature/")
# # path <- paste0(simu_repo, "/output/PLOTS/FINAL_Ne_Htz_mature_short/") # not used anymore
# PATHOUT  = "output/PLOTS/KEY_TRAITS_MATURE/"
# if (!dir.exists(PATHOUT)) dir.create(PATHOUT, recursive = TRUE)
# pathoutputRDATA = paste0(PATHOUT, "/RDATA_MATURE")
# if(!dir.exists(pathoutputRDATA)) dir.create(pathoutputRDATA, recursive = T)
# base_dir <- paste0(path, "INIT")
# folders  <- list.dirs(base_dir, recursive = F) # files INIT2 with raw Ne datasets
# # folders_long <- list.dirs(paste0(path_long, "INIT"), recursive = F)
# 
# folders_traits_y <- paste0(simu_repo, "output/PLOTS/Trait_variables_by_year_BYSP/RData/")
# folders_traits_y <- list.dirs(folders_traits_y, recursive = F)
# folders_traits_c <- paste0(simu_repo, "output/PLOTS/Trait_variables_by_cohort/RData/")
# folders_traits_c  <- list.dirs(folders_traits_c, recursive = F)
# folders_traits_y_a <- paste0(simu_repo, "output/PLOTS/Trait_variables_by_year_by_age2/RData/")
# folders_traits_y_a <- list.dirs(folders_traits_y_a, recursive = F)
# folders_traits_c_a <- paste0(simu_repo, "output/PLOTS/Trait_variables_by_cohort_by_age2/RData/")
# folders_traits_c_a  <- list.dirs(folders_traits_c_a, recursive = F)
# 
# #### DATASETS load or format dataset ------------------------
# Time1=150
# Time2=160
# species_df=species
# 
# 
# ## RAW Ne & associated variables, dimensions : cohort, age, replicate, Fstart, Fmax, Species
# ### produce dataset reading multiple .RData to make 1 dataframe, stored as a .RData file ----------
# sp_i = "Hake"
# FORCE=T
# full_dataset = list()
# full_dataset_age = list()
# dataset="data_Ne"
# name_out = ifelse(dataset=="data_Ne", "Complete_INIT2", "Complete_INIT2_raw")
# folders <- list(folders)
# path = list(path)
# path = path[[ifelse(dataset=="data_Ne", 1, 2)]]
# for (f_i in folders[[ifelse(dataset=="data_Ne", 1, 2)]]){
#   cat(f_i, "\n")
#   
#   Fmax   <- as.numeric(stringr::str_split(stringr::str_split(f_i, pattern="Fmax")[[1]][2], pattern="_Fstart")[[1]][1])
#   Fstart <- as.numeric(stringr::str_split(stringr::str_split(f_i, pattern="Fstart")[[1]][2], pattern="_sp")[[1]][1])
#   
#   if (Fstart != 1) next
#   
#   load(list.files(f_i, full.names = T))
#   dataNe <- INIT2[[dataset]]
#   rm(INIT2)
#   
#   for (sp in names(dataNe)[1]){
#     if (file.exists(paste0(path, "/", name_out, sp_i, ".RData")) && FORCE !=T){
#       next
#     }else{
#       
#       cat(sp, "\n")
#       
#       if(dataset == "data_Ne") {
#         f_i_sp <- dataNe[[sp]]
#         f_i_sp[["Ne_N_alpha/Numerator_N_alpha"]] <- f_i_sp$Ne_N_alpha/f_i_sp$Numerator_N_alpha
#       }
#       if(dataset != "data_Ne") f_i_sp <- dataNe
#       
#       
#       for (dim in names(f_i_sp)){
#         if (dim == "cohort_time_step_corresponding"){next}
#         f_i_sp_d <- f_i_sp[[dim]]
#         if(dataset != "data_Ne") {
#           sp = names(f_i_sp[[dim]])[1]
#           f_i_sp_d <- list()
#           for (rep in 1:dim(f_i_sp[[dim]][[sp]])[3]){ # all the dimensions variables are considered from age 1, thus Ne is estimated by variables from age 1, except for Numerator_N_alpha
#             f_i_sp_d[[paste0("rep_", rep)]] <- f_i_sp[[dim]][[sp]][which(as.numeric(rownames(f_i_sp[[dim]][[sp]]))>=1), ,paste0("rep_", rep)]
#           }
#         }
#         
#         if (is.data.frame(f_i_sp_d) || is.matrix(f_i_sp_d)){
#           
#           if (is.matrix(f_i_sp_d)) f_i_sp_d <- as.data.frame(f_i_sp_d)
#           f_i_sp_d$cohort <- rownames(f_i_sp_d)
#           f_i_sp_d$cohort = stringr::str_remove(as.matrix(f_i_sp_d[, "cohort"]), pattern="cohort_")
#           f_i_sp_d <- pivot_longer(f_i_sp_d, names_to = "replicate", values_to = "value", cols=-cohort)
#           f_i_sp_d$replicate = stringr::str_remove(as.matrix(f_i_sp_d[, "replicate"]), pattern="rep_")
#           f_i_sp_d <- f_i_sp_d[which(f_i_sp_d$cohort %in% c(Time1:Time2)), ]
#           f_i_sp_d$Species <- sp
#           f_i_sp_d$variable <- dim
#           f_i_sp_d$Fmax <- Fmax
#           f_i_sp_d$Fstart <- Fstart
#           
#           # harmonic mean
#           data_Ne_Nalpha = NULL
#           vector_harmonic_mean = c("Ne_N_alpha", "Numerator_N_alpha", 
#                                    "Ne_N_alpha/Numerator_N_alpha")
#           vector_harmonic_mean <- c(vector_harmonic_mean, improve_vriable_names(vector_harmonic_mean))
#           if (any(f_i_sp_d$variable %in% vector_harmonic_mean & "cohort" %in% colnames(f_i_sp_d))){
#             data_Ne_Nalpha <- f_i_sp_d [which(f_i_sp_d$variable %in% vector_harmonic_mean), ]
#             f_i_sp_d   <- f_i_sp_d [-which(f_i_sp_d$variable %in% vector_harmonic_mean), ]
#             data_Ne_Nalpha <- data_Ne_Nalpha %>%
#               group_by(Fmax, Fstart, variable, replicate, Species) %>%
#               summarise(value = psych::harmonic.mean(value, na.rm = TRUE))
#           }
#           # arithmetic mean mean
#           f_i_sp_d <- f_i_sp_d %>%
#             group_by(Fmax, Fstart, variable, replicate, Species) %>%
#             summarise(value = mean(value, na.rm = TRUE))
#           # bind arithmetic and harmonic mean dataset
#           if (!is.null(data_Ne_Nalpha)) {
#             f_i_sp_d = rbind(data_Ne_Nalpha,  f_i_sp_d)
#           }
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
#           f_i_sp_d <- f_i_sp_d[which(f_i_sp_d$cohort %in% c(Time1:Time2)), ]
#           f_i_sp_d$Species <- sp
#           f_i_sp_d$variable <- dim
#           f_i_sp_d$Fmax <- stringr::str_split(stringr::str_split(f_i, pattern="Fmax")[[1]][2], pattern="_Fstart")[[1]][1]
#           f_i_sp_d$Fstart <- stringr::str_split(stringr::str_split(f_i, pattern="Fstart")[[1]][2], pattern="_sp")[[1]][1]
#           
#           
#           # mean of the cohorts : arithmetic or harmonic
#           # harmonic mean
#           data_Ne_Nalpha = NULL
#           vector_harmonic_mean = c("Ne_N_alpha", "Numerator_N_alpha", 
#                                    "Ne_N_alpha/Numerator_N_alpha")
#           vector_harmonic_mean <- c(vector_harmonic_mean, improve_vriable_names(vector_harmonic_mean))
#           if (any(f_i_sp_d$variable %in% vector_harmonic_mean & "cohort" %in% colnames(f_i_sp_d))){
#             data_Ne_Nalpha <- f_i_sp_d [which(f_i_sp_d$variable %in% vector_harmonic_mean), ]
#             f_i_sp_d   <- f_i_sp_d [-which(f_i_sp_d$variable %in% vector_harmonic_mean), ]
#             data_Ne_Nalpha <- data_Ne_Nalpha %>%
#               group_by(age, Fmax, Fstart, variable, replicate, Species) %>%
#               summarise(value = psych::harmonic.mean(value, na.rm = TRUE))
#           }
#           # bind arithmetic and harmonic mean dataset
#           if (!is.null(data_Ne_Nalpha)) {
#             data_Ne_Nalpha$cohort <- 1
#             df_trait_sp = rbind(data_Ne_Nalpha, df_trait_sp)
#           }
#           # arithmetic mean mean
#           f_i_sp_d <- f_i_sp_d %>%
#             group_by(age, Fmax, Fstart, variable, replicate, Species) %>%
#             summarise(value = mean(value, na.rm = TRUE))
#           
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
#   if (file.exists(paste0(path, "/", name_out, sp_i, ".RData")) && FORCE !=T){
#     next
#   }else{
#     
#     # produce a dataset without age dimension : weighted mean
#     # weight_df <- full_dataset_age[[sp_i]] %>%
#     #   dplyr::filter(variable == "abundanceByAge") %>%
#     #   dplyr::select(Species, replicate, Fstart, Fmax, cohort, age, weight_value = value)
#     # data_with_w <- full_dataset_age[[sp_i]] %>%
#     #   dplyr::left_join(weight_df, by = c("Species", "replicate", "cohort", "age", "Fmax", "Fstart"), relationship = "many-to-many")
#     # data_with_w <- data_with_w[-which(data_with_w$weight_value == 0), ]
#     # data_with_w$age <- as.numeric(data_with_w$age)
#     data_out <- full_dataset_age[[sp_i]] %>%
#       dplyr::group_by(Species, replicate, Fstart, Fmax, cohort, variable) %>%
#       dplyr::summarise(
#         value = mean(value, na.rm=T), #stats::weighted.mean(value, w = weight_value, na.rm = TRUE),
#         .groups = "drop"
#       )
#     final_dataset <- rbind(full_dataset[[sp_i]], data_out)
#     
#     # produce a dataste with the age dimension
#     # full_dataset_add_age <- lapply(names(table(full_dataset_age[[sp_i]]$age)), FUN = function(x){
#     #   out = as.data.frame(full_dataset_age[[sp_i]])
#     #   out$age = x
#     #   return(out)
#     # })
#     # full_dataset_add_age <- do.call(rbind, full_dataset_age)
#     
#     final_dataset_age <- (full_dataset_age[[sp_i]])
#     
#     # save all
#     complete = list(final_dataset, final_dataset_age)
#     
#     save(complete, file=paste0(path, "/", name_out, sp_i, ".RData"))
#   }
# }
# 
