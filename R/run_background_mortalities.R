
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
path <- c("C:/Users/mbeneat/Documents/osmose/Chap2/osmose-ns/osmose-ns-v4.4.1_clim/output/PLOTS/FINAL_Ne_Htz_mature/")
PATHOUT  = "osmose-ns-v4.4.1_clim/output/PLOTS/KEY_TRAITS_MATURE/"
if (!dir.exists(PATHOUT)) dir.create(PATHOUT, recursive = TRUE)
pathoutputRDATA = paste0(PATHOUT, "/RDATA_MATURE")
if(!dir.exists(pathoutputRDATA)) dir.create(pathoutputRDATA, recursive = T)
base_dir <- paste0(path, "INIT")
folders  <- list.dirs(base_dir, recursive = F) # files INIT2 with raw Ne datasets

folders_traits_y <- "osmose-ns-v4.4.1_clim/output/PLOTS/Trait_variables_by_year_BYSP//RData/"
folders_traits_y <- list.dirs(folders_traits_y, recursive = F)
folders_traits_c <- "osmose-ns-v4.4.1_clim/output/PLOTS/Trait_variables_by_cohort/RData/"
folders_traits_c  <- list.dirs(folders_traits_c, recursive = F)
folders_traits_y_a <- "osmose-ns-v4.4.1_clim/output/PLOTS/Trait_variables_by_year_by_age2/RData/"
folders_traits_y_a <- list.dirs(folders_traits_y_a, recursive = F)
folders_traits_c_a <- "osmose-ns-v4.4.1_clim/output/PLOTS/Trait_variables_by_cohort_by_age2/RData/"
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
sp_i = "Hake"
FORCE=F
full_dataset = list()
full_dataset_age = list()
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
base_path = "osmose-ns-v4.4.1_clim/output/PLOTS/Trait_variables_by_cohort_by_age2"
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


