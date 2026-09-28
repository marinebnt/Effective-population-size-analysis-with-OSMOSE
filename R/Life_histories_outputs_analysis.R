# the idea is to store separately all the key parameters per species into plots
# and to produce one pdf sheet per species with all the key variables to have it all summarized at one glance




#### Initialisation --------------
library(legendry)
library(ggpubr)

# functions and environment
main = "New_gpu_conf.R" # config file
source("R/functions/source.R") # load config and key parameters and functions
source("R/functions/function_data_traits_summary.R") # datasets preparation & loading
source("R/functions/functions_key_traits.R") # plots 

# paths
# path <- paste0(simu_repo, "/output/PLOTS/FINAL_Ne_Htz_mature_short/")
PATHOUT  = "output/PLOTS/KEY_TRAITS_MATURE/"
if (!dir.exists(PATHOUT)) dir.create(PATHOUT, recursive = TRUE)
if (!dir.exists(paste0(PATHOUT, "/Clean_for_paper"))) dir.create(paste0(PATHOUT, "/Clean_for_paper"), recursive = TRUE)
pathoutputRDATA = paste0(PATHOUT, "RDATA_MATURE")
if(!dir.exists(pathoutputRDATA)) dir.create(pathoutputRDATA, recursive = T)
base_path = paste0("output/DATASETS/Trait_variables_by_cohort_by_age2")
base_path = paste0("output/DATASETS/Trait_variables_by_cohort_by_age_withWbyAge/")
# folders_long <- list.dirs(paste0(path_long, "INIT"), recursive = F)

folders_traits_c <- paste0("output/DATASETS/Trait_variables_by_cohort/RData/")
folders_traits_c  <- list.dirs(folders_traits_c, recursive = F)
folders_traits_c_a <- paste0("output/DATASETS/Trait_variables_by_cohort_by_age/RData/")
folders_traits_c_a <- paste0("output/DATASETS/Trait_variables_by_cohort_by_age_withWbyAge/RData/")
folders_traits_c_a  <- list.dirs(folders_traits_c_a, recursive = F)

#### DATASETS load or format dataset ------------------------
Time1=150
Time2=160
species_df=species
sp_i = "Hake"
FORCE=T
dataset="data_Ne"

## % activate if needed : will produce .RData files in the folder paste0(simu_repo, "/output/PLOTS/FINAL_Ne_Htz_mature/")
## % 1 file per species, with the dataset at the scale of the cohorts, starting at age = alpha (first mature age)
# dataset_Ne_mean_cohorts(species, FORCE, pathdata=paste0("output/DATASETS/FINAL_Ne_Htz_mature/"),
#         pathout=pathoutputRDATA, Time1, Time2)

## loading the dataset of traits averaged over cohorts Time1 to Time2, and reorganised to remove age=0, add reassign natural mortality as add mortality, and estimate natural mortality
data_complete_c <- load_trait_df(base_path = base_path, pathout=pathoutputRDATA)


# dataset by species 
sp="AtlanticHerring"
## load and filter data
load(file=paste0(pathoutputRDATA, "/Complete_INIT2", sp, ".RData"))
final_dataset_age = complete[[2]]
min_age <- min(as.numeric(complete[[2]]$age), na.rm=T)
rm(complete)










