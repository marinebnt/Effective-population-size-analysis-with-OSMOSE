
main = "New_gpu_conf.R"
source("osmose-ns-v4.4.1_clim/R/functions/function_plot_trait_age_or_time.R")
source("osmose-ns-v4.4.1_clim/R/functions/source.R")
source("osmose-ns-v4.4.1_clim/R/functions/function_data_traits_summary.R")

base_dir <- "osmose-ns-v4.4.1_clim/output/PLOTS/All_clean_genetic_drift_indicators/INIT/"
folders  <- list.dirs(base_dir, recursive = F)
repo        = "Simu_LONG"
multip_name = "y"
multipliers = c(1,0.5,0.75,1.25,1.5)
ideal_order_factor <- paste0("Fmax=", rep(c(0.5,0.75,1,1.25,1.5), 5),
                             " Fstart=", c(rep(0.5, 5),rep(0.75, 5),rep(1, 5),rep(1.25, 5),rep(1.5, 5)))



trait_dim= c("eggNumber", "abundanceByAge", "biomassByAge", "mortalityByAge", "schoolNumber", "Mpred", "Mstar", "Mnat", "F", "Z", "Mfor", "Dis")


# path_in <- c("C:/Users/mbeneat/Documents/osmose/Chap2/osmose-ns/osmose-ns-v4.4.1_clim/output/test")
# path_in <- c("/home1/datawork/mbeneat/SIMU_CHAP2_FIN")
# list_d <- list.dirs(path = path_in, recursive = F, full.names = T)

F_vec <- list()



# BY COHORT -----------------------------------

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
# data_complete_c[grep(data_complete_c$TRAIT, pattern="fec"), "TRAIT"] <- "mean_of_eggNumber_per_age"
unique(data_complete_c$TRAIT)

# add the non raw traits 
Vk_surv <- get(load(file=paste0(configDir, "/output/RData/Vkx_surviv.RData")))
Vk_surv$var_F_for_sp <- NULL
data_complete_c <- rbind(data_complete_c, Vk_surv)
unique(data_complete_c$TRAIT)

# produce plots (saved directly in pdf format)
function_plot_trait_age_or_time(SP=species, age_or_cohort="age", paste0(pathout, "/Vk"), Fmax_constrain = NULL, Fstart_constrain = c(1), 
                                list_traits = list(c("mean_survivorship_curve"), c("mean_Vx")), 
                                data_complete=data_complete_c)
function_plot_trait_age_or_time(SP=species, age_or_cohort="age", paste0(pathout, "/Vk"), Fmax_constrain = NULL, Fstart_constrain = c(1), 
                                list_traits = list(c("sum_of_eggNumberByAge_per_age"), c("mean_Vx")), 
                                data_complete=data_complete_c)
function_plot_trait_age_or_time(SP=species, age_or_cohort="age", paste0(pathout, "/Vk"), Fmax_constrain = NULL, Fstart_constrain = c(1), 
                                list_traits = list(c("mean_of_fecundityByAge_per_age"), c("mean_Vx")), 
                                data_complete=data_complete_c)
function_plot_trait_age_or_time(SP=species, age_or_cohort="age", paste0(pathout, "/Vk"), Fmax_constrain = NULL, Fstart_constrain = c(1), 
                                list_traits = list(c("mean_of_fecundityByAge_per_age"), c("mean_Vx")), 
                                data_complete=data_complete_c)
function_plot_trait_age_or_time(SP=species, age_or_cohort="age", paste0(pathout, "/fec"), Fmax_constrain = NULL, Fstart_constrain = c(1), 
                                list_traits = list(c("mean_survivorship_curve"), c("mean_of_fecundityByAge_per_age")), 
                                data_complete=data_complete_c)
function_plot_trait_age_or_time(SP=species, age_or_cohort="age", paste0(pathout, "/fec"), Fmax_constrain = NULL, Fstart_constrain = c(1), 
                                list_traits = list(c("mean_of_Mnat_per_age"), c("mean_of_fecundityByAge_per_age")), 
                                data_complete=data_complete_c)
function_plot_trait_age_or_time(SP=species, age_or_cohort="age", paste0(pathout, "/fec"), Fmax_constrain = NULL, Fstart_constrain = c(1), 
                                list_traits = list(c("mean_of_Mstar_per_age"), c("mean_of_fecundityByAge_per_age")), 
                                data_complete=data_complete_c)
function_plot_trait_age_or_time(SP=species, age_or_cohort="age", paste0(pathout, "/fec"), Fmax_constrain = NULL, Fstart_constrain = c(1), 
                                list_traits = list(c("logsum_of_eggNumberByAge_per_age"), c("mean_of_fecundityByAge_per_age", "mean_of_Mnat_per_age", "mean_of_Mpred_per_age", "mean_survivorship_curve")), 
                                data_complete=data_complete_c)
function_plot_trait_age_or_time(SP=species, age_or_cohort="age", paste0(pathout, "/fec"), Fmax_constrain = NULL, Fstart_constrain = c(1), 
                                list_traits = list(c("logsum_of_eggNumberByAge_per_age"), c("mean_of_fecundityByAge_per_age", "mean_of_Mnat_per_age", "mean_of_Mpred_per_age", "mean_survivorship_curve")), 
                                data_complete=data_complete_c, relative_to_ref = T)
function_plot_trait_age_or_time(SP=species, age_or_cohort="age", paste0(pathout, "/fec"), Fmax_constrain = NULL, Fstart_constrain = c(1), 
                                list_traits = list(c("mean_of_fecundityByAge_per_age"), c("logsum_of_eggNumberByAge_per_age")), 
                                data_complete=data_complete_c, relative_to_ref = T)
function_plot_trait_age_or_time(SP=species, age_or_cohort="age", paste0(pathout, "/fec"), Fmax_constrain = NULL, Fstart_constrain = c(1), 
                                list_traits = list(c("logsum_of_eggNumberByAge_per_age"), c("mean_of_Mnat_per_age", "mean_of_Mpred_per_age", "mean_survivorship_curve")), 
                                data_complete=data_complete_c, relative_to_ref = T)
function_plot_trait_age_or_time(SP=species, age_or_cohort="age", paste0(pathout, "/fec"), Fmax_constrain = NULL, Fstart_constrain = c(1), 
                                list_traits = list(c("mean_of_sizeMature_per_age"), c("logsum_of_eggNumberByAge_per_age")), 
                                data_complete=data_complete_c, relative_to_ref = T)
function_plot_trait_age_or_time(SP=species, age_or_cohort="age", paste0(pathout, "/fec"), Fmax_constrain = NULL, Fstart_constrain = c(1), 
                                list_traits = list(c("sum_of_schoolNumber_per_age"), c("logsum_of_eggNumberByAge_per_age")), 
                                data_complete=data_complete_c, relative_to_ref = T)
function_plot_trait_age_or_time(SP=species, age_or_cohort="age", paste0(pathout, "/fec"), Fmax_constrain = NULL, Fstart_constrain = c(1), 
                                list_traits = list(c("sum_of_schoolNumber_per_age"), c("logsum_of_eggNumberByAge_per_age")), 
                                data_complete=data_complete_c, relative_to_ref = F)
function_plot_trait_age_or_time(SP=species, age_or_cohort="age", paste0(pathout, "/mort"), Fmax_constrain = NULL, Fstart_constrain = c(1), 
                                list_traits = list(c("sum_of_schoolNumber_per_age"), c("mean_of_Mnat_per_age")), 
                                data_complete=data_complete_c, relative_to_ref = F)
function_plot_trait_age_or_time(SP=species, age_or_cohort="age", paste0(pathout, "/mort"), Fmax_constrain = NULL, Fstart_constrain = c(1), 
                                list_traits = list(c("mean_of_Mpred_per_age"), c("sum_of_schoolNumber_per_age")), 
                                data_complete=data_complete_c, relative_to_ref = F)
function_plot_trait_age_or_time(SP=species, age_or_cohort="age", paste0(pathout), Fmax_constrain = NULL, Fstart_constrain = c(1), 
                                list_traits = list(c("mean_of_Mfor_per_age")), 
                                data_complete=data_complete_c, relative_to_ref = F)

function_plot_trait_age_or_time(SP=species, age_or_cohort="age", paste0(pathout), Fmax_constrain = NULL, Fstart_constrain = c(1), 
                                list_traits = list(c("mean_of_Mfor_per_age")), 
                                data_complete=data_complete_c, relative_to_ref = F)


