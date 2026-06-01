

main = "New_gpu_conf.R" # config file
source("osmose-ns-v4.4.1_clim/R/functions/source.R") # load config and key parameters and functions
source("osmose-ns-v4.4.1_clim/R/functions/functions_key_traits.R") # plots against Fmax
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
data_complete_c <- data_complete_c[-which(data_complete_c$age == 0), ]

PATHOUT  = "osmose-ns-v4.4.1_clim/output/PLOTS/KEY_TRAITS/"


# adding variables or transforming variables
data_complete_supp = NULL
for (sp in species){
  cat(sp)
  
  df_trait_sp <- data_complete_c %>% filter(Species == sp)
  df_trait_sp <- df_trait_sp %>% filter(Fstart == 1)
  
  df_trait_t    <- estimateSum(trait_data_full_i = df_trait_sp, var=paste0("mean_of_F_per_age", "\nsum", "mean_of_Mnat_per_age"))
  df_trait_t$variable    <- "sum_mortality_rates"
  df_trait_sp   <- rbind(df_trait_sp, df_trait_t)
  
  df_trait_t    <- estimateSurv(trait_data_full_i = df_trait_sp)
  df_trait_sp   <- rbind(df_trait_sp, df_trait_t)
  
  df_trait_t    <- estimateDiv(trait_data_full_i = df_trait_sp, var="mean_of_eggNumber_per_age\nDivmean_of_abundanceByAge_per_age")
  df_trait_sp   <- rbind(df_trait_sp, df_trait_t)
  
  df_trait_t    <- estimateDiv(trait_data_full_i = df_trait_sp, var="mean_of_Mpred_per_age\nDivsum_mortality_rates")
  df_trait_sp   <- rbind(df_trait_sp, df_trait_t)
  
  var_i    <- unique(df_trait_sp$variable)
  # compare all mortality rates to their value *abundance to have the number of individuals to understand where is the increase coming from
  for (v in var_i){
    if (any(grepl(pattern="mean_of_M", v) | grepl(pattern="mean_of_F", v) | grepl(pattern="sum_mortality_rates", v))) {
      df_trait_t    <- estimateMnumber(trait_data_full_i = df_trait_sp, var=paste0(v))
      df_trait_sp   <- rbind(df_trait_sp, df_trait_t)
    }
  }
  var_i    <- unique(df_trait_sp$variable)
  for (v in var_i){
    if (any(grepl(pattern="deaths", v))) {
      df_trait_t    <- estimatecolSum(trait_data_full_i = df_trait_sp, var=paste0("colSum", v))
      df_trait_sp   <- rbind(df_trait_sp, df_trait_t)
    }
  }
  
  if (is.null(data_complete_supp)) {
    data_complete_supp <- df_trait_sp
  } else(data_complete_supp   <- rbind(data_complete_supp, df_trait_sp))
}

# Build a parallel dataset with slopes of each variable against Fmax.
# One slope is estimated for each Species x age x time x replicate x Fstart x variable group.
safe_lm_coef <- function(x, y, which_coef = c("slope", "intercept")) {
  which_coef <- match.arg(which_coef)
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]
  y <- y[ok]
  if (length(x) < 2 || length(unique(x)) < 2) return(NA_real_)
  fit <- lm(y ~ x)
  if (which_coef == "slope") return(unname(coef(fit)[2]))
  return(unname(coef(fit)[1]))
}

# dataset_slope_long <- data_complete_supp %>%
#   group_by(Species, age, replicate, Fstart, variable) %>%
#   summarize(
#     slope_Fmax = safe_lm_coef(Fmax, value, "slope"),
#     intercept_Fmax = safe_lm_coef(Fmax, value, "intercept"),
#     n_Fmax = n_distinct(Fmax[is.finite(value)]),
#     value_range = {
#       x <- value[is.finite(value)]
#       if (length(x) == 0) NA_real_ else diff(range(x))
#     },
#     .groups = "drop"
#   ) %>%
#   filter(n_Fmax > 1, is.finite(slope_Fmax), value_range > 0)
# 
# dataset_slope <- dataset_slope_long %>%
#   select(Species, age, replicate, Fstart, variable, slope_Fmax) %>%
#   pivot_wider(names_from = "variable", values_from = "slope_Fmax")

dataset_slope$`Ne/N` <- ifelse(dataset_slope$Species %in% 
                               c("Hake", "Whiting", "AtlanticHerring", "AtlanticCod"), "increasing",
                             ifelse(dataset_slope$Species %in% 
                                      c("EuropeanSprat"), "increasingCST",
                                    ifelse(dataset_slope$Species %in% 
                                             c("AtlanticMackerel", "Saithe", "CommonDab", "GreyGurnard", "Sandeel", 
                                               "CommonSole"), "decreasing",
                                           ifelse(dataset_slope$Species %in% 
                                                    c("EuropeanPlaice", "NorwayPout"), "decreasingCST", "CST"))))
dataset_slope_long$`Ne/N` <- ifelse(dataset_slope_long$Species %in% 
                               c("Hake", "Whiting", "AtlanticHerring", "AtlanticCod"), "increasing",
                             ifelse(dataset_slope_long$Species %in% 
                                      c("EuropeanSprat"), "increasingCST",
                                    ifelse(dataset_slope_long$Species %in% 
                                             c("AtlanticMackerel", "Saithe", "CommonDab", "GreyGurnard", "Sandeel", 
                                               "CommonSole"), "decreasing",
                                           ifelse(dataset_slope_long$Species %in% 
                                                    c("EuropeanPlaice", "NorwayPout"), "decreasingCST", "CST"))))

base_path = "osmose-ns-v4.4.1_clim/output/PLOTS/RDATA/"
data_slope <- list(dataset_slope_long, dataset_slope)
# save(data_slope, file=paste0(base_path, "/traitSLOPE.RData"))
load(file=paste0(base_path, "/traitSLOPE.RData"))

dataset_slope_long <- data_slope[[1]]
dataset_slope      <- data_slope[[1]]


# finish to clean dataset
data_w <- pivot_wider(data_complete_supp, names_from = "variable", values_from = "value")
cat("removing the NANs, don't know if good move")
#normalize
# data_w_norm <- data_w[, c("mean_of_F_per_age", "mean_of_Mpred_per_age", "mean_of_Mfor_per_age", "mean_of_Mstar_per_age")]
data_w_norm <- data_w
data_w_norm <- data_w_norm[, -which(colnames(data_w_norm) == c("colSumMpred\nDivsum_mortality_rates_deaths"))]
colSums(is.na(data_w_norm))
# data_w_norm[is.na(data_w_norm)] <- 0
data_w_norm <- na.omit(data_w_norm)
# data_w_norm[, -c(1:5)] =  scale(data_w_norm[, -c(1:5)]) #not scaling anymore
data_w_norm[is.na(data_w_norm)] <- 0
# data_w[, c("mean_of_F_per_age", "mean_of_Mpred_per_age", "mean_of_Mfor_per_age", "mean_of_Mstar_per_age")] =
#   scale(data_w[, c("mean_of_F_per_age", "mean_of_Mpred_per_age", "mean_of_Mfor_per_age", "mean_of_Mstar_per_age")])
data_w_norm$`Ne/N` <- ifelse(data_w_norm$Species %in% 
                               c("Hake", "Whiting", "AtlanticHerring", "AtlanticCod"), "increasing",
                             ifelse(data_w_norm$Species %in% 
                                      c("EuropeanSprat"), "increasingCST",
                                    ifelse(data_w_norm$Species %in% 
                                             c("AtlanticMackerel", "Saithe", "CommonDab", "GreyGurnard", "Sandeel", 
                                               "CommonSole"), "decreasing",
                                           ifelse(data_w_norm$Species %in% 
                                                    c("EuropeanPlaice", "NorwayPout"), "decreasingCST", "CST"))))

# keeping only Fmax == 1 and Fstart ==1, because we are trying to characterize the species, not to see the consequences of fishing 
# data_w_norm <- data_w_norm[data_w_norm$Fmax==1, ]
# data_w_norm <- data_w_norm[data_w_norm$Fstart==1, ]
# data_w_norm <- data_w_norm[, -which(colnames(data_w_norm) %in% c("Fstart", "Fmax"))]








ggplot(data_w_norm, aes(x=mean_of_biomassByAge_per_age, y=surv_curve, color=.data[["Ne/N"]]=="increasingCST"))+
  geom_point(alpha=ifelse(data_w_norm[,"Ne/N"]=="increasingCST", 1, 0.01))+
  facet_wrap(~Fmax)
ggplot(data_w_norm, aes(x=mean_of_Mpred_per_age, y=sum_mortality_rates, color=.data[["Ne/N"]]=="CST"))+
  geom_point(alpha=ifelse(data_w_norm[,"Ne/N"]=="CST", 1, 0.01))+
  facet_wrap(~Fmax)
ggplot(data_w_norm, aes(x=mean_of_Mpred_per_age, y=sum_mortality_rates, color=.data[["Ne/N"]]))+
  geom_point()+
  facet_wrap(~Fmax)
ggplot(data_w_norm, aes(x=mean_of_Mpred_per_age, y=all_deaths, color=.data[["Ne/N"]]))+
  geom_point(alpha=ifelse(data_w_norm[,"Ne/N"]=="increasingCST", 1, 0.01))+
  facet_wrap(~age)
ggplot(data_w_norm, aes(x=mean_of_Mpred_per_age, y=sum_mortality_rates, color=.data[["Ne/N"]]))+
  geom_point()+
  facet_wrap(~Fmax)
ggplot(data_w_norm, aes(x=.data[["mean_of_Mpred_per_age"]]/.data[["sum_mortality_rates"]], y=mean_of_Mpred_per_age, color=.data[["Ne/N"]]))+
  geom_point(alpha=ifelse(data_w_norm[,"Ne/N"]=="increasingCST", 1, 0.01))+
  facet_wrap(~Fmax)
ggplot(data_w_norm, aes(x=.data[["mean_of_F_per_age"]]/.data[["sum_mortality_rates"]], y=mean_of_Mpred_per_age, color=.data[["Ne/N"]]))+
  geom_point(alpha=ifelse(data_w_norm[,"Ne/N"]=="increasingCST", 1, 0.01))+
  facet_wrap(~Fmax)
ggplot(data_w_norm, aes(x=log(.data[["colSumMpred_deaths"]]/.data[["colSumall_deaths"]]), y=mean_of_Mpred_per_age, color=.data[["Ne/N"]]))+
  geom_point()+
  facet_wrap(~Fmax)

ggplot(data_w_norm, aes(x=log((.data[["sum_mortality_rates"]])/(.data[["mean_of_Mpred_per_age"]])), y=, color=.data[["Ne/N"]]))+
  geom_point()+
  facet_wrap(~Fmax)

data_w_norm$`logMpred/Mtot` <- log(data_w_norm[["mean_of_Madd_per_age"]]/data_w_norm[["sum_mortality_rates"]])
ggplot(data_w_norm, aes(x=`logMpred/Mtot`, y=mean_of_abundanceByAge_per_age, color=.data[["Ne/N"]]))+
  geom_point()+
  facet_wrap(~Fmax)


data_w_norm$`logMadd/Mtot` <- log(data_w_norm[["mean_of_Madd_per_age"]]/data_w_norm[["sum_mortality_rates"]])
data_w_norm$`logMnat/Mtot` <- log(data_w_norm[["mean_of_Mnat_per_age"]]/data_w_norm[["sum_mortality_rates"]])
data_w_norm$`logMpred/Mtot` <- log(data_w_norm[["mean_of_Mpred_per_age"]]/data_w_norm[["sum_mortality_rates"]])
data_w_norm$`logF/Mtot` <- log(data_w_norm[["mean_of_F_per_age"]]/data_w_norm[["sum_mortality_rates"]])
data_w_norm$`logMstar/Mtot` <- log(data_w_norm[["mean_of_Mstar_per_age"]]/data_w_norm[["sum_mortality_rates"]])
data_w_norm$`logMpred/F` <- log(data_w_norm[["mean_of_Mpred_per_age"]]/data_w_norm[["mean_of_F_per_age"]])
data_w_norm$`logMpred/Mnat` <- log(data_w_norm[["mean_of_Mpred_per_age"]]/data_w_norm[["mean_of_Mnat_per_age"]])
data_w_norm$`logMadd/Mtot` <- log(data_w_norm[["mean_of_Madd_per_age"]]/data_w_norm[["sum_mortality_rates"]])
data_w_norm$`logMadd/Mnat` <- log(data_w_norm[["mean_of_Madd_per_age"]]/data_w_norm[["mean_of_Mnat_per_age"]])
data_w_norm$`logMadd/F` <- log(data_w_norm[["mean_of_Madd_per_age"]]/data_w_norm[["mean_of_F_per_age"]])
data_w_norm$sp_shape <- ifelse(data_w_norm[,"Ne/N"]=="CST",(data_w_norm$Species), "")
ggplot(data_w_norm, aes(y=`logMnat/Mtot`, x=log(.data[["mean_of_eggNumber_per_age
Divmean_of_abundanceByAge_per_age"]]), color=.data[["Ne/N"]]))+
  geom_point(alpha=ifelse(data_w_norm[,"Ne/N"]=="CST", 1, 0.1), 
             aes(shape=sp_shape))+
  facet_wrap(~age)

ggplot(data_w_norm, aes(y=, x=`logMnat/Mtot`, color=.data[["Ne/N"]]))+
  geom_point(#alpha=ifelse(data_w_norm[,"Ne/N"]=="CST", 1, 0.1),
             aes(shape=sp_shape))+
  facet_wrap(~age)

ggplot(data_w_norm,aes(y=(.data[["mean_of_biomassByAge_per_age"]]), x=mean_of_abundanceByAge_per_age, color=.data[["Ne/N"]]))+
  geom_point(alpha=ifelse(data_w_norm[,"Ne/N"]=="increasingCST", 1, 0.1),
    aes(shape=sp_shape))+
  facet_wrap(~age)

#WITH THE SLOPE DATASET ------------------

dataset_slope$`logMadd/Mtot` <- log(dataset_slope[["mean_of_Madd_per_age"]]/dataset_slope[["sum_mortality_rates"]])
dataset_slope$`logMnat/Mtot` <- log(dataset_slope[["mean_of_Mnat_per_age"]]/dataset_slope[["sum_mortality_rates"]])
dataset_slope$`logMpred/Mtot` <- log(dataset_slope[["mean_of_Mpred_per_age"]]/dataset_slope[["sum_mortality_rates"]])
dataset_slope$`logF/Mtot` <- log(dataset_slope[["mean_of_F_per_age"]]/dataset_slope[["sum_mortality_rates"]])
dataset_slope$`logMstar/Mtot` <- log(dataset_slope[["mean_of_Mstar_per_age"]]/dataset_slope[["sum_mortality_rates"]])
dataset_slope$`logMfor/Mtot` <- log(dataset_slope[["mean_of_Mfor_per_age"]]/dataset_slope[["sum_mortality_rates"]])


ggplot(dataset_slope, aes(y=mean_of_Mstar_per_age, x=surv_curve, color=.data[["Ne/N"]]))+
  geom_point(#alpha=ifelse(data_w_norm[,"Ne/N"]=="CST", 1, 0.1),
             #  aes(shape=sp_shape)
             )+
  facet_wrap(~age)

ggplot(dataset_slope, aes(y=(.data[["mean_of_biomassByAge_per_age"]]), x=mean_of_abundanceByAge_per_age, color=.data[["Ne/N"]]))+
  geom_point(#alpha=ifelse(data_w_norm[,"Ne/N"]=="CST", 1, 0.1),
    #  aes(shape=sp_shape)
  )+
  facet_wrap(~age)


ggplot(dataset_slope, aes(y=(.data[["mean_of_biomassByAge_per_age"]]), x=mean_of_abundanceByAge_per_age, color=.data[["Ne/N"]]))+
  geom_point(#alpha=ifelse(data_w_norm[,"Ne/N"]=="CST", 1, 0.1),
    #  aes(shape=sp_shape)
  )+
  facet_wrap(~age)







## SUMMARY plot


# incresing constant : highest young age Mstar
ggplot(data_w_norm,aes(x=Fmax, y=mean_of_Mstar_per_age, 
                          color=.data[["Ne/N"]]))+
  geom_point(alpha=ifelse(data_w_norm[,"Ne/N"]=="increasingCST", 1, 0.1))+
  facet_wrap(~age)

# decreasing constant : highest log(Mpred/Mtot) + longest lifespan
ggplot(data_w_norm,aes(x=age, y=`logMpred/Mtot`, 
                          color=.data[["Ne/N"]]))+
  geom_point(alpha=ifelse(data_w_norm[,"Ne/N"]=="decreasingCST", 1, 0.1))+
  facet_wrap(~Fmax)
ggplot(data_w_norm,aes(x=age, y=`logF/Mtot`, 
                       color=.data[["Ne/N"]]))+
  geom_point(alpha=ifelse(data_w_norm[,"Ne/N"]=="decreasingCST", 1, 0.1))+
  facet_wrap(~Fmax)
ggplot(data_w_norm,aes(x=age, y=`logMstar/Mtot`, 
                       color=.data[["Ne/N"]]))+
  geom_point(alpha=ifelse(data_w_norm[,"Ne/N"]=="decreasingCST", 1, 0.1))+
  facet_wrap(~Fmax)
ggplot(data_w_norm,aes(x=age, y=`logMpred/F`, 
                       color=.data[["Ne/N"]]))+
  geom_point(alpha=ifelse(data_w_norm[,"Ne/N"]=="decreasingCST", 1, 0.1))+
  facet_wrap(~Fmax)
ggplot(data_w_norm,aes(x=age, y=`logMpred/Mnat`, 
                       color=.data[["Ne/N"]]))+
  geom_point(alpha=ifelse(data_w_norm[,"Ne/N"]=="decreasingCST", 1, 0.1))+
  facet_wrap(~Fmax)


# constant : highest logMadd/Mtot
ggplot(data_w_norm,aes(x=age, y=`logMadd/Mtot`, 
                       color=.data[["Ne/N"]]))+
  geom_point(alpha=ifelse(data_w_norm[,"Ne/N"]=="CST", 1, 0.1))+
  facet_wrap(~Fmax)
ggplot(data_w_norm,aes(x=age, y=`logMnat/Mtot`, 
                       color=.data[["Ne/N"]]))+
  geom_point(alpha=ifelse(data_w_norm[,"Ne/N"]=="CST", 1, 0.1))+
  facet_wrap(~Fmax)

# increasing : lowest Mpred
ggplot(data_w_norm,aes(x=Fmax, y=1/(.data[["mean_of_Mpred_per_age"]]+0.001), 
                       color=.data[["Ne/N"]]))+
  geom_point(alpha=ifelse(data_w_norm[,"Ne/N"]=="increasing", 1, 0.1)
             )+
  facet_wrap(~age)




