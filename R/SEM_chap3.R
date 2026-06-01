##  Marine Beneat
##  10/05/26
## here we are trying to explain the specific NE/NT ratio in response to fishing pressure
## species are described by their functional traits 
## SEM models : Lefcheck 2016 (piecewise VS lavaan : Grace 2022)
## 


# scientific question 
# can we explain the Ne/N variability between species by using the species functional traits ? 
# species functional traits are depending on age, replicates and Fmax
# it could also be seen as a multivariate problem, when considering Ne and N simultaneously vary with the same functional traits


# PRACTICE ON 1 SP -----------------
# Collect dataset 

library("osmose")
library("piecewiseSEM")

# functions and environment
main = "New_gpu_conf.R" # config file
source("R/functions/source.R") # load config and key parameters and functions
source("R/functions/function_data_traits_summary.R") # datasets preparation & loading
source("R/functions/functions_key_traits.R") # plots 

# paths
simu_repo <- "C:/Users/mbeneat/Documents/osmose/Chap2/osmose-ns/osmose-ns-v4.4.1_clim/"
path <- paste0(simu_repo, "/output/PLOTS/FINAL_Ne_Htz_mature/")
# path <- paste0(simu_repo, "/output/PLOTS/FINAL_Ne_Htz_mature_short/")
PATHOUT  = "output/PLOTS/KEY_TRAITS_MATURE/"
if (!dir.exists(PATHOUT)) dir.create(PATHOUT, recursive = TRUE)
pathoutputRDATA = paste0(PATHOUT, "/RDATA_MATURE")
if(!dir.exists(pathoutputRDATA)) dir.create(pathoutputRDATA, recursive = T)
base_dir <- paste0(path, "INIT")
base_path = paste0(simu_repo, "output/PLOTS/Trait_variables_by_cohort_by_age2")
folders  <- list.dirs(base_dir, recursive = F) # files INIT2 with raw Ne datasets
# folders_long <- list.dirs(paste0(path_long, "INIT"), recursive = F)

#### DATASETS load or format dataset ------------------------
Time1=150
Time2=160
species_df=species
sp_i = "Hake"
FORCE=T
dataset="data_Ne"
folders <- list(folders) #, folders_long)
path = list(path) #, path_long)
path = path[[ifelse(dataset=="data_Ne", 1, 2)]]

## % activate if needed : will produce .RData files in the folder paste0(simu_repo, "/output/PLOTS/FINAL_Ne_Htz_mature/")
## % 1 file per species, with the dataset at the scale of the cohorts, starting at age = alpha (first mature age)
#dataset_Ne_mean_cohorts(species, FORCE, folders, path, Time1, Time2)

## loading the dataset of traits averaged over cohorts Time1 to Time2, and reorganised to remove age=0, add reassign natural mortality as add mortality, and estimate natural mortality
pathout   = paste0(base_path, "/otherplots")
if (!dir.exists(pathout)) dir.create(pathout, recursive = T)
data_complete_c <- load_trait_df(base_path)
  
sp = "AtlanticHerring"
sp2 = "AtlanticMackerel"
complete2 <- get(load(file=paste0(path, "/Complete_INIT2", sp2, ".RData")))[["complete_with_age"]]
complete <- get(load(file=paste0(path, "/Complete_INIT2", sp, ".RData")))[["complete_with_age"]]
min_age <- min(as.numeric(complete$age), na.rm=T)
complete$age[which(is.na(complete$age))] <- min_age
min_age2 <- min(as.numeric(complete2$age), na.rm=T)
complete2$age[which(is.na(complete2$age))] <- min_age2
complete <- rbind(complete, complete2)
for (var in unique(complete$variable)){
  if (length(table(complete$age[which(complete$variable == var)])) ==1 ){
    df_var_i <- complete[which(complete$variable == var), ]
    df_else  <- complete[which(complete$variable != var), ]
    df_var = NULL
    for (age in unique(df_else$age)) {
      df_var_i$age <- age
      df_var <- rbind(df_var, df_var_i)
    }
    complete = rbind(df_var, df_else)
  }
}
complete$age <- as.numeric(complete$age)
complete <- pivot_wider(complete, names_from="variable", values_from="value")
complete[, "Ne_N_alphadivNT"] <- complete$Ne_N_alpha/complete$NT
complete[, "logDeltasquare"] <- log(complete$Deltasquare)
complete[, "logDead"] <- log(complete$Dead)
complete[, "logNe_N_alphadivNT"] <- log(complete$Ne_N_alphadivNT)
complete[, "logsurvivorship_curve"] <- log(complete$survivorship_curve)
complete$Species <- as.factor(complete$Species)
complete$age <- as.factor(complete$age)
complete$replicate <- as.factor(complete$replicate)
# Create SEM using `psem`
complete <- na.omit(complete)

completei <- complete[, c("age", "replicate", "Species", "b_rate_sc", "Dead", "survivorship_curve", "Deltasquare", "Ne_N_alphadivNT")]
modelList <- psem(
  glm(Dead ~ age + survivorship_curve + Species, "poisson", data=complete),
  glm(Deltasquare ~ age + Species + b_rate_sc, "poisson", data=complete),
  lm(Ne_N_alphadivNT ~ Dead + Deltasquare + Species, data=complete)
)
modelList <- psem(
   lm(logDead ~ age + logsurvivorship_curve + Species, data=complete),
   lm(Deltasquare ~ age + Species + b_rate_sc, data=complete),
   lm(Ne_N_alphadivNT ~ logDead+Deltasquare + Species, data=complete)
)
modelList <- psem(
  lm(logDead ~ age + logsurvivorship_curve + Species + b_rate_sc + Fmax, data=complete),
  lm(Deltasquare ~ age + logDead + Species + logsurvivorship_curve + b_rate_sc + Fmax, data=complete),
  lm(Ne_N_alphadivNT ~ age + logDead+Deltasquare + Species + b_rate_sc + logsurvivorship_curve + Fmax,  data=complete)
)
dSep(modelList) # https://jonlefcheck.net/2014/07/06/piecewise-structural-equation-modeling-in-ecological-research/ : what paths to add? 
residuals(modelList)
print(modelList)
anova(modelList) # can also be used to compare models
AIC_psem(modelList) 
plot(modelList)
# plot.sem(modelList)
# Run summary
summary(modelList)
# Address conflict using conserve = T
summary(modelList, conserve = T)
# Address conflict using direction = c()
summary(modelList, direction = c("y2 <- y1"))
# Address conflict using correlated errors
modelList2 <- update(modelList, y2 %~~% y1)
summary(modelList2)




# EXAMPLES -----------
# https://jslefche.github.io/piecewiseSEM/
# Load library
library(piecewiseSEM)

# Create fake data
set.seed(1)

data <- data.frame(
  x = runif(100),
  y1 = runif(100),
  y2 = rpois(100, 1),
  y3 = runif(100)
)

# Create SEM using `psem`
modelList <- psem(
  lm(y1 ~ x, data),
  glm(y2 ~ x, "poisson", data),
  lm(y3 ~ y1 + y2, data),
  data
)

# Run summary
summary(modelList)

# Address conflict using conserve = T
summary(modelList, conserve = T)

# Address conflict using direction = c()
summary(modelList, direction = c("y2 <- y1"))

# Address conflict using correlated errors
modelList2 <- update(modelList, y2 %~~% y1)

summary(modelList2)
