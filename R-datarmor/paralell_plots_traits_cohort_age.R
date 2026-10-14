library(future)
library(future.apply)
library(progressr)

plan(multisession, workers = 15) #parallel::detectCores() - 1)
handlers("txtprogressbar") 
SP = c("AtlanticHerring",  "AtlanticMackerel", "Sandeel",  "EuropeanSprat", "NorwayPout",       
       "EuropeanPlaice", "CommonSole", "Saithe", "AtlanticCod", "Haddock",  "HorseMackerel",
       "Whiting",  "CommonDab", "GreyGurnard", "Hake")
main = "New_gpu_conf.R"
source("R-datarmor/PLOT_trait_replicates_cohorts_age.R")
source("R/functions/source.R")


# path_in <- c("C:/Users/mbeneat/Documents/osmose/Ne_on_NT_with_F/osmose-ns/osmose-ns-v4.4.1_clim/output/old/test")
# path_in <- c("Y:/outputs-meso-marine")
path_in <- c("/home/datawork-marbec-scenlab/OSMOSE/Marine/SIMU_Ne_on_NT_with_climate_change_FIN")
path_in = c("/home1/datawork/mbeneat/SIMU_Ne_on_NT_with_F_FIN")
simu_name = "Ne_on_NT_with_climate_change_enet"
simu_name = "Ne_on_NT_with_F_enet"
with_progress({
  p <- progressor(along = SP)    # create a progressor function with length = length(SP)
  res <- future_lapply(SP, function(sp) {
    PLOT_trait_replicates2(sp, path_in=path_in, simu_name=simu_name)    # your function
  })
})
