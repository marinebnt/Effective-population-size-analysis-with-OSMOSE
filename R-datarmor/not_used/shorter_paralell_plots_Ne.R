library(future)
library(future.apply)
library(progressr)


plan(multisession, workers = 15) #parallel::detectCores() - 1)
handlers("txtprogressbar") 
SP = c("AtlanticHerring",  "AtlanticMackerel", "Sandeel",  "EuropeanSprat", "NorwayPout",       
       "EuropeanPlaice", "CommonSole", "Saithe", "AtlanticCod", "Haddock",  "HorseMackerel",
       "Whiting",  "CommonDab", "GreyGurnard", "Hake")
main = "New_gpu_conf.R"
source("osmose-ns-v4.4.1_clim/R-datarmor/plot_shorter_all_Ne2.R")
source("osmose-ns-v4.4.1_clim/R/functions/source.R")

with_progress({
  p <- progressor(along = SP)    # create a progressor function with length = length(SP)
  res <- future_lapply(SP, function(sp) {
    Ne_variables_calculations_and_plots(sp)             # your function
  })
})
