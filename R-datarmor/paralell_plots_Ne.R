library(future)
library(future.apply)
library(progressr)


plan(multisession, workers = 15) #parallel::detectCores() - 1)
handlers("txtprogressbar") 
SP = c("AtlanticHerring",  "AtlanticMackerel", "Sandeel",  "EuropeanSprat", "NorwayPout",       
       "EuropeanPlaice", "CommonSole", "Saithe", "AtlanticCod", "Haddock",  "HorseMackerel",
       "Whiting",  "CommonDab", "GreyGurnard", "Hake")
main = "New_gpu_conf_RCP.R"
source("R-datarmor/plot_all_Ne2.R")
source("R/functions/source.R")

with_progress({
  p <- progressor(along = SP)    # create a progressor function with length = length(SP)
  res <- future_lapply(SP, function(sp) {
    plot_all_Ne2(sp)             # your function
  })
})
