
rm(list=ls())
require("osmose")
require("mgcv")

setwd("C:/R/osmose-ns/")

filename = "Eco-OSMOSE-NS/Eco_Bioen_npv.R"

restart = TRUE

source("C:/R/osmose-ns/Eco-OSMOSE-NS/osmose-auxiliar.R")
source("C:/R/osmose-ns/Eco-OSMOSE-NS/osmose-fnpv.R")

output = F_npv(input.file = filename, Fmin = 0, Fmax=2, 
               simu=10, version="4.3.3", restart=restart)
  
simu = 10
title = seq.int(0, simu)
npvpersimu = read.csv(paste0("Eco-OSMOSE-NS/npv", "_after", simu, "replicat.csv"), col.names = title, check.names=FALSE)

# Find which simulation has the highest NPV
maxnpv = colnames(npvpersimu)[apply(npvpersimu,1,which.max)]
maxnpv

# Extract alpha for each species for the max NPV
maxalpha = read.csv(paste0("Eco-OSMOSE-NS/Fnpv/NPV_alpha_simu_",maxnpv, ".csv"))
maxalpha = maxalpha[,2:17]

maxalpha2 = t(maxalpha)
