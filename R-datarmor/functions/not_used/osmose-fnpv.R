# a <- readOsmoseConfiguration("C:/Users/amorell/Desktop/osmose-ns/Bioen-OSMOSE-NS - 4.3.3 - Morell et al. 2023/calib_config.csv")


#' Finds the F_npv value for a given species
#'
#' @param sp Index of the species from which F_npv is analysed  // I need to do on all species
#' @param input.file Path of the configuration file
#' @param restart True if Fnpv should start from restart
#' @param Fmin Minimum value of alpha
#' @param Fmax Maximum value of alpha
#' @param simu number of simulation needed
#' @param ...  Additional arguments of the run_osmose function
#'
#' @export
F_npv <- function(input.file, restart=TRUE, 
                  Fmin=0, Fmax=2, simu = 1, ...) {
  
  # Initial values
  input.folder = "C:/R/osmose-ns"
  jar.file = "C:/R/osmose-4.3.3/inst/java/osmose_4.3.3-jar-with-dependencies.jar"
  
  Param = readOsmoseConfiguration(input.file)
  nspecies = as.numeric(Param$simulation$nspecies)
  
  savenpv = list() #create an empty list
  
  for (s in 1: simu) {
    npvalpha =  list() #create an empty list
  
    for (i in 1 : nspecies) {
      # For each species, call the fishing mortality file
      index = paste0("sp",i-1)
      selectivity.file = paste0(input.folder,"/Eco-OSMOSE-NS/fishing_rate/mortality.fishing.rate.bydt.bysize.file.",index,".csv")
      F.rate = read.csv(selectivity.file, sep=",", header = F)
      fishing.folder = "C:/R/osmose-ns/Eco-OSMOSE-NS/fishing_rate"
  
      # Random draw of a parameter that will modify the fishing mortality F
      alpha <- runif(1, Fmin, Fmax)
      npvalpha[[i]] = alpha
      cF=ncol(F.rate)
      rF=nrow(F.rate)
      NewF <- F.rate[2:rF,2:cF] * alpha    
    
      write.csv(NewF, file=paste0("Eco-OSMOSE-NS/Fnpv/Fnpv_sp", i, "_simu", s, ".csv"), row.names = F)
      Param$mortality$fishing$rate$bydt$bysize$file[i] = paste0("Eco-OSMOSE-NS/Fnpv/Fnpv_sp", i,"_simu", s, ".csv")
   
    }

    write.csv(npvalpha, file= paste0("Eco-OSMOSE-NS/Fnpv/Npv_alpha_simu_", s, ".csv"))
  
    simulation = paste0("npv_simu_", s)
  
    configFile = file.path(input.folder, input.file) # path to main configuration file
    outputDir  = file.path(input.folder, "Eco-OSMOSE-NS/output", simulation)
  
    run_osmose(input = configFile, output = outputDir,
              osmose = jar.file, version="4.3.3",  parameters="-Poutput.netcdf.format=netcdf3")

  ### Get Net Present Value ---------------------------------------------------
  
    replicat = Param$simulation$nsimulation
    replicat = as.numeric(replicat)-1
    npv =  list() #create an empty list
  
    for (r in 0: replicat) {
      mylist <- list() #create an empty list
      nreplicat <- paste0("_simu", r)
    
      for (i in 1: nspecies) {
        spname = Param$species$name[i]
        file = paste0(outputDir,"/Econ/ns_HarvestingCosts", spname, nreplicat, ".csv")
        HarvestCost = read.csv(file)
        mylist[[i]] <- HarvestCost[,2] #put all vectors in the list
      }
    
      HCosts = as.data.frame(mylist)
      totCosts = rowSums(HCosts)
    
      timestart = as.numeric(Param$output$start$year)
      timeend = as.numeric(Param$simulation$time$nyear)
      time = timeend - timestart
    
      file = paste0(outputDir,"/Econ/ns_Utility", nreplicat, ".csv")
      utility = read.csv(file)
      socialoptimum <- list() #create an empty list
      eta = as.numeric(Param$elasticity$demand$fish)
      gamma = as.numeric(Param$weight$fish$consumption)
    
      discountfactor = exp(log(0.05) / time)
      
      for (t in 1: time) {
        utilityt = utility[t, 2]
        totalutility = gamma * (eta / (eta - 1)) * utilityt ^ ((eta -1)/ eta)
        socialoptimum[[t]] = (discountfactor ^ t) *(totalutility - totCosts[t])
      }
      socialoptimum = as.data.frame(socialoptimum)
      netpresentvalue = rowSums(socialoptimum)
    
      npv[[r+1]] = netpresentvalue
    }
    
    npv = as.numeric(npv)
    mean.npv = mean(npv)
    
    savenpv[[s]] = mean.npv
  }
  write.csv(savenpv, file=paste0("Eco-OSMOSE-NS/npv", "_after", s, "replicat.csv"))
} # end of Fnpv
