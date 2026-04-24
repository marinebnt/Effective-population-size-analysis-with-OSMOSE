# equations to estimate temporal Ne --------------------
# based on Marty et al 2014 and Waples 1989
# marine beneat
# 11/06/25


library(osmose)
library(stringr)
library(reshape2)
library(ggplot2)
library(tidyr)


estimate_temporalNe <- function(outputDir_Simu, data_Ne, species, conf, seeding_y){
  temporalNe <- list()
  Ne = c()
  for (sp in species){
    cat("\n", sp)
    k=0
    Ne = c()
    cohort_time_corr <- data_Ne[["cohort_time_step_corresponding"]]
    co_to_include <- which(colSums(data_Ne[[sp]][["longevity"]][[paste0("rep_", rep)]] )>0)
    
    # allele freq data 
    data <- ncdf4::nc_open(paste0(outputDir_Simu, "/Genetic/ns_alleleFrequency-", sp, "Simu0.nc"))
    data_nc <- ncdf4::ncvar_get(data, "allele_occurrence_frequency")
    trait1_data <- data_nc[, , 1, ] # trait 1 is the neutral trait 
    ncdf4::nc_close(data)
    trait1_data <- trait1_data/100
    
    for (t1 in seeding_y:(years-20)){
      cat(t1)
      k=k+1
      t2 = t1+20
      
      lambdaa = data.frame(data_Ne$survivorship_curve[[sp]][[rep]])
      # lambdaa <- as.data.frame(rowMeans(lambdaa))
      lambdaa$age <- rownames(lambdaa)
      ba = data.frame(data_Ne$b_rate_sc[[sp]][[1]])
      # ba <- as.data.frame(rowMeans(ba))
      ba$age <- as.numeric(rownames(ba))
      a = as.numeric(lambdaa$age)
      generation_length <-c()
      for (co in 1:ncohort){
        generation_length[co] <- sum(a*lambdaa[,co]*ba[,co])/ sum(lambdaa[,co]*ba[,co])
      }
      
        # abundance
      abund <- ns1[["abundanceByAge"]][sp][[rep]]
      # cat(generation_time, "\n")
      abund <- rowSums(abund[,,1])
      abund_1 <- abund[t1]
      abund_2 <- abund[t2]
      
      # F : the variance in allele frequency change between two samples 
      # According to Nei and Tajima 1981
       #number of  allelic states
       L = as.numeric(get_var(conf, "neutral.trait.nval.sp0"))
       #number of loci
       K = as.numeric(get_var(conf, "neutral.trait.nlocus.sp0"))
       #frequency of allelic state l at locus k
       Sum_on_L = 0
       for (li in 1:L){
         Sum_on_K = 0
         for (ki in 1:K){
           numerator   <- (trait1_data[li,ki,t1]-trait1_data[li,ki,t2])^2
           denominator <- (trait1_data[li,ki,t1]+trait1_data[li,ki,t2])/2 -
             (trait1_data[li,ki,t1]*trait1_data[li,ki,t2])
           Sum_on_K = Sum_on_K + numerator/denominator
         }
         Sum_on_L = Sum_on_L + (1/K)*Sum_on_K
       }
       F_Ne = (1/L)*Sum_on_L
      # Ne : the effective population size
       # longevity <- data_Ne[["longevity"]][[sp]][[paste0("rep_", rep)]]
       # generation_length <- mean(longevity[1,2:20])
       
       numerator = t2-t1+1
       denominator = 2*mean(generation_length, na.rm = T)*(F_Ne+1/(2*abund_1)-1/(2*abund_2))
       Ne[k] = numerator/denominator
    }
    names(Ne) <- seeding_y:(years-20)
    temporalNe[[sp]] = Ne
  }
  save(temporalNe, file=paste0(outputDir_Ne, "/temporalNe.RData"))
  return(temporalNe)
}