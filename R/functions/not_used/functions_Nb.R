# equations to estimate Ne --------------------
# MB
# 11/06/25
# based on Waples et al 2025

library(osmose)
library(stringr)
library(reshape2)
library(ggplot2)
library(tidyr)
library(ncdf4)



estimate_Nb <- function(outputDir_Simu, outputDir_Ne, species, nb_rep, input_dir, OUTPUT_data, data_Ne_output){
  
  list_output <- list()
  out_plots <- paste0(outputDir_Ne, "/plots")
  if(!dir.exists(out_plots)) dir.create(out_plots, recursive = T)
  new_age_one <- list()
  
  for (sp in species){
    cat("\n", sp, " Longevity : ")
    new_age_one <- data_Ne_output$new_first_age[[sp]]
    if (any(is.na(unlist(new_age_one[[sp]])))) {new_age_one[[sp]][[which(is.na(unlist(new_age_one[[sp]])))]] <- mean(unlist(new_age_one[[sp]]), na.rm=T)}
    for (rep in 1:nb_rep){
      
     ## same variables as in Ne_theoretical estimation process, but by year
      #abundance by age : sum of all biomasses at all ages 
      abundanceByAge           <- t(data.frame((OUTPUT_data[["abundanceByAge"]][[sp]][,,rep])))
      rownames(abundanceByAge) <- colnames(OUTPUT_data[["abundanceByAge"]][[sp]][,,rep])
      N_age_1                  <- abundanceByAge[which(rownames(abundanceByAge)==1), ]
      abundanceByAge           <- abundanceByAge[which(rownames(abundanceByAge)>=1), ]
      if(rownames(data_Ne_output$abundanceByAge[[sp]][[paste0("rep_", rep)]])[1] != rownames(abundanceByAge)[1]){
        a_nrow    <- which(rownames(abundanceByAge) %in% rownames(data_Ne_output$abundanceByAge[[sp]][[paste0("rep_", rep)]]))
        abundanceByAge       <- abundanceByAge[a_nrow, ]
      }
      NT                       <- colSums(abundanceByAge)
      list_output[["NT"]][[sp]][[paste0("rep_", rep)]] <- NT
      
      #dead
      dead_mat <- matrix(c(-99), nrow = nrow(abundanceByAge), ncol=ncol(abundanceByAge))
      colnames(dead_mat) = colnames(abundanceByAge)
      rownames(dead_mat) <- rownames(abundanceByAge)
      pop_size <- t(matrix(rep(abundanceByAge[1,], nrow(abundanceByAge)), ncol=nrow(abundanceByAge), nrow=ncol(abundanceByAge)))
      first_cohort = as.numeric(colnames(abundanceByAge)[1])
      dimensions_outputs = dim(OUTPUT_data$biomass)
      for(y in (ny-dimensions_outputs[1]+1):ny){
        colID =which(colnames(abundanceByAge)==y)
        dead_mat[,colID] <- c(abundanceByAge[1:c(nrow(abundanceByAge)-1),colID]-
                               abundanceByAge[2:nrow(abundanceByAge),colID],0)
        dead_mat[which(dead_mat[,colID]==0)[1],colID] <- pop_size[1,colID]-sum(dead_mat[,colID]) 
      }
      list_output[["Dead"]][[sp]][[paste0("rep_", rep)]] <- dead_mat
      
      # produce the survival indicator # OK
      dead <- dead_mat
      years <- nrow(dead)
      ncohort <- ncol(dead)
      # sx
      alive_frac <- rbind(abundanceByAge[c(2:years),]/abundanceByAge[c(1:c(years-1)),], rep(0, ncohort))
      alive_frac[which(alive_frac=="NaN", arr.ind = T)] <- 0
      # lx
      alive_cumsum <- dead # just giving it the size of a matrix as it is expected
      i=0
      for (y in (ny-dimensions_outputs[1]+1):ny){
        i=i+1
        alive_cumsum[,i] <- abundanceByAge[,which(colnames(abundanceByAge)==y)]/abundanceByAge[1,which(colnames(abundanceByAge)==y)]
      }
      alive_cumsum[which(alive_cumsum=="NaN", arr.ind = T)] <- 0
      if(rownames(data_Ne_output$abundanceByAge[[sp]][[paste0("rep_", rep)]])[1] != rownames(alive_frac)[1]){
        a_nrow    <- which(rownames(abundanceByAge) %in% rownames(data_Ne_output$abundanceByAge[[sp]][[paste0("rep_", rep)]]))
        abundanceByAge       <- abundanceByAge[a_nrow, ]
      }
      list_output[["survival_rate"]][[sp]][[paste0("rep_", rep)]] <- alive_frac
      list_output[["surv_curve"]][[sp]][[paste0("rep_", rep)]] <- alive_cumsum
      
     ## new variables 
      #kx and k_point
      eggNb_raw           <- read.csv(paste0(outputDir_Simu, "/AgeIndicators/ns_fecundityDistribByAge_Simu", c(rep-1),".csv"), skip=c(1))
          # from the first age at reproduction
      eggNb_raw           <- eggNb_raw[as.numeric(eggNb_raw$Age)>=1, ]
      eggNb_raw           <- eggNb_raw[, colnames(eggNb_raw) %in% c("Age","Time",sp)]
      eggNb               <- eggNb_raw %>%
        select(-Age) %>%         # remove age
        group_by(Time) %>%       # group by time
        summarise(across(everything(), sum, na.rm = TRUE))
      # rownames(eggNb) <- eggNb$Time
      list_output[["eggNb"]][[sp]][[paste0("rep_", rep)]] <- eggNb
      k_point <- as.data.frame(eggNb)[,sp]/NT # eggNb from first reproduction age
      k_point[which(k_point=="NaN", arr.ind = T)] <- 0
      list_output[["k_point"]][[sp]][[paste0("rep_", rep)]] <- k_point
      
      #Vkx and Vk_point
      Vk = c()
      sum_test = 0
      for (t in names(table(eggNb_raw$Time))){
        
        # scaling birth rate and egg number
        abundanceByAge_t <- abundanceByAge[,which(colnames(abundanceByAge)==as.numeric(t))]
        eggNb_t          <- eggNb_raw[which(eggNb_raw$Time==as.numeric(t)),]
        eggNb_t          <- eggNb_t[which(eggNb_t$Age %in% rownames(data_Ne_output$abundanceByAge[[sp]][[paste0("rep_", rep)]])),sp]
        birthrate        <- eggNb_t/abundanceByAge_t
        if (any(birthrate=="NaN", arr.ind = T)) {birthrate[which(birthrate=="NaN", arr.ind = T)] <- 0}
        alive_cumsum_t   <- alive_cumsum[,which(colnames(alive_cumsum)==as.numeric(t))]
        rescaling_f      <- 0.5*sum(as.matrix(birthrate*alive_cumsum_t), na.rm = T)
        birthrate_scaled <- birthrate/rescaling_f # rescaling method 1 
        eggNb_sc         <- birthrate_scaled*abundanceByAge_t
        k_point_t        <- k_point[which(names(k_point)==as.numeric(t))]
        # longevity
        # Longevity = sum(as.numeric(names(eggNb_sc))*eggNb_sc, na.rm = T)/N_age_1[which(names(N_age_1)==as.numeric(t))]
        Longevity = sum(as.numeric(names(birthrate_scaled))*birthrate_scaled*alive_cumsum_t)/sum(birthrate_scaled*alive_cumsum_t)
        # kx and Vkx are the same = birth rate*2 because hypo = hermaphrodite diploid population
        Vkx = 2*birthrate_scaled
        
        # SSDIx SSDGx
        SSDIx = 0.5*birthrate_scaled*dead[,which(colnames(dead)==as.numeric(t))]
        SSDGx = 0.5*dead[,which(colnames(dead)==as.numeric(t))]*
          (birthrate_scaled-2*N_age_1[which(names(N_age_1)==as.numeric(t))]/
             NT[which(names(N_age_1)==as.numeric(t))])^2
        
        D   = dead_mat[,which(colnames(dead_mat) == t)]
        SSDxi = Vkx*D
        SSDxg = D*(Vkx - k_point[t])
        SSDx  = SSDIx + SSDGx
        Vk[t] = sum(SSDx)/(0.5*N_age_1[which(names(N_age_1)==as.numeric(t))]) # Vkx[1] : number of individuals at first age at reproduction
        
        cat("TESTS : FALSE = Failed, TRUE = Successed")
        #test : sum birthrate_scaled*alivecumsum = 1
        sum_test = sum(birthrate_scaled*alive_cumsum_t)
        cat(sum_test == 1)
        # test : it is expected that mean(k_bar) = k_point
        cat((weighted.mean(x = Vkx, w = abundanceByAge_t)/k_point_t)==1)
        
        
      }
      cat("   ", Longevity)
      list_output[["Vk_point"]][[sp]][[paste0("rep_", rep)]] <- Vk

      
      Nb <- (k_point*NT-2)/(k_point-1+Vk/k_point)
      list_output[["Nb"]][[sp]][[paste0("rep_", rep)]] <- Nb
    }
  }
  
  data_Nb <- list_output
  save(data_Nb, file=paste0(outputDir_Ne, "/Nb.RData"))
  return(data_Nb)
}
