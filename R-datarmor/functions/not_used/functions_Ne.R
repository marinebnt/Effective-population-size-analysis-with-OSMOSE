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


estimate_Ne <- function(outputDir_Ne, data_Ne_shape, species, nb_rep, input_dir){
  
  out_plots <- paste0(outputDir_Ne, "/plots")
  dir.create(out_plots, recursive = T)
  new_age_one <- list()
  # if (!file.exists(paste0(out_plots, "/abundance_mature.RData"))){
  # table_mature <- abund_mature_by_sp(input_dir = input_dir, outputDir_Ne=out_plots)
  # } else (table_mature <- get(load(paste0(out_plots, "/abundance_mature.RData"))))
  for (sp in species){
    new_age_one[[sp]] <- floor(min(ns1$ageMature[,sp,rep], na.rm = T))
    for (nms in names(data_Ne_shape)){
      cat(sp, "changing dim data_Ne_shape\n")
      if (nms == "cohort_time_step_corresponding"){next}
      if (is.matrix(data_Ne_shape[[nms]])){
        cat("dim [[nms]]", nms, "\n")
        data_Ne_shape[[nms]] <- data_Ne_shape[[nms]][which(rownames(data_Ne_shape[[nms]][[sp]])>=new_age_one[[sp]]), ]
      } 
      else {
        cat("dim [[nms]][[sp]]", nms, "\n")
        data_Ne_shape[[nms]][[sp]] <- data_Ne_shape[[nms]][[sp]][,,1][which(rownames(data_Ne_shape[[nms]][[sp]])>=new_age_one[[sp]]), ]
      } 
    }
  }
  
    
  data_Ne_shape[["new_first_age"]] <- new_age_one
  # data_Ne_shape[["table_mature"]]  <- table_mature
  
  
  # Dead individuals # OK
  for (sp in species){
    for (rep in 1:nb_rep){
      abundanceByAge <- data_Ne_shape[["abundanceByAge"]][[sp]]
      years <- nrow(abundanceByAge)
      ncohort <- ncol(abundanceByAge)
      dead_mat <- matrix(c(-99), nrow = years, ncol=ncohort)
      pop_size <- t(matrix(rep(abundanceByAge[1,], nrow(abundanceByAge)), ncol=nrow(abundanceByAge), nrow=ncol(abundanceByAge)))
      for(coho in 1:ncohort){
        dead_mat[,coho] <- c(abundanceByAge[1:c(years-1),coho]-abundanceByAge[2:years,coho],0)
        dead_mat[which(dead_mat[,coho]==0)[1],coho] <- pop_size[1,coho]-sum(dead_mat[,coho]) 
      }
      colnames(dead_mat) <- colnames(abundanceByAge)
      rownames(dead_mat) <- rownames(abundanceByAge)
      data_Ne_shape[["Dead"]][[sp]][[paste0("rep_", rep)]] <- dead_mat
    }
  }
  
  # produce the survival indicator # OK
  for (sp in species){
    for (rep in 1:nb_rep){
      abundanceByAge <- data_Ne_shape[["abundanceByAge"]][[sp]]
      dead <- data_Ne_shape[["Dead"]][[sp]][[paste0("rep_", rep)]]
      years <- nrow(dead)
      ncohort <- ncol(dead)
      
      # sx
      alive_frac <- rbind(abundanceByAge[c(2:years),]/abundanceByAge[c(1:c(years-1)),], rep(0, ncohort))
      alive_frac[which(alive_frac=="NaN", arr.ind = T)] <- 0
      
      # lx
      alive_cumsum <- dead # just giving it the size of a matrix as it is expected
      for (co in 1:ncohort){
        alive_cumsum[,co] <- abundanceByAge[,co]/abundanceByAge[1,co]
      }
      alive_cumsum[which(alive_cumsum=="NaN", arr.ind = T)] <- 0
      
      colnames(alive_frac) <- colnames(dead)
      colnames(alive_cumsum) <- colnames(dead)
      rownames(alive_frac) <- rownames(dead)
      rownames(alive_cumsum) <- rownames(dead)
      data_Ne_shape[["survival_rate"]][[sp]][[paste0("rep_", rep)]] <- alive_frac
      data_Ne_shape[["survivorship_curve"]][[sp]][[paste0("rep_", rep)]] <- alive_cumsum
    }
  }
  
  
  # rescale birth rate and and offspring number (constant pop size of N1) #OK (issue of the not respected colSums(birthrate_scaled*surviv_curve, na.rm = T) always = 1 )
  
  for (sp in species){
    for (rep in 1:nb_rep){
      surviv_curve = data_Ne_shape[["survivorship_curve"]][[sp]][[paste0("rep_", rep)]]
      # sex_r <- sex_ratio[[stringr::str_replace(species_i[which(species == sp)], pattern = "species.name.", replacement = "species.sexratio.")]]
      eggNb <- data_Ne_shape[["eggNumber"]][[sp]] #/sex_r # calculations are different if sex ratio isn't 0.5 
      eggNb <- t(apply(data_Ne_shape[["eggNumber"]][[sp]], MARGIN = 1, function(x) {x*(data_Ne_shape[["abundanceByAge"]][[sp]][1,]/colSums(eggNb))}))
      abund <- data_Ne_shape[["abundanceByAge"]][[sp]]
      years <- nrow(abund)
      ncohort <- ncol(abund)
      pop_size <- t(matrix(rep(abund[1,], nrow(abund)), ncol=nrow(abund), nrow=ncol(abund))) # individuals of the population at first mature age
      
      birthrate <- as.data.frame(eggNb)/abund #bx: mean number of newborns produced by 1 individual at age x
      birthrate[which(birthrate=="NaN", arr.ind = T)] <- 0
      rescaling_f <- t(matrix(apply(birthrate*surviv_curve, 2, sum), ncol= nrow(birthrate), nrow=ncol(birthrate)))
      birthrate_scaled <- birthrate/rescaling_f # rescaling method 1 
      
      # birthrate_scaled <- N1*birthrate/(apply(birthrate*abundanceByAge, 2, cumsum)) # rescaling method 2 
      # birthrate_scaled[which(birthrate_scaled=="NaN", arr.ind = T)] <- 0
      birthrate_scaled[1,] <- 0
      # eggNb_scaled <- birthrate_scaled*abundanceByAge #Bx: number of births by cohort of age x (bx*Nx)
      # eggNb_scaled[which(eggNb_scaled=="NaN", arr.ind = T)] <- 0
      colnames(birthrate_scaled) <- colnames(eggNb)
      rownames(birthrate_scaled) <- rownames(eggNb)
      # colnames(eggNb_scaled) <- colnames(eggNb)
      # rownames(eggNb_scaled) <- rownames(eggNb)
      # data_Ne_shape[["eggNb_sc"]][[sp]][[paste0("rep_", rep)]] <- eggNb_scaled
      data_Ne_shape[["b_rate_sc"]][[sp]][[paste0("rep_", rep)]] <- birthrate_scaled
      data_Ne_shape[["b_rate"]][[sp]][[paste0("rep_", rep)]] <- birthrate
      
      # save also the egg number associated with the scaled birthrate
      data_Ne_shape[["eggNumber_sc"]][[sp]][[rep]]  <- abund*birthrate_scaled 
    }
  }
  # verifying that the egg production = as many adult N1 of this cohort after rescaling
  # eggNb_scaled_calculations <- eggNb_scaled
  # eggNb_scaled_calculations[which(eggNb_scaled_calculations == 0, arr.ind = T)] <- NA
  # should_be_popsize <- colSums(eggNb_scaled_calculations, na.rm = T) # should be N1
  # should_be_one <- !any(c(pop_size[1,] - should_be_popsize) >1)
  
  # birth rate * survivorship curve = 1
  test <- apply(birthrate_scaled*surviv_curve, 2, cumsum) 
  test <- apply(test, 2, function(x) max(x, na.rm = T)) 
    
  cat("\nThe sum of survival curve * birth rate is ", test, "\n")
  # check : colSums(birthrate_scaled*surviv_curve, na.rm = T) always = 1 except last cohort 
  
  # estimate generation length
  
  for (sp in species){
    for (rep in 1:nb_rep){
      eggNb <- data_Ne_shape[["eggNumber_sc"]][[sp]][[rep]]
      abund <- data_Ne_shape[["abundanceByAge"]][[sp]]
      cnames <- colnames(eggNb) 
      rnames <- as.numeric(rownames(eggNb))
      eggNb  <- as.data.frame(eggNb)
      pop_size <-  t(matrix(rep(abund[1,], nrow(abund)), ncol=nrow(abund), nrow=ncol(abund)))
      years    <- nrow(eggNb)
      ncohort  <- ncol(eggNb)
      Longevity <- matrix(c(-99), nrow = years, ncol=ncohort)
      for(coho in 1:ncohort){
        Longevity[,coho] <- sum(eggNb[,coho]*c(1:years))/pop_size[,coho] # start from age 0 or 1 ??
      }
      colnames(Longevity) <- cnames
      rownames(Longevity) <- rnames
      data_Ne_shape[["longevity"]][[sp]][[paste0("rep_", rep)]] <- Longevity
    }
  }
  
  # estimate mean and variance lifetime reproductive success
  
  for (sp in species){
    for (rep in 1:nb_rep){
      birthrate_scaled <- as.matrix(data_Ne_shape[["b_rate_sc"]][[sp]][[rep]])
      ncohort <- ncol(birthrate_scaled)
      data_Ne_shape[["kx_bar"]][[sp]][[paste0("rep_", rep)]] <- apply(ifelse(is.nan(birthrate_scaled), 0, birthrate_scaled), 2, cumsum)
      data_Ne_shape[["Vx"]][[sp]][[paste0("rep_", rep)]]     <- data_Ne_shape[["kx_bar"]][[sp]][[rep]]
      colnames(data_Ne_shape[["kx_bar"]][[sp]][[paste0("rep_", rep)]]) <- colnames(birthrate_scaled)
      rownames(data_Ne_shape[["kx_bar"]][[sp]][[paste0("rep_", rep)]]) <- rownames(birthrate_scaled)
      colnames(data_Ne_shape[["Vx"]][[sp]][[paste0("rep_", rep)]]) <- colnames(birthrate_scaled)
      rownames(data_Ne_shape[["Vx"]][[sp]][[paste0("rep_", rep)]]) <- rownames(birthrate_scaled)
    }
  }
  
  
  # Mean lifetime repro success overall
  for (sp in species){
    for (rep in 1:nb_rep){
      dead <- data_Ne_shape[["Dead"]][[sp]][[rep]]
      kx_bar <- data_Ne_shape[["kx_bar"]][[sp]][[rep]]
      pop_size <- data_Ne_shape[["abundanceByAge"]][[sp]][1,]
      years <- nrow(dead)
      ncohort <- ncol(dead)
      data_Ne_shape[["k_point"]][[sp]][[paste0("rep_", rep)]] <- t(matrix(rep(colSums(dead*kx_bar)/pop_size, years), ncol=years, nrow=ncohort))
      colnames(data_Ne_shape[["k_point"]][[sp]][[paste0("rep_", rep)]]) <- colnames(dead)
      rownames(data_Ne_shape[["k_point"]][[sp]][[paste0("rep_", rep)]]) <- rownames(dead)
    }
  }
  
  # Mean square deviation lifetime repro success overall 
  for (sp in species){
    for (rep in 1:nb_rep){
      dead <- data_Ne_shape[["Dead"]][[sp]][[rep]]
      var_x <- data_Ne_shape[["Vx"]][[sp]][[rep]]
      k_point <- data_Ne_shape[["k_point"]][[sp]][[rep]]
      kx_bar <- data_Ne_shape[["kx_bar"]][[sp]][[rep]]
      pop_size <-  t(matrix(rep(data_Ne_shape[["abundanceByAge"]][[sp]][1,], nrow(dead)), ncol=nrow(dead), nrow=ncol(dead)))
      years <- nrow(dead)
      ncohort <- ncol(dead)
      SSDix <- var_x*dead
      Delta <- kx_bar - k_point
      SSDgx <- dead*Delta^2
      SSD <- SSDix+SSDgx
      Vk_point <- t(matrix(rep(colSums(SSD), years), ncohort, years)) /pop_size
      Vk_point[which(SSD == 0, arr.ind = T)] <- 0
      data_Ne_shape[["Vk_point"]][[sp]][[paste0("rep_", rep)]] <- Vk_point #range is wrong
      colnames(data_Ne_shape[["Vk_point"]][[sp]][[paste0("rep_", rep)]]) <- colnames(dead)
      rownames(data_Ne_shape[["Vk_point"]][[sp]][[paste0("rep_", rep)]]) <- rownames(dead)
    }
  }
  
  
  # Effective population size # OK
  for (sp in species){
    for (rep in 1:nb_rep){
      longevity <- data_Ne_shape[["longevity"]][[sp]][[rep]]
      k_point <- data_Ne_shape[["k_point"]][[sp]][[rep]]
      Vk_point <- data_Ne_shape[["Vk_point"]][[sp]][[rep]]
      pop_size <-  t(matrix(rep(data_Ne_shape[["abundanceByAge"]][[sp]][1,], nrow(k_point)), ncol=nrow(k_point), nrow=ncol(k_point)))
      years <- nrow(longevity)
      ncohort <- ncol(longevity)
      Ne <- (4*pop_size*longevity)/(Vk_point+2) 
      Ne[which(Vk_point == 0, arr.ind = T)] <- 0
      data_Ne_shape[["Ne"]][[sp]][[paste0("rep_", rep)]] <- Ne
      colnames(data_Ne_shape[["Ne"]][[sp]][[paste0("rep_", rep)]]) <- colnames(longevity)
      rownames(data_Ne_shape[["Ne"]][[sp]][[paste0("rep_", rep)]]) <- rownames(longevity)
    }
  }
  data_Ne <- data_Ne_shape
  save(data_Ne, file=paste0(outputDir_Ne, "/Ne.RData"))
  return(data_Ne)
}
