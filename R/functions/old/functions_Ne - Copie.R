# equations to estimate Ne --------------------
# MB
# 11/06/25
# based on Waples et al 2025

library(osmose)
library(stringr)
library(reshape2)
library(ggplot2)
library(tidyr)


estimate_Ne <- function(data_Ne, species, nb_replicates){
  
  # rescale abundance, to have a constant pop size at birth of 1000 # OK
  for (sp in species){
    for (rep in 1:nb_rep){
      abund <- data_Ne[["abundanceByAge"]][[sp]]
      if(dim(abund)[3]>1) abund <- abund[[rep]]
      years <- nrow(abund)
      ncohort <- ncol(abund)
      abund_scaled <- matrix(c(-99), nrow = years, ncol=ncohort)
      for(coho in 1:ncohort){
        abund_scaled[,coho] <- 1000*abund[,coho,rep]/abund[1,coho,rep] # abundance 2 ???????***********
      }
      colnames(abund_scaled) <- colnames(abund)
      rownames(abund_scaled) <- rownames(abund)
      data_Ne[["abund_scaled"]][[sp]][[paste0("rep_", rep)]] <- (abund_scaled)
    }
  }
  
  # Dead individuals # OK
  for (sp in species){
    for (rep in 1:nb_rep){
      abund_scaled <- data_Ne[["abund_scaled"]][[sp]][[rep]]
      years <- nrow(abund_scaled)
      ncohort <- ncol(abund_scaled)
      dead_mat <- matrix(c(-99), nrow = years, ncol=ncohort)
      for(coho in 1:ncohort){
        dead_mat[,coho] <- c(abund_scaled[1:c(years-1),coho]-abund_scaled[2:years,coho], 0)
        dead_mat[years,coho] <- 1000-sum(dead_mat[,coho])
      }
      colnames(dead_mat) <- colnames(abund_scaled)
      rownames(dead_mat) <- rownames(abund_scaled)
      data_Ne[["Dead"]][[sp]][[paste0("rep_", rep)]] <- dead_mat
    }
  }
  
  # produce the survival indicator # OK
  for (sp in species){
    for (rep in 1:nb_rep){
      dead <- data_Ne[["Dead"]][[sp]][[paste0("rep_", rep)]]
      years <- nrow(dead)
      ncohort <- ncol(dead)
      dead_sum <- apply(dead, 2, cumsum)
      alive_sum <- 1000-dead_sum
      alive_frac <- matrix(-99, nrow = years, ncol=ncohort)
      alive_frac <- rbind(abund_scaled[c(2:years),]/abund_scaled[c(1:c(years-1)),], rep(0, ncohort))
      alive_frac[which(alive_frac=="NaN", arr.ind = T)] <- 0
      alive_cumsum <- matrix(NA, nrow = years, ncol=ncohort)
      alive_cumsum[1,] <- rep(1, ncohort)
      for (i in 2:years){
        alive_cumsum[i,] <- alive_cumsum[c(i-1),]*alive_frac[c(i-1),]
      }
          #wrong
         # alive_cumsum <- apply(alive_sum, 2, cumsum)
         # alive_cumsum_frac <- alive_sum*1/1000
      colnames(alive_frac) <- colnames(dead)
      colnames(alive_cumsum) <- colnames(dead)
      rownames(alive_frac) <- rownames(dead)
      rownames(alive_cumsum) <- rownames(dead)
      data_Ne[["survival_rate"]][[sp]][[paste0("rep_", rep)]] <- alive_frac
      data_Ne[["survivorship_curve"]][[sp]][[paste0("rep_", rep)]] <- alive_cumsum
    }
  }
  
  
  # rescale birth rate and and offspring number (constant pop size of 1000) #OK (issue of the not respected colSums(birthrate_scaled*surviv_curve, na.rm = T) always = 1 )
  
  for (sp in species){
    for (rep in 1:nb_rep){
      
      # sex_r <- sex_ratio[[stringr::str_replace(species_i[which(species == sp)], pattern = "species.name.", replacement = "species.sexratio.")]]
      eggNb <- data_Ne[["eggNumber"]][[sp]] #/sex_r # calculations are different if sex ratio isn't 0.5 
      surviv_curve <- data_Ne[["survivorship_curve"]][[sp]][[paste0("rep_", rep)]]
      if(dim(eggNb)[3]>1) eggNb <- eggNb[[rep]]
      abund <- data_Ne[["abundanceByAge"]][[sp]]
      if(dim(abund)[3]>1) abund <- abund[[rep]]
      abund_scaled <- data_Ne[["abund_scaled"]][[sp]][[rep]]
      years <- nrow(abund)
      ncohort <- ncol(abund)
      
      birthrate <- as.data.frame(eggNb)/abund #bx: mean number of newborns produced by 1 individual at age x
      birthrate[which(birthrate=="NaN", arr.ind = T)] <- 0
      birthrate_scaled <- birthrate/(apply(birthrate*surviv_curve, 2, cumsum)) # rescaling method 1 
      birthrate_scaled <- 1000*birthrate/(apply(birthrate*abund_scaled, 2, cumsum)) # rescaling method 2 
      # birthrate_scaled[which(birthrate_scaled=="NaN", arr.ind = T)] <- 0
      birthrate_scaled[1,] <- 0
      eggNb_scaled <- birthrate_scaled*abund_scaled #Bx: number of births by cohort of age x (bx*Nx)
      eggNb_scaled[which(eggNb_scaled=="NaN", arr.ind = T)] <- 0
      colnames(birthrate_scaled) <- colnames(eggNb)
      rownames(birthrate_scaled) <- rownames(eggNb)
      colnames(eggNb_scaled) <- colnames(eggNb)
      rownames(eggNb_scaled) <- rownames(eggNb)
      data_Ne[["eggNb_sc"]][[sp]][[paste0("rep_", rep)]] <- eggNb_scaled
      data_Ne[["b_rate_sc"]][[sp]][[paste0("rep_", rep)]] <- birthrate_scaled
      data_Ne[["b_rate"]][[sp]][[paste0("rep_", rep)]] <- birthrate
    }
  }
  eggNb_scaled_calculations <- eggNb_scaled
  eggNb_scaled_calculations[which(eggNb_scaled_calculations == 0, arr.ind = T)] <- NA
  should_be_popsize <- colMeans(eggNb_scaled_calculations, na.rm = T) # should be 1000
  birthrate_scaled_calculation <- birthrate_scaled
  birthrate_scaled_calculation[which(birthrate_scaled_calculation == 0, arr.ind = T)] <- NA
  should_be_one <-colMeans(birthrate_scaled_calculation*surviv_curve, na.rm = T)
  
  boolean <- all(round(should_be_one) == 1, na.rm = T)
  cat("\nThe sum of survival curve * birth rate == 1, is ", boolean, "\n")
  # check : colSums(birthrate_scaled*surviv_curve, na.rm = T) always = 1 except last cohort 
  
  # estimate generation length
  
  for (sp in species){
    for (rep in 1:nb_rep){
      eggNb_scaled <- data_Ne[["eggNb_sc"]][[sp]]
      if(length(eggNb_scaled)>0) eggNb_scaled <- eggNb_scaled[[rep]]
      cnames <- colnames(eggNb_scaled) 
      rnames <- rownames(eggNb_scaled) 
      eggNb_scaled <- as.data.frame(eggNb_scaled)
      pop_size <- 1000 # if it was scaled
      years <- nrow(eggNb_scaled)
      ncohort <- ncol(eggNb_scaled)
      Longevity <- matrix(c(-99), nrow = years, ncol=ncohort)
      for(coho in 1:ncohort){
        Longevity[,coho] <- sum(eggNb_scaled[,coho]*c(1:years)/pop_size) # start from age 0 or 1 ??
      }
      colnames(Longevity) <- cnames
      rownames(Longevity) <- rnames
      data_Ne[["longevity"]][[sp]][[paste0("rep_", rep)]] <- Longevity
    }
  }
  
  # estimate mean and variance lifetime reproductive success
  
  for (sp in species){
    for (rep in 1:nb_rep){
      birthrate_scaled <- as.matrix(data_Ne[["b_rate_sc"]][[sp]][[rep]])
      years <- nrow(birthrate_scaled)
      ncohort <- ncol(birthrate_scaled)
      data_Ne[["kx_bar"]][[sp]][[paste0("rep_", rep)]] <- apply(ifelse(is.nan(birthrate_scaled), 0, birthrate_scaled), 2, cumsum)
      data_Ne[["Vx"]][[sp]][[paste0("rep_", rep)]] <- data_Ne[["kx_bar"]][[sp]][[rep]]
      colnames(data_Ne[["kx_bar"]][[sp]][[paste0("rep_", rep)]]) <- colnames(birthrate_scaled)
      rownames(data_Ne[["kx_bar"]][[sp]][[paste0("rep_", rep)]]) <- rownames(birthrate_scaled)
      colnames(data_Ne[["Vx"]][[sp]][[paste0("rep_", rep)]]) <- colnames(birthrate_scaled)
      rownames(data_Ne[["Vx"]][[sp]][[paste0("rep_", rep)]]) <- rownames(birthrate_scaled)
    }
  }
  
  
  # Mean lifetime repro success overall
  for (sp in species){
    for (rep in 1:nb_rep){
      dead <- data_Ne[["Dead"]][[sp]][[rep]]
      kx_bar <- data_Ne[["kx_bar"]][[sp]][[rep]]
      pop_size <- 1000
      years <- nrow(dead)
      ncohort <- ncol(dead)
      # data_Ne[["k_point"]][[sp]][[paste0("rep_", rep)]] <- apply(dead*kx_bar, 2, cumsum)/pop_size
      data_Ne[["k_point"]][[sp]][[paste0("rep_", rep)]] <- t(matrix(rep(colSums(dead*kx_bar)/pop_size, years), ncol=years, nrow = ncohort))
      colnames(data_Ne[["k_point"]][[sp]][[paste0("rep_", rep)]]) <- colnames(dead)
      rownames(data_Ne[["k_point"]][[sp]][[paste0("rep_", rep)]]) <- rownames(dead)
    }
  }
  
  # Mean square deviation lifetime repro success overall 
  for (sp in species){
    for (rep in 1:nb_rep){
      dead <- data_Ne[["Dead"]][[sp]][[rep]]
      var_x <- data_Ne[["Vx"]][[sp]][[rep]]
      k_point <- data_Ne[["k_point"]][[sp]][[rep]]
      kx_bar <- data_Ne[["kx_bar"]][[sp]][[rep]]
      pop_size <- 1000
      years <- nrow(dead)
      ncohort <- ncol(dead)
      SSDix <- var_x*dead
      Delta <- kx_bar - k_point
      SSDgx <- dead*Delta^2
      SSD <- SSDix+SSDgx
      Vk_point <- t(matrix(rep(colSums(SSD), years), ncohort, years)) /pop_size
      Vk_point[which(SSD == 0, arr.ind = T)] <- 0
      data_Ne[["Vk_point"]][[sp]][[paste0("rep_", rep)]] <- Vk_point #range is wrong
      
      colnames(data_Ne[["Vk_point"]][[sp]][[paste0("rep_", rep)]]) <- colnames(dead)
      rownames(data_Ne[["Vk_point"]][[sp]][[paste0("rep_", rep)]]) <- rownames(dead)
    }
  }
  
  
  # Effective population size # OK
  for (sp in species){
    for (rep in 1:nb_rep){
      longevity <- data_Ne[["longevity"]][[sp]][[rep]]
      k_point <- data_Ne[["k_point"]][[sp]][[rep]]
      Vk_point <- data_Ne[["Vk_point"]][[sp]][[rep]]
      pop_size <- 1000
      years <- nrow(dead)
      ncohort <- ncol(dead)
      data_Ne[["Ne"]][[sp]][[paste0("rep_", rep)]] <-
        (4*pop_size*longevity)/(Vk_point+2)      
      data_Ne[["Ne"]][[sp]][[paste0("rep_", rep)]][which(Vk_point == 0, arr.ind = T)] <- 0
      colnames(data_Ne[["Ne"]][[sp]][[paste0("rep_", rep)]]) <- colnames(longevity)
      rownames(data_Ne[["Ne"]][[sp]][[paste0("rep_", rep)]]) <- rownames(longevity)
    }
  }
  
  return(data_Ne)
}
