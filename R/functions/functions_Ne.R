# equations to estimate Ne --------------------
# MB
# 11/06/25
# based on Waples et al 2025

# WARNING : all variables estimated based on haploid life-histories. 
#        The cases haploid and diploid are described in Waples et al., 2011. 
#        We have adapted the equation of the haploid populations to the diploid case, 
#        while employing the parameters estimated for the haploid case.
#        Full description is give :  DATA/Haploid_Ne_for_Diploid.pdf



estimate_Ne <- function(outputDir_Ne, data_Ne_shape, species, nb_rep, input_dir, OUTPUT_data, removing_senescence=F){

  cutoff_age <- get_par(conf, "output.cutoff.age")
  data_Ne_shape2 <- list()
  out_plots <- paste0(outputDir_Ne, "/plots")

  new_age_one <- list()

  # quick fix : Mortality total is the sum of natural mortality and fishing mortaity. 
  # natural mortality is the sum of Mpred, Madd, Mstar
  # Mnat in the actual output is in fact Madd. Let's correct that and creat Mtot
  names(data_Ne_shape)[which(names(data_Ne_shape) == "Mnat")] <- "Madd"
  for (sp in species){
    data_Ne_shape[["Mnat"]][[sp]] <- data_Ne_shape[["Madd"]][[sp]] + data_Ne_shape[["Mpred"]][[sp]] + data_Ne_shape[["Mstar"]][[sp]]
    data_Ne_shape[["Mtot"]][[sp]] <- data_Ne_shape[["Mnat"]][[sp]] + data_Ne_shape[["F"]][[sp]]
  }

  for (sp in species){
    
    # # we are removing the larval stage of our species
    # sp_id <- c(0:c(nb_sp-1))[which(species == sp)]
    # cutoff_age[grep(sp_id, names(cutoff_age))][1]
    # 
    # for (rep in 1:nb_rep){
    #   row_fec_not_null <- which(data_Ne_shape$eggNumber[[sp]][,,paste0("rep_", rep)]>0, arr.ind = T)
    #   row_fec_not_null.min <- as.data.frame(row_fec_not_null) %>%
    #     group_by(col) %>%
    #     summarize(minrow = min(row))
    #   new_age_one[[sp]][[rep]] = as.numeric(rownames(data_Ne_shape$eggNumber[[sp]][,,paste0("rep_", rep)])[floor(mean(row_fec_not_null.min$minrow, na.rm = T))])
    # }
    # if (any(is.na(unlist(new_age_one[[sp]])))) {new_age_one[[sp]][[which(is.na(unlist(new_age_one[[sp]])))]] <- mean(unlist(new_age_one[[sp]]), na.rm=T)}
    # 
    
    # reloading data_Ne but reversing the order of the outputs to facilitate manipulation
    for (rep in 1:nb_rep){
     for (nms in names(data_Ne_shape)){ # all the dimensions variables are considered from age 1, thus Ne is estimated by variables from age 1, except for Numerator_N_alpha
        if (nms == "cohort_time_step_corresponding"){next}
        if (is.matrix(data_Ne_shape[[nms]])){
          data_Ne_shape2[[sp]][[nms]][[paste0("rep_", rep)]] <- data_Ne_shape[[nms]][which(as.numeric(rownames(data_Ne_shape[[nms]][[sp]]))>=1), ]
        }
        else {
          data_Ne_shape2[[sp]][[nms]][[paste0("rep_", rep)]] <- data_Ne_shape[[nms]][[sp]][,,paste0("rep_", rep)][which(as.numeric(rownames(data_Ne_shape[[nms]][[sp]]))>=1), ]
        }
      }
    }

    # remove abundance before maturity
    for (rep in 1:nb_rep){
      abund <- data_Ne_shape2[[sp]][["abundanceByAge"]][[paste0("rep_", rep)]]
      egg   <- data_Ne_shape2[[sp]][["eggNumber"]][[paste0("rep_", rep)]]
      if (any(abund>0 & egg!=0)) {        
        id_not_remove  <- which(abund>0 & egg!=0, arr.ind = T)
        mean_age_abund <- mean(id_not_remove[,"row"], na.rm=T)
      } else mean_age_abund = 5
      id_remove      <- which((abund>0 & egg==0), arr.ind = T)
      if(length(id_remove)>1){
        id_remove  <- id_remove[which(id_remove[,"row"]<mean_age_abund), "row"]
        if(length(id_remove)>1){ 
          row_remove <- unique(id_remove)
          # cat(row_remove, "\n")
          data_Ne_shape2[[sp]] <- lapply(data_Ne_shape2[[sp]], FUN = function(list) {
            if(any(class(list) == "list")) list[[rep]] <- list[[rep]][-row_remove,]
            return(list)
          })
        }
      }
    }
    # get_abund_age <- function(age_get, rep_i){
    #   if (age_get == "maturation"){age_get = as.numeric(new_age_one[[sp]][[rep_i]])}
    #   matrix <- data_Ne_shape[["abundanceByAge"]][[sp]][which(as.numeric(rownames(data_Ne_shape[["abundanceByAge"]][[sp]][,,rep_i]))==age_get), ,paste0("rep_", rep_i)]
    #   return(matrix)
    # }
    # 
    # # if new age 1 is maturation age, or age 1 = total population size at age 1 or mature
    # numerator_n_alpha <- sapply(1:nb_rep, FUN= function(x) get_abund_age("maturation", x))
    # colnames(numerator_n_alpha) <- paste0("rep_", 1:nb_rep)
    # data_Ne_shape2[[sp]][["Numerator_N_alpha"]]<- numerator_n_alpha
    # numerator_n_1    <- sapply(1:nb_rep, FUN= function(x) get_abund_age(1, x))
    # colnames(numerator_n_1) <- paste0("rep_", 1:nb_rep)
    data_Ne_shape2[[sp]][["Numerator_N_alpha"]]<- as.data.frame(t(do.call(rbind, lapply(data_Ne_shape2[[sp]][["abundanceByAge"]], FUN = function(matrix_rep) {
        matrix_rep[1,]
      }
    ))))
    
  }
  


  
  # estimate Dead before removing senescence (removing data at last age), otherwise huge mortality at last step
  for (sp in species){
    for (rep in 1:nb_rep){
      # Dead individuals # OK
      abundanceByAge <- data_Ne_shape2[[sp]][["abundanceByAge"]][[paste0("rep_", rep)]]
      years_i  <- nrow(abundanceByAge)
      ncohort  <- ncol(abundanceByAge)
      dead_mat <- matrix(c(-99), nrow = years_i, ncol=ncohort)
      pop_size <- data_Ne_shape2[[sp]][["Numerator_N_alpha"]][,paste0("rep_", rep)]
      for(coho in 1:ncohort){
        dead_mat[,coho] <- c(abundanceByAge[1:c(years_i-1),coho]-abundanceByAge[2:years_i,coho],0)
        dead_mat[which(dead_mat[,coho]==0)[1],coho] <- pop_size[coho]-sum(dead_mat[,coho]) 
      }
      for (Mort in c("Mpred", 'Mstar', "Mnat", 'F', "Mfor", "Mtot", "Madd")){
        data_Ne_shape2[[sp]][[paste0("Dead_", Mort)]][[paste0("rep_", rep)]] <- 
          data_Ne_shape2[[sp]][[paste0(Mort)]][[paste0("rep_", rep)]]/data_Ne_shape2[[sp]][["Mtot"]][[paste0("rep_", rep)]]*dead_mat
      }
      colnames(dead_mat) <- colnames(abundanceByAge)
      rownames(dead_mat) <- rownames(abundanceByAge)
      data_Ne_shape2[[sp]][["Dead"]][[paste0("rep_", rep)]] <- dead_mat
    }
  }
  # removing senescence
  if (removing_senescence) {
    data_Ne_shape2 <- lapply(data_Ne_shape2, function(splist){
      lapply(splist, function(traitlist){
        if(is.matrix(traitlist)) {
          return(traitlist)}
        lapply(traitlist, function(replist){
          matrix <- as.matrix(data.frame(row=(as.numeric(apply(replist, 2, function(x) min(which(x==0))))-1), col=1:ncol(replist)))
          replist[matrix] <- 0
          return(replist)
        })
      })
    })
  }  
  # add data new first age
  data_Ne_shape2[["new_first_age"]] <- new_age_one
  # data_Ne_shape2[["table_mature"]]  <- table_mature
  
 # output dimension and data frame
  dimensions_outputs = dim(OUTPUT_data$biomass)
  first_y_Ne_tempo <- max(seeding_y, c(years-dimensions_outputs[1]+1))
  generation_length = data.frame(matrix(ncol=nb_rep, 
              nrow=length(first_y_Ne_tempo:(years))))
  rownames(generation_length) = first_y_Ne_tempo:(years)
  colnames(generation_length) = paste0("rep_", 1:nb_rep)
  Vk_point_df = generation_length
  Ne_N_alpha_df = generation_length
  Ne_N_1_df = generation_length
  k_point_df = generation_length

  for (sp in species){
    cat("\n", "sp", sp, "\n  Sum survival curve * birth rate [50]= ")
    for (rep in 1:nb_rep){
      
      # major input
      lifespan <- get_par(conf, par = paste0("lifespan.sp", unlist(stringr::str_split(names(get_par(conf, paste0("species.name")))[which(get_par(conf, paste0("species.name"))==sp)], pattern="sp"))[3]))[1]
      abundanceByAge <- data_Ne_shape2[[sp]][["abundanceByAge"]][[paste0("rep_", rep)]]
      dead           <- data_Ne_shape2[[sp]][["Dead"]][[paste0("rep_", rep)]]
      dead_Mpred     <- data_Ne_shape2[[sp]][["Dead_Mpred"]][[paste0("rep_", rep)]]
      dead_Mstar     <- data_Ne_shape2[[sp]][["Dead_Mstar"]][[paste0("rep_", rep)]]
      dead_Mnat      <- data_Ne_shape2[[sp]][["Dead_Mnat"]][[paste0("rep_", rep)]]
      dead_F         <- data_Ne_shape2[[sp]][["Dead_F"]][[paste0("rep_", rep)]]
      dead_Mfor      <- data_Ne_shape2[[sp]][["Dead_Mfor"]][[paste0("rep_", rep)]]
      dead_Madd      <- data_Ne_shape2[[sp]][["Dead_Madd"]][[paste0("rep_", rep)]]
      eggNb          <- data_Ne_shape2[[sp]][["eggNumber"]][[paste0("rep_", rep)]] #/sex_r # calculations are different if sex ratio isn't 0.5 
      pop_size       <- data_Ne_shape2[[sp]][["Numerator_N_alpha"]][,paste0("rep_", rep)]
      
      # scale of the outputs
      years_i <- nrow(dead)
      ncohort <- ncol(dead)
      
      # birth rate
      birthrate <- as.data.frame(eggNb)/abundanceByAge #bx: mean number of newborns produced by 1 individual at age x
      birthrate[which(birthrate=="NaN", arr.ind = T)] <- 0
      birthrate[which(birthrate=="Inf", arr.ind = T)] <- 0
      
      # produce the survival indicator # OK
      # sx
      alive_frac <- rbind(abundanceByAge[c(2:years_i),]/abundanceByAge[c(1:c(years_i-1)),], 0)
      rownames(alive_frac) <- rownames(dead)
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
      data_Ne_shape2[[sp]][["survival_rate"]][[paste0("rep_", rep)]] <- alive_frac
      data_Ne_shape2[[sp]][["survivorship_curve"]][[paste0("rep_", rep)]] <- alive_cumsum

      # rescale birth rate and and offspring number (constant pop size of N1) #OK (issue of the not respected colSums(birthrate_scaled*surviv_curve, na.rm = T) always = 1 )
      surviv_curve = alive_cumsum
      # sex_r <- sex_ratio[[stringr::str_replace(species_i[which(species == sp)], pattern = "species.name.", replacement = "species.sexratio.")]]
    
      #scaling
      rescaling_f <- t(matrix(apply(birthrate*surviv_curve, 2, sum, na.omit=T), ncol= nrow(birthrate), nrow=ncol(birthrate)))
      birthrate_scaled <- birthrate/rescaling_f # rescaling method 1 
      
      colnames(birthrate_scaled) <- colnames(dead)
      rownames(birthrate_scaled) <- rownames(dead)      
      colnames(birthrate) <- colnames(dead)
      rownames(birthrate) <- rownames(dead)
      data_Ne_shape2[[sp]][["b_rate_sc"]][[paste0("rep_", rep)]] <- birthrate_scaled
      data_Ne_shape2[[sp]][["b_rate"]][[paste0("rep_", rep)]] <- birthrate
      
      # save also the egg number associated with the scaled birthrate
      eggNb <- abundanceByAge*birthrate_scaled 
      eggNb[which(eggNb=="NaN", arr.ind = T)] <- 0
      colnames(eggNb) <- colnames(dead)
      rownames(eggNb) <- rownames(dead)
      data_Ne_shape2[[sp]][["eggNumber_sc"]][[paste0("rep_", rep)]]  <- eggNb

      # verifying that the egg production = as many adult N1 of this cohort after rescaling
      # eggNb_scaled_calculations <- eggNb_scaled
      # eggNb_scaled_calculations[which(eggNb_scaled_calculations == 0, arr.ind = T)] <- NA
      # should_be_popsize <- colSums(eggNb_scaled_calculations, na.rm = T) # should be N1
      # should_be_one <- !any(c(pop_size[1,] - should_be_popsize) >1)
      
      # birth rate * survivorship curve = 1
      test <- apply(birthrate_scaled*surviv_curve, 2, cumsum) 
      # test <- apply(test, 2, function(x) max(x, na.rm = T)) 
      cat(max(test[,50]))
      # check : colSums(birthrate_scaled*surviv_curve, na.rm = T) always = 1 except last cohort 

      # estimate generation length
      cnames <- colnames(eggNb) 
      rnames <- as.numeric(rownames(eggNb))
      eggNb  <- as.data.frame(eggNb)
      ages     <- as.numeric(rownames(eggNb))
      ncohort  <- ncol(eggNb)
      Longevity <- c()
      for(coho in 1:ncohort){
        Longevity[coho] <- sum(eggNb[,coho]*ages)/pop_size[coho] # start from age 0 or 1 ??
      }
      Longevity <- t(as.data.frame(Longevity))
      rownames(Longevity) <- NULL
      colnames(Longevity) <- cnames
      generation_length[, paste0("rep_", rep)] <- t(Longevity)

      # estimate mean and variance lifetime reproductive success
      birthrate_scaled <- as.matrix(birthrate_scaled)
      ncohort <- ncol(birthrate_scaled)
      data_Ne_shape2[[sp]][["kx_bar"]][[paste0("rep_", rep)]] <- apply(ifelse(is.nan(birthrate_scaled), 0, birthrate_scaled), 2, cumsum)
      data_Ne_shape2[[sp]][["Vx"]][[paste0("rep_", rep)]]     <- data_Ne_shape2[[sp]][["kx_bar"]][[paste0("rep_", rep)]]
      colnames(data_Ne_shape2[[sp]][["kx_bar"]][[paste0("rep_", rep)]]) <- colnames(birthrate_scaled)
      rownames(data_Ne_shape2[[sp]][["kx_bar"]][[paste0("rep_", rep)]]) <- rownames(birthrate_scaled)
      colnames(data_Ne_shape2[[sp]][["Vx"]][[paste0("rep_", rep)]])     <- colnames(birthrate_scaled)
      rownames(data_Ne_shape2[[sp]][["Vx"]][[paste0("rep_", rep)]])     <- rownames(birthrate_scaled)

      # Mean lifetime repro success overall
      kx_bar   <- data_Ne_shape2[[sp]][["kx_bar"]][[paste0("rep_", rep)]]
      k_point <- t(as.matrix(colSums(dead*kx_bar)/pop_size))
      if (any(is.nan(k_point))) k_point[which(is.nan(k_point), arr.ind = T)] <- 0
      k_point_df[, paste0("rep_", rep)] <- t(k_point)
      if(is.null( data_Ne_shape2[[sp]][["k_point"]])) { 
        data_Ne_shape2[[sp]][["k_point"]] <- data.frame(matrix(ncol=nb_rep, nrow=c(dimensions_outputs[1])))
        colnames(data_Ne_shape2[[sp]][["k_point"]]) <- paste0("rep_", 1:nb_rep)
        rownames(data_Ne_shape2[[sp]][["k_point"]]) <- colnames(dead)
      }
      data_Ne_shape2[[sp]][["k_point"]][, rep] <- k_point_df[, paste0("rep_", rep)]

      # Mean square deviation lifetime repro success overall 
      var_x    <- data_Ne_shape2[[sp]][["Vx"]][[paste0("rep_", rep)]]
      SSDix    <- var_x*dead
      Delta    <- kx_bar - t(matrix(rep(k_point, nrow=nrow(kx_bar)), ncol=nrow(kx_bar), nrow=ncol(kx_bar)))
      SSDgx    <- dead*Delta^2
      SSDgx_F      <- dead_F*Delta^2
      SSDgx_Mpred  <- dead_Mpred*Delta^2
      SSDgx_Mfor   <- dead_Mfor*Delta^2
      SSDgx_Mstar  <- dead_Mstar*Delta^2
      SSDgx_Madd   <- dead_Madd*Delta^2
      SSDgx_Mnat   <- dead_Mnat*Delta^2
      SSD      <- SSDix+SSDgx
      Vk_point <- colSums(SSD, na.rm = T)/pop_size
      if(any(is.nan(Vk_point))) Vk_point[which(is.nan(Vk_point), arr.ind = T)] <- 0
      Vk_point_df[, paste0("rep_", rep)] <- Vk_point
      
      data_Ne_shape2[[sp]][["Delta"]][[paste0("rep_", rep)]] <- Delta
      colnames(data_Ne_shape2[[sp]][["Delta"]][[paste0("rep_", rep)]]) <- colnames(birthrate_scaled)
      rownames(data_Ne_shape2[[sp]][["Delta"]][[paste0("rep_", rep)]]) <- rownames(birthrate_scaled)           
      data_Ne_shape2[[sp]][["Deltasquare"]][[paste0("rep_", rep)]] <- Delta^2
      colnames(data_Ne_shape2[[sp]][["Deltasquare"]][[paste0("rep_", rep)]]) <- colnames(birthrate_scaled)
      rownames(data_Ne_shape2[[sp]][["Deltasquare"]][[paste0("rep_", rep)]]) <- rownames(birthrate_scaled)     
      data_Ne_shape2[[sp]][["SSDix"]][[paste0("rep_", rep)]] <- SSDix
      colnames(data_Ne_shape2[[sp]][["SSDix"]][[paste0("rep_", rep)]]) <- colnames(birthrate_scaled)
      rownames(data_Ne_shape2[[sp]][["SSDix"]][[paste0("rep_", rep)]]) <- rownames(birthrate_scaled)      
      data_Ne_shape2[[sp]][["SSDgx"]][[paste0("rep_", rep)]] <- SSDgx
      colnames(data_Ne_shape2[[sp]][["SSDgx"]][[paste0("rep_", rep)]]) <- colnames(birthrate_scaled)
      rownames(data_Ne_shape2[[sp]][["SSDgx"]][[paste0("rep_", rep)]]) <- rownames(birthrate_scaled)
      data_Ne_shape2[[sp]][["SSDgx_F"]][[paste0("rep_", rep)]] <- SSDgx_F
      colnames(data_Ne_shape2[[sp]][["SSDgx_F"]][[paste0("rep_", rep)]]) <- colnames(birthrate_scaled)
      rownames(data_Ne_shape2[[sp]][["SSDgx_F"]][[paste0("rep_", rep)]]) <- rownames(birthrate_scaled)       
      data_Ne_shape2[[sp]][["SSDgx_Mpred"]][[paste0("rep_", rep)]] <- SSDgx_Mpred
      colnames(data_Ne_shape2[[sp]][["SSDgx_Mpred"]][[paste0("rep_", rep)]]) <- colnames(birthrate_scaled)
      rownames(data_Ne_shape2[[sp]][["SSDgx_Mpred"]][[paste0("rep_", rep)]]) <- rownames(birthrate_scaled)       
      data_Ne_shape2[[sp]][["SSDgx_Mfor"]][[paste0("rep_", rep)]] <- SSDgx_Mfor
      colnames(data_Ne_shape2[[sp]][["SSDgx_Mfor"]][[paste0("rep_", rep)]]) <- colnames(birthrate_scaled)
      rownames(data_Ne_shape2[[sp]][["SSDgx_Mfor"]][[paste0("rep_", rep)]]) <- rownames(birthrate_scaled)       
      data_Ne_shape2[[sp]][["SSDgx_Mstar"]][[paste0("rep_", rep)]] <- SSDgx_Mstar
      colnames(data_Ne_shape2[[sp]][["SSDgx_Mstar"]][[paste0("rep_", rep)]]) <- colnames(birthrate_scaled)
      rownames(data_Ne_shape2[[sp]][["SSDgx_Mstar"]][[paste0("rep_", rep)]]) <- rownames(birthrate_scaled)       
      data_Ne_shape2[[sp]][["SSDgx_Madd"]][[paste0("rep_", rep)]] <- SSDgx_Madd
      colnames(data_Ne_shape2[[sp]][["SSDgx_Madd"]][[paste0("rep_", rep)]]) <- colnames(birthrate_scaled)
      rownames(data_Ne_shape2[[sp]][["SSDgx_Madd"]][[paste0("rep_", rep)]]) <- rownames(birthrate_scaled)       
      data_Ne_shape2[[sp]][["SSDgx_Mnat"]][[paste0("rep_", rep)]] <- SSDgx_Mnat
      colnames(data_Ne_shape2[[sp]][["SSDgx_Mnat"]][[paste0("rep_", rep)]]) <- colnames(birthrate_scaled)
      rownames(data_Ne_shape2[[sp]][["SSDgx_Mnat"]][[paste0("rep_", rep)]]) <- rownames(birthrate_scaled)       


      # Effective population size # OK
      longevity <- Longevity
      years_i   <- nrow(k_point)
      ncohort   <- ncol(k_point)

      # Ne <- (4*pop_size*longevity)/(Vk_point+2) # diploid, haploid
      Ne_N_alpha        <- (4*pop_size*longevity)/(2*Vk_point + 2 + (2*colSums(SSDgx))/pop_size) # diploid and 2 sexes (diploid)
      # here : all variables estimated based on haploid life-histories. 
      #        The cases haploid and diploid are described in Waples et al., 2011. 
      #        We have adapted the equation of the haploid populations to the diploid case, 
      #        while employing the parameters estimated for the haploid case.
      #        Full description is give :  DATA/Haploid_Ne_for_Diploid.pdf
      Ne_N_alpha[which(Vk_point == 0, arr.ind = T)] <- 0  
      Ne_N_alpha_df[, paste0("rep_", rep)] <- t(Ne_N_alpha)
    }# rep

    data_Ne_shape2[[sp]][["DeadSum"]] <- (sapply(1:nb_rep, function(x) colSums(data_Ne_shape2[[sp]][["Dead"]][[x]])))
    data_Ne_shape2[[sp]][["generation_length"]] <- generation_length
    data_Ne_shape2[[sp]][["k_point"]] <- k_point_df
    data_Ne_shape2[[sp]][["Ne_N_alpha"]] <- Ne_N_alpha_df
    data_Ne_shape2[[sp]][["Vk_point"]] <- Vk_point_df

    colnames(data_Ne_shape2[[sp]][["DeadSum"]]) <- colnames(Vk_point_df)
    for (dims in names(data_Ne_shape2[[sp]])){
      if (is.list(data_Ne_shape2[[sp]][[dims]]) && length(data_Ne_shape2[[sp]][[dims]])==1){
        data_Ne_shape2[[sp]][[dims]] <- reshape_data(data_Ne_shape2[[sp]][[dims]], as.variable = T, varName = dims)
      }
    }
  }# sp
  data_Ne <- data_Ne_shape2
  return(data_Ne)
}

