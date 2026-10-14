# equations to estimate Ne --------------------
# MB
# 11/06/25
# based on Waples et al 2025


estimate_Nb <- function(outputDir_Simu, outputDir_Ne, species, nb_rep, OUTPUT_data){
  
  list_output <- list()
  new_age_one <- list()
  dimensions_outputs = dim(OUTPUT_data$biomass)
  first_y_Ne_tempo <- max(seeding_y, c(years-dimensions_outputs[1]+1))
  
  for (sp in species){
    NT = data.frame(matrix(ncol=nb_rep, nrow=length(first_y_Ne_tempo:(years))))
    rownames(NT) = first_y_Ne_tempo:(years)
    colnames(NT) = paste0("rep_", 1:nb_rep)
    k_point <- NT
    V_k_point <- NT
    Nb <- NT
    
    for (rep in 1:nb_rep){
      cat("\n", sp)
      
      #abundance by age : sum of all biomasses at all ages
      abundanceByAge           <- t(data.frame((OUTPUT_data[["abundanceByAge"]][[sp]][,,rep])))
      rownames(abundanceByAge) <- colnames(OUTPUT_data[["abundanceByAge"]][[sp]][,,rep])
      N_age_1                  <- abundanceByAge[which(rownames(abundanceByAge)==1), ]
      abundanceByAge           <- abundanceByAge[which(rownames(abundanceByAge)>=1), ]
      NT[, paste0("rep_", rep)]<- colSums(abundanceByAge)

      
      # produce the survival indicator # OK
      # sx
      age_max = nrow(abundanceByAge)
      alive_frac <- rbind(abundanceByAge[c(2:age_max),]/abundanceByAge[c(1:c(age_max-1)),],0)
      rownames(alive_frac) <- rownames(abundanceByAge)
      alive_frac[which(alive_frac=="NaN", arr.ind = T)] <- 0
      # lx
      alive_cumsum <- alive_frac # just giving it the size of a matrix as it is expected
      i=0
      for (y in colnames(alive_cumsum)){
        i=i+1
        alive_cumsum[,i] <- abundanceByAge[,which(colnames(abundanceByAge)==y)]/abundanceByAge[1,which(colnames(abundanceByAge)==y)]
      }
      alive_cumsum[which(alive_cumsum=="NaN", arr.ind = T)] <- 0
      list_output[[sp]][["survival_rate"]][[paste0("rep_", rep)]] <- alive_frac
      list_output[[sp]][["surv_curve"]][[paste0("rep_", rep)]] <- alive_cumsum
      
     
      # k_point
      k_point[, paste0("rep_", rep)] = 2*N_age_1/NT[, paste0("rep_", rep)]
      
      #birth rate
      eggNb_raw           <- read.csv(paste0(outputDir_Simu, "/AgeIndicators/ns_fecundityDistribByAge_Simu", c(rep-1),".csv"), skip=c(1))
      # from the first age at reproduction
      eggNb_raw           <- eggNb_raw[as.numeric(eggNb_raw$Age)>=1, ]
      eggNb_raw           <- eggNb_raw[, colnames(eggNb_raw) %in% c("Age","Time",sp)]
      eggNb               <- eggNb_raw %>%
        select(-Age) %>%         # remove age
        group_by(Time) %>%       # group by time
        summarise(across(everything(), sum, na.rm = TRUE))
      # rownames(eggNb) <- eggNb$Time
      list_output[[sp]][["eggNb"]][[paste0("rep_", rep)]] <- eggNb
      birthrate_scaled = data.frame("time"=names(table(eggNb_raw$Age)))
      for (t in names(table(eggNb_raw$Time))){
        # scaling birth rate and egg number
        abundanceByAge_t <- abundanceByAge[,which(colnames(abundanceByAge)==as.numeric(t))]
        eggNb_t          <- eggNb_raw[which(eggNb_raw$Time==as.numeric(t)),]
        birthrate        <- eggNb_t/abundanceByAge_t
        if (any(birthrate=="NaN", arr.ind = T)) {birthrate[which(birthrate=="NaN", arr.ind = T)] <- 0}
        alive_cumsum_t   <- alive_cumsum[,which(colnames(alive_cumsum)==as.numeric(t))]
        rescaling_f      <- 0.5*sum(as.matrix(birthrate*alive_cumsum_t), na.rm = T)
        birthrate_scaled[,t] <- birthrate[, sp]/rescaling_f # rescaling method 1 
      }
      birthrate_scaled <- birthrate_scaled[, -which(colnames(birthrate_scaled) == c("time"))]
      list_output[[sp]][["birthrate_scaled"]][[paste0("rep_", rep)]] = birthrate_scaled
      
      
      #test : we should have mean birth rate = k_point 
      #cat("    TESTS : FALSE = Failed, TRUE = Successed ===>   ")
      #cat(any((N_age_1/NT[, paste0("rep_", rep)] * colSums(birthrate_scaled*alive_cumsum, na.rm = T)) == k_point[, paste0("rep_", rep)], na.rm = T))
      
      
      Vk_point_mat = abundanceByAge*(birthrate_scaled+(birthrate_scaled-k_point[, paste0("rep_", rep)])^2)
      V_k_point[, paste0("rep_", rep)] = sapply(1:ncol(Vk_point_mat), function(x){sum(Vk_point_mat[,x]/NT[x, paste0("rep_", rep)], na.rm = T)})

      Nb[, paste0("rep_", rep)] <- (k_point[, paste0("rep_", rep)]*NT[, paste0("rep_", rep)]-2)/(k_point[, paste0("rep_", rep)]-1+V_k_point[, paste0("rep_", rep)]/k_point[, paste0("rep_", rep)])
    }
    list_output[[sp]][["NT"]] <- NT
    list_output[[sp]][["k_point"]] <- k_point
    list_output[[sp]][["V_k_point"]]<- V_k_point
    list_output[[sp]][["Nb"]] <- Nb
  }
  
  data_Nb <- list_output
  return(data_Nb)
}
