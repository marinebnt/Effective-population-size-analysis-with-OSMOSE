# equations to estimate temporal Ne --------------------
# based on Marty et al 2014 and Waples 1989
# MB
# 11/06/25


estimate_temporalNe <- function(outputDir_Simu, outputDir_Ne, species, conf, seeding_y, OUTPUT_data){
  temporalNe <- list()
  Ne = c()
  for (sp in species){
    cat("\n", sp)
    sp_id <- which(sp %in% sp_names)

    # output dimension and data frame
    dimensions_outputs = dim(OUTPUT_data$biomass)
    first_y_Ne_tempo <- max(seeding_y, c(years-dimensions_outputs[1]+1))
    Ne = data.frame(matrix(ncol=nb_rep, nrow=length(first_y_Ne_tempo:(years-20))))
    rownames(Ne) = first_y_Ne_tempo:(years-20)
    colnames(Ne) = paste0("rep_", 1:nb_rep)
    gene_length <- data.frame(matrix(ncol=nb_rep, nrow=years))
    rownames(gene_length) = 1:years
    colnames(gene_length) = paste0("rep_", 1:nb_rep)

    for (rep in 1:nb_rep){

      # allele freq data 
      data <- ncdf4::nc_open(paste0(outputDir_Simu, "/Genetic/ns_alleleFrequencyDistribByAge-", sp, "Simu", c(rep-1), ".nc"))
      data_nc_raw <- ncdf4::ncvar_get(data, "allele_occurrence_frequency")
      data_nc_raw <- data_nc_raw[ , ,-1, ] # removing age 0 
      data_nc <- apply(data_nc_raw, c(1,2,4), function(x) mean(x, na.rm = TRUE)) # mean over ages
      if (length(dim(data_nc))==4) {trait1_data <- data_nc[, , 1, ] # trait 1 is the neutral trait 
        }else {trait1_data <- data_nc} 
      ncdf4::nc_close(data)
      trait1_data <- trait1_data/100 # we want something <1

      # life history parameters
    # abundance by age
      abund_tot <- OUTPUT_data[["abundanceByAge"]][[sp]][,,rep]
    # survivorship curve
      dead_mat <- (abund_tot[1:c(nrow(abund_tot)-1),-1]-abund_tot[2:nrow(abund_tot),-1]) #removing age 0
      dead_mat <- rbind(dead_mat, 0) # to have the same dimension as abund_tot
      rownames(dead_mat) <- rownames(abund_tot)
      a <- as.numeric(colnames(dead_mat)) #ages
    # sx
      alive_frac <- rbind(abund_tot[c(2:nrow(abund_tot)),-1]/abund_tot[c(1:c(nrow(abund_tot)-1)),-1], 0)
      alive_frac[which(alive_frac=="NaN", arr.ind = T)] <- 0
      rownames(alive_frac) <- rownames(abund_tot)
    # lx
      alive_cumsum <- abund_tot/abund_tot[,2]
      alive_cumsum <- alive_cumsum[, -1] #removing age 0
      alive_cumsum[which(alive_cumsum=="NaN", arr.ind = T)] <- 0
      lambdaa = t(alive_cumsum)
    #  raw birth rate
      eggNb_full = read.csv(paste0(outputDir_Simu, "/AgeIndicators/ns_fecundityDistribByAge_Simu", c(rep-1), ".csv"), skip=c(1))
      eggNb <- dcast(eggNb_full, Time ~ Age, value.var = sp) # does not include the replicates
      eggNb <- eggNb[ , -1] # removing time column
      eggNb <- eggNb[ , -1] # removing age 0
      birthrate <- eggNb/abund_tot[, which(colnames(abund_tot) %in% colnames(eggNb))] #bx: mean number of newborns produced by 1 individual at age x
      birthrate[which(birthrate=="NaN", arr.ind = T)] <- 0
      birthrate[which(birthrate=="Inf", arr.ind = T)] <- 0
    # rescaling birth rate 
      rescaling_f <- t(matrix(apply(birthrate*alive_cumsum, 2, sum), ncol= nrow(birthrate), nrow=ncol(birthrate)))
      birthrate_scaled <- birthrate/rescaling_f # rescaling method 1 
      birthrate_scaled[which(birthrate_scaled=="NaN", arr.ind = T)] <- 0
      birthrate_scaled[which(birthrate_scaled=="Inf", arr.ind = T)] <- 0
      ba = t(birthrate_scaled)
      generation_length <- c() #rep(NA, years)
      for (y in c(years-dimensions_outputs[1]+1):years){
        ii <- grep(pattern=y, colnames(lambdaa))[1]
        generation_length[y] <- sum(a*(lambdaa[,ii]*ba[,ii]/sum(lambdaa[,ii]*ba[,ii])), na.rm = T)
      }
      names(generation_length) = 1:years
      gene_length[, paste0("rep_", rep)] <- generation_length

      # save info in the output list
      temporalNe[[sp]][["alive_cumsum"]][[paste0("rep_", rep)]] = lambdaa
      temporalNe[[sp]][["birthrate_scaled"]][[paste0("rep_", rep)]] = ba

      # loop over years to estimate Nb
      t1=0
      Ne_tempo <- c()
      allelet1 <- c()
      allelet2 <- c()
      for (year_i in (first_y_Ne_tempo:(years-20))){
        t1 = t1+1
        t2 = t1+20

        abund <- rowSums(abund_tot)
        abund_1 <- abund[t1]
        abund_2 <- abund[t2]

        # F : the variance in allele frequency change between two samples 
        # According to Nei and Tajima 1981
        #number of  allelic states
        L = as.numeric(get_var(conf, paste0("neutral.trait.nval.sp", sp_id)))
        #number of loci
        K = as.numeric(get_var(conf, paste0("neutral.trait.nlocus.sp", sp_id)))
        #frequency of allelic state l at locus k
        Sum_on_L = 0
        for (li in 1:L){
          Sum_on_K = 0
          for (ki in 1:K){
            numerator   <- (trait1_data[li,ki,(t1)]-trait1_data[li,ki,(t2)])^2
            denominator <- (trait1_data[li,ki,(t1)]+trait1_data[li,ki,(t2)])/2 - (trait1_data[li,ki,t1]*trait1_data[li,ki,t2])
            Sum_on_K = Sum_on_K + numerator/denominator
          }
          Sum_on_L = Sum_on_L + (1/K)*Sum_on_K
        }
        F_Ne = (1/L)*Sum_on_L
        # Ne : the effective population size
        # longevity <- data_Ne[["longevity"]][[sp]][[paste0("rep_", rep)]]
        # generation_length <- mean(longevity[1,2:20])
        
        numerator = t2-t1+1
        denominator = 2*mean(generation_length[which(names(generation_length) == year_i)], na.rm = T)*(F_Ne+1/(2*abund_1)-1/(2*abund_2))
        Ne_tempo[t1] = numerator/denominator
        allelet1[t1] = trait1_data[li,ki,(t1)]
        allelet2[t1] = trait1_data[li,ki,(t2)]
      }
      Ne[, paste0("rep_", rep)] = Ne_tempo
    }#rep
    temporalNe[[sp]][["allelet2"]] = allelet2
    temporalNe[[sp]][["allelet1"]] = allelet1
    temporalNe[[sp]][["Ne_temporal"]] = Ne
    temporalNe[[sp]][["generation_length"]] = gene_length[which(rownames(gene_length) %in% rownames(Ne)), ]
  }#sp
  return(temporalNe)
}
# end of file ----------------------------------------------------------------
