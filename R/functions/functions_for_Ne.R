# data conversion, to further input in the process to estimate Ne --------------------
# marine beneat
# 11/06/25


estimate_for_Ne <- function(outputDir_Simu, outputDir_Ne, dim_vec, species, nb_rep, OUTPUT_data, rep){
  VAR_for_NE <- list() # this is the final list with 1 list per trait dimension, and each time an element of the list per species. And a species contains an array with 3Dimensions : per time step and per cohort and per replicate
  for (TRAIT in dim_vec){
    cat("\n", TRAIT)
    tempo_list_sp <- list()
    tempo_list_sp_survivor <- list()

    for (sp in species){

      expected_array_dim <- dim(as.data.frame(OUTPUT_data[["abundanceByAge"]][[sp]][,,1]))
      expected_cohort_number <- expected_array_dim[1]
      expected_max_age       <- expected_array_dim[2]
      vec_ages               <- round(as.numeric(colnames(as.data.frame(OUTPUT_data[["abundanceByAge"]][[sp]][,,1]))))
      vec_year               <- round(as.numeric(rownames(as.data.frame(OUTPUT_data[["abundanceByAge"]][[sp]][,,1]))))
      tempo_array <- array(data=c(NA), dim=c(expected_max_age, expected_cohort_number, nb_rep))
      tempo_array <- provideDimnames(tempo_array, sep = "_", base = list('time','cohort','rep'))

      for (rep in (0:c(nb_rep-1))){

        # dimensions of the dataset for that species and that replicate
        if(any(grepl("egg", TRAIT) || grepl("school", TRAIT) || grepl("mature", TRAIT))){
          if (grepl("egg", TRAIT)) path_i = "/AgeIndicators/ns_fecundityDistribByAge_Simu"
          if (grepl("mature", TRAIT)) path_i = "/Indicators/ns_matureAbundanceDistribByAge_Simu"
          if (grepl("school", TRAIT)) path_i = "/ns_nschoolDistribByAge_Simu"
          d_i <- read.csv(paste0(outputDir_Simu, path_i, rep, ".csv"), skip=c(1))
          data <- dcast(d_i, Time ~ Age, value.var = sp) # does not include the replicates
          rownames(data) <- data[,1]
          data <- data[,-1]
        }
        if(grepl("mortality", TRAIT)){
          tempo_array_survivor <- tempo_array
          mortality <- OUTPUT_data[[TRAIT]][[sp]]
          mortality <- Reduce("+", lapply(mortality, function(x) x[,,rep+1]))
          # data      <- exp(-mortality) # survival rate
          data      <- mortality
          if (ncol(data) < 30) data <- cbind(data, matrix(0, nrow=nrow(data), ncol=(31 - ncol(data)))) # to have at least 30 ages
        }
        if (TRAIT %in% c("abundanceMature", "sizeMature")){
          age_Mature = OUTPUT_data[["ageMature"]]
          age_Mature = lapply(seq(dim(age_Mature)[3]), function(x) as.data.frame(age_Mature[ , , x]))
          age_mat_sp_rep = age_Mature[[rep+1]][,sp]
          if (TRAIT == "abundanceMature") vec_y = OUTPUT_data[["abundanceByAge"]][[sp]][,,rep+1]
          if (TRAIT == "sizeMature")      vec_y = OUTPUT_data[["meanSizeByAge"]][[sp]][,,rep+1]
          for (y in 1:length(age_mat_sp_rep)){
              if (TRAIT == "abundanceMature") vec_y[,which(as.numeric(colnames(OUTPUT_data$abundanceByAge[[sp]]))< ceiling(age_mat_sp_rep[y]))] <- 0
              if (TRAIT == "sizeMature")      vec_y[,which(as.numeric(colnames(OUTPUT_data$meanSizeByAge[[sp]]))< ceiling(age_mat_sp_rep[y]))] <- 0
          }
          data = cbind(vec_y, matrix(rep(NA, nrow(vec_y)*(31-ncol(vec_y))), nrow=nrow(vec_y), ncol=31-ncol(vec_y))) # to have at least 30 ages
        }
        if (TRAIT == ("fecundityByAge")) {
          eggNumber <- read.csv(paste0(outputDir_Simu, "/AgeIndicators/ns_fecundityDistribByAge_Simu", rep, ".csv"), skip=c(1))
          eggNumber <- eggNumber[, which(colnames(eggNumber)%in% c("Age", sp, "Time"))]
          abundance <- as.data.frame(OUTPUT_data[["abundanceByAge"]][[sp]][,,rep+1])
          abundance$Time <- as.numeric(rownames(abundance))
          abundance <- pivot_longer(abundance, cols =  -Time, names_to = "Age", values_to = "abundance")
          abundance$Age <- as.numeric(abundance$Age)
          full_join(eggNumber, abundance, by = c("Time", "Age")) %>%
            mutate(fecundity = get(sp)/abundance) -> fecundity_data
          fecundity_data$fecundity[-(grep('^-?[0-9.]+$', fecundity_data$fecundity))] <- 0
          data_tot <- fecundity_data[, which(colnames(fecundity_data) %in% c("Age", "Time", "fecundity"))]
          colnames(data_tot) <- c("Time", "Age", sp)
          data_tot <- data_tot[, which(colnames(data_tot) %in% c("Time", "Age", sp))]
          data <- pivot_wider(data_tot, names_from = Age, values_from = colnames(data_tot)[-c(1,2)])
          rownames(data) <- as.matrix(data[,1])
          colnames(data) <- as.numeric(colnames(data))
          data <- data[,-1]
        }

        if(!((TRAIT == "mature") || grepl("mortality", TRAIT) || grepl("egg", TRAIT)  || grepl("school", TRAIT) || 
             grepl("sizeMature", TRAIT) || grepl(TRAIT, pattern="abundanceMature") || TRAIT == ("fecundityByAge"))) {
            data <- OUTPUT_data[[TRAIT]][[sp]]
          if(!(TRAIT %in% names(OUTPUT_data[["mortalityByAge"]][[sp]]))) {data <- data[,,rep+1]
          }else {
            data <- OUTPUT_data[["mortalityByAge"]][[sp]][[TRAIT]][,,rep+1]
          }
          if (ncol(data) < 30) data <- cbind(data, matrix(0, nrow=nrow(data), ncol=(31 - ncol(data)))) # to have at least 30 ages
        }

        # error message dataset missing
        if (any(dim(data)==c(0,0))) {
          cat("ERROR : Dimension data trait ", TRAIT, "are NULL")
          next
        }
        
        # identify the cohorts in the dataset
        tempo <- as.data.frame(data)
        ID_matrix <- row(tempo) - col(tempo)
        colnames(ID_matrix) <- vec_ages
        rownames(ID_matrix) <- vec_year 
        ID_cohort <- table(ID_matrix)
        
        # order data by cohort
        by_cohort <- list()
        for (i in 1:length(ID_cohort)){
          ID <- names(ID_cohort)[i]
          by_cohort[[ID]] <- tempo[which(ID_matrix == ID, arr.ind = T)]
        }
        ID_cohort_positive <- which(as.numeric(names(by_cohort))>=0)
        nb_cohort <- length(ID_cohort_positive)
        
        # species, trait data
        mat_tempo <- matrix(NA, nrow = expected_max_age, ncol=nb_cohort)
        dimensions_outputs = dim(OUTPUT_data$biomass)
        colnames(mat_tempo) <- paste0("cohort_", c(ny-dimensions_outputs[1]+1):ny)
        rownames(mat_tempo) <- vec_ages
        k=0
        # if(grepl("mortality", TRAIT)){mat_tempo_survivor <- mat_tempo}
        for (cohort in 0:c(nb_cohort-1)){
          k=k+1
          vec_tempo     <- rep(0, expected_max_age)
          arr_id <- which(ID_matrix == cohort, arr.ind = T)
          
          arr_id_cohort <- as.numeric(arr_id[,2])
          vec_tempo[arr_id_cohort] <- as.data.frame(data)[arr_id] #unlist(by_cohort[[as.character(cohort)]])
          mat_tempo[,k] <- vec_tempo
          # if(grepl("mortality", TRAIT)){ # estimate the survivorship curve from the survival rate
          #   conditional_survival <- rep(NA, length(vec_tempo) - 1)
          #   for (i in 1:(length(vec_tempo)-1)) {
          #     if (vec_tempo[i] > 0) {
          #       conditional_survival[i] <- vec_tempo[i + 1] / vec_tempo[i]
          #     } else {
          #       conditional_survival[i] <- 0  # or NA, depending on how you want to treat total extinction
          #     }
          #   }
          #   mat_tempo_survivor[,k] <- c(1,conditional_survival)
          # }
        }#Cohort
        tempo_array[,,(rep+1)] <- mat_tempo
        # if(grepl("mortality", TRAIT)){tempo_array_survivor[,,rep] <- mat_tempo_survivor}
      }#Replicate
      dimnames(tempo_array) <- list(rownames(mat_tempo), colnames(mat_tempo), paste0("rep_", 1:nb_rep))
      tempo_list_sp[[sp]] <- tempo_array
      # if(grepl("mortality", TRAIT)){
      #   dimnames(tempo_array_survivor) <- list(rownames(mat_tempo), colnames(mat_tempo), paste0("rep_", 1:nb_rep))
      #   tempo_list_sp_survivor[[sp]] <- tempo_array_survivor
      # }
    }#Species
    VAR_for_NE[[TRAIT]] <- tempo_list_sp
    # if(grepl("mortality", TRAIT)){
    #   VAR_for_NE[["Survivorship"]] <- tempo_list_sp_survivor
    # }
  }#Trait
  
  VAR_for_NE[["cohort_time_step_corresponding"]] <- ID_matrix
  #if(!dir.exists(outputDir_Ne)) dir.create(outputDir_Ne, recursive = T)
  #save(VAR_for_NE, file=paste0(outputDir_Ne, "/Data_for_Ne.RData"))
  return(VAR_for_NE)
}
