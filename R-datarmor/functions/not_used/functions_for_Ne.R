# data conversion, to further input in the process to estimate Ne --------------------
# marine beneat
# 11/06/25

library(osmose)
library(stringr)
library(reshape2)


estimate_for_Ne <- function(outputDir_Simu, outputDir_Ne, dim_vec, species, nb_rep, OUTPUT_data){
  VAR_for_NE <- list() # this is the final list with 1 list per trait dimension, and each time an element of the list per species. And a species contains an array with 3Dimensions : per time step and per cohort and per replicate
  for (TRAIT in dim_vec){
    cat("\n", TRAIT)
    tempo_list_sp <- list()
    tempo_list_sp_survivor <- list()
    for (sp in species){
      expected_array_dim <- dim(as.data.frame(OUTPUT_data[[dim_vec[1]]][sp][[1]]))
      expected_cohort_number <- expected_array_dim[1]
      tempo_array <- array(data=c(-99), dim=c(length(time_vec), expected_cohort_number, nb_rep))
      tempo_array <- provideDimnames(tempo_array, sep = "_", base = list('time','cohort','rep'))
      if(grepl("mortality", TRAIT)){tempo_array_survivor <- tempo_array}
      for (rep in 1:nb_rep){
        # dimensions of the dataset for that species and that replicate
        if(grepl("egg", TRAIT)){
          fecundity <- read.csv(paste0(outputDir_Simu, "/AgeIndicators/ns_fecundityDistribByAge_Simu0.csv"), skip=c(1))
          data <- dcast(fecundity, Time ~ Age, value.var = sp) # does not include the replicates
          rownames(data) <- data[,1]
          data <- data[,-1]
        }
        if(grepl("school", TRAIT)){
          school <- read.csv(paste0(outputDir_Simu, "/ns_nschoolDistribByAge_Simu0.csv"), skip=c(1))
          data <- dcast(school, Time ~ Age, value.var = sp) # does not include the replicates
          rownames(data) <- data[,1]
          data <- data[,-1]
        }
        if(grepl("mortality", TRAIT)){
          mortality <- OUTPUT_data[[TRAIT]][sp][[rep]]
          mortality <- as.data.frame(Reduce('+', mortality))
          data      <- exp(-mortality) # survival
        }
        if(!(grepl("mortality", TRAIT) || grepl("egg", TRAIT)  || grepl("school", TRAIT))) data <- OUTPUT_data[[TRAIT]][sp][[rep]]
        
        # error message dataset missing
        if (any(dim(data)==c(0,0))) {
          cat("ERROR : Dimension data trait ", TRAIT, "are NULL")
          next
        }
        
        # identify the cohorts in the dataset
        tempo <- as.data.frame(data)
        ID_matrix <- row(tempo) - col(tempo)
        rownames(ID_matrix) <- 0:c(years-1)
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
        mat_tempo <- matrix(-99, nrow = length(time_vec), ncol=nb_cohort)
        colnames(mat_tempo) <- paste0("cohort_", 0:c(nb_cohort-1))
        rownames(mat_tempo) <- rownames(ID_matrix)
        k=0
        if(grepl("mortality", TRAIT)){mat_tempo_survivor <- mat_tempo}
        for (cohort in 0:c(nb_cohort-1)){
          k=k+1
          vec_tempo     <- rep(0, length(time_vec))
          arr_id <- which(ID_matrix == cohort, arr.ind = T)
          
          arr_id_cohort <- as.numeric(arr_id[,2])
          vec_tempo[arr_id_cohort] <- as.data.frame(data)[arr_id] #unlist(by_cohort[[as.character(cohort)]])
          mat_tempo[,k] <- vec_tempo
          if(grepl("mortality", TRAIT)){ # estimate the survivorship curve from the survival rate
            conditional_survival <- rep(NA, length(vec_tempo) - 1)
            for (i in 1:(length(vec_tempo)-1)) {
              if (vec_tempo[i] > 0) {
                conditional_survival[i] <- vec_tempo[i + 1] / vec_tempo[i]
              } else {
                conditional_survival[i] <- 0  # or NA, depending on how you want to treat total extinction
              }
            }
            mat_tempo_survivor[,k] <- c(1,conditional_survival)
          }
        }#Cohort
        tempo_array[,,rep] <- mat_tempo
        if(grepl("mortality", TRAIT)){tempo_array_survivor[,,rep] <- mat_tempo_survivor}
      }#Replicate
      dimnames(tempo_array) <- list(rownames(mat_tempo), colnames(mat_tempo), paste0("rep_", 1:nb_rep))
      tempo_list_sp[[sp]] <- tempo_array
      if(grepl("mortality", TRAIT)){
        dimnames(tempo_array_survivor) <- list(rownames(mat_tempo), colnames(mat_tempo), paste0("rep_", 1:nb_rep))
        tempo_list_sp_survivor[[sp]] <- tempo_array_survivor
      }
    }#Species
    VAR_for_NE[[TRAIT]] <- tempo_list_sp
    if(grepl("mortality", TRAIT)){VAR_for_NE[["Survivorship"]] <- tempo_list_sp}
  }#Trait
  
  VAR_for_NE[["cohort_time_step_corresponding"]] <- ID_matrix
  dir.create(outputDir_Ne, recursive = T)
  save(VAR_for_NE, file=paste0(outputDir_Ne, "/Data_for_Ne.RData"))
  return(VAR_for_NE)
}
