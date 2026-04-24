##########
## Nc has multiple ways to be estimated 
###########


# Nc = sum all individuals from the first age at maturity 
Nc_from_maturity <- function(DATA_NE, species){
  Nc_from_maturity <- list()
  raw_ab_c         <- data_Ne$abundanceByAge
  for (sp in species){
    raw_ab_c_i <- as.data.frame(raw_ab_c[[sp]])
    Nc <- as.numeric(colSums(raw_ab_c_i))
    Nc_from_maturity[[sp]] <- Nc
  }
  return(Nc_from_maturity)
}


# Nc = sum of only mature individuals from the first age at maturity 
Nc_only_mature <- function(DATA_NE, species){
  tab_Nc_Nr_all <- DATA_NE$table_mature
  Nc_only_mature <- list()
  cohort_time_corr <- DATA_NE$cohort_time_step_corresponding
  for (sp in species){
    tab_Nc_Nr <- NULL
    tab_Nc_Nr <- tab_Nc_Nr_all[,,sp]
    Nc <- c()
    for (co in names(table(cohort_time_corr[which(cohort_time_corr>=0, arr.ind = T)]))){
      coho_right_dt        <- cohort_time_corr[which(rownames(cohort_time_corr) %in% rownames(tab_Nc_Nr)),]
      ID_coho_abund_mature <- which(coho_right_dt==co, arr.ind = T)
      if(length(ID_coho_abund_mature)>0) {
        if(any(as.numeric(ID_coho_abund_mature[,"row"]) <= nrow(tab_Nc_Nr) & 
               as.numeric(ID_coho_abund_mature[,"col"]) <= ncol(tab_Nc_Nr))){
          id_keep <- which(as.numeric(ID_coho_abund_mature[,"row"]) <= nrow(tab_Nc_Nr) & 
                             as.numeric(ID_coho_abund_mature[,"col"]) <= ncol(tab_Nc_Nr))
          Nc[co] <- sum(tab_Nc_Nr[ID_coho_abund_mature[as.numeric(id_keep),]])
        }else { Nc[co] <- 0}
      }else {Nc[co] <- 0}
    }
    Nc <- Nc[which(names(Nc) %in% rownames(tab_Nc_Nr))]
    Nc_only_mature[[sp]] <- Nc
  }
  return(Nc_only_mature)
}


# Nc = sum of mature individuals from first age at maturity * survival at each time step 
Nc_mature_times_survival <- function(DATA_NE, species){
  tab_Nc_Nr_all <- DATA_NE$table_mature
  surviv        <- DATA_NE$survivorship_curve
  new_first_age <- DATA_NE$new_first_age
  Nc_mature_times_survival <- list()
  for (sp in species){
    N_first_age  <- tab_Nc_Nr_all[,which(colnames(tab_Nc_Nr_all) == new_first_age[sp]), sp]
    surv_i       <- surviv[[sp]][[1]]
    cohorts_surv <- surv_i[,which(stringr::str_remove(pattern = "cohort_", string = colnames(surv_i)) %in% names(N_first_age))]
    cohorts_surv_sum <- colSums(cohorts_surv)
    if (any(N_first_age<0)) { N_first_age[which(N_first_age<0)] <- 0}
    Nc_mature_times_survival[[sp]] <- N_first_age*cohorts_surv_sum
  }
  return(Nc_mature_times_survival)
}

# Nc = like Kuparinen et al 2016 = sum of mature individuals at first age at maturity * survival at each time step 
Nc_tot_agemat_times_survival <- function(DATA_NE, species){
  tab_Nc_Nr_all <- DATA_NE$abundanceByAge
  surviv        <- DATA_NE$survivorship_curve
  new_first_age <- DATA_NE$new_first_age
  Nc_tot_agemat_times_survival <- list()
  for (sp in species){
    N_first_age  <- tab_Nc_Nr_all[[sp]][1,]
    surv_i       <- surviv[[sp]][[1]]
    cohorts_surv <- surv_i
    cohorts_surv_sum <- colSums(cohorts_surv)
    if (any(N_first_age<0)) { N_first_age[which(N_first_age<0)] <- 0}
    Nc_tot_agemat_times_survival[[sp]] <- N_first_age*cohorts_surv_sum
  }
  return(Nc_tot_agemat_times_survival)
}

