#' data_traits_summary
#' 
#' This function loads and combines all trait data files from the OSMOSE output directory.
#' It binds together all RData files from output/PLOTS/Trait_variables_by_cohort_by_age/RDATA/
#' and keeps the age dimension intact.
#'
#' @param base_path Path to the Trait_variables_by_cohort_by_age directory 
#'        (e.g., "osmose-ns-v4.4.1_clim/output/PLOTS/Trait_variables_by_cohort_by_age")
#'
#' @return A single dataframe combining all species and trait data with columns:
#'   - age: age class
#'   - variable: trait variable name
#'   - value: trait value (averaged over cohorts Time1-Time2, during the data extraction)
#'   - replicate: replicate number
#'   - Fish_press: fishing pressure (format: FmaxX_FstartY)
#'   - Species: species name
#'   - Fmax: Fmax value
#'   - Fstart: Fstart value
#'
#' @examples
#' \dontrun{
#'   data_summary <- data_traits_summary(
#'     base_path = "osmose-ns-v4.4.1_clim/output/PLOTS/Trait_variables_by_cohort_by_age"
#'   )
#' }

data_traits_summary <- function(base_path, pathout) {
  
  # Initialize combined dataframe
  data_combined <- NULL
  
  # Get path to RDATA subfolder
  rdata_path <- file.path(pathout)
  
  # Check if directory exists
  if (!dir.exists(rdata_path)) {
    stop(paste("Directory does not exist:", rdata_path))
  }
  
  # List all subdirectories in RDATA folder (e.g., Fmax1.5_Fstart1.5_spAtlanticCod)
  subdirs <- list.dirs(file.path(here::here(), base_path, "RDATA"), recursive = FALSE, full.names = TRUE)
  
  if (length(subdirs) == 0) {
    warning(paste("No subdirectories found in:", file.path(here::here(), base_path, "RDATA")))
    return(NULL)
  }
  
  cat("Found", length(subdirs), "files/configuration directories\n")
  
  # Loop through each subdirectory
  for (subdir in subdirs) {
    
    # Extract species name from folder name
    folder_name <- basename(subdir)
    spsuffix    <- stringr::str_locate(pattern = "sp", folder_name)
    sp_modelled <- stringr::str_sub(start = (spsuffix[2]+1), end=nchar(folder_name), string = folder_name)  # keep underscores
    
    for (fi in list.files(subdir, full.names = T)){

      if (!file.exists(fi)) {
        warning(paste("RData file not found:", fi))
        next
      }

      # Load the data
      cat("Loading:", folder_name, "\n")
      data_temp <- get(load(fi))
      data_temp <- data_temp[which(data_temp$Species == sp_modelled),]
      
      if (!("Fmax" %in% colnames(data_temp))) {
        fmax <- unlist(stringr::str_split(pattern="Fmax", fi))
        fmax <- fmax[grep("summary", fmax)]
        fmax <- unlist(stringr::str_split(pattern="_Fstart", fmax))[1]
        data_temp$Fmax <- fmax
      }      
      if (!("Fstart" %in% colnames(data_temp))) {
        fstart <- unlist(stringr::str_split(pattern="Fstart", fi))
        fstart <- fstart[grep("summary", fstart)]
        fstart <- unlist(stringr::str_split(pattern="_", fstart))[1]
        data_temp$Fstart <- fstart
      }
      
      # Check if data is NULL or empty
      if (is.null(data_temp)) {
        warning(paste("Empty data in:", folder_name))
        next
      }
      
      # Ensure data is a dataframe (not a list)
      if (is.list(data_temp) && !is.data.frame(data_temp)) {
        # If it's a list of dataframes, combine them
        data_temp <- do.call(rbind, lapply(data_temp, function(x) {
          if (is.data.frame(x)) x else NULL
        }))
      }
      
      # Check if we have the expected columns
      if (!is.data.frame(data_temp)) {
        warning(paste("Data is not a dataframe in:", folder_name))
        next
      }
      
      # Combine with main dataset
      if (is.null(data_combined)) {
        data_combined <- data_temp
      } else {
        # Use bind_rows to handle potential column mismatches
        data_combined <- dplyr::bind_rows(data_combined, data_temp)
      }
    }
  }
  
  
  # duplica <- data_combined  %>%
  #   dplyr::group_by(rep, Fish_press, age, Species, TRAIT) %>%
  #   dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
  #   dplyr::filter(n > 1L) 
  # if(nrow(duplica)>0) {
  #   cat("Error: duplicates in the dataset")
  #   stop()
  # }
  
  if (is.null(data_combined)) {
    warning("No data was successfully loaded")
    return(NULL)
  }
  
  # Ensure Fish_press column exists
  if (!("Fish_press" %in% colnames(data_combined))) {
    if (all(c("Fmax", "Fstart") %in% colnames(data_combined))) {
      data_combined$Fish_press <- paste0("Fmax", data_combined$Fmax, "_Fstart", data_combined$Fstart)
    }
  }
  
  # Select and order columns (keep only relevant ones)
  desired_cols <- c("age", "variable", "value", "replicate", "Fish_press", "Species", "Fmax", "Fstart", "cohort", "time", "replicate", "rep", "TRAIT")
  available_cols <- intersect(desired_cols, colnames(data_combined))
  data_combined <- data_combined[, available_cols]
  
  cat("\nSuccessfully combined data:\n")
  cat("  - Total rows:", nrow(data_combined), "\n")
  cat("  - Unique species:", length(unique(data_combined$Species)), "\n")
  cat("  - Unique variables:", length(unique(data_combined$variable)), "\n")
  
  # Save combined data to RDATA file
  # Extract folder type from base_path (e.g., "Trait_variables_by_cohort_by_age" -> "cohort_age")
  folder_name <- basename(base_path)
  suffix <- gsub("Trait_variables_by_", "", folder_name)
  suffix <- gsub("_", "_", suffix)  # keep underscores
  
  # Build save path
  if (!dir.exists(pathout)) {
    dir.create(pathout, recursive = TRUE)
  }
  
  save_file <- file.path(pathout, paste0("summary_trait_", suffix, ".RData"))
  save(data_combined, file = save_file)
  cat("\nData saved to:", pathout, "\n")
  
  return(data_combined)
}




## RAW Ne & associated variables, dimensions : cohort, age, replicate, Fstart, Fmax, Species
### produce dataset reading multiple .RData to make 1 dataframe, stored as a .RData file ----------

dataset_Ne_mean_cohorts <- function(species, FORCE, pathdata, pathout, Time1, Time2){
  
  base_dir <- paste0(pathdata, "INIT")
  folders  <- list.dirs(paste0(here::here(), "/", base_dir), recursive = F) # files INIT2 with raw Ne datasets
  full_dataset = list()
  full_dataset_age = list()
  name_out = ifelse(dataset=="data_Ne", "Complete_INIT2", "Complete_INIT2_raw")
  
  for (f_i in folders){ # 374,362
    cat(f_i, "\n")
    
    Fmax   <- as.numeric(stringr::str_split(stringr::str_split(f_i, pattern="Fmax")[[1]][2], pattern="_Fstart")[[1]][1])
    Fstart <- as.numeric(stringr::str_split(stringr::str_split(f_i, pattern="Fstart")[[1]][2], pattern="_sp")[[1]][1])
    
    if (Fstart != 1) next
    
    cat("Loading   ",  f_i)
    
    load(list.files(f_i, full.names = T))
    dataNe <- INIT2[[dataset]]
    rm(INIT2)
    
    for (sp in names(dataNe)[1]){
      if (file.exists(paste0(pathout, "/", name_out, sp_i, ".RData")) && FORCE !=T){
        next
      }else{
        
        cat(sp, "\n")
        
        if(dataset == "data_Ne") {
          f_i_sp <- dataNe[[sp]]
          f_i_sp[["NT"]] <- do.call(cbind, (lapply(f_i_sp$abundanceByAge, FUN=function(x) colSums(x))))
        }
        if(dataset != "data_Ne") f_i_sp <- dataNe
        
        for (dim in names(f_i_sp)){
          if (dim == "cohort_time_step_corresponding"){next}
          f_i_sp_d <- f_i_sp[[dim]]
          
          if(dataset != "data_Ne") {
            sp = names(f_i_sp[[dim]])[1]
            f_i_sp_d <- list()
            for (rep in 1:dim(f_i_sp[[dim]][[sp]])[3]){ # all the dimensions variables are considered from age 1, thus Ne is estimated by variables from age 1, except for Numerator_N_alpha
              f_i_sp_d[[paste0("rep_", rep)]] <- f_i_sp[[dim]][[sp]][which(as.numeric(rownames(f_i_sp[[dim]][[sp]]))>=1), ,paste0("rep_", rep)]
            }
          }
          
          if (is.data.frame(f_i_sp_d) || is.matrix(f_i_sp_d)){
            
            if (is.matrix(f_i_sp_d)) f_i_sp_d <- as.data.frame(f_i_sp_d)
            f_i_sp_d$cohort <- rownames(f_i_sp_d)
            f_i_sp_d$cohort = stringr::str_remove(as.matrix(f_i_sp_d[, "cohort"]), pattern="cohort_")
            f_i_sp_d <- pivot_longer(f_i_sp_d, names_to = "replicate", values_to = "value", cols=-cohort)
            f_i_sp_d$replicate = stringr::str_remove(as.matrix(f_i_sp_d[, "replicate"]), pattern="rep_")
            f_i_sp_d <- f_i_sp_d[which(f_i_sp_d$cohort %in% c(Time1:Time2)), ]
            f_i_sp_d$Species <- sp
            f_i_sp_d$variable <- dim
            f_i_sp_d$Fmax <- as.numeric(Fmax)
            f_i_sp_d$Fstart <- as.numeric(Fstart)
            
            if (any(na.omit(f_i_sp_d$value)==0)) f_i_sp_d[(which(f_i_sp_d$value==0)), "value"] <- NA
            
            # harmonic mean
            data_Ne_Nalpha = NULL
            vector_harmonic_mean = c("Ne_N_alpha", "Numerator_N_alpha", "NT")#, "generation_length", "Vk_point")
            vector_harmonic_mean <- c(vector_harmonic_mean, improve_vriable_names(vector_harmonic_mean))
            if (any(f_i_sp_d$variable %in% vector_harmonic_mean & "cohort" %in% colnames(f_i_sp_d))){
              data_Ne_Nalpha <- f_i_sp_d[which(f_i_sp_d$variable %in% vector_harmonic_mean), ]
              f_i_sp_d   <- f_i_sp_d[-which(f_i_sp_d$variable %in% vector_harmonic_mean), ]
              data_Ne_Nalpha <- data_Ne_Nalpha %>%
                group_by(Fmax, Fstart, variable, replicate, Species) %>%
                summarise(value = psych::harmonic.mean(value, na.rm = TRUE))%>%
                ungroup()
            }
            f_i_sp_d <- rbind(f_i_sp_d, data_Ne_Nalpha)
            # arithmetic mean mean
            f_i_sp_d <- f_i_sp_d %>%
              group_by(Fmax, Fstart, variable, replicate, Species) %>%
              summarise(value = mean(value, na.rm = TRUE))%>%
              ungroup()
            f_i_sp_d$age <- NA
              
            if (is.null(full_dataset_age[[sp]])){
              full_dataset_age[[sp]] <- f_i_sp_d
            } else {
              full_dataset_age[[sp]] <- rbind(full_dataset_age[[sp]], f_i_sp_d)
            }
            
          } else { # (!is.data.frame(f_i_sp_d) & !is.matrix(f_i_sp_d))
            
            f_i_sp_d = lapply(names(f_i_sp_d), FUN = function(x) {
              f_i_sp_d[[x]] <- as.data.frame(f_i_sp_d[[x]][,which(colnames(f_i_sp_d[[x]]) %in% paste0("cohort_", Time1:Time2))])
              f_i_sp_d[[x]][, "replicate"] = stringr::str_remove(x, pattern="rep_")
              f_i_sp_d[[x]][, "age"] = rownames(f_i_sp_d[[x]])
              return(f_i_sp_d[[x]])
            })
            f_i_sp_d = do.call(rbind, f_i_sp_d)
            
            f_i_sp_d <- pivot_longer(f_i_sp_d, names_to = "cohort", values_to = "value", cols=-c(replicate, age))
            f_i_sp_d[, "cohort"] = stringr::str_remove(as.matrix(f_i_sp_d[, "cohort"]), pattern="cohort_")
            f_i_sp_d$Species <- sp
            f_i_sp_d$variable <- dim
            f_i_sp_d$Fmax <- as.numeric(stringr::str_split(stringr::str_split(f_i, pattern="Fmax")[[1]][2], pattern="_Fstart")[[1]][1])
            f_i_sp_d$Fstart <- as.numeric(stringr::str_split(stringr::str_split(f_i, pattern="Fstart")[[1]][2], pattern="_sp")[[1]][1])
            
            if (any(na.omit(f_i_sp_d$value)==0)) f_i_sp_d[(which(f_i_sp_d$value==0)), "value"] <- NA
            
            # mean of the cohorts : arithmetic or harmonic
            # harmonic mean
            data_Ne_Nalpha = NULL
            vector_harmonic_mean = c("Ne_N_alpha", "Numerator_N_alpha", "NT")#, "generation_length", "Vk_point")
            vector_harmonic_mean <- c(vector_harmonic_mean, improve_vriable_names(vector_harmonic_mean))
            if (any(f_i_sp_d$variable %in% vector_harmonic_mean & "cohort" %in% colnames(f_i_sp_d))){
              data_Ne_Nalpha <- f_i_sp_d [which(f_i_sp_d$variable %in% vector_harmonic_mean), ]
              f_i_sp_d   <- f_i_sp_d [-which(f_i_sp_d$variable %in% vector_harmonic_mean), ]
              data_Ne_Nalpha <- data_Ne_Nalpha %>%
                group_by(age, Fmax, Fstart, variable, replicate, Species) %>%
                summarise(value = psych::harmonic.mean(value, na.rm = TRUE))%>%
                ungroup()
            }
            f_i_sp_d <- rbind(f_i_sp_d, data_Ne_Nalpha)
            # arithmetic mean mean
            f_i_sp_d <- f_i_sp_d %>%
              group_by(age, Fmax, Fstart, variable, replicate, Species) %>%
              summarise(value = mean(value, na.rm = TRUE))%>%
              ungroup()
            
            
            if (is.null(full_dataset_age[[sp]])){
              full_dataset_age[[sp]] <- f_i_sp_d
            } else {
              full_dataset_age[[sp]] <- rbind(full_dataset_age[[sp]], f_i_sp_d)
            }
            
          } # if list
        }# variable dimension
      }# only if file does not exist yet
    }# species dimension
  }# folder dimension (Fmax*Fstart)
  #----- MERGE raw variables without age dimension AND variables with age dimension ---#
  for (sp_i in names(full_dataset_age)){
    if (file.exists(paste0(pathout, "/", name_out, sp_i, ".RData")) && FORCE !=T){
      next
    }else{
      data_out <- full_dataset_age[[sp_i]] %>%
        dplyr::group_by(Species, replicate, Fstart, Fmax, variable) %>%
        dplyr::summarise(
          value = mean(value, na.rm=T), #stats::weighted.mean(value, w = weight_value, na.rm = TRUE),
          .groups = "drop"
        )
      data_out$Fmax <- as.numeric(data_out$Fmax)
      data_out$Fstart <- as.numeric(data_out$Fstart)
      final_dataset <- rbind(full_dataset[[sp_i]], data_out)
      final_dataset_age <- (full_dataset_age[[sp_i]])
      
      # save all
      complete = list("complete_without_age"=final_dataset, "complete_with_age"=final_dataset_age)
      
      cat("\nNB var: ", length(table(complete$complete_with_age$variable)))
      
      save(complete, file=paste0(pathout, "/", name_out, sp_i, ".RData"))
      
      cat("\nSuccessfully combined data related to Ne by age:\n")
      cat("  - Total rows:", nrow(final_dataset_age), "\n")
      cat("  - Unique species:", length(unique(final_dataset_age$Species)), "\n")
      cat("  - Unique variables:", length(unique(final_dataset_age$variable)), "\n")
      
      cat("\nData saved to:", pathout, "/", name_out, sp_i, ".RData", "\n")
    }
  }
}



# RAW TRAITS, dimensions :  age, replicate, Fish_press, Species, TRAIT

load_trait_df <- function(base_path, pathout, FORCE=F, removeage0=T){
   # load dataset
   data_complete_sp = NULL
   data_complete_c = NULL
   folder_name <- basename(base_path)
   suffix <- gsub("Trait_variables_by_", "", folder_name)
   suffix <- gsub("_", "_", suffix)  # keep underscores
   save_dir  <- file.path(pathout)
   save_file <- file.path(save_dir, paste0("summary_trait_", suffix, ".RData"))
   if (!dir.exists(save_dir)){dir.create(save_dir, recursive = T)}

   if (!file.exists(save_file) || FORCE==T){
     data_complete_c <- data_traits_summary(base_path, pathout)
   } else {
     data_complete_c = get(load(file=save_file))
   }
   
   # removing age 0
   if (removeage0 == T) data_complete_c <- data_complete_c[-which(data_complete_c$age == 0), ]
   # correction
   if (any(grepl(data_complete_c$TRAIT, pattern="fec"))) data_complete_c[grep(data_complete_c$TRAIT, pattern="fec"), "TRAIT"] <- "mean_of_eggNumber_per_age"
   
   # adding natural mortality dimension
   if (any((data_complete_c$TRAIT) == "mean_of_Mnat_per_age")){
     data_complete_c[which((data_complete_c$TRAIT) == "mean_of_Mnat_per_age"), "TRAIT"] <- "mean_of_Madd_per_age"
     data_complete_c <- pivot_wider(data = data_complete_c, names_from = "TRAIT", values_from = "value")
     data_complete_c$mean_of_Mnat_per_age <- data_complete_c$mean_of_Madd_per_age + data_complete_c$mean_of_Mpred_per_age + data_complete_c$mean_of_Mstar_per_age + data_complete_c$mean_of_Mfor_per_age
   }
   if("Fish_press" %in% colnames(data_complete_c)) {
     data_complete_c$Fish_press <- stringr::str_replace(data_complete_c$Fish_press, pattern = "_", replacement = " ")
     data_complete_c$Fish_press <- stringr::str_remove(data_complete_c$Fish_press, pattern = "=")
     data_complete_c$Fish_press <- stringr::str_remove(data_complete_c$Fish_press, pattern = "=")
     data_complete_c$Fmax   <- as.numeric(unlist(stringr::str_remove(unlist(stringr::str_split(data_complete_c$Fish_press, " Fstart"))[seq(1, 2*length(data_complete_c$Fish_press), by=2)], "Fmax")))
     data_complete_c$Fstart <- as.numeric(unlist(stringr::str_split(data_complete_c$Fish_press, " Fstart"))[seq(2, 2*length(data_complete_c$Fish_press), by=2)])
     data_complete_c <- data_complete_c[,-which(colnames(data_complete_c) == "Fish_press")]
   }
   if (!("TRAIT" %in% colnames(data_complete_c))) 
     {data_complete_c <- pivot_longer(data = data_complete_c, names_to = "TRAIT", values_to = "value", 
                                   cols=-colnames(data_complete_c)[which(colnames(data_complete_c) %in% 
                                                                c("Fstart", "Fmax", "Species", "age", "replicate"))])}
   
   # changing colnames to match expectations
   colnames(data_complete_c)[which(colnames(data_complete_c) == "TRAIT")] <- "variable"

   return(data_complete_c)
}

# FUNCTIONS TO CLEA DATASETS
when_variables_missing_age_dim <- function(dataset, add_missing_ages = T){
  nvarinit <- length(table(dataset$variable))
  dataset_age <- (as.numeric(as.matrix(dataset[which(dataset$variable %in% c("eggNumber", "mean_of_eggNumberByAge_per_age")), "age"])))
  min_age <- min(dataset_age)
  if (any(as.numeric(as.matrix(dataset$age))<min_age) || any(is.na(as.matrix(dataset$age)))) {
    if(any(is.na(as.matrix(dataset$age))) & add_missing_ages == F) {
      cat("\nRemplacing NA in age by min_age only")
      data_na <- dataset[which(is.na(as.matrix(dataset$age))), ]
      data_na$age <- min_age
      data_not_na <- dataset[-which(is.na(as.matrix(dataset$age))), ]
      dataset <- rbind(data_na, data_not_na)
    }
    if(any(is.na(as.matrix(dataset$age))) & add_missing_ages == T) {
      cat("\nRemplacing NA in age by a looped vector with all the ages from maturation to mortality")
      data_rep <- dataset[which(is.na(as.matrix(dataset$age))), ]
      data_not_rep <- dataset[-which(is.na(as.matrix(dataset$age))), ]
      dataset_repet_i = NULL
      for (x in min_age:max(as.numeric(as.matrix(dataset$age)), na.rm=T)) {
        tempo <- data_rep
        tempo$age <- x
        dataset_repet_i <- rbind(dataset_repet_i, tempo)
      }
      dataset <- rbind(dataset_repet_i, data_not_rep)
    }
    if(any(as.numeric(as.matrix(dataset$age))<min_age)) {
      cat("\nCorrecting error in the dimension of the input dataset : certain outputs have the age < maturation age\n")
      data_wrong_range_age  <- dataset[which(as.numeric(as.matrix(dataset$age))<min_age), "variable"]
      data_good_range_age   <- dataset[-which(dataset$variable %in% as.matrix(unique(data_wrong_range_age))), ]
      data_good_range_age_i <- as.matrix(data_good_range_age$variable)[1]
      data_good_range_age_i <- dataset[which(dataset$variable == data_good_range_age_i), ]
      ages_i_ok <- as.numeric(as.matrix(data_good_range_age_i$age))
      for(vi in names(table(data_wrong_range_age))){
        data_wrong_range_age_i <- dataset[which(dataset$variable == vi), ]
        ages_i <- as.numeric(as.matrix(data_wrong_range_age_i$age))
        if ((max(ages_i)-min(ages_i)) == (max(ages_i_ok)-min(ages_i_ok))){
          data_wrong_range_age_i$age = as.numeric(as.matrix(data_wrong_range_age_i$age))-min(ages_i)+min(ages_i_ok)
        } else cat("=====> cannot be fixed, normal if add_missing_ages == F <=============")
        data_good_range_age <- rbind(data_good_range_age, data_wrong_range_age_i)
      }
      dataset <- data_good_range_age
    }
  }
  nvarend <- length(table(dataset$variable))
  if(nvarend!=nvarend) stop("Dimensions have changed: lost variable(s)")
  return(dataset)
}
when_null_or_low_abund <- function(dataset){
  tab <- table(dataset$variable)
  var_less_rep <- names(tab)[which(tab<max(tab))]
  dataset_t <- pivot_wider(dataset, names_from = "variable", values_from = "value")
  if (any(is.na(as.matrix(dataset_t[,grep(colnames(dataset_t), pattern="abundanceByAge")])))) {
    dataset_t <- dataset_t[-which(is.na(as.matrix(dataset_t[,grep(colnames(dataset_t), pattern="abundanceByAge")]))),]
    dataset   <- pivot_longer(dataset_t, cols = -c("age", "replicate", "Fmax", "Fstart", 'Species'), 
                              names_to = "variable", values_to = "value")
  }
  if (!(length(var_less_rep)==0)) {
    data_t  <- dataset[-which(dataset$variable %in% var_less_rep), ]
  } else {
    data_t <- dataset
  }
  dataset <- rbind(data_t, na.omit(dataset[which(dataset$variable %in% var_less_rep), ]))
  if (any(as.matrix(dataset_t[,grep(colnames(dataset_t), pattern="abundanceByAge")])<100, na.rm=T)){
    dataset_t[which(as.matrix(dataset_t[,grep(colnames(dataset_t), pattern="abundanceByAge")])<100), -which(colnames(dataset_t) %in%
                                                                                                              c("age", "replicate", "Fmax", "Fstart", 'Species'))] <- NA
    dataset <- pivot_longer(dataset_t, cols = -c("age", "replicate", "Fmax", "Fstart", 'Species'), names_to = "variable", values_to = "value")
  }
  if (!(length(var_less_rep)==0)) {
    data_t  <- dataset[-which(dataset$variable %in% var_less_rep), ]
  } else {
    data_t <- dataset
  }
  dataset <- rbind(data_t, na.omit(dataset[which(dataset$variable %in% var_less_rep), ]))
  return(dataset)
}
what_ages_have_null_dead_nb <- function(dataset){
  # remove ages that are not included
  dataset$age <- as.numeric(dataset$age)
  ages_to_remove <- table(dataset[which(dataset$variable == "Dead" & dataset$value ==0),"age"])
  if (dim(ages_to_remove)>0) ages_to_remove <- ages_to_remove[which(ages_to_remove == max(ages_to_remove))]
  if (any(dataset$age %in% as.numeric(names(ages_to_remove)))) dataset <- dataset[-which(dataset$age %in% as.numeric(names(ages_to_remove))), ]
  if (any(which(is.na(dataset$age)))) dataset$age[which(is.na(dataset$age))] <- min_age
  return(dataset)
}
         