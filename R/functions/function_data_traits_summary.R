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

data_traits_summary <- function(base_path) {
  
  # Initialize combined dataframe
  data_combined <- NULL
  
  # Get path to RDATA subfolder
  rdata_path <- file.path(base_path, "RDATA")
  
  # Check if directory exists
  if (!dir.exists(rdata_path)) {
    stop(paste("Directory does not exist:", rdata_path))
  }
  
  # List all subdirectories in RDATA folder (e.g., Fmax1.5_Fstart1.5_spAtlanticCod)
  subdirs <- list.dirs(rdata_path, recursive = FALSE, full.names = TRUE)
  
  if (length(subdirs) == 0) {
    warning(paste("No subdirectories found in:", rdata_path))
    return(NULL)
  }
  
  cat("Found", length(subdirs), "species/configuration directories\n")
  
  # Loop through each subdirectory
  for (subdir in subdirs) {
    
    # Look for data_traits_summary.RData file
    rdata_file <- file.path(subdir, "data_traits_summary.RData")
    
    if (!file.exists(rdata_file)) {
      warning(paste("RData file not found:", rdata_file))
      next
    }
    
    # Extract species name from folder name
    folder_name <- basename(subdir)
    spsuffix <- stringr::str_locate(pattern = "sp", folder_name)
    sp_modelled <- stringr::str_sub(start = (spsuffix[2]+1), end=nchar(folder_name), string = folder_name)  # keep underscores
    
    # Load the data
    cat("Loading:", folder_name, "\n")
    data_temp <- get(load(rdata_file))
    data_temp <- data_temp[which(data_temp$Species == sp_modelled),]
    
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
  save_dir <- file.path(dirname(dirname(base_path)), "RDATA")
  if (!dir.exists(save_dir)) {
    dir.create(save_dir, recursive = TRUE)
  }
  
  save_file <- file.path(save_dir, paste0("summary_trait_", suffix, ".RData"))
  save(data_combined, file = save_file)
  cat("\nData saved to:", save_file, "\n")
  
  return(data_combined)
}




## RAW Ne & associated variables, dimensions : cohort, age, replicate, Fstart, Fmax, Species
### produce dataset reading multiple .RData to make 1 dataframe, stored as a .RData file ----------

dataset_Ne_mean_cohorts <- function(species, FORCE, folders, path, Time1, Time2){
  
  full_dataset = list()
  full_dataset_age = list()
  name_out = ifelse(dataset=="data_Ne", "Complete_INIT2", "Complete_INIT2_raw")
  
  for (f_i in folders[[ifelse(dataset=="data_Ne", 1, 2)]]){
    cat(f_i, "\n")
    
    Fmax   <- as.numeric(stringr::str_split(stringr::str_split(f_i, pattern="Fmax")[[1]][2], pattern="_Fstart")[[1]][1])
    Fstart <- as.numeric(stringr::str_split(stringr::str_split(f_i, pattern="Fstart")[[1]][2], pattern="_sp")[[1]][1])
    
    if (Fstart != 1) next
    
    load(list.files(f_i, full.names = T))
    dataNe <- INIT2[[dataset]]
    rm(INIT2)
    
    for (sp in names(dataNe)[1]){
      if (file.exists(paste0(path, "/", name_out, sp_i, ".RData")) && FORCE !=T){
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
            
            # harmonic mean
            data_Ne_Nalpha = NULL
            vector_harmonic_mean = c("Ne_N_alpha", "Numerator_N_alpha", "NT")
            vector_harmonic_mean <- c(vector_harmonic_mean, improve_vriable_names(vector_harmonic_mean))
            if (any(f_i_sp_d$variable %in% vector_harmonic_mean & "cohort" %in% colnames(f_i_sp_d))){
              data_Ne_Nalpha <- f_i_sp_d [which(f_i_sp_d$variable %in% vector_harmonic_mean), ]
              f_i_sp_d   <- f_i_sp_d[-which(f_i_sp_d$variable %in% vector_harmonic_mean), ]
              data_Ne_Nalpha <- data_Ne_Nalpha %>%
                group_by(Fmax, Fstart, variable, replicate, Species) %>%
                summarise(value = psych::harmonic.mean(value, na.rm = TRUE))
            }
            f_i_sp_d <- rbind(f_i_sp_d, data_Ne_Nalpha)
            # arithmetic mean mean
            f_i_sp_d <- f_i_sp_d %>%
              group_by(Fmax, Fstart, variable, replicate, Species) %>%
              summarise(value = mean(value, na.rm = TRUE))
            
            if (is.null(full_dataset_age[[sp]])){
              full_dataset_age[[sp]] <- f_i_sp_d
            } else {
              full_dataset_age[[sp]] <- rbind(full_dataset_age[[sp]], f_i_sp_d)
            }
            
          } # if matrix or df
          if (!is.data.frame(f_i_sp_d) & !is.matrix(f_i_sp_d)){
            
            f_i_sp_d = lapply(names(f_i_sp_d), FUN = function(x) {
              f_i_sp_d[[x]] <- as.data.frame(f_i_sp_d[[x]])
              f_i_sp_d[[x]][, "replicate"] = stringr::str_remove(x, pattern="rep_")
              f_i_sp_d[[x]][, "age"] = rownames(f_i_sp_d[[x]])
              return(f_i_sp_d[[x]])
            })
            f_i_sp_d = do.call(rbind, f_i_sp_d)
            
            f_i_sp_d <- pivot_longer(f_i_sp_d, names_to = "cohort", values_to = "value", cols=-c(replicate, age))
            f_i_sp_d[, "cohort"] = stringr::str_remove(as.matrix(f_i_sp_d[, "cohort"]), pattern="cohort_")
            f_i_sp_d <- f_i_sp_d[which(f_i_sp_d$cohort %in% c(Time1:Time2)), ]
            f_i_sp_d$Species <- sp
            f_i_sp_d$variable <- dim
            f_i_sp_d$Fmax <- as.numeric(stringr::str_split(stringr::str_split(f_i, pattern="Fmax")[[1]][2], pattern="_Fstart")[[1]][1])
            f_i_sp_d$Fstart <- as.numeric(stringr::str_split(stringr::str_split(f_i, pattern="Fstart")[[1]][2], pattern="_sp")[[1]][1])
            
            
            # mean of the cohorts : arithmetic or harmonic
            # harmonic mean
            data_Ne_Nalpha = NULL
            vector_harmonic_mean = c("Ne_N_alpha", "Numerator_N_alpha", "NT",
                                     "Ne_N_alpha/Numerator_N_alpha")
            vector_harmonic_mean <- c(vector_harmonic_mean, improve_vriable_names(vector_harmonic_mean))
            if (any(f_i_sp_d$variable %in% vector_harmonic_mean & "cohort" %in% colnames(f_i_sp_d))){
              data_Ne_Nalpha <- f_i_sp_d [which(f_i_sp_d$variable %in% vector_harmonic_mean), ]
              f_i_sp_d   <- f_i_sp_d [-which(f_i_sp_d$variable %in% vector_harmonic_mean), ]
              data_Ne_Nalpha <- data_Ne_Nalpha %>%
                group_by(age, Fmax, Fstart, variable, replicate, Species) %>%
                summarise(value = psych::harmonic.mean(value, na.rm = TRUE))
            }
            f_i_sp_d <- rbind(f_i_sp_d, data_Ne_Nalpha)
            # arithmetic mean mean
            f_i_sp_d <- f_i_sp_d %>%
              group_by(age, Fmax, Fstart, variable, replicate, Species) %>%
              summarise(value = mean(value, na.rm = TRUE))
            
            
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
    if (file.exists(paste0(path, "/", name_out, sp_i, ".RData")) && FORCE !=T){
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
      
      save(complete, file=paste0(path, "/", name_out, sp_i, ".RData"))
    }
  }
}



# RAW TRAITS, dimensions :  age, replicate, Fish_press, Species, TRAIT

load_trait_df <- function(base_path){
   # load dataset
   data_complete_sp = NULL
   data_complete_c = NULL
   folder_name <- basename(base_path)
   suffix <- gsub("Trait_variables_by_", "", folder_name)
   suffix <- gsub("_", "_", suffix)  # keep underscores
   save_dir <- file.path(dirname(dirname(base_path)), "RDATA")
   save_file <- file.path(save_dir, paste0("summary_trait_", suffix, ".RData"))
   if (!dir.exists(save_dir)){dir.create(save_dir, recursive = T)}
   FORCE=F
   if (!file.exists(save_file) || FORCE==T){
     data_complete_c <- data_traits_summary(base_path)
   } else {
     data_complete_c = get(load(file=save_file))
   }
   
   # removing age 0
   data_complete_c <- data_complete_c[-which(data_complete_c$age == 0), ]
   # correction
   data_complete_c[grep(data_complete_c$TRAIT, pattern="fec"), "TRAIT"] <- "mean_of_eggNumber_per_age"
   
   # adding natural mortality dimension
   data_complete_c[which((data_complete_c$TRAIT) == "mean_of_Mnat_per_age"), "TRAIT"] <- "mean_of_Madd_per_age"
   data_complete_c <- pivot_wider(data = data_complete_c, names_from = "TRAIT", values_from = "value")
   data_complete_c <- data_complete_c[,-grep(colnames(data_complete_c), pattern="mean_of_Dis_per_age")]
   data_complete_c <- data_complete_c[,-grep(colnames(data_complete_c), pattern="mean_of_Z_per_age")]
   data_complete_c$mean_of_Mnat_per_age <- data_complete_c$mean_of_Madd_per_age + data_complete_c$mean_of_Mpred_per_age + data_complete_c$mean_of_Mstar_per_age
   data_complete_c <- pivot_longer(data = data_complete_c, names_to = "TRAIT", values_to = "value", cols=-colnames(data_complete_c)[1:4])
   
   # changing colnames to match expectations
   colnames(data_complete_c)[which(colnames(data_complete_c) == "TRAIT")] <- "variable"
   data_complete_c$Fmax   <- as.numeric(unlist(stringr::str_remove(unlist(stringr::str_split(data_complete_c$Fish_press, " Fstart"))[seq(1, 2*length(data_complete_c$Fish_press), by=2)], "Fmax=")))
   data_complete_c$Fstart <- as.numeric(unlist(stringr::str_split(data_complete_c$Fish_press, " Fstart="))[seq(2, 2*length(data_complete_c$Fish_press), by=2)])
   data_complete_c <- (data_complete_c[,-which(colnames(data_complete_c) == "Fish_press")])

   return(data_complete_c)
}

         
         