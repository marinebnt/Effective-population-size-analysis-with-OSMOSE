#' data_traits_summary
#' 
#' This function loads and combines all trait data files from the OSMOSE output directory.
#' It binds together all RData files from output/PLOTS/Trait_variables_by_cohort_by_age/RDATA/
#' and keeps the age dimension intact.
#'
#' @param base_path Path to the Trait_variables_by_cohort_by_age directory 
#'        (e.g., "osmose-ns-v4.4.1_clim/output/PLOTS/Trait_variables_by_cohort_by_age")
#' @param Time1 Start of cohort range (default: 150)
#' @param Time2 End of cohort range (default: 160)
#'
#' @return A single dataframe combining all species and trait data with columns:
#'   - age: age class
#'   - variable: trait variable name
#'   - value: trait value (averaged over cohorts Time1-Time2)
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
    
    # duplica <- data_temp  %>%
    #   dplyr::group_by(rep, Fish_press, age, Species, TRAIT) %>%
    #   dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
    #   dplyr::filter(n > 1L) 
    # if(any(row(duplica)>0)) {
    #   cat("Error: duplicates in the dataset")
    #   stop()
    # }
    
    
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
