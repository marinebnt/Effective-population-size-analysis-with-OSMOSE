
function_INIT2 <- function(outputDir_Simu, outputDir_Ne, dim_vec, species, nb_rep, OUTPUT_data, conf, seeding_y, REPO, TEMPORAL=TRUE, HTZ=FALSE, 
      removing_senescence_param=F){

  # effective population sizes
  cat("\n=> prepare osmose data to fit the by-cohort requisites")
  data_Ne_raw <- estimate_for_Ne(outputDir_Simu, outputDir_Ne, dim_vec, species, nb_rep, OUTPUT_data, rep)
  # save.image(file=paste0("env_for_data_Ne.RData"))
  # load(paste0("env_for_data_Ne.RData"))
  cat("\n=> producing Ne data\n")
  data_Ne <- estimate_Ne(outputDir_Ne = outputDir_Ne, data_Ne_shape = data_Ne_raw, species = species, 
                          nb_rep, input_dir=outputDir_Simu, OUTPUT_data, removing_senescence=removing_senescence_param)
  tempoNe <- NA
  if(TEMPORAL) {
    cat("\n=> producing Nb data\n")
    data_Nb <- estimate_Nb(outputDir_Simu, outputDir_Ne, species, nb_rep, OUTPUT_data)
    cat("\n=> producing temporal Ne data\n")
    tempoNe <- estimate_temporalNe(outputDir_Simu, outputDir_Ne, species, conf, seeding_y, OUTPUT_data)}
  
  htz_list <- list()
  if(HTZ) {
    dimensions_outputs = dim(OUTPUT_data$biomass)
    first_y_Ne_tempo <- max(seeding_y, c(years-dimensions_outputs[1]+1))
    cat("\n=> producing htz data\n")
    for (z in species){        
      cat(z)
      htz <- NA
      htz <- data.frame(matrix(ncol=nb_rep, nrow=c(dimensions_outputs[1])))
      colnames(htz) = paste0("rep_", 1:nb_rep)
      for(rep in 1:nb_rep){
      # allele freq data by age
        data <- ncdf4::nc_open(paste0(outputDir_Simu, "/Genetic/ns_observedHtzDistribByAge-", z, "Simu", c(rep-1), ".nc"))
        data_nc_raw <- ncdf4::ncvar_get(data, "observed_heterozygosity")
        data_nc <- data_nc_raw[ ,2, ] # keeping only age 1 (i=1 is the age 0)
        # data_nc <- apply(data_nc_raw, c(1,3), function(x) mean(x, na.rm = TRUE)) # mean over ages 
        # if ( length(dim(data_nc))==4) {trait1_data <- data_nc[, , 1, ] # trait 1 is the neutral trait 
        #   }else {trait1_data <- data_nc} 
        ncdf4::nc_close(data)

        if (length(dim(data_nc))==3) {trait1_data <- data_nc[, 1, ] # trait 1 is the neutral trait 
        }else {trait1_data <- data_nc} 

        mean_htz <- colMeans(trait1_data)
        time <- 1:c(dim(trait1_data)[2])
        htz[,paste0("rep_", rep)] <- mean_htz
        rownames(htz) =  first_y_Ne_tempo:(years)
      }
      htz_list[[z]] <- htz
    }
  }
  if (TEMPORAL)  list_output <- list("data_Ne_raw"=data_Ne_raw, "data_Ne"=data_Ne, "data_Nb"=data_Nb, "tempoNe"=tempoNe, "htz_list"=htz_list)
  if (!TEMPORAL) list_output <- list("data_Ne_raw"=data_Ne_raw, "data_Ne"=data_Ne, "htz_list"=htz_list)
  return(list_output)
}

