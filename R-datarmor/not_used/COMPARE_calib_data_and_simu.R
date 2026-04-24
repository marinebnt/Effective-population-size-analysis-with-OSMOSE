

main="New_gpu_conf.R"
source(file.path("osmose-ns-v4.4.1_clim/R/functions/source.R"))
repo        = "RESTART"
multip_name <- c("y")
multipliers = c("0.5", "1", "1.5")
multipliers = 200

# prepare environment for this particular plot
rep=1
incr_multiplier <- 0
cat("load simulation output")
simulation = paste0(repo, "/output_", multip_name, multipliers)
Ne_version = paste0("plot_biomass", multip_name, multipliers)

ns1 = NULL
outputDir_Simu  = file.path(configDir, "output", simulation)
outputDir_Ne =paste0(configDir, "/output/", "Ne","/PLOT_REPLICATES/")
ns1  = read_osmose(path = outputDir_Simu)
conf = read_osmose(path = paste0(outputDir_Simu, "/ns-configuration.osm"))
nrep = conf$simulation.nsimulation
# multipliers = "output_etaequalsone_bioen_neutral"
trait_dim <- c("biomass", "landings", "meanSize")

fi              <-list.files("../osmose-ns_calibration_data/OSMOSE_NS_calibration/data/", full.names = T)
numbers         <- gregexpr("[0-9]+", fi)
result          <- regmatches(fi, numbers) 
numeric_result  <- as.numeric(unlist(result))
max_age         <- numeric_result[which.max(numeric_result)]

list_comparison <- list()
for (rep in 1:nrep){
  cat("\nREP ", rep, "\n")
  mean_fr <- as.data.frame(matrix(-99, ncol=length(species), nrow=c(max_age+length(trait_dim)-1)))
  colnames(mean_fr) <- species
  rownames(mean_fr) <- c(trait_dim[1], "yield", paste0(trait_dim[3], "ByAge", 0:c(max_age-1)))  
  sd_fr <- as.data.frame(matrix(-99, ncol=length(species), nrow=c(max_age+length(trait_dim)-1)))
  colnames(sd_fr) <- species
  rownames(sd_fr) <- c(trait_dim[1], "yield", paste0(trait_dim[3], "ByAge", 0:c(max_age-1)))
  for (sp in species){
    for (trait in trait_dim){
      cat("\n\t--->\t", sp, "\t", trait, '\n')
      f     <- fi[grep(sp, fi)][grep(trait, fi[grep(sp, fi)])]
      if (length(f)==0){next}
      if (trait == "landings") trait = "yield"
      
      # extract the corresponding data from the simulation output
      if (length(f)>1){
        trait <- paste0(trait, "ByAge")
        data  <- ns1[[trait]][[sp]][,,rep]
        for (fa in f){
          # calibration data
          fr    <- read.csv(fa)
          
          # extract that age data from simulation
          age     <- as.numeric(paste(str_extract_all(string = fa, pattern = '\\(\\d{2}\\)|\\d')[[1]], collapse=""))-1
          trait_i <- paste0(trait, age)
          data_a  <- data[, which(colnames(data) == age)]
          
          # compare simulation and calibration data
          data_trunc <- as.data.frame(data_a)
          data_trunc$time <- as.numeric(rownames(data_trunc))
          if (any(colnames(data_trunc) == "Time")) {
            data_trunc$time <- data_trunc$Time
            data_trunc      <- data_trunc[,-which(colnames(data_trunc) == c("Time"))]}
          # keep only the data corresponding to the calibration data
          data_trunc      <- data_trunc[which(data_trunc$time %in% time), ]
          
          # Compare calibration goal and simulated
          mean_fr[trait_i, sp] <- mean((data_trunc[, -which(colnames(data_trunc) =="time")] - fr[, -which(colnames(fr) =="Time")])/fr[, -which(colnames(fr) =="Time")], na.rm = T)
          sd_fr[trait_i, sp]   <- sd((data_trunc[, -which(colnames(data_trunc) =="time")] - fr[, -which(colnames(fr) =="Time")])/fr[, -which(colnames(fr) =="Time")], na.rm = T)
        }
      }else {
        fr    <- read.csv(f)
        time  <- fr$Time
        data <- ns1[[trait]][,sp, rep]
        if(is.null(data)){
          data     <- list()
          files    <- list.files(outputDir_Simu, full.names = T, recursive = T)
          files    <- files[c(grep(".csv", files))]
          files_rep0    <- files[c(grep("Simu0", files))]
          ID_files <- grep(trait, files_rep0)
          dim_data_ID_files <- unlist(sapply(X = files_rep0[ID_files], 
                                             FUN=function(x) {dim(read.csv(x, skip = c(1)))}))
          ID_file_selected  <- which.min(dim_data_ID_files[1,]) # chose the file the smallest : exclude the data by age or by size ...
          if (nrep>1){
            str_a <- stringr::str_remove(pattern = "0.csv", string = files_rep0[ID_files[ID_file_selected]])
            str_b <- stringr::str_remove(pattern = "0.nc", string = files_rep0[ID_files[ID_file_selected]])
            ID_file_selected_Simu <- grep(c(str_a, str_b)[which.min(c(nchar(str_a), nchar(str_b)))], files)
          } else {ID_file_selected_Simu <- which(files == files_rep0[ID_files[ID_file_selected]])}
          for (ID_i in 1:length(ID_file_selected_Simu)){
            if (length(ID_i)>1){cat("error dim selected outputs")
              stop()}
            ID <- ID_file_selected_Simu[ID_i]
            if(grepl("ByAge", files[ID]) || grepl("BySize", files[ID]) ){
              data_A <- as.data.frame(read.csv(files[ID], skip=c(1)))
              data_A <- data_A[, -c(which(colnames(data_A)=="Age"), which(colnames(data_A)=="Size"))]
              avg_by_species_age <- data_A %>%
                group_by(Time) %>%
                summarise(across(where(is.numeric), \(x) mean(x, na.rm = TRUE)), .groups = "drop")
              data <- avg_by_species_age[,sp]
              if (grepl("yieldN", files[ID])){
                data[1:c(nrow(data)-1), sp] <- data[1:c(nrow(data)-1), sp] - data[2:nrow(data), sp]
              }
            }else {data <- as.data.frame(read.csv(files[ID], skip=c(1)))}
          }
        }
        data_trunc <- as.data.frame(data)
        data_trunc$time <- as.numeric(rownames(data_trunc))
        if (any(colnames(data_trunc) == "Time")) {
          data_trunc$time <- data_trunc$Time
          data_trunc      <- data_trunc[,-which(colnames(data_trunc) == c("Time"))]}
        # keep only the data corresponding to the calibration data
        data_trunc         <- data_trunc[which(data_trunc$time %in% time), ]
        mean_fr[trait, sp] <- mean(((data_trunc[,1] - fr[, -which(colnames(fr) =="Time")])/fr[, -which(colnames(fr) =="Time")]), na.rm = T)
        sd_fr[trait, sp]   <- sd(((data_trunc[,1] - fr[, -which(colnames(fr) =="Time")])/ fr[, -which(colnames(fr) =="Time")]), na.rm = T)
      }
    }
  } 
  list_comparison[[rep]] <- list(MEAN = mean_fr, SD = sd_fr) 
}


save(file="COMPARE_RESTART_REP_WITH_CALIB.RData", list_comparison)


