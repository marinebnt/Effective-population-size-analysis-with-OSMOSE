# GOAL  : create new configuration files each associated to a different repo


####### CHANGING VARIABLES WITH CONF CHOSEN ############
output_repo = "new_config_Fmax_RCP"
main        = "New_gpu_conf_RCP.R" # "New_gpu_conf.R" # "New_param_Eco_Bioen_neutral.R"#_neutral.R"
  
### inputs
source(file.path("R/functions/source.R")) 
repo        = "Simu_Fisheries"
multip_name = "Fmax_"
multipliers = c(0.5,0.75,1,1.25,1.5)

### FISHING RATE DATA ############

path_raw    <- c("../osmose-ns/osmose-ns-v4.4.1_clim/fishing_rate") # raw fishing rate
nspecies    <- get_par(conf, "simulation.nspecies")
species     <- get_par(conf, "species.name")[1:nspecies]
namessp     <- stringr::str_remove(names(species), pattern = "species.name.sp")





### RUN SIMULATION FOR VARIABLE FMAX #####
if(!dir.exists(paste0(output_repo))){dir.create(paste0(output_repo))}
list_fish_rate <- list.files(path_raw, pattern = ".csv")
list_fish_rate <- list_fish_rate[order(as.numeric(unlist(stringr::str_split(unlist(stringr::str_split(list_fish_rate,
                                 pattern = "\\."))[seq(7, length(unlist(stringr::str_split(list_fish_rate, 
                                 pattern = "\\."))), by=8)],
                                 pattern="sp"))[seq(2,16*2, by=2)]))]
for (sp in 1:nspecies){
  sp_i   <- as.character(species[sp])
  sp_id  <- namessp[sp]
  if (sp_i == "Shrimp") next
  cat("species =", sp_i, "\n")
  
  file <- grep(paste0(sp_id, ".csv"), list_fish_rate)
  fish_rate_file <- read.csv(paste0(path_raw, "/", list_fish_rate[file[1]]))
  colnames(fish_rate_file) <- stringr::str_remove(colnames(fish_rate_file), "X")
  
  for (Fmax_m in multipliers){
    # what parameter value change
    cat("Fmax_m =", Fmax_m, "\n")
    # F_max_m_val  =  Fmax_i*Fmax_m #[grep(sp_id, names(Fmax))[1]]
    new_fish_rate_matFMAX <- fish_rate_file
    new_fish_rate_matFMAX[, -1] <- Fmax_m*new_fish_rate_matFMAX[, -1]
    # if (Fmax_m == 1) {new_fish_rate_mat[, -1] <- new_fish_rate_mat[, -1]/Fmax_i}
    
    
    for (Fstart_m in 1){
      name_infile <- paste0("F_Fstart", Fstart_m, "_Fmax", Fmax_m, "_sp", sp_id, ".csv")
      # what parameter value change
      cat("Fstart_m =", Fstart_m, "\n")
      # F_max_m_val  =  Fmax_i*Fmax_m #[grep(sp_id, names(Fmax))[1]]
      
      # increase or decrease (translation) of the first age with fishing mortality
      ID_no_zero     <- which(new_fish_rate_matFMAX[1,] != 0)
      DIVIDER        <- 1*(ncol(fish_rate_file) != c(as.numeric(colnames(fish_rate_file)[ncol(fish_rate_file)])+2))+1
      FracDIV        <- 1/DIVIDER
      first_length   <- (ID_no_zero[1]-2)/DIVIDER
      first_length_m <- first_length*Fstart_m
      dec <- first_length_m - floor(first_length_m)
      if(c(c(dec < 0.5 & dec > 0) || c(dec < 1 & dec > 0.5)) && c(DIVIDER != 1)){first_length_m = first_length_m + 0.25}
      if (DIVIDER == 1 && dec>0){first_length_m <- round(first_length_m + 0.1)}
      SIGN                <- sign(as.numeric(rev(colnames(new_fish_rate_matFMAX[, -1]))[1]) - first_length_m)
      names_1             <- seq(0,first_length_m-FracDIV, by=FracDIV)
      if (SIGN<1){
        names_2           <- seq(names_1[length(names_1)]+1,  first_length_m, by=FracDIV)
      }else{
        names_2           <- seq(first_length_m, as.numeric(rev(colnames(new_fish_rate_matFMAX[, -1]))[1]), by = SIGN*FracDIV)
      }
      df_1                <- matrix(rep(0, nrow(new_fish_rate_matFMAX)*length(names_1)), ncol=length(names_1), nrow=nrow(new_fish_rate_matFMAX))
      df_2                <- as.matrix(new_fish_rate_matFMAX[, ID_no_zero])
      missing_dimensions  <- length(names_2)-ncol(df_2)
      if(missing_dimensions>0) {df_2 <- cbind(df_2, matrix(rep(df_2[,ncol(df_2)], missing_dimensions), ncol=missing_dimensions))}
      if(missing_dimensions<0) {df_2 <- df_2[, -c((ncol(df_2)+missing_dimensions+1) : ncol(df_2))]}
      
      new_fish_rate_matFSTART           <- cbind(df_1, df_2)
      colnames(new_fish_rate_matFSTART) <- c(names_1, names_2)
    
      if(file.exists(paste0(output_repo, "/", name_infile))){file.remove(paste0(output_repo, "/", name_infile))}
      write.csv(file=paste0(output_repo, "/", name_infile), new_fish_rate_matFSTART, row.names = F)
      
      name_file <- paste0(output_repo, "/conf_Fmax", Fmax_m, "_Fstart", Fstart_m, "_sp", sp_id, ".R")    
      file_path <- paste0(name_file)
      file.copy(paste0(main), 
                file_path, overwrite = T)
      lines <- readLines(file_path)
      ID_to_change <- grep(paste0("mortality.fishing.rate.byDt.bySize.file.sp", sp_id), lines)[1]
      lines[ID_to_change]
      lines[ID_to_change] <- paste0("mortality.fishing.rate.byDt.bySize.file.sp", sp_id, " = ", name_infile)
      writeLines(lines, file_path)
      readLines(file_path, n = ID_to_change)
    }
  }
}
