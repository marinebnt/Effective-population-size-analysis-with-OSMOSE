# GOAL  : create new configuration files each associated to a different repo

main       = "New_param_Eco_Bioen_neutral.R"#_neutral.R"
source(file.path("osmose-ns-v4.4.1_clim/R/functions/source.R"))
repo        = "Simu_Fisheries"
multip_name = "Fmax_"
multipliers = c(0.5,0.75,1,1.25,1.5)



### FISHING RATE DATA ############

path_raw <- c("osmose-ns-v4.4.1_clim/fishing_rate")

# fish_raw    <- osmose::read_osmose(paste0(path_raw, "/Fishing_rate_raw.R"))
# Fisheries   <- get_par(fish_raw, "fisheries")
# Fmax        <- as.numeric(get_par(fish_raw, "fisheries.rate.base"))/12
# selectivity <- get_par(fish_raw, "fisheries.selectivity.values")
# vec_size    <- get_par(fish_raw, "fisheries.selectivity.threshold")
nspecies    <- get_par(conf, "simulation.nspecies")
species     <- get_par(conf, "species.name")[1:nspecies]
namessp     <- stringr::str_remove(names(species), pattern = "species.name.sp")
# max_age     <- get_par(conf, "age.max")

# 
# paths_fish <- get_par(conf, "fishing.rate")[grep(x = get_par(conf, "fishing.rate"),pattern = ".csv")]
# list_nls <- list()
# k=0
# plots <- list()
# for (namessp_i in namessp){
#   p <- NULL
#   k=k+1
#   table_fish_rate <- read.csv(paste0("osmose-ns-v4.4.1_clim/", as.character(paths_fish[grep(x = names(paths_fish),pattern =  namessp_i)[1]])), row.names = c(1))
#   tab_nls         <- data.frame("size"=as.numeric(stringr::str_remove(colnames(table_fish_rate), pattern = "X")), "m_F"=as.numeric(table_fish_rate[1,]))
#   colnames(table_fish_rate) <- tab_nls$size
#   speci <- species[[as.numeric(namessp_i) +1]]
#   
#   plot(tab_nls$m_F~tab_nls$size)
#   
#   if (namessp_i %in% c("3", "15")){
#     selectiv_ft <- ifelse(tab_nls$size < max(tab_nls$size), 0, 1)
#     Fmax_ft = max(tab_nls$m_F)
#     list_nls[[namessp_i]] <- list("selectiv_ft"=selectiv_ft, "Fmax_ft"=Fmax_ft)
#     tab <- data.frame("size"=tab_nls$size, "estimate"=selectiv_ft*Fmax_ft, 
#                       "m_F"=tab_nls$m_F)
#     p <- ggplot(tab, aes(x = size, y = m_F)) +
#       geom_point(shape = 1) +
#       geom_line(aes(x = size, y = estimate), color = "red", size = 1) +
#       labs(x = "Size", y = "m_F", title =paste0(speci, "Binary fit")) +
#       theme_minimal(base_size = 14)
#   }
#   if (namessp_i == "2"){
#     selectiv_ft <- ifelse(tab_nls$size < mean(tab_nls$size), 0, 1)
#     Fmax_ft = max(tab_nls$m_F)
#     list_nls[[namessp_i]] <- list("selectiv_ft"=selectiv_ft, "Fmax_ft"=Fmax_ft)
#     tab <- data.frame("size"=tab_nls$size, "estimate"=selectiv_ft*Fmax, 
#                       "m_F"=tab_nls$m_F)
#     p <- ggplot(tab, aes(x = size, y = m_F)) +
#       geom_point(shape = 1) +
#       geom_line(aes(x = size, y = estimate), color = "red", size = 1) +
#       labs(x = "Size", y = "m_F", title = paste0(speci, "Binary fit")) +
#       theme_minimal(base_size = 14)
#   }
#   
#   if (namessp_i %in% c("0", "1")){
#     fit1 <- nlsLM(
#       m_F ~ Fmax /(1 + exp(-alpha*(size-halfsize))),
#       data = tab_nls,
#       start = list(
#         Fmax = max(tab_nls$m_F),
#         halfsize = 20,   # guess the "half-rise" size from your plot
#         alpha = 3
#       ),
#       lower = c(Fmax = 0, halfsize = 0, alpha = 0.1),
#       upper = c(Fmax = 1, halfsize = max(tab_nls$size), alpha = 50),
#       control = nls.lm.control(maxiter = 500)
#     )
#     # prediction grid
#     sz <- seq(min(tab_nls$size), max(tab_nls$size), length.out = 200)
#     pred_df <- data.frame(
#       size = sz,
#       m_F  = predict(fit1, newdata = data.frame(size = sz))
#     )
#     # plot with ggplot
#     p <- ggplot(tab_nls, aes(x = size, y = m_F)) +
#       geom_point(shape = 1) +
#       geom_line(data = pred_df, aes(x = size, y = m_F), color = "red", size = 1) +
#       labs(x = "Size", y = "m_F", title = paste0(speci,"Logistic fit")) +
#       theme_minimal(base_size = 14)
#   }
#   if (namessp_i %in% c("4", "5", "6")){
#     
#     tab_nls2 <- subset(tab_nls, size > 0) # to avoid having zeros which stops nls function loop
#     fit1 <- nlsLM(
#       m_F ~ Fmax*(1/((size/sizemax)^(sizemax/p)*exp((sizemax-size)/p))),
#       data = tab_nls2,
#       start = list(
#         Fmax=max(tab_nls$m_F),
#         sizemax = max(tab_nls$size),   # guess the "half-rise" size from your plot
#         p = 1
#       ),
#       lower = c(Fmax = 0, sizemax = 0, p = max(tab_nls$size, na.rm = T) / 5),
#       upper = c(Fmax = 1, sizemax = 2*max(tab_nls$size), p = max(tab_nls$size)),
#       control = nls.lm.control(maxiter = 500)
#     )
#     # prediction grid
#     sz <- seq(min(tab_nls$size), max(tab_nls$size), length.out = 200)
#     pred_df <- data.frame(
#       size = sz,
#       m_F  = predict(fit1, newdata = data.frame(size = sz))
#     )
#     # plot with ggplot
#     p <- ggplot(tab_nls, aes(x = size, y = m_F)) +
#       geom_point(shape = 1) +
#       geom_line(data = pred_df, aes(x = size, y = m_F), color = "red", size = 1) +
#       labs(x = "Size", y = "m_F", title = paste0(speci, "Gamma fit")) +
#       theme_minimal(base_size = 14)
#   }
#   plots[[k]] <-p
#   
#   # 
#   # list_nls[[namessp_i]] <- nls(m_F ~ Fmax * (1 + exp(alpha * (size - halfsize))), 
#   #                              data=tab_nls, 
#   #                              start = c(Fmax=max(tab_nls$m_F), alpha=10, halfsize=mean(tab_nls$size))) # logistic & gamma
#   # list_nls[[namessp_i]] <- nls(m_F ~ Fmax*(1/((size/sizemax)^(sizemax/p)*exp((sizemax-size)/p))), 
#   #                              data=tab_nls, 
#   #                              start = c(Fmax=max(tab_nls$m_F), alpha=0.1, halfsize=mean(tab_nls$size)/2, p=5, sizemax=max(tab_nls$size))) # logistic & gamma
#   # list_nls[[namessp_i]] <- nls(m_F ~ Fmax*(1/(1+exp(alpha*(size-halfsize))) + (size/sizemax)^(sizemax/p)*exp((sizemax-size)/p)), 
#   #                              data=tab_nls, 
#   #                              start = c(Fmax=max(tab_nls$m_F), alpha=0.1, halfsize=mean(tab_nls$size)/2, p=5, sizemax=max(tab_nls$size))) # logistic & gamma
#   # assign(paste0("t_fish", namessp_i), table_fish_rate)
# }
# 
# cowplot::plot_grid(plots[[1]], plots[[2]], plots[[3]], plots[[4]],plots[[5]],plots[[6]], plots[[15]], ncol=3)



# ### RUN SIMULATION FOR VARIABLE FMAX #####
# if(!dir.exists("osmose-ns-v4.4.1_clim/new_config")){dir.create("osmose-ns-v4.4.1_clim/new_config")}
# list_fish_rate <- list.files(path_raw, pattern = ".csv")
# for (sp in 1:nspecies){
#   sp_i   <- as.character(species[sp])
#   sp_id  <- namessp[sp]
#   cat("species =", sp_i, "\n")
#   # Fmax_i = Fmax[sp]
#   file <- grep(paste0(sp_id, ".csv"), list_fish_rate)
#   fish_rate_file <- read.csv(paste0(path_raw, "/", list_fish_rate[file[1]]))
#   colnames(fish_rate_file) <- stringr::str_remove(colnames(fish_rate_file), "X")
#   
#   for (Fmax_m in multipliers){
#     name_infile <- paste0("new_config/F_Fmax", Fmax_m, "_sp", sp_id, ".csv")
#     # what parameter value change
#     cat("Fmax_m =", Fmax_m, "\n")
#     # F_max_m_val  =  Fmax_i*Fmax_m #[grep(sp_id, names(Fmax))[1]]
#     new_fish_rate_mat <- fish_rate_file
#     new_fish_rate_mat[, -1] <- Fmax_m*new_fish_rate_mat[, -1]
#     # if (Fmax_m == 1) {new_fish_rate_mat[, -1] <- new_fish_rate_mat[, -1]/Fmax_i}
#     if(file.exists(paste0("osmose-ns-v4.4.1_clim/", name_infile))){file.remove(paste0("osmose-ns-v4.4.1_clim/", name_infile))}
#     write.csv(file=paste0("osmose-ns-v4.4.1_clim/", name_infile), new_fish_rate_mat, row.names = F)
#     
#     name_file <- paste0("new_config/conf_Fmax", Fmax_m, "_sp", sp_id, ".R")    
#     file_path <- paste0("osmose-ns-v4.4.1_clim/", name_file)
#     file.copy("osmose-ns-v4.4.1_clim/New_gpu_conf.R", 
#               file_path, overwrite = T)
#     lines <- readLines(file_path)
#     ID_to_change <- grep(paste0("mortality.fishing.rate.byDt.bySize.file.sp", sp_id), lines)[1]
#     lines[ID_to_change]
#     lines[ID_to_change] <- paste0("mortality.fishing.rate.byDt.bySize.file.sp", sp_id, " = ", name_infile)
#     writeLines(lines, file_path)
#     readLines(file_path, n = 1)
#     
#     # run simulation with updated parameter
#   #   to launch the simulation with the corresponding new parameters
#   #   parameters <- paste0("-Pmortality.fishing.rate.byDt.bySize.file.sp", sp_id, "=osmose-ns-v4.4.1_clim/fishing_rate/tempo.csv")
#   #   simulation = paste0("output_Fmax_", Fmax_m)
#   #   outputDir  = file.path(configDir, "Simu_Fisheries/output", simulation)
#   #   parameters = paste0("-Pfisheries.rate.base.sp", sp_id,"=", F_max_val)
#   #   options = paste0("-Xmx3g -Xms1g")
#   #   run_osmose(input = configFile, 
#   #              output = paste0(outputDir, "/", sp_i), 
#   #              osmose = jarFile, 
#   #              version=version, 
#   #              options=options, 
#   #              parameters = parameters)
#   # }
#   }
# }
# 
# 
# ### RUN SIMULATION FOR VARIABLE FIRST AGE FISHED #####
# list_fish_rate <- list.files(path_raw, pattern = ".csv")
# for (sp in 1:nspecies){#1:nspecies){
#   sp_i   <- as.character(species[sp])
#   sp_id  <- namessp[sp]
#   cat("species =", sp_i, "\n")
#   # Fmax_i = Fmax[sp]
#   file <- grep(paste0(sp_id, ".csv"), list_fish_rate)
#   fish_rate_file <- read.csv(paste0(path_raw, "/", list_fish_rate[file[1]]))
#   colnames(fish_rate_file) <- stringr::str_remove(colnames(fish_rate_file), "X")
#   for (Fstart_m in multipliers){
#     name_infile <- paste0("new_config/F_Fstart", Fstart_m, "_sp", sp_id, ".csv")
#     # what parameter value change
#     cat("Fstart_m =", Fstart_m, "\n")
#     # F_max_m_val  =  Fmax_i*Fmax_m #[grep(sp_id, names(Fmax))[1]]
#     
#     # increase or decrease (translation) of the first age with fishing mortality
#     ID_no_zero     <- which(fish_rate_file[1,] != 0)
#     first_length   <- (ID_no_zero[1]/2 - 1)
#     first_length_m <- first_length*Fstart_m
#     dec <- first_length_m - floor(first_length_m)
#     if (dec < 0.5 & dec > 0){first_length_m = first_length_m + 0.25}
#     if (dec < 1 & dec > 0.5){first_length_m = first_length_m + 0.25}
#     names_1             <- seq(0,first_length_m-0.5, by=0.5)
#     names_2             <- seq(first_length_m, as.numeric(rev(colnames(fish_rate_file[, -1]))[1]), by = 0.5)
#     df_1                <- matrix(rep(0, nrow(fish_rate_file)*length(names_1)), ncol=length(names_1), nrow=nrow(fish_rate_file))
#     df_2                <- as.matrix(fish_rate_file[, ID_no_zero])
#     missing_dimensions  <- length(names_2)-ncol(df_2)
#     if(missing_dimensions>0) {df_2 <- cbind(df_2, matrix(rep(df_2[,ncol(df_2)], missing_dimensions), ncol=missing_dimensions))}
#     if(missing_dimensions<0) {df_2 <- df_2[, -c((ncol(df_2)+missing_dimensions+1) : ncol(df_2))]}
#     
#     new_fish_rate_mat           <- cbind(df_1, df_2)
#     colnames(new_fish_rate_mat) <- c(names_1, names_2)
# 
#     
#     if(file.exists(paste0("osmose-ns-v4.4.1_clim/", name_infile))){file.remove(paste0("osmose-ns-v4.4.1_clim/", name_infile))}
#     write.csv(file=paste0("osmose-ns-v4.4.1_clim/", name_infile), new_fish_rate_mat, row.names = F)
#     
#     # 2 options : save the updated config in a separate repo OR launch simu with the new F mortality
#     name_file <- paste0("new_config/conf_Fstart", Fstart_m, "_sp", sp_id, ".R")    
#     file_path <- paste0("osmose-ns-v4.4.1_clim/", name_file)
#     file.copy("osmose-ns-v4.4.1_clim/New_gpu_conf.R", 
#               file_path, overwrite = T)
#     lines <- readLines(file_path)
#     ID_to_change <- grep(paste0("mortality.fishing.rate.byDt.bySize.file.sp", sp_id), lines)[1]
#     lines[ID_to_change]
#     lines[ID_to_change] <- paste0("mortality.fishing.rate.byDt.bySize.file.sp", sp_id, " = ", name_infile)
#     writeLines(lines, file_path)
#     readLines(file_path, n = 1)
#     
#     # run simulation with updated parameter
#     # simulation = paste0("output_Fstart_", Fstart_m)
#     # outputDir  = file.path(configDir, "Simu_Fisheries/output", simulation)
#     # parameters <- paste0("-Pmortality.fishing.rate.byDt.bySize.file.sp", sp_id, "=osmose-ns-v4.4.1_clim/fishing_rate/tempo2.csv")
#     # parameters = paste0("-Pfisheries.rate.base.sp", sp_id,"=", F_max_val)
#     # options = paste0("-Xmx3g -Xms1g")
#     # run_osmose(input = configFile,
#     #            output = paste0(outputDir, "/", sp_i),
#     #            osmose = jarFile,
#     #            version=version,
#     #            options=options,
#     #            parameters = parameters)
#   }
# }




# crossed changes of Fmax and Fstart

### RUN SIMULATION FOR VARIABLE FMAX #####
if(!dir.exists("osmose-ns-v4.4.1_clim/new_config")){dir.create("osmose-ns-v4.4.1_clim/new_config")}
list_fish_rate <- list.files(path_raw, pattern = ".csv")
list_fish_rate <- list_fish_rate[order(as.numeric(unlist(stringr::str_split(unlist(stringr::str_split(list_fish_rate,
                                 pattern = "\\."))[seq(7, length(unlist(stringr::str_split(list_fish_rate, 
                                 pattern = "\\."))), by=8)],
                                 pattern="sp"))[seq(2,16*2, by=2)]))]
for (sp in 1:nspecies){
  sp_i   <- as.character(species[sp])
  sp_id  <- namessp[sp]
  cat("species =", sp_i, "\n")
  # Fmax_i = Fmax[sp]
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
    
    
    for (Fstart_m in multipliers){
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
    
      if(file.exists(paste0("osmose-ns-v4.4.1_clim/new_config/", name_infile))){file.remove(paste0("osmose-ns-v4.4.1_clim/new_config/", name_infile))}
      write.csv(file=paste0("osmose-ns-v4.4.1_clim/new_config/", name_infile), new_fish_rate_matFSTART, row.names = F)
      
      name_file <- paste0("new_config/conf_Fmax", Fmax_m, "_Fstart", Fstart_m, "_sp", sp_id, ".R")    
      file_path <- paste0("osmose-ns-v4.4.1_clim/", name_file)
      file.copy("osmose-ns-v4.4.1_clim/New_gpu_conf.R", 
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



# 
# 
# ### RUN SIMULATION FOR VARIABLE SELECTIVITES #####
# selectivity_m <- c(0.5,1,2)
# for (sp in 1:nspecies){
#   
#   sp_i   <- as.character(species[sp])
#   sp_id  <- namessp[sp]
#   cat("species =", sp_i, "\n")
#   
#   for (selectivity_m in Fmax_m_v){
#     
#     # what parameter value change
#     cat("Fmax_m =", Fmax_m, "\n")
#     F_max_val  = ifelse(Fmax_m>1, 1, Fmax*Fmax_m) #[grep(sp_id, names(Fmax))[1]]
#     simulation = paste0("output_Fmax_", Fmax_m)
#     outputDir  = file.path(configDir, "output", simulation)
#     
#     # run simulation with updated parameter
#     parameters = paste0("-Pfisheries.rate.base.sp", sp_id,"=", )
#     options = paste0("-Xmx3g -Xms1g")
#     run_osmose(input = configFile, 
#                output = outputDir, 
#                osmose = jarFile, 
#                version=version, 
#                options=options, 
#                parameters = parameters)
#     
#   }
# }
