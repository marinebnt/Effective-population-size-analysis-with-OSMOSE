
main="New_gpu_conf.R"
source(file.path("osmose-ns-v4.4.1_clim/R/functions/source.R"))
repo        = "RESTART"
multip_name <- c("y")
multipliers = c("0.5", "1", "1.5")
multipliers = 1
multipliers = get_par(conf, "simulation.time.nyear")

# function to prepare environment : launch simulation if repo missing
function_INIT1(MULTIP=multipliers, REPO=repo)

# prepare environment for this particular plot
rep=1
incr_multiplier <- 0
vecdim= c("biomass", "school", "observedHtz", "yield", "meanSize")
# multipliers = "output_etaequalsone_bioen_neutral"

for (trait_dim in vecdim){
  for (m in multipliers){
    
    cat("load simulation output")
    simulation = paste0(repo, "/output_", multip_name, m)
    Ne_version = paste0("plot_biomass", multip_name, m)
    cat("\n\t---> folder ", simulation, "\n", trait_dim, "\n")
    ns1 = NULL
    outputDir_Simu  = file.path(configDir, "output", simulation)
    outputDir_Ne =paste0(configDir, "/output/", "Ne","/PLOT_REPLICATES/")
    if (any(list.files(outputDir_Simu) %in% species)){
      outputDir_Simu    <- paste0(outputDir_Simu, "/", sp_focus)
      outputDir_Ne <- paste0(outputDir_Ne, "Fmax_on_", sp_focus, "/")
    }
    if(!dir.exists(outputDir_Ne)) {dir.create(outputDir_Ne)}
    ns1  = read_osmose(path = outputDir_Simu)
    conf = read_osmose(path = paste0(outputDir_Simu, "/ns-configuration.osm"))
    nrep = conf$simulation.nsimulation
    data <- ns1[[trait_dim]]
    if(is.null(data)){
      if(!(grepl("Htz", trait_dim))){
        data     <- list()
        files    <- list.files(outputDir_Simu, full.names = T, recursive = T)
        files    <- files[c(grep(".csv", files), grep(".nc", files))]
        files_rep0    <- files[c(grep("Simu0", files))]
        ID_files <- grep(trait_dim, files_rep0)
        dim_data_ID_files <- unlist(sapply(X = files_rep0[ID_files], 
              FUN=function(x) {dim(read.csv(x, skip = c(1)))}))
        ID_file_selected  <- which.min(dim_data_ID_files[1,]) # chose the file the smallest : exclude the data by age or by size ...
        if (nrep>1){
          str_a <- stringr::str_remove(pattern = "0.csv", string = files_rep0[ID_files[ID_file_selected]])
          str_b <- stringr::str_remove(pattern = "0.nc", string = files_rep0[ID_files[ID_file_selected]])
          ID_file_selected_Simu <- grep(c(str_a, str_b)[which.min(c(nchar(str_a), nchar(str_b)))], files)
        } else {ID_file_selected_Simu <- which(files == files_rep0[ID_files[ID_file_selected]])}
        for (ID_i in 1:length(ID_file_selected_Simu)){
          ID <- ID_file_selected_Simu[ID_i]
          if(grepl("ByAge", files[ID])){
            data_A <- as.data.frame(read.csv(files[ID], skip=c(1)))
            data_A <- data_A[, -which(colnames(data_A)=="Age")]
            avg_by_species_age <- data_A %>%
              group_by(Time) %>%
              summarise(across(where(is.numeric), \(x) mean(x, na.rm = TRUE)), .groups = "drop")
            data[[ID_i]] <- avg_by_species_age
          }else {data[[ID_i]] <- as.data.frame(read.csv(files[ID], skip=c(1)))}
        }
      } else {
        data_r <- NULL
        data_r_list <- array(data=NA, dim=c(years, length(species)+1, nrep))
        for (r in 0:c(nrep-1)){
          for (sp in species){
            data <- ncdf4::nc_open(paste0(outputDir_Simu, "/Genetic/ns_observedHtz-", sp, "Simu", r,".nc"))
            data_nc <- ncdf4::ncvar_get(data, "observed_heterozygosity")
  cat("\n", (paste0(outputDir_Simu, "/Genetic/ns_observedHtz-", sp, "Simu", r,".nc")))
  cat("\nDIM HTZ DATA ", dim(data_nc))
            trait1_data <- data_nc[, ]
            ncdf4::nc_close(data)
            mean_htz <- colMeans(trait1_data)
            time <- 0:c(dim(data_nc)[2]-1)
            data_htz <- data.frame(Time=time, HTZ=mean_htz)
            if (sp == species[1]) {
              data_r <- data_htz
            }else {data_r <- cbind(data_r, data_htz[,2])}
          }
          colnames(data_r) <- c("time", species)
          data_r_list[,,r+1] = as.matrix(data_r)
        }
        dimnames(data_r_list) <- list(time, c("time", species), paste0("rep", 0:c(nrep-1)))
        data = data_r_list
      }
    }
    
    p <- list()
    for (i in 1:nrep){
      data_t <- data[[i]]
      if(is.null(dim(data_t))){data_t <- data[,,i]}
      dataset_plot <- as.data.frame(data_t)
      dataset_plot$time <- as.numeric(rownames(dataset_plot))
      if (any(colnames(dataset_plot) == "Time")) {
        dataset_plot$time <- dataset_plot$Time
        dataset_plot <- dataset_plot[,-which(colnames(dataset_plot) == c("Time"))]}
      dataset_plot <- dataset_plot %>% 
        pivot_longer(!time, values_to = trait_dim, names_to = "species")
      # dataset_plot[,trait_dim] <- as.numeric(dataset_plot[,trait_dim])
      if (trait_dim == "biomass"){dataset_plot[,trait_dim] <- dataset_plot[,trait_dim]/10^3}
      dataset_plot <- dataset_plot[which(dataset_plot$time>100),]
      p[[i]] <- ggplot(dataset_plot, aes_string(x='time', y=trait_dim))+ 
          geom_line()+
          facet_wrap(~ species, ncol=4, nrow=4)+
          theme_minimal()+
          labs(title=paste0("Replicate ", i, "(x 10^3 if biomass)"))
    }
    
    q <- list()
    dataset_plot <- NULL
    for (sp in species){
      for (i in 1:nrep){
        data_t <- data[[i]]
        if(is.null(dim(data_t))){data_t <- data[,,i]}
        if (i==1) {dataset_plot <- data_t[,sp]
        }else {dataset_plot <- as.data.frame(cbind(dataset_plot, data_t[,sp]))}
      }
      if(!is.null(dim(dataset_plot))){colnames(dataset_plot) <- paste0("rep", paste0("rep", 0:c(nrep-1)))
      }else{dataset_plot <- data.frame("rep1"=dataset_plot)}
      if (any(colnames(data_t) == "Time")) {
        dataset_plot$time <- data_t$Time}else {dataset_plot$time <- as.numeric(rownames(data))}
      dataset_plot <- dataset_plot %>% 
        pivot_longer(!time, values_to = trait_dim, names_to = "replicate")
      # dataset_plot[,trait_dim] <- as.numeric(dataset_plot[,trait_dim])
      if (trait_dim == "biomass"){dataset_plot[,trait_dim] <- dataset_plot[,trait_dim]/10^3}
      dataset_plot <- dataset_plot[which(dataset_plot$time>100),]
      q[[sp]] <- ggplot(dataset_plot, aes_string(x='time', y=trait_dim))+ 
        geom_line(aes(colour = replicate))+
        theme_minimal()+
        labs(title=paste0(sp," ", trait_dim, " for ", nrep, " replicates"))
    }
    
    pdf(file = paste0(outputDir_Ne, trait_dim, "_BY_REP_ALLrep.pdf"))
    for (i in 1:ns1$model$simus){
      print(p[[i]])
    }
    dev.off()
    
    pdf(file = paste0(outputDir_Ne, trait_dim, "_BY_SP_ALLrep.pdf"))
    for (sp in species){
      print(q[[sp]])
    }
    dev.off()
  }
}
