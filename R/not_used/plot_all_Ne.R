
setwd("X:/osmose-ns")

main = "New_gpu_conf.csv"
source(file.path("osmose-ns-v4.4.1_clim/R/functions/source.R"))
repo        = "Simu_LONG"
multip_name = "y"
multipliers = c(1,0.5,0.75,1.25,1.5)


path_in <- c("Y:/outputs-meso-marine")
path_in <- c("/home1/datawork/outputs-meso-marine")
list_d <- list.dirs(path = path_in, recursive = F, full.names = T)


for (sp in seq(0, 15)){
  
  # plots lists
  plot_Neth <- c()
  plot_Nete <- c()
  plot_Nb   <- c()
  plot_htz  <- c()
  
  # open one dir at a time
  sp_dir <- list_d[grep(pattern = paste0("sp", sp), list_d)]
  for (Fmax in multipliers){
    Fmax_dir <- sp_dir[grep(pattern = paste0("Fmax", Fmax, "_"), sp_dir)]
    for (Fstart in multipliers){
      dir <- Fmax_dir[grep(pattern = paste0("Fstart", Fstart, "_"), Fmax_dir)]
        
      # initialisation
      simulation = dir
      cat("\n\t---> folder Fmax", Fmax, " Fstart", Fstart, " sp", sp, "\n")
      
      ns1 = NULL
      outputDir_Simu  = simulation
      outputDir_Ne    = file.path(configDir, "output", "PLOTS", stringr::str_split(dir, pattern = "/")[[1]][3])
      dir.create(outputDir_Ne, recursive = T)
      dim_vec  <- c("abundanceByAge", "eggNumber", "mortalityByAge", "schoolNumber")
      ns1 = read_osmose(path = outputDir_Simu)
      
      INIT2 <- function_INIT2(outputDir_Simu, outputDir_Ne, dim_vec, species, nb_rep, OUTPUT_data=ns1, conf, seeding_y, REPO=repo, TEMPORAL=FALSE)
      data_Ne_raw  <- INIT2[[1]]
      data_Ne      <- INIT2[[2]]
      
      Neth <- data_Ne[[(as.numeric(sp)+1)]][[]]
      Yea  <- 
      
      # produce data 
      if (Fmax == 1 & Fstart == 1){
        # data_f for plots
        Neth_df <- data.frame(reference=rep("reference", length(Neth)), 
                              year = Yea, 
                              Ne_th = Neth, 
                              Fish_press = rep(paste0("Fmax", Fmax, " Fstart", Fstart), length(Neth)))
        Nete_df <- NULL
        Nb_df   <- NULL
        htz_df  <- NULL
      }
      
      Neth_df <- rbind(data.frame(reference=rep("reference", length(Neth)), 
                            year = Yea, 
                            Ne_th = Neth, 
                            Fish_press = rep(paste0("Fmax", Fmax, " Fstart", Fstart), length(Neth))))
      Nete_df <- NULL
      Nb_df   <- NULL
      htz_df  <- NULL
     
        
    }#Fstart
  } #Fmax
  
  # i want this subset one to appear on every plot
  reference <- subset(Neth_df, Fish_press == paste0("Fmax", 1, " Fstart", 1))[, c("Ne_th", "year")]
  plot_Neth <- ggplot(data = Neth_df, aes(x = year, y = Ne_th, group = Fish_press)) +
    # 1) Layer: each panel gets its own group (via facet_wrap)
    geom_line() +
    
    # 2) Layer: reference group A, drawn in all panels
    geom_line(data = reference,
              aes(x = year, y = Ne_th, group = 1),
              linetype = "dashed", linewidth = 0.7) +
    
    facet_wrap(~ Fish_press)
  
  
  
} #sp
