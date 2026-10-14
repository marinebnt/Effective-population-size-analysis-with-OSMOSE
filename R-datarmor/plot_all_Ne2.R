# setwd("X:/osmose-ns")

plot_all_Ne2 <- function(SP){
  
  # INITIALISATION
  repo        = "Simu_LONG"
  multip_name = "y"
  multipliers = c(1,0.5,0.75,1.25,1.5)
  removing_senescence_param = F #***************
  configDir=getwd()
  
  ideal_order_factor <- paste0("Fmax=", rep(c(0.5,0.75,1,1.25,1.5), 5),
                               " Fstart=", c(rep(0.5, 5),rep(0.75, 5),rep(1, 5),rep(1.25, 5),rep(1.5, 5)))
  
  
  path_in <- c("output/DATASETS/old")
  path_in <- c("C:/Users/mbeneat/Documents/ZAutres/simu_test")
  path_in <- c("/home1/datawork/mbeneat/SIMU_Ne_on_NT_with_F_FIN")
  list_d <- list.dirs(path = path_in, recursive = F, full.names = T)[1]
  # for the test : has to be erased when real launch
  # nb_sp = 7
  # species = species[1:7]
  # sp_names = sp_names[1:7]
  # species_i = species_i[1:7]
  # years = 120
  # ny       <- 120-get_par(conf, "output.start.years")
  # time_vec <- seq(0,ny, by=1/dt)[-1]
  # nb_rep = 1
  # ncohort = years
  
  pathout <- paste0(configDir, "/output/PLOTS/All_clean_genetic_drift_indicators")
  pathout <- paste0("/home1/scratch/mbeneat/output/PLOTS/Ne_NT_calculated_variables")
  
  if (removing_senescence_param)  pathout <- paste0(pathout, "_noSenescence")

plot_Neth_Nalpha    <- list()
plot_Neth_Nalpha_NT <- list()
plot_Neth_Nalpha_Nalpha <- list()
plot_Neth_N1        <- list()
plot_Neth_N1_NT     <- list()
plot_Neth_N1_Nalpha <- list()
plot_Nete <- list()
plot_Nete_NT <- list()
plot_Nb   <- list()
plot_Nb_NT   <- list()
plot_htz  <- list()
F_vec     <- list()


spec_names <- species


for (sp in SP){
  # progressr::handlers(global = TRUE)
  sp_i <- which(spec_names == sp)-1
  species = spec_names[c(c(sp_i+1), ifelse(sp_i<14, sp_i+2, 1))]
  # open one dir at a time
  sp_dir <- list_d[grep(pattern = paste0("sp", sp_i), list_d)]
  
  if (length(sp_dir)<1) next
  
  Neth_N_alpha_df = list()
  Neth_N_1_df = list()
  Netp_df = list()
  Nb_df = list()
  htz_df = list()
  
  plot_Neth_N_alpha_NT_sp     <- list()
  plot_Neth_N_1_NT_sp         <- list()  
  plot_Neth_N_alpha_Nalpha_sp <- list()
  plot_Neth_N_1_Nalpha_sp     <- list()
  plot_Nete_NT_sp             <- list()
  plot_Nb_NT_sp               <- list()
  plot_Neth_Nalpha_sp <- list()
  plot_Neth_N1_sp     <- list()
  plot_Nete_sp <- list()
  plot_Nb_sp   <- list()
  plot_htz_sp  <- list()
  Fmax_1 = -99
  Fstart_1 = -99
  
#  folder_metadata <- tibble(folder = folders) %>%
#    separate(folder, into = c("Fmax", "Fstart"), 
#             regex = "Fmax([^_]+)_Fstart([^_]+)", remove = FALSE) %>%
#    select(folder, Fmax, Fstart)
  folder_metadata <- tibble(folder = sp_dir)
  folder_metadata$Fmax <- sapply(sp_dir, FUN = function(x) as.numeric(stringr::str_split(stringr::str_split(x, pattern="Fmax")[[1]][2], pattern="_Fstart")[[1]][1]))
  folder_metadata$Fstart <- sapply(sp_dir, FUN = function(x) as.numeric(stringr::str_split(stringr::str_split(x, pattern="Fstart")[[1]][2], pattern="_sp")[[1]][1]))


  for (i in seq_len(nrow(folder_metadata))) {
    f_i <- folder_metadata$folder[i]
    Fmax <- folder_metadata$Fmax[i]
    Fstart <- folder_metadata$Fstart[i]
    
    if (Fmax_1==-99 && Fstart_1 ==-99){
      Fmax_1 = Fmax
      Fstart_1 = Fstart
      F_vec[[sp]] <- c(Fmax_1, Fstart_1)
    }
    
    # initialisation
    simulation = f_i
    cat("\n\t---> folder Fmax", Fmax, " Fstart", Fstart, " sp", sp, " dir ",f_i, "\n")
    
    ns1 = NULL
    outputDir_Simu  = f_i 
    outputDir_Ne    = file.path(configDir, "output", "PLOTS", paste0("Fmax", Fmax, "Fstart", Fstart, "sp", sp_i))
    if(!dir.exists(outputDir_Ne)) {dir.create(outputDir_Ne, recursive = T)}
    dim_vec  <- c( "Mpred", "Mstar", "Mnat", "F", "Mfor","schoolNumber", "abundanceByAge", "eggNumber", "mortalityByAge")
    ns1 = read_osmose(path = paste0(outputDir_Simu, "/"))
      
      
    INIT2 <- function_INIT2(outputDir_Simu=outputDir_Simu, outputDir_Ne, dim_vec, species, nb_rep,
     OUTPUT_data=ns1, conf, seeding_y, REPO=dir,
     TEMPORAL=F, HTZ=TRUE, removing_senescence_param=removing_senescence_param)
    pathoutputINIT <-paste0(pathout, "/INIT/Fmax", Fmax, "_Fstart", Fstart, "_sp", sp)
    if (!dir.exists(pathoutputINIT)){dir.create(pathoutputINIT, recursive = T)}
    if (file.exists(paste0(pathoutputINIT, "/INIT2.RData"))) file.remove(paste0(pathoutputINIT, "/INIT2.RData"))
    save(INIT2, file=paste0(pathoutputINIT, "/INIT2.RData"))
    #load(file=paste0(pathoutputINIT, "/INIT2.RData"))
      
    data_Nete      <- INIT2[["data_Ne"]]
#    data_Nb        <- INIT2[["data_Nb"]]
#    data_Netp      <- INIT2[["tempoNe"]]
    data_htz       <- INIT2[["htz_list"]]
rm(INIT2) 
    for (sp_plotted in species){
      
      for (rep in 1:nb_rep){
        #Ne theoretical and associated variables
        first_maturation_age <- data_Nete$new_first_age[[sp_plotted]][[rep]]
        Sum_abund_by_age     <- colSums(data_Nete[[sp_plotted]][["abundanceByAge"]][[paste0("rep_", rep)]])
        abund_mature         <- as.data.frame(data_Nete[[sp_plotted]][["Numerator_N_alpha"]])
        Neth_Nalpha <- data_Nete[[sp_plotted]][["Ne_N_alpha"]][,paste0("rep_", rep)]
        # Neth_N1     <- data_Nete[[sp_plotted]][["Ne_N_1"]][,paste0("rep_", rep)]
        
        #Nb and associated variables
 #       NT                   <- data_Nb[[sp_plotted]][["NT"]][,paste0("rep_", rep)]
 #       Nb                   <- data_Nb[[sp_plotted]][["Nb"]][,paste0("rep_", rep)]
        
        #Ne temporal 
 #       Netp <- data_Netp[[sp_plotted]][["Ne_temporal"]][,paste0("rep_", rep)]
        
        # heterozygosity
        Htz  <- data_htz[[sp_plotted]][,paste0("rep_", rep)]
        
        # else
        dimensions_outputs = years-80
        first_y_Ne_tempo <- max(seeding_y, c(years-dimensions_outputs[1]+1))
        Yea  <- first_y_Ne_tempo:(years)
        Yea_tp  <- first_y_Ne_tempo:(years-20) #20 :2* generation length max : gap between t1 and t2 to estimate Netp
        
        # keep only data for the corresponding years
        Neth_Nalpha <- Neth_Nalpha[which(rownames(Neth_Nalpha) %in% paste0("cohort_", Yea_tp))] 
        # Neth_N1 <- Neth_N1[which(rownames(Neth_N1) %in% paste0("cohort_", Yea_tp))] 
  #      Nb   <- Nb[which(rownames(data_Nb[[sp_plotted]][["Nb"]]) %in% Yea_tp)]   
        Htz  <- Htz[which(rownames(data_htz[[sp_plotted]]) %in% Yea_tp)]   
  #      Netp <- Netp[which(rownames(data_Netp[[sp_plotted]][["Ne_temporal"]]) %in% (Yea_tp))]   
        abund_mature <- abund_mature[which(rownames(abund_mature) %in% paste0("cohort_", Yea_tp)), paste0("rep_", rep)] 
        Sum_abund_by_age     <- Sum_abund_by_age[which(names(Sum_abund_by_age) %in% paste0("cohort_", Yea_tp))] 
  #      NT_tp                   <- NT[which(rownames(data_Nb[[sp_plotted]][["NT"]]) %in% paste0(Yea_tp))] 
  #      NT                   <- NT[which(rownames(data_Nb[[sp_plotted]][["NT"]]) %in% paste0(Yea_tp))] 
        
        # produce data 
        if (Fmax == Fmax_1 && Fstart == Fstart_1 && rep==1){
            # data_f for plots
            Neth_N_alpha_df[[sp_plotted]] <- data.frame(years = Yea_tp, 
                                          Ne_th = (Neth_Nalpha), 
                                          Ne_th_NT = (Neth_Nalpha)/Sum_abund_by_age, 
                                          Ne_th_Nalpha = (Neth_Nalpha)/abund_mature, 
                                          Fish_press = rep(paste0("Fmax=", Fmax, " Fstart=", Fstart), length(Neth_Nalpha)),
                                          replicate = rep-1)              
            # Neth_N_1_df[[sp_plotted]] <- data.frame(years = Yea_tp, 
            #                               Ne_th = (Neth_N1), 
            #                               Ne_th_NT = (Neth_N1)/Sum_abund_by_age, 
            #                               Ne_th_Nalpha = (Neth_N1)/abund_mature, 
            #                               Fish_press = rep(paste0("Fmax=", Fmax, " Fstart=", Fstart), length(Neth_N1)),
            #                               replicate = rep-1)
 #           Netp_df[[sp_plotted]] <- data.frame(years = Yea_tp, 
 #                                         Ne_tp = (Netp), 
 #                                         Netp_NT = (Netp)/NT_tp, 
 #                                         Fish_press = rep(paste0("Fmax=", Fmax, " Fstart=", Fstart), length(Netp)),
 #                                         replicate = rep-1)
 #           Nb_df[[sp_plotted]]   <- data.frame(data.frame(years = Yea_tp, 
 #                                         Nb = Nb, 
 #                                         Nb_NT = (Nb)/NT, 
 #                                         Fish_press = rep(paste0("Fmax=", Fmax, " Fstart=", Fstart), length(Nb))),
 #                                         replicate = rep-1)
            htz_df[[sp_plotted]]  <- data.frame(data.frame(years = Yea_tp, 
                                          Htz = Htz, 
                                          Fish_press = rep(paste0("Fmax=", Fmax, " Fstart=", Fstart), length(Yea_tp))),
                                          replicate = rep-1)

        } else {
          Neth_N_alpha_df[[sp_plotted]] <- rbind(Neth_N_alpha_df[[sp_plotted]],
                           data.frame(years = Yea_tp,
                                        Ne_th = (Neth_Nalpha), 
                                        Ne_th_NT = (Neth_Nalpha)/Sum_abund_by_age, 
                                        Ne_th_Nalpha = (Neth_Nalpha)/abund_mature, 
                                        Fish_press = rep(paste0("Fmax=", Fmax, " Fstart=", Fstart), length(Neth_Nalpha)),
                                        replicate = rep-1))            
          # Neth_N_1_df[[sp_plotted]] <- rbind(Neth_N_1_df[[sp_plotted]],
          #                  data.frame(years = Yea_tp,
          #                               Ne_th = (Neth_N1), 
          #                               Ne_th_NT = (Neth_N1)/Sum_abund_by_age, 
          #                               Ne_th_Nalpha = (Neth_N1)/abund_mature, 
          #                               Fish_press = rep(paste0("Fmax=", Fmax, " Fstart=", Fstart), length(Neth_N1)),
          #                               replicate = rep-1))
  #        Netp_df[[sp_plotted]] <- rbind(Netp_df[[sp_plotted]],
  #                         data.frame(years = Yea_tp, 
  #                                      Ne_tp = (Netp), 
  #                                      Netp_NT = (Netp)/NT_tp, 
  #                                      Fish_press = rep(paste0("Fmax=", Fmax, " Fstart=", Fstart), length(Netp)),
  #                                      replicate = rep-1))
  #        Nb_df[[sp_plotted]]   <- rbind(Nb_df[[sp_plotted]], 
  #                         data.frame(years = Yea_tp, 
  #                                      Nb = Nb, 
  #                                      Nb_NT = (Nb)/NT, 
  #                                      Fish_press = rep(paste0("Fmax=", Fmax, " Fstart=", Fstart), length(Nb)),
  #                                      replicate = rep-1))
          htz_df[[sp_plotted]]  <-  rbind(htz_df[[sp_plotted]], data.frame(years = Yea_tp, 
                                        Htz = (Htz), 
                                        Fish_press = rep(paste0("Fmax=", Fmax, " Fstart=", Fstart), length(Yea_tp)),
                                        replicate = rep-1))
           
        }
      }#rep_plot
    }#sp_plotted
    # p(paste0("\n\t--->done: folder Fmax", Fmax, " Fstart", Fstart, " sp", sp, "\n"))
  }#Folders Fmax&Fstart
  
  Neth_N_alpha_df_mean = NULL
  Neth_N_1_df_mean = NULL
  Neth_N_alpha_NT_df_mean = NULL
  Neth_N_1_NT_df_mean = NULL
  Neth_N_alpha_Nalpha_df_mean = NULL
  Neth_N_1_Nalpha_df_mean = NULL
  Netp_df_mean = NULL
  Netp_NT_df_mean = NULL
  Nb_df_mean = NULL
  Nb_NT_df_mean = NULL
  htz_df_mean = NULL
  
  for (sp_plotted in species){
    # outputs Ne theoretical
    Neth_N_alpha_df_mean <- Neth_N_alpha_df[[sp_plotted]] %>%
      group_by(Fish_press, years) %>%
      summarise(Ne_th_mean = mean(Ne_th, na.rm = TRUE),
                Ne_th_sd   = sd(Ne_th, na.rm = TRUE))        
    # Neth_N_1_df_mean <- Neth_N_1_df[[sp_plotted]] %>%
    #   group_by(Fish_press, years) %>%
    #   summarise(Ne_th_mean = mean(Ne_th, na.rm = TRUE),
    #             Ne_th_sd   = sd(Ne_th, na.rm = TRUE))
    Neth_N_alpha_NT_df_mean <- Neth_N_alpha_df[[sp_plotted]] %>%
      group_by(Fish_press, years) %>%
      summarise(Ne_th_mean = mean(Ne_th_NT, na.rm = TRUE),
                Ne_th_sd   = sd(Ne_th_NT, na.rm = TRUE))        
    # Neth_N_1_NT_df_mean <- Neth_N_1_df[[sp_plotted]] %>%
    #   group_by(Fish_press, years) %>%
    #   summarise(Ne_th_mean = mean(Ne_th_NT, na.rm = TRUE),
    #             Ne_th_sd   = sd(Ne_th_NT, na.rm = TRUE))
    Neth_N_alpha_Nalpha_df_mean <- Neth_N_alpha_df[[sp_plotted]] %>%
      group_by(Fish_press, years) %>%
      summarise(Ne_th_mean = mean(Ne_th_Nalpha, na.rm = TRUE),
                Ne_th_sd   = sd(Ne_th_Nalpha, na.rm = TRUE))        
    # Neth_N_1_Nalpha_df_mean <- Neth_N_1_df[[sp_plotted]] %>%
    #   group_by(Fish_press, years) %>%
    #   summarise(Ne_th_mean = mean(Ne_th_Nalpha, na.rm = TRUE),
    #             Ne_th_sd   = sd(Ne_th_Nalpha, na.rm = TRUE))
    
    #outputs Ne tempo
 #   Netp_df_mean <- Netp_df[[sp_plotted]] %>%
 #     group_by(Fish_press, years) %>%
 #     summarise(Ne_tp_mean = mean(Ne_tp, na.rm = TRUE),
 #               Ne_tp_sd   = sd(Ne_tp, na.rm = TRUE))
 #   Netp_NT_df_mean <- Netp_df[[sp_plotted]] %>%
 #     group_by(Fish_press, years) %>%
 #     summarise(Ne_tp_mean = mean(Netp_NT, na.rm = TRUE),
 #               Ne_tp_sd   = sd(Netp_NT, na.rm = TRUE))
    
    # outputs Nb
 #   Nb_df_mean <- Nb_df[[sp_plotted]] %>%
 #     group_by(Fish_press, years) %>%
 #     summarise(Nb_mean = mean(Nb, na.rm = TRUE),
 #               Nb_sd   = sd(Nb, na.rm = TRUE))        
 #   Nb_NT_df_mean <- Nb_df[[sp_plotted]] %>%
 #     group_by(Fish_press, years) %>%
 #     summarise(Nb_mean = mean(Nb_NT, na.rm = TRUE),
 #               Nb_sd   = sd(Nb_NT, na.rm = TRUE))
    
    # outputs heterozygosity
    htz_df_mean <- htz_df[[sp_plotted]] %>%
      group_by(Fish_press, years) %>%
      summarise(htz_mean = mean(Htz, na.rm = TRUE),
                htz_sd   = sd(Htz, na.rm = TRUE))
    
    plot_Neth_Nalpha_sp[[sp_plotted]] <- Neth_N_alpha_df_mean
    # plot_Neth_N1_sp[[sp_plotted]]     <- Neth_N_1_df_mean        
    plot_Neth_N_alpha_NT_sp[[sp_plotted]] <- Neth_N_alpha_NT_df_mean
    # plot_Neth_N_1_NT_sp[[sp_plotted]]     <- Neth_N_1_NT_df_mean        
    plot_Neth_N_alpha_Nalpha_sp[[sp_plotted]] <- Neth_N_alpha_Nalpha_df_mean
    # plot_Neth_N_1_Nalpha_sp[[sp_plotted]]     <- Neth_N_1_Nalpha_df_mean
  #  plot_Nete_sp[[sp_plotted]]    <- Netp_df_mean
  #  plot_Nete_NT_sp[[sp_plotted]] <- Netp_NT_df_mean
  #  plot_Nb_sp[[sp_plotted]]      <- Nb_df_mean
  #  plot_Nb_NT_sp[[sp_plotted]]   <- Nb_NT_df_mean
    plot_htz_sp[[sp_plotted]]     <- htz_df_mean        
    # plot_Neth_sp[[sp_plotted]] <- Neth_df_sd
    # plot_Nete_sp[[sp_plotted]] <- Netp_df_sd
    # plot_Nb_sp[[sp_plotted]]   <- Nb_df_sd
    # plot_htz_sp[[sp_plotted]]  <- htz_df_sd
  }
  
  pathoutputRDATA <-paste0(pathout, "/RDATA/Fmax", Fmax, "_Fstart", Fstart, "_sp", sp)
  if (!dir.exists(pathoutputRDATA)){dir.create(pathoutputRDATA, recursive = T)}
  
  save(plot_Neth_Nalpha_sp, file=paste0(pathoutputRDATA, "/plot_Neth_Nalpha_sp.RData"))
  # save(plot_Neth_N1_sp, file=paste0(pathoutputRDATA, "/plot_Neth_N1_sp.RData"))  
  save(plot_Neth_N_alpha_NT_sp, file=paste0(pathoutputRDATA, "/plot_Neth_N_alpha_NT_sp.RData"))
  # save(plot_Neth_N_1_NT_sp, file=paste0(pathoutputRDATA, "/plot_Neth_N_1_NT_sp.RData"))  
  save(plot_Neth_N_alpha_Nalpha_sp, file=paste0(pathoutputRDATA, "/plot_Neth_N_alpha_Nalpha_sp.RData"))
  # save(plot_Neth_N_1_Nalpha_sp, file=paste0(pathoutputRDATA, "/plot_Neth_N_1_Nalpha_sp.RData"))
 # save(plot_Nete_sp, file=paste0(pathoutputRDATA, "/plot_Nete_sp.RData"))
 # save(plot_Nete_NT_sp, file=paste0(pathoutputRDATA, "/plot_Nete_NT_sp.RData"))
 # save(plot_Nb_sp, file=paste0(pathoutputRDATA,"/plot_Nb_sp.RData"))
 # save(plot_Nb_NT_sp, file=paste0(pathoutputRDATA,"/plot_Nb_NT_sp.RData"))
  save(plot_htz_sp, file=paste0(pathoutputRDATA, "/plot_htz_sp.RData"))
  

  for (sp_plotted in species){  

    # reorder Fmax/Fstart factor for prettier plots
    plot_Neth_Nalpha_sp[[sp_plotted]]  <- as.data.frame(plot_Neth_Nalpha_sp[[sp_plotted]])
    # plot_Neth_N1_sp[[sp_plotted]]<- as.data.frame(plot_Neth_N1_sp[[sp_plotted]])
    plot_Neth_N_alpha_NT_sp[[sp_plotted]]<- as.data.frame(plot_Neth_N_alpha_NT_sp[[sp_plotted]])
    # plot_Neth_N_1_NT_sp[[sp_plotted]]<- as.data.frame(plot_Neth_N_1_NT_sp[[sp_plotted]])
    plot_Neth_N_alpha_Nalpha_sp[[sp_plotted]]<- as.data.frame(plot_Neth_N_alpha_Nalpha_sp[[sp_plotted]])
    # plot_Neth_N_1_Nalpha_sp[[sp_plotted]]<- as.data.frame(plot_Neth_N_1_Nalpha_sp[[sp_plotted]])
 #   plot_Nete_NT_sp[[sp_plotted]]<- as.data.frame(plot_Nete_NT_sp[[sp_plotted]])
 #   plot_Nete_sp[[sp_plotted]]<- as.data.frame(plot_Nete_sp[[sp_plotted]])
 #   plot_Nb_sp[[sp_plotted]]<- as.data.frame(plot_Nb_sp[[sp_plotted]])
 #   plot_Nb_NT_sp[[sp_plotted]]<- as.data.frame(plot_Nb_NT_sp[[sp_plotted]])
    plot_htz_sp[[sp_plotted]]<- as.data.frame(plot_htz_sp[[sp_plotted]])
    
    
    # reorder Fmax/Fstart factor for prettier plots
    plot_Neth_Nalpha_sp[[sp_plotted]]$Fish_press        <- factor(as.factor(plot_Neth_Nalpha_sp[[sp_plotted]]$Fish_press), levels = ideal_order_factor) 
    # plot_Neth_N1_sp[[sp_plotted]]$Fish_press            <- factor(as.factor(plot_Neth_N1_sp[[sp_plotted]]$Fish_press), levels = ideal_order_factor) 
    plot_Neth_N_alpha_NT_sp[[sp_plotted]]$Fish_press    <- factor(as.factor(plot_Neth_N_alpha_NT_sp[[sp_plotted]]$Fish_press), levels = ideal_order_factor) 
    # plot_Neth_N_1_NT_sp[[sp_plotted]]$Fish_press        <- factor(as.factor(plot_Neth_N_1_NT_sp[[sp_plotted]]$Fish_press), levels = ideal_order_factor) 
    plot_Neth_N_alpha_Nalpha_sp[[sp_plotted]]$Fish_press<- factor(as.factor(plot_Neth_N_alpha_Nalpha_sp[[sp_plotted]]$Fish_press), levels = ideal_order_factor) 
    # plot_Neth_N_1_Nalpha_sp[[sp_plotted]]$Fish_press    <- factor(as.factor(plot_Neth_N_1_Nalpha_sp[[sp_plotted]]$Fish_press), levels = ideal_order_factor) 
 #   plot_Nete_NT_sp[[sp_plotted]]$Fish_press            <- factor(as.factor(plot_Nete_NT_sp[[sp_plotted]]$Fish_press), levels = ideal_order_factor) 
 #   plot_Nete_sp[[sp_plotted]]$Fish_press               <- factor(as.factor(plot_Nete_sp[[sp_plotted]]$Fish_press), levels = ideal_order_factor) 
 #   plot_Nb_sp[[sp_plotted]]$Fish_press                 <- factor(as.factor(plot_Nb_sp[[sp_plotted]]$Fish_press), levels = ideal_order_factor) 
 #   plot_Nb_NT_sp[[sp_plotted]]$Fish_press              <- factor(as.factor(plot_Nb_NT_sp[[sp_plotted]]$Fish_press), levels = ideal_order_factor) 
    plot_htz_sp[[sp_plotted]]$Fish_press                <- factor(as.factor(plot_htz_sp[[sp_plotted]]$Fish_press), levels = ideal_order_factor) 
    
    
 ######### PLOTS NE temporal 
    reference <- subset(na.omit(plot_Neth_Nalpha_sp[[sp_plotted]]), Fish_press == paste0("Fmax=", F_vec[[sp]][1], " Fstart=", 
                                                                         F_vec[[sp]][2]))[, c("Ne_th_mean", "Ne_th_sd", "years")]
    plot_Neth_Nalpha[[sp]][[sp_plotted]] <- ggplot(data = na.omit(plot_Neth_Nalpha_sp[[sp_plotted]]), aes(x = years, group = Fish_press)) +
      geom_ribbon(aes(ymin = Ne_th_mean-Ne_th_sd, ymax = Ne_th_mean + Ne_th_sd), 
                  alpha = 0.5, fill = "grey70")+
      # 1) Layer: each panel gets its own group (via facet_wrap)
      geom_line(aes(y = Ne_th_mean)) +
      # 2) Layer: reference group A, drawn in all panels
      geom_ribbon(data = reference, 
                  fill = "grey70", alpha = 0.5,
                  aes(x = years, ymin = Ne_th_mean-Ne_th_sd, ymax = Ne_th_mean + Ne_th_sd, group = 1)) +
      geom_line(data = reference,
                aes(x = years, y = Ne_th_mean, group = 1),
                linetype = "dashed", linewidth = 0.7) +
      theme_bw()+
      labs(title = paste0(sp_plotted, " Ne theoretical (numerator: Nalpha),\nref is ", sp, F_vec))+
      facet_wrap(~ Fish_press)
    
    # reference <- subset(na.omit(plot_Neth_N1_sp[[sp_plotted]]), Fish_press == paste0("Fmax=", F_vec[[sp]][1], " Fstart=", 
    #                                                                      F_vec[[sp]][2]))[, c("Ne_th_mean", "Ne_th_sd", "years")]
    # plot_Neth_N1[[sp]][[sp_plotted]] <- ggplot(data = na.omit(plot_Neth_N1_sp[[sp_plotted]]), aes(x = years, group = Fish_press)) +
    #   theme_bw()+
    #   labs(title = paste0(sp_plotted, " Ne theoretical (numerator: N1),\nref is ", sp, F_vec))+
    #   geom_ribbon(aes(ymin = Ne_th_mean-Ne_th_sd, ymax = Ne_th_mean + Ne_th_sd), 
    #               alpha = 0.5, fill = "grey70")+
    #   # 1) Layer: each panel gets its own group (via facet_wrap)
    #   geom_line(aes(y = Ne_th_mean)) +
    #   # 2) Layer: reference group A, drawn in all panels
    #   geom_ribbon(data = reference, 
    #               fill = "grey70", alpha = 0.5,
    #               aes(x = years, ymin = Ne_th_mean-Ne_th_sd, ymax = Ne_th_mean + Ne_th_sd, group = 1)) +
    #   geom_line(data = reference,
    #             aes(x = years, y = Ne_th_mean, group = 1),
    #             linetype = "dashed", linewidth = 0.7) +
    #   facet_wrap(~ Fish_press)
    
    reference <- subset(na.omit(plot_Neth_N_alpha_NT_sp[[sp_plotted]]), Fish_press == paste0("Fmax=", F_vec[[sp]][1], " Fstart=", 
                                                                            F_vec[[sp]][2]))[, c("Ne_th_mean", "Ne_th_sd", "years")]
    plot_Neth_Nalpha_NT[[sp]][[sp_plotted]] <- ggplot(data = na.omit(plot_Neth_N_alpha_NT_sp[[sp_plotted]]), aes(x = years, group = Fish_press)) +
      theme_bw()+
      labs(title = paste0(sp_plotted, " Ne theoretical (numerator: Nalpha)\ndivided by NT, ref is ", sp, F_vec))+
      geom_ribbon(aes(ymin = Ne_th_mean-Ne_th_sd, ymax = Ne_th_mean + Ne_th_sd), 
                  alpha = 0.5, fill = "grey70")+
      # 1) Layer: each panel gets its own group (via facet_wrap)
      geom_line(aes(y = Ne_th_mean)) +
      # 2) Layer: reference group A, drawn in all panels
      geom_ribbon(data = reference, 
                  fill = "grey70", alpha = 0.5,
                  aes(x = years, ymin = Ne_th_mean-Ne_th_sd, ymax = Ne_th_mean + Ne_th_sd, group = 1)) +
      geom_line(data = reference,
                aes(x = years, y = Ne_th_mean, group = 1),
                linetype = "dashed", linewidth = 0.7) +
      facet_wrap(~ Fish_press)
    
    # reference <- subset(na.omit(plot_Neth_N_1_NT_sp[[sp_plotted]]), Fish_press == paste0("Fmax=", F_vec[[sp]][1], " Fstart=", 
    #                                                                             F_vec[[sp]][2]))[, c("Ne_th_mean", "Ne_th_sd", "years")]
    # plot_Neth_N1_NT[[sp]][[sp_plotted]] <- ggplot(data = na.omit(plot_Neth_N_1_NT_sp[[sp_plotted]]), aes(x = years, group = Fish_press)) +
    #   theme_bw()+
    #   labs(title = paste0(sp_plotted, " Ne theoretical (numerator: N1)\ndivided by NT, ref is ", sp, F_vec))+
    #   geom_ribbon(aes(ymin = Ne_th_mean-Ne_th_sd, ymax = Ne_th_mean + Ne_th_sd), 
    #               alpha = 0.5, fill = "grey70")+
    #   # 1) Layer: each panel gets its own group (via facet_wrap)
    #   geom_line(aes(y = Ne_th_mean)) +
    #   # 2) Layer: reference group A, drawn in all panels
    #   geom_ribbon(data = reference, 
    #               fill = "grey70", alpha = 0.5,
    #               aes(x = years, ymin = Ne_th_mean-Ne_th_sd, ymax = Ne_th_mean + Ne_th_sd, group = 1)) +
    #   geom_line(data = reference,
    #             aes(x = years, y = Ne_th_mean, group = 1),
    #             linetype = "dashed", linewidth = 0.7) +
    #   facet_wrap(~ Fish_press)
    
    reference <- subset(na.omit(plot_Neth_N_alpha_Nalpha_sp[[sp_plotted]]), Fish_press == paste0("Fmax=", F_vec[[sp]][1], " Fstart=", 
                                                                                F_vec[[sp]][2]))[, c("Ne_th_mean", "Ne_th_sd", "years")]
    plot_Neth_Nalpha_Nalpha[[sp]][[sp_plotted]] <- ggplot(data = na.omit(plot_Neth_N_alpha_Nalpha_sp[[sp_plotted]]), aes(x = years, group = Fish_press)) +
      theme_bw()+
      labs(title = paste0(sp_plotted, " Ne theoretical (numerator: Nalpha)\ndivided by Nalpha, ref is ", sp, F_vec))+
      geom_ribbon(aes(ymin = Ne_th_mean-Ne_th_sd, ymax = Ne_th_mean + Ne_th_sd), 
                  alpha = 0.5, fill = "grey70")+
      # 1) Layer: each panel gets its own group (via facet_wrap)
      geom_line(aes(y = Ne_th_mean)) +
      # 2) Layer: reference group A, drawn in all panels
      geom_ribbon(data = reference, 
                  fill = "grey70", alpha = 0.5,
                  aes(x = years, ymin = Ne_th_mean-Ne_th_sd, ymax = Ne_th_mean + Ne_th_sd, group = 1)) +
      geom_line(data = reference,
                aes(x = years, y = Ne_th_mean, group = 1),
                linetype = "dashed", linewidth = 0.7) +
      facet_wrap(~ Fish_press)
    # 
    # reference <- subset(na.omit(plot_Neth_N_1_Nalpha_sp[[sp_plotted]]), Fish_press == paste0("Fmax=", F_vec[[sp]][1], " Fstart=", 
    #                                                                             F_vec[[sp]][2]))[, c("Ne_th_mean", "Ne_th_sd", "years")]
    # plot_Neth_N1_Nalpha[[sp]][[sp_plotted]] <- ggplot(data = na.omit(plot_Neth_N_1_Nalpha_sp[[sp_plotted]]), aes(x = years, group = Fish_press)) +
    #   theme_bw()+
    #   labs(title = paste0(sp_plotted, " Ne theoretical (numerator: N1)\ndivided by Nalpha, ref is ", sp, F_vec))+
    #   geom_ribbon(aes(ymin = Ne_th_mean-Ne_th_sd, ymax = Ne_th_mean + Ne_th_sd), 
    #               alpha = 0.5, fill = "grey70")+
    #   # 1) Layer: each panel gets its own group (via facet_wrap)
    #   geom_line(aes(y = Ne_th_mean)) +
    #   # 2) Layer: reference group A, drawn in all panels
    #   geom_ribbon(data = reference, 
    #               fill = "grey70", alpha = 0.5,
    #               aes(x = years, ymin = Ne_th_mean-Ne_th_sd, ymax = Ne_th_mean + Ne_th_sd, group = 1)) +
    #   geom_line(data = reference,
    #             aes(x = years, y = Ne_th_mean, group = 1),
    #             linetype = "dashed", linewidth = 0.7) +
    #   facet_wrap(~ Fish_press)
    # 
    
    
    
  ############ PLOTS Ne temporal  
 #   reference <- subset(na.omit(plot_Nete_sp[[sp_plotted]]), Fish_press == paste0("Fmax=", F_vec[[sp]][1], " Fstart=", 
 #                                                                        F_vec[[sp]][2]))[, c("Ne_tp_mean", "Ne_tp_sd", "years")]
 #   plot_Nete[[sp]][[sp_plotted]] <- ggplot(data = na.omit(plot_Nete_sp[[sp_plotted]]), aes(x = years, group = Fish_press)) +
 #     theme_bw()+
 #     labs(title = paste0(sp_plotted, " Ne temporal, ref is ", sp, F_vec))+
 #     geom_ribbon(aes(ymin = Ne_tp_mean-Ne_tp_sd, ymax = Ne_tp_mean + Ne_tp_sd), 
 #                 alpha = 0.5, fill = "grey70")+
 #     # 1) Layer: each panel gets its own group (via facet_wrap)
 #     geom_line(aes(y = Ne_tp_mean)) +
 #     # 2) Layer: reference group A, drawn in all panels
 #     geom_ribbon(data = reference, 
 #                 fill = "grey70", alpha = 0.5,
 #                 aes(x = years, ymin = Ne_tp_mean-Ne_tp_sd, ymax = Ne_tp_mean + Ne_tp_sd, group = 1)) +
 #     geom_line(data = reference,
 #               aes(x = years, y = Ne_tp_mean, group = 1),
 #               linetype = "dashed", linewidth = 0.7) +
 #     facet_wrap(~ Fish_press)
 #   reference <- subset(na.omit(plot_Nete_NT_sp[[sp_plotted]]), Fish_press == paste0("Fmax=", F_vec[[sp]][1], " Fstart=", 
 #                                                                        F_vec[[sp]][2]))[, c("Ne_tp_mean", "Ne_tp_sd", "years")]
 #   plot_Nete_NT[[sp]][[sp_plotted]] <- ggplot(data = na.omit(plot_Nete_NT_sp[[sp_plotted]]), aes(x = years, group = Fish_press)) +
 #     theme_bw()+
 #     labs(title = paste0(sp_plotted, " Ne temporal divided by NT,\nref is ", sp, F_vec))+
 #     geom_ribbon(aes(ymin = Ne_tp_mean-Ne_tp_sd, ymax = Ne_tp_mean + Ne_tp_sd), 
 #                 alpha = 0.5, fill = "grey70")+
 #     # 1) Layer: each panel gets its own group (via facet_wrap)
 #     geom_line(aes(y = Ne_tp_mean)) +
 #     # 2) Layer: reference group A, drawn in all panels
 #     geom_ribbon(data = reference, 
 #                 fill = "grey70", alpha = 0.5,
 #                 aes(x = years, ymin = Ne_tp_mean-Ne_tp_sd, ymax = Ne_tp_mean + Ne_tp_sd, group = 1)) +
 #     geom_line(data = reference,
 #               aes(x = years, y = Ne_tp_mean, group = 1),
 #               linetype = "dashed", linewidth = 0.7) +
 #     facet_wrap(~ Fish_press)
 #   
 # 
 # ########### PLOTS Nb 
 #   reference <- subset(na.omit(plot_Nb_sp[[sp_plotted]]), Fish_press == paste0("Fmax=", F_vec[[sp]][1], " Fstart=", 
 #                                                                        F_vec[[sp]][2]))[, c("Nb_mean", "Nb_sd", "years")]
 #   plot_Nb[[sp]][[sp_plotted]] <- ggplot(data = na.omit(plot_Nb_sp[[sp_plotted]]), aes(x = years, group = Fish_press)) +
 #     theme_bw()+
 #     labs(title = paste0(sp_plotted, " Nb, ref is ", sp, F_vec))+
 #     geom_ribbon(aes(ymin = Nb_mean-Nb_sd, ymax = Nb_mean + Nb_sd), 
 #                 alpha = 0.5, fill = "grey70")+
 #     # 1) Layer: each panel gets its own group (via facet_wrap)
 #     geom_line(aes(y = Nb_mean)) +
 #     # 2) Layer: reference group A, drawn in all panels
 #     geom_ribbon(data = reference, 
 #                 fill = "grey70", alpha = 0.5,
 #                 aes(x = years, ymin = Nb_mean-Nb_sd, ymax = Nb_mean + Nb_sd, group = 1)) +
 #     geom_line(data = reference,
 #               aes(x = years, y = Nb_mean, group = 1),
 #               linetype = "dashed", linewidth = 0.7) +
 #     facet_wrap(~ Fish_press)    
 #   reference <- subset(na.omit(plot_Nb_NT_sp[[sp_plotted]]), Fish_press == paste0("Fmax=", F_vec[[sp]][1], " Fstart=", 
 #                                                                        F_vec[[sp]][2]))[, c("Nb_mean", "Nb_sd", "years")]
 #   plot_Nb_NT[[sp]][[sp_plotted]] <- ggplot(data = na.omit(plot_Nb_NT_sp[[sp_plotted]]), aes(x = years, group = Fish_press)) +
 #     theme_bw()+
 #     labs(title = paste0(sp_plotted, " Nb divided by NT,\nref is ", sp, F_vec))+
 #     geom_ribbon(aes(ymin = Nb_mean-Nb_sd, ymax = Nb_mean + Nb_sd), 
 #                 alpha = 0.5, fill = "grey70")+
 #     # 1) Layer: each panel gets its own group (via facet_wrap)
 #     geom_line(aes(y = Nb_mean)) +
 #     # 2) Layer: reference group A, drawn in all panels
 #     geom_ribbon(data = reference, 
 #                 fill = "grey70", alpha = 0.5,
 #                 aes(x = years, ymin = Nb_mean-Nb_sd, ymax = Nb_mean + Nb_sd, group = 1)) +
 #     geom_line(data = reference,
 #               aes(x = years, y = Nb_mean, group = 1),
 #               linetype = "dashed", linewidth = 0.7) +
 #     facet_wrap(~ Fish_press)
 #   
    
    
 ############# PLOTS heterozygosity
    reference <- subset(na.omit(plot_htz_sp[[sp_plotted]]), Fish_press == paste0("Fmax=", F_vec[[sp]][1], " Fstart=", 
                                                                         F_vec[[sp]][2]))[, c("htz_mean", "htz_sd", "years")]
    plot_htz[[sp]][[sp_plotted]] <- ggplot(data = na.omit(plot_htz_sp[[sp_plotted]]), aes(x = years, group = Fish_press)) +
      theme_bw()+
      labs(title = paste0(sp_plotted, " HTZ theoretical, ref is ", sp, F_vec))+
      geom_ribbon(aes(ymin = htz_mean-htz_sd, ymax = htz_mean + htz_sd), 
                  alpha = 0.5, fill = "grey70")+
      # 1) Layer: each panel gets its own group (via facet_wrap)
      geom_line(aes(y = htz_mean)) +
      # 2) Layer: reference group A, drawn in all panels
      geom_ribbon(data = reference, 
                  fill = "grey70", alpha = 0.5,
                  aes(x = years, ymin = htz_mean-htz_sd, ymax = htz_mean + htz_sd, group = 1)) +
      geom_line(data = reference,
                aes(x = years, y = htz_mean, group = 1),
                linetype = "dashed", linewidth = 0.7) +
      facet_wrap(~ Fish_press)
    
  }
  
  if(!dir.exists(pathout)) {dir.create(pathout, recursive = T) }
  
  # pdf(file = paste0(pathout, "/plot_Neth_N1_", sp, ".pdf"), width = 10, height = 8)
  # par(mfrow=c(1,1))
  # print(plot_Neth_N1[[sp]])
  # dev.off()
  
  pdf(file = paste0(pathout, "/plot_Neth_Nalpha_", sp, ".pdf"), width = 10, height = 8)
  par(mfrow=c(1,1))
  print(plot_Neth_Nalpha[[sp]])
  dev.off()
  # 
  # pdf(file = paste0(pathout, "/plot_Neth_N1_NT_", sp, ".pdf"), width = 10, height = 8)
  # par(mfrow=c(1,1))
  # print(plot_Neth_N1_NT[[sp]])
  # dev.off()
  
  pdf(file = paste0(pathout, "/plot_Neth_Nalpha_NT_", sp, ".pdf"), width = 10, height = 8)
  par(mfrow=c(1,1))
  print(plot_Neth_Nalpha_NT[[sp]])
  dev.off()
  
  pdf(file = paste0(pathout, "/plot_Neth_Nalpha_Nalpha_", sp, ".pdf"), width = 10, height = 8)
  par(mfrow=c(1,1))
  print(plot_Neth_Nalpha_Nalpha[[sp]])
  dev.off()
  
  # pdf(file = paste0(pathout, "/plot_Neth_N1_Nalpha_", sp, ".pdf"), width = 10, height = 8)
  # par(mfrow=c(1,1))
  # print(plot_Neth_N1_Nalpha[[sp]])
  # dev.off()
  
 # pdf(file = paste0(pathout, "/plot_Nete_", sp, ".pdf"), width = 10, height = 8)
 # par(mfrow=c(1,1))
 # print(plot_Nete[[sp]])
 # dev.off()
 # 
 # pdf(file = paste0(pathout, "/plot_Nete_NT_", sp, ".pdf"), width = 10, height = 8)
 # par(mfrow=c(1,1))
 # print(plot_Nete_NT[[sp]])
 # dev.off()
 # 
 # pdf(file = paste0(pathout, "/plot_Nb_", sp, ".pdf"), width = 10, height = 8)
 # par(mfrow=c(1,1))
 # print(plot_Nb[[sp]])
 # dev.off()
 # 
 # pdf(file = paste0(pathout, "/plot_Nb_NT_", sp, ".pdf"), width = 10, height = 8)
 # par(mfrow=c(1,1))
 # print(plot_Nb_NT[[sp]])
 # dev.off()
 # 
  pdf(file = paste0(pathout, "/plot_htz_", sp, ".pdf"), width = 10, height = 8)
  par(mfrow=c(1,1))
  print(plot_htz[[sp]])
  dev.off()
  
  # save(plot_Neth_Nalpha, file=paste0(pathoutputRDATA, "/plot_Neth_Nalpha.RData"))
  # save(plot_Neth_N1, file=paste0(pathoutputRDATA, "/plot_Neth_N1.RData"))
  # save(plot_Neth_Nalpha_NT, file=paste0(pathoutputRDATA, "/plot_Neth_Nalpha_NT.RData"))
  # save(plot_Neth_N1_NT, file=paste0(pathoutputRDATA, "/plot_Neth_N1_NT.RData"))
  # save(plot_Neth_Nalpha_Nalpha, file=paste0(pathoutputRDATA, "/plot_Neth_Nalpha_Nalpha.RData"))
  # save(plot_Neth_N1_Nalpha, file=paste0(pathoutputRDATA, "/plot_Neth_N1_Nalpha.RData"))
  # save(plot_Nete, file=paste0(pathoutputRDATA, "/plot_Nete.RData"))
  # save(plot_Nete_NT, file=paste0(pathoutputRDATA, "/plot_Nete_NT.RData"))
  # save(plot_Nb, file=paste0(pathoutputRDATA, "/plot_Nb.RData"))
  # save(plot_Nb_NT, file=paste0(pathoutputRDATA, "/plot_Nb_NT.RData"))
  # save(plot_htz, file=paste0(pathoutputRDATA, "/plot_htz.RData"))
  
} #sp
}# loop
# for (sp in species){
#   pathout <- paste0(configDir, "/output/PLOTS/All_genetic_drift_indicators")
#   dir.create(pathout, recursive = T)
#   
#   pdf(file = paste0(pathout, "/plot_Neth_", sp, ".pdf"), width = 10, height = 8)
#   par(mfrow=c(1,1))
#   print(plot_Neth[[sp]])
#   dev.off()
#   
#   pdf(file = paste0(pathout, "/plot_Nete_", sp, ".pdf"), width = 10, height = 8)
#   par(mfrow=c(1,1))
#   print(plot_Nete[[sp]])
#   dev.off()
#   
#   pdf(file = paste0(pathout, "/plot_Nb_", sp, ".pdf"), width = 10, height = 8)
#   par(mfrow=c(1,1))
#   print(plot_Nb[[sp]])
#   dev.off()
#   
#   pdf(file = paste0(pathout, "/plot_htz_", sp, ".pdf"), width = 10, height = 8)
#   par(mfrow=c(1,1))
#   print(plot_htz[[sp]])
#   dev.off()
# }







