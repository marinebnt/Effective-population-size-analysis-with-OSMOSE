

PLOT_trait_replicates2 <- function(SP){

# setwd("X:/osmose-ns")

repo        = "Simu_LONG"
multip_name = "y"
multipliers = c(1,0.5,0.75,1.25,1.5)

trait_dim= c("SSB", "sizeMature", "abundanceMature", "fecundityByAge", "meanSizeByAge", 
     #"MpredbyY", "MstarbyY", "MnatbyY", "FbyY", "ZbyY", "MforbyY", "DisbyY", "eggNumberByAge", 
     "abundanceByAge", "biomass", "mortality", "school")
#trait_dim = "mortality"


ideal_order_factor <- paste0("Fmax=", rep(c(0.5,0.75,1,1.25,1.5), 5),
                             " Fstart=", c(rep(0.5, 5),rep(0.75, 5),rep(1, 5),rep(1.25, 5),rep(1.5, 5)))



path_in <- c("C:/Users/mbeneat/Documents/osmose/Chap2/osmose-ns/osmose-ns-v4.4.1_clim/output/test")
# path_in <- c("Y:/outputs-meso-marine")
# path_in <- c("/home1/datawork/mbeneat/SIMU_CHAP2_FIN")
list_d <- list.dirs(path = path_in, recursive = F, full.names = T)

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

pathout <- paste0(configDir, "/output/PLOTS/Trait_variables_by_year")

plot_Neth_Nalpha    <- list()
F_vec     <- list()


spec_names <- species


for (sp in SP){
  data_complete <- NULL
  data_all <- NULL

  # Plot only target species + one additional species (to reduce plot size)
  sp_i <- which(spec_names == sp)-1
  species = spec_names[c(c(sp_i+1), ifelse(sp_i<14, sp_i+2, 1))]  
  # Locate directories for the focal species
  sp_dir <- list_d[grep(pattern = paste0("sp", sp_i), list_d)]
  
  plot_Neth_N_alpha_NT_sp     <- list()
  Fmax_1 = -99
  Fstart_1 = -99
  age_Mature_df_all = NULL
  for (Fmax in multipliers){
    Fmax_dir <- sp_dir[grep(pattern = paste0("Fmax", Fmax, "_"), sp_dir)]
    for (Fstart in multipliers){
      dir <- Fmax_dir[grep(pattern = paste0("Fstart", Fstart, "_sp", sp_i, ".R"), Fmax_dir)]

      if (length(list.files(dir))<1){
        cat("\ndata missing", dir, "\n")
        next}
      if (Fmax_1==-99 && Fstart_1 ==-99){
        Fmax_1 = Fmax
        Fstart_1 = Fstart
        F_vec[[sp]] <- c(Fmax_1, Fstart_1)
      }

      # initialisation
      simulation = dir
      cat("\n\t---> folder Fmax", Fmax, " Fstart", Fstart, " sp", sp, "\n")

      ns1 = NULL
      outputDir_Simu  = simulation
      outputDir_Ne    = file.path(configDir, "output", "outputdirNe")
      #if(!dir.exists(outputDir_Ne)) {dir.create(outputDir_Ne, recursive = T)}
      ns1 = read_osmose(path = paste0(outputDir_Simu, "/"))

      # age mature
      age_Mature = ns1[["ageMature"]]
      age_Mature = lapply(seq(dim(age_Mature)[3]), function(x) as.data.frame(age_Mature[ , , x]))
      for (i in 1:length(age_Mature)) {
        age_Mature[[i]]$rep <- i
        age_Mature[[i]]$time <- rownames(age_Mature[[i]])
      }
      age_Mature_df = do.call(rbind, age_Mature)
      age_Mature_df <- as.data.frame(pivot_longer(age_Mature_df, cols = -c(rep, time), names_to = "Species", values_to = "value"))
      age_Mature_df$TRAIT <- "age_Mature"
      age_Mature_df$Fish_press = paste0("Fmax=", Fmax, " Fstart=", Fstart)
      if (exists("age_Mature_df_all")){
        age_Mature_df_all <- rbind(age_Mature_df_all, age_Mature_df)
      } else {
        age_Mature_df_all <- age_Mature_df
      }

      for (TRAIT in trait_dim){
        cat("Processing trait:", TRAIT, "\n")
        data_tot = NULL
        data <- ns1[[TRAIT]]
        for (i in 0:c(nb_rep-1)){ 
          # get the abundance by age to be able to weight the variables by the abundance of each age class (for meanSizeByAge for example)
          data_tot_abund <- list()
          data_abund_ab <- ns1[["abundanceByAge"]]
          data_tot_abund_t <- as.tibble(data_abund_ab[[species[1]]][,,i+1])
          # data_tot_abund_t$Time <- rownames(data_tot_abund_t)
          data_tot_abund[[species[1]]] <- data_tot_abund_t #pivot_longer(data_tot_abund_t, cols = -Time, names_to = "Age", values_to = "abundance")
          for (spi in species[-1]){
            data_tot_abund_t <- as.tibble(data_abund_ab[[spi]][,,i+1])
            # data_tot_abund_t$Time <- rownames(data_tot_abund_t)
            data_tot_abund[[spi]] <- data_tot_abund_t #pivot_longer(data_tot_abund_t, cols = -Time, names_to = "Age", values_to = "abundance")
          }

          if (TRAIT == "mortality"){
            data_tot <- data.frame("Time"= rownames(ns1$mortalityByAge[[species[1]]]$Mpred))
            rownames(data_tot) = data_tot$Time
            names_list <- names(ns1$mortalityByAge[[species[1]]])
            data_tot[, species[1]] = rowSums(sapply(names_list, FUN = function(x) {
              data_t = NULL
              for (y in 1:nrow(data_tot)){
                if (x == "Mnat") {
                  data_t[y] <- 0
                } else {
                  data_t[y] <- weighted.mean(x=ns1$mortalityByAge[[species[1]]][[x]][y,,(i+1)], 
                                                              w=as.matrix(data_tot_abund[[species[1]]][y,]), na.rm = TRUE)
                }
              }
              # if(x == "Mnat") rep(0, nrow(data_tot)) else rowWeightedMeans(x=ns1$mortalityByAge[[species[1]]][[x]][,,(i+1)], 
              #                                                              w=data_tot_abund[[species[1]]], na.rm = TRUE)
              return(data_t)
              }))
            for (spi in species[-1]){
              data_tot[, spi] = rowSums(sapply(names_list, FUN = function(x) {
                data_t = NULL
                for (y in 1:nrow(data_tot)){
                  if (x == "Mnat") {
                    data_t[y] <- 0
                  } else {
                    data_t[y] <- weighted.mean(x=ns1$mortalityByAge[[spi]][[x]][y,,(i+1)], 
                                                                w=as.matrix(data_tot_abund[[spi]][y,]), na.rm = TRUE)
                  }
                }
                return(data_t)
              }))
            }
          }
          if (TRAIT %in% c("abundanceMature", "sizeMature")){
            if (TRAIT == "abundanceMature") data_tot <- data.frame("Time"= rownames(ns1$abundanceByAge[[species[1]]]))
            if (TRAIT == "sizeMature")      data_tot <- data.frame("Time"= rownames(ns1$meanSizeByAge[[species[1]]]))
            rownames(data_tot) = data_tot$Time
            age_mat_sp_rep = age_Mature[[i+1]][, species[1]]
            vec_y = c()
            for (y in 1:length(age_mat_sp_rep)){
               if (TRAIT == "abundanceMature") {
                  vec_tempo <- ns1$abundanceByAge[[species[1]]][y,which(as.numeric(colnames(ns1$abundanceByAge[[species[1]]]))== ceiling(age_mat_sp_rep[y])), i+1]
                  vec_y[y] <- ifelse(length(vec_tempo)==0, NA, vec_tempo)
               }
               if (TRAIT == "sizeMature")      {
                  vec_tempo <- ns1$meanSizeByAge[[species[1]]][y,which(as.numeric(colnames(ns1$meanSizeByAge[[species[1]]]))== ceiling(age_mat_sp_rep[y])), i+1]
                  vec_y[y] <- ifelse(length(vec_tempo)==0, NA, vec_tempo)
               }
            }
            data_tot[,species[1]] <- vec_y
            for (spi in species[-1]){
              age_mat_sp_rep = age_Mature[[i+1]][, spi]
              vec_y = c()
              for (y in 1:length(age_mat_sp_rep)){
                if (TRAIT == "abundanceMature") {
                  vec_tempo <- ns1$abundanceByAge[[spi]][y,which(as.numeric(colnames(ns1$abundanceByAge[[spi]]))== ceiling(age_mat_sp_rep[y])), i+1]
                  vec_y[y] <- ifelse(length(vec_tempo)==0, NA, vec_tempo) 
                }
                if (TRAIT == "sizeMature")      {
                  tempo_vec_i <- ns1$meanSizeByAge[[spi]][y,which(as.numeric(colnames(ns1$meanSizeByAge[[spi]]))== max(c(ceiling(age_mat_sp_rep[y]), mean(na.omit(age_mat_sp_rep))), na.rm=T)), i+1]
                  vec_y[y] <- ifelse(length(tempo_vec_i)==0, NA, tempo_vec_i)
                }
              }
              data_tot[,spi] <- vec_y
            }
          }
          if (TRAIT %in% c("meanSizeByAge", "abundanceByAge")) {
            data_tot_t <- data[[species[1]]]
            data_tot <- data.frame("Time"= rownames(data_tot_t[,,i+1]))
            rownames(data_tot) = data_tot$Time
            if (TRAIT == ("meanSizeByAge"))  data_tot[,species[1]] <- rowMeans(data_tot_t[,,i+1],   na.rm = TRUE)
            if (TRAIT == ("abundanceByAge")) data_tot[,species[1]] <- data_tot_t[,2,i+1]
            for (spi in species[-1]){
              data_tot_t <- data[[spi]]
              if (TRAIT == ("meanSizeByAge")) data_tot[,spi] <- rowMeans(data_tot_t[,,i+1],   na.rm = TRUE)
              if (TRAIT == ("abundanceByAge")) data_tot[,spi] <- data_tot_t[,2,i+1]
            }
          }
          if (TRAIT %in% c("MpredbyY", "MstarbyY", "MnatbyY", "FbyY", "ZbyY", "MforbyY", "DisbyY")) {
            TRAIT_t <- stringr::str_remove(TRAIT, "byY")
            data_tot_all <- read.csv(paste0(outputDir_Simu, "/Mortality/ns_mortalityRateDistribByAge-", species[1], "_Simu", i, ".csv"), skip=c(1))
            data_tot_all <- data_tot_all[, -ncol(data_tot_all)] #remove useless column, also called Age
            data_tot   <- data.frame("Time"= unique(data_tot_all$Time))
            for (y in unique(data_tot_all$Time)){
              data_tot_t <- pivot_wider(data_tot_all[which(data_tot_all$Time == y), which(colnames(data_tot_all) %in% c(TRAIT_t, "Time", "Age"))], 
                                        names_from = "Age", values_from = TRAIT_t) %>%
                               select(-Time)
              data_tot[which(data_tot$Time == y),species[1]] <- weighted.mean(x=data_tot_t,
                                                        w = as.matrix(data_tot_abund[[species[1]]][which(rownames(data_tot_abund[[species[1]]]) == c(y-80)),]), na.rm = TRUE)
            }
            for (spi in species[-1]){
              data_tot_all <- read.csv(paste0(outputDir_Simu, "/Mortality/ns_mortalityRateDistribByAge-", spi, "_Simu", i, ".csv"), skip=c(1))
              for (y in unique(data_tot_all$Time)){
                data_tot_t <- pivot_wider(data_tot_all[which(data_tot_all$Time == y), which(colnames(data_tot_all) %in% c(TRAIT_t, "Time", "Age"))], 
                                          names_from = "Age", values_from = TRAIT_t) %>%
                                select(-Time)
                data_tot[which(data_tot$Time == y),spi] <- weighted.mean(x=data_tot_t,
                                                          w = as.matrix(data_tot_abund[[spi]][which(rownames(data_tot_abund[[spi]]) == c(y-80)),]), na.rm = TRUE)
              }
            }
          }
          if (TRAIT == ("fecundityByAge")) {
            eggNumber <- read.csv(paste0(outputDir_Simu, "/AgeIndicators/ns_fecundityDistribByAge_Simu", i, ".csv"), skip=c(1))
            eggNumber <- eggNumber[, which(colnames(eggNumber)%in% c("Age", species[1], "Time"))]
            abundance <- as.data.frame(ns1[["abundanceByAge"]][[species[1]]][,,i+1])
            abundance$Time <- as.numeric(rownames(abundance))
            abundance <- pivot_longer(abundance, cols =  -Time, names_to = "Age", values_to = "abundance")
            abundance$Age <- as.numeric(abundance$Age)
            full_join(eggNumber, abundance, by = c("Time", "Age")) %>%
              mutate(fecundity = get(species[1])/abundance) -> fecundity_data
            fecundity_data$fecundity[-(grep('^-?[0-9.]+$', fecundity_data$fecundity))] <- 0
            data_tot <- fecundity_data[, which(colnames(fecundity_data) %in% c("Age", "Time", "fecundity"))]
            colnames(data_tot) <- c("Time", "Age", species[1])
            for (spi in species[-1]){
              eggNumber <- read.csv(paste0(outputDir_Simu, "/AgeIndicators/ns_fecundityDistribByAge_Simu", i, ".csv"), skip=c(1))
              eggNumber <- eggNumber[, which(colnames(eggNumber)%in% c("Age", spi, "Time"))]
              abundance <- as.data.frame(ns1[["abundanceByAge"]][[spi]][,,i+1])
              abundance$Time <- as.numeric(rownames(abundance))
              abundance <- pivot_longer(abundance, cols =  -Time, names_to = "Age", values_to = "abundance")
              abundance$Age <- as.numeric(abundance$Age)
              full_join(eggNumber, abundance, by = c("Time", "Age")) %>%
                mutate(fecundity = get(spi)/abundance) -> fecundity_data
              fecundity_data$fecundity[-(grep('^-?[0-9.]+$', fecundity_data$fecundity))] <- 0
              data_tot <- full_join(data_tot, fecundity_data[, which(colnames(fecundity_data) %in% c("Age", "Time", "fecundity"))],
                                   by = c("Age", "Time"))
              colnames(data_tot)[ncol(data_tot)] <- spi
            }
            data_tot_w <- data.frame("Time"= unique(data_tot$Time))
            for (spi in species){
              for (y in unique(data_tot$Time)){
                  data_tot_t <- pivot_wider(data_tot[which(data_tot$Time == y), which(colnames(data_tot) %in% c(spi, "Time", "Age"))], 
                                            names_from = "Age", values_from = spi) %>%
                                  select(-Time)
                  data_tot_w[which(data_tot_w$Time == y),spi] <- weighted.mean(x=data_tot_t,
                                                            w = as.matrix(data_tot_abund[[spi]][which(rownames(data_tot_abund[[spi]]) == c(y-80)),]), na.rm = TRUE)
              }
            }
            data_tot <- data_tot_w
          }
          if (TRAIT %in% c("SSB")) {
            data_tot_all <- read.csv(paste0(outputDir_Simu, "/ns_SSB_Simu", i, ".csv"), skip=c(1))
            data_tot <- data_tot_all[, which(colnames(data_tot_all) %in% c("Time", species))]
          }
          if (TRAIT %in% c("school", 'eggNumberByAge')) {
            if (TRAIT == c("eggNumberByAge")) {
              data_tot_all <- read.csv(paste0(outputDir_Simu, "/AgeIndicators/ns_fecundityDistribByAge_Simu", i, ".csv"), skip=c(1))
              data_tot_sum <- data_tot_all %>%
                group_by(Time) %>%
                summarise(across(all_of(species), ~ sum(.x, na.rm = TRUE))) %>%
                as.data.frame()
            }
            if (TRAIT == c("school")) {
              data_tot_all <- read.csv(paste0(outputDir_Simu, "/ns_nschoolDistribByAge_Simu", i, ".csv"), skip=c(1))
              data_tot_sum <- data_tot_all[which(data_tot_all$Age == 1), ] # only age 1
            }
            data_tot <- as.data.frame(data_tot_sum[, which(colnames(data_tot_sum) %in% c(species, "Time"))])
          }
          if (!(TRAIT %in% c("mortality", "school", 'eggNumber', 'abundanceByAge', "eggNumberByAge",
                              "abundanceMature", "fecundityByAge", "meanSizeByAge", "SSB", "sizeMature",
                              "MpredbyY", "MstarbyY", "MnatbyY", "FbyY", "ZbyY", "MforbyY", "DisbyY"))){
            data_tot <- as.data.frame(data[,,i+1])
            data_tot$Time <- rownames(data_tot)
            data_tot <- as.data.frame(data_tot[, which(colnames(data_tot) %in% c(species, "Time"))])
          }
          data_tot$time <- data_tot$Time
          data_tot$Time <- NULL
          data_tot <- data_tot%>%
            pivot_longer(!time, names_to = "Species", values_to = "value")
          data_tot$rep <- i+1
          data_tot$Fish_press = paste0("Fmax=", Fmax, " Fstart=", Fstart)
          data_tot$TRAIT = TRAIT

          if(TRAIT == trait_dim[1] && i==0){
            data_all <- rbind(data_tot, age_Mature_df_all)
          } else {
            data_all <- rbind(data_all, data_tot, age_Mature_df_all)
          }
        }#rep
      }#TRAIT

      if (Fmax == Fmax_1 & Fstart == Fstart_1) {
        data_complete <- (data_all)
      }else {
        data_complete <- rbind(data_complete, data_all)
      }
      # p(paste0("\n\t--->done: folder Fmax", Fmax, " Fstart", Fstart, " sp", sp, "\n"))
    }#Fstart
  }#Fmax
  # PLOT the mean and the std deviation 
     # adding the sum of mortalities by age 
  
  data_complete <- (data_complete)

  pathoutputRDATA <-paste0(pathout, "/RDATA/Fmax", Fmax, "_Fstart", Fstart, "_sp", sp)
  if (!dir.exists(pathoutputRDATA)){dir.create(pathoutputRDATA, recursive = T)}  
  # if (!file.exists(paste0(pathoutputRDATA, "/data_traits_summary.RData"))){
    save(data_complete, file=paste0(pathoutputRDATA, "/data_traits_summary.RData"))
  # } else {
    # load(file=paste0(pathoutputRDATA, "/data_traits_summary.RData"))
    #   data_save <- rbind(loaded_data, data_complete)
    #   save(data_save, file=paste0(pathoutputRDATA, "/data_traits_summary.RData"))
  # }  

  # coherently with the transformations applied earlier
  new_trait_dim <- trait_dim
  for (nm in trait_dim){
    ifelse_name <- ifelse(nm %in% c("eggNumberByAge"), "sum", 
                   ifelse(nm %in% c("mortality", "meanSizeByAge", "fecundityByAge", "MpredbyY", "MstarbyY", "MnatbyY", "FbyY", "ZbyY", "MforbyY", "DisbyY"), "mean", 
                   ifelse(nm %in% c("sizeMature", "abundanceMature"), "age_mat", 
                   ifelse(nm %in% c("school", "abundanceByAge"), "age_1", "total"))))
    new_nm <- paste0(ifelse_name, "_of_", nm, "_per_y")
    new_trait_dim[which(trait_dim == nm)] <- new_nm
  }
  new_trait_dim <- c(new_trait_dim, "age_Mature")
  # if (all(c("MpredbyY", "MstarbyY", "FbyY", "ZbyY", "MforbyY", "DisbyY") %in% trait_dim)){
  #   data_tempo_mortbyage <-
  #     data_complete[which(data_complete$TRAIT %in% c("MpredbyY", "MstarbyY", "FbyY", "ZbyY", "MforbyY", "DisbyY")),]
  #   data_tempo_mortbyage <- data_tempo_mortbyage %>%
  #                 group_by(Fish_press, Species, time, rep) %>%
  #                 summarise(sum = sum(value, na.rm = TRUE))
  #   data_tempo_mortbyage$TRAIT = c("mortality")
  #   data_complete <- rbind(data_complete,
  #                              data.frame("time"= data_tempo_mortbyage$time,
  #                                         "Species"= data_tempo_mortbyage$Species,
  #                                         "value"= data_tempo_mortbyage$sum,
  #                                         "rep"= data_tempo_mortbyage$rep,
  #                                         "Fish_press"= data_tempo_mortbyage$Fish_press,
  #                                         "TRAIT"= data_tempo_mortbyage$TRAIT))
  #   trait_dim = c(trait_dim, "mortalityByAge")
  # }

  #data_complete$value[grep(data_complete$TRAIT, pattern="eggNumber")] <- log(data_complete$value[grep(data_complete$TRAIT, pattern= "eggNumber")]+10^{-8}) # easier to plot
  data_tot_summary <- data_complete %>%
    group_by(time, Species, Fish_press, TRAIT) %>%
    summarise(mean = mean(value, na.rm = TRUE),
              sd   = sd(value, na.rm = TRUE))   
  data_tot_summary$time <- as.numeric(data_tot_summary$time)
  data_tot_summary <- data_tot_summary[which(data_tot_summary$time %in% (min(data_tot_summary$time)+19):(max(data_tot_summary$time)-20)), ]

    
  plot_trait = list()
  
  for (TRAIT in c(trait_dim, "age_Mature")){
    for (sp_plotted in sp){

      TRAIT_new <- new_trait_dim[grep(pattern=TRAIT, new_trait_dim)]

      data_tot_summary_tr    <- data_tot_summary[which(data_tot_summary$TRAIT==TRAIT),]
      data_tot_summary_tr_sp <- data_tot_summary_tr[which(data_tot_summary_tr$Species==sp_plotted),] 

      # reorder Fmax/Fstart factor for prettier plots
      data_tot_summary_tr_sp$Fish_press    <- factor(as.factor(data_tot_summary_tr_sp$Fish_press), levels = ideal_order_factor)

      
      if(nrow(data_tot_summary_tr_sp)==0){
        plot_trait[[TRAIT]][[sp_plotted]] <- ggplot(data = data_tot_summary_tr_sp, aes(x = time, group = Fish_press)) 
        next}
      
      reference <- subset(na.omit(data_tot_summary_tr_sp), Fish_press == paste0("Fmax=", F_vec[[sp]][1], " Fstart=", F_vec[[sp]][2]))


      x_axis <- "Time"
      #tr = ifelse(TRAIT == "fecundityByAge", "log egg number", TRAIT)
      title = paste0(sp_plotted, " ranges of ", TRAIT, "\nF varies for ", sp, ", reference is", F_vec[[sp]])
      plot_trait[[TRAIT]][[sp_plotted]] <- ggplot(data = data_tot_summary_tr_sp, aes(x = time, group = Fish_press)) +
        geom_ribbon(aes(ymin = mean-sd, ymax = mean + sd),
                    alpha = 0.5, fill = "grey70")+
        # 1) Layer: each panel gets its own group (via facet_wrap)
        geom_line(aes(y = mean)) +
        # 2) Layer: reference group A, drawn in all panels
        geom_ribbon(data = reference,
                    fill = "grey70", alpha = 0.5,
                    aes(x = time, ymin = mean-sd, ymax = mean+sd, group = 1)) +
        geom_line(data = reference,
                  aes(x = time, y = mean, group = 1),
                  linetype = "dashed", linewidth = 0.7) +
        theme_bw()+
        ylab(label = TRAIT_new)+
        xlab(label = x_axis)+
        labs(title = title)+
        facet_wrap(~ Fish_press)
    }
    
    pdf(file = paste0(pathout, "/plot_", TRAIT_new, "_Sp_fish_", sp, ".pdf"), width = 10, height = 8)
    par(mfrow=c(1,1))
    print(plot_trait[[TRAIT]])
    dev.off()
    
  } #trait plot
} #sp plot
}#function



# else{
#             if (!(TRAIT %in% c("mortality", "school", 'fecundity', "abundanceByAge", "fecundityByAge",
#                                "MpredbyY", "MstarbyY", "MnatbyY", "FbyY", "ZbyY", "MforbyY", "DisbyY"))){
#               data_t <- as.data.frame(data[,,i+1])
#             }
#             if (TRAIT == "mortality"){
#               data_t <- data.frame("time"= rownames(ns1$mortalityByAge[[species[1]]]$Mpred))
#               rownames(data_t) = data_t$time
#               names_list <- names(ns1$mortalityByAge[[species[1]]])
#               data_t[, species[1]] = rowSums(sapply(names_list, FUN = function(x) {
#                 if(x == "Mnat") rep(0, nrow(data_t)) else rowMeans(ns1$mortalityByAge[[species[1]]][[x]][,,(i+1)])
#               }))
#               for (spi in species[-1]){
#                 data_t[, spi] = rowSums(sapply(names_list, FUN = function(x) {
#                   if(x == "Mnat") rep(0, nrow(data_t)) else rowMeans(ns1$mortalityByAge[[spi]][[x]][,,(i+1)])
#                 }))
#               }
#             }
#             if (TRAIT == "eggNumberByAge"){
#               data_t_all <- read.csv(paste0(outputDir_Simu, "/AgeIndicators/ns_fecundityDistribByAge_Simu", i, ".csv"), skip=c(1))
#               data_t_t <- data_t_all[, which(colnames(data_t_all)%in% c("Time", species[1]))]
#               data_t   <- data.frame("time"= names(table(data_t_t$Time)))
#               data_t_t_t <- data_t_all %>%
#                 group_by(Time) %>%
#                 summarise(mean = mean(get(species[1]), na.rm = TRUE))
#               data_t[,species[1]] <- data_t_t_t$mean
#               for (spi in species[-1]){
#                 data_t_t <- data_t_all[, which(colnames(data_t_all)%in% c("Time", spi))]
#                 data_t_t_t <- data_t_all %>%
#                   group_by(Time) %>%
#                   summarise(mean = mean(get(spi), na.rm = TRUE))
#                 data_t[,spi] <- data_t_t_t$mean
#               }
#             }
#             if (TRAIT %in% c("MpredbyY", "MstarbyY", "MnatbyY", "FbyY", "ZbyY", "MforbyY", "DisbyY")) {
#               TRAIT_t <- stringr::str_remove(TRAIT, "byY")
#               data_t_all <- read.csv(paste0(outputDir_Simu, "/Mortality/ns_mortalityRateDistribByAge-", species[1], "_Simu", i, ".csv"), skip=c(1))
#               data_t_all <- data_t_all[, -ncol(data_t_all)] #remove useless column, also called Age
#               data_t_t <- data_t_all[, which(colnames(data_t_all)%in% c("Time", species[1]))]
#               data_t   <- data.frame("time"= unique(data_t_all$Time))
#               data_t_t_t <- data_t_all %>%
#                 group_by(Time) %>%
#                 summarise(mean = mean(get(TRAIT_t), na.rm = TRUE))
#               data_t[,species[1]] <- data_t_t_t$mean
#               for (spi in species[-1]){
#                 data_t_t <- read.csv(paste0(outputDir_Simu, "/Mortality/ns_mortalityRateDistribByAge-", spi, "_Simu", i, ".csv"), skip=c(1))
#                 data_t_t_t <- data_t_t %>%
#                   group_by(Time) %>%
#                   summarise(mean = mean(get(TRAIT_t), na.rm = TRUE))
#                 data_t[,spi] <- data_t_t_t$mean
#               }
#             }
#             if (TRAIT == "abundanceByAge"){
#               data_t_t <- data[[species[1]]]
#               data_t <- data.frame("time"= colnames(data_t_t[,,i+1]))
#               rownames(data_t) = data_t$time
#               data_t[,species[1]] <- colSums(data_t_t[,,i+1])
#               for (spi in species[-1]){
#                 data_t_t <- data[[spi]]
#                 data_t[,spi] <- colSums(data_t_t[,,i+1])
#               }
#             }
#             if (TRAIT %in% c("school", 'fecundity')) {
#               if (TRAIT == c("fecundity")) data_t_all <- read.csv(paste0(outputDir_Simu, "/AgeIndicators/ns_fecundityDistribByAge_Simu", c(i), ".csv"), skip=c(1))
#               if (TRAIT == c("school")) data_t_all <- read.csv(paste0(outputDir_Simu, "/ns_nschoolDistribByAge_Simu", c(i), ".csv"), skip=c(1))
#               data_t_sum <- t(sapply(names(table(data_t_all$Time)), FUN = function(x) {colSums(data_t_all[which(data_t_all$Time == x), ])}))
#               data_t <- data_t_sum[, -which(colnames(data_t_sum) %in% c("Time"))]
#               data_t <- as.data.frame(data_t[, which(colnames(data_t) %in% species)])
#             }

          # if (!grepl(TRAIT, pattern="byY")) data_t$time <- rownames(data_t)
          # data_t <- data_t%>%
          #     pivot_longer(!time, names_to = "Species", values_to = "value")
          # data_t$rep <- i
          # data_t$Fish_press = paste0("Fmax=", Fmax, " Fstart=", Fstart)
          # data_t$TRAIT = TRAIT



