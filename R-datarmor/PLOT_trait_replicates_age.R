

PLOT_trait_replicates2 <- function(SP){

# setwd("X:/osmose-ns")

repo        = "Simu_LONG"
multip_name = "y"
multipliers = c(1,0.5,0.75,1.25,1.5)


trait_dim= c("biomassByAge", "sizeMature", "abundanceMature", "fecundityByAge", "meanSizeByAge", 
     "MpredbyY", "MstarbyY", "MnatbyY", "FbyY", "ZbyY", "MforbyY", "DisbyY", "eggNumberByAge", 
     "abundanceByAge", "mortality", "school")
#trait_dim=c("abundanceMature", "sizeMature")

ideal_order_factor <- paste0("Fmax=", rep(c(0.5,0.75,1,1.25,1.5), 5),
                             " Fstart=", c(rep(0.5, 5),rep(0.75, 5),rep(1, 5),rep(1.25, 5),rep(1.5, 5)))



path_in <- c("C:/Users/mbeneat/Documents/osmose/Ne_on_NT_with_F/osmose-ns/osmose-ns-v4.4.1_clim/output/test")
# path_in <- c("Y:/outputs-meso-marine")
path_in <- c("/home1/datawork/mbeneat/SIMU_Ne_on_NT_with_F_FIN")
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

pathout <- paste0(configDir, "/output/PLOTS/Trait_variables_by_year_by_age150-160")

plot_Neth_Nalpha    <- list()
F_vec     <- list()

for (sp in SP){
  data_all <- NULL
  sp_i <- which(species == sp)-1
  # open one dir at a time
  sp_dir <- list_d[grep(pattern = paste0("sp", sp_i), list_d)]
  
  plot_Neth_N_alpha_NT_sp     <- list()
  Fmax_1 = -99
  Fstart_1 = -99
  
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


      age_Mature = ns1[["ageMature"]]
      age_Mature = lapply(seq(dim(age_Mature)[3]), function(x) as.data.frame(age_Mature[ , , x]))
      for (TRAIT in trait_dim[1]){
        cat("Processing trait:", TRAIT, "\n")
        data_tot = NULL
        data <- ns1[[TRAIT]]
        
        for (i in 0:c(nb_rep-1)){
          if (TRAIT == "mortality"){
            data_tot <- data.frame("Age"= colnames(ns1$mortalityByAge[[species[1]]]$Mpred))
            rownames(data_tot) = data_tot$Age
            names_list <- names(ns1$mortalityByAge[[species[1]]])
            data_tot[, species[1]] = rowSums(sapply(names_list, FUN = function(x) {
              dat <- ns1$mortalityByAge[[species[1]]][[x]][,,(i+1)]
              dat <- dat[which(as.numeric(rownames(dat)) %in% 150:160),]
              if(x == "Mnat") rep(0, nrow(data_tot)) else colMeans(dat)
            }))
            for (spi in species[-1]){
              data_tot[, spi] = rowSums(sapply(names_list, FUN = function(x) {
                dat <- ns1$mortalityByAge[[spi]][[x]][,,(i+1)]
                dat <- dat[which(as.numeric(rownames(dat)) %in% 150:160),]
                if(x == "Mnat") rep(0, nrow(data_tot)) else colMeans(dat)
              }))
            }
          }
          if (TRAIT %in% c("abundanceMature", "sizeMature")){
            data_tot <- data.frame("Age"=0:30)
            rownames(data_tot) = data_tot$Age
            age_mat_sp_rep = age_Mature[[i+1]][, species[1]]
            if (TRAIT == "abundanceMature") vec_y = ns1$abundanceByAge[[species[1]]][,,i+1]
            if (TRAIT == "sizeMature")      vec_y = ns1$meanSizeByAge[[species[1]]][,,i+1]
            for (y in 1:length(age_mat_sp_rep)){
                if (TRAIT == "abundanceMature") vec_y[,which(as.numeric(colnames(ns1$abundanceByAge[[species[1]]]))< ceiling(age_mat_sp_rep[y]))] <- 0
                if (TRAIT == "sizeMature")      vec_y[,which(as.numeric(colnames(ns1$meanSizeByAge[[species[1]]]))< ceiling(age_mat_sp_rep[y]))] <- 0
            }
            vec_y <- vec_y[which(as.numeric(rownames(vec_y)) %in% 150:160),]
            if (TRAIT == "abundanceMature") data_tot[,species[1]] <- c(colMeans(vec_y), rep(0, 31 - ncol(vec_y)))
            if (TRAIT == "sizeMature")      data_tot[,species[1]] <- c(colMeans(vec_y), rep(0, 31 - ncol(vec_y)))
            for (spi in species[-1]){
              age_mat_sp_rep = age_Mature[[i+1]][, spi]
              if (TRAIT == "abundanceMature") vec_y = ns1$abundanceByAge[[spi]][,,i+1]
              if (TRAIT == "sizeMature")      vec_y = ns1$meanSizeByAge[[spi]][,,i+1]
              for (y in 1:length(age_mat_sp_rep)){
                if (TRAIT == "abundanceMature") vec_y[,which(as.numeric(colnames(ns1$abundanceByAge[[spi]]))< ceiling(age_mat_sp_rep[y]))] <- 0
                if (TRAIT == "sizeMature")      vec_y[,which(as.numeric(colnames(ns1$meanSizeByAge[[spi]]))< ceiling(age_mat_sp_rep[y]))] <- 0
              }
              vec_y <- vec_y[which(as.numeric(rownames(vec_y)) %in% 150:160),]
              if (TRAIT == "abundanceMature") data_tot[,spi] <- c(colMeans(vec_y), rep(0, 31 - ncol(vec_y)))
              if (TRAIT == "sizeMature")      data_tot[,spi] <- c(colMeans(vec_y), rep(0, 31 - ncol(vec_y)))
            }

          }
          if (TRAIT %in% c("abundanceByAge", "meanSizeByAge")){
            data_tot_t <- data[[species[1]]]
            data_tot <- data.frame("Age"= 0:30)
            rownames(data_tot) = data_tot$Age
            if (TRAIT == ("abundanceByAge")) data_tot[1:length(colSums(data_tot_t[which(rownames(data_tot_t) %in% 150:160),,i+1])),species[1]] <- colSums(data_tot_t[which(rownames(data_tot_t) %in% 150:160),,i+1],   na.rm = TRUE)
            if (TRAIT == ("meanSizeByAge")) data_tot[1:length(colMeans(data_tot_t[which(rownames(data_tot_t) %in% 150:160),,i+1])),species[1]] <- colMeans(data_tot_t[which(rownames(data_tot_t) %in% 150:160),,i+1],   na.rm = TRUE)
            for (spi in species[-1]){
              data_tot_t <- data[[spi]]
              if (TRAIT == ("abundanceByAge")) data_tot[1:length(colSums(data_tot_t[which(rownames(data_tot_t) %in% 150:160),,i+1])),spi] <- colSums(data_tot_t[which(rownames(data_tot_t) %in% 150:160),,i+1],   na.rm = TRUE)
              if (TRAIT == ("meanSizeByAge"))  data_tot[1:length(colMeans(data_tot_t[which(rownames(data_tot_t) %in% 150:160),,i+1])),spi] <- colMeans(data_tot_t[which(rownames(data_tot_t) %in% 150:160),,i+1],  na.rm = TRUE)
            }
          }
          if (TRAIT %in% c("MpredbyY", "MstarbyY", "MnatbyY", "FbyY", "ZbyY", "MforbyY", "DisbyY")) {
            TRAIT_t <- stringr::str_remove(TRAIT, "byY")
            data_tot_all <- read.csv(paste0(outputDir_Simu, "/Mortality/ns_mortalityRateDistribByAge-", species[1], "_Simu", i, ".csv"), skip=c(1))
            data_tot_all <- data_tot_all[, -ncol(data_tot_all)] #remove useless column, also called Age
            data_tot_t <- data_tot_all[, which(colnames(data_tot_all)%in% c("Age", species[1]))]
            data_tot   <- data.frame("Age"= unique(data_tot_all$Age))
            data_tot_all <- data_tot_all[which(data_tot_all$Time %in% 150:160),]
            data_tot_t_t <- data_tot_all %>%
              group_by(Age) %>%
              summarise(mean = mean(get(TRAIT_t), na.rm = TRUE))
            data_tot[,species[1]] <- as.data.frame(data_tot_t_t$mean)
            for (spi in species[-1]){
              data_tot_t <- read.csv(paste0(outputDir_Simu, "/Mortality/ns_mortalityRateDistribByAge-", spi, "_Simu", i, ".csv"), skip=c(1))
              data_tot_t_t <- data_tot_t %>%
                group_by(Age) %>%
                summarise(mean = mean(get(TRAIT_t), na.rm = TRUE))
              data_tot[,spi] <- data_tot_t_t$mean
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
            fecundity_data <- fecundity_data[which(fecundity_data$Time %in% 150:160),]
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
              fecundity_data <- fecundity_data[which(fecundity_data$Time %in% 150:160),]
              data_tot <- full_join(data_tot, fecundity_data[, which(colnames(fecundity_data) %in% c("Age", "Time", "fecundity"))],
                                   by = c("Age", "Time"))
              colnames(data_tot)[ncol(data_tot)] <- spi
            }
            data_tot <- data_tot %>%
              group_by(Age) %>%
              summarise(across(all_of(species), ~ mean(.x, na.rm = TRUE)))
          }
          if (TRAIT %in% c("school", 'eggNumberByAge')) {
            if (TRAIT == c("eggNumberByAge")) data_tot_all <- read.csv(paste0(outputDir_Simu, "/AgeIndicators/ns_fecundityDistribByAge_Simu", c(i), ".csv"), skip=c(1))
            if (TRAIT == c("school")) data_tot_all <- read.csv(paste0(outputDir_Simu, "/ns_nschoolDistribByAge_Simu", c(i), ".csv"), skip=c(1))
            data_tot_all <- data_tot_all[which(data_tot_all$Time %in% 150:160),]
            data_tot_sum <- data_tot_all %>%
              group_by(Age) %>%
              summarise(across(all_of(species), ~ sum(.x, na.rm = TRUE))) %>%
              as.data.frame()
            data_tot <- as.data.frame(data_tot_sum[, which(colnames(data_tot_sum) %in% c(species, "Age"))])
          }
          if (!(TRAIT %in% c("mortality", "school", 'eggNumber', 'abundanceByAge', "eggNumberByAge",
                              "MpredbyY", "MstarbyY", "MnatbyY", "FbyY", "ZbyY", "MforbyY", "DisbyY", 
                              "abundanceMature", "sizeMature", "meanSizeByAge", "fecundityByAge"))){
            if (is.list(data)) {data_tot <- sapply(species, FUN = function(x) colSums(as.data.frame(data[[x]][,,i+1])))
              data_tot <- as.data.frame(data_tot)
            } else {data_tot <- as.data.frame(data[,,i+1])}
            data_tot$Age <- 0:30
            data_tot <- as.data.frame(data_tot[, which(colnames(data_tot) %in% c(species, "Age"))])
          }
          data_tot$age <- data_tot$Age
          data_tot$Age <- NULL
          data_tot <- data_tot%>%
            pivot_longer(!age, names_to = "Species", values_to = "value")
          data_tot$rep <- i
          data_tot$Fish_press = paste0("Fmax=", Fmax, " Fstart=", Fstart)
          data_tot$TRAIT = TRAIT
        
          if(TRAIT == trait_dim[1] && i==0){
            data_all <- data_tot
          } else {
            data_all <- rbind(data_all, data_tot)
          }

        }#rep
      }#TRAIT
      if (Fmax == Fmax_1 & Fstart == Fstart_1) {
        data_complete <- data_all
      }else {
        data_complete <- rbind(data_complete, data_all)
      }
      # p(paste0("\n\t--->done: folder Fmax", Fmax, " Fstart", Fstart, " sp", sp, "\n"))
    }#Fstart
  }#Fmax
  
  # PLOT the mean and the std deviation 
     # adding the sum of mortalities by age 
  
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
    ifelse_name <- ifelse(nm %in% c("school", "biomassByAge", "eggNumberByAge", "abundanceByAge", "abundanceMature"), "sum", 
                   ifelse(nm %in% c("sizeMature", "mortality", "meanSizeByAge", "fecundityByAge", "MpredbyY", "MstarbyY", "MnatbyY", "FbyY", "ZbyY", "MforbyY", "DisbyY"), "mean", 
                  #  ifelse(nm %in% c(), "age_mat", 
                  #  ifelse(nm %in% c(), "age_1", "total")
                   "total"))
    new_nm <- paste0(ifelse_name, "_of_", nm, "_per_age")
    new_trait_dim[which(trait_dim == nm)] <- new_nm
  }

  # if (all(c("MpredbyY", "MstarbyY", "MnatbyY", "FbyY", "ZbyY", "MforbyY", "DisbyY") %in% trait_dim)){
  #   data_tempo_mortbyage <-
  #     data_complete[which(data_complete$TRAIT %in% c("MpredbyY", "MstarbyY", "MnatbyY", "FbyY", "ZbyY", "MforbyY", "DisbyY")),]
  #   data_tempo_mortbyage <- data_tempo_mortbyage %>%
  #                 group_by(Fish_press, Species, age, rep) %>%
  #                 summarise(sum = sum(value, na.rm = TRUE))
  #   data_tempo_mortbyage$TRAIT = c("mortalityByY")
  #   data_complete <- rbind(data_complete,
  #                              data.frame("age"= data_tempo_mortbyage$age,
  #                                         "Species"= data_tempo_mortbyage$Species,
  #                                         "value"= data_tempo_mortbyage$sum,
  #                                         "rep"= data_tempo_mortbyage$rep,
  #                                         "Fish_press"= data_tempo_mortbyage$Fish_press,
  #                                         "TRAIT"= data_tempo_mortbyage$TRAIT))
  #   trait_dim = c(trait_dim, "mortalityByAge")
  # }
  #data_complete$value[which(data_complete$TRAIT == "fecundityByAge")] <- log(data_complete$value[which(data_complete$TRAIT == "fecundityByAge")]+10^{-8}) # easier to plot
  data_tot_summary <- data_complete %>%
    group_by(age, Species, Fish_press, TRAIT) %>%
    summarise(mean = mean(value, na.rm = TRUE),
              sd   = sd(value, na.rm = TRUE))   
  data_tot_summary$age <- as.numeric(data_tot_summary$age)
  data_tot_summary <- data_tot_summary[which(data_tot_summary$age %in% c(1:30)), ]
    
  plot_trait = list()
  
  for (TRAIT in unique(data_tot_summary$TRAIT)){      
    data_tot_summary_tr    <- data_tot_summary[which(data_tot_summary$TRAIT==TRAIT),]
    for (sp_plotted in sp){

      #TRAIT_new <- new_trait_dim[grep(pattern=TRAIT, new_trait_dim)]


      data_tot_summary_tr_sp <- data_tot_summary_tr[which(data_tot_summary_tr$Species==sp_plotted),] 

      # reorder Fmax/Fstart factor for prettier plots
      data_tot_summary_tr_sp$Fish_press    <- factor(as.factor(data_tot_summary_tr_sp$Fish_press), levels = ideal_order_factor)

      
      if(nrow(data_tot_summary_tr_sp)==0){
        plot_trait[[TRAIT]][[sp_plotted]] <- ggplot(data = data_tot_summary_tr_sp, aes(x = age, group = Fish_press)) 
        next}
      
      reference <- subset(na.omit(data_tot_summary_tr_sp), Fish_press == paste0("Fmax=", F_vec[[sp]][1], " Fstart=", F_vec[[sp]][2]))


      x_axis <- ifelse(TRAIT %in% c("abundanceByAge", "fecundityByAge"),
                        "Age", "Time")
      #tr = ifelse(TRAIT == "fecundityByAge", "log egg number", TRAIT)
      title = paste0(sp_plotted, " ranges of ", TRAIT, "\nF varies for", sp, ", reference is", F_vec[[sp]])
      plot_trait[[TRAIT]][[sp_plotted]] <- ggplot(data = data_tot_summary_tr_sp, aes(x = age, group = Fish_press)) +
        geom_ribbon(aes(ymin = mean-sd, ymax = mean + sd),
                    alpha = 0.5, fill = "grey70")+
        # 1) Layer: each panel gets its own group (via facet_wrap)
        geom_line(aes(y = mean)) +
        # 2) Layer: reference group A, drawn in all panels
        geom_ribbon(data = reference,
                    fill = "grey70", alpha = 0.5,
                    aes(x = age, ymin = mean-sd, ymax = mean+sd, group = 1)) +
        geom_line(data = reference,
                  aes(x = age, y = mean, group = 1),
                  linetype = "dashed", linewidth = 0.7) +
        theme_bw()+
        ylab(label = TRAIT)+
        xlab(label = x_axis)+
        labs(title = title)+
        facet_wrap(~ Fish_press)
    }
    
    pdf(file = paste0(pathout, "/plot_", TRAIT, "_Sp_fish_", sp, ".pdf"), width = 10, height = 8)
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
#               data_t <- data.frame("age"= colnames(ns1$mortalityByAge[[species[1]]]$Mpred))
#               rownames(data_t) = data_t$age
#               names_list <- names(ns1$mortalityByAge[[species[1]]])
#               data_t[, species[1]] = rowSums(sapply(names_list, FUN = function(x) {
#                 if(x == "Mnat") rep(0, nrow(data_t)) else colMeans(ns1$mortalityByAge[[species[1]]][[x]][,,(i+1)])
#               }))
#               for (spi in species[-1]){
#                 data_t[, spi] = rowSums(sapply(names_list, FUN = function(x) {
#                   if(x == "Mnat") rep(0, nrow(data_t)) else colMeans(ns1$mortalityByAge[[spi]][[x]][,,(i+1)])
#                 }))
#               }
#             }
#             if (TRAIT %in% c("MpredbyY", "MstarbyY", "MnatbyY", "FbyY", "ZbyY", "MforbyY", "DisbyY")) {
#               TRAIT_t <- stringr::str_remove(TRAIT, "byY")
#               data_t_all <- read.csv(paste0(outputDir_Simu, "/Mortality/ns_mortalityRateDistribByAge-", species[1], "_Simu", i, ".csv"), skip=c(1))
#               data_t_all <- data_t_all[, -ncol(data_t_all)] #remove useless column, also called Age
#               data_t_t <- data_t_all[, which(colnames(data_t_all)%in% c("Age", species[1]))]
#               data_t   <- data.frame("age"= unique(data_t_all$Age))
#               data_t_t_t <- data_t_all %>%
#                 group_by(Age) %>%
#                 summarise(mean = mean(get(TRAIT_t), na.rm = TRUE))
#               data_t[,species[1]] <- data_t_t_t$mean
#               for (spi in species[-1]){
#                 data_t_t <- read.csv(paste0(outputDir_Simu, "/Mortality/ns_mortalityRateDistribByAge-", spi, "_Simu", i, ".csv"), skip=c(1))
#                 data_t_t_t <- data_t_t %>%
#                   group_by(Age) %>%
#                   summarise(mean = mean(get(TRAIT_t), na.rm = TRUE))
#                 data_t[,spi] <- data_t_t_t$mean
#               }
#             }
#             if (TRAIT == "abundanceByAge"){
#               data_t_t <- data[[species[1]]]
#               data_t <- data.frame("age"= colnames(data_t_t[,,i+1]))
#               rownames(data_t) = data_t$age
#               data_t[,species[1]] <- colSums(data_t_t[,,i+1])
#               for (spi in species[-1]){
#                 data_t_t <- data[[spi]]
#                 data_t[,spi] <- colSums(data_t_t[,,i+1])
#               }
#             }
#             if (TRAIT %in% c("school", 'fecundity')) {
#               if (TRAIT == c("fecundity")) data_t_all <- read.csv(paste0(outputDir_Simu, "/AgeIndicators/ns_fecundityDistribByAge_Simu", c(i-1), ".csv"), skip=c(1))
#               if (TRAIT == c("school")) data_t_all <- read.csv(paste0(outputDir_Simu, "/ns_nschoolDistribByAge_Simu", c(i-1), ".csv"), skip=c(1))
#               data_t_sum <- t(sapply(names(table(data_t_all$Age)), FUN = function(x) {colSums(data_t_all[which(data_t_all$Age == x), ])}))
#               data_t <- data_t_sum[, -which(colnames(data_t_sum) %in% c("Time", "Age"))]
#               data_t <- as.data.frame(data_t[, which(colnames(data_t) %in% species)])
#             }
#             if (!grepl(TRAIT, pattern="byY")) data_t$age <- rownames(data_t)
#             data_t <- data_t%>%
#                 pivot_longer(!age, names_to = "Species", values_to = "value")
#             data_t$rep <- i
#             data_t$Fish_press = paste0("Fmax=", Fmax, " Fstart=", Fstart)
#             data_t$TRAIT = TRAIT
