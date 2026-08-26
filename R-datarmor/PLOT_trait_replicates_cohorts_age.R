
PLOT_trait_replicates2 <- function(SP, Time1=150, Time2=160, path_in, simu_name){

# setwd("X:/osmose-ns")
if (grepl(simu_name, pattern= "Chap3")) {
    main = "New_gpu_conf_RCP.R"
    term = "$"} 
if (grepl(simu_name, pattern= "Chap2")) {
    main = "New_gpu_conf.R"
    term = ".R"}

repo        = "Simu_LONG"
multip_name = "y"
multipliers = c(1,0.5,0.75,1.25,1.5)

trait_dim= c("GonadWeight", "meanWeightByAge", "abundanceMature", 
     "meanSizeByAge", "fecundityByAge",
     "Mpred", "Mstar", "Mnat", "F", "Mfor",
     "abundanceByAge", "biomassByAge", 
     "mortalityByAge", 
     "sizeMature", "eggNumberByAge")

trait_dim = c("meanWeightByAge", "enetByAge")

ideal_order_factor <- paste0("Fmax=", rep(c(0.5,0.75,1,1.25,1.5), 5),
                             " Fstart=", c(rep(0.5, 5),rep(0.75, 5),rep(1, 5),rep(1.25, 5),rep(1.5, 5)))


configDir=getwd()

list_d <- list.dirs(path = path_in, recursive = F, full.names = T)

# for the test : has to be erased when real launch
# nb_sp = 7
# species = species[1:7]
# sp_names = sp_names[1:7]
# species_i = species_i[1:7]
# years = 120
# ny       <- 120-get_par(conf, "output.start.years")
# age_vec <- seq(0,ny, by=1/dt)[-1]
# nb_rep = 1
# nage = years

pathout <- paste0("output/PLOTS/Trait_variables_by_cohort_by_age_withWbyAge")
pathout <- paste0("/home1/scratch/mbeneat/output/SIMU/", simu_name, "/Trait_variables_by_cohort_by_age_withWbyAge")
if (!dir.exists(pathout)) dir.create(pathout, recursive=T)

F_vec <- list()
#species = species[1:2]
spec_names <- species
cat(species, "\n\n\n")

for (sp in SP){  

  cat("\nProcessing species ", sp, "\n")

  # only two sp plotted out of the 15, to gain space
  sp_i <- which(spec_names == sp)-1
  species = spec_names[c(c(sp_i+1), ifelse(sp_i<14, sp_i+2, 1))]  
  sp_dir <- list_d[grep(pattern = paste0("sp", sp_i), list_d)]

  cat("sp_dir = ", sp_dir)
  
  # initialisation
  data_tot <- NULL
  plot_Neth_N_alpha_NT_sp     <- list()
  Fmax_1 = -99
  Fstart_1 = -99
  
  for (Fmax in multipliers){
    Fmax_dir <- sp_dir[grep(pattern = paste0("Fmax", Fmax, "_"), sp_dir)]
    for (Fstart in 1){ #multipliers){
      dir <- Fmax_dir[grep(pattern = paste0("Fstart", Fstart, "_sp", sp_i, term), Fmax_dir)]
     
      cat(dir)
 
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
      data_Ne_raw <- list()
      data_Ne_raw[["data_Ne"]] <- estimate_for_Ne(outputDir_Simu, outputDir_Ne, trait_dim, species, nb_rep, OUTPUT_data=ns1)
      data_Ne_raw[["data_Ne"]][[which(names(data_Ne_raw[["data_Ne"]]) == "cohort_time_step_corresponding")]] <- NULL
  
  
      reshape_data_Ne = lapply(data_Ne_raw, FUN = myf2)
      reshape_data_Ne = lapply(reshape_data_Ne, FUN = reshape_list) # have the name of the species before the trait name in the list order
      reshape_data_Ne = reshape_list(reshape_data_Ne)
  
      new_trait_dim <- trait_dim
      # bind replicate dimensions and remove year or cohort dimension
      for (rep in 0:c(nb_rep-1)){
       age_Mature = ns1[["ageMature"]]
       age_Mature = lapply(seq(dim(age_Mature)[3]), function(x) as.data.frame(age_Mature[ , , x]))

       for (nm in names(reshape_data_Ne[[sp]][["data_Ne"]])){
           # should we sum or average the cohorts dimensions ? 
          fct_i = ifelse(nm %in% c("biomassByAge", "schoolNumber", "abundanceMature", "eggNumberByAge", "abundanceByAge", "GonadWeight"), "mean", 
                   ifelse(nm %in% c("sizeMature", "mortalityByAge", "meanSizeByAge", "fecundityByAge", "Mpred", "Mstar", "Mnat", "F", "Z", "Mfor", "Dis", "meanWeightByAge", "enetByAge"), "mean", "total"))
          new_nm <- paste0(fct_i, "_of_", nm, "_per_age") 
          new_trait_dim[which(trait_dim == nm)] <- new_nm
          data_t_t <- reshape_data_Ne[[species[1]]][["data_Ne"]][[nm]]
          data_t_t <- data_t_t[,which(colnames(data_t_t) %in% paste0("cohort_", Time1:Time2)),]
          
          if (length(trait_dim)>2){ # average between cohorts if too any dimensions
            data_t <- data.frame(matrix(nrow=nrow(data_t_t)))
            data_t[,"age"] <- as.numeric(rownames(data_t_t[,,rep+1]))
            if (fct_i %in% c("sum", "mean")) data_t[,"value"] <- apply(data_t_t[,,rep+1], 1, function(x) if(fct_i == "sum") sum(x, na.rm=T) else mean(x, na.rm=T))
            data_t[,"value"] <- apply(data_t_t[,,rep+1], 1, fct_i)
            data_t[,"Species"] <- species[1]
            for (sp_j in species[-1]){
              data_t_t <- reshape_data_Ne[[sp_j]][["data_Ne"]][[nm]]
              data_t_t <- data_t_t[,which(colnames(data_t_t) %in% paste0("cohort_", Time1:Time2)),]
              data_tj <- data_t
              data_tj[,"age"] <- as.numeric(rownames(data_t_t[,,rep+1]))
              if (fct_i %in% c("sum", "mean")) data_tj[,"value"] <- apply(data_t_t[,,rep+1], 1, function(x) if(fct_i == "sum") sum(x, na.rm=T) else mean(x, na.rm=T))
              data_tj[,"Species"] <- sp_j
              data_t <- rbind(data_t, data_tj)
            }#species
          } else { # keep cohort dimension if few dimensions
            data_t <- data.frame(matrix(nrow=nrow(data_t_t)))
            data_t[,"age"] <- as.numeric(rownames(data_t_t[,,rep+1]))
            pivot_longer(data_t_t[,,rep+1], cols = colnames(data_t_t[,,rep+1]), names_to = "cohort", values_to = "value")
            data_t[,"Species"] <- species[1]
            for (sp_j in species[-1]){
              data_t_t <- reshape_data_Ne[[sp_j]][["data_Ne"]][[nm]]
              data_t_t <- data_t_t[,which(colnames(data_t_t) %in% paste0("cohort_", Time1:Time2)),]
              data_tj <- data_t
              data_tj[,"age"] <- as.numeric(rownames(data_t_t[,,rep+1]))
              pivot_longer(data_t_t[,,rep+1], cols = colnames(data_t_t[,,rep+1]), names_to = "cohort", values_to = "value")
              data_tj[,"Species"] <- sp_j
              data_t <- rbind(data_t, data_tj)
            }#species
          }
          
          data_t$TRAIT <- new_nm
          if (nm == names(reshape_data_Ne[[species[1]]][["data_Ne"]])[1]){
            data_rep <- data_t
          } else {
            data_rep <- rbind(data_rep, data_t)
          }
        }#trait
  
        data_t <- data_rep[,-1]
        data_t$replicate <- rep
        data_t$var_F_for_sp <- sp
        
        #Replicate : saving one replicate at a time instead of all of them together
        # if (rep == 0){
        #   dat <- data_t
        # } else {
        #   dat <- rbind(dat, data_t)
        # }
        
        pathoutputRDATA <-paste0(pathout, "/RDATA/Fmax", Fmax, "_Fstart", Fstart, "_sp", sp)
        if (!dir.exists(pathoutputRDATA)){dir.create(pathoutputRDATA, recursive = T)}
        save(data_t, file=paste0(pathoutputRDATA, "/data_traits_summary_Fmax_", Fmax, "_Fstart_", Fstart, "_replicate_", rep, ".RData"))
        
      }#replicate

      #Fmax & Fstart : saving one at a time instead of all together
      # dat$Fish_press <- paste0("Fmax=", Fmax, " Fstart=", Fstart)
      # if (Fmax == Fmax_1 & Fstart == Fstart_1) {
      #   data_complete <- dat
      # } else { data_complete <- rbind(data_complete, dat) }
  
    }#Fstart
  }#Fmax


  
  #save dataset with Fmax and Fstart and replicates varying
  # pathoutputRDATA <-paste0(pathout, "/RDATA/Fmax", Fmax, "_Fstart", Fstart, "_sp", sp)
  # if (!dir.exists(pathoutputRDATA)){dir.create(pathoutputRDATA, recursive = T)}
  # if (!file.exists(paste0(pathoutputRDATA, "/data_traits_summary.RData"))){
  #   save(data_complete, file=paste0(pathoutputRDATA, "/data_traits_summary.RData"))
  # } else {
  #   load(file=paste0(pathoutputRDATA, "/data_traits_summary.RData"))
  #   data_save <- rbind(loaded_data, data_complete)
  #   save(data_save, file=paste0(pathoutputRDATA, "/data_traits_summary.RData"))
  # }
    
  
  #PLOT the mean and the std deviation
  # data_tot_summary <- data_complete %>%
  #   group_by(age, Species, Fish_press, TRAIT) %>%
  #   summarise(mean = mean(value, na.rm = TRUE),
  #             sd   = sd(value, na.rm = TRUE))
  # data_tot_summary$age <- as.numeric(data_tot_summary$age)
  # 
  # plot_trait = list()
  # for (TRAIT in new_trait_dim){
  #   for (sp_plotted in sp){
  # 
  #     data_tot_summary_tr    <- data_tot_summary[which(data_tot_summary$TRAIT==TRAIT),]
  #     data_tot_summary_tr_sp <- data_tot_summary_tr[which(data_tot_summary_tr$Species==sp_plotted),] 
  # 
  #     # reorder Fmax/Fstart factor for prettier plots
  #     data_tot_summary_tr_sp$Fish_press    <- factor(as.factor(data_tot_summary_tr_sp$Fish_press), levels = ideal_order_factor)
  # 
  #     
  #     if(nrow(data_tot_summary_tr_sp)==0){
  #       plot_trait[[TRAIT]][[sp_plotted]] <- ggplot(data = data_tot_summary_tr_sp, aes(x = age, group = Fish_press)) 
  #       next}
  #     
  #     reference <- subset(na.omit(data_tot_summary_tr_sp), Fish_press == paste0("Fmax=", F_vec[[sp]][1], " Fstart=", F_vec[[sp]][2]))
  # 
  # 
  #     x_axis <- "Age"
  #     #tr = ifelse(grepl(TRAIT, pattern= "fecundityByAge"), "log egg number", TRAIT)
  #     title = paste0("By cohort ", sp_plotted, " ranges of ", TRAIT, "\nF varies for", sp, ", reference is", F_vec[[sp]])
  #     plot_trait[[TRAIT]][[sp_plotted]] <- ggplot(data = data_tot_summary_tr_sp, aes(x = age, group = Fish_press)) +
  #       geom_ribbon(aes(ymin = mean-sd, ymax = mean + sd),
  #                   alpha = 0.5, fill = "grey70")+
  #       # 1) Layer: each panel gets its own group (via facet_wrap)
  #       geom_line(aes(y = mean)) +
  #       # 2) Layer: reference group A, drawn in all panels
  #       geom_ribbon(data = reference,
  #                   fill = "grey70", alpha = 0.5,
  #                   aes(x = age, ymin = mean-sd, ymax = mean+sd, group = 1)) +
  #       geom_line(data = reference,
  #                 aes(x = age, y = mean, group = 1),
  #                 linetype = "dashed", linewidth = 0.7) +
  #       theme_bw()+
  #       ylab(label = TRAIT)+
  #       xlab(label = x_axis)+
  #       labs(title = title)+
  #       facet_wrap(~ Fish_press)
  #   }
  #   
  #   pdf(file = paste0(pathout, "/plot_", TRAIT, "_Sp_fish_", sp, ".pdf"), width = 10, height = 8)
  #   par(mfrow=c(1,1))
  #   print(plot_trait[[TRAIT]])
  #   dev.off()
  #   
  # } #trait plot
} #sp plot
}#function







