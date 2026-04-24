

  
# initialisation
main = "New_gpu_conf.R"
source("osmose-ns-v4.4.1_clim/R/functions/source.R")
base_dir <- "osmose-ns-v4.4.1_clim/output/PLOTS/FINAL_Ne_Htz/INIT"
folders  <- list.dirs(base_dir, recursive = F)
multipliers = c(1,0.5,0.75,1.25,1.5)
pathout <- paste0(configDir, "/output/PLOTS/Trait_variables_by_cohort_by_age2")
F_vec <- list()
spec_names <- species
ideal_order_factor <- paste0("Fmax=", rep(c(0.5,0.75,1,1.25,1.5), 5),
                             " Fstart=", c(rep(0.5, 5),rep(0.75, 5),rep(1, 5),rep(1.25, 5),rep(1.5, 5)))

dimensions <- c("survivorship_curve", "Vx")

for (sp in spec_names){  
  
  cat("\nProcessing species ", sp, "\n")
  
  # only two sp plotted out of the 15, to gain space
  sp_i <- which(spec_names == sp)-1
  species = spec_names[c(c(sp_i+1), ifelse(sp_i<14, sp_i+2, 1))]  
  sp_dir <- folders[grep(pattern = paste0("sp", sp), folders)]
  
  # initialisation
  data_tot <- NULL
  plot_Neth_N_alpha_NT_sp     <- list()
  Fmax_1 = -99
  Fstart_1 = -99
  
  for (Fmax in multipliers){
    Fmax_dir <- sp_dir[grep(pattern = paste0("Fmax", Fmax, "_"), sp_dir)]
    for (Fstart in multipliers){
      dir <- Fmax_dir[grep(pattern = paste0("Fstart", Fstart, "_sp", sp), Fmax_dir)]
      
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
      #cat("\n\t---> folder Fmax", Fmax, " Fstart", Fstart, " sp", sp, "\n")
      
      data_Ne_raw <- get(load(paste0(simulation, "/INIT2.RData")))
      data_Ne <- data_Ne_raw$data_Ne[[sp]]
      data_Ne <- data_Ne[which(names(data_Ne) %in% dimensions)]
      

      for (rep in 1:c(nb_rep)){
        for (nm in names(data_Ne)){
          #  fct_i = ifelse(nm %in% c("mortalityByAge", "abundanceByAge"), sum, mean)
          fct_i = ifelse(nm %in% c("biomassByAge", "schoolNumber", "abundanceMature", "eggNumberByAge", "abundanceByAge"), "sum", 
                         ifelse(nm %in% c(dimensions, "sizeMature", "mortalityByAge", "meanSizeByAge", "fecundityByAge", "Mpred", "Mstar", "Mnat", "F", "Z", "Mfor", "Dis"), "mean", "total"))
          data_t_t <- data_Ne[[nm]][[rep]]
          data_t_t <- data_t_t[,which(colnames(data_t_t) %in% paste0("cohort_", 100:180))]
          data_t <- data.frame(matrix(nrow=nrow(data_t_t)))
          data_t[,"age"] <- (rownames(data_t_t))
          if (fct_i %in% c("sum", "mean")) data_t[,"value"] <- apply(data_t_t, 1, function(x) if(fct_i == "sum") sum(x, na.rm=T) else mean(x, na.rm=T))
          data_t[,"value"] <- apply(data_t_t[,], 1, fct_i)
          data_t[,"Species"] <- species[1]
          data_t$TRAIT <- paste0(fct_i, "_",nm)
          if (nm == names(data_Ne)[1]){
            data_rep <- data_t
          } else {
            data_rep <- rbind(data_rep, data_t)
          }
        }
        data_t <- data_rep[,-1]
        data_t$replicate <- rep
        data_t$var_F_for_sp <- sp
        if (rep == 1){
          dat <- data_t
        } else {
          dat <- rbind(dat, data_t)
        }
      }
      dat$Fish_press <- paste0("Fmax=", Fmax, " Fstart=", Fstart)
      
      if (Fmax == Fmax_1 & Fstart == Fstart_1) {
        data_complete <- dat
      } else { data_complete <- rbind(data_complete, dat) }
      
    }#Fstart
  }#Fmax
  
  plot_trait = list()
  
  for (TRAIT in dimensions){
    for (sp_plotted in sp){
      
      data_tot_summary_tr    <- data_complete[grep(data_complete$TRAIT, pattern=TRAIT),]
      data_tot_summary_tr_sp <- data_tot_summary_tr[which(data_tot_summary_tr$Species==sp_plotted),] 
      data_tot_summary_tr_sp <- data_tot_summary_tr_sp %>%
        group_by(age, Species, TRAIT, Fish_press) %>%
        summarise(mean=mean(value, na.rm=T),
                  sd=sd(value, na.rm=T))
      data_tot_summary_tr_sp$mean <- as.numeric(data_tot_summary_tr_sp$mean)
      data_tot_summary_tr_sp$sd <- as.numeric(data_tot_summary_tr_sp$sd)
      data_tot_summary_tr_sp$age <- as.numeric(data_tot_summary_tr_sp$age)
      # reorder Fmax/Fstart factor for prettier plots
      data_tot_summary_tr_sp$Fish_press    <- factor(as.factor(data_tot_summary_tr_sp$Fish_press), levels = ideal_order_factor)
      
      
      if(nrow(data_tot_summary_tr_sp)==0){
        plot_trait[[TRAIT]][[sp_plotted]] <- ggplot(data = data_tot_summary_tr_sp, aes(x = age, group = Fish_press)) 
        next}
      
      reference <- subset(na.omit(data_tot_summary_tr_sp), Fish_press == paste0("Fmax=", F_vec[[sp]][1], " Fstart=", F_vec[[sp]][2]))
      
      
      x_axis <- "Age"
      #tr = ifelse(grepl(TRAIT, pattern= "fecundityByAge"), "log egg number", TRAIT)
      title = paste0("By cohort ", sp_plotted, " ranges of ", TRAIT, "\nF varies for", sp, ", reference is", F_vec[[sp]])
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
  
  if (sp == spec_names[1])  {
    data_complete_all_sp <- data_complete
  } else  data_complete_all_sp <- rbind(data_complete_all_sp, data_complete)
  
} #sp plot

save(data_complete_all_sp, file=paste0(configDir, "/output/RData/Vkx_surviv.RData"))





