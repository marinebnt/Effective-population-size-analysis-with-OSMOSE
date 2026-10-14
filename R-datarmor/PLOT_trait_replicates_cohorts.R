
# Plot trait trajectories by cohort for a given species (and a neighbor species for context)
# across fishing pressure scenarios (Fmax/Fstart). Saves a summary RData and PDF plots.
PLOT_trait_replicates2 <- function(SP){

# setwd("X:/osmose-ns")
# Configuration file used by the current run (kept for reference)
main = "New_gpu_conf.R"

# Multipliers for the fishing scenarios (currently fixed further below)
repo        = "Simu_LONG"
multip_name = "y"
multipliers = c(1,0.5,0.75,1.25,1.5)

# Trait list to extract from the simulation outputs
trait_dim= c("sizeMature", "abundanceMature", "meanSizeByAge", "fecundityByAge",
     "Mpred", "Mstar", "Mnat", "F", "Mfor",
     "abundanceByAge", "biomassByAge", "mortalityByAge", "schoolNumber", "eggNumberByAge")



# Ordering of Fmax/Fstart in facets for consistent plotting
ideal_order_factor <- paste0("Fmax=", rep(c(0.5,0.75,1,1.25,1.5), 5),
                             " Fstart=", c(rep(0.5, 5),rep(0.75, 5),rep(1, 5),rep(1.25, 5),rep(1.5, 5)))



# Input paths (last assignment is the active one)
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
# cohort_vec <- seq(0,ny, by=1/dt)[-1]
# nb_rep = 1
# ncohort = years

# Output folder for plots and summary data
pathout <- paste0(configDir, "/output/PLOTS/Trait_variables_by_cohort")

# Store the reference F scenario per species (used for dashed line)
F_vec <- list()
# Species list from the current configuration
spec_names <- species


for (sp in SP){
  cat("\nProcessing species ", sp, "\n")

  # Plot only target species + one additional species (to reduce plot size)
  sp_i <- which(spec_names == sp)-1
  species = spec_names[c(c(sp_i+1), ifelse(sp_i<14, sp_i+2, 1))]  
  # Locate directories for the focal species
  sp_dir <- list_d[grep(pattern = paste0("sp", sp_i), list_d)]

  # Initialization
  data_tot <- NULL
  plot_Neth_N_alpha_NT_sp     <- list()
  Fmax_1 = -99
  Fstart_1 = -99
  
  # Default/reference values (overwritten by loop below)
  Fmax=1.5
  Fstart = 1
  
  
   # Loop through fishing pressure scenarios (currently fixed to 0.5)
   for (Fmax in multipliers){
     Fmax_dir <- sp_dir[grep(pattern = paste0("Fmax", Fmax, "_"), sp_dir)]
     for (Fstart in multipliers){
       dir <- Fmax_dir[grep(pattern = paste0("Fstart", Fstart, "_sp", sp_i, ".R"), Fmax_dir)]
       
       # Skip if directory is empty
       if (length(list.files(dir))<1){
         cat("\ndata missing", dir, "\n")
         next}
       # First scenario becomes the reference for dashed line plotting
       if (Fmax_1==-99 && Fstart_1 ==-99){
         Fmax_1 = Fmax
         Fstart_1 = Fstart
         F_vec[[sp]] <- c(Fmax_1, Fstart_1)
       }
       
      # Read simulation outputs and compute Ne-related traits
       simulation = dir
       cat("\n\t---> folder Fmax", Fmax, " Fstart", Fstart, " sp", sp, "\n")
       
       ns1 = NULL
       outputDir_Simu  = simulation
       outputDir_Ne    = file.path(configDir, "output", "outputdirNe")
       #if(!dir.exists(outputDir_Ne)) {dir.create(outputDir_Ne, recursive = T)}
       ns1 = read_osmose(path = paste0(outputDir_Simu, "/"))
       data_Ne_raw <- list()
       data_Ne_raw[["data_Ne"]] <- estimate_for_Ne(outputDir_Simu, outputDir_Ne, dim_vec=trait_dim, species, nb_rep, OUTPUT_data=ns1)
       data_Ne_raw[["data_Ne"]][[which(names(data_Ne_raw[["data_Ne"]]) == "cohort_time_step_corresponding")]] <- NULL
   
   
      # Reshape outputs into list structure: species -> trait -> array
       reshape_data_Ne = lapply(data_Ne_raw, FUN = myf2)
       reshape_data_Ne = lapply(reshape_data_Ne, FUN = reshape_list) # have the name of the species before the trait name in the list order
       reshape_data_Ne = reshape_list(reshape_data_Ne)

       # Aggregate per cohort and species for each replicate
       new_trait_dim <- trait_dim

       for (rep in 0:c(nb_rep-1)){
         age_Mature = ns1[["ageMature"]]
         age_Mature = lapply(seq(dim(age_Mature)[3]), function(x) as.data.frame(age_Mature[ , , x]))

         for (nm in names(reshape_data_Ne[[sp]][["data_Ne"]])){
            cat("\nProcessing trait ", nm)
            #  fct_i = ifelse(nm %in% c("mortalityByAge", "abundanceByAge"), sum, mean)
            fct_i = ifelse(nm %in% c("eggNumberByAge"), "sum", 
                    ifelse(nm %in% c("mortalityByAge", "meanSizeByAge", "fecundityByAge", "Mpred", "Mstar", "Mnat", "F", "Z", "Mfor", "Dis"), "mean", 
                    ifelse(nm %in% c("sizeMature", "abundanceMature"), "mean", 
                    ifelse(nm %in% c("schoolNumber", "biomassByAge", "abundanceByAge"), "age_1", "total"))))
            new_nm <- paste0(fct_i, "_of_", nm, "_per_y")
            new_trait_dim[which(trait_dim == nm)] <- new_nm
            data_t_t <- reshape_data_Ne[[species[1]]][["data_Ne"]][[nm]]   
            data_t <- data.frame("cohort" = stringr::str_replace(colnames(data_t_t[,,rep+1]), pattern="cohort_", replacement=""))
            if (fct_i %in% c("sum")) data_t[,"value"] <- apply(data_t_t[,,rep+1], 2, function(x) sum(x, na.rm=T))
            if (fct_i %in% c("mean")) {
              k=80
              for (y in colnames(data_t_t[,,rep+1])){
                k=k+1
                tempo <-  weighted.mean(x=data_t_t[,y,rep+1],
                                        w=reshape_data_Ne[[sp_j]][["data_Ne"]][["abundanceByAge"]][,y,rep+1], 
                                        na.rm = T)
                rownames(tempo) <- NULL
                data_t[which(data_t$cohort == k),"value"] <- tempo
              }
            }
            if (fct_i %in% c("age_1"))   data_t[,"value"] <- apply(data_t_t[,,rep+1], 2, FUN=function(x) x[2])
            if (fct_i %in% c("age_mat")) {
              age_mat_sp_rep = age_Mature[[rep+1]][,species[1]]
              vec_y = c()
              for (y in 1:length(age_mat_sp_rep)){
                vec_tempo <- data_t_t[which(rownames(data_t_t) == ceiling(age_mat_sp_rep[y])),y, rep+1]
                vec_y[y] <- ifelse(length(vec_tempo)==0, NA, vec_tempo)
              }
              data_t[,"value"] <- vec_y
            }    
            data_t[,"Species"] <- species[1]

            for (sp_j in species[-1]){
            data_t_t <- reshape_data_Ne[[sp_j]][["data_Ne"]][[nm]]
            data_tj <- data_t
            data_tj[,"cohort"] <- stringr::str_replace(colnames(data_t_t[,,rep+1]), pattern="cohort_", replacement="")
            if (fct_i %in% c("sum")) data_tj[,"value"] <- apply(data_t_t[,,rep+1], 2, function(x) sum(x, na.rm=T))
            if (fct_i %in% c("mean")) {
              k=80
              data_tj = data.frame("cohort"=81:200)
              for (y in colnames(data_t_t[,,rep+1])){
                k=k+1
                tempo <-  weighted.mean(x=data_t_t[,y,rep+1],
                                        w=reshape_data_Ne[[sp_j]][["data_Ne"]][["abundanceByAge"]][,y,rep+1], na.rm = T)
                rownames(tempo) <- NULL
                data_tj[which(data_tj$cohort == k),"value"] <- tempo
              }
            }
            if (fct_i %in% c("age_1"))       data_tj[,"value"] <- apply(data_t_t[,,rep+1], 2, FUN=function(x) x[2])
            if (fct_i %in% c("age_mat")) {
                age_mat_sp_rep = age_Mature[[rep+1]][,sp_j]
                for (y in 1:length(age_mat_sp_rep)){
                  vec_tempo <- data_t_t[which(rownames(data_t_t) == ceiling(age_mat_sp_rep[y])),y, rep+1]
                  vec_y[y] <- ifelse(length(vec_tempo)==0, NA, vec_tempo)
                }
                data_tj[,"value"] <- vec_y
            }    
            data_tj[,"Species"] <- sp_j
            data_t <- rbind(data_t, data_tj)
            }#species
            data_t$TRAIT <- new_nm
            
            if (nm == names(reshape_data_Ne[[species[1]]][["data_Ne"]])[1]){
            data_rep <- data_t
            } else {
            data_rep <- rbind(data_rep, data_t)
            }
            #cat(colnames(data_t))
         }#trait 
         data_rep$replicate <- rep
         data_rep$var_F_for_sp <- sp
         if (rep == 0){
           dat <- data_rep
         } else {
           dat <- rbind(dat, data_rep)
         }
       }#replicate
       dat$Fish_press <- paste0("Fmax=", Fmax, " Fstart=", Fstart)
       # Build dataset across F scenarios
       if (Fmax == Fmax_1 & Fstart == Fstart_1) {
         data_complete <- dat
       } else { data_complete <- rbind(data_complete, dat) }
   
     }#Fstart
   }#Fmax
  data_complete$cohort <- as.numeric(data_complete$cohort)
  
  # Save dataset for later reuse
  pathoutputRDATA <-paste0(pathout, "/RDATA/Fmax", Fmax, "_Fstart", Fstart, "_sp", sp)
  if (!dir.exists(pathoutputRDATA)){dir.create(pathoutputRDATA, recursive = T)}
  # if (!file.exists(paste0(pathoutputRDATA, "/data_traits_summary.RData"))){
    save(data_complete, file=paste0(pathoutputRDATA, "/data_traits_summary.RData"))
  # } else {
    # load(file=paste0(pathoutputRDATA, "/data_traits_summary.RData"))
    #   data_save <- rbind(loaded_data, data_complete)
    #   save(data_save, file=paste0(pathoutputRDATA, "/data_traits_summary.RData"))
  # }
  

  # Summarize per cohort and plot mean +/- sd across replicates
  #data_complete$value[which(data_complete$TRAIT == "fecundityByAge")] <- log(data_complete$value[which(data_complete$TRAIT == "fecundityByAge")]+10^{-8}) # easier to plot
  data_tot_summary <- data_complete %>%
    group_by(cohort, Species, TRAIT, Fish_press) %>%
    summarise(mean = mean(value, na.rm = TRUE),
              sd   = sd(value, na.rm = TRUE))   
  data_tot_summary$cohort <- as.numeric(data_tot_summary$cohort)
  # data_tot_summary_t <- data_tot_summary[which(data_tot_summary$TRAIT %in% 
  #                c("abundanceByAge", "fecundityByAge", "mortalityByAge", "Mpred", "Mstar", "Mnat", "F", "Z", "Mfor", "Dis")),]
  data_tot_summary <- data_tot_summary[which(data_tot_summary$cohort %in% (min(data_tot_summary$cohort)+19):(max(data_tot_summary$cohort)-20)), ]
  # data_tot_summary <- data_tot_summary[-which(data_tot_summary$TRAIT %in%  
  #                c("abundanceByAge", "fecundityByAge", "mortalityByAge", "Mpred", "Mstar", "Mnat", "F", "Z", "Mfor", "Dis")), ]
  # data_tot_summary <- rbind(data_tot_summary_t, data_tot_summary)
  
  
    
  # Prepare plots per trait
  plot_trait = list()
  
  for (TRAIT in new_trait_dim){
    data_tot_summary_tr    <- data_tot_summary[grep(data_tot_summary$TRAIT, pattern=TRAIT),]
    for (sp_plotted in sp){

      data_tot_summary_tr_sp <- data_tot_summary_tr[which(data_tot_summary_tr$Species==sp_plotted),] 
      # reorder Fmax/Fstart factor for prettier plots
      data_tot_summary_tr_sp$Fish_press    <- factor(as.factor(data_tot_summary_tr_sp$Fish_press), levels = ideal_order_factor)
      data_tot_summary_tr_sp$cohort <- as.numeric(data_tot_summary_tr_sp$cohort)
      # If no data, return an empty plot placeholder
      if(nrow(data_tot_summary_tr_sp)==0){
        plot_trait[[TRAIT]][[sp_plotted]] <- ggplot(data = data_tot_summary_tr_sp, aes(x = cohort, group = Fish_press)) 
        next}
      
      # Reference scenario (first encountered) plotted as dashed line
      reference <- subset(na.omit(data_tot_summary_tr_sp), Fish_press == paste0("Fmax=", F_vec[[sp]][1], " Fstart=", F_vec[[sp]][2]))


      x_axis <- "Cohort"
      # tr = ifelse(grepl(TRAIT, pattern= "fecundityByAge"), "log egg number", TRAIT)
      title = paste0("By cohort ", sp_plotted, " ranges of ", TRAIT, "\nF varies for", sp, ", reference is", F_vec[[sp]])
      plot_trait[[TRAIT]][[sp_plotted]] <- ggplot(data = data_tot_summary_tr_sp, aes(x = cohort, group = Fish_press)) +
        geom_ribbon(aes(ymin = mean-sd, ymax = mean + sd),
                    alpha = 0.5, fill = "grey70")+
        # 1) Layer: each panel gets its own group (via facet_wrap)
        geom_line(aes(y = mean)) +
        # 2) Layer: reference group A, drawn in all panels
        geom_ribbon(data = reference,
                    fill = "grey70", alpha = 0.5,
                    aes(x = cohort, ymin = mean-sd, ymax = mean+sd, group = 1)) +
        geom_line(data = reference,
                  aes(x = cohort, y = mean, group = 1),
                  linetype = "dashed", linewidth = 0.7) +
        theme_bw()+
        ylab(label = TRAIT)+
        xlab(label = x_axis)+
        labs(title = title)+
        facet_wrap(~ Fish_press)
    }
    
    # Save plot per trait
    pdf(file = paste0(pathout, "/plot_", TRAIT, "_Sp_fish_", sp, ".pdf"), width = 10, height = 8)
    par(mfrow=c(1,1))
    print(plot_trait[[TRAIT]])
    dev.off()
    
  } #trait plot
} #sp plot
}#function
