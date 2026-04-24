

PLOT_trait_replicates2 <- function(SP){

# setwd("X:/osmose-ns")

repo        = "Simu_LONG"
multip_name = "y"
multipliers = c(1,0.5,0.75,1.25,1.5)

trait_dim= c("abundanceByAge", "biomass", "school", "yield", "fecundity", "mortality")


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

pathout <- paste0(configDir, "/output/PLOTS/Trait_variables")

plot_Neth_Nalpha    <- list()
F_vec     <- list()

for (sp in SP){
  sp_i <- which(species == sp)-1
  # open one dir at a time
  sp_dir <- list_d[grep(pattern = paste0("sp", sp_i), list_d)]
  
  plot_Neth_N_alpha_NT_sp     <- list()
  Fmax_1 = -99
  Fstart_1 = -99
  
  for (Fmax in c(1,0.5)){#multipliers){
    Fmax_dir <- sp_dir[grep(pattern = paste0("Fmax", Fmax, "_"), sp_dir)]
    for (Fstart in 0.5){#multipliers){
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
      outputDir_Ne    = file.path(configDir, "output", "PLOTS", paste0("Fmax", Fmax, "Fstart", Fstart, "sp", sp_i))
      if(!dir.exists(outputDir_Ne)) {dir.create(outputDir_Ne, recursive = T)}
      ns1 = read_osmose(path = paste0(outputDir_Simu, "/"))
      
      for (TRAIT in trait_dim[1]){
        data_t = NULL
        data <- ns1[[TRAIT]]
        for (i in 1:c(nb_rep-1)){            
          # cat(TRAIT, "\n")
          if(TRAIT == trait_dim[1] && i==1){
            if (TRAIT == "mortality"){
              data_tot <- data.frame("time"= rownames(data$AtlanticCod$eggs))
              rownames(data_tot) = data_tot$time
              data_tot[, species[1]] <- rowSums(ns1$mortalityByAge[[species[1]]]$Mpred[,,i])
              for (spi in species[-1]){
                data_tot[, spi] <- rowSums(ns1$mortalityByAge[[spi]]$Mpred[,,i])
              }
            }
            if (TRAIT %in% c("school", 'fecundity')) {
              if (TRAIT == c("fecundity")) {
                data_tot_all <- read.csv(paste0(outputDir_Simu, "/AgeIndicators/ns_fecundityDistribByAge_Simu", c(0), ".csv"), skip=c(1))
                data_tot_sum <- t(sapply(names(table(data_tot_all$Age)), FUN = function(x) {colSums(data_tot_all[which(data_tot_all$Age == x), ])}))
              }
              if (TRAIT == c("school")) {
                data_tot_all <- read.csv(paste0(outputDir_Simu, "/ns_nschoolDistribByAge_Simu", c(0), ".csv"), skip=c(1))
                data_tot_sum <- NULL
              }
              data_tot_sum <- t(sapply(names(table(data_tot_all$Time)), FUN = function(x) {colSums(data_tot_all[which(data_tot_all$Time == x), ])}))
              data_tot <- data_tot_sum[, -which(colnames(data_tot_sum) %in% c("Time", "Age"))]
              data_tot <- as.data.frame(data_tot[, which(colnames(data_tot) %in% species)])
            }
            if (TRAIT %in% c("abundanceByAge")){
              data_tot <- data[[1]][,,i]
              species_1<- species[1]
              data_tot <- data.frame(time = colnames(data_tot))
              data_tot[, species_1] <- colSums(data[[species_1]][,,i])
              for (spi in species[-1]){
                data_tot[[spi]] <- colSums(data[[spi]][,,i])
              }
            }
            if (!(TRAIT %in% c("mortality", "school", 'fecundity', 'abundanceByAge'))){
              data_tot <- as.data.frame(data[,,i]) 
            }
            if(!(TRAIT %in% c("abundanceByAge"))) {
              data_tot$time <- rownames(data_tot)
              data_tot <- data_tot%>%
                pivot_longer(!time, names_to = "Species", values_to = "value")
              data_tot$rep <- 1
              data_tot$Fish_press = paste0("Fmax=", Fmax, " Fstart=", Fstart)
              data_tot$TRAIT = TRAIT
            }
            if((TRAIT %in% c("abundanceByAge"))) {
              data_tot <- data_tot%>%
                pivot_longer(!time, names_to = "Species", values_to = "value")
              data_tot$rep <- 1
              data_tot$Fish_press = paste0("Fmax=", Fmax, " Fstart=", Fstart)
              data_tot$TRAIT = TRAIT
            }
            
          } else{
            if (TRAIT %in% c("abundanceByAge")){
              data_t <- data[[1]][,,i]
              species_1<- species[1]
              data_t <- data.frame(time = colnames(data_t))
              data_t[, species_1] <- colSums(data[[species_1]][,,i])
              for (spi in species[-1]){
                data_t[[spi]] <- colSums(data[[spi]][,,i])
              }
            }
            if (!(TRAIT %in% c("mortality", "school", 'fecundity', "abundanceByAge"))){
              data_t <- as.data.frame(data[,,i]) 
            }
            if (TRAIT == "mortality"){
              data_t <- data.frame("time"= rownames(data$AtlanticCod$eggs))
              rownames(data_t) = data_t$time
              data_t[, species[1]] <- rowSums(ns1$mortalityByAge[[species[1]]]$Mpred[,,i])
              for (spi in species[-1]){
                data_t[, spi] <- rowSums(ns1$mortalityByAge[[spi]]$Mpred[,,i])
              }
            }
            if (TRAIT %in% c("school", 'fecundity')) {
              if (TRAIT == c("fecundity")) data_t_all <- read.csv(paste0(outputDir_Simu, "/AgeIndicators/ns_fecundityDistribByAge_Simu", c(i-1), ".csv"), skip=c(1))
              if (TRAIT == c("school")) data_t_all <- read.csv(paste0(outputDir_Simu, "/ns_nschoolDistribByAge_Simu", c(i-1), ".csv"), skip=c(1))
              data_t_sum <- t(sapply(names(table(data_t_all$Time)), FUN = function(x) {colSums(data_t_all[which(data_t_all$Time == x), ])}))
              data_t <- data_t_sum[, -which(colnames(data_t_sum) %in% c("Time", "Age"))]
              data_t <- as.data.frame(data_t[, which(colnames(data_t) %in% species)])
            }
            if(!(TRAIT %in% c("abundanceByAge"))) {
              data_t$time <- rownames(data_t)
              data_t <- data_t%>%
                pivot_longer(!time, names_to = "Species", values_to = "value")
              data_t$rep <- i
              data_t$Fish_press = paste0("Fmax=", Fmax, " Fstart=", Fstart)
              data_t$TRAIT = TRAIT         
            }
            if((TRAIT %in% c("abundanceByAge"))) {
              data_t <- data_t%>%
                pivot_longer(!time, names_to = "Species", values_to = "value")
              data_t$rep <- 1
              data_t$Fish_press = paste0("Fmax=", Fmax, " Fstart=", Fstart)
              data_t$TRAIT = TRAIT
            }
            data_tot <- rbind(data_tot, data_t)  
          }#rep_cdtion
          # cat(TRAIT, "\n")
          # print(head(data_t))
        }#rep
      }#TRAIT
      if (Fmax == Fmax_1 & Fstart == Fstart_1) {
        data_complete <- data_tot
      }else {
        data_complete <- rbind(data_complete, data_tot)
      }
      # p(paste0("\n\t--->done: folder Fmax", Fmax, " Fstart", Fstart, " sp", sp, "\n"))
    }#Fstart
  }#Fmax
  
  # PLOT the mean and the std deviation 
  
  data_tot_summary <- data_complete %>%
    group_by(time, Species, Fish_press, TRAIT) %>%
    summarise(mean = mean(value, na.rm = TRUE),
              sd   = sd(value, na.rm = TRUE))   
  data_tot_summary$time <- as.numeric(data_tot_summary$time)
  if(!(TRAIT %in% c("abundanceByAge"))) { 
    data_tot_summary <- data_tot_summary[which(data_tot_summary$time %in% (min(data_tot_summary$time)+19):(max(data_tot_summary$time)-20)), ]
  }
  


  
  pathoutputRDATA <-paste0(pathout, "/RDATA/Fmax", Fmax, "_Fstart", Fstart, "_sp", sp)
  if (!dir.exists(pathoutputRDATA)){dir.create(pathoutputRDATA, recursive = T)}
  save(data_tot_summary, file=paste0(pathoutputRDATA, "/data_traits_summary.RData"))
    
  plot_trait = list()
  
  for (TRAIT in trait_dim){
    for (sp_plotted in species){  
      
      x_axis <- ifelse(TRAIT=="abundanceByAge", "Age", "Time")
      data_tot_summary_tr    <- data_tot_summary[which(data_tot_summary$TRAIT==TRAIT),]
      data_tot_summary_tr_sp <- data_tot_summary_tr[which(data_tot_summary_tr$Species==sp_plotted),]
      
      reference <- subset(na.omit(data_tot_summary_tr_sp), Fish_press == paste0("Fmax=", F_vec[[sp]][1], " Fstart=", 
                                                                                           F_vec[[sp]][2]))[, c("mean", "sd", "time")]
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
        ylab(label = TRAIT)+
        xlab(label = x_axis)+
        labs(title = paste0(sp_plotted, " ranges of ", TRAIT, "\nF varies for", sp, ", reference is", F_vec[[sp]]))+
        facet_wrap(~ Fish_press)
    }
    
    pdf(file = paste0(pathout, "/plot_", TRAIT, "_Sp_fish_", sp, ".pdf"), width = 10, height = 8)
    par(mfrow=c(1,1))
    print(plot_trait[[TRAIT]])
    dev.off()
    
  }
} #sp
}







