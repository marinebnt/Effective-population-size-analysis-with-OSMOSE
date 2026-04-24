
library(osmose)
library(ggplot2)

repo
# initialisation
configDir <- "osmose-ns-v4.4.1_clim"
outputDir_Ne =file.path(configDir, "output", "Ne","Fisheries")
dir.create(outputDir_Ne)
source(paste0(configDir, "/R/functions/functions_for_Ne.R"))
source(file.path(configDir, "R/functions/functions_Ne.R"))
source(file.path(configDir, "R/functions/save_plots_ne.R"))
source(file.path(configDir, "R/functions/functions_temporalNe.R"))
source(file.path(configDir, "R/functions/save_plots_ne_tempo.R"))
source(file.path(configDir, "R/functions/Data_abund_mat_by_age.R"))

# load config
conf <- osmose::read_osmose(paste0(configDir, "/New_param_Eco_Bioen_neutral.R"))

# store usefull parameters of the config 
dt        <- get_par(conf, "simulation.time.ndtPerYear")/get_par(conf, "output.recordfrequency.ndt")
# ny        <- get_par(conf, "simulation.time.nyear")-get_par(conf, "output.start.year")
time_vec  <- seq(0,ny, by=1/dt)[-1] # if year 0 is removed
years     <- time_vec[length(time_vec)]
nb_rep    <- get_par(conf, "simulation.nsimulation")
nb_sp     <- get_par(conf, "simulation.nspecies")
nspecies  <- nb_sp
species   <- unlist(as.matrix(get_par(conf, "species.name.")[1:nb_sp]))
sp_names  <- unlist(as.matrix(get_par(conf, "species.name.")[1:nb_sp]))
species_i <- rownames(as.matrix(get_par(conf, "species.name.")[1:nb_sp]))
sex_ratio <- get_par(conf, "sexratio")
ncohort   <- years
dim_vec   <- c("abundanceByAge", "eggNumber", "schoolNumber")
max_age   <- get_par(conf, "age.max")

# list the folders where the outputs = f(Fmax) are stored
dir_list <- list.dirs(path = "osmose-ns-v4.4.1_clim/output", full.names = T)
dir_list <- dir_list[grep("Fmax", dir_list)]


for (sp in c("Sandeel", "AtlanticHerring")){#c("Sandeel")){#,  #"Saithe",
             #"CommonSole", "NorwayPout",
             #"Sandeel")){ #, "Saithe", "Hake" "AtlanticMackerel" "AtlanticHerring" "CommonSole" "EuropeanPlaice" "NorwayPout"
# for (sp in c("Sandeel")){
  q <- list()
  p <- list()
  mort_htz_l <-  list()
  mort_Ne_l <- list()
  multip <- c()
  cat("\n", sp)
  dir_sp <- dir_list[grep(sp, dir_list)]    
  dir_sp <- dir_sp[grep("Fisheries", dir_sp)]
  seq_multipl <- c("Fmax2_0.5/","Fmax2_1/","Fmax2_1.5/")
  df <- data.frame("Time"=0:c(years-1))
  k=0
  i=0
  for (m in c(0.5,1,1.5)){
    if(k==0) k=k+2
    else k=k+7
    i=i+1
    long_multi <- seq_multipl[i]
    # load outputs
    dir <- dir_sp[grep(long_multi, dir_sp)][1]
    if (m == 1){dir <- c("osmose-ns-v4.4.1_clim/output/Simu_MULTI/output_m1/")}
    ns1 <- osmose::read_osmose(dir)
    out_plots <- paste0(outputDir_Ne, "/fisheries_", sp) 
    
    # effective population sizes
    cat("\nprepare osmose data to fit the by-cohort requisites\n")
    cat("producing raw by cohort data for Ne and temporal Ne estimation\n")
    data_Ne_raw <- estimate_for_Ne(outputDir_Simu=dir, outputDir_Ne = out_plots, dim_vec = dim_vec, species=sp, nb_rep = nb_rep, OUTPUT_data=ns1)
    cat("producing Ne data")
    data_Ne <- estimate_Ne(outputDir_Ne=out_plots, data_Ne_raw, species=sp, nb_rep, input_dir=dir)
    #cat("producing temporal Ne data\n")
    #tempoNe <- estimate_temporalNe(outputDir_Simu=out_plots, data_Ne, species, conf, seeding_y)
   
      cat("saving htz data\n")
      data <- ncdf4::nc_open(paste0(dir, "/Genetic/ns_observedHtz-", sp, "Simu0.nc"))
      data_nc <- ncdf4::ncvar_get(data, "observed_heterozygosity")
      trait1_data <- data_nc[, 1, ]
      ncdf4::nc_close(data)
      mean_htz <- colMeans(trait1_data)
      time <- 0:c(dim(data_nc)[3]-1)
      data_htz <- data.frame(Time=time, HTZ=mean_htz)
      save(data_htz, file=paste0(outputDir_Ne, "/htz", sp, ".RData"))
      data_htz <- get(load(paste0(outputDir_Ne, "/htz", sp, ".RData")))
    
    df[,k] <- data_htz$HTZ
    df[,k+1] <- as.numeric(as.matrix(data_Ne$Ne[[sp]]$rep_1[1,]))
    df[,k+2] <- rowSums(ns1$mortalityByAge[[sp]]$Mpred)
    df[,k+3] <- as.numeric(as.matrix(data_Ne$Vk_point[[sp]]$rep_1[1,]))
    df[,k+4] <- as.numeric(as.matrix(colSums(data_Ne$kx_bar[[sp]]$rep_1)))
    df[,k+5] <- as.numeric(as.matrix(data_Ne$longevity[[sp]]$rep_1[1,]))
    df[,k+6] <- as.numeric(as.matrix(data_Ne$abundanceByAge[[sp]][1,]))
    colnames(df)[k] <- paste0("htzFmax=",m)
    colnames(df)[k+1] <- paste0("NeFmax=", m, " ")
    colnames(df)[k+2] <- paste0("Predation_mortalityFmax=", m, "  ")
    colnames(df)[k+3] <- paste0("Vk_pointFmax=", m, "  ")
    colnames(df)[k+4] <- paste0("kx_barFmax=", m, "  ")
    colnames(df)[k+5] <- paste0("longevityFmax=", m, "  ")
    colnames(df)[k+6] <- paste0("abundanceByAgeFmax=", m, "  ")
  }
  
  rownames(df) <- rownames(as.matrix(data_Ne$Ne[[sp]]$rep_1[1,]))
  df_long_ne <- df[, c(1,grep(pattern = "Ne", colnames(df)))] %>%
    pivot_longer(-Time, names_to = "Fmax_multiplier", values_to = "Effective_population_size_Ne")
  if(any(is.na(df_long_ne$Effective_population_size_Ne))){df_long_ne <- df_long_ne[-which(is.na(df_long_ne$Effective_population_size_Ne)),]}
  df_long_htz<-df[, c(1,grep(pattern = "htz", colnames(df)))] %>%
    pivot_longer(-Time, names_to = "Fmax_multiplier", values_to = "Heterozygosity")
  if(any(is.na(df_long_htz$Heterozygosity))){df_long_htz <- df_long_htz[-which(is.na(df_long_htz$Heterozygosity)),]}    
  
  df_long_ne$Fmax_multiplier <- (stringr::str_remove(df_long_ne$Fmax_multiplier, pattern = "Ne"))    
  df_long_htz$Fmax_multiplier <- (stringr::str_remove(df_long_htz$Fmax_multiplier, pattern = "htz"))
  
  for (column in c("Vk_point", "Predation_mortality", "kx_bar", "longevity", "abundanceByAge")){ # loop to plot all variables

    cat("plot ", column, "\n")
    
    df_long_trait<- df[, c(1,grep(pattern = column, colnames(df)))] %>%
      pivot_longer(-Time, names_to = "Fmax_multiplier", values_to = column)
    
    df_long_trait$Fmax_multiplier <- (stringr::str_remove(df_long_trait$Fmax_multiplier, pattern = column))
    # df_long_ne$Fmax_multiplier <- stringr::str_remove(df_long_ne$Fmax_multiplier, pattern = "\\.2")
    # if(any(is.na(df_long_trait$Mortality_tot))){df_long_trait <- df_long_trait[-which(is.na(df_long_htz$Mortality_tot)),]}
    # df_long <- full_join(df_long_ne, df_long_htz)
    # assign(x = paste0("df_long", sp), value = df_long_ne)
    
    mort_htz <-  max(df_long_trait[,3], na.rm = T)/max(df_long_htz$Heterozygosity, na.rm = T)
    mort_Ne <-  as.numeric(median(as.numeric(as.matrix(df_long_trait[,column])), na.rm = T)/max(df_long_ne$Effective_population_size_Ne, na.rm = T))
    mort_htz_l <- c(mort_htz_l, mort_htz)
    mort_Ne_l <- c(mort_Ne_l, mort_Ne)
    df_long_trait[,3] <- df_long_trait[,3]/mort_Ne
    m_vec = max(df_long_htz[,"Heterozygosity"], na.rm = T)/max(df_long_ne[-which(is.nan(df_long_ne$Effective_population_size_Ne)),"Effective_population_size_Ne"], na.rm = T)
    multip <- c(multip, m_vec)
    p[[sp]] <- ggplot()+
      geom_line(data = df_long_ne[which(df_long_ne$Time %in% c(70:100)),], aes(x=Time, 
                  y=Effective_population_size_Ne, color=Fmax_multiplier),
                size=2, alpha=0.5)+
      geom_line(data=df_long_trait[which(df_long_ne$Time %in% c(70:100)),],
                aes(x=Time, y=get(column), linetype=Fmax_multiplier),
                size=1, alpha=0.5)+
      theme(legend.position = "none")+  
      theme_bw()+
      labs(title=paste0(sp," Ne\ according to Fmax"))+
      scale_fill_manual(values = c("Fmax=0.5" = "darkgreen",
                                   "Fmax=1" = "firebrick",
                                   "Fmax=1.5" = "royalblue"))+
      scale_linetype_manual(values = c("Fmax=0.5  " = "dotted",
                                       "Fmax=1  " = "solid",
                                       "Fmax=1.5  " = "twodash"))+
      scale_y_continuous(
        name = "Ne",
        sec.axis = sec_axis(~.*mort_Ne,
                            name = column))
    
    q[[sp]] <- ggplot(df_long_htz[which(df_long_htz$Time %in% c(70:100)),], aes(x=Time))+
      geom_line(aes(y=Heterozygosity, color=Fmax_multiplier), 
                size=1, alpha=0.5)+
      theme_bw()+
      labs(title=paste0(sp," heterozygosity according to Fmax"))+
      scale_fill_manual(values = c("Fmax=0.5" = "darkgreen",
                                   "Fmax=1" = "firebrick",
                                   "Fmax=1.5" = "royalblue"))
    
    pdf(file = paste0("osmose-ns-v4.4.1_clim/output/Ne/Fisheries/2Ne_Htz", sp, column, ".pdf"), width=15, height=6)
    # pdf(file = paste0(paste0(pathNCDF, "_FINAL_plot/", "tot", ".pdf")), width=25, height=15)
    par(mfrow=c(1,1))
    for (k in 1:length(p)) {
      multip <- m_vec[k]
      mort_Ne <- mort_Ne[[k]]
      print(cowplot::plot_grid(q[[k]],p[[k]]))
    }
    dev.off()
    # print(cowplot::plot_grid(q[[sp]],p[[sp]]))
  }
}





