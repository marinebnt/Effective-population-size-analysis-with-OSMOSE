library(ggplot2)
library(osmose)
source("R/functions/source.R")

path <- "new_config_Fmax_Fstart"
path_init <- 'output/DATASETS/Chap2/FINAL_Ne_Htz_mature/INIT'
nfiles <- list.files(path, pattern = ".csv")
inifiles <- list.files(path_init, pattern = ".RData", recursive = TRUE, full.names = TRUE)
inifiles <- inifiles[grep(pattern = "INIT2.RData", inifiles)]
multip <- c(1, 1.25, 1.5, 0.5, 0.75)
ideal_order_factor <- paste0("Fmax=", rep(c(0.5,0.75,1,1.25,1.5), 5),
                               "_Fstart=", c(rep(0.5, 5),rep(0.75, 5),rep(1, 5),rep(1.25, 5),rep(1.5, 5)))
  

sp_vec = c("AtlanticHerring",  "AtlanticMackerel", "Sandeel",  "EuropeanSprat", "NorwayPout",       
       "EuropeanPlaice", "CommonSole", "Saithe", "AtlanticCod", "Haddock",  "HorseMackerel",
       "Whiting",  "CommonDab", "GreyGurnard", "Hake")

plot_sp <- list()
gene_l_sp_coho <- c()
gene_l_sp_year <- c()
data_surv_curve <- list()
for (species in sp_vec){

  colnames_df <- c("length", "F", "Fmax_Fstart")
  dataset_plot <- data.frame(matrix(nrow=0, ncol = length(colnames_df)))
  colnames(dataset_plot) = colnames_df
  id_sp <- grep(pattern = species, sp_vec)-1
  id_sp_files <- grep(pattern = paste0("sp", id_sp, ".csv"), nfiles)
  max_l <- as.numeric(get_par(conf, paste0("linf.sp", id_sp))[1])
  for (fmax in 1){
    for (fstart in 1){ #multip

      # find the right F pressure to plot 
      id_mx_files <- grep(pattern = paste0("Fmax", fmax, "_"), nfiles)
      id_st_files <- grep(pattern = paste0("Fstart", fstart, "_"), nfiles)
      files_id <- table(c(id_sp_files, id_mx_files, id_st_files))
      file_chosen <- nfiles[as.numeric(names(files_id)[which(files_id == 3)])]
      tab <- read.csv(file.path(path, file_chosen))
      colnames(tab) <- as.numeric(stringr::str_remove(colnames(tab), pattern = "X"))
      by=as.numeric(colnames(tab)[2])-as.numeric(colnames(tab)[1])
      add_cols <- round(max_l-max(as.numeric(colnames(tab)),na.rm=T))/by
      colname_max <- max(as.numeric(colnames(tab)), na.rm=T)
      if (add_cols > 0) {
        tab <- cbind(tab, matrix(rep(tab[,ncol(tab)], add_cols), nrow=nrow(tab)))
        colnames(tab) <- seq(as.numeric(colnames(tab)[1]), colname_max+add_cols*by, by=by)
      }
      
      # save F pressure in the dataset to plot
      add_mat <-  matrix(rep(NA, ncol(tab)*ncol(dataset_plot)), nrow=ncol(tab), ncol=ncol(dataset_plot))
      colnames(add_mat) <- colnames(dataset_plot)
      dataset_plot2 <- rbind(dataset_plot, add_mat)
      dataset_plot2$F <- c(dataset_plot$F[], as.numeric(tab[1,]))
      dataset_plot2$Fmax_Fstart <- c(dataset_plot$Fmax_Fstart, rep(paste0("Fmax=", fmax, "_Fstart=", fstart), ncol(tab)))
      dataset_plot2$length <- c(dataset_plot$length, as.numeric(colnames(tab)))
      dataset_plot <- dataset_plot2

      # find the right generation length to plot
      inifiles_id_sp <- inifiles[grep(pattern = species, inifiles)]
      inifiles_id_fmax <- inifiles[grep(pattern = paste0("Fmax", fmax, "_"), inifiles)]
      inifiles_id_fstart <- inifiles[grep(pattern = paste0("Fstart", fstart, "_"), inifiles)]
      inifiles_id <- table(c(inifiles_id_sp, inifiles_id_fmax, inifiles_id_fstart))
      inifiles_id <- names(inifiles_id)[which(inifiles_id == 3)]
      # inifiles_id <- "osmose-ns-v4.4.1_clim/output/PLOTS/All_genetic_drift_indicators_test/INIT/Fmax0.5_Fstart0.5_spAtlanticHerring/INIT2.RData" #test
      gene_length_data <- get(load(inifiles_id))
      # gene_length_sp <- gene_length_data[["tempoNe"]][[species]]$generation_length
      # mean_gene_length_sp <- mean(apply(gene_length_sp, 2, mean, na.rm = T))
      # gene_l_sp_year[species] <- mean_gene_length_sp
      gene_length_sp <- gene_length_data[["data_Ne"]][[species]]$generation_length
      mean_gene_length_sp <- mean(apply(gene_length_sp, 2, mean, na.rm = T))
      gene_l_sp_coho[species] <- mean_gene_length_sp

    }# fstart
  }# fmax

 # clean dataset to plot 
  # normalize F by max F at Fmax=1 and keep length up to the max length at F=1
  max_length <- max(dataset_plot$length[which(dataset_plot$Fmax_Fstart == paste0("Fmax=", 1, "_Fstart=", 1))], na.rm=T)
  max_F <- max(dataset_plot$F[which(dataset_plot$Fmax_Fstart == paste0("Fmax=", 1, "_Fstart=", 1))])
  #normalize F
  dataset_plot$F_sc <- dataset_plot$F/max_F
  # keep only length up to the max length at F=1 
  dataset_plot <- dataset_plot[which(dataset_plot$length <= max_length), ]
  # order factor levels
  dataset_plot$Fmax_Fstart        <- factor(as.factor(dataset_plot$Fmax_Fstart), levels = ideal_order_factor) 
  dataset_plot$Fmax <- unlist(stringr::str_split(dataset_plot$Fmax_Fstart, pattern = "_Fstart"))[seq(1,(2*length(dataset_plot$Fmax_Fstart)), by=2)]
 #produce plots
  colnames(dataset_plot)[which( colnames(dataset_plot) == "Fmax")] <- "\u03B8F"
  dataset_plot[, "θF"] <- stringr::str_replace(dataset_plot$θF, pattern = "Fmax", replacement = "θF")
  plot_sp[[species]] <- ggplot(dataset_plot, aes(x=length, y=F_sc))+ 
        geom_line(lwd=1.2, col="blue3")+
        theme_test()+
        theme(text = element_text(size = 13),
          plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5))+
        ylab("")+
        # facet_wrap(~ Fmax, ncol=5, nrow=1)+
        labs(title=species)+ #,        "\nMean generation length by year: ", round(gene_l_sp_year[species], 2)
        theme(legend.title = element_blank())
  data_surv_curve[[species]] <- dataset_plot
  
}# species

#save
pathout <- paste0("output/PLOTS/selectivity_curves_per_Fmax")
save(data_surv_curve, file=paste0(pathout, "dataset_surv-curve.RData"))

if(!dir.exists(pathout)) {dir.create(pathout, recursive = T)}
# for (species in sp_vec){
#   png(file = paste0(pathout, "/plot_Fpressure_",species, ".png"), width = 1500, height = 350)
#   par(mfrow=c(1,1))
#   print(plot_sp[[species]])
#   dev.off()
# }
library(grDevices)
par(mfrow=c(3,2))
windowsFonts(myCustomWindowsFontName=windowsFont("DejaVu Sans"))
cairo_pdf(file = paste0(pathout, "/plot_Fpressure_Fmax1",species, ".pdf"),
          width = 14, height =15, onefile=T)
ggarrange(plotlist=plot_sp, ncol=3, nrow=5)
dev.off()

  