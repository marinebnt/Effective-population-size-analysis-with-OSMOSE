
main="New_gpu_conf.R"
source(file.path("osmose-ns-v4.4.1_clim/R/functions/source.R"))
repo        = "RESTART"
multip_name <- c("y")
multipliers = c("0.5", "1", "1.5")
multipliers = get_par(conf, "simulation.time.nyear") 

# function to prepare environment : launch simulation if repo missing
function_INIT1(MULTIP=multipliers, REPO=repo)

# prepare environment for this particular plot
rep=1
incr_multiplier <- 0


# multipliers = "output_etaequalsone_bioen_neutral"

sp_focus <- "Sandeel"
for (m in multipliers){
  
  cat("load simulation output")
  simulation = paste0(repo, "/output_", multip_name, m)
  Ne_version = paste0("plot_biomass", multip_name, m)
  cat("\n\t---> folder ", simulation, "\n")
  ns1 = NULL
  outputDir  = file.path(configDir, "output", simulation)
  outputDir_Ne =paste0(configDir, "/output/", "Ne","/PLOT_Biomass/")
  if (any(list.files(outputDir) %in% species)){
    outputDir    <- paste0(outputDir, "/", sp_focus)
    outputDir_Ne <- paste0(outputDir_Ne, "Fmax_on_", sp_focus, "/")
  }
  if(!dir.exists(outputDir_Ne)) {dir.create(outputDir_Ne)}
  ns1 = read_osmose(path = outputDir)
  
  p <- list()
  for (i in 1:ns1$model$simus){
    dataset_plot <- as.data.frame(ns1$biomass[,,i])
    dataset_plot$time <- as.numeric(rownames(dataset_plot))
    dataset_plot <- dataset_plot %>% 
      pivot_longer(!time, values_to = "biomass", names_to = "species")
    dataset_plot$biomass <- as.numeric(dataset_plot$biomass)/10^3
    p[[i]] <- ggplot(dataset_plot, aes(x=time, y=biomass))+ 
        geom_line()+
        facet_wrap(~ species, ncol=4, nrow=4)+
        theme_minimal()+
        labs(title=paste0("Replicate ", i, "x 10^3"))
  }
  
  pdf(file = paste0(outputDir_Ne, Ne_version,"ALLrep.pdf"))
  for (i in 1:ns1$model$simus){
    print(p[[i]])
  }
  dev.off()
}
