# compare the Ne/Nc estimates depending on different Nc shapes



library(ggpubr)


source(file.path("osmose-ns-v4.4.1_clim/R/functions/source.R"))
repo        = "Simu_MULTI"
multip_name = "m"
multipliers = na.omit(as.numeric(stringr::str_remove(list.files("osmose-ns-v4.4.1_clim/output/Simu_MULTI"), "output_m")))


# function to prepare environment : launch simulation if repo missing
function_INIT1(MULTIP=multipliers, REPO=repo)


# prepare environment for this particular plot
rep=1
m=1


cat("load simulation output")
simulation = paste0(repo, "/output_", multip_name, m)
Ne_version = paste0("School_", multip_name, m)
if (!(m %in% c(1,2,3,4,5,6,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9))){simulation = paste0(m)}
cat("\n\t---> folder ", simulation, "\n")
ns1 = NULL
outputDir_Simu  = file.path(configDir, "output", simulation)
outputDir_Ne =file.path(configDir, "output", "Ne","by_multiplier", Ne_version)
ns1 = read_osmose(path = outputDir_Simu)
if(!dir.exists(outputDir_Ne)) {dir.create(outputDir_Ne)}

# function to prepare environment : 
INIT2 <- function_INIT2(outputDir_Simu, outputDir_Ne, dim_vec, species, nb_rep, OUTPUT_data=ns1, conf, seeding_y, REPO=repo, TEMPORAL = F)
data_Ne_raw  <- INIT2[[1]]
data_Ne      <- INIT2[[2]]
tempoNe      <- INIT2[[3]]

Nc_from_maturity         <- Nc_from_maturity(data_Ne, species)
Nc_only_mature           <- Nc_only_mature(data_Ne, species)
Nc_mature_times_survival <- Nc_mature_times_survival(DATA_NE = data_Ne, species)

g <- list()
q <- list()
p <- list()
m1 <- c()
m2 <- c()
m3 <- c()
for (sp in species){
  dt_Nc <- as.numeric(names(Nc_only_mature[[sp]]))+1 #(we want the position in the vector [1:120], not the time step [0:119])
  f_Ne_Nc <- data.frame("Ne" = data_Ne$Ne[[sp]][[1]][1,dt_Nc], 
                        "Nc_from_maturity" = Nc_from_maturity[[sp]][dt_Nc], 
                        "Nc_only_mature" = Nc_only_mature[[sp]], 
                        "Nc_mature_times_survival" = Nc_mature_times_survival[[sp]],
                        "Ne/Nc_from_maturity" = data_Ne$Ne[[sp]][[1]][1,dt_Nc]/Nc_from_maturity[[sp]][dt_Nc], 
                        "Ne/Nc_only_mature"   = data_Ne$Ne[[sp]][[1]][1,dt_Nc]/Nc_only_mature[[sp]],
                        "Ne/Nc_mature_times_survival" = data_Ne$Ne[[sp]][[1]][1,dt_Nc]/Nc_mature_times_survival[[sp]],
                        "Time" = dt_Nc)
  f_Ne_Nc <- na.omit(f_Ne_Nc)
  f_Ne_Nc <- f_Ne_Nc[which(f_Ne_Nc$Time <101),]
  f_Ne_Nc <- pivot_longer(f_Ne_Nc, cols = c("Ne.Nc_from_maturity" , "Ne.Nc_only_mature", "Ne.Nc_mature_times_survival"))
  if(any(is.infinite(f_Ne_Nc$value))) {f_Ne_Nc <- f_Ne_Nc[-which(f_Ne_Nc$value == Inf),]}
  
  multip <- mean(c(f_Ne_Nc$value), na.rm = T)/(max(f_Ne_Nc$Ne, na.rm = T))
  m1 <- c(m1, multip)
  
  p[[sp]] <- ggplot(f_Ne_Nc, aes(x=Time))+
    # ylim(0, max(df_var[, "max_Ne"], na.rm = T))+
    geom_line(aes(y=Ne), col="red", size=3)+
    geom_line(aes(y=value/multip, col=name), size=3, alpha = 0.5)+
    scale_y_continuous(
      # Features of the first axis
      name = "Ne (red)", #"mean temporal Ne (blue)",
      # Add a second axis and specify its features
      sec.axis = sec_axis(~.*multip, name=paste0("Ne/Nc (shades of blue)"))
    )+
    theme(legend.position="top")+
    # theme(axis.line = element_line(),
    #       plot.margin = margin(10, 10, 10, 30))+
    labs(title = paste0(sp))
  
  multip2 <- median(c(f_Ne_Nc$value[-which(f_Ne_Nc$name=="Ne.Nc_mature_times_survival")]), na.rm = T)/(max(f_Ne_Nc$Ne[-which(f_Ne_Nc$name=="Ne.Nc_mature_times_survival")], na.rm = T))
  m2 <- c(m2, multip2)
  
  ID_remove <- f_Ne_Nc$name[which(f_Ne_Nc$value == max(f_Ne_Nc$value))]
  q[[sp]] <- ggplot(f_Ne_Nc[-which(f_Ne_Nc$name==ID_remove),], aes(x=Time))+
    # ylim(0, max(df_var[, "max_Ne"], na.rm = T))+
    geom_line(aes(y=Ne), col="red", size=3)+
    geom_line(aes(y=value/multip2, col=name), size=3, alpha = 0.5)+
    scale_y_continuous(
      # Features of the first axis
      name = "Ne (red)", #"mean temporal Ne (blue)",
      # Add a second axis and specify its features
      sec.axis = sec_axis(~.*multip2, name=paste0("Ne/Nc (see legend)"))
    )+
    theme(legend.position="top")+
    # theme(axis.line = element_line(),
    #       plot.margin = margin(10, 10, 10, 30))+
    labs(title = paste0(sp))
  
  p[[sp]] <- ggplot(f_Ne_Nc, aes(x=Time))+
    # ylim(0, max(df_var[, "max_Ne"], na.rm = T))+
    geom_line(aes(y=Ne), col="red", size=3)+
    geom_line(aes(y=value/multip, col=name), size=3, alpha = 0.5)+
    scale_y_continuous(
      # Features of the first axis
      name = "Ne (red)", #"mean temporal Ne (blue)",
      # Add a second axis and specify its features
      sec.axis = sec_axis(~.*multip, name=paste0("Ne/Nc (shades of blue)"))
    )+
    theme(legend.position="top")+
    # theme(axis.line = element_line(),
    #       plot.margin = margin(10, 10, 10, 30))+
    labs(title = paste0(sp))
  

  
  ID_keep <- f_Ne_Nc$name[which(f_Ne_Nc$value == min(f_Ne_Nc$value))]
  multip3 <- median(c(f_Ne_Nc$value[which(f_Ne_Nc$name==ID_keep)]), na.rm = T)/(max(f_Ne_Nc$Ne[which(f_Ne_Nc$name==ID_keep)], na.rm = T))
  m3 <- c(m3, multip3)
  g[[sp]] <- ggplot(f_Ne_Nc[which(f_Ne_Nc$name==ID_keep),], aes(x=Time))+
    geom_line(aes(y=value), col="red", size=3)+
    theme(legend.position="top")+
    # theme(axis.line = element_line(),
    #       plot.margin = margin(10, 10, 10, 30))+
    labs(title = paste0(sp))
}



# STORE THE PLOTS
dir.create(paste0(configDir, "/output/Ne"))
pdf(file = paste0(configDir, "/output/Ne/plot_NeonNc_by_Nc.pdf"), width = 15, height = 7)
# pdf(file = paste0(paste0(pathNCDF, "_FINAL_plot/", "tot", ".pdf")), width=25, height=15)
par(mfrow=c(1,1))
for (k in 1:length(p)){
  # multiplier2 = coef_vec[k]
  multip = m1[k]
  multip2 = m2[k]
  print(ggarrange(p[[k]], q[[k]]))
}
dev.off()


# STORE THE PLOTS
dir.create(paste0(configDir, "/output/Ne"))
pdf(file = paste0(configDir, "/output/Ne/plot_MIN_NeonNc_by_Nc.pdf"), width = 5, height = 4)
par(mfrow=c(1,1))
for (k in 1:length(p)){
  # multiplier2 = coef_vec[k]
  multip = m1[k]
  multip2 = m2[k]
  print(g[[k]])
}
dev.off()
