# Test multiple school size multipliers and then explore the effect on Ne & else
# 
# marine beneat
# 02/06/25
# 


source(file.path("osmose-ns-v4.4.1_clim/R/functions/source.R"))
repo        = "Simu_MULTI"
multip_name = "m"
multipliers = na.omit(as.numeric(stringr::str_remove(list.files("osmose-ns-v4.4.1_clim/output/Simu_MULTI"), "output_m")))


# function to prepare environment : launch simulation if repo missing
function_INIT1(MULTIP=multipliers, REPO=repo)


# prepare environment for this particular plot
rep=1
incr_multiplier <- 0
var_to_plot <- c("mean_Ne", "median_Ne", "min_Ne", "max_Ne", "Ne/Nc", 
                 "mean_Ne_tempo", "mean_htz", "quantile95_Ne", "quantile5_Ne", 
                 "quantile95_Ne_tempo", "quantile5_Ne_tempo", "quantile95_htz","quantile5_htz")
var_list <- array(dim = c(length(var_to_plot), length(multipliers) , 16), 
                  dimnames = list(var_to_plot,
                             multipliers,
                             unlist(get_par(conf, "species.name"))[1:16]))


for (m in multipliers){ # loop over the school number multipliers
  
  cat("load simulation output")
  simulation = paste0(repo, "/output_", multip_name, m)
  Ne_version = paste0("School_", multip_name, m)
  if (!(m %in% c(1,2,3,4,5,6,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9))){simulation = paste0(m)}
  cat("\n\t---> folder ", simulation, "\n")
  ns1 = NULL
  outputDir_Simu  = file.path(configDir, "output", simulation)
  outputDir_Ne    = file.path(configDir, "output", "Ne","by_multiplier", Ne_version)
  ns1 = read_osmose(path = outputDir_Simu)
  if(!dir.exists(outputDir_Ne)) {dir.create(outputDir_Ne)}
  
  # function to prepare environment : 
  INIT2 <- function_INIT2(outputDir_Simu, outputDir_Ne, dim_vec, species, nb_rep, OUTPUT_data=ns1, conf, seeding_y, REPO=repo)
  data_Ne_raw  <- INIT2[[1]]
  data_Ne      <- INIT2[[2]]
  tempoNe      <- INIT2[[3]]
  
  incr_multiplier <- incr_multiplier + 1
  
  for(sp in species){
    cat(sp)
    time1 =c(80)
    time2 =c(90)
    
    cohort_time_corr <- data_Ne[["cohort_time_step_corresponding"]]
    co_to_include <- which(colSums(data_Ne[["longevity"]][[sp]][[paste0("rep_", rep)]] )>0)
    
    # what is the first age of each species ? 
    raw_ab_c_i     <- as.data.frame(data_Ne$abundanceByAge[[sp]])
    raw_ab_c_i$age <- as.numeric(rownames(raw_ab_c_i))
    new_age_one    <- as.numeric(data_Ne$new_first_age[sp])
    tab_Nc_Nr      <- data_Ne$table_mature[[sp]]
    
    # VAR/MEAN Ne
    # Nc
    # Nc <- as.data.frame(data_Ne[["abundanceByAge"]][[sp]])
    # Nc   <-c()
    # Nc_h <- c()
    # Nc <- Nc_from_maturity(DATA_NE=data_Ne, species)
    # # Add time as a column
    # Nc$age <- as.numeric(rownames(Nc))
    # # Reshape data to long format
    # df_nc_long_i <- pivot_longer(Nc, cols = -age, names_to = "Cohort", values_to = "Nc")

    
    
    
    # Ne
    Ne <- as.data.frame(data_Ne[["Ne"]][[sp]][[paste0("rep_", nb_rep)]])
    # Add time as a column
    Ne$age <- as.numeric(rownames(Ne))
    # Reshape data to long format
    df_ne_long_i <- pivot_longer(Ne, cols = -age, names_to = "Cohort", values_to = "Ne")
    # match the cohort age with the time step to which it corresponds
    for (co in 0:c(ncohort-1)){
      ID_co <- which(cohort_time_corr==co, arr.ind = T)[,1]
      df_ne_long_i[grep(paste0("\\bcohort_", co, "\\b"), df_ne_long_i$Cohort)[1:length(ID_co)], "Time"] <- ID_co
    }
    # merge Nc and Ne dataset
    data_tot <- full_join(df_nc_long_i, df_ne_long_i)
    
    # Ne tempo data
    df      <- data.frame(tempoNe[sp])
    df$Time <- as.numeric(rownames(df))
    colnames(df) <- c("Ne_tempo", "Time")
    
    # merge Ne, Nc, weighted Ne, tempo Ne
    data_tot <-full_join(data_tot, df)
      
    # VAR/MEAN htz
    if (!file.exists(paste0(outputDir_Ne, "/htz", sp, ".RData"))){
      data <- ncdf4::nc_open(paste0(outputDir, "/Genetic/ns_observedHtz-", sp, "Simu0.nc"))
      data_nc <- ncdf4::ncvar_get(data, "observed_heterozygosity")
      trait1_data <- data_nc[, 1, ]
      ncdf4::nc_close(data)
      mean_htz <- colMeans(trait1_data)
      time <- 0:c(dim(data_nc)[3]-1)
      data_htz <- data.frame(Time=time, HTZ=mean_htz)
      save(data_htz, file=paste0(outputDir_Ne, "/htz", sp, ".RData"))
    }  else {
      data_htz <- get(load(paste0(outputDir_Ne, "/htz", sp, ".RData")))
      cat("=> loading already existing htz info\n")}
    var_htz <- sd(data_htz$HTZ[-which(is.nan(data_htz$HTZ))])
    
    data_htz <- data_htz[which(data_htz$Time <= time2 & data_htz$Time >= time1), ]
    data_tot <- data_tot[which(data_tot$Time <= time2 & data_tot$Time >= time1), ]
    data_tot <- data_tot[which(data_tot$age==min(data_tot$age)),]
    
    var_list["mean_Ne",incr_multiplier,sp] <- mean(data_tot$Ne) #, w=data_tot$Nc, na.rm = T)
    var_list["median_Ne",incr_multiplier,sp] <- median(data_tot$Ne) #, w=data_tot$Nc, na.rm = T)
    var_list["max_Ne",incr_multiplier,sp] <- max(data_tot$Ne) #, w=data_tot$Nc, na.rm = T)
    var_list["min_Ne",incr_multiplier,sp] <- min(data_tot$Ne) #, w=data_tot$Nc, na.rm = T)
    var_list["Ne/Nc",incr_multiplier,sp] <- mean(data_tot$Ne/data_tot$Nc, na.rm = T) #, w=data_tot$Nc, na.rm = T)
    var_list["mean_Ne_tempo",incr_multiplier,sp] <- mean(data_tot$Ne_tempo, na.rm = T)
    var_list["mean_htz",incr_multiplier,sp] <- mean(data_htz$HTZ, na.rm = T)
    var_list["quantile95_Ne", incr_multiplier, sp] <- quantile(data_tot$Ne, 0.95, na.rm = T)
    var_list["quantile5_Ne", incr_multiplier, sp] <- quantile(data_tot$Ne, 0.05, na.rm = T)
    var_list["quantile95_Ne_tempo", incr_multiplier, sp] <- quantile(data_tot$Ne_tempo, 0.95, na.rm = T)
    var_list["quantile5_Ne_tempo", incr_multiplier, sp] <- quantile(data_tot$Ne_tempo, 0.05, na.rm = T)
    var_list["quantile95_htz", incr_multiplier, sp] <- quantile(data_htz$HTZ, 0.95, na.rm = T)
    var_list["quantile5_htz", incr_multiplier, sp] <- quantile(data_htz$HTZ, 0.05, na.rm = T)
  }
}

save(var_list, file= file.path(configDir, "output", "Ne","by_multiplier/var_list.RData"))
var_list<-get(load(file= file.path(configDir, "output", "Ne","by_multiplier/var_list.RData")))

mu_vec <- c()
plot <- list()
coef_vec <- c()
for(sp in species){
  cat(sp)
  par(mfrow=c(1,2))
  df_var <- as.data.frame(var_list[,,sp])
  # colnames(df_var) <- factor(colnames(df_var), levels=c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2,3,4,5))
  df_var <- t(rbind("School_number_multiplier"=multipliers, df_var))
  df_var <- df_var[,order(colnames(df_var))]
  # coef <- (max(df_var[, ""], na.rm = T)/max(df_var[,""], na.rm = T)) #which(df_long_tot$age==1 & 
  # coef_vec <- c(coef_vec, coef)
  
  if((ncol(df_var[,])-1)*nrow(df_var[,]) == sum(is.na(df_var))){
    cat("table empty for", sp, "\n")
    next
  }
  
  multiplier = mean(df_var[,"Ne/Nc"], na.rm = T)/(max(df_var[,"mean_Ne"], na.rm = T))
  multiplier2 = mean(df_var[,"mean_htz"], na.rm = T)/mean(df_var[,"mean_Ne_tempo"], na.rm = T)
  
  p3 <- ggplot(df_var, aes(x=School_number_multiplier))+
    ylim(0, max(df_var[, "max_Ne"], na.rm = T))+
    geom_point(aes(y=mean_Ne), col="green", size=3)+
    geom_point(aes(y=median_Ne), col="blue", size=3, alpha = 0.5)+
    geom_point(aes(y=max_Ne), col="red", size=3, alpha=0.5)+
    geom_point(aes(y=min_Ne), col="red", size=3, alpha=0.5)+
    geom_errorbar(aes(ymin = quantile5_Ne, ymax = quantile95_Ne), col="red")+
    geom_point(aes(y=`Ne/Nc`/multiplier), col="black")+
    # geom_errorbar(aes(ymin = quantile5_Ne/multiplier, ymax = quantile95_Ne/multiplier))+
    # geom_errorbar(aes(ymin = quantile5_Ne_tempo, ymax = quantile95_Ne_tempo), col="blue")+
    scale_y_continuous(
      # Features of the first axis
      name = "Ne (green : mean, blue : median, red : min/max/IC)", #"mean temporal Ne (blue)",
      # Add a second axis and specify its features
      sec.axis = sec_axis(~.*multiplier, name=paste0("Ne/Nc (black)"))
    )+
    # theme(axis.line = element_line(),
    #       plot.margin = margin(10, 10, 10, 30))+
    labs(title = paste0("Mean, median, min, max and 95% confidence intervalle of Ne\nover ",  sp," cohorts, according to the \nnumber of school multiplier"))
  
  # plot[[sp]] <- wrap_elements(get_plot_component(p1, "ylab-l")) +
  #   wrap_elements(get_y_axis(p1)) +
  #   wrap_elements(get_plot_component(p2, "ylab-l")) +
  #   wrap_elements(get_y_axis(p2)) +
  #   p3+
  #   plot_layout(widths = c(3, 1, 3, 1, 40))
  
  plot[[sp]] <- p3 
  
  mu_vec <- c(mu_vec, multiplier)
  coef_vec <- c(coef_vec, multiplier2)
}


# STORE THE PLOTS
dir.create(paste0(configDir, "/output/Ne/PLOT_NeVSmultiplier"))
pdf(file = paste0(configDir, "/output/Ne/PLOT_NeVSmultiplier/plot_Ne_multiplier.pdf"), width = 10)
# pdf(file = paste0(paste0(pathNCDF, "_FINAL_plot/", "tot", ".pdf")), width=25, height=15)
par(mfrow=c(1,1))
for (k in 1:length(plot)){
  # multiplier2 = coef_vec[k]
  multiplier =   mu_vec[k]
  print(plot[[k]])
}
dev.off()



# mu_vec <- c()
# plot <- list()
# coef_vec <- c()
# for(sp in species){
#   cat(sp)
#   par(mfrow=c(1,2))
#   df_var <- as.data.frame(var_list[,,sp])
#   colnames(df_var) <- factor(colnames(df_var), levels=c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2,3,4,5))
#   df_var <- t(rbind("School_number_multiplier"=multipliers, df_var))
#   df_var <- df_var[,order(colnames(df_var))]
#   # coef <- (max(df_var[, ""], na.rm = T)/max(df_var[,""], na.rm = T)) #which(df_long_tot$age==1 & 
#   # coef_vec <- c(coef_vec, coef)
#   
#   if((ncol(df_var[,])-1)*nrow(df_var[,]) == sum(is.na(df_var))){
#     cat("table empty for", sp, "\n")
#     next
#   }
#   
#   multiplier = max(df_var[,"mean_Ne"], na.rm = T)/max(df_var[,"mean_htz"], na.rm = T)
#   multiplier2 = max(df_var[,"mean_htz"], na.rm = T)/max(df_var[,"mean_Ne_tempo"], na.rm = T)
#   
#   
#   p1 <- ggplot(df_var, aes(x=School_number_multiplier, y=mean_Ne)) + 
#     scale_y_continuous(
#     # Features of the first axis
#     name = "mean Ne weighted by Nc (black) and IC of Ne")
#   p2 <- ggplot(df_var, aes(x=School_number_multiplier, y=mean_htz))+ 
#     scale_y_continuous(
#       # Features of the first axis
#       name = "mean Htz (red)")
#   p3 <- ggplot(df_var, aes(x=School_number_multiplier, y=mean_htz))+
#     # geom_point(aes(y=mean_Ne_tempo), col="blue", size=3)+
#     geom_point(aes(y=mean_htz), col="red", size=3)+
#     geom_errorbar(aes(ymin = quantile5_htz, ymax = quantile95_htz), col="red")+
#     geom_point(aes(y=mean_Ne/multiplier))+
#     geom_errorbar(aes(ymin = quantile5_Ne/multiplier, ymax = quantile95_Ne/multiplier))+
#     # geom_errorbar(aes(ymin = quantile5_Ne_tempo, ymax = quantile95_Ne_tempo), col="blue")+
#     scale_y_continuous(
#       # Features of the first axis
#       name = "mean Htz (red)", #"mean temporal Ne (blue)",
#       # Add a second axis and specify its features
#       sec.axis = sec_axis(~.*multiplier, name=paste0("Ne (black)"))
#       )+
#     theme(axis.line = element_line(),
#           plot.margin = margin(10, 10, 10, 30))+
#     labs(title = paste0("Mean and 95% confidence intervalle of Ne and Htz \nover ",  sp," cohorts, according to the \nnumber of school multiplier"))
#   
#   # plot[[sp]] <- wrap_elements(get_plot_component(p1, "ylab-l")) +
#   #   wrap_elements(get_y_axis(p1)) +
#   #   wrap_elements(get_plot_component(p2, "ylab-l")) +
#   #   wrap_elements(get_y_axis(p2)) +
#   #   p3+
#   #   plot_layout(widths = c(3, 1, 3, 1, 40))
#   
#   plot[[sp]] <- p3 
#   
#   mu_vec <- c(mu_vec, multiplier)
#   coef_vec <- c(coef_vec, multiplier2)
# }





