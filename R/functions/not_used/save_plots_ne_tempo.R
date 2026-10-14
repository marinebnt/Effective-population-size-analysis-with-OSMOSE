# monitoring Ne variability --------------------
# MB
# 11/06/25

library(tidyr)
library(dplyr)


save_plots_ne_tempo <- function(data_Ne_t, data_Ne, config, outputDir_Ne, outputDir_Simu, seeding_y, print=FALSE){
  

  
  # Plot 1 ------------------------------------------------
  # variability_over_time
  
  Ne_var <- list()
  k=0
  for(sp in species){
    
    # Ne dataset
    df      <- as.data.frame(data_Ne_t[sp])
    df$Time <- rownames(df)
    
    k=k+1
    
    a <- var(df[, sp], na.rm = T)
    Ne_var[[k]] <- a
    
    cat(sp, ": var :", a, "\n")
  }
  
  # Plot 2 -----------------------------------------------------
  # Plotting the cohorts Ne (2 plots combined)
  
  # Assuming your data is in a data frame called `df`
  # Rows: time steps, Columns: cohort abundances
  # Example: df <- data.frame(Cohort1 = ..., Cohort2 = ..., ...)
  
  # 1- Ne // Nc
  
  par(mfrow = c(3,3))
  k=0
  Nc_Ne <- list()
  multi_vec <- c()
  for(sp in species){
    cat(sp, "\n")
   
    # Nc data 
    data <- config[["abundance"]][,,1][,sp]

    # Ne data
    df      <- data.frame(data_Ne_t[sp])
    df$Time <- as.numeric(rownames(df))
    colnames(df) <- c("Ne", "Time")
    
    # Nc in Ne shape 
    Nc <- as.data.frame(data[which(names(data) %in% df$Time)])
    
    # bind datasets
    df_long_tot <- cbind(df, Nc)
    colnames(df_long_tot)[3] <- "Nc"
    if(any(is.na(df_long_tot$Ne))){
      df_long_tot[-which(is.na(df_long_tot$Ne)), ] -> df_long_tot
    }
    if(!length(df_long_tot$Time > seeding_y)>0){
      cat("no data in the required time and age categories\n")
      next
    }
    multip <- (max(df_long_tot$Nc[which(df_long_tot$Time > seeding_y)], na.rm = T)/max(df_long_tot$Ne, na.rm = T))
    na_Nc <- which(is.na(df_long_tot$Nc))
    df_long_tot$Nc[na_Nc] <- 0
    
    # Plot using ggplot2
    
    a <- ggplot(df_long_tot[which(df_long_tot$Time > seeding_y),], aes(x = Time)) +
      geom_point(aes(y = Ne), col="red", pch=24) +
      geom_smooth(method = loess, aes(y = Ne), lwd=0.7, col="red")+
      geom_point(aes(y = Nc/multip), size=3) +
      geom_smooth(method = loess, aes(y=Nc/multip), lwd=0.7)+
      labs(title = paste0(sp, ": Comparing Ne and Nc"), x = "Time Step", y = "Ne") +
      theme_minimal()+ 
      theme(legend.position="none")+
      scale_y_continuous(
        # Features of the first axis
        name = "Ne (red triangles)",
        # Add a second axis and specify its features
        sec.axis = sec_axis(~.*multip, name="Total species Nc (black)")
      ) 
    par(mfrow = c(3,3))
    
    if (print) {
      print(a)
      readline(prompt="Press [enter] to continue")
    }
    k=k+1
    multi_vec <- c(multi_vec, multip)
    Nc_Ne[[k]] <- a
  }
  
  
  # 2- Ne // HTZ
  
  par(mfrow = c(3,3))
  k=0
  HTZ_Ne <- list()
  coef_vec <- c()
  for(sp in species){
    cat(sp, "\n")
    
    # HTZ data 
    if (!file.exists(paste0(outputDir_Ne, "/htz", sp, ".RData"))){
      cat("saving htz data\n")
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
    
    # Ne data
    df      <- data.frame(data_Ne_t[sp])
    df$Time <- as.numeric(rownames(df))
    colnames(df) <- c("Ne", "Time")
    
    # Nc in Ne shape 
    data_htz <- as.data.frame(data_htz[which(data_htz$Time %in% df$Time),])
    
    # bind datasets
    df_long_tot <- dplyr::left_join(df, data_htz, by = "Time")
    if(any(is.na(df_long_tot$Ne))){
      df_long_tot[-which(is.na(df_long_tot$Ne)), ] -> df_long_tot
    }
    if(!length(df_long_tot$Time > seeding_y)>0){
      cat("no data in the required time and age categories\n")
      next
    }
    coef <- (max(df_long_tot$HTZ, na.rm = T)/max(df_long_tot$Ne[which(df_long_tot$Time > seeding_y)], na.rm = T)) #which(df_long_tot$age==1 & 
    lm_htz <- lm(data_htz$HTZ~data_htz$Time) 
    slope <- lm_htz$coefficients[2]
  
  # na_HTZ <- which(is.na(df_long_tot$HTZ))
  # df_long_tot$HTZ[na_HTZ] <- 0
  
  # Plot using ggplot2
  a1 <- ggplot(df_long_tot[which(df_long_tot$Time > seeding_y),], aes(x = Time)) + #which(df_long_tot$age==1 & 
    geom_point(aes(x=Time, y = Ne), size=2, col="red") +
    geom_smooth(method = loess, aes(y = Ne), lwd=0.7, col="red")+
    geom_point(aes(y = HTZ/coef), lwd = 0.75)+
    geom_smooth(method = loess, aes(y = HTZ/coef), lwd=0.7)+
    labs(title = "Comparing observed HTZ and Ne" , x = "Time Step", y = "Ne") +
    theme_minimal()+ 
    theme(legend.position="none")+
    scale_y_continuous(
      # Features of the first axis
      name = "Ne (color by cohort) (red)",
      # Add a second axis and specify its features
      sec.axis = sec_axis(~.*coef, name=paste0("Mean species HTZ (black), slope lm model (blue)= ", slope))
    ) 
  par(mfrow = c(3,3))
  
  if (print) {
    print(a1)
    readline(prompt="Press [enter] to continue")
  }
  k=k+1
  coef_vec <- c(coef_vec, coef)
  HTZ_Ne[[k]] <- a1
  }
  
  # Save as .pdf
  pdf(file = paste0(outputDir_Ne, "/plot_Ne_shorter.pdf"), width = 10)
  # pdf(file = paste0(paste0(pathNCDF, "_FINAL_plot/", "tot", ".pdf")), width=25, height=15)
  par(mfrow=c(1,1))
  for (k in 1:length(Nc_Ne)) {
    multip <- multi_vec[k]
    coef <- coef_vec[k]
    print(cowplot::plot_grid(Nc_Ne[[k]], HTZ_Ne[[k]], ncol=2, nrow=1))
  }
  dev.off()  
  
  # # Plot 3 -----------------------------------------------------
  # # Plotting the cohorts Ne/Nc
  # 
  # # Assuming your data is in a data frame called `df`
  # # Rows: time steps, Columns: cohort abundances
  # # Example: df <- data.frame(Cohort1 = ..., Cohort2 = ..., ...)
  # 
  # par(mfrow = c(3,3))
  # k=0
  # NebyNc <- list()
  # multi_vec_nebync <- c()
  # for(sp in species){
  #   # Nc data 
  #   k=k+1
  #   data <- config[["abundance"]][,,1][,sp]
  #   
  #   # Ne data
  #   df      <- data.frame(data_Ne_t[sp])
  #   df$Time <- as.numeric(rownames(df))
  #   colnames(df) <- c("Ne", "Time")
  #   
  #   # Nc in Ne shape 
  #   Nc <- as.data.frame(data[which(names(data) %in% df$Time)])
  #   
  #   # bind datasets
  #   df_long_tot <- cbind(df, Nc)
  #   colnames(df_long_tot)[3] <- "Nc"
  #   if(any(which(df_long_tot$Ne == 0))){
  #     df_long_tot[-which(df_long_tot$Ne == 0),] -> df_long_tot
  #   }
  #   if(any(is.na(df_long_tot$Ne))){
  #     df_long_tot[-which(is.na(df_long_tot$Ne)), ] -> df_long_tot
  #   }
  #   df_long_tot$NebyNc <- df_long_tot$Ne/df_long_tot$Nc
  #   multip <- (max(df_long_tot$Nc, na.rm = T)/max(df_long_tot$Ne, na.rm = T))
  #   
  #   
  #   # Plot using ggplot2
  #   a <- ggplot(df_long_tot[-1,][which(df_long_tot$Time > seeding_y),], aes(x = Time)) +
  #     geom_point(aes(y = NebyNc), size=2) +
  #     labs(title = paste0(sp, ": Ne/Nc"), x = "Time Step", y = "Ne/Nc") +
  #     theme_minimal()+ 
  #     theme(legend.position="none")+
  #     scale_y_continuous(
  #       # Features of the first axis
  #       name = "Ne/Nc (color by cohort)",
  #       # Add a second axis and specify its features
  #       sec.axis = sec_axis(~.*multip, name="Total species Nc (black)")
  #     ) 
  #   par(mfrow = c(3,3))
  #   
  #   if (print) {
  #     print(a)
  #     readline(prompt="Press [enter] to continue")
  #   }
  #   
  #   multi_vec_nebync <- c(multi_vec_nebync, multip)
  #   NebyNc[[k]] <- a
  # }
  # 
  # # STORE THE PLOTS
  # pdf(file = paste0(outputDir_Ne, "/tempo_NebyNc.pdf"), width = 10)
  # # pdf(file = paste0(paste0(pathNCDF, "_FINAL_plot/", "tot", ".pdf")), width=25, height=15)
  # par(mfrow=c(1,1))
  # for (k in 1:length(NebyNc)) {
  #   multip <- multi_vec_nebync[k]
  #   print(NebyNc[[k]])
  # }
  # dev.off()
  
  #* Plot 4 -----------------------------------------------------
  #* 
  #* Weighted plots
  #* 
  #* 1- Ne Weighted by school nb
  
  par(mfrow = c(3,3))
  k=0
  NebySchool <- list()
  multi_vec_nebyschool <- c()
  for(sp in species){
    k=k+1
    cohort_time_corr <- data_Ne[["cohort_time_step_corresponding"]]
    co_to_include <- which(colSums(data_Ne[["longevity"]][[sp]][[paste0("rep_", rep)]] )>0)
    
    # Nc data 
    data <- config[["abundance"]][,,1][,sp]
    
    # Ne data
    df      <- data.frame(data_Ne_t[sp])
    df$Time <- as.numeric(rownames(df))
    colnames(df) <- c("Ne", "Time")
    
    # Nc in Ne shape 
    Nc <- as.data.frame(data[which(names(data) %in% df$Time)])
    
    # bind datasets
    df_long_tot <- cbind(df, Nc)
    colnames(df_long_tot) <- c("Ne", "Time", "Nc")
    
    # school number by cohort 
    df_ab<- as.data.frame(data_Ne[["schoolNumber"]][[sp]])
    # Add time as a column
    df_ab$age <- as.numeric(rownames(df_ab))
    colnames(df_ab) <- stringr::str_remove(colnames(df_ab), pattern = paste0(".rep_", rep))
    # Reshape data to long format
    df_ab_long_i <- pivot_longer(df_ab, cols = -age, names_to = "Cohort", values_to = "nb_schools")
    
    # match the school number to each time step
    for (co in 0:c(ncohort-3)){ # removing the last 3 cohorts that are incomplete
      ID_co <- which(cohort_time_corr==co, arr.ind = T)[,1]
      df_ab_long_i[grep(paste0("\\bcohort_", co, "\\b"), df_ab_long_i$Cohort)[1:length(ID_co)], "Time"] <- ID_co
    }
    schools = aggregate(data=df_ab_long_i, nb_schools~Time, sum)
    schools = as.data.frame(schools[which(schools$Time %in% df$Time),])
    
    # Ne + nb school data
    df_long_tot <- dplyr::full_join(df_long_tot, schools)
    
    # weighted mean per time step
    mat_ne_plot <- NULL
    for (ti in df_long_tot$Time){
      mat_ne_plot$Time[ti] <- ti
      mat_ne_plot$Ne_wei_school[ti] <- weighted.mean(x=df_long_tot$Ne[which(df_long_tot$Time == ti)], w=df_long_tot$nb_schools[which(df_long_tot$Time == ti)], na.rm = T)
      mat_ne_plot$Ne_wei_school_byNc[ti] <- mat_ne_plot$Ne_wei_school[ti]/df_long_tot$Nc[which(df_long_tot$Time == ti)]
    }
    mat_ne_plot <- as.data.frame(mat_ne_plot)
    mat_ne_plot <- mat_ne_plot[-which(is.na(mat_ne_plot$Time)),]
    
    # clean dataset
    if(any(which(df_long_tot$Ne == 0))){
      df_long_tot[-which(df_long_tot$Ne == 0),] -> df_long_tot
    }
    if(any(is.na(df_long_tot$Ne))){
      df_long_tot[-which(is.na(df_long_tot$Ne)), ] -> df_long_tot
    }
    if(!length(which(df_long_tot$Time > seeding_y))>0){
      cat("no data in the required time and age categories\n")
      next
    }    
    
    # df_long_tot$NebyNc <- df_long_tot$Ne/df_long_tot$Nc
    multip_sc <- (max(mat_ne_plot$Ne_wei_school_byNc, na.rm = T)/max(mat_ne_plot$Ne_wei_school, na.rm = T))
    
    
    # Plot using ggplot2
    a <- ggplot(mat_ne_plot, aes(x = Time)) +
      geom_line(aes(y = Ne_wei_school)) +
      geom_line(aes(y = Ne_wei_school_byNc/multip_sc), col="orange", lty="longdash", size=0.75)+
      labs(title = paste0(sp, ": Mean Ne weighted by school \nnumber at each time step \nand weighted Ne/Nc"), x = "Time Step", y = "Ne_wei") +
      theme_minimal()+ 
      theme(legend.position="none")+
      scale_y_continuous(
        # Features of the first axis
        name = "weighted N (black)",
        # Add a second axis and specify its features
        sec.axis = sec_axis(~.*multip_sc, name="weighted Ne by school nb by Nc (orange)")
      ) 
    par(mfrow = c(3,3))
    
    if (print) {
      print(a)
      readline(prompt="Press [enter] to continue")
    }
    
    multi_vec_nebyschool <- c(multi_vec_nebyschool, multip_sc)
    NebySchool[[k]] <- a
  }
  
  
  #* 2- Ne Weighted by school Nc
  #* 
  par(mfrow = c(3,3))
  k=0
  NebyNc <- list()
  multi_vec_nebync <- c()
  for(sp in species){
    k=k+1
    cohort_time_corr <- data_Ne[["cohort_time_step_corresponding"]]
    co_to_include <- which(colSums(data_Ne[["longevity"]][[sp]][[paste0("rep_", rep)]] )>0)
    
    # Nc data 
    data <- config[["abundance"]][,,1][,sp]
    
    # Ne data
    df      <- data.frame(data_Ne_t[sp])
    df$Time <- as.numeric(rownames(df))
    colnames(df) <- c("Ne", "Time")
    
    # Nc in Ne shape 
    Nc <- as.data.frame(data[which(names(data) %in% df$Time)])
    
    # bind datasets
    df_long_tot <- cbind(df, Nc)
    colnames(df_long_tot) <- c("Ne", "Time", "Nc")
  
    
    # abundance by cohort 
    df_ab<- as.data.frame(data_Ne[["abundanceByAge"]][[sp]][,,rep])
    # Add time as a column
    df_ab$age <- as.numeric(rownames(df_ab))
    # Reshape data to long format
    df_ab_long_i <- pivot_longer(df_ab, cols = -age, names_to = "Cohort", values_to = "abund_scaled")
    
    # match the abund number to each time step
    for (co in 0:c(ncohort-3)){ # removing the last 3 cohorts that are incomplete
      ID_co <- which(cohort_time_corr==co, arr.ind = T)[,1]
      df_ab_long_i[grep(paste0("\\bcohort_", co, "\\b"), df_ab_long_i$Cohort)[1:length(ID_co)], "Time"] <- ID_co
    }
    abund = aggregate(data=df_ab_long_i, abund_scaled~Time, sum)
    abund = as.data.frame(abund[which(abund$Time %in% df$Time),])
    
    
    # Ne + abund data
    df_long_tot <- dplyr::full_join(df_long_tot, abund)
    
    # weighted mean per time step
    mat_ne_plot <- NULL
    for (ti in df_long_tot$Time){
      mat_ne_plot$Time[ti] <- ti
      mat_ne_plot$Ne_wei_Nc[ti] <- weighted.mean(x=df_long_tot$Ne[which(df_long_tot$Time == ti)], w=df_long_tot$abund_scaled[which(df_long_tot$Time == ti)], na.rm = T)
      mat_ne_plot$Ne_wei_Nc_byNc[ti] <- mat_ne_plot$Ne_wei_Nc[ti]/df_long_tot$Nc[which(df_long_tot$Time == ti)]
    }
    mat_ne_plot <- as.data.frame(mat_ne_plot)
    mat_ne_plot <- mat_ne_plot[-which(is.na(mat_ne_plot$Time)),]
    
    # clean dataset
    if(any(which(df_long_tot$Ne == 0))){
      df_long_tot[-which(df_long_tot$Ne == 0),] -> df_long_tot
    }
    if(any(is.na(df_long_tot$Ne))){
      df_long_tot[-which(is.na(df_long_tot$Ne)), ] -> df_long_tot
    }
    if(!length(df_long_tot$Time > seeding_y)>0){
      cat("no data in the required time and age categories\n")
      next
    }    
    
    # df_long_tot$NebyNc <- df_long_tot$Ne/df_long_tot$Nc
    multip_nc <- (max(mat_ne_plot$Ne_wei_Nc_byNc, na.rm = T)/max(mat_ne_plot$Ne_wei_Nc, na.rm = T))
    
    
    # Plot using ggplot2
    a <- ggplot(mat_ne_plot, aes(x = Time)) +
      geom_line(aes(y = Ne_wei_Nc)) +
      geom_line(aes(y = Ne_wei_Nc_byNc/multip_nc), col="orange", lty="longdash", size=0.75)+
      labs(title = paste0(sp, ": Mean Ne weighted by cohort \nabundance at each time step and \nweighted Ne/Nc"), x = "Time Step", y = "Ne_wei") +
      theme_minimal()+ 
      theme(legend.position="none")+
      scale_y_continuous(
        # Features of the first axis
        name = "weighted N (black)",
        # Add a second axis and specify its features
        sec.axis = sec_axis(~.*multip_nc, name="weighted Ne by school nb by Nc (orange)")
      ) 
    par(mfrow = c(3,3))
    
    if (print) {
      print(a)
      readline(prompt="Press [enter] to continue")
    }
    
    multi_vec_nebync <- c(multi_vec_nebync, multip_nc)
    NebyNc[[k]] <- a
  }
  
  
  # STORE THE PLOTS
  pdf(file = paste0(outputDir_Ne, "/tempo_plot_Neweighted.pdf"), width = 10)
  # pdf(file = paste0(paste0(pathNCDF, "_FINAL_plot/", "tot", ".pdf")), width=25, height=15)
  par(mfrow=c(1,1))
  for (k in 1:length(NebyNc)) {
    multip_nc <- multi_vec_nebync[k]
    multip_sc <- multi_vec_nebyschool[k]
    print(cowplot::plot_grid(NebySchool[[k]], NebyNc[[k]], ncol=2, nrow=1))
  }
  dev.off()
  
}
