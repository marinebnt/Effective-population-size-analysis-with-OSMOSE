library(legendry)
library(ggpubr)

function_plot_trait_age_or_time <- function(SP, age_or_cohort, pathout, Fmax_constrain = NULL, Fstart_constrain = NULL, 
                                            list_traits = NULL, data_complete, relative_to_ref = FALSE){
  if (age_or_cohort == "age" & !grepl(pattern="age", pathout)) {
    cat("dataset has no age data in it")
    stop()
  }
  if (age_or_cohort == "cohort" & grepl(pattern="age", pathout)) {
    cat("dataset has age data in it and shouldn't have")
    stop()
  }
  
  plot_trait = list()
  scale_v = list()
  
  # clean dataset
  if (age_or_cohort == "age" & all(!grepl(pattern="age", colnames(data_complete)))) {
    colnames(data_complete)[which(colnames(data_complete) == "time")] <- "age"
  }
  if (any(unique(data_complete$TRAIT) %in% c("mean_of_Z_per_age", "mean_of_Dis_per_age"))) data_complete <- data_complete[-which((data_complete$TRAIT) %in% c("mean_of_Z_per_age", "mean_of_Dis_per_age")),]
  data_complete[,age_or_cohort] <- as.numeric(as.matrix(data_complete[, age_or_cohort]))
  if (age_or_cohort == "age")    data_complete <- subset(data_complete, age <= 30 & age >=1)
  if (age_or_cohort == "cohort") data_complete <- subset(data_complete, cohort >= (min(data_complete[, age_or_cohort])+19) & cohort <= (max(data_complete[,age_or_cohort])-20))
  # add fecunditybyage dimension
  if (!(any(grepl("fec",  data_complete$TRAIT)))) {
    egg = NULL
    abund = NULL
    egg = pivot_wider(data_complete[grep(data_complete$TRAIT, pattern="egg"),], values_from = "value", names_from = "TRAIT")
    abund = pivot_wider(data_complete[c(grep(data_complete$TRAIT, pattern="abundance")),], values_from = "value", names_from = "TRAIT")
    fe_data = full_join(egg, abund)
    fe_data[, "fecundityByAge"] = fe_data[,c(grep(colnames(fe_data), pattern="egg"))]/fe_data[,c(grep(colnames(fe_data), pattern="abundance"))]
    fe_data <- pivot_longer(fe_data[, -c(grep(colnames(fe_data), pattern="egg"), grep(colnames(fe_data), pattern="abund"))], names_to = "TRAIT", values_to = "value", cols = "fecundityByAge")
    data_complete <- rbind(data_complete, fe_data)
  }
  # select the traits to plot 
  TRAITS_to_plot = unique(data_complete$TRAIT) # plot all the traits of the dataset
  other_TRAITS_to_plot = NULL
  if (!is.null(list_traits)){ # or only the ones selected
    TRAITS_to_plot = list_traits[[1]]
    TRAITS_to_plot2 = TRAITS_to_plot
    if (any(grepl(TRAITS_to_plot, pattern="log"))) TRAITS_to_plot2 <- stringr::str_remove(pattern = "log", TRAITS_to_plot) 
    if (any(grepl(TRAITS_to_plot, pattern="exp")) & !any(grepl(TRAITS_to_plot, pattern="exp-"))) TRAITS_to_plot2 <- stringr::str_remove(pattern = "exp", TRAITS_to_plot) 
    if (any(grepl(TRAITS_to_plot, pattern="exp-"))) TRAITS_to_plot2 <- stringr::str_remove(pattern = "exp-", TRAITS_to_plot) 
    if (!any(grepl(pattern=TRAITS_to_plot2, unique(data_complete$TRAIT))) & (TRAITS_to_plot != "fecundityByAge")){
      cat("traits to plot are not in the dataset. Dataset contains : \n", unique(data_complete$TRAIT))
      stop()
    }
    if (length(list_traits)>1) {            # one trait at a time, or multiple
      other_TRAITS_to_plot = list_traits[-1]
    }
  }
  
  for (sp in SP){
    
  # PLOT the mean and the std deviation 
    k=0
    other_TRAIT = NULL
    for (TRAIT in TRAITS_to_plot){
      k=k+1
      for (sp_plotted in sp){
        traits_df <- c(TRAIT, other_TRAITS_to_plot[[k]])
        rm("id_log", "id_exp_neg", "id_exp_pos")
        if (any(grepl(c(TRAIT, other_TRAITS_to_plot[[k]]), pattern="log"))) {
          traits_df <- stringr::str_remove(pattern = "log", c(TRAIT, other_TRAITS_to_plot[[k]])) 
          id_log = grep(c(TRAIT, other_TRAITS_to_plot[[k]]), pattern="log")
        }
        if (any(grepl(c(TRAIT, other_TRAITS_to_plot[[k]]), pattern="exp-"))) {
          traits_df <- stringr::str_remove(pattern = "exp-", traits_df)
          id_exp_neg = which(grepl(c(TRAIT, other_TRAITS_to_plot[[k]]), pattern = "exp-"))
        }
        if (any(grepl(c(TRAIT, other_TRAITS_to_plot[[k]]), pattern="exp")) & !exists("id_exp_neg")) {
          traits_df <- stringr::str_remove(pattern = "exp", traits_df)
          id_exp_pos = which(grepl(c(TRAIT, other_TRAITS_to_plot[[k]]), pattern = "exp") &
                               !grepl(c(TRAIT, other_TRAITS_to_plot[[k]]), pattern = "exp-"))
        }
        data_tot_summary_tr    <- data_complete[which(data_complete$TRAIT %in% traits_df),]
        if (exists("id_log")){
          id_log <- which(data_complete$TRAIT %in% traits_df[id_log])
          log_val  <- log(as.matrix(data_tot_summary_tr[id_log, "value"]))
          log_val[which(is.infinite(log_val))] <- 0
          df_log <- data_tot_summary_tr[id_log,]
          data_tot_summary_tr <- data_tot_summary_tr[-id_log, ]
          df_log[, "value"] <- log_val
          data_tot_summary_tr <- rbind(data_tot_summary_tr, df_log)
        }
        if (exists("id_exp_neg")){
          id_exp_neg <- which(data_complete$TRAIT %in% traits_df[id_exp_neg])
          exp_neg_val  <- exp(-as.matrix(data_tot_summary_tr[id_exp_neg, "value"]))
          df_exp_neg <- data_tot_summary_tr[id_exp_neg,]
          data_tot_summary_tr <- data_tot_summary_tr[-id_exp_neg, ]
          df_exp_neg[, "value"] <- exp_neg_val
          data_tot_summary_tr <- rbind(data_tot_summary_tr, df_exp_neg)
        }
        if (exists("id_exp_pos")){
          id_exp_pos <- which(data_complete$TRAIT %in% traits_df[id_exp_pos])
          exp_pos_val  <- exp(as.matrix(data_tot_summary_tr[id_exp_pos, "value"]))
          df_exp_pos <- data_tot_summary_tr[id_exp_pos,]
          data_tot_summary_tr <- data_tot_summary_tr[-id_exp_pos, ]
          df_exp_pos[, "value"] <- exp_pos_val
          data_tot_summary_tr <- rbind(data_tot_summary_tr, df_exp_pos)
        }
        data_tot_summary_tr_sp <- data_tot_summary_tr[which(data_tot_summary_tr$Species==sp_plotted),] 
        data_tot_summary_tr_sp$value <- as.numeric(data_tot_summary_tr_sp$value)
        data_tot_summary_tr_sp <- data_tot_summary_tr_sp %>%
          group_by(get(age_or_cohort), Species, Fish_press, TRAIT) %>%
          summarise(mean = mean(value, na.rm = TRUE),
                    sd   = sd(value, na.rm = TRUE))
        colnames(data_tot_summary_tr_sp)[which(colnames(data_tot_summary_tr_sp) == "get(age_or_cohort)")] <- age_or_cohort
        # reorder Fmax/Fstart factor for prettier plots
        data_tot_summary_tr_sp$Fish_press    <- factor(as.factor(data_tot_summary_tr_sp$Fish_press), levels = ideal_order_factor)
        
        # select the Fishing pressures I want to keep
        if (!is.null(Fmax_constrain) | !is.null(Fstart_constrain)){
          F_keep = NULL
          F_keep = NULL
          if (!is.null(Fmax_constrain))    F_keep = (grep(data_tot_summary_tr_sp$Fish_press, pattern=paste0("Fmax=", Fmax_constrain, " Fstart")))
          if (!is.null(Fstart_constrain))  F_keep = c(F_keep, grep(data_tot_summary_tr_sp$Fish_press, pattern=paste0("Fstart=", Fstart_constrain, "$")))
          if(any(duplicated(F_keep))) F_keep <- F_keep[duplicated(F_keep)]
          data_tot_summary_tr_sp <- data_tot_summary_tr_sp[F_keep, ]
        }
        reference <- subset((data_tot_summary_tr_sp), Fish_press == paste0("Fmax=", 1, " Fstart=", 1))
        reference <- reference[, -which(colnames(reference) == 'Fish_press')]
        if (relative_to_ref && nrow(reference) > 0) {
          reference <- reference %>%
            dplyr::select(all_of(age_or_cohort), Species, TRAIT, mean) %>%
            dplyr::rename(ref_mean = mean)
          data_tot_summary_tr_sp <- dplyr::left_join(
            data_tot_summary_tr_sp,
            reference,
            by = c(age_or_cohort, "Species", "TRAIT")
          )
          data_tot_summary_tr_sp$mean <- 100*(data_tot_summary_tr_sp$mean - data_tot_summary_tr_sp$ref_mean)/data_tot_summary_tr_sp$ref_mean
          data_tot_summary_tr_sp$ref_mean <- NULL
        }
        
        if(nrow(data_tot_summary_tr_sp)==0){
          plot_trait[[TRAIT]][[sp_plotted]] <- ggplot(data = data_tot_summary_tr_sp, aes(x = age_or_cohort, group = Fish_press)) 
          next
        }
        
        x_axis <- age_or_cohort
        title = paste0("By cohort ", sp_plotted, ifelse(relative_to_ref==T, " % relative (to 1,1) changes of ", " ranges of "), 
                       paste(c(TRAIT, unlist(other_TRAITS_to_plot)), collapse=" "), "\nF varies for ", sp)
        vec_linetype <- c("solid", "dotted", "dashed", "D3")[1:length(unique(data_tot_summary_tr_sp$TRAIT))]
        # vec_pointtype <- c("solid", "dotted", "dashed", "D3")[1:length(unique(data_tot_summary_tr_sp$TRAIT))]
        data_tot_summary_tr_sp$vec_linetype = sapply(data_tot_summary_tr_sp$TRAIT, function(x) vec_linetype[which(unique(data_tot_summary_tr_sp$TRAIT) == x)])
        data_tot_summary_tr_sp$mean  = as.numeric(data_tot_summary_tr_sp$mean)
        data_tot_summary_tr_sp$sd    = as.numeric(data_tot_summary_tr_sp$sd)
        data_tot_summary_tr_sp$TRAIT = as.factor(data_tot_summary_tr_sp$TRAIT)
        data_tot_summary_tr_sp$Fish_press = as.factor(data_tot_summary_tr_sp$Fish_press)
        NAN_remove <- sapply(unique(data_tot_summary_tr_sp$TRAIT), FUN = function(x) all(is.nan(as.matrix(data_tot_summary_tr_sp[which(data_tot_summary_tr_sp$TRAIT == x),"mean"]))))
        NAN_remove <- sapply(unique(data_tot_summary_tr_sp$TRAIT), FUN = function(x) all(is.nan(as.matrix(data_tot_summary_tr_sp[which(data_tot_summary_tr_sp$TRAIT == x),"mean"]))))
        if(any(NAN_remove)) {
          data_tot_summary_tr_sp <- data_tot_summary_tr_sp[-which(data_tot_summary_tr_sp$TRAIT == unique(data_tot_summary_tr_sp$TRAIT)[which(NAN_remove)]),]
          title <- paste0(title, "\n removed ", unique(data_tot_summary_tr_sp$TRAIT)[which(NAN_remove)], " due to NAN")
        }
        
        
        scale = 1
        var_scale_regular = TRAIT
        name_second_axis = ""
        if(!is.null(other_TRAITS_to_plot)){
          data_key <- data_tot_summary_tr_sp[which(data_tot_summary_tr_sp$TRAIT %in% as.matrix(other_TRAITS_to_plot[[1]])),]
          key <- key_auto()(
            scale_y_continuous(
              limits = range(data_key$mean, na.rm=T),
              breaks = c(quantile(data_key$mean, probs = 0.25, na.rm = T), quantile(data_key$mean, probs = 0.5, na.rm = T), 
                         quantile(data_key$mean, probs = 0.75, na.rm = T)),
              labels = \(x) x / 1000
            )
          )
          for_scale <- data_tot_summary_tr_sp[!is.infinite(data_tot_summary_tr_sp$mean), ]
          mean_var  <- aggregate(for_scale[, "mean"], list(for_scale$TRAIT), FUN = function(x) mean(x, na.rm=T))
          median_var<- aggregate(for_scale[, "mean"], list(for_scale$TRAIT), median)
          var_scale_regular = mean_var[which(mean_var[,2] > 0.8*(max(mean_var[,2]))),1]
          scale = mean(mean_var[-which.max(mean_var[,2]),"mean" ], na.rm=T)/
            quantile(mean(mean_var[which.max(mean_var[,2]),"mean" ], na.rm=T), 0.95)
          data_tot_summary_tr_sp[-which(data_tot_summary_tr_sp$TRAIT %in%  var_scale_regular), "mean"] <- 
            data_tot_summary_tr_sp[-which(data_tot_summary_tr_sp$TRAIT %in% var_scale_regular), "mean"]/scale
          name_second_axis = paste(mean_var[-which(mean_var[,2] > 0.8*(max(mean_var[,2]))),1], collapse = "\n")
        }
        scale <- ifelse(scale ==0, 1, scale)
        scale_v[[TRAIT]][[sp]] <- scale

        if (length(unique(data_tot_summary_tr_sp$TRAIT))>3) {
          col_fac = "TRAIT"
          shape_fac = "Fish_press"
        } else {
          col_fac = "Fish_press"
          shape_fac = "TRAIT"
        }
        
        
        plot_trait[[TRAIT]][[sp]] <- ggplot(data = data_tot_summary_tr_sp, aes(x = .data[[age_or_cohort]], 
                                                                               group=interaction(Fish_press, TRAIT))) +
          # geom_ribbon(aes(ymin = mean-sd, ymax = mean + sd),
            # alpha = 0.5, fill = "grey70")+
          geom_line(aes(y = mean, col=.data[[col_fac]], linetype=.data[[shape_fac]]), lwd=1, alpha=0.6) +
          geom_point(aes(y = mean, col=.data[[col_fac]], pch=.data[[shape_fac]]), alpha=0.5, col="black") + #linetype=data_tot_summary_tr_sp$TRAIT
          scale_color_brewer(palette = "Spectral", direction = -1) +
          theme_bw()+
          ylab(label = paste(var_scale_regular, collapse="\n"))+
          xlab(label = x_axis)+
          labs(title = title)+
          # facet_wrap(~ Fish_press)+ 
          # scale_linetype_manual(#values = data_tot_summary_tr_sp$vec_linetype, 
          #   name="Functional trait")+ 
          # scale_shape_manual(#(values = data_tot_summary_tr_sp$vec_pointtype,
          #   name="Functional trait")+ 
          theme(legend.key.width = unit(2,"cm"))+
          scale_y_continuous( 
            sec.axis = sec_axis(~. *scale, 
                    name= name_second_axis))
        
        result <- try(plot_trait[[TRAIT]][[sp]])
        if (inherits(result, 'try-error')) {
          NAN_remove <- sapply(unique(data_tot_summary_tr_sp$TRAIT), FUN = function(x) all(is.nan(as.matrix(data_tot_summary_tr_sp[which(data_tot_summary_tr_sp$TRAIT == x),"mean"]))))
          if(any(NAN_remove)) {
            data_tot_summary_tr_sp <- data_tot_summary_tr_sp[-which(data_tot_summary_tr_sp$TRAIT == unique(data_tot_summary_tr_sp$TRAIT)[which(NAN_remove)]),]
            title <- paste0(title, "\n removed ", unique(data_tot_summary_tr_sp$TRAIT)[which(NAN_remove)], " due to NAN")
          }
          plot_trait[[TRAIT]][[sp]] <- ggplot(data = data_tot_summary_tr_sp, aes(x = .data[[age_or_cohort]], 
                                                                                 group=interaction(Fish_press, TRAIT))) +
            geom_line(aes(y = mean, col=Fish_press, linetype=TRAIT), lwd=1, alpha=0.6) +
            scale_color_brewer(palette = "Spectral", direction = -1) +
            theme_bw()+
            ylab(label = paste(var_scale_regular, collapse="\n"))+
            xlab(label = x_axis)+
            labs(title = title)+
            # facet_wrap(~ Fish_press)+ 
            scale_linetype_manual(values = data_tot_summary_tr_sp$vec_linetype, 
                                  name="Functional trait")+ 
            theme(legend.key.width = unit(2,"cm"))+
            scale_y_continuous( 
              sec.axis = sec_axis(~. *scale, 
                                  name= name_second_axis))
          
        }# try error loop
      }
    } #trait plot
  } #sp plot
  for (Tr in names(plot_trait)){
    trait_combined <- ifelse(length(list_traits)>1, paste0(Tr, "_WITH_"), Tr)
    suffix_relative <- ifelse(relative_to_ref, "RELATIVE_", "")
    if(!dir.exists(pathout)) dir.create(pathout, recursive = T)
    pdf(file = paste0(pathout, "/", suffix_relative, "plot_", trait_combined, ".pdf"), width = 10, height = 8)
    par(mfrow=c(1,1))
    for (sp_plot in names(plot_trait[[Tr]])){
       scale = scale_v[[Tr]][[sp_plot]]
       print(plot_trait[[Tr]][[sp_plot]])
    }
    dev.off()
  }
}#function



function_1_trait_combined_sp<- function(SP, age_or_cohort, pathout, Fmax_constrain = NULL, Fstart_constrain = NULL, 
                                            list_traits = NULL, data_complete, relative_to_ref = FALSE){
  if (age_or_cohort == "age" & !grepl(pattern="age", pathout)) {
    cat("dataset has no age data in it")
    stop()
  }
  if (age_or_cohort == "cohort" & grepl(pattern="age", pathout)) {
    cat("dataset has age data in it and shouldn't have")
    stop()
  }
  
  plot_trait = list()
  scale_v = list()
  
  # clean dataset
  if (age_or_cohort == "age" & all(!grepl(pattern="age", colnames(data_complete)))) {
    colnames(data_complete)[which(colnames(data_complete) == "time")] <- "age"
  }
  if (any(unique(data_complete$TRAIT) %in% c("mean_of_Z_per_age", "mean_of_Dis_per_age"))) data_complete <- data_complete[-which((data_complete$TRAIT) %in% c("mean_of_Z_per_age", "mean_of_Dis_per_age")),]
  data_complete[,age_or_cohort] <- as.numeric(as.matrix(data_complete[, age_or_cohort]))
  if (age_or_cohort == "age")    data_complete <- subset(data_complete, age <= 30 & age >=1)
  if (age_or_cohort == "cohort") data_complete <- subset(data_complete, cohort >= (min(data_complete[, age_or_cohort])+19) & cohort <= (max(data_complete[,age_or_cohort])-20))
  # add fecunditybyage dimension
  if (!(any(grepl("fec",  data_complete$TRAIT)))) {
    egg = NULL
    abund = NULL
    egg = pivot_wider(data_complete[grep(data_complete$TRAIT, pattern="egg"),], values_from = "value", names_from = "TRAIT")
    abund = pivot_wider(data_complete[c(grep(data_complete$TRAIT, pattern="abundance")),], values_from = "value", names_from = "TRAIT")
    fe_data = full_join(egg, abund)
    fe_data[, "fecundityByAge"] = fe_data[,c(grep(colnames(fe_data), pattern="egg"))]/fe_data[,c(grep(colnames(fe_data), pattern="abundance"))]
    fe_data <- pivot_longer(fe_data[, -c(grep(colnames(fe_data), pattern="egg"), grep(colnames(fe_data), pattern="abund"))], names_to = "TRAIT", values_to = "value", cols = "fecundityByAge")
    data_complete <- rbind(data_complete, fe_data)
  }
  # select the traits to plot 
  TRAITS_to_plot = unique(data_complete$TRAIT) # plot all the traits of the dataset
  other_TRAITS_to_plot = NULL
  if (!is.null(list_traits)){ # or only the ones selected
    TRAITS_to_plot = list_traits[[1]]
    TRAITS_to_plot2 = TRAITS_to_plot
    if (any(grepl(TRAITS_to_plot, pattern="log"))) TRAITS_to_plot2 <- stringr::str_remove(pattern = "log", TRAITS_to_plot) 
    if (any(grepl(TRAITS_to_plot, pattern="exp")) & !any(grepl(TRAITS_to_plot, pattern="exp-"))) TRAITS_to_plot2 <- stringr::str_remove(pattern = "exp", TRAITS_to_plot) 
    if (any(grepl(TRAITS_to_plot, pattern="exp-"))) TRAITS_to_plot2 <- stringr::str_remove(pattern = "exp-", TRAITS_to_plot) 
    if (!any(grepl(pattern=TRAITS_to_plot2, unique(data_complete$TRAIT))) & (TRAITS_to_plot != "fecundityByAge")){
      cat("traits to plot are not in the dataset. Dataset contains : \n", unique(data_complete$TRAIT))
      stop()
    }
    if (length(list_traits)>1) {            # one trait at a time, or multiple
      other_TRAITS_to_plot = list_traits[-1]
    }
  }
  
  for (sp in SP){
    
    # PLOT the mean and the std deviation 
    k=0
    other_TRAIT = NULL
    for (TRAIT in TRAITS_to_plot){
      k=k+1
      for (sp_plotted in sp){
        traits_df <- c(TRAIT, other_TRAITS_to_plot[[k]])
        rm("id_log", "id_exp_neg", "id_exp_pos")
        if (any(grepl(c(TRAIT, other_TRAITS_to_plot[[k]]), pattern="log"))) {
          traits_df <- stringr::str_remove(pattern = "log", c(TRAIT, other_TRAITS_to_plot[[k]])) 
          id_log = grep(c(TRAIT, other_TRAITS_to_plot[[k]]), pattern="log")
        }
        if (any(grepl(c(TRAIT, other_TRAITS_to_plot[[k]]), pattern="exp-"))) {
          traits_df <- stringr::str_remove(pattern = "exp-", traits_df)
          id_exp_neg = which(grepl(c(TRAIT, other_TRAITS_to_plot[[k]]), pattern = "exp-"))
        }
        if (any(grepl(c(TRAIT, other_TRAITS_to_plot[[k]]), pattern="exp")) & !exists("id_exp_neg")) {
          traits_df <- stringr::str_remove(pattern = "exp", traits_df)
          id_exp_pos = which(grepl(c(TRAIT, other_TRAITS_to_plot[[k]]), pattern = "exp") &
                               !grepl(c(TRAIT, other_TRAITS_to_plot[[k]]), pattern = "exp-"))
        }
        data_tot_summary_tr    <- data_complete[which(data_complete$TRAIT %in% traits_df),]
        if (exists("id_log")){
          id_log <- which(data_complete$TRAIT %in% traits_df[id_log])
          log_val  <- log(as.matrix(data_tot_summary_tr[id_log, "value"]))
          log_val[which(is.infinite(log_val))] <- 0
          df_log <- data_tot_summary_tr[id_log,]
          data_tot_summary_tr <- data_tot_summary_tr[-id_log, ]
          df_log[, "value"] <- log_val
          data_tot_summary_tr <- rbind(data_tot_summary_tr, df_log)
        }
        if (exists("id_exp_neg")){
          id_exp_neg <- which(data_complete$TRAIT %in% traits_df[id_exp_neg])
          exp_neg_val  <- exp(-as.matrix(data_tot_summary_tr[id_exp_neg, "value"]))
          df_exp_neg <- data_tot_summary_tr[id_exp_neg,]
          data_tot_summary_tr <- data_tot_summary_tr[-id_exp_neg, ]
          df_exp_neg[, "value"] <- exp_neg_val
          data_tot_summary_tr <- rbind(data_tot_summary_tr, df_exp_neg)
        }
        if (exists("id_exp_pos")){
          id_exp_pos <- which(data_complete$TRAIT %in% traits_df[id_exp_pos])
          exp_pos_val  <- exp(as.matrix(data_tot_summary_tr[id_exp_pos, "value"]))
          df_exp_pos <- data_tot_summary_tr[id_exp_pos,]
          data_tot_summary_tr <- data_tot_summary_tr[-id_exp_pos, ]
          df_exp_pos[, "value"] <- exp_pos_val
          data_tot_summary_tr <- rbind(data_tot_summary_tr, df_exp_pos)
        }
        data_tot_summary_tr_sp <- data_tot_summary_tr[which(data_tot_summary_tr$Species==sp_plotted),] 
        data_tot_summary_tr_sp$value <- as.numeric(data_tot_summary_tr_sp$value)
        data_tot_summary_tr_sp <- data_tot_summary_tr_sp %>%
          group_by(get(age_or_cohort), Species, Fish_press, TRAIT) %>%
          summarise(mean = mean(value, na.rm = TRUE),
                    sd   = sd(value, na.rm = TRUE))
        colnames(data_tot_summary_tr_sp)[which(colnames(data_tot_summary_tr_sp) == "get(age_or_cohort)")] <- age_or_cohort
        # reorder Fmax/Fstart factor for prettier plots
        data_tot_summary_tr_sp$Fish_press    <- factor(as.factor(data_tot_summary_tr_sp$Fish_press), levels = ideal_order_factor)
        
        # select the Fishing pressures I want to keep
        if (!is.null(Fmax_constrain) | !is.null(Fstart_constrain)){
          F_keep = NULL
          F_keep = NULL
          if (!is.null(Fmax_constrain))    F_keep = (grep(data_tot_summary_tr_sp$Fish_press, pattern=paste0("Fmax=", Fmax_constrain, " Fstart")))
          if (!is.null(Fstart_constrain))  F_keep = c(F_keep, grep(data_tot_summary_tr_sp$Fish_press, pattern=paste0("Fstart=", Fstart_constrain, "$")))
          if(any(duplicated(F_keep))) F_keep <- F_keep[duplicated(F_keep)]
          data_tot_summary_tr_sp <- data_tot_summary_tr_sp[F_keep, ]
        }
        reference <- subset((data_tot_summary_tr_sp), Fish_press == paste0("Fmax=", 1, " Fstart=", 1))
        reference <- reference[, -which(colnames(reference) == 'Fish_press')]
        if (relative_to_ref && nrow(reference) > 0) {
          reference <- reference %>%
            dplyr::select(all_of(age_or_cohort), Species, TRAIT, mean) %>%
            dplyr::rename(ref_mean = mean)
          data_tot_summary_tr_sp <- dplyr::left_join(
            data_tot_summary_tr_sp,
            reference,
            by = c(age_or_cohort, "Species", "TRAIT")
          )
          data_tot_summary_tr_sp$mean <- 100*(data_tot_summary_tr_sp$mean - data_tot_summary_tr_sp$ref_mean)/data_tot_summary_tr_sp$ref_mean
          data_tot_summary_tr_sp$ref_mean <- NULL
        }
        
        if(nrow(data_tot_summary_tr_sp)==0){
          plot_trait[[TRAIT]][[sp_plotted]] <- ggplot(data = data_tot_summary_tr_sp, aes(x = age_or_cohort, group = Fish_press)) 
          next
        }
        
        x_axis <- age_or_cohort
        title = paste0("By cohort ", sp_plotted, ifelse(relative_to_ref==T, " % relative (to 1,1) changes of ", " ranges of "), 
                       paste(c(TRAIT, unlist(other_TRAITS_to_plot)), collapse=" "), "\nF varies for ", sp)
        #vec_linetype <- c("solid", "dotted", "dashed", "D3")[1:length(unique(data_tot_summary_tr_sp$TRAIT))]
        # vec_pointtype <- c("solid", "dotted", "dashed", "D3")[1:length(unique(data_tot_summary_tr_sp$TRAIT))]
        data_tot_summary_tr_sp$vec_linetype = sapply(data_tot_summary_tr_sp$TRAIT, function(x) vec_linetype[which(unique(data_tot_summary_tr_sp$TRAIT) == x)])
        data_tot_summary_tr_sp$mean  = as.numeric(data_tot_summary_tr_sp$mean)
        data_tot_summary_tr_sp$sd    = as.numeric(data_tot_summary_tr_sp$sd)
        data_tot_summary_tr_sp$TRAIT = as.factor(data_tot_summary_tr_sp$TRAIT)
        data_tot_summary_tr_sp$Fish_press = as.factor(data_tot_summary_tr_sp$Fish_press)
        NAN_remove <- sapply(unique(data_tot_summary_tr_sp$TRAIT), FUN = function(x) all(is.nan(as.matrix(data_tot_summary_tr_sp[which(data_tot_summary_tr_sp$TRAIT == x),"mean"]))))
        NAN_remove <- sapply(unique(data_tot_summary_tr_sp$TRAIT), FUN = function(x) all(is.nan(as.matrix(data_tot_summary_tr_sp[which(data_tot_summary_tr_sp$TRAIT == x),"mean"]))))
        if(any(NAN_remove)) {
          data_tot_summary_tr_sp <- data_tot_summary_tr_sp[-which(data_tot_summary_tr_sp$TRAIT == unique(data_tot_summary_tr_sp$TRAIT)[which(NAN_remove)]),]
          title <- paste0(title, "\n removed ", unique(data_tot_summary_tr_sp$TRAIT)[which(NAN_remove)], " due to NAN")
        }
        
        
        scale = 1
        var_scale_regular = TRAIT
        name_second_axis = ""
        if(!is.null(other_TRAITS_to_plot)){
          data_key <- data_tot_summary_tr_sp[which(data_tot_summary_tr_sp$TRAIT %in% as.matrix(other_TRAITS_to_plot[[1]])),]
          key <- key_auto()(
            scale_y_continuous(
              limits = range(data_key$mean, na.rm=T),
              breaks = c(quantile(data_key$mean, probs = 0.25, na.rm = T), quantile(data_key$mean, probs = 0.5, na.rm = T), 
                         quantile(data_key$mean, probs = 0.75, na.rm = T)),
              labels = \(x) x / 1000
            )
          )
          for_scale <- data_tot_summary_tr_sp[!is.infinite(data_tot_summary_tr_sp$mean), ]
          mean_var  <- aggregate(for_scale[, "mean"], list(for_scale$TRAIT), FUN = function(x) mean(x, na.rm=T))
          median_var<- aggregate(for_scale[, "mean"], list(for_scale$TRAIT), median)
          var_scale_regular = mean_var[which(mean_var[,2] > 0.8*(max(mean_var[,2]))),1]
          scale = mean(mean_var[-which.max(mean_var[,2]),"mean" ], na.rm=T)/
            quantile(mean(mean_var[which.max(mean_var[,2]),"mean" ], na.rm=T), 0.95)
          data_tot_summary_tr_sp[-which(data_tot_summary_tr_sp$TRAIT %in%  var_scale_regular), "mean"] <- 
            data_tot_summary_tr_sp[-which(data_tot_summary_tr_sp$TRAIT %in% var_scale_regular), "mean"]/scale
          name_second_axis = paste(mean_var[-which(mean_var[,2] > 0.8*(max(mean_var[,2]))),1], collapse = "\n")
        }
        scale <- ifelse(scale ==0, 1, scale)
        scale_v[[TRAIT]][[sp]] <- scale
        
        if (length(unique(data_tot_summary_tr_sp$TRAIT))>3) {
          col_fac = "TRAIT"
          shape_fac = "Fish_press"
        } else {
          col_fac = "Fish_press"
          shape_fac = "TRAIT"
        }
        
        
        plot_trait[[TRAIT]][[sp]] <- ggplot(data = data_tot_summary_tr_sp, aes(x = .data[[age_or_cohort]], 
                                                                               group=interaction(Fish_press, TRAIT))) +
          # geom_ribbon(aes(ymin = mean-sd, ymax = mean + sd),
          # alpha = 0.5, fill = "grey70")+
          geom_line(aes(y = mean, col=.data[[col_fac]], linetype=.data[[shape_fac]]), lwd=0.5, alpha=0.6) +
          # geom_point(aes(y = mean, col=.data[[col_fac]], pch=.data[[shape_fac]]), alpha=0.5, col="black") + #linetype=data_tot_summary_tr_sp$TRAIT
          scale_color_brewer(palette = "Spectral", direction = -1) +
          theme_bw()+
          ylab(label = paste(var_scale_regular, collapse="\n"))+
          xlab(label = x_axis)+
          labs(title = title)+
          # facet_wrap(~ Fish_press)+ 
          # scale_linetype_manual(#values = data_tot_summary_tr_sp$vec_linetype, 
          #   name="Functional trait")+ 
          # scale_shape_manual(#(values = data_tot_summary_tr_sp$vec_pointtype,
          #   name="Functional trait")+ 
          theme(legend.key.width = unit(2,"cm"))+
          scale_y_continuous( 
            sec.axis = sec_axis(~. *scale, 
                                name= name_second_axis))
        
        result <- try(plot_trait[[TRAIT]][[sp]])
        if (inherits(result, 'try-error')) {
          NAN_remove <- sapply(unique(data_tot_summary_tr_sp$TRAIT), FUN = function(x) all(is.nan(as.matrix(data_tot_summary_tr_sp[which(data_tot_summary_tr_sp$TRAIT == x),"mean"]))))
          if(any(NAN_remove)) {
            data_tot_summary_tr_sp <- data_tot_summary_tr_sp[-which(data_tot_summary_tr_sp$TRAIT == unique(data_tot_summary_tr_sp$TRAIT)[which(NAN_remove)]),]
            title <- paste0(title, "\n removed ", unique(data_tot_summary_tr_sp$TRAIT)[which(NAN_remove)], " due to NAN")
          }
          plot_trait[[TRAIT]][[sp]] <- ggplot(data = data_tot_summary_tr_sp, aes(x = .data[[age_or_cohort]], 
                                                                                 group=interaction(Fish_press, TRAIT))) +
            geom_line(aes(y = mean, col=Fish_press, linetype=TRAIT), lwd=1, alpha=0.6) +
            scale_color_brewer(palette = "Spectral", direction = -1) +
            theme_bw()+
            ylab(label = paste(var_scale_regular, collapse="\n"))+
            xlab(label = x_axis)+
            labs(title = title)+
            # facet_wrap(~ Fish_press)+ 
            scale_linetype_manual(values = data_tot_summary_tr_sp$vec_linetype, 
                                  name="Functional trait")+ 
            theme(legend.key.width = unit(2,"cm"))+
            scale_y_continuous( 
              sec.axis = sec_axis(~. *scale, 
                                  name= name_second_axis))
          
        }# try error loop
      }
    } #trait plot
  } #sp plot
  for (Tr in names(plot_trait)){
    trait_combined <- ifelse(length(list_traits)>1, paste0(Tr, "_WITH_"), Tr)
    suffix_relative <- ifelse(relative_to_ref, "RELATIVE_", "")
    if(!dir.exists(pathout)) dir.create(pathout, recursive = T)
    pdf(file = paste0(pathout, "/", suffix_relative, "plot_", trait_combined, ".pdf"), width = 14, height = 10)
    par(mfrow=c(1,1))
    # for (sp_plot in names(plot_trait[[Tr]])){
    #   scale = scale_v[[Tr]][[sp_plot]]
    combined <- ggarrange(plotlist=(plot_trait[[Tr]]), ncol=4, nrow=2, common.legend = T)
    plot(combined[[1]])
    plot(combined[[2]])
    # }
    dev.off()
  }
}#function











