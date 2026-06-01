
main = "New_gpu_conf.R"
source("osmose-ns-v4.4.1_clim/R/functions/source.R")

base_dir <- "osmose-ns-v4.4.1_clim/output/PLOTS/All_clean_genetic_drift_indicators/INIT/"
folders  <- list.dirs(base_dir, recursive = F)
INIT2 <- get(load(list.files(folders[1], full.names = T)))

folders_traits_y <- "osmose-ns-v4.4.1_clim/output/PLOTS/Trait_variables_by_year/RData/"
folders_traits_y <- list.dirs(folders_traits_y, recursive = F)

folders_traits_c <- "osmose-ns-v4.4.1_clim/output/PLOTS/Trait_variables_by_cohort/RData/"
folders_traits_c  <- list.dirs(folders_traits_c, recursive = F)

# Auxiliar functions ------------------------------------------------------

# load
source("osmose-ns-v4.4.1_clim/R/functions/auxiliary_for_lists.R")
source("osmose-ns-v4.4.1_clim/R/functions/functions_output_SimuRData_into_df_for_plots.R")


# Processing one list -----------------------------------------------------

# INIT2b = INIT2
# INIT2b[2:3] = lapply(INIT2b[2:3], FUN = myf2)
# INIT2b[1:3] = lapply(INIT2b[1:3], FUN = reshape_list)
# INIT2b[[5]] = lapply(INIT2b[[5]], FUN = as.list)
# INIT2b[[5]] = reshape_list(INIT2b[[5]])
# INIT2b[[5]] = lapply(INIT2b[[5]], FUN = myf)
# INIT2b[[5]]$age = NULL
# INIT2b[[1]] = NULL
# INIT2b = reshape_list(INIT2b)


INIT2b = INIT2
INIT2b = lapply(INIT2b, FUN = myf2)
INIT2b = reshape_list(INIT2b)



# longeH = reshape_data(INIT2b$Haddock$data_Ne$generation_length, "generation_length", as.variable = TRUE)
# Ne_N1H = reshape_data(INIT2b$Haddock$data_Ne$Ne_N_1, "Ne_N1", as.variable = TRUE)
# dat = rbind(longeH, Ne_N1H)
# dat$fishing = 0


# Loop for Fishing Pressure -----------------------------------------------


pathout <- paste0(configDir, "/output/PLOTS/Htz_vs_Ne_variables")

# from RData output of the loops producing the plots of Trait by time and Genetic drift indicators, to df
# produce a dataframe with the right format : Ne variables and htz variables together
Time1=150
Time2=160
species_df=species
list_data_traits <- data_trait_and_genetic_into_df(Time1=150, Time2=160, pathoutput=paste0(configDir, "/output/PLOTS/RDATA"), 
                                                   FORCE = F, species = species_df)

# produce a dataframe with the right format : Ne variables and htz variables together
list_data_traits_ratio <- data_trait_and_genetic_into_df(Time1=150, Time2=160, pathoutput=paste0(configDir, "/output/PLOTS/RDATA"), 
                                                              FORCE = T, species = species_df, RATIO=T)


# list of outputs 
list_outputs <- list("traits"=list_data_traits, "trait_ratio"=list_data_traits_ratio)




# # Bind together the trait data coming from similar dimensions : cohort or year data together (by genetic drift dimension)

pathoutputRDATA = paste0(configDir, "/output/PLOTS/RDATA")
data_full_merged <- dataset_attach_all_traits_genetic_drift(pathoutputRDATA, list_outputs=list("traits"=list_data_traits), 
                                                            trait_cv = list("folders_traits_c"=folders_traits_c, "folders_traits_y"=folders_traits_y))
dataRATIO_full_merged <- dataset_attach_all_traits_genetic_drift(pathoutputRDATA, list_outputs=list("ratio"=list_data_traits_ratio), 
                                                                 trait_cv = list("folders_traits_c"=folders_traits_c, "folders_traits_y"=folders_traits_y))


data_full_merged_age <- dataset_attach_all_traits_genetic_drift(pathoutputRDATA, list_outputs = list("traits"=list_data_traits), s_pecies = species,
                                                                trait_cv = list("folders_traits_c_a"=folder_traits_c_a, "folders_traits_y_a"=folder_traits_c_a),
                                                                keep_cohort_or_age_dimension = T, FORCE=F)
# 
# # Dataset ready for modeling --------------------------------------------
# # scaling as a standardized loop
data_mean_sc = scale_dataset_list(pathoutputRDATA, data_mean=data_full_merged[[1]], FORCE=FALSE)
dataRATIO_mean_sc = scale_dataset_list(pathoutputRDATA, data_mean=dataRATIO_full_merged[[1]], FORCE=FALSE)
# 
# # new dataset of the difference between Fishing pressure 1,1 and Else--------------------
diff_data = relative_dataset(data_mean_sc, pathoutputRDATA, FORCE=F)
diffRATIO_data = relative_dataset(dataRATIO_mean_sc, pathoutputRDATA, FORCE=FALSE, RATIO=T)


# plot TRAITs against GENETIC drift var

p =NULL
for (SP in species){
  for (species_plotted in SP){
    cat(SP)
    for (i in names(data_mean_sc)){
      cat(i)
      dat_plot_w <- data_mean_sc[[i]][which(data_mean_sc[[i]]$Species == SP), ]
      i_i <- colnames(dat_plot_w)[grep(colnames(dat_plot_w), pattern = i)]
      dat_plot_w$Fmax = as.factor(dat_plot_w$Fmax)
      dat_plot_w$Fstart = as.factor(dat_plot_w$Fstart)

      for (trait_var in colnames(dat_plot_w)[-c(which(colnames(dat_plot_w) %in% c("replicate", "Fstart", "Fmax", "Species", "F_p")),
                                                grepl(i, pattern="htz_slope"),grepl(i, pattern="DisbyY"),grepl(i, pattern="ZbyY"))]){
        cat(trait_var)
        # removing when Ne =0 and trait =0
        # if (any(dat_plot_w[[ncol(dat_plot_w)]] == 0 & dat_plot_w[[ncol(dat_plot_w)-1]] == 0)) {
          # dat_plot_w <- dat_plot_w[-which(dat_plot_w[[ncol(dat_plot_w)]] == 0 & dat_plot_w[[ncol(dat_plot_w)-1]] == 0),] }

        if(range(dat_plot_w[, trait_var], na.rm = T)[1]>1){
          p[[i]][[trait_var]]<-
            ggplot(dat_plot_w, aes_string(x=trait_var,  y=i_i, color="Fstart"))+
            geom_point(alpha=0.6)+
            geom_smooth(method = "lm", formula = y ~ x, se = TRUE, alpha=0.4, aes(
              group = Fstart), col="black")+
            facet_grid(~Fmax, labeller = labeller(.rows = label_both, .cols = label_both))+
            labs(title = paste0(species_plotted, " ", i, " against ", trait_var, "\nReference species is ", SP, "+ lm"))+
            ylab(label = i)+
            #geom_abline(intercept = 0, slope = 1, col="red")+
            scale_x_log10() +
            scale_y_log10()+
            theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
        } else {
          p[[i]][[trait_var]] <-
            ggplot(dat_plot_w, aes_string(x=trait_var,  y=i_i, color="Fstart"))+
            geom_point(alpha=0.6)+
            geom_smooth(method = "lm", formula = y ~ x, se = TRUE, alpha=0.4, aes(
              group = dat_plot_w$Fstart), col="black")+
            facet_grid(~Fmax, labeller = labeller(.rows = label_both, .cols = label_both))+
            labs(title = paste0(species_plotted, " ", i, " against ", trait_var, "\nReference species is ", SP, "+ lm"))+
            ylab(label = i)+
            theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

        }
      }
      if(!is.null(p[[i]])){
        cat("Saved for species ", SP, "\n")
        pathout = paste0("osmose-ns-v4.4.1_clim/output/PLOTS/Var_VS_Var/",  Time1, "_", Time2 ,"/NEW",i)
        if(!dir.exists(pathout)) dir.create(pathout, recursive = T)
        pdf(file = paste0(pathout, "/plot_", i, "_Fvarforspecies_", SP, ".pdf"), width = 10, height = 8)
        par(mfrow=c(1,1))
        print(p[[i]])
        dev.off()

      } # save plots
    }
  }# sp plotted
}




## SAME BUT WITH THE RATIOS ---------------------



# plot TRAITs against GENETIC drift var

p =NULL
for (SP in species[8:length(species)]){
  for (species_plotted in SP){
    cat(SP)
    for (i in names(dataRATIO_mean_sc)){  
      cat(i)
      dat_plot_w <- dataRATIO_mean_sc[[i]][which(dataRATIO_mean_sc[[i]]$Species == SP), ]
      i_i <- colnames(dat_plot_w)[grep(colnames(dat_plot_w), pattern = i)]
      dat_plot_w$Fmax = as.factor(dat_plot_w$Fmax)
      dat_plot_w$Fstart = as.factor(dat_plot_w$Fstart)
      
      for (trait_var in colnames(dat_plot_w)[-c(which(colnames(dat_plot_w) %in% c("replicate", "Fstart", "Fmax", "Species", "F_p")), 
                                                grepl(i, pattern="htz_slope"),grepl(i, pattern="DisbyY"),grepl(i, pattern="ZbyY"))]){
        cat(trait_var)
        # removing when Ne =0 and trait =0
        # if (any(dat_plot_w[[ncol(dat_plot_w)]] == 0 & dat_plot_w[[ncol(dat_plot_w)-1]] == 0)) {
        # dat_plot_w <- dat_plot_w[-which(dat_plot_w[[ncol(dat_plot_w)]] == 0 & dat_plot_w[[ncol(dat_plot_w)-1]] == 0),] }
        
        if(range(dat_plot_w[, trait_var], na.rm = T)[1]>0 & range(dat_plot_w[, i_i], na.rm = T)[1]>0){
          p[[i]][[trait_var]] <-
            ggplot(dat_plot_w, aes_string(x=trait_var,  y=i_i, color="Fstart"))+
            geom_point(alpha=0.6)+
            geom_smooth(method = "lm", formula = y ~ x, se = TRUE, alpha=0.4, aes(
              group = Fstart), col="black")+            
            facet_grid(~Fmax, labeller = labeller(.rows = label_both, .cols = label_both))+
            labs(title = paste0(species_plotted, " ", i, " against ", trait_var, "\nReference species is ", SP, "+ lm"))+
            ylab(label = i)+
            #geom_abline(intercept = 0, slope = 1, col="red")+
            scale_x_log10() +
            scale_y_log10()+
            theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
        } else {
          p[[i]][[trait_var]] <-
            ggplot(dat_plot_w, aes_string(x=trait_var,  y=i_i, color="Fstart"))+
            geom_point(alpha=0.6)+
            geom_smooth(method = "lm", formula = y ~ x, se = TRUE, alpha=0.4, aes(
              group = dat_plot_w$Fstart), col="black")+
            facet_grid(~Fmax, labeller = labeller(.rows = label_both, .cols = label_both))+
            labs(title = paste0(species_plotted, " ", i, " against ", trait_var, "\nReference species is ", SP, "+ lm"))+
            ylab(label = i)+
            theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
          
        }
      }
      if(!is.null(p[[i]])){
        cat("Saved for species ", SP, "\n")   
        pathout = paste0("osmose-ns-v4.4.1_clim/output/PLOTS/Var_VS_Var/",  Time1, "_", Time2 ,"/NEW",i)
        if(!dir.exists(pathout)) dir.create(pathout, recursive = T)
        pdf(file = paste0(pathout, "/plotRATIO_", i, "_Fvarforspecies_", SP, ".pdf"), width = 10, height = 8)
        par(mfrow=c(1,1))
        print(p[[i]])
        dev.off()  
        
      } # save plots
    } 
  }# sp plotted
}









## SAME BUT WITH THE DIFFERENCE (relative) ---------------------

# plot TRAITs against GENETIC drift var

p =NULL
for (SP in species){
  for (species_plotted in SP){
    cat(SP)
    for (i in names(diff_data)){  
      cat(i)
      dat_plot_w <- diff_data[[i]][which(diff_data[[i]]$Species == SP), ]
      dat_plot_w <- pivot_wider(data=dat_plot_w, names_from = "variable", values_from = "value")
      i_i <- colnames(dat_plot_w)[grep(colnames(dat_plot_w), pattern = i)]
      dat_plot_w$Fmax = as.factor(dat_plot_w$Fmax)
      dat_plot_w$Fstart = as.factor(dat_plot_w$Fstart)
      
      for (trait_var in colnames(dat_plot_w)[-c(which(colnames(dat_plot_w) %in% c("replicate", "Fstart", "Fmax", "Species", "F_p")), 
                                                grepl(i, pattern="htz_slope"),grepl(i, pattern="DisbyY"),grepl(i, pattern="ZbyY"))]){
        cat(trait_var)
        # removing when Ne =0 and trait =0
        # if (any(dat_plot_w[[ncol(dat_plot_w)]] == 0 & dat_plot_w[[ncol(dat_plot_w)-1]] == 0)) {
        # dat_plot_w <- dat_plot_w[-which(dat_plot_w[[ncol(dat_plot_w)]] == 0 & dat_plot_w[[ncol(dat_plot_w)-1]] == 0),] }
        
        if(range(dat_plot_w[, trait_var], na.rm = T)[1]>0 & range(dat_plot_w[, i_i], na.rm = T)[1]>0){
          p[[i]][[trait_var]] <-
            ggplot(dat_plot_w, aes_string(x=trait_var,  y=i_i, color="Fstart"))+
            geom_point(alpha=0.6)+
            geom_smooth(method = "lm", formula = y ~ x, se = TRUE, alpha=0.4, aes(
              group = Fstart), col="black")+            
            facet_grid(~Fmax, labeller = labeller(.rows = label_both, .cols = label_both))+
            labs(title = paste0(species_plotted, " ", i, " against ", trait_var, "\nReference species is ", SP, "+ lm"))+
            ylab(label = i)+
            #geom_abline(intercept = 0, slope = 1, col="red")+
            scale_x_log10() +
            scale_y_log10()+
            theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
        } else {
          p[[i]][[trait_var]] <-
            ggplot(dat_plot_w, aes_string(x=trait_var,  y=i_i, color="Fstart"))+
            geom_point(alpha=0.6)+
            geom_smooth(method = "lm", formula = y ~ x, se = TRUE, alpha=0.4, aes(
              group = dat_plot_w$Fstart), col="black")+
            facet_grid(~Fmax, labeller = labeller(.rows = label_both, .cols = label_both))+
            labs(title = paste0(species_plotted, " ", i, " against ", trait_var, "\nReference species is ", SP, "+ lm"))+
            ylab(label = i)+
            theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
          
        }
      }
      if(!is.null(p[[i]])){
        cat("Saved for species ", SP, "\n")   
        pathout = paste0("osmose-ns-v4.4.1_clim/output/PLOTS/Var_VS_Var/",  Time1, "_", Time2 ,"/NEW",i)
        if(!dir.exists(pathout)) dir.create(pathout, recursive = T)
        pdf(file = paste0(pathout, "/plotDIFF_", i, "_Fvarforspecies_", SP, ".pdf"), width = 10, height = 8)
        par(mfrow=c(1,1))
        print(p[[i]])
        dev.off()  
        
      } # save plots
    } 
  }# sp plotted
}










# folders
# for (SP in species){
#   dat_tot = NULL
#   dat_plot = NULL
#   p = list()
#   ID_fold <- grep(SP, folders)
  
#   trait_file <- list.files(folders_traits[grep(folders_traits, pattern=SP)], full.names = T)
#   trait_var <- get(load(trait_file))
  
#   mortalityH <- trait_var[which(trait_var$TRAIT == "mortality"), ]
#   mortalityH <- mortalityH[which(mortalityH$Species == SP), ]
#   mortalityH$Fmax <- unlist(stringr::str_split(unlist(stringr::str_split(mortalityH$Fish_press, "Fmax="))[seq(by=2, 2, 2*length(mortalityH$Fish_press))], " "))[seq(by=2, 1, 2*length(mortalityH$Fish_press))]
#   mortalityH$Fstart <- unlist(stringr::str_split(mortalityH$Fish_press, "Fstart="))[seq(by=2, 2, 2*length(mortalityH$Fish_press))]
#   mortalityH <- mortalityH[, which(colnames(mortalityH) %in% c("time", "mean", "Fmax", "Fstart"))]
#   colnames(mortalityH) <- c("cohort", "value", "Fmax", "Fstart")
#   mortalityH$variable = "mean_natural_mortality"
#   mortalityH_l <-bind_rows(replicate(20,mortalityH, simplify=FALSE)) 
#   mortalityH_l$replicate <- as.vector(sapply(1:20, FUN=function(x) rep(x, times=nrow(mortalityH))))

  
#   for (f_sp in ID_fold){
#     files <- list.files(folders[f_sp], full.names = T)
#     data_sp <- get(load(files))
#     Fmax <- unlist(stringr::str_split(unlist(stringr::str_split(folders[f_sp], "Fmax"))[2], "_"))[1]
#     Fstart <- unlist(stringr::str_split(unlist(stringr::str_split(folders[f_sp], "Fstart"))[2], "_"))[1]
  
    
#     INIT2b = data_sp
#     INIT2b[2:3] = lapply(INIT2b[2:3], FUN = myf2)
#     INIT2b[1:3] = lapply(INIT2b[1:3], FUN = reshape_list)
#     INIT2b[[5]] = lapply(INIT2b[[5]], FUN = as.list)
#     INIT2b[[5]] = reshape_list(INIT2b[[5]])
#     INIT2b[[5]] = lapply(INIT2b[[5]], FUN = myf)
#     INIT2b[[5]]$age = NULL
#     INIT2b[[1]] = NULL
#     INIT2b = reshape_list(INIT2b)

    
#     longeH = reshape_data(INIT2b[[SP]]$data_Ne$longevity, "longevity", as.variable = TRUE)
#     Ne_N1H = reshape_data(INIT2b[[SP]]$data_Ne$Ne_N_1, "Ne_N1", as.variable = TRUE)
#     Ne_NalphaH = reshape_data(INIT2b[[SP]]$data_Ne$Ne_N_alpha, "Ne_Nalpha", as.variable = TRUE)
#     Vk_pointH = reshape_data(INIT2b[[SP]]$data_Ne$Vk_point, "Vk_point", as.variable = TRUE)
  
    
#     Nb_H_l = expand.grid(cohort=(1:nrow(INIT2b[[SP]]$data_Nb$Nb))+80,replicate=as.numeric(gsub(colnames(INIT2b[[SP]]$data_Nb$Nb), pattern="rep_", replacement="")))
#     Nb_H_l$value = as.vector(INIT2b[[SP]]$data_Nb$Nb)
#     Nb_H_l$variable = "Nb"
#     Ne_t = myf(INIT2b[[SP]]$tempoNe)
#     Ne_t_l = expand.grid(cohort=(1:nrow(Ne_t))+80,replicate=as.numeric(gsub(colnames(Ne_t), pattern="rep_", replacement="")))
#     Ne_t_l$value = as.vector(Ne_t)
#     Ne_t_l$variable = "Ne_t"
#     htzH_l = expand.grid(cohort=(1:nrow(INIT2b[[SP]]$htz_list))+80,replicate=as.numeric(gsub(colnames(INIT2b[[SP]]$htz_list), pattern="rep_", replacement="")))
#     htzH_l$value = as.vector(INIT2b[[SP]]$htz_list)
#     htzH_l$variable = "htz"

    
#     dat = rbind(longeH, Ne_N1H, Ne_NalphaH, htzH_l, Ne_t_l, Nb_H_l, Vk_pointH)
#     dat$Fmax = Fmax
#     dat$Fstart = Fstart
#     # dat$Species = SP
    
#     if(f_sp != ID_fold[1]) dat_tot = rbind(dat_tot, dat)
#     if(f_sp == ID_fold[1]) dat_tot = dat
#   }
  
#   dat_tot <- rbind(dat_tot, mortalityH_l)
  
#   dat_plot <- dat_tot
  
#   Y_TRAIT = "htz"#c("Ne_N1")
  
#   # dat_plot <- dat_tot[which(dat_tot$Species==SP),]
#   dat_plot$y <- c(dat_plot$variable %in% Y_TRAIT)
#   for (trait in names(table(dat_plot$variable[-which(dat_plot$y)]))){
#     dat_plot_t <- dat_plot[which((dat_plot$cohort %in% c((min(dat_plot$cohort)+20):(max(dat_plot$cohort)-20)))),]
#     dat_plot_t <- dat_plot_t[which((dat_plot_t$variable %in% c(trait)) | dat_plot_t$y),]
#     dat_plot_t <- dat_plot_t[, -which(colnames(dat_plot_t) == "y")]
#     dat_plot_w <- dat_plot_t %>%
#       pivot_wider(
#         names_from  = variable,
#         values_from = value
#       )
#     ID_y <- which(colnames(dat_plot_w) == Y_TRAIT)
#     ID_x = ifelse(ID_y==5, c(6),c(5))
#     y_var = colnames(dat_plot_w)[ID_y]
#     if (Y_TRAIT %in% c("htz")) {
#       dat_plot_w[, Y_TRAIT] =  log(dat_plot_w[,Y_TRAIT])
#       y_var <- paste0("log ",y_var)
#     }
#     if (all(is.na(dat_plot_w[, trait]))) next
#     p[[SP]][[trait]] <- ggplot(dat_plot_w, aes_string(x=colnames(dat_plot_w)[ID_x], y=colnames(dat_plot_w)[ID_y]))+
#       geom_point()+
#       geom_smooth(method = loess, formula = y ~ x, se = FALSE)+
#       facet_grid(Fstart~Fmax,labeller = labeller(.rows = label_both, .cols = label_both))+
#       labs(title = paste0(SP, " ", y_var, " against ", colnames(dat_plot_w)[ID_x]))+
#       ylab(label = y_var)+
#       geom_abline(intercept = 0, slope = 1, col="red")
#   }
#   pathout = paste0("osmose-ns-v4.4.1_clim/output/PLOTS/Var_VS_Var/log",Y_TRAIT)
#   if(!dir.exists(pathout)) dir.create(pathout)
#   for (n in names(p)){
#     pdf(file = paste0(pathout, "/plot_", n, ".pdf"), width = 10, height = 8)
#     par(mfrow=c(1,1))
#     print(p[[SP]])
#     dev.off()
#   }
# }



# Plots -------------------------------------------------------------------
