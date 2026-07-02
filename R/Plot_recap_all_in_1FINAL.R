# the idea is to store separately all the key parameters per species into plots
# and to produce one pdf sheet per species with all the key variables to have it all summarized at one glance




#### Initialisation --------------
library(legendry)
library(ggpubr)

# functions and environment
main = "New_gpu_conf.R" # config file
source("R/functions/source.R") # load config and key parameters and functions
source("R/functions/function_data_traits_summary.R") # datasets preparation & loading
source("R/functions/functions_key_traits.R") # plots 

# paths
# path <- paste0(simu_repo, "/output/PLOTS/FINAL_Ne_Htz_mature_short/")
PATHOUT  = "output/PLOTS/KEY_TRAITS_MATURE/"
if (!dir.exists(PATHOUT)) dir.create(PATHOUT, recursive = TRUE)
if (!dir.exists(paste0(PATHOUT, "/Clean_for_paper"))) dir.create(paste0(PATHOUT, "/Clean_for_paper"), recursive = TRUE)
pathoutputRDATA = paste0(PATHOUT, "RDATA_MATURE")
if(!dir.exists(pathoutputRDATA)) dir.create(pathoutputRDATA, recursive = T)
base_path = paste0("output/DATASETS/Chap2/Trait_variables_by_cohort_by_age2")
base_path = paste0("output/DATASETS/Chap2/Trait_variables_by_cohort_by_age_withWbyAge/")
# folders_long <- list.dirs(paste0(path_long, "INIT"), recursive = F)

folders_traits_c <- paste0("output/DATASETS/Chap2/Trait_variables_by_cohort/RData/")
folders_traits_c  <- list.dirs(folders_traits_c, recursive = F)
folders_traits_c_a <- paste0("output/DATASETS/Chap2/Trait_variables_by_cohort_by_age/RData/")
folders_traits_c_a <- paste0("output/DATASETS/Chap2/Trait_variables_by_cohort_by_age_withWbyAge/RData/")
folders_traits_c_a  <- list.dirs(folders_traits_c_a, recursive = F)

#### DATASETS load or format dataset ------------------------
Time1=150
Time2=160
species_df=species
sp_i = "Hake"
FORCE=T
dataset="data_Ne"

## % activate if needed : will produce .RData files in the folder paste0(simu_repo, "/output/PLOTS/FINAL_Ne_Htz_mature/")
## % 1 file per species, with the dataset at the scale of the cohorts, starting at age = alpha (first mature age)
# dataset_Ne_mean_cohorts(species, FORCE, pathdata=paste0("output/DATASETS/FINAL_Ne_Htz_mature/"),
#         pathout=pathoutputRDATA, Time1, Time2)

## loading the dataset of traits averaged over cohorts Time1 to Time2, and reorganised to remove age=0, add reassign natural mortality as add mortality, and estimate natural mortality
data_complete_c <- load_trait_df(base_path = base_path, pathout=pathoutputRDATA, FORCE=F)




# _  ------------------



######## PLOT -0 selectivity curves --------------
path_surv_curve <- paste0("output/PLOTS/selectivity_curves_per_Fmax")
load(file=paste0(path_surv_curve, "dataset_surv-curve.RData"))
surv_curve = list()
for (sp in species){
  data_surv_curve[[sp]]$Fmax <- stringr::str_remove(data_surv_curve[[sp]]$Fmax, pattern = "Fmax=")
  data_surv_curve[[sp]]$Fmax = as.factor(data_surv_curve[[sp]]$Fmax)
  surv_curve[[sp]] <- ggplot(data_surv_curve[[sp]], aes(x=length, y=F_sc, group=Fmax)) + 
    geom_line(aes(linetype = data_surv_curve[[sp]]$Fmax), linewidth=1.3)+
    labs(title="Fisheries selectivity curve", y="", x="Body length", linetype="Fmax")+ #, y="SSD\u0262"
    theme(strip.background = element_rect(fill="white"),
          panel.background = element_rect(fill = "white", colour = "grey50"),
          plot.title = element_text(hjust = 0.5, size=15),
          axis.text.x = element_text(angle = 60, hjust = 1),
          axis.text = element_text(colour = 1, size = 12),
          panel.border = element_rect(colour = "black", fill=NA),
          strip.text = element_text(size = 12), 
          legend.key.size = unit(1.5, "cm"),
          legend.text=element_text(size=12))
}


path_surv_curve <- paste0("output/PLOTS/selectivity_curves_per_Fmax")
load(file=paste0(path_surv_curve, "dataset_surv-curve.RData"))
surv_curve = list()
k=0
for (sp in species){
  k=k+1
  data_surv_curve[[sp]]$Fmax <- stringr::str_remove(data_surv_curve[[sp]]$Fmax, pattern = "Fmax=")
  data_surv_curve[[sp]]$Fmax = as.factor(data_surv_curve[[sp]]$Fmax)
  surv_curve[[k]] <- ggplot(data_surv_curve[[sp]], aes(x=length, y=F_sc, group=Fmax)) + 
    geom_line(aes(linetype = Fmax), linewidth=1.3)+
    labs(title=sp, y="", x="Body length", linetype="Fmax")+ #, y="SSD\u0262"
    theme(strip.background = element_rect(fill="white"),
          panel.background = element_rect(fill = "white", colour = "grey50"),
          plot.title = element_text(hjust = 0.5, size=15),
          axis.text.x = element_text(angle = 60, hjust = 1),
          axis.text = element_text(colour = 1, size = 12),
          panel.border = element_rect(colour = "black", fill=NA),
          strip.text = element_text(size = 12), 
          legend.key.size = unit(1.5, "cm"),
          legend.text=element_text(size=12))
}


windowsFonts(myCustomWindowsFontName=windowsFont("DejaVu Sans"))
cairo_pdf(file = paste0(PATHOUT, "/Clean_for_paper/TOT_surv_curve.pdf"),
          width = 10, height =10, onefile=T)
for (sp in  1:ceiling(length(species)/4)){
  if (sp == 4) {plot_list <- c(surv_curve[[(sp-1)*4+1]], 
                              surv_curve[[(sp-1)*4+2]],
                              surv_curve[[(sp-1)*4+3]])
  }else {
    plot_list <- c(surv_curve[[(sp-1)*4+1]], 
                   surv_curve[[(sp-1)*4+2]],
                   surv_curve[[(sp-1)*4+3]], 
                   surv_curve[[(sp-1)*4+4]])
  }
  combined_plot <- ggarrange(plotlist = plot_list,
                             ncol = 2,
                             nrow = 2, common.legend = T)
  combined_plot <- annotate_figure(combined_plot)
  print(combined_plot)
}
dev.off()  

#### PLOT 0- Biomass vs Fmax & keep it on the side -------------------
biomassplot = NULL
rel_biomassplot = NULL
list_relative_biomass <- list()
abundanceplot = list()
biomassplot_curve = list()
abundanceplot_curve = list()
for (nm in "data_Ne"){
  for (sp in species){
    # d <- data_complete_c
    # d <- d[d$Fstart==1, ]
    # id_biom <- c(grep(pattern='biomass', d$variable), grep(pattern='abundance', d$variable))
    complete <- get(load(file=paste0(pathoutputRDATA, "/Complete_INIT2", sp, ".RData")))
    d<- complete$complete_with_age
    # we want census abundance = NT
    # min_age <- min(as.numeric(complete[[2]]$age), na.rm=T)
    # rm(complete)
    # d <- d[which(d$age == min_age),]    
    id_biom <- c(grep(pattern='biomass', d$variable), grep(pattern='abundance', d$variable))
    d_i_i <- d[id_biom, ]
    d_i_i <- d_i_i[which(d_i_i$Species == sp), ]
    d_i_i <- pivot_wider(d_i_i, values_from = "value", names_from = "variable")
    if (any(is.na(d_i_i$mean_of_biomassByAge_per_age))) d_i_i <- d_i_i[-which(is.na(d_i_i[, "mean_of_biomassByAge_per_age"])), ]

    # biomassplot[[sp]] <- ggplot(d_i_i, aes(y=mean_of_biomassByAge_per_age, x=Fmax)) +
    #   geom_point(aes(colour = Iszero))+
    #   labs(title = "Biomass per fishing pressure", subtitle = sp)
    
    abundanceplot[[sp]] <- ggplot(d_i_i, aes(y=abundanceByAge, x=Fmax)) +
      geom_point(aes(colour = Iszero))+
      labs(title = "Abundance per fishing pressure", subtitle = sp)
    
    
    abundanceplot_curve[[sp]] <- ggplot(d_i_i, aes(y=abundanceByAge, x=Fmax)) +
      geom_smooth(aes(colour = Iszero), se = F)+
      labs(title = "Abundance per fishing pressure", subtitle = sp)+ 
      theme_bw()+
      ylab("Abundance")+
      theme(legend.position = "none")
    
    d_i_i <- d_i_i %>%
      group_by(Fmax, Fstart, replicate, Species)%>%
      summarise(abundanceByAge = sum(abundanceByAge))
    d_i_i$Iszero <- ifelse(d_i_i$abundanceByAge < 100, "Abundance <100", "Abundance >100")
    
    ref_data_percent <- d_i_i %>%
      filter(Fmax == 0.5 & Fstart == 1) %>%
      dplyr::select(replicate, abundanceByAge) %>%
      dplyr::rename(ref_value_percent_abund = abundanceByAge)
    # Subtract reference from all means
    rel_d_i_i_perc <- d_i_i %>%
      left_join(ref_data_percent[, -which(colnames(ref_data_percent) %in% c("Fmax", "Fstart"))], by = c("replicate"), relationship = "many-to-one") %>%
      mutate(relative_abundance_percent = 100*abundanceByAge/ref_value_percent_abund-100) %>%
      dplyr::select(-ref_value_percent_abund)
    list_relative_biomass[[sp]] <- rel_d_i_i_perc
    rel_biomassplot[[sp]] <- ggplot(rel_d_i_i_perc, aes(y=relative_abundance_percent, x=Fmax)) +
      geom_point(aes(colour = Iszero))+
      labs(title = "Relative biomass per fishing pressure", subtitle = sp)
  }
  
}


## Estimate the % in decline (in biomass and abundance, at first age at reproduction) from Fmax=1.5 to Fmax = 0.5
df_relative_N <- do.call(rbind, list_relative_biomass)
summary_df_relative_abunance <- df_relative_N %>%
  group_by(Fmax, Species)%>%
  summarize(percent_difference_withFmax0.5 = mean(relative_abundance_percent, na.rm=T))
summary_df_relative_abunance <- summary_df_relative_abunance %>%
  dplyr::mutate(percent_difference_withFmax0.5=round(percent_difference_withFmax0.5, 2))
write.csv(summary_df_relative_abunance, file=paste0(PATHOUT, "/Clean_for_paper/Abund_declinetrends.csv"))

## Estimate the % in decline (in Ne) from Fmax=1.5 to Fmax = 0.5
list_relative_Ne <- list()
for (nm in "data_Ne"){
  for (sp in species){
    cat(sp)
    ## load and filter data
    load(file=paste0(pathoutputRDATA, "/Complete_INIT2", sp, ".RData"))
    min_age <- min(as.numeric(complete[[2]]$age), na.rm=T)
    final_dataset_age = complete[[2]]
    rm(complete)

    # remove ages that are not included
    final_dataset_age$age <- as.numeric(final_dataset_age$age)
    final_dataset_age$age[which(is.na(final_dataset_age$age))] <- min_age
    final_dataset_age <- final_dataset_age[which(final_dataset_age$age == min_age),]
    final_dataset_age <- pivot_wider(final_dataset_age, names_from = "variable", values_from = "value")
    ref_data_percent <- final_dataset_age %>%
      filter(Fmax == 0.5 & Fstart == 1) %>%
      dplyr::select(Ne_N_alpha, replicate, age) %>%
      dplyr::rename(ref_value_percent = Ne_N_alpha)
    # Subtract reference from all means
    rel_d_i_i_perc <- final_dataset_age %>%
      dplyr::filter(age %in% c(min_age, NA))%>%
      left_join(ref_data_percent[, -which(colnames(ref_data_percent) %in% c("Fmax", "Fstart"))], by = c("replicate", "age"), relationship = "many-to-one") %>%
      mutate(relative_Ne_percent = 100*Ne_N_alpha/ref_value_percent-100) %>%
      dplyr::select(-ref_value_percent)
    list_relative_Ne[[sp]] <- rel_d_i_i_perc
  }
}
df_relative_Ne <- do.call(rbind, list_relative_Ne)
summary_df_relative_abunance <- df_relative_Ne %>%
  group_by(Fmax, Species)%>%
  summarize(percent_difference_withFmax0.5 = mean(relative_Ne_percent, na.rm=T))
summary_df_relative_abunance <- summary_df_relative_abunance %>%
  dplyr::mutate(percent_difference_withFmax0.5=round(percent_difference_withFmax0.5, 2))
write.csv(summary_df_relative_abunance, file=paste0(PATHOUT, "/Clean_for_paper/Ne_declinetrends.csv"), )





#### PLOT 1- NE/N VS FMAX + Biomass VS FMAX ------------------

traits_to_plot <- list(c("generation_length", "Vk_point", "T/Vk."), 
                       c("Ne_N_alpha", "NT", "Ne_N_alpha\nDivNT"), 
                       c("Ne_N_alpha", "Numerator_N_alpha", "NT"), 
                       c("k_point"))


Ne_plot_list = list()
for (sp in species){
  cat(sp)
  
  ## load and filter data
  load(file=paste0(pathoutputRDATA, "/Complete_INIT2", sp, ".RData"))
  final_dataset_age = complete[[2]]
  min_age <- min(as.numeric(complete[[2]]$age), na.rm=T)
  rm(complete)
  
  # remove ages that are not included
  final_dataset_age$age <- as.numeric(final_dataset_age$age)
  ages_to_remove <- table(final_dataset_age[which(final_dataset_age$variable == "Dead" & final_dataset_age$value ==0),"age"])
  ages_to_remove <- ages_to_remove[which(ages_to_remove == max(ages_to_remove))]
  final_dataset_age <- final_dataset_age[-which(final_dataset_age$age %in% as.numeric(names(ages_to_remove))), ]
  if (any(which(is.na(final_dataset_age$age)))) final_dataset_age$age[which(is.na(final_dataset_age$age))] <- min_age
  
  # # avoid the estimation bias of very populations of very low size
  # if (any(((final_dataset_age[which(final_dataset_age$variable == "abundanceByAge"), "value"])<10))){
  #   final_dataset_age <- pivot_wider(final_dataset_age, names_from = "variable", values_from = "value")
  #   final_dataset_age <- final_dataset_age[-which((final_dataset_age$abundanceByAge)<10), ]
  #   final_dataset_age <- pivot_longer(final_dataset_age, cols = -c("age", "replicate", "Fmax", "Fstart", 'Species', "cohort"), names_to = "variable", values_to = "value")
  # }
  
  ## dataset
  data <-  final_dataset_age %>% filter(Species == sp, Fstart == 1)
  
  if (is.data.frame(data)) data <- list("data_Ne"=data)
  # for the plot
  if (!any(grepl(names(data), pattern="Ne"))) {
    tempo <- c(names(data[[1]])[1], names(data[[1]])[2], names(data[[1]])[3])
    if (length(data)==3) data <- list(do.call(rbind, list(data[[1]][[names(data[[1]])[1]]], data[[2]][[names(data[[1]])[1]]])),
                                      do.call(rbind, list(data[[1]][[names(data[[1]])[2]]], data[[2]][[names(data[[1]])[2]]])),
                                      do.call(rbind, list(data[[1]][[names(data[[1]])[3]]], data[[2]][[names(data[[1]])[3]]])))
    if (length(data)==1)data <- list(do.call(rbind, list(data[[1]][[names(data[[1]])[1]]], data[[2]][[names(data[[1]])[1]]])))
    names(data) <- tempo[1:length(data)]
  }
  scale_v = list()
  for (genetic_drift_data in (1:length(data))[1]){
    
    # select dataset
    summary = data[[genetic_drift_data]]
    if (names(data)[genetic_drift_data] == "tempoNe") {
      if (ncol(data[["data_Nb"]])>6)  summary <- full_join(summary, data[["data_Nb"]][,c(1:4,grep(colnames(data[["data_Nb"]]), pattern="mean_V"))])
      if (ncol(data[["data_Nb"]])==6) summary <- rbind(summary, data[["data_Nb"]][grep(data[["data_Nb"]]$variable, pattern="mean_V"),])
    }
    if (names(data)[genetic_drift_data] == "data_Nb") {
      if (ncol(data[["data_Nb"]])>6)  summary <- full_join(summary, data[["tempoNe"]][,c(1:4,grep(colnames(data[["tempoNe"]]), pattern="generation"))])
      if (ncol(data[["data_Nb"]])==6) summary <- rbind(summary, data[["tempoNe"]][grep(data[["tempoNe"]]$variable, pattern="generation"),])
    }
    
    summary <- as_long_if_needed(summary, columns_keep_pivot = colnames(summary)[-which(colnames(summary) %in% c("replicate", "Fmax", "Fstart", "Species", "age"))])
    df_trait_t <- estimateT_Vk(trait_data_full_i = summary)
    summary <- rbind(summary, df_trait_t)
    
    # add the Ne/NT
    df_trait_t    <- estimateDiv(trait_data_full_i = summary, var=paste0("Ne_N_alpha\nDivNT"))
    summary       <- rbind(summary, df_trait_t)
    
    # add the Ne/Nalpha
    df_trait_t    <- estimateDiv(trait_data_full_i = summary, var=paste0("Ne_N_alpha\nDivNumerator_N_alpha"))
    summary       <- rbind(summary, df_trait_t)
    
    # WEIGHT = summary[which(summary$variable == "abundanceByAge"),]
    # colnames(WEIGHT)[which(colnames(WEIGHT) == "value")] <- "weight"
    
    
    if (nrow(summary)<1) next
    df_trait_sp <- summary %>% filter(Species == sp, Fstart == 1, variable %in% unlist(traits_to_plot))
    df_trait_sp$Fmax <- as.factor(df_trait_sp$Fmax)
    for (relative_to_ref in c(TRUE, FALSE, "proportional")){  
      
      proportional = FALSE
      if (relative_to_ref == "proportional") {
        relative_to_ref = FALSE
        proportional = TRUE
      }
      relative_to_ref = as.logical(relative_to_ref)
      
      for (traits_selected in traits_to_plot) {
        
        df_trait <- df_trait_sp[which(df_trait_sp$variable %in% traits_selected),]
        
        if (sum(is.nan(df_trait$value))>0) df_trait <- df_trait[-which(is.nan(df_trait$value)), ]
        if (sum(is.infinite(df_trait$value))>0) df_trait <- df_trait[-which(is.infinite(df_trait$value)), ]
        
        # conditions to stop
        if (nrow(df_trait)==0) next
        if (all(is.na(df_trait$value))) {
          Ne_plot_list[[paste0(as.character(relative_to_ref), as.character(proportional))]][[sp]][[traits_selected]] <- NULL
          next}
        
        # Use the function
        # if (length(var_i)>2){
        k=1
        id_to_plot = NULL
        # for (i in 1:c((length(var_i)+1)/3)){
        # id_to_plot = k:(k+2)
        # if (traits_selected == "N&T/Vk") id_to_plot=c(which(var_i == "T/Vk"),
        #                                               ifelse(i==1, which(var_i != "T/Vk")[1], which(var_i != "T/Vk")[2]))
        df_trait$variable <- improve_vriable_names(df_trait$variable)
        t <- improve_vriable_names(traits_selected)
        if (all(is.na(df_trait$value))) next
        result <- plot_trait_vs_fmax(df_trait_sp = df_trait, trait_name = t, age_dim = F,#weight = WEIGHT,
                                     Time1, Time2, relative_to_ref = relative_to_ref, spec=sp, proportional=proportional,
                                     position_legend = "right", facet = F)
        
        if (all(t=="Nα")) result <- plot_trait_vs_fmax(df_trait_sp = df_trait, trait_name = t, age_dim = F,#weight = WEIGHT,
                                                       Time1, Time2, relative_to_ref = relative_to_ref, spec=sp, proportional=proportional,
                                                       facet = F, position_legend = "right")
        
        Ne_plot_list[[paste0(as.character(relative_to_ref), as.character(proportional))]][[sp]][[paste(traits_selected, collapse="")]]<- result$plot
        # scale_v[[sp]][[paste(t, collapse="")]] <- result$scale
        k=k+2
        
        # } else {
        #   df_trait_sp_i <- df_trait_sp
        #   df_trait_sp_i$variable <- improve_vriable_names(df_trait_sp_i$variable)
        #   t <- improve_vriable_names(var_i)
        #   result <- plot_trait_vs_fmax(df_trait_sp = df_trait_sp, trait_name = t, age_dim = F, weight = WEIGHT,
        #                                Time1, Time2, relative_to_ref = relative_to_ref, spec=sp, proportional = proportional)
        #   Ne_plot_list[[paste0(as.character(relative_to_ref), as.character(proportional))]][[sp]][[t]]       <- result$plot
        #   scale_v[[sp]][[t]] <- result$scale
        # }
      } # trait selected in the list
    }#sp
  } # genetic drift dim
}# relative
save(Ne_plot_list, file=paste0(pathoutputRDATA, "/Ne_plots.RData"))
load(file=paste0(pathoutputRDATA, "/Ne_plots.RData"))


#### PLOT 2- SSDg and SSDg/N1 and by mortality : only the area under the curve ------------


combined_plot =list()
Time1=150
Time2=160
plot_Deathsdelta_age = list()
plot_Abund_age = list()
plots_by_age_delta= list()
inset.tot.F.Nat = list()
plot_Deathsdelta_age2 = list()
plot_ssdg_age <- list()
for (sp in species){
  cat(sp)
  
  load(file=paste0(pathoutputRDATA, "/Complete_INIT2", sp, ".RData"))
  final_dataset_age = complete[[2]]
  min_age <- min(as.numeric(complete[[2]]$age), na.rm=T)
  rm(complete)
  
  # remove ages that are not included
  final_dataset_age$age <- as.numeric(final_dataset_age$age)
  ages_to_remove <- table(final_dataset_age[which(final_dataset_age$variable == "Dead" & final_dataset_age$value ==0),"age"])
  ages_to_remove <- ages_to_remove[which(ages_to_remove == max(ages_to_remove))]
  final_dataset_age <- final_dataset_age[-which(final_dataset_age$age %in% as.numeric(names(ages_to_remove))), ]
  if (any(which(is.na(df_trait_sp$age)))) df_trait_sp$age[which(is.na(df_trait_sp$age))] <- min_age
  
  # avoid the estimation bias of very populations of very low size
  # if (any(((final_dataset_age[which(final_dataset_age$variable == "abundanceByAge"), "value"])<10))){
  #   final_dataset_age <- pivot_wider(final_dataset_age, names_from = "variable", values_from = "value")
  #   final_dataset_age <- final_dataset_age[-which((final_dataset_age$abundanceByAge)<10), ]
  #   final_dataset_age <- pivot_longer(final_dataset_age, cols = -c("age", "replicate", "Fmax", "Fstart", 'Species', "cohort"), names_to = "variable", values_to = "value")
  # }
  # 
  
  ########### SSDG ----
  df_trait_sp <- final_dataset_age %>% filter(Species == sp, Fstart == 1)
  
  # keep abundance by age variable on the side
  #WEIGHT = df_trait_sp[which(df_trait_sp$variable=="abundanceByAge"),]
  
  # add the SSDgx/Nalpha
  # df_trait_t    <- estimateDiv(trait_data_full_i = df_trait_sp, var=paste0("SSDgx\nDiv", "Numerator_N_alpha"))
  # df_trait_sp   <- rbind(df_trait_sp, df_trait_t)
  
  # all \nDiv Numerator_N_alpha
  var_i <- unique(df_trait_sp$variable)
  for (v in var_i){
    if (v != "Numerator_N_alpha") {
      df_trait_t    <- estimateDiv(trait_data_full_i = df_trait_sp, var=paste0(v, "\nDivNumerator_N_alpha")) # div by population size
      df_trait_sp   <- rbind(df_trait_sp, df_trait_t)
    }
  }
  df_trait_sp$Fmax <- as.factor(df_trait_sp$Fmax)
  
  
  # remove low abundances
  df_trait_sp <- pivot_wider(data = df_trait_sp, names_from = "variable", values_from = "value")
  ID_low_abund <- which(df_trait_sp$abundanceByAge < 10)
  if(length(ID_low_abund)>1) df_trait_sp <- df_trait_sp[-ID_low_abund,]
  df_trait_sp <- pivot_longer(data = df_trait_sp, cols = colnames(df_trait_sp)[-c(1:5)], names_to = "variable", values_to = "value")
  
  
  # produce and store plots
  for (relative_to_ref in c(FALSE)){
    proportional = FALSE
    if (relative_to_ref == "proportional") {
      relative_to_ref = FALSE
      proportional = TRUE
    }
    relative_to_ref = as.logical(relative_to_ref)
    plot_name_ch=c()
    plot_name <- paste0(as.character(relative_to_ref), as.character(proportional))
    if (plot_name=="TRUEFALSE") plot_name_ch <- c("Relative")
    if (plot_name=="FALSEFALSE") plot_name_ch <- c("Absolute")
    if (plot_name=="FALSETRUE") plot_name_ch <- c("Proportional")
    # for the plot
    plot_biom = ifelse(relative_to_ref, "rel_biomassplot", "biomassplot")
    plot_biom = ifelse(relative_to_ref, "biomassplot", "biomassplot")
    
    v= c("SSDgx")
    df_trait_sp_plot <- df_trait_sp[which(df_trait_sp$variable %in% v),]
    ssdg <-  add_SSDg(df=df_trait_sp_plot)#, weight=WEIGHT)
    ssdg$variable <- improve_vriable_names(ssdg$variable)
    df_trait_sp_plot$variable  <- improve_vriable_names(df_trait_sp_plot$variable)
    t <- improve_vriable_names(v)
    ssdg$SSDg <- round(ssdg$SSDg, 1)
    
    inset.tot.F.Nat[[plot_name]][[sp]][[paste(t, collapse="_")]] <- ggplot(ssdg, group = (variable)) +
      geom_boxplot(aes(y=SSDg, x=Fmax)) +
      facet_wrap(~variable, scales="free", ncol=2) + #
      labs(title=sp, y="")+ #, y="SSD\u0262"
      theme(strip.background = element_rect(fill="white"),
            panel.background = element_rect(fill = "white", colour = "grey50"),
            plot.title = element_text(hjust = 0.5),
            axis.text.x = element_text(angle = 60, hjust = 1),
            axis.text = element_text(colour = 1, size = 9),
            panel.border = element_rect(colour = "black", fill=NA),
            strip.text = element_text(size = 10))
    
    result <- plot_trait_vs_fmax(df_trait_sp_plot, trait_name = t, relative_to_ref = relative_to_ref, 
                                 Time1 = Time1, Time2=Time2, #weight=WEIGHT, 
                                 spec = sp, age_dim = T, sp_dim = F, reference_input = NULL, 
                                 FORCElm = FALSE, nrow = 1, facet = T, position_legend = "right")
    plot_ssdg_age[[plot_name]][[sp]][[paste(t, collapse="_")]] <- result
    
    v= c("SSDgx\nDivNumerator_N_alpha")
    df_trait_sp_plot <- df_trait_sp[which(df_trait_sp$variable %in% v),]
    ssdg <-  add_SSDg(df=df_trait_sp_plot)#, weight=WEIGHT)
    ssdg$variable <- improve_vriable_names(ssdg$variable)
    df_trait_sp_plot$variable  <- improve_vriable_names(df_trait_sp_plot$variable)
    t <- improve_vriable_names(v)
    ssdg$SSDg <- round(ssdg$SSDg, 1)
    
    inset.tot.F.Nat[[plot_name]][[sp]][[paste(t, collapse="_")]] <- ggplot(ssdg, group = (variable)) +
      geom_boxplot(aes(y=SSDg, x=Fmax)) +
      facet_wrap(~variable, scales="free", ncol=2) + #
      labs(title=sp, y="")+ #, y="SSD\u0262"
      theme(strip.background = element_rect(fill="white"),
            panel.background = element_rect(fill = "white", colour = "grey50"),
            plot.title = element_text(hjust = 0.5),
            axis.text.x = element_text(angle = 60, hjust = 1),
            axis.text = element_text(colour = 1, size = 9),
            panel.border = element_rect(colour = "black", fill=NA),
            strip.text = element_text(size = 10))
    
    result <- plot_trait_vs_fmax(df_trait_sp_plot, trait_name = t, relative_to_ref = relative_to_ref, 
                                 Time1 = Time1, Time2=Time2, #weight=WEIGHT, 
                                 spec = sp, age_dim = T, sp_dim = F, reference_input = NULL, 
                                 FORCElm = FALSE, nrow = 1, facet = T, position_legend = "right")
    plot_ssdg_age[[plot_name]][[sp]][[paste(t, collapse="_")]] <- result
    
    v= c("SSDgx_F", "SSDgx_Mnat", "SSDgx_Mpred", "SSDgx_Mstar")
    df_trait_sp_plot <- df_trait_sp[which(df_trait_sp$variable %in% v),]
    df_trait_sp_plot$variable <- factor(df_trait_sp_plot$variable, levels = v)
    ssdg <-  add_SSDg(df=df_trait_sp_plot)#, weight=WEIGHT)
    ssdg$variable <- improve_vriable_names(ssdg$variable)
    t <- improve_vriable_names(v)
    ssdg$variable <- factor(ssdg$variable, levels = t)
    ssdg$SSDg <- round(ssdg$SSDg, 1)
    
    inset.tot.F.Nat[[plot_name]][[sp]][[paste(t, collapse="_")]] <- ggplot(ssdg, group = (variable)) +
      geom_boxplot(aes(y=SSDg, x=Fmax)) +
      facet_wrap(~variable, ncol=2) + #scales="free", 
      labs(title=sp, y="")+ #, y="SSD\u0262"
      theme(strip.background = element_rect(fill="white"),
            panel.background = element_rect(fill = "white", colour = "grey50"),
            plot.title = element_text(hjust = 0.5),
            axis.text.x = element_text(angle = 60, hjust = 1),
            axis.text = element_text(colour = 1, size = 9),
            panel.border = element_rect(colour = "black", fill=NA),
            strip.text = element_text(size = 10))    
    
    v= c("SSDgx_F\nDivNumerator_N_alpha", "SSDgx_Mnat\nDivNumerator_N_alpha",
         "SSDgx_Mpred\nDivNumerator_N_alpha", "SSDgx_Mstar\nDivNumerator_N_alpha")
    df_trait_sp_plot <- df_trait_sp[which(df_trait_sp$variable %in% v),]
    ssdg <-  add_SSDg(df=df_trait_sp_plot)#, weight=WEIGHT)
    ssdg$variable <- improve_vriable_names(ssdg$variable)
    t <- improve_vriable_names(v)
    ssdg$variable <- factor(ssdg$variable, levels = t)
    ssdg$SSDg <- round(ssdg$SSDg, 1)
    
    inset.tot.F.Nat[[plot_name]][[sp]][[paste(t, collapse="_")]] <- ggplot(ssdg, group = (variable)) +
      geom_boxplot(aes(y=SSDg, x=Fmax)) +
      facet_wrap(~variable, ncol=2) + #scales="free", 
      labs(title=sp, y="")+ #, y="SSD\u0262"
      theme(strip.background = element_rect(fill="white"),
            panel.background = element_rect(fill = "white", colour = "grey50"),
            plot.title = element_text(hjust = 0.5),
            axis.text.x = element_text(angle = 60, hjust = 1),
            axis.text = element_text(colour = 1, size = 9),
            panel.border = element_rect(colour = "black", fill=NA),
            strip.text = element_text(size = 10))   
    
  }
}#sp
save(inset.tot.F.Nat, file=paste0(pathoutputRDATA, "/inset.tot.F.Nat.RData"))
load(file=paste0(pathoutputRDATA, "/inset.tot.F.Nat.RData"))
save(plot_ssdg_age, file=paste0(pathoutputRDATA, "/plot_ssdg_age.RData"))
load(file=paste0(pathoutputRDATA, "/plot_ssdg_age.RData"))



#### PLOT 3- Species functional traits ? Delta² Dead and Ne VS FMAX and AGE ------------------

Traitplot = list()
Traitplot_age = list()
for (sp in species){
  cat(sp)
  
  load(file=paste0(pathoutputRDATA, "/Complete_INIT2", sp, ".RData"))
  final_dataset_age = complete[[2]]
  min_age <- min(as.numeric(complete[[2]]$age), na.rm=T)
  rm(complete)
  
  # remove ages that are not included
  final_dataset_age$age <- as.numeric(final_dataset_age$age)
  ages_to_remove <- table(final_dataset_age[which(final_dataset_age$variable == "Dead" & final_dataset_age$value ==0),"age"])
  ages_to_remove <- ages_to_remove[which(ages_to_remove == max(ages_to_remove))]
  final_dataset_age <- final_dataset_age[-which(final_dataset_age$age %in% as.numeric(names(ages_to_remove))), ]
  
  #### traits to plot ####
  traits_sp_to_plot = list(c("Numerator_N_alpha","Deltasquare", "mean_of_abundanceMature_per_age\nDivabundanceByAge"), 
                           c("Dead","Deltasquare"), 
                           c("Numerator_N_alpha"),
                           c("mean_of_meanWeightByAge_per_age", "Deltasquare", "b_rate_sc"),
                           c("mean_of_GonadWeight_per_age", "Deltasquare", "b_rate_sc"), 
                           c("Deltasquare", "eggNumber", "mean_of_abundanceMature_per_age"))
  df_trait_sp <- final_dataset_age %>% filter(Species == sp, Fstart == 1)
  rm(final_dataset_age)
  
  #Devide by Nalpha
  for (v in c("eggNumber_sc\nDivNumerator_N_alpha",  "Dead\nDivNumerator_N_alpha")){
    df_trait_t    <- estimateDiv(trait_data_full_i = df_trait_sp, var=v)
    df_trait_sp   <- rbind(df_trait_sp, df_trait_t)
  }
  
  # individual weight data
  data_complete_c ->df_trait_spw
  df_trait_spw <- df_trait_spw[which(df_trait_spw$Species == sp), ]
  df_trait_spw <- df_trait_spw[which(df_trait_spw$Fstart == 1), ]
  df_trait_spw <- df_trait_spw[which(df_trait_spw$variable %in% c("mean_of_meanWeightByAge_per_age", "mean_of_GonadWeight_per_age", "mean_of_abundanceMature_per_age")), ]
  df_trait_spw$replicate <- as.character(df_trait_spw$replicate)
  df_trait_sp <- rbind(df_trait_sp, df_trait_spw)
    
  for (v in c("mean_of_abundanceMature_per_age")){
    df_trait_t    <- estimateDiv(trait_data_full_i = df_trait_sp, var=paste0(v, "\nDivabundanceByAge"))
    df_trait_sp   <- rbind(df_trait_sp, df_trait_t)
  }
  
  # remove low abundances
  df_trait_sp <- pivot_wider(data = df_trait_sp, names_from = "variable", values_from = "value")
  ID_low_abund <- which(df_trait_sp$abundanceByAge < 10)
  if(length(ID_low_abund)>1) df_trait_sp <- df_trait_sp[-ID_low_abund,]
  df_trait_sp <- pivot_longer(data = df_trait_sp, cols = colnames(df_trait_sp)[-c(1:5)], names_to = "variable", values_to = "value")
  
  for (trait_i in traits_sp_to_plot){
    
    cat(trait_i)
    
    # reduce dataset dimensions
    df_trait_sp_plot  <- df_trait_sp[which(df_trait_sp$variable %in% c(trait_i)), ]
    
    for (relative_to_ref in c(FALSE, TRUE, "proportional")){
      
      cat(relative_to_ref)
      
      proportional = FALSE
      if (relative_to_ref == "proportional") {
        relative_to_ref = FALSE
        proportional = TRUE
      }
      relative_to_ref = as.logical(relative_to_ref)
      # for the plot
      plot_biom = ifelse(relative_to_ref, "rel_biomassplot", "biomassplot")
      # plot_biom = ifelse(relative_to_ref, "biomassplot", "biomassplot")
      #
      # ref_data <- df_trait_sp_plot %>%
      #   filter(Fmax == 1 & Fstart == 1) %>%
      #   dplyr::select(variable, value, replicate, cohort, age) %>%
      #   dplyr::rename(ref_value = value)
      
      df_trait_sp_plot <- clean_values(df_trait_sp_plot)
      
      # conditions to stop
      if (nrow(df_trait_sp)==0) next
      if (all(is.na(df_trait_sp_plot$value) | df_trait_sp_plot$value == 0)) {
        next
      }
      
      var_ii                        <- improve_vriable_names(trait_i)
      df_trait_sp_plot$variable     <- improve_vriable_names(df_trait_sp_plot$variable)
      legend = ifelse(length(trait_i)==2, "right", "right")
      df_trait_sp_plot$Fmax <- as.factor(df_trait_sp_plot$Fmax)
      # ref_data$variable             <- improve_vriable_names(ref_data$variable)
      
      df_trait_sp_plot$age[which(is.na(df_trait_sp_plot$age))] <- min(df_trait_sp_plot$age, na.rm = T)
      ssdg <-  add_SSDg(df=df_trait_sp_plot, relative=relative_to_ref)#, weight=WEIGHT)
      ssdg$variable <- improve_vriable_names(ssdg$variable)      
      df_trait_sp_plot$variable  <- improve_vriable_names(df_trait_sp_plot$variable)
      ssdg$variable <- factor(ssdg$variable, levels = var_ii)
      
      Traitplot[[paste0(as.character(relative_to_ref), as.character(proportional))]][[sp]][[paste(var_ii, collapse="")]]  <- 
        ggplot(ssdg) +
        geom_boxplot(aes(y=SSDg, x=Fmax)) +
        facet_wrap(~variable, scales="free", ncol=length(var_ii)) + #
        labs(title=sp, y="")+ #, y="SSD\u0262"
        theme(strip.background = element_rect(fill="white"),
              panel.background = element_rect(fill = "white", colour = "grey50"),
              plot.title = element_text(hjust = 0.5),
              axis.text.x = element_text(angle = 60, hjust = 1),
              axis.text = element_text(colour = 1, size = 9),
              panel.border = element_rect(colour = "black", fill=NA),
              strip.text = element_text(size = 10))
      
      
      if (all(is.na(df_trait_sp_plot$value))) next
      result <- plot_trait_vs_fmax(df_trait_sp = df_trait_sp_plot, trait_name = var_ii,
                                   Time1, Time2, relative_to_ref = relative_to_ref, spec=sp,
                                   ncol=3, nrow=3, 
                                   facet=F, subtitle=F, scale_z = F,
                                   age_dim=T, FORCElm=F, #reference_input=ref_data,
                                   proportional = proportional, #weight = weight_df, 
                                   position_legend="right")
      Traitplot_age[[paste0(as.character(relative_to_ref), as.character(proportional))]][[sp]][[paste(var_ii, collapse="")]] <- result$plot
      
      result <- plot_trait_vs_fmax(df_trait_sp = df_trait_sp_plot, trait_name = var_ii,
                                   Time1, Time2, relative_to_ref = relative_to_ref, spec=sp,
                                   ncol=2, nrow=3, 
                                   legend_justif="right",
                                   facet=T, subtitle=F, scale_z = F,
                                   age_dim=T, FORCElm=F, #reference_input=ref_data,
                                   proportional = proportional) #weight = weight_df)
      Traitplot_age[[paste0(as.character(relative_to_ref), as.character(proportional))]][[sp]][[paste(var_ii, collapse="facet")]] <- result$plot
      
      
      # result <- plot_trait_vs_fmax(df_trait_sp = df_trait_sp_plot, trait_name = var_ii,
      #                              Time1, Time2, relative_to_ref = relative_to_ref, spec=sp,
      #                              ncol=2, nrow=1,
      #                              facet=F, scale_z=F, scale_facet = "free", 
      #                              age_dim=F,  FORCElm=F, #reference_input=ref_data,
      #                              proportional = proportional, #weight = weight_df, 
      #                              position_legend="right")
      
    }# relative to ref
  }# trait
}#sp


save(Traitplot, file=paste0(pathoutputRDATA, "/Traitplot.RData"))
load(file=paste0(pathoutputRDATA, "/Traitplot.RData"))
save(Traitplot_age, file=paste0(pathoutputRDATA, "/Traitplot_age.RData"))
load(file=paste0(pathoutputRDATA, "/Traitplot_age.RData"))



#### PLOT 5- Survivors area under the curve (Nalpha by M) + Weight at age------------

Time1=150
Time2=160
deadinteger <- list()
deadonNinteger <- list()
matrix_compare <- NULL
aliveinteger <- list()
for (sp in species){
  cat(sp)
  
  data_complete_c ->df_trait_sp
  df_trait_sp <- df_trait_sp[which(df_trait_sp$Species == sp), ]
  df_trait_sp <- df_trait_sp[which(df_trait_sp$Fstart == 1), ]
  # add deaths estimates per mortaltiy source
  
  load(file=paste0(pathoutputRDATA, "/Complete_INIT2", sp, ".RData"))
  alpha_age <- min(as.numeric(complete[[2]]$age), na.rm=T)  
  rm(complete)
  
  alpha_agei <- ifelse(alpha_age==1, 2, alpha_age)
  # df_trait_sp <- df_trait_sp[which(df_trait_sp$age<=alpha_agei), ] # from first age to first age at maturation
  
  df_trait_t    <- estimateSum(trait_data_full_i = df_trait_sp, var=paste0("mean_of_F_per_age", "\nsum", "mean_of_Mnat_per_age"))
  df_trait_t$variable    <- "total_M"
  df_trait_sp   <- rbind(df_trait_sp, df_trait_t)
  var_i <- unique(df_trait_sp$variable)
  
  df_trait_t    <- estimateDead(trait_data_full_i = df_trait_sp)
  df_trait_sp   <- rbind(df_trait_sp, df_trait_t)
  
  for (v in var_i){
    if (any(grepl(pattern="mean_of_M", v) | grepl(pattern="mean_of_F", v) | grepl(pattern="total_M", v))) {
      df_trait_t    <- estimateMnumber(trait_data_full_i = df_trait_sp, var=paste0(v))
      df_trait_sp   <- rbind(df_trait_sp, df_trait_t)
    }
  }
  df_trait_t    <- estimatecolSum(trait_data_full_i = df_trait_sp, var=paste0("colSumtot_deaths"))
  df_trait_sp   <- rbind(df_trait_sp, df_trait_t)
  var_i <- unique(df_trait_sp$variable)
  for (v in var_i){
    if (any(grepl(pattern="death", v) & c(grepl(pattern="M", v)) | (grepl(pattern="F_death", v)) | (v=="all_deaths")| (v=="all_deaths_rate") | grepl(pattern="Dead", v))) {
      df_trait_t    <- estimateDiv(trait_data_full_i = df_trait_sp, var=paste0(v, "\nDiv", "colSumtot_deaths"))
      df_trait_t2   <- estimatecolSum(trait_data_full_i = df_trait_t, var=paste0("colSum", (df_trait_t$variable)[1]))
      df_trait_sp   <- rbind(df_trait_sp, df_trait_t)
      df_trait_sp   <- rbind(df_trait_sp, df_trait_t2)
    }
  }
  
  df_trait_t <- estimateSum(df_trait_sp, var =  paste0("F_deaths\nDivcolSumtot_deaths\nsumMnat_deaths\nDivcolSumtot_deaths"), soustraction=T)
  df_trait_sp <- rbind(df_trait_sp, df_trait_t)
  
  df_trait_sp  %>% # replicates ? 
    dplyr::group_by(age, replicate, Species, Fmax, Fstart, variable) %>%
    dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
    dplyr::filter(n > 1L) 
  

  # extract the age at age 1 to substract the dead number from that N1 population size
  df_trait_sp_abund_1 <- df_trait_sp[which(df_trait_sp$variable == "mean_of_abundanceByAge_per_age" & df_trait_sp$age == 1), ]
  df_trait_sp_abund_1 <- df_trait_sp_abund_1 %>%
    group_by(Fmax, Fstart, Species) %>%
    summarise(value = mean(value))
  df_trait_sp_abund_1$Fmax <- as.factor(df_trait_sp_abund_1$Fmax)
  
  # extract Fmax=1 to estimate the proprotional change
  df_trait_sp_w_1 <- df_trait_sp[which(df_trait_sp$variable == "mean_of_meanWeightByAge_per_age" & df_trait_sp$age == 1), ]
  df_trait_sp_w_1 <- df_trait_sp_w_1 %>%
    group_by(Fmax, Fstart, Species) %>%
    summarise(value = mean(value)) %>%
    filter(Fmax==1)
  df_trait_sp_w_1$Fmax <- as.factor(df_trait_sp_w_1$Fmax)
  

  # remove low abundances
  df_trait_sp <- pivot_wider(data = df_trait_sp, names_from = "variable", values_from = "value")
  ID_low_abund <- which(df_trait_sp$mean_of_abundanceByAge_per_age < 10)
  if(length(ID_low_abund)>1) df_trait_sp <- df_trait_sp[-ID_low_abund,]
  df_trait_sp <- pivot_longer(data = df_trait_sp, cols = colnames(df_trait_sp)[-c(1:5)], names_to = "variable", values_to = "value")
  
  
  # produce and store plots
  for (relative_to_ref in c(FALSE)){
    proportional = FALSE
    if (relative_to_ref == "proportional") {
      relative_to_ref = FALSE
      proportional = TRUE
    }
    relative_to_ref = as.logical(relative_to_ref)
    plot_name_ch=c()
    plot_name <- paste0(as.character(relative_to_ref), as.character(proportional))
    if (plot_name=="TRUEFALSE") plot_name_ch <- c("Relative")
    if (plot_name=="FALSEFALSE") plot_name_ch <- c("Absolute")
    if (plot_name=="FALSETRUE") plot_name_ch <- c("Proportional")
    # for the plot
    plot_biom = ifelse(relative_to_ref, "rel_biomassplot", "biomassplot")
    plot_biom = ifelse(relative_to_ref, "biomassplot", "biomassplot")
    
    v= c("F_deaths", "Mnat_deaths","Mpred_deaths", "Mstar_deaths")
    t <- improve_vriable_names(v)
    df_trait_sp_plot <- df_trait_sp[which(df_trait_sp$variable %in% v),]
    df_trait_sp_plot$Fmax <- as.factor(df_trait_sp_plot$Fmax)
    # consider only dead that are happening at ages < age=alpha (maturation)
    df_trait_sp_plot <- df_trait_sp_plot[which(df_trait_sp_plot$age<alpha_age), ] # from first age to first age at maturation -1
    ssdg <-  add_SSDg(df=df_trait_sp_plot)#, weight=WEIGHT)
    ssdg$variable <- improve_vriable_names(ssdg$variable)
    df_trait_sp_plot$variable  <- improve_vriable_names(df_trait_sp_plot$variable)
    ssdg <- full_join(ssdg, df_trait_sp_abund_1)
    ssdg <- na.omit(ssdg)
    ssdg <- ssdg %>% #  = abund_age_1 - dead_by_mortalities_until_age_i
      group_by(variable, Fmax, Fstart, Species) %>%
      summarise(SSDg = value-SSDg) # abundance age 1 - dead number before age a maturity
    t= stringr::str_replace(t, "deaths", "N\u03B1")
    ssdg$variable <- stringr::str_replace(ssdg$variable, 
                                          "deaths", "N\u03B1")
    ssdg$variable <- factor(ssdg$variable, levels = t)
    ssdg$SSDg <- round(ssdg$SSDg, 1)
    
    surv <- ssdg
    
    slope = NULL
    # add slope to plot 
    for (tr in names(table(df_trait_sp_plot$variable))){
      tempo <- lm(df_trait_sp_plot$value[which(df_trait_sp_plot$variable == tr)] ~ as.numeric(df_trait_sp_plot$Fmax[which(df_trait_sp_plot$variable == tr)]))
      slope[tr] <- sign(coef(tempo)[2])
    }
    label <- NULL
    for (l in levels(ssdg$variable)){
      slopel <- as.numeric(slope[l])
      label[l] <- paste0(l, "\n", ifelse(slopel>0, "Positive slope", ifelse(slopel<0, "Negative slope", "No slope")))
    }
    
    aliveinteger[[plot_name]][[sp]][[paste(t, collapse="_")]] <- ggplot(ssdg) +
      geom_boxplot(aes(y=SSDg, x=Fmax)) +
      labs(title=sp, y="")+ #, y="SSD\u0262"
      facet_wrap(~variable, ncol=4, labeller=as_labeller(label)) +
      theme(strip.background = element_rect(fill="white"),
            panel.background = element_rect(fill = "white", colour = "grey50"),
            plot.title = element_text(hjust = 0.5),
            axis.text.x = element_text(angle = 60, hjust = 1),
            axis.text = element_text(colour = 1, size = 9),
            panel.border = element_rect(colour = "black", fill=NA),
            strip.text = element_text(size = 10))
    
    
    
    v= c("Dead", "mean_of_F_per_age", "mean_of_Mnat_per_age","mean_of_Mpred_per_age")
    t <- improve_vriable_names(v)
    df_trait_sp_plot <- df_trait_sp[which(df_trait_sp$variable %in% v),]
    df_trait_sp_plot$Fmax <- as.factor(df_trait_sp_plot$Fmax)
    # consider only dead that are happening at ages < age=alpha (maturation)
    df_trait_sp_plot <- df_trait_sp_plot[which(df_trait_sp_plot$age>=alpha_age), ] # from first age to first age at maturation -1
    slope <- NULL
    ssdg <-  add_SSDg(df=df_trait_sp_plot)#, weight=WEIGHT)
    ssdg$variable <- improve_vriable_names(ssdg$variable)
    df_trait_sp_plot$variable  <- improve_vriable_names(df_trait_sp_plot$variable)
    ssdg$variable <- factor(ssdg$variable, levels = t)
    ssdg$SSDg <- round(ssdg$SSDg, 1)
    ssdg <- na.omit(ssdg)
    
    # ssdg <- rbind(surv, ssdg)
    
    sign = NULL
    # add slope to plot 
    for (tr in names(table(df_trait_sp_plot$variable))){
      tempo <- lm(df_trait_sp_plot$value[which(df_trait_sp_plot$variable == tr)] ~ as.numeric(df_trait_sp_plot$Fmax[which(df_trait_sp_plot$variable == tr)]))
      slope[tr] <- sign(coef(tempo)[2])
    }
    label <- NULL
    for (l in levels(ssdg$variable)){
      slopel <- as.numeric(slope[l])
      label[l] <- paste0(l, "\n", ifelse(slopel>0, "Positive slope", ifelse(slopel<0, "Negative slope", "No slope")))
    }
    
    aliveinteger[[plot_name]][[sp]][[paste(t, collapse="_")]] <- ggplot(ssdg) +
      geom_boxplot(aes(y=SSDg, x=Fmax)) +
      labs(title=sp, y="")+ #, y="SSD\u0262"
      facet_wrap(~variable, ncol=4, scales="free")+#, labeller=as_labeller(label)) +
      theme(strip.background = element_rect(fill="white"),
            panel.background = element_rect(fill = "white", colour = "grey50"),
            plot.title = element_text(hjust = 0.5),
            axis.text.x = element_text(angle = 60, hjust = 1),
            axis.text = element_text(colour = 1, size = 9),
            panel.border = element_rect(colour = "black", fill=NA),
            strip.text = element_text(size = 10))
    
    
    
    v= c("mean_of_meanWeightByAge_per_age")
    t <- improve_vriable_names(v)
    df_trait_sp_plot <- df_trait_sp[which(df_trait_sp$variable %in% v),]
    df_trait_sp_plot$Fmax <- as.factor(df_trait_sp_plot$Fmax)
    ssdg <-  add_SSDg(df=df_trait_sp_plot)#, weight=WEIGHT)
    ssdg$variable <- improve_vriable_names(ssdg$variable)
    df_trait_sp_plot$variable  <- improve_vriable_names(df_trait_sp_plot$variable)
    ssdg <- full_join(ssdg, df_trait_sp_w_1[, -which(colnames(df_trait_sp_w_1) %in% c("Fmax", "Fstart"))])
    ssdg <- ssdg %>% #  = abund_age_1 - dead_by_mortalities_until_age_i
      group_by(variable, Fmax, Fstart, Species) %>%
      summarise(SSDg = (SSDg-value)/(value)) # relative change in weight
    t= stringr::str_replace(t, "deaths", "N\u03B1")
    ssdg$variable <- stringr::str_replace(ssdg$variable, 
                                          "deaths", "N\u03B1")
    ssdg$variable <- factor(ssdg$variable, levels = t)
    
    aliveinteger[[plot_name]][[sp]][[paste(t, collapse="_")]] <- ggplot(ssdg) +
      geom_boxplot(aes(y=SSDg, x=Fmax)) +
      labs(title=sp, y="")+ #, y="SSD\u0262"
      facet_wrap(~variable, ncol=4) +
      theme(strip.background = element_rect(fill="white"),
            panel.background = element_rect(fill = "white", colour = "grey50"),
            plot.title = element_text(hjust = 0.5),
            axis.text.x = element_text(angle = 60, hjust = 1),
            axis.text = element_text(colour = 1, size = 9),
            panel.border = element_rect(colour = "black", fill=NA),
            strip.text = element_text(size = 10))
  }
}#sp
save(aliveinteger, file=paste0(pathoutputRDATA, "/aliveinteger.RData"))
load(file=paste0(pathoutputRDATA, "/aliveinteger.RData"))

# save(deadinteger, file=paste0(pathoutputRDATA, "/deadinteger.RData"))
# load(file=paste0(pathoutputRDATA, "/deadinteger.RData"))
# write.csv(matrix_compare, file=paste0(PATHOUT, "/csv_compare_mortality_proportion.csv"))

#### PLOTS : combine all and SAVE ----------------


path_surv_curve <- paste0("output/PLOTS/selectivity_curves_per_Fmax")
load(file=paste0(path_surv_curve, "dataset_surv-curve.RData"))


## ALL plots RELATIVE -----------------
for (relative_to_ref in c(TRUE)){
  if (!dir.exists(paste0(PATHOUT, "/Clean_for_paper"))) dir.create(paste0(PATHOUT, "/Clean_for_paper"))
  proportional = FALSE
  if (relative_to_ref == "proportional") {
    relative_to_ref = FALSE
    proportional = TRUE
  }
  relative_to_ref <- as.logical(relative_to_ref)
  
  plot_name_ch=c()
  plot_name <- paste0(as.character(relative_to_ref), as.character(proportional))
  if (plot_name=="TRUEFALSE") plot_name_ch <- c("relative")
  if (plot_name=="FALSEFALSE") plot_name_ch <- c("")
  if (plot_name=="FALSETRUE") plot_name_ch <- c("proportional")
  
  plot_biom = ifelse(relative_to_ref, "rel_biomassplot", "biomassplot")
  plot_biom = ifelse(relative_to_ref, "biomassplot", "biomassplot")
  suffix_relative <- plot_name_ch
  
  
  windowsFonts(myCustomWindowsFontName=windowsFont("DejaVu Sans"))
  cairo_pdf(file = paste0(PATHOUT, "/Clean_for_paper/Ne_NT_NeonNT", plot_name_ch, ".pdf"),
            width = 10, height =8, onefile=T)
  inset_no_title <- Ne_plot_list[[plot_name]][["Whiting"]][2] |>
    lapply(\(x) x + labs(title = "Whiting", subtitle="")+ theme(legend.position="none")) |>
    ggarrange(plotlist = _)  
  inset_no_title2 <- Ne_plot_list[[plot_name]][["GreyGurnard"]][2] |>
    lapply(\(x) x + labs(title = "Grey Gurnard", subtitle="")) |>
    ggarrange(plotlist = _)  
  inset_no_title3 <- Ne_plot_list[[plot_name]][["NorwayPout"]][2] |>
    lapply(\(x) x + labs(title = "Norway Pout", subtitle="")+ theme(legend.position="none")) |>
    ggarrange(plotlist = _)  
  inset_no_title4 <- Ne_plot_list[[plot_name]][["HorseMackerel"]][2] |>
    lapply(\(x) x + labs(title = "Horse Mackerel", subtitle="")+ theme(legend.position="none")) |>
    ggarrange(plotlist = _)
  plot_list <- c(inset_no_title, inset_no_title2, inset_no_title3, inset_no_title4)
  # plot_list <- lapply(plot_list, function(x) x + theme_minimal(base_size = 11) + theme(legend.position = "bottom"))
  combined_plot <- ggarrange(plotlist = plot_list,
                             ncol = 2,
                             nrow = 2, common.legend = T)
  combined_plot <- annotate_figure(combined_plot)
  print(combined_plot)
  dev.off()  
  
  windowsFonts(myCustomWindowsFontName=windowsFont("DejaVu Sans"))
  cairo_pdf(file = paste0(PATHOUT, "/Clean_for_paper/TOTNe_Nalpha", plot_name_ch, ".pdf"),
            width = 10, height =8, onefile=T)
  for (sp in  1:ceiling(length(species)/4)){
    if (sp == 4) plot_list <- c(Ne_plot_list[[plot_name]][[(sp-1)*4+1]][3], 
                                Ne_plot_list[[plot_name]][[(sp-1)*4+2]][3],
                                Ne_plot_list[[plot_name]][[(sp-1)*4+3]][3])
    else {
      plot_list <- c(Ne_plot_list[[plot_name]][[(sp-1)*4+1]][3], 
                     Ne_plot_list[[plot_name]][[(sp-1)*4+2]][3],
                     Ne_plot_list[[plot_name]][[(sp-1)*4+3]][3],
                     Ne_plot_list[[plot_name]][[(sp-1)*4+4]][3])
    }
    combined_plot <- ggarrange(plotlist = plot_list,
                               ncol = 2,
                               nrow = 2, common.legend = T)
    combined_plot <- annotate_figure(combined_plot)
    print(combined_plot)
  }
  dev.off()  
  
  windowsFonts(myCustomWindowsFontName=windowsFont("DejaVu Sans"))
  cairo_pdf(file = paste0(PATHOUT, "/Clean_for_paper/Vk_T_TonVk", plot_name_ch, ".pdf"),
            width = 10, height =8, onefile=T)
  inset_no_title <- Ne_plot_list[[plot_name]][["Whiting"]][1] |>
    lapply(\(x) x + labs(title = "Whiting", subtitle="")+ theme(legend.position="none")) |>
    ggarrange(plotlist = _)  
  inset_no_title2 <- Ne_plot_list[[plot_name]][["GreyGurnard"]][1] |>
    lapply(\(x) x + labs(title = "Grey Gurnard", subtitle="")) |>
    ggarrange(plotlist = _)  
  inset_no_title3 <- Ne_plot_list[[plot_name]][["NorwayPout"]][1] |>
    lapply(\(x) x + labs(title = "Norway Pout", subtitle="")+ theme(legend.position="none")) |>
    ggarrange(plotlist = _)  
  inset_no_title4 <- Ne_plot_list[[plot_name]][["HorseMackerel"]][1] |>
    lapply(\(x) x + labs(title = "Horse Mackerel", subtitle="")+ theme(legend.position="none")) |>
    ggarrange(plotlist = _)
  plot_list <- c(inset_no_title, inset_no_title2, inset_no_title3, inset_no_title4)
  # plot_list <- lapply(plot_list, function(x) x + theme_minimal(base_size = 11) + theme(legend.position = "bottom"))
  combined_plot <- ggarrange(plotlist = plot_list,
                             ncol = 2,
                             nrow = 2, common.legend = T)
  combined_plot <- annotate_figure(combined_plot)
  print(combined_plot)
  dev.off()  
  
  
  windowsFonts(myCustomWindowsFontName=windowsFont("DejaVu Sans"))
  cairo_pdf(file = paste0(PATHOUT, "/Clean_for_paper/TOTVk_T_TonVk", plot_name_ch, ".pdf"),
            width = 10, height =8, onefile=T)
  for (sp in  1:ceiling(length(species)/4)){
    if (sp == 4) plot_list <- c(Ne_plot_list[[plot_name]][[(sp-1)*4+1]][1], 
                                Ne_plot_list[[plot_name]][[(sp-1)*4+2]][1],
                                Ne_plot_list[[plot_name]][[(sp-1)*4+3]][1])
    else {
      plot_list <- c(Ne_plot_list[[plot_name]][[(sp-1)*4+1]][1], 
                     Ne_plot_list[[plot_name]][[(sp-1)*4+2]][1],
                     Ne_plot_list[[plot_name]][[(sp-1)*4+3]][1],
                     Ne_plot_list[[plot_name]][[(sp-1)*4+4]][1])
    }
    combined_plot <- ggarrange(plotlist = plot_list,
                               ncol = 2,
                               nrow = 2, common.legend = T)
    combined_plot <- annotate_figure(combined_plot)
    print(combined_plot)
  }
  dev.off()  
  
  windowsFonts(myCustomWindowsFontName=windowsFont("DejaVu Sans"))
  cairo_pdf(file = paste0(PATHOUT, "/Clean_for_paper/TOTNe_NT_NeonNT", plot_name_ch, ".pdf"),
            width = 10, height =8, onefile=T)
  for (sp in  1:ceiling(length(species)/4)){
    if (sp == 4) plot_list <- c(Ne_plot_list[[plot_name]][[(sp-1)*4+1]][2], 
                                Ne_plot_list[[plot_name]][[(sp-1)*4+2]][2],
                                Ne_plot_list[[plot_name]][[(sp-1)*4+3]][2])
    else {
      plot_list <- c(Ne_plot_list[[plot_name]][[(sp-1)*4+1]][2], 
                     Ne_plot_list[[plot_name]][[(sp-1)*4+2]][2],
                     Ne_plot_list[[plot_name]][[(sp-1)*4+3]][2],
                     Ne_plot_list[[plot_name]][[(sp-1)*4+4]][2])
    }
    combined_plot <- ggarrange(plotlist = plot_list,
                               ncol = 2,
                               nrow = 2, common.legend = T)
    combined_plot <- annotate_figure(combined_plot)
    print(combined_plot)
  }
  dev.off() 
  
}
  
  
## ALL PLOTS absolute ----------------
for (relative_to_ref in c(FALSE)){
  if (!dir.exists(paste0(PATHOUT, "/Clean_for_paper"))) dir.create(paste0(PATHOUT, "/Clean_for_paper"))
  proportional = FALSE
  if (relative_to_ref == "proportional") {
    relative_to_ref = FALSE
    proportional = TRUE
  }
  relative_to_ref <- as.logical(relative_to_ref)
  
  plot_name_ch=c()
  plot_name <- paste0(as.character(relative_to_ref), as.character(proportional))
  if (plot_name=="TRUEFALSE") plot_name_ch <- c("relative")
  if (plot_name=="FALSEFALSE") plot_name_ch <- c("")
  if (plot_name=="FALSETRUE") plot_name_ch <- c("proportional")
  
  plot_biom = ifelse(relative_to_ref, "rel_biomassplot", "biomassplot")
  plot_biom = ifelse(relative_to_ref, "biomassplot", "biomassplot")
  suffix_relative <- plot_name_ch
  
  
### SSD
  
  for (sp in c("Whiting", "GreyGurnard", "CommonDab", "EuropeanPlaice")){
    inset_no_title <- inset.tot.F.Nat[[plot_name]][[sp]][1] |>
      lapply(\(x) x + labs(title = "")) |>
      ggarrange(plotlist = _)  
    inset_no_title2 <- inset.tot.F.Nat[[plot_name]][[sp]][2] |>
      lapply(\(x) x + labs(title = "")) |>
      ggarrange(plotlist = _)  
    # inset_no_title3 <- inset.tot.F.Nat[[plot_name]][[sp]][3] |>
    #   lapply(\(x) x + labs(title = "")) |>
    #   ggarrange(plotlist = _)  
    # inset_no_title4 <- inset.tot.F.Nat[[plot_name]][[sp]][4] |>
    #   lapply(\(x) x + labs(title = "")) |>
    #   ggarrange(plotlist = _)
    plot_list <- c(inset_no_title, inset_no_title2)#, inset_no_title3, inset_no_title4)
    # plot_list <- lapply(plot_list, function(x) x + theme_minimal(base_size = 11) + theme(legend.position = "bottom"))
    if (sp=="CommonDab") sp <- "Common Dab"
    if (sp=="EuropeanPlaice") sp <- "European Plaice"
    if (sp=="GreyGurnard") sp <- "Grey Gurnard"
    combined_plot <- ggarrange(plotlist = plot_list,
                               ncol = 2,
                               nrow = 1, common.legend = T)
    assign(x = paste0("plot", sp), annotate_figure(combined_plot,  text_grob(sp, color = "black", face = "bold", size = 13)))
  }  
  windowsFonts(myCustomWindowsFontName=windowsFont("DejaVu Sans"))
  cairo_pdf(file = paste0(PATHOUT, "/Clean_for_paper/SSDgonN_final_", plot_name_ch, ".pdf"),
            width = 10.5, height =6, onefile=T)
  par(mfrow=c(2,2))
  combined_plot <- ggarrange(plotlist = list(get(paste0("plotWhiting")), get(paste0("plotGrey Gurnard")), get(paste0("plotCommon Dab")), get(paste0("plotEuropean Plaice"))),
            ncol = 2,
            nrow = 2, common.legend = T) 
  print(combined_plot)
  dev.off()  
    
  for (sp in species){
    inset_no_title <- inset.tot.F.Nat[[plot_name]][[sp]][1] |>
      lapply(\(x) x + labs(title = "")) |>
      ggarrange(plotlist = _)  
    inset_no_title2 <- inset.tot.F.Nat[[plot_name]][[sp]][2] |>
      lapply(\(x) x + labs(title = "")) |>
      ggarrange(plotlist = _)  
    plot_list <- c(inset_no_title, inset_no_title2)#, inset_no_title3, inset_no_title4)
    combined_plot <- ggarrange(plotlist = plot_list,
                               ncol = 2,
                               nrow = 1, common.legend = T)
    assign(x = paste0("plot", sp), annotate_figure(combined_plot,  text_grob(sp, color = "black", face = "bold", size = 13)))
  }  
  windowsFonts(myCustomWindowsFontName=windowsFont("DejaVu Sans"))
  cairo_pdf(file = paste0(PATHOUT, "/Clean_for_paper/TOTSSDgonN_final_", plot_name_ch, ".pdf"),
            width = 10.5, height =6, onefile=T)
  par(mfrow=c(2,2))
  combined_plot <- ggarrange(plotlist = sapply(species, function(x) get(paste0("plot", x))),
                             ncol = 2,
                             nrow = 2, common.legend = T) 
  print(combined_plot)
  dev.off()  
  
  
  
  
### DELTA SQUARE AND EGG AND MATURE
  
  
  windowsFonts(myCustomWindowsFontName=windowsFont("DejaVu Sans"))
  cairo_pdf(file = paste0(PATHOUT, "/Clean_for_paper/TOT_Deltasquare_&mat", plot_name_ch, ".pdf"),
            width = 13, height =6, onefile=T)
  for (sp in  1:ceiling(length(species)/4)){
    if (sp == 4) plot_list <- c(Traitplot$TRUEFALSE[[(sp-1)*4+1]]$`Δ²eggNumberAbundance mature`, 
                                Traitplot$TRUEFALSE[[(sp-1)*4+2]]$`Δ²eggNumberAbundance mature`,
                                Traitplot$TRUEFALSE[[(sp-1)*4+3]]$`Δ²eggNumberAbundance mature`)
    else {
      plot_list <- c(Traitplot$TRUEFALSE[[(sp-1)*4+1]]$`Δ²eggNumberAbundance mature`, 
                     Traitplot$TRUEFALSE[[(sp-1)*4+2]]$`Δ²eggNumberAbundance mature`,
                     Traitplot$TRUEFALSE[[(sp-1)*4+3]]$`Δ²eggNumberAbundance mature`, 
                     Traitplot$TRUEFALSE[[(sp-1)*4+4]]$`Δ²eggNumberAbundance mature`)
    }
    combined_plot <- ggarrange(plotlist = plot_list,
                               ncol = 2,
                               nrow = 2, common.legend = T)
    combined_plot <- annotate_figure(combined_plot)
    print(combined_plot)
  }
  dev.off()  
  
  
  windowsFonts(myCustomWindowsFontName=windowsFont("DejaVu Sans"))
  cairo_pdf(file = paste0(PATHOUT, "/Clean_for_paper/Deltasquare_&mat.pdf"),
            width = 13, height =6, onefile=T)
  inset_no_title <- list(Traitplot$TRUEFALSE[["Whiting"]]$`Δ²eggNumberAbundance mature`)  |>
    lapply(\(x) x + labs(subtitle = "Whiting", title="") + 
             theme(legend.position="none", 
                   plot.title = element_text(hjust = -0.05, vjust=1, face = "bold", size=15))) |>
    ggarrange(plotlist = _)  
  inset_no_title2 <- list(Traitplot$TRUEFALSE[["GreyGurnard"]]$`Δ²eggNumberAbundance mature`) |>
    lapply(\(x) x + labs(title = "", subtitle = "Grey Gurnard") + theme(legend.position="none")) |>
    ggarrange(plotlist = _)  
  inset_no_title3 <- list(Traitplot$TRUEFALSE[["CommonDab"]]$`Δ²eggNumberAbundance mature`) |>
    lapply(\(x) x + labs(title = "", subtitle = "Common Dab") + theme(legend.position="none")) |>
    ggarrange(plotlist = _)  
  inset_no_title4 <- list(Traitplot$TRUEFALSE[["EuropeanPlaice"]]$`Δ²eggNumberAbundance mature`)  |>
    lapply(\(x) x + labs(title = "", subtitle = "European Plaice") + theme(legend.position="none")) |>
    ggarrange(plotlist = _)  
  plot_list <- c(inset_no_title, inset_no_title2, inset_no_title3, inset_no_title4)
  combined_plot <- ggarrange(plotlist = plot_list,
                             ncol = 2,
                             nrow = 2, common.legend = T)
  combined_plot <- annotate_figure(combined_plot)
  print(combined_plot)
  dev.off()  
  
 ### DEAD & MORTALITY SOURCES
  
  windowsFonts(myCustomWindowsFontName=windowsFont("DejaVu Sans"))
  cairo_pdf(file = paste0(PATHOUT, "/Clean_for_paper/TOT_Deaths&M", plot_name_ch, ".pdf"),
            width = 15, height =7, onefile=T)
  for (sp in  1:ceiling(length(species)/4)){
    if (sp == 4) {plot_list <- c(aliveinteger[[plot_name]][[(sp-1)*4+1]]$`Dead_Fishing mortality rate_Natural mortality rate_Predation mortality rate`, 
                                aliveinteger[[plot_name]][[(sp-1)*4+2]]$`Dead_Fishing mortality rate_Natural mortality rate_Predation mortality rate`,
                                aliveinteger[[plot_name]][[(sp-1)*4+3]]$`Dead_Fishing mortality rate_Natural mortality rate_Predation mortality rate`)
    }else {
      plot_list <- c(aliveinteger[[plot_name]][[(sp-1)*4+1]]$`Dead_Fishing mortality rate_Natural mortality rate_Predation mortality rate`, 
                     aliveinteger[[plot_name]][[(sp-1)*4+2]]$`Dead_Fishing mortality rate_Natural mortality rate_Predation mortality rate`,
                     aliveinteger[[plot_name]][[(sp-1)*4+3]]$`Dead_Fishing mortality rate_Natural mortality rate_Predation mortality rate`,
                     aliveinteger[[plot_name]][[(sp-1)*4+4]]$`Dead_Fishing mortality rate_Natural mortality rate_Predation mortality rate`)
    }
    combined_plot <- ggarrange(plotlist = plot_list,
                               ncol = 2,
                               nrow = 2, common.legend = T)
    combined_plot <- annotate_figure(combined_plot)
    print(combined_plot)
  }
  dev.off()  
  

  windowsFonts(myCustomWindowsFontName=windowsFont("DejaVu Sans"))
  cairo_pdf(file = paste0(PATHOUT, "/Clean_for_paper/Deaths&M&mat.pdf"),
            width = 15, height =7, onefile=T)
  inset_no_title <- list(aliveinteger[[plot_name]][["Whiting"]]$`Dead_Fishing mortality rate_Natural mortality rate_Predation mortality rate`)  |>
    lapply(\(x) x + labs(subtitle = "Whiting", title="") + 
             theme(legend.position="none", 
                   plot.title = element_text(hjust = -0.05, vjust=1, face = "bold", size=15))) |>
    ggarrange(plotlist = _)  
  inset_no_title2 <- list(aliveinteger[[plot_name]][["GreyGurnard"]]$`Dead_Fishing mortality rate_Natural mortality rate_Predation mortality rate`) |>
    lapply(\(x) x + labs(title = "", subtitle = "Grey Gurnard") + theme(legend.position="none")) |>
    ggarrange(plotlist = _)  
  inset_no_title3 <- list(aliveinteger[[plot_name]][["CommonDab"]]$`Dead_Fishing mortality rate_Natural mortality rate_Predation mortality rate`) |>
    lapply(\(x) x + labs(title = "", subtitle = "Common Dab") + theme(legend.position="none")) |>
    ggarrange(plotlist = _)  
  inset_no_title4 <- list(aliveinteger[[plot_name]][["EuropeanPlaice"]]$`Dead_Fishing mortality rate_Natural mortality rate_Predation mortality rate`)  |>
    lapply(\(x) x + labs(title = "", subtitle = "European Plaice") + theme(legend.position="none")) |>
    ggarrange(plotlist = _)  
  plot_list <- c(inset_no_title, inset_no_title2, inset_no_title3, inset_no_title4)
  combined_plot <- ggarrange(plotlist = plot_list,
                             ncol = 2,
                             nrow = 2, common.legend = T)
  combined_plot <- annotate_figure(combined_plot)
  print(combined_plot)
  dev.off()  
  
 ### WEIGHT
  
  
  windowsFonts(myCustomWindowsFontName=windowsFont("DejaVu Sans"))
  cairo_pdf(file = paste0(PATHOUT, "/Clean_for_paper/TOTWeight", plot_name_ch, ".pdf"),
            width = 10, height =8, onefile=T)
  for (sp in  1:ceiling(length(species)/4)){
    if (sp == 4) plot_list <- c(aliveinteger[[plot_name]][[(sp-1)*4+1]]$`Individual weight`, 
                                aliveinteger[[plot_name]][[(sp-1)*4+2]]$`Individual weight`,
                                aliveinteger[[plot_name]][[(sp-1)*4+3]]$`Individual weight`)
    else {
      plot_list <- c(aliveinteger[[plot_name]][[(sp-1)*4+1]]$`Individual weight`, 
                     aliveinteger[[plot_name]][[(sp-1)*4+2]]$`Individual weight`,
                     aliveinteger[[plot_name]][[(sp-1)*4+3]]$`Individual weight`,
                     aliveinteger[[plot_name]][[(sp-1)*4+4]]$`Individual weight`)
    }
    combined_plot <- ggarrange(plotlist = plot_list,
                               ncol = 2,
                               nrow = 2, common.legend = T)
    combined_plot <- annotate_figure(combined_plot)
    print(combined_plot)
  }
  dev.off()  
}

