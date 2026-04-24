

#### Initialisation --------------
library(legendry)
library(ggpubr)

# functions and environment
main = "New_gpu_conf.R" # config file
source("osmose-ns-v4.4.1_clim/R/functions/source.R") # load config and key parameters and functions
source("osmose-ns-v4.4.1_clim/R/functions/functions_key_traits.R") # plots against Fmax
source("osmose-ns-v4.4.1_clim/R/functions/auxiliary_for_lists.R") # from list to df
source("osmose-ns-v4.4.1_clim/R/functions/functions_output_SimuRData_into_df_for_plots.R") # to transform dataset from RData to long df
source("osmose-ns-v4.4.1_clim/R/functions/function_data_traits_summary.R")

# paths
path <- c("C:/Users/mbeneat/Documents/osmose/Chap2/osmose-ns/osmose-ns-v4.4.1_clim/output/PLOTS/FINAL_Ne_Htz_mature/")
PATHOUT  = "osmose-ns-v4.4.1_clim/output/PLOTS/KEY_TRAITS_MATURE/"
if (!dir.exists(PATHOUT)) dir.create(PATHOUT, recursive = TRUE)
pathoutputRDATA = paste0(PATHOUT, "/RDATA_MATURE")
if(!dir.exists(pathoutputRDATA)) dir.create(pathoutputRDATA, recursive = T)
base_dir <- paste0(path, "INIT")
folders  <- list.dirs(base_dir, recursive = F) # files INIT2 with raw Ne datasets

folders_traits_y <- "osmose-ns-v4.4.1_clim/output/PLOTS/Trait_variables_by_year_BYSP//RData/"
folders_traits_y <- list.dirs(folders_traits_y, recursive = F)
folders_traits_c <- "osmose-ns-v4.4.1_clim/output/PLOTS/Trait_variables_by_cohort/RData/"
folders_traits_c  <- list.dirs(folders_traits_c, recursive = F)
folders_traits_y_a <- "osmose-ns-v4.4.1_clim/output/PLOTS/Trait_variables_by_year_by_age2/RData/"
folders_traits_y_a <- list.dirs(folders_traits_y_a, recursive = F)
folders_traits_c_a <- "osmose-ns-v4.4.1_clim/output/PLOTS/Trait_variables_by_cohort_by_age2/RData/"
folders_traits_c_a  <- list.dirs(folders_traits_c_a, recursive = F)

#### DATASETS load or format dataset ------------------------
Time1=150
Time2=160
species_df=species
# DATA NE and related traits, dimensions:, replicate, Fstart, Fmax, Species
### df by loading the trait datasets by Fish_press (by year, by cohort, by year by age, by cohort by age) and binding them together
# list_data_traits_ratio <- data_trait_and_genetic_into_df(folders_in = folders, Time1=150, Time2=160, pathoutput=pathoutputRDATA,
#                                                          FORCE = F, species = species_df, RATIO=T)
# 
# # DATA NE + TRAITS dimensions : (age), replicate, Fstart, Fmax, Species
# dataRATIO_full_merged_age <- dataset_attach_all_traits_genetic_drift(pathoutput=pathoutputRDATA, list_outputs=list("ratio"=list_data_traits_ratio),
#                                                                  trait_cv = list("folders_traits_c_a"=folders_traits_c_a, "folders_traits_y_a"=folders_traits_y_a),
#                                                                  FORCE=F, keep_cohort_or_age_dimension = T)
# dataRATIO_full_merged <- dataset_attach_all_traits_genetic_drift(pathoutput=pathoutputRDATA, list_outputs=list("ratio"=list_data_traits_ratio),
#                                                                  trait_cv = list("folders_traits_c"=folders_traits_c, "folders_traits_y"=folders_traits_y),
#                                                                  FORCE=F, keep_cohort_or_age_dimension = F)
# RAW Ne & associated variables, dimensions : cohort, age, replicate, Fstart, Fmax, Species
# ### produce dataset reading multiple .RData to make 1 dataframe, stored as a .RData file ----------
sp_i = "Hake"
FORCE=F
full_dataset = list()
full_dataset_age = list()
# for (f_i in folders){
#   cat(f_i, "\n")
# 
#   Fmax   <- as.numeric(stringr::str_split(stringr::str_split(f_i, pattern="Fmax")[[1]][2], pattern="_Fstart")[[1]][1])
#   Fstart <- as.numeric(stringr::str_split(stringr::str_split(f_i, pattern="Fstart")[[1]][2], pattern="_sp")[[1]][1])
# 
#   if (Fstart != 1) next
# 
#   load(list.files(f_i, full.names = T))
#   dataNe <- INIT2$data_Ne
# 
#   for (sp in names(dataNe)[1]){
#     if (file.exists(paste0(path, "/Complete_INIT2", sp_i, ".RData")) && FORCE !=T){
#       next
#     }else{
# 
#       cat(sp, "\n")
# 
#       f_i_sp <- dataNe[[sp]]
# 
#       for (dim in names(f_i_sp)){
# 
#         f_i_sp_d <- f_i_sp[[dim]]
# 
#         if (is.data.frame(f_i_sp_d) || is.matrix(f_i_sp_d)){
# 
#           if (is.matrix(f_i_sp_d)) f_i_sp_d <- as.data.frame(f_i_sp_d)
#           f_i_sp_d$cohort <- rownames(f_i_sp_d)
#           f_i_sp_d <- pivot_longer(f_i_sp_d, names_to = "replicate", values_to = "value", cols=-cohort)
#           f_i_sp_d$Species <- sp
#           f_i_sp_d$variable <- dim
#           f_i_sp_d$Fmax <- Fmax
#           f_i_sp_d$Fstart <- Fstart
# 
#           if (is.null(full_dataset[[sp]])){
#             full_dataset[[sp]] <- f_i_sp_d
#           } else {
#             full_dataset[[sp]] <- rbind(full_dataset[[sp]], f_i_sp_d)
#           }
# 
#         } # if matrix or df
#         if (!is.data.frame(f_i_sp_d) & !is.matrix(f_i_sp_d)){
# 
#           f_i_sp_d = lapply(names(f_i_sp_d), FUN = function(x) {
#             f_i_sp_d[[x]] <- as.data.frame(f_i_sp_d[[x]])
#             f_i_sp_d[[x]][, "replicate"] = stringr::str_remove(x, pattern="rep_")
#             f_i_sp_d[[x]][, "age"] = rownames(f_i_sp_d[[x]])
#             return(f_i_sp_d[[x]])
#           })
#           f_i_sp_d = do.call(rbind, f_i_sp_d)
# 
#           f_i_sp_d <- pivot_longer(f_i_sp_d, names_to = "cohort", values_to = "value", cols=-c(replicate, age))
#           f_i_sp_d[, "cohort"] = stringr::str_remove(as.matrix(f_i_sp_d[, "cohort"]), pattern="cohort_")
#           f_i_sp_d$Species <- sp
#           f_i_sp_d$variable <- dim
#           f_i_sp_d$Fmax <- stringr::str_split(stringr::str_split(f_i, pattern="Fmax")[[1]][2], pattern="_Fstart")[[1]][1]
#           f_i_sp_d$Fstart <- stringr::str_split(stringr::str_split(f_i, pattern="Fstart")[[1]][2], pattern="_sp")[[1]][1]
# 
#           if (is.null(full_dataset_age[[sp]])){
#             full_dataset_age[[sp]] <- f_i_sp_d
#           } else {
#             full_dataset_age[[sp]] <- rbind(full_dataset_age[[sp]], f_i_sp_d)
#           }
# 
#         } # if list
#       }# variable dimension
#     }# only if file does not exist yet
#   }# species dimension
# }# folder dimension (Fmax*Fstart)
# #----- MERGE raw variables without age dimension AND variables with age dimension ---#
# for (sp_i in names(full_dataset_age)){
#   if (file.exists(paste0(path, "/Complete_INIT2", sp_i, ".RData")) && FORCE !=T){
#     next
#   }else{
# 
#   # produce a dataset without age dimension : weighted mean
#   weight_df <- full_dataset_age[[sp_i]] %>%
#     dplyr::filter(variable == "abundanceByAge") %>%
#     dplyr::select(Species, replicate, Fstart, Fmax, cohort, age, weight_value = value)
#   data_with_w <- full_dataset_age[[sp_i]] %>%
#     dplyr::left_join(weight_df, by = c("Species", "replicate", "cohort", "age", "Fmax", "Fstart"), relationship = "many-to-many")
#   data_with_w <- data_with_w[-which(data_with_w$weight_value == 0), ]
#   data_with_w$age <- as.numeric(data_with_w$age)
#   data_out <- data_with_w %>%
#     dplyr::group_by(Species, replicate, Fstart, Fmax, cohort, variable) %>%
#     dplyr::summarise(
#       value = stats::weighted.mean(value, w = weight_value, na.rm = TRUE),
#       .groups = "drop"
#     )
#   final_dataset <- rbind(full_dataset[[sp_i]], data_out)
# 
#   # produce a dataste with the age dimension
#   full_dataset_add_age <- lapply(names(table(full_dataset_age[[sp_i]]$age)), FUN = function(x){
#     out = as.data.frame(full_dataset[[sp_i]])
#     out$age = x
#     return(out)
#   })
#   full_dataset_add_age <- do.call(rbind, full_dataset_add_age)
# 
#   final_dataset_age <- rbind(full_dataset_add_age, full_dataset_age[[sp_i]])
# 
#   # save all
#   complete = list(final_dataset, final_dataset_age)
# 
#   save(complete, file=paste0(path, "/Complete_INIT2", sp_i, ".RData"))
#   }
# }



# RAW TRAITS, dimensions :  age, replicate, Fish_press, Species, TRAIT
base_path = "osmose-ns-v4.4.1_clim/output/PLOTS/Trait_variables_by_cohort_by_age2"
pathout   = paste0(base_path, "/otherplots")
if (!dir.exists(pathout)) dir.create(pathout, recursive = T)

# load dataset
data_complete_sp = NULL
data_complete_c = NULL
folder_name <- basename(base_path)
suffix <- gsub("Trait_variables_by_", "", folder_name)
suffix <- gsub("_", "_", suffix)  # keep underscores
save_dir <- file.path(dirname(dirname(base_path)), "RDATA")
save_file <- file.path(save_dir, paste0("summary_trait_", suffix, ".RData"))
if (!dir.exists(save_dir)){dir.create(save_dir, recursive = T)}
FORCE=F
if (!file.exists(save_file) || FORCE==T){
  data_complete_c <- data_traits_summary(base_path)
} else {
  data_complete_c = get(load(file=save_file))
}

# avoid the estimation bias of very populations of very low size
# if (any(((data_complete_c[which(data_complete_c$TRAIT == "mean_of_abundanceByAge_per_age"), "value"])<10))){
#   data_complete_c <- pivot_wider(data_complete_c, names_from = "TRAIT", values_from = "value")
#   data_complete_c <- data_complete_c[-which((data_complete_c$mean_of_abundanceByAge_per_age)<10), ]
#   data_complete_c <- pivot_longer(data_complete_c, cols = -c("age", "replicate", "Fish_press", 'Species'), 
#                                   names_to = "TRAIT", values_to = "value")
# }

# removing age 0
data_complete_c <- data_complete_c[-which(data_complete_c$age == 0), ]
# correction
data_complete_c[grep(data_complete_c$TRAIT, pattern="fec"), "TRAIT"] <- "mean_of_eggNumber_per_age"

# adding natural mortality dimension
data_complete_c[which((data_complete_c$TRAIT) == "mean_of_Mnat_per_age"), "TRAIT"] <- "mean_of_Madd_per_age"
data_complete_c <- pivot_wider(data = data_complete_c, names_from = "TRAIT", values_from = "value")
data_complete_c <- data_complete_c[,-grep(colnames(data_complete_c), pattern="mean_of_Dis_per_age")]
data_complete_c <- data_complete_c[,-grep(colnames(data_complete_c), pattern="mean_of_Z_per_age")]
data_complete_c$mean_of_Mnat_per_age <- data_complete_c$mean_of_Madd_per_age + data_complete_c$mean_of_Mpred_per_age + data_complete_c$mean_of_Mstar_per_age
data_complete_c <- pivot_longer(data = data_complete_c, names_to = "TRAIT", values_to = "value", cols=-colnames(data_complete_c)[1:4])

# changing colnames to match expectations
colnames(data_complete_c)[which(colnames(data_complete_c) == "TRAIT")] <- "variable"
data_complete_c$Fmax   <- as.numeric(unlist(stringr::str_remove(unlist(stringr::str_split(data_complete_c$Fish_press, " Fstart"))[seq(1, 2*length(data_complete_c$Fish_press), by=2)], "Fmax=")))
data_complete_c$Fstart <- as.numeric(unlist(stringr::str_split(data_complete_c$Fish_press, " Fstart="))[seq(2, 2*length(data_complete_c$Fish_press), by=2)])
data_complete_c <- (data_complete_c[,-which(colnames(data_complete_c) == "Fish_press")])




#### PLOT 1- NE/N VS FMAX + Biomass VS FMAX ------------------

traits_to_plot <- list(c("mean_generation_length", "mean_Vk_point", "T/Vk."), c("Ne_N_alpha", "mean_Numerator_N_alpha", "Ne_N_alpha\nDivmean_Numerator_N_alpha"), 
                       c("Ne_N_alpha", "mean_Numerator_N_alpha"), c("mean_Numerator_N_alpha", "T/Vk."), 
                       c("Ne_N_alpha\nDivmean_Numerator_N_alpha", "Ne_N_alpha\nDivmean_NT"))


Ne_plot_list = list()
for (sp in species){
  cat(sp)
  
  ## load and filter data
  load(file=paste0(path, "/Complete_INIT2", sp, ".RData"))
  final_dataset_age = complete[[2]]
  final_dataset_age <- final_dataset_age %>%
    dplyr::mutate(dplyr::across(c(cohort, replicate), ~ stringr::str_remove(.x, "^(cohort_|rep_)")))
  final_dataset_age <- final_dataset_age[which(final_dataset_age$cohort %in% Time1:Time2),]
  rm(complete)
  
  # remove ages that are not included
  final_dataset_age$age <- as.numeric(final_dataset_age$age)
  ages_to_remove <- table(final_dataset_age[which(final_dataset_age$variable == "Dead" & final_dataset_age$value ==0),"age"])
  ages_to_remove <- ages_to_remove[which(ages_to_remove == max(ages_to_remove))]
  final_dataset_age <- final_dataset_age[-which(final_dataset_age$age %in% as.numeric(names(ages_to_remove))), ]
  
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
    
    summary <- as_long_if_needed(summary, columns_keep_pivot = colnames(summary)[-which(colnames(summary) %in% c("replicate", "Fmax", "Fstart", "Species", "age", "cohort"))])
    df_trait_t <- estimateT_Vk(trait_data_full_i = summary)
    summary <- rbind(summary, df_trait_t)
    
    # add NT
    df_trait_t    <- estimateX(trait_data_full_i = summary, var=paste0("abundanceByAgexsurvivorship_curve"))
    df_trait_t$variable <- "NT"
    summary       <- rbind(summary, df_trait_t)      
    
    # add the Ne/NT
    df_trait_t    <- estimateDiv(trait_data_full_i = summary, var=paste0("Ne_N_alpha\nDivNT"))
    summary       <- rbind(summary, df_trait_t)
    
    # add the Ne/Nalpha
    df_trait_t    <- estimateDiv(trait_data_full_i = summary, var=paste0("Ne_N_alpha\nDivNumerator_N_alpha"))
    summary       <- rbind(summary, df_trait_t)
    
    WEIGHT = summary[which(summary$variable == "abundanceByAge"),]
    colnames(WEIGHT)[which(colnames(WEIGHT) == "value")] <- "weight"
    
    
    if (nrow(summary)<1) next
    df_trait_sp <- summary %>% filter(Species == sp, Fstart == 1, variable %in% unlist(traits_to_plot))
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
        result <- plot_trait_vs_fmax(df_trait_sp = df_trait, trait_name = t, age_dim = F,weight = WEIGHT,
                                     Time1, Time2, relative_to_ref = relative_to_ref, spec=sp, proportional=proportional,
                                     position_legend = "right")
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
for (sp in species[2]){
  cat(sp)
  
  load(file=paste0(path, "/Complete_INIT2", sp, ".RData"))
  final_dataset_age = complete[[2]]
  final_dataset_age <- final_dataset_age %>%
    dplyr::mutate(dplyr::across(c(cohort, replicate), ~ stringr::str_remove(.x, "^(cohort_|rep_)")))
  final_dataset_age <- final_dataset_age[which(final_dataset_age$cohort %in% Time1:Time2),]
  rm(complete)
  
  # remove ages that are not included
  final_dataset_age$age <- as.numeric(final_dataset_age$age)
  ages_to_remove <- table(final_dataset_age[which(final_dataset_age$variable == "Dead" & final_dataset_age$value ==0),"age"])
  ages_to_remove <- ages_to_remove[which(ages_to_remove == max(ages_to_remove))]
  final_dataset_age <- final_dataset_age[-which(final_dataset_age$age %in% as.numeric(names(ages_to_remove))), ]
  
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
  WEIGHT = df_trait_sp[which(df_trait_sp$variable=="abundanceByAge"),]
  
  # add the SSDgx/Nalpha
  df_trait_t    <- estimateDiv(trait_data_full_i = df_trait_sp, var=paste0("SSDgx\nDiv", "Numerator_N_alpha"))
  df_trait_sp   <- rbind(df_trait_sp, df_trait_t)
  
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
    
    v= c("SSDgx", "SSDgx\nDivNumerator_N_alpha")
    df_trait_sp_plot <- df_trait_sp[which(df_trait_sp$variable %in% v),]
    ssdg <-  add_SSDg(df=df_trait_sp_plot, weight=WEIGHT)
    ssdg$variable <- improve_vriable_names(ssdg$variable)
    t <- improve_vriable_names(v)
    
    inset.tot.F.Nat[[plot_name]][[sp]][[paste(t, collapse="_")]] <- ggplot(ssdg, group = (variable)) +
      geom_boxplot(aes(y=SSDg, x=Fmax)) +
      facet_wrap(~variable, scales="free", ncol=2) +
      labs(title=sp, y="")+ #, y="SSD\u0262"
      theme(strip.background = element_rect(fill="white"),
            panel.background = element_rect(fill = "white", colour = "grey50"),
            plot.title = element_text(hjust = 0.5),
            axis.text.x = element_text(angle = 60, hjust = 1),
            axis.text = element_text(colour = 1, size = 9),
            panel.border = element_rect(colour = "black", fill=NA),
            strip.text = element_text(size = 10))
  }
  
  
  ######### Delta²xMortalitiesDead area under the curves -----------
  
  df_trait_sp <- final_dataset_age %>% filter(Species == sp, Fstart == 1)
  # remove cohort dimension
  df_trait_sp <- df_trait_sp %>%
    mutate(age = as.numeric(age))%>%
    mutate(replicate = as.numeric(replicate)-1)%>% # to make the datasets match : unify the replicate range
    mutate(Fmax = as.numeric(Fmax))%>%
    mutate(Fstart = as.numeric(Fstart))
  weight_df <- build_weight_df(df_trait_sp, weight_var = "abundanceByAge", key_cols = c("age", "Fmax", "Fstart", "Species", "replicate"))
  df_trait_sp <- weighted_reduce(
    df = df_trait_sp,
    weight_df = weight_df,
    join_cols = c("Species", "replicate", "Fmax", "Fstart", "age"),
    group_cols = c("Species", "replicate", "Fstart", "Fmax", "variable", "age"),
    drop_age_zero = TRUE
  )
  
  # keep abundance by age variable on the side
  WEIGHT = df_trait_sp[which(df_trait_sp$variable=="abundanceByAge"),]
  
  # combine all traits: Mortality + Ne traits : both need to be pivot_wider to be combined
  df_trait_sp <- pivot_wider(df_trait_sp, names_from = "variable", values_from = "value")
  data_complete_c_i <- data_complete_c %>%
    filter(Species == sp)
  data_complete_c_i <- pivot_wider(data_complete_c_i, names_from = "variable", values_from = "value")
  df_trait_sp <- left_join(data_complete_c_i,df_trait_sp, by=c("age", "Fmax", "Fstart", "Species", "replicate"))
  if (any(df_trait_sp$age == 0, na.rm=T)) df_trait_sp <- df_trait_sp[-which(df_trait_sp$age == 0), ]
  df_trait_sp <- pivot_longer(data = df_trait_sp, cols = -which(colnames(df_trait_sp) %in% c("age", "Fmax", "Fstart", "Species", "replicate")),
                              names_to = "variable", values_to = "value")
  df_trait_sp <- na.omit(df_trait_sp)
  
  # total mortality rate 
  df_trait_t    <- estimateSum(trait_data_full_i = df_trait_sp, var=("mean_of_F_per_age\nsummean_of_Mnat_per_age"))
  df_trait_t$variable <- "total_M"
  df_trait_sp   <- rbind(df_trait_sp, df_trait_t)
  
  # add deaths estimates per mortaltiy source
  var_i <- unique(df_trait_sp$variable)
  for (v in var_i){
    if (any(grepl(pattern="mean_of_M", v) | grepl(pattern="mean_of_F", v) | grepl(pattern="total_M", v))) {
      df_trait_t    <- estimateMnumber(trait_data_full_i = df_trait_sp, var=paste0(v))
      df_trait_sp   <- rbind(df_trait_sp, df_trait_t)
    }
  }
  
  df_trait_t <- estimateSum(trait_data_full_i = df_trait_sp, var = paste0("Mnat_deaths\nsumF_deaths"))
  df_trait_t$variable <- "all_deaths"
  df_trait_sp   <- rbind(df_trait_sp, df_trait_t)
  
  df_trait_t    <- estimateX(trait_data_full_i = df_trait_sp, var=paste0("survival_ratexabundanceByAge"))
  #add delta square
  df_trait_t    <- estimateSquare(trait_data_full_i = df_trait_sp, var=paste0("Deltasquare"))
  df_trait_sp   <- rbind(df_trait_sp, df_trait_t)
  
  #add delta square times dead number
  var_i <- unique(df_trait_sp$variable)
  for (v in var_i){
    if (any(grepl(pattern="deaths", v) | (v=="Dead"))) {
      cat(v)
      df_trait_t    <- estimateX(trait_data_full_i = df_trait_sp, var=paste0("Deltasquarex", v))
      df_trait_sp   <- rbind(df_trait_sp, df_trait_t)
    }
  }
  # divide all by Nalpha
  var_i <- unique(df_trait_sp$variable)
  for (v in var_i){
    if (v != "Numerator_N_alpha") {
      cat(v)
      df_trait_t    <- estimateDiv(trait_data_full_i = df_trait_sp, var=paste0(v, "\nDivNumerator_N_alpha")) # div by population size
      df_trait_sp   <- rbind(df_trait_sp, df_trait_t)
    }
  }
  
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
    
    v= names(table(df_trait_sp$variable[which((grepl(df_trait_sp$variable, pattern="death") | grepl(df_trait_sp$variable, pattern="Dead")) &
                                                !grepl(df_trait_sp$variable, pattern="\nDiv") &
                                                grepl(df_trait_sp$variable, pattern="Deltasquare"))]))
    v <- v[c(grep("Dead", v), grep("nat", v), grep("F", v), grep("pred", v))]
    df_trait_sp_plot <- df_trait_sp[which(df_trait_sp$variable %in% c(v)),]#, "SSDgx", "SSDgx\nDivNumerator_N_alpha")),]
    df_trait_sp_plot$Fmax <- as.factor(df_trait_sp_plot$Fmax)
    df_trait_sp_plot$Fstart <- as.factor(df_trait_sp_plot$Fstart)
    WEIGHT$Fmax <- as.factor(WEIGHT$Fmax)
    WEIGHT$Fstart <- as.factor(WEIGHT$Fstart)
    WEIGHT$replicate <- as.numeric(WEIGHT$replicate)
    ssdg <-  add_SSDg(df=df_trait_sp_plot, weight=WEIGHT)
    ssdg$variable <- improve_vriable_names(ssdg$variable)
    t <- improve_vriable_names(v)
    
    inset.tot.F.Nat[[plot_name]][[sp]][[paste(t, collapse="_")]] <- ggplot(ssdg, group = (variable)) +
      geom_boxplot(aes(y=SSDg, x=Fmax)) +
      facet_wrap(~variable, scales="free", ncol=2, nrow=2) +
      labs(title=sp, y="")+ #, y="SSD\u0262"
      theme(strip.background = element_rect(fill="white"),
            panel.background = element_rect(fill = "white", colour = "grey50"),
            plot.title = element_text(hjust = 0.5),
            axis.text.x = element_text(angle = 60, hjust = 1),
            axis.text = element_text(colour = 1, size = 9),
            panel.border = element_rect(colour = "black", fill=NA),
            strip.text = element_text(size = 10))
    
    
    
    v= names(table(df_trait_sp$variable[which((grepl(df_trait_sp$variable, pattern="death") | grepl(df_trait_sp$variable, pattern="Dead")) &
                                                grepl(df_trait_sp$variable, pattern="\nDiv") &
                                                grepl(df_trait_sp$variable, pattern="Deltasquare"))]))
    v <- v[c(grep("Dead", v), grep("nat", v), grep("F", v), grep("pred", v))]
    df_trait_sp_plot <- df_trait_sp[which(df_trait_sp$variable %in% v),]
    df_trait_sp_plot$Fmax <- as.factor(df_trait_sp_plot$Fmax)
    df_trait_sp_plot$Fstart <- as.factor(df_trait_sp_plot$Fstart)
    WEIGHT$Fmax <- as.factor(WEIGHT$Fmax)
    WEIGHT$Fstart <- as.factor(WEIGHT$Fstart)
    WEIGHT$replicate <- as.numeric(WEIGHT$replicate)
    ssdg <-  add_SSDg(df=df_trait_sp_plot, weight=WEIGHT)
    ssdg$variable <- improve_vriable_names(ssdg$variable)
    t <- improve_vriable_names(v)
    
    inset.tot.F.Nat[[plot_name]][[sp]][[paste(t, collapse="_")]] <- ggplot(ssdg, group = (variable)) +
      geom_boxplot(aes(y=SSDg, x=Fmax)) +
      facet_wrap(~variable, scales="free", ncol=2, nrow=2) +
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

