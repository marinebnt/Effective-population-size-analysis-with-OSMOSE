# the idea is to store separately all the key parameters per species into plots
# and to produce one pdf sheet per species with all the key variables to have it all summarized at one glance




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
PATHOUT  = "osmose-ns-v4.4.1_clim/output/PLOTS/KEY_TRAITS/"
if (!dir.exists(PATHOUT)) dir.create(PATHOUT, recursive = TRUE)
base_dir <- "osmose-ns-v4.4.1_clim/output/PLOTS/FINAL_Ne_Htz2/INIT"
folders  <- list.dirs(base_dir, recursive = F) # files INIT2 with raw Ne datasets

folders_traits_y <- "osmose-ns-v4.4.1_clim/output/PLOTS/Trait_variables_by_year_BYSP//RData/"
folders_traits_y <- list.dirs(folders_traits_y, recursive = F)
folders_traits_c <- "osmose-ns-v4.4.1_clim/output/PLOTS/Trait_variables_by_cohort/RData/"
folders_traits_c  <- list.dirs(folders_traits_c, recursive = F)
folders_traits_y_a <- "osmose-ns-v4.4.1_clim/output/PLOTS/Trait_variables_by_year_by_age2/RData/"
folders_traits_y_a <- list.dirs(folders_traits_y_a, recursive = F)
folders_traits_c_a <- "osmose-ns-v4.4.1_clim/output/PLOTS/Trait_variables_by_cohort_by_age2/RData/"
folders_traits_c_a  <- list.dirs(folders_traits_c_a, recursive = F)

# load or format dataset
pathoutputRDATA = paste0(configDir, "/output/PLOTS/KEY_TRAITS/RDATA")
if(!dir.exists(pathoutputRDATA)) dir.create(pathoutputRDATA, recursive = T)
Time1=150
Time2=160
species_df=species
# df by loading the trait datasets by Fish_press (by year, by cohort, by year by age, by cohort by age) and binding them together
list_data_traits_ratio <- data_trait_and_genetic_into_df(folders_in = folders, Time1=150, Time2=160, pathoutput=pathoutputRDATA,
                                                         FORCE = F, species = species_df, RATIO=T)
dataRATIO_full_merged <- dataset_attach_all_traits_genetic_drift(pathoutput=pathoutputRDATA, list_outputs=list("ratio"=list_data_traits_ratio),
                                                                 trait_cv = list("folders_traits_c"=folders_traits_c, "folders_traits_y"=folders_traits_y),
                                                                 FORCE=F)





Time1=150
Time2=160

PATHOUT  = "osmose-ns-v4.4.1_clim/output/PLOTS/KEY_TRAITS/"
path <- c("C:/Users/mbeneat/Documents/osmose/Chap2/osmose-ns/osmose-ns-v4.4.1_clim/output/PLOTS/FINAL_Ne_Htz2")


### produce dataset reading multiple .RData to make 1 dataframe, stored as a .RData file ----------
sp_i = "Hake"
FORCE=F
if (!file.exists(paste0(path, "/Complete_INIT2", sp_i, ".RData")) || FORCE ==T){
  full_dataset = list()
  full_dataset_age = list()
  for (f_i in folders){
    cat(f_i, "\n")
    
    Fmax   <- stringr::str_split(stringr::str_split(f_i, pattern="Fmax")[[1]][2], pattern="_Fstart")[[1]][1]
    Fstart <- stringr::str_split(stringr::str_split(f_i, pattern="Fstart")[[1]][2], pattern="_sp")[[1]][1]
    
    if (Fstart != 1) next
    
    load(list.files(f_i, full.names = T))
    dataNe <- INIT2$data_Ne
    
    for (sp in names(dataNe)){
      
      cat(sp, "\n")
      
      f_i_sp <- dataNe[[sp]]
      
      for (dim in names(f_i_sp)){
        
        f_i_sp_d <- f_i_sp[[dim]]
        
        if (is.data.frame(f_i_sp_d) || is.matrix(f_i_sp_d)){
          
          if (is.matrix(f_i_sp_d)) f_i_sp_d <- as.data.frame(f_i_sp_d)
          f_i_sp_d$cohort <- rownames(f_i_sp_d)
          f_i_sp_d <- pivot_longer(f_i_sp_d, names_to = "replicate", values_to = "value", cols=-cohort)
          f_i_sp_d$Species <- sp
          f_i_sp_d$variable <- dim
          f_i_sp_d$Fmax <- Fmax
          f_i_sp_d$Fstart <- Fstart
          
          if (is.null(full_dataset[[sp]])){
            full_dataset[[sp]] <- f_i_sp_d
          } else {
            full_dataset[[sp]] <- rbind(full_dataset[[sp]], f_i_sp_d)
          }
          
        } # if matrix or df
        if (!is.data.frame(f_i_sp_d) & !is.matrix(f_i_sp_d)){
          
          f_i_sp_d = lapply(names(f_i_sp_d), FUN = function(x) {
            f_i_sp_d[[x]] <- as.data.frame(f_i_sp_d[[x]])
            f_i_sp_d[[x]][, "replicate"] = stringr::str_remove(x, pattern="rep_")
            f_i_sp_d[[x]][, "age"] = rownames(f_i_sp_d[[x]])
            return(f_i_sp_d[[x]])
          })
          f_i_sp_d = do.call(rbind, f_i_sp_d)
          
          f_i_sp_d <- pivot_longer(f_i_sp_d, names_to = "cohort", values_to = "value", cols=-c(replicate, age))
          f_i_sp_d[, "cohort"] = stringr::str_remove(as.matrix(f_i_sp_d[, "cohort"]), pattern="cohort_")
          f_i_sp_d$Species <- sp
          f_i_sp_d$variable <- dim
          f_i_sp_d$Fmax <- stringr::str_split(stringr::str_split(f_i, pattern="Fmax")[[1]][2], pattern="_Fstart")[[1]][1]
          f_i_sp_d$Fstart <- stringr::str_split(stringr::str_split(f_i, pattern="Fstart")[[1]][2], pattern="_sp")[[1]][1]
          
          if (is.null(full_dataset_age[[sp]])){
            full_dataset_age[[sp]] <- f_i_sp_d
          } else {
            full_dataset_age[[sp]] <- rbind(full_dataset_age[[sp]], f_i_sp_d)
          }
          
        } # if list
      }# variable dimension
    }# species dimension
  }# folder dimension (Fmax*Fstart)
  
  for (sp_i in names(full_dataset_age)){
    # MERGE raw variables without age dimension AND variables with age dimension --------------
    
    # produce a dataset without age dimension : weighted mean
    weight_df <- full_dataset_age[[sp_i]] %>%
      dplyr::filter(variable == "abundanceByAge") %>%
      dplyr::select(Species, replicate, Fstart, Fmax, cohort, age, weight_value = value)
    data_with_w <- full_dataset_age[[sp_i]] %>%
      dplyr::left_join(weight_df, by = c("Species", "replicate", "cohort", "age", "Fmax", "Fstart"), relationship = "many-to-many")
    data_with_w <- data_with_w[-which(data_with_w$weight_value == 0), ]
    data_with_w$age <- as.numeric(data_with_w$age)
    data_out <- data_with_w %>%
      dplyr::group_by(Species, replicate, Fstart, Fmax, cohort, variable) %>%
      dplyr::summarise(
        value = stats::weighted.mean(value, w = weight_value, na.rm = TRUE),
        .groups = "drop"
      )
    final_dataset <- rbind(full_dataset[[sp_i]], data_out)
    
    # produce a dataste with the age dimension
    full_dataset_add_age <- lapply(names(table(full_dataset_age[[sp_i]]$age)), FUN = function(x){
      out = as.data.frame(full_dataset[[sp_i]])
      out$age = x
      return(out)
    })
    full_dataset_add_age <- do.call(rbind, full_dataset_add_age)
    
    final_dataset_age <- rbind(full_dataset_add_age, full_dataset_age[[sp_i]])
    
    # save all
    complete = list(final_dataset, final_dataset_age)
    
    save(complete, file=paste0(path, "/Complete_INIT2", sp_i, ".RData"))
  }
} else {
  cat("you have to load the dataset species by species")
}

traits_sp_to_plot = c("colSumSSDix\nDivNumerator_N_1", "colSumSSDgx\nDivNumerator_N_1", "Vk_point",
                      "colSumDeltasquare\nDivNumerator_N_1", "DeadSum", "colSumDeltasquare\nDivNumerator_N_1\nxDead")
SSDplot = list()
for (sp in species){
  cat(sp)

  load(file=paste0(path, "/Complete_INIT2", sp, ".RData"))
  final_dataset_age = complete[[2]]
  final_dataset_age <- final_dataset_age %>%
    dplyr::mutate(dplyr::across(c(cohort, replicate), ~ stringr::str_remove(.x, "^(cohort_|rep_)")))
  final_dataset_age <- final_dataset_age[which(final_dataset_age$cohort %in% Time1:Time2),]
  rm(complete)

  df_trait_sp <- final_dataset_age %>% filter(Species == sp, Fstart == 1)

  #add delta square
  df_trait_t    <- estimateSquare(trait_data_full_i = df_trait_sp, var=paste0("Deltasquare"))
  df_trait_sp   <- rbind(df_trait_sp, df_trait_t)

  #add sum by age of SSDgx and SSDix
  df_trait_t    <- estimatecolSum(trait_data_full_i = df_trait_sp, var=paste0("colSumSSDix"))
  df_trait_sp   <- rbind(df_trait_sp, df_trait_t)
  df_trait_t    <- estimatecolSum(trait_data_full_i = df_trait_sp, var=paste0("colSumSSDgx"))
  df_trait_sp   <- rbind(df_trait_sp, df_trait_t)
  df_trait_t    <- estimatecolSum(trait_data_full_i = df_trait_sp, var=paste0("colSumDeltasquare"))
  df_trait_sp   <- rbind(df_trait_sp, df_trait_t)

  var_i <- unique(df_trait_sp$variable)
  # compare all traits to their value /N1 to understand where is the increase coming from
  for (v in var_i){
    if (v != "Numerator_N_1") {
      df_trait_t    <- estimateDiv(trait_data_full_i = df_trait_sp, var=paste0(v, "\nDivNumerator_N_1"))
      df_trait_sp   <- rbind(df_trait_sp, df_trait_t)
    }
  }

  #add the SSDgx estimation by hand
  df_trait_t    <- estimateX(trait_data_full_i = df_trait_sp, var=paste0("Deltasquare\nDivNumerator_N_1\nxDead"))
  df_trait_sp   <- rbind(df_trait_sp, df_trait_t)
  df_trait_t    <- estimatecolSum(trait_data_full_i = df_trait_sp, var=paste0("colSumDeltasquare\nDivNumerator_N_1\nxDead"))
  df_trait_sp   <- rbind(df_trait_sp, df_trait_t)

  # reduce dataset dimensions
  var_i = traits_sp_to_plot
  df_trait_sp  <- df_trait_sp[which(df_trait_sp$variable %in% c(var_i, "abundanceByAge")), ]

  # reduce dimension number by averaging (weighted mean)
  weight_df <- build_weight_df(df = df_trait_sp, weight_var = "abundanceByAge",
                               key_cols = c("replicate", "Species", "Fmax", "Fstart"))
  df_trait_sp <- weighted_reduce(
    df = df_trait_sp,
    weight_df = weight_df,
    join_cols = c("Species", "replicate", "Fmax", "Fstart"),
    group_cols = c("Species", "replicate", "Fstart", "Fmax", "variable"),
    drop_age_zero = TRUE
  )

  for (relative_to_ref in c(FALSE, TRUE, "proportional")){
    proportional = FALSE
    if (relative_to_ref == "proportional") {
      relative_to_ref = FALSE
      proportional = TRUE
    }
    relative_to_ref = as.logical(relative_to_ref)
    # for the plot
    plot_biom = ifelse(relative_to_ref, "rel_biomassplot", "biomassplot")
    plot_biom = ifelse(relative_to_ref, "biomassplot", "biomassplot")

    ref_data <- df_trait_sp %>%
      filter(Fmax == 1 & Fstart == 1) %>%
      dplyr::select(variable, value, replicate) %>%
      dplyr::rename(ref_value = value)

    df_trait_sp <- clean_values(df_trait_sp)

    # conditions to stop
    if (nrow(df_trait_sp)==0) next
    if (all(is.na(df_trait_sp$value) | df_trait_sp$value == 0)) {
      next
    }

    # Use the function
    k=1
    id_to_plot = NULL
    for (i in seq_len(2)){
      if (i==1 ) id_to_plot = 1:3
      if (i==2 ) id_to_plot = 4:5

      var_ii <- improve_vriable_names(var_i[id_to_plot])
      # produce plot
      data_split <- df_trait_sp[which(df_trait_sp$variable %in% var_i[id_to_plot] ),]
      data_split$variable <- improve_vriable_names(data_split$variable)
      # data_split <- data_split[, -which(colnames(data_split) == "Fstart")]
      ref_data_split <- ref_data[which(ref_data$variable %in% var_i[id_to_plot]),]
      ref_data_split$variable <- improve_vriable_names(ref_data_split$variable)
      # ref_data_split <- ref_data_split[, -which(colnames(ref_data_split) == "Fstart")]
      # p_l <- list(c(0.35, 0.35), c(0.35, 0.70))[[i]]
      if (all(is.na(data_split$value))) next
      result <- plot_trait_vs_fmax(df_trait_sp = data_split, trait_name = var_ii,
                                   Time1, Time2, relative_to_ref = relative_to_ref, spec=sp,
                                   age_dim=F, reference_input=ref_data_split, FORCElm=F,
                                   proportional = proportional, weight = weight_df, position_legend="auto")

      SSDplot[[paste0(as.character(relative_to_ref), as.character(proportional))]][[sp]][[paste(var_ii, collapse="")]]   <- result$plot
      k=k+2
    }
  } # relative to ref
  # # save plot
  # suffix_relative <- paste0(ifelse(relative_to_ref, "_RELATIVE", ""))
  # if (!dir.exists(paste0(PATHOUT, "/byage", sp))) dir.create(paste0(PATHOUT, "/byage", sp))
  # if (is.null(p[[sp]])) next
  # pdf(file = paste0(PATHOUT, "/byage", sp, "/VxRAW2plot_", "data_Ne", sp, "vsFmax", Fmax_i, suffix_relative, ".pdf"),
  #     width = 15, height = 8)
  # for (trait_plotted in names(p[[sp]])){
  #   plot_list <- c(SSDplot[[as.character(relative_to_ref)]][[sp]][[trait_plotted]], list(get(plot_biom)[[sp]]))
  #   plot_list <- lapply(plot_list, function(x) x + theme_minimal(base_size = 11) + theme(legend.position = "bottom"))
  #   combined_plot <- ggarrange(plotlist = plot_list, ncol = 3,
  #                              nrow = ceiling((length(plot_list))/3),
  #                              common.legend = TRUE, legend = "bottom")
  #   print(combined_plot)
  # }
  # dev.off()
}#sp
save(SSDplot, file=paste0(pathoutputRDATA, "/SSDplots.RData"))
load(file=paste0(pathoutputRDATA, "/SSDplots.RData"))