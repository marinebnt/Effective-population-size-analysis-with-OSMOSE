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

when_variables_missing_age_dim <- function(dataset, add_missing_ages = T){
  nvarinit <- length(table(dataset$variable))
  dataset_age <- (as.numeric(as.matrix(dataset[which(dataset$variable %in% c("eggNumber", "mean_of_eggNumberByAge_per_age")), "age"])))
  min_age <- min(dataset_age)
  if (any(as.numeric(as.matrix(dataset$age))<min_age) || any(is.na(as.matrix(dataset$age)))) {
    if(any(is.na(as.matrix(dataset$age))) & add_missing_ages == F) {
      cat("\nRemplacing NA in age by min_age only")
      data_na <- dataset[which(is.na(as.matrix(dataset$age))), ]
      data_na$age <- min_age
      data_not_na <- dataset[-which(is.na(as.matrix(dataset$age))), ]
      dataset <- rbind(data_na, data_not_na)
    }
    if(any(is.na(as.matrix(dataset$age))) & add_missing_ages == T) {
      cat("\nRemplacing NA in age by a looped vector with all the ages from maturation to mortality")
      data_rep <- dataset[which(is.na(as.matrix(dataset$age))), ]
      data_not_rep <- dataset[-which(is.na(as.matrix(dataset$age))), ]
      dataset_repet_i = NULL
      for (x in min_age:max(as.numeric(as.matrix(dataset$age)), na.rm=T)) {
        tempo <- data_rep
        tempo$age <- x
        dataset_repet_i <- rbind(dataset_repet_i, tempo)
      }
      dataset <- rbind(dataset_repet_i, data_not_rep)
    }
    if(any(as.numeric(as.matrix(dataset$age))<min_age)) {
      cat("\nCorrecting error in the dimension of the input dataset : certain outputs have the age < maturation age\n")
      data_wrong_range_age  <- dataset[which(as.numeric(as.matrix(dataset$age))<min_age), "variable"]
      data_good_range_age   <- dataset[-which(dataset$variable %in% as.matrix(unique(data_wrong_range_age))), ]
      data_good_range_age_i <- as.matrix(data_good_range_age$variable)[1]
      data_good_range_age_i <- dataset[which(dataset$variable == data_good_range_age_i), ]
      ages_i_ok <- as.numeric(as.matrix(data_good_range_age_i$age))
      for(vi in names(table(data_wrong_range_age))){
        data_wrong_range_age_i <- dataset[which(dataset$variable == vi), ]
        ages_i <- as.numeric(as.matrix(data_wrong_range_age_i$age))
        if ((max(ages_i)-min(ages_i)) == (max(ages_i_ok)-min(ages_i_ok))){
          data_wrong_range_age_i$age = as.numeric(as.matrix(data_wrong_range_age_i$age))-min(ages_i)+min(ages_i_ok)
        } else cat("=====> cannot be fixed, normal if add_missing_ages == F <=============")
        data_good_range_age <- rbind(data_good_range_age, data_wrong_range_age_i)
      }
      dataset <- data_good_range_age
    }
  }
  nvarend <- length(table(dataset$variable))
  if(nvarend!=nvarend) stop("Dimensions have changed: lost variable(s)")
  return(dataset)
}
when_null_or_low_abund <- function(dataset){
  tab <- table(dataset$variable)
  var_less_rep <- names(tab)[which(tab<max(tab))]
  dataset_t <- pivot_wider(dataset, names_from = "variable", values_from = "value")
  if (any(is.na(as.matrix(dataset_t[,grep(colnames(dataset_t), pattern="abundanceByAge")])))) {
    dataset_t <- dataset_t[-which(is.na(as.matrix(dataset_t[,grep(colnames(dataset_t), pattern="abundanceByAge")]))),]
    dataset   <- pivot_longer(dataset_t, cols = -c("age", "replicate", "Fmax", "Fstart", 'Species'), 
                              names_to = "variable", values_to = "value")
  }
  if (!(length(var_less_rep)==0)) {
    data_t  <- dataset[-which(dataset$variable %in% var_less_rep), ]
  } else {
    data_t <- dataset
  }
  dataset <- rbind(data_t, na.omit(dataset[which(dataset$variable %in% var_less_rep), ]))
  if (any(as.matrix(dataset_t[,grep(colnames(dataset_t), pattern="abundanceByAge")])<100, na.rm=T)){
    dataset_t[which(as.matrix(dataset_t[,grep(colnames(dataset_t), pattern="abundanceByAge")])<100), -which(colnames(dataset_t) %in%
                                                                                                              c("age", "replicate", "Fmax", "Fstart", 'Species'))] <- NA
    dataset <- pivot_longer(dataset_t, cols = -c("age", "replicate", "Fmax", "Fstart", 'Species'), names_to = "variable", values_to = "value")
  }
  if (!(length(var_less_rep)==0)) {
    data_t  <- dataset[-which(dataset$variable %in% var_less_rep), ]
  } else {
    data_t <- dataset
  }
  dataset <- rbind(data_t, na.omit(dataset[which(dataset$variable %in% var_less_rep), ]))
  return(dataset)
}
what_ages_have_null_dead_nb <- function(dataset){
  # remove ages that are not included
  dataset$age <- as.numeric(dataset$age)
  ages_to_remove <- table(dataset[which(dataset$variable == "Dead" & dataset$value ==0),"age"])
  if (dim(ages_to_remove)>0) ages_to_remove <- ages_to_remove[which(ages_to_remove == max(ages_to_remove))]
  if (any(dataset$age %in% as.numeric(names(ages_to_remove)))) dataset <- dataset[-which(dataset$age %in% as.numeric(names(ages_to_remove))), ]
  if (any(which(is.na(dataset$age)))) dataset$age[which(is.na(dataset$age))] <- min_age
  return(dataset)
}


## LOOP TO EXTRACT DATA BY SCENARIO ---------
for (Scenario in c("Chap2", "Chap3")[1]){

  cat("=======> ", Scenario, "\n")
  PATHOUT  = paste0("output/PLOTS/", toupper(Scenario), "/KEY_TRAITS_MATURE_removed0inNedata/")
  if (!dir.exists(PATHOUT)) dir.create(PATHOUT, recursive = TRUE)
  if (!dir.exists(paste0(PATHOUT, "/Clean_for_paper"))) dir.create(paste0(PATHOUT, "/Clean_for_paper"), recursive = TRUE)
  pathoutputRDATA = paste0(PATHOUT, "RDATA")
  if(!dir.exists(pathoutputRDATA)) dir.create(pathoutputRDATA, recursive = T)
  base_path = paste0("output/DATASETS/", Scenario, "/Trait_variables_by_cohort_by_age_withWbyAge/")
  
  ## loading the dataset of traits averaged over cohorts Time1 to Time2, and reorganised to remove age=0, add reassign natural mortality as add mortality, and estimate natural mortality
  data_complete_c <- load_trait_df(base_path = base_path, pathout=pathoutputRDATA, FORCE=F)
  
  #### DATASETS clean ------------------------
  data_complete_c <- when_variables_missing_age_dim(data_complete_c)
  data_complete_c <- when_null_or_low_abund(dataset=data_complete_c)
  

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
    data_complete_c_sp <- data_complete_c[which(data_complete_c$Species == sp),] %>%
      filter(Fstart == 1)
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
      
      
      data_complete_c_sp$variable <- improve_vriable_names(data_complete_c_sp$variable)
      for (dim in names(table(data_complete_c_sp$variable))){
        data_plot <- data_complete_c_sp[which(data_complete_c_sp$variable == dim),]
        result <- plot_trait_vs_fmax(data_plot, title = "",
                                     trait_name = dim, 
                                     relative_to_ref = relative_to_ref, 
                                     Time1 = Time1, 
                                     Time2=Time2, #weight=WEIGHT, 
                                     spec = sp,
                                     age_dim = T, 
                                     sp_dim = F,
                                     reference_input = NULL, 
                                     FORCElm = FALSE, nrow = 1, facet = T, position_legend = "right")
        plot_ssdg_age[[plot_name]][[sp]][[paste(dim, collapse="_")]] <- result$plot
      }
    }
  }
  
  
  for (sp in species) {
    plot = c(plot_ssdg_age[[plot_name]][[sp]])
    combined_plot <- ggarrange(plotlist = plot,
                             ncol = 5,
                             nrow = 4, common.legend = T)
    ggsave(plot = combined_plot, file=paste0(PATHOUT, "Clean_for_paper/TF/", sp, ".png"), create.dir = T, 
           width=25, height=20)
  }
}



## LOOP TO EXTRACT DATA BY SCENARIO ---------
plot_tot <- list()
for (Scenario in c("Chap2", "Chap3")[1]){
  
  cat("=======> ", Scenario, "\n")
  PATHOUT  = paste0("output/PLOTS/", toupper(Scenario), "/KEY_TRAITS_MATURE_removed0inNedata/")
  if (!dir.exists(PATHOUT)) dir.create(PATHOUT, recursive = TRUE)
  if (!dir.exists(paste0(PATHOUT, "/Clean_for_paper"))) dir.create(paste0(PATHOUT, "/Clean_for_paper"), recursive = TRUE)
  pathoutputRDATA = paste0(PATHOUT, "RDATA")
  if(!dir.exists(pathoutputRDATA)) dir.create(pathoutputRDATA, recursive = T)
  base_path = paste0("output/DATASETS/", Scenario, "/Trait_variables_by_cohort_by_age_withWbyAge/")
  
  
  traits_to_plot <- list(c("b_rate", "survivorship_curve"))
  plot_list = list()
  for (sp in species){
    cat(sp)
    
    ## load and filter data
    load(file=paste0(pathoutputRDATA, "/Complete_INIT2", sp, ".RData"))
    final_dataset_age = complete$complete_with_age
    min_age <- min(as.numeric(complete$complete_with_age$age), na.rm=T)
    rm(complete)
    final_dataset_age <- when_variables_missing_age_dim(final_dataset_age)
    final_dataset_age <- what_ages_have_null_dead_nb(final_dataset_age)
    final_dataset_age <- when_null_or_low_abund(final_dataset_age)
    
    ## dataset
    data <-  final_dataset_age %>% filter(Species == sp, Fstart == 1)
    df_trait_sp <- data %>% filter(Species == sp, Fstart == 1, variable %in% unlist(traits_to_plot))
    df_trait_sp$Fmax <- as.factor(df_trait_sp$Fmax)
    
    for (relative_to_ref in c(FALSE)){
      
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
          plot_list[[paste0(as.character(relative_to_ref), as.character(proportional))]][[sp]][[traits_selected]] <- NULL
          next}
        df_trait$variable <- improve_vriable_names(df_trait$variable)
        t <- improve_vriable_names(traits_selected)
        result <- plot_trait_vs_fmax(df_trait_sp = df_trait, trait_name = t, age_dim = T,#weight = WEIGHT,
                                     Time1, Time2, relative_to_ref = relative_to_ref, spec=sp, proportional=proportional,
                                     position_legend = "right", facet = F, subtitle = F, chosen_title = "")
        plot_list[[paste0(as.character(relative_to_ref), as.character(proportional))]][[sp]][[paste(traits_selected, collapse="")]]<- result$plot
      }
    }# sp
  }
  plot_tot[[Scenario]] <- plot_list
}
b <- (plot_list$FALSEFALSE$EuropeanPlaice)
a <- (plot_list$FALSEFALSE$Whiting)
ggarrange(plotlist=c(a, b), nrow=2, common.legend = T, legend="right")

