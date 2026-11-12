# the idea is to store separately all the key parameters per species into plots
# and to produce one pdf sheet per species with all the key variables to have it all summarized at one glance


#### Initialisation --------------
library(legendry)
library(ggpubr)
library(ggrepel)
library(ggforce)
library(ggnewscale)
library(patchwork)
library(ggbreak) 

# functions and environment
main = "New_gpu_conf.R" # config file
source("R/functions/source.R") # load config and key parameters and functions
source("R/functions/function_data_traits_summary.R") # datasets preparation & loading
source("R/functions/functions_plots.R") # plots 

# paths
PATHOUT  = "output/PLOTS/Ne_on_NT_with_F/KEY_TRAITS_MATURE_removed0inNedata/"#KEY_TRAITS_MATURE_150-150/"
if (!dir.exists(PATHOUT)) dir.create(PATHOUT, recursive = TRUE)
if (!dir.exists(paste0(PATHOUT, "/Clean_for_paper"))) dir.create(paste0(PATHOUT, "/Clean_for_paper"), recursive = TRUE)
pathoutputRDATA = paste0(PATHOUT, "RDATA")
if(!dir.exists(pathoutputRDATA)) dir.create(pathoutputRDATA, recursive = T)
base_path = paste0("output/DATASETS/Ne_on_NT_with_F/Trait_variables_by_cohort_by_age2")
base_path = paste0("output/DATASETS/Ne_on_NT_with_F/Other_variables/")

folders_traits_c <- paste0("output/DATASETS/Ne_on_NT_with_F/Trait_variables_by_cohort/RData/")
folders_traits_c  <- list.dirs(folders_traits_c, recursive = F)
folders_traits_c_a <- paste0("output/DATASETS/Ne_on_NT_with_F/Trait_variables_by_cohort_by_age/RData/")
folders_traits_c_a <- paste0("output/DATASETS/Ne_on_NT_with_F/Other_variables/RData/")
folders_traits_c_a  <- list.dirs(folders_traits_c_a, recursive = F)

#### DATASETS load or format dataset ------------------------
Time1=150
Time2=160
species_df=species
species_imp <- improve_sp_names(species)
sp_i = "Hake"
FORCE=F #***************
dataset="data_parameters"
sp_main = (c("Whiting", "EuropeanPlaice", "AtlanticMackerel", "HorseMackerel"))
sp_main_imp = improve_sp_names(sp_main)

## % produce .RData files in the folder paste0(simu_repo, "/output/PLOTS/Ne_NT_calculated_variables/")
## % 1 file per species, with the dataset at the scale of the cohorts, starting at age = alpha (first mature age)
dataset_Ne_mean_cohorts(species, FORCE, pathdata=paste0("output/DATASETS/Ne_on_NT_with_F/Ne_NT_calculated_variables/"),
                        pathout=pathoutputRDATA, Time1, Time2, dataset = "data_Ne") 

## loading the dataset of traits averaged over cohorts Time1 to Time2, and reorganised to remove age=0, 
# add reassign natural mortality as add mortality, and estimate natural mortality
data_variables <- load_trait_df(base_path = base_path, pathout=pathoutputRDATA, FORCE)

CLEAN = T # if TRUE : remove low abundances 
          # if FALSE : remove data points with low number of replicates (< 5 replicates )
#### DATASETS clean ------------------------
# CLEAN DATASET OF TRAITS VARIABLES
data_variables <- when_variables_missing_age_dim(data_variables)
data_variables <- what_ages_have_null_dead_nb(data_variables)
if (CLEAN){
  data_variables <- when_null_or_low_abund(dataset=data_variables)
  data_variables <- data_variables[which(data_variables$Fstart == 1),]
}


# CLEAN DATASET OF TRAITS PARAMETERS
data_parameters <- NULL
min_age = c()
for (sp in species){
  cat(sp)
  
  ## load and filter data
  load(file=paste0(pathoutputRDATA, "/Complete_INIT2", sp, ".RData"))
  data_t = complete$complete_with_age
  min_age[sp] <- min(as.numeric(as.matrix(data_t[which(data_t$variable=="survivorship_curve"), "age"])), na.rm=T)
  rm(complete)
  data_t$age <- as.numeric(data_t$age)
  data_t <- when_variables_missing_age_dim(dataset = data_t)
  if (CLEAN){
    # data_t <- what_ages_have_null_dead_nb(data_t)
    data_t <- when_null_or_low_abund(data_t)
  }
  data_t$Fmax <- as.factor(data_t$Fmax) 
  data_t$min_age <- min_age[[sp]]
  data_t = what_ages_have_null_dead_nb(data_t)
  
  ## dataset
  data <-  data_t %>% filter(Species == sp, Fstart == 1)
  if (is.null(data_parameters)){ data_parameters <- data 
  } else data_parameters <- rbind(data_parameters, data)
}




Age_plot_list = list()
traits_to_plot <- list(c("b_rate_sc"))

# select data_parameters
data_plot_1b <- as_long_if_needed(data_parameters, columns_keep_pivot = 
                                    colnames(data_parameters)[-which(colnames(data_parameters) %in% 
                                                                       c("replicate", "Fmax", "Fstart", "Species", "age"))])
tab_rep <- table(data_plot_1b$Fmax, data_plot_1b$age, data_plot_1b$Species)
tab_rep_all <-NULL
# if (!CLEAN){
  # Identify the age & Fmax combinations that do not have > 5 replicates 
  for (i in dimnames(tab_rep)[[3]]){
    tempo <- data.frame(tab_rep[,,i])
    colnames(tempo) <- c("Fmax", "age", "replicates") 
    tempo[, "Species"] <- i
    if(is.null(tab_rep_all)) {tab_rep_all <- tempo
    } else tab_rep_all <- rbind(tab_rep_all, tempo)
  }
  tab_remove <- tab_rep_all[which(tab_rep_all$replicates<6 & tab_rep_all$replicates>0),]
  for (s in unique(tab_remove$Species)){
    tempo <- tab_remove[tab_remove$Species == s, ]
    id_rm <- which(data_plot_1b$Species == s &
                     paste(data_plot_1b$Fmax, data_plot_1b$age) %in% paste(tempo$Fmax, tempo$age))   # exact (Fmax, age) pairs of this species
    data_plot_1b[id_rm, "value"] <- NA                                                              
  }
  
  # Identify the age & Fmax combinations that are not realistic => NA
  data_plot_1b_test <- pivot_wider(data_plot_1b, values_from = "value", names_from = "variable")
  data_plot_1b_test <- data_plot_1b_test %>% 
    ungroup()
  data_plot_1b_abund <- data_plot_1b_test[, c("age", "Fmax", "abundanceByAge", "Species", "replicate")]
  data_plot_1b_abund <- na.omit(data_plot_1b_abund)
  tempo_abund_mean <- data_plot_1b_abund %>%
    group_by(age, Fmax, Species) %>%
    summarise(mean_abund=mean(abundanceByAge, na.rm=T)) %>%
    ungroup()
  tempo_abund_mean <- pivot_wider(tempo_abund_mean, values_from = "mean_abund", names_from = "Fmax")
  tempo_abund_mean <- tempo_abund_mean[order(tempo_abund_mean$Species),]
  tempo_abund_mean[,c("0.5", "0.75", "1", "1.25", "1.5")] <- log10(tempo_abund_mean[,c("0.5", "0.75", "1", "1.25", "1.5")])
  
  age_remove <- list()
  for (sp in species[-which(species == "Sandeel")]){
    for (fmax in c("0.5", "0.75", "1", "1.25", "1.5")){
      data_sp_fmax <- as.matrix(tempo_abund_mean[which(tempo_abund_mean$Species == sp), which(colnames(tempo_abund_mean) == fmax)])   # abs() removed
      age_sp_fmax <- as.matrix(tempo_abund_mean[which(tempo_abund_mean$Species == sp), "age"])
      data_sp_fmax
      # if (all(!is.na(data_sp_fmax))) next
      if (all(data_sp_fmax>2, na.rm = T)) next
      age_remove[[sp]][[fmax]] <- age_sp_fmax[which(data_sp_fmax<2)]
      age_remove[[sp]][[fmax]] <- age_sp_fmax[which(data_sp_fmax<2)]
      data_plot_1b[which((data_plot_1b$Species == sp & 
                            data_plot_1b$Fmax == fmax) & 
                           data_plot_1b$age %in% age_sp_fmax[which(data_sp_fmax<2)]), "value"] <- NA
    }
  }
# }

df_trait_t <- estimateX(trait_data_full_i = data_plot_1b, var = paste0("b_rate", "x", "survivorship_curve"))
data_plot_1b <- rbind(data_plot_1b, df_trait_t)

df_trait_t <- estimateDiv(trait_data_full_i = data_plot_1b, var = paste0("Dead", "\nDiv", "Numerator_N_alpha"))
data_plot_1b <- rbind(data_plot_1b, df_trait_t)

df_trait_t <- estimateX(trait_data_full_i = data_plot_1b, var = paste0("Deltasquare\nxDead\nDivNumerator_N_alpha"))
data_plot_1b <- rbind(data_plot_1b, df_trait_t)

# clean dataset
data_plot_1b <- data_plot_1b %>% filter(variable %in% c(unlist(traits_to_plot), "SSDgx", "survivorship_curve", "Deltasquare", "Dead", "Mtot"))

# condition to stop
if (nrow(data_plot_1b)<1) next

#loop to plot relative and absolute data
slopes_list <- list()
for (relative_to_ref in c(F, T, "proportional")){
  
  cat("\n======>", relative_to_ref)
  # variable for relative or not relative plots
  proportional = FALSE
  if (relative_to_ref == "proportional") {
    relative_to_ref = FALSE
    proportional = TRUE
  }
  relative_to_ref = as.logical(relative_to_ref)
  
  # loop to plot each set of traits 
  for (traits_selected in traits_to_plot) {
    
    cat("\n======>", traits_selected)
    # trait selected
    df_trait_t <- data_plot_1b[which(data_plot_1b$variable %in% traits_selected),]
    
    
    ### *** remove age corresponding to the max age *** ----------------
    # makes more readable plots
    lifespan <- as.data.frame(as.matrix(get_par(conf, "lifespan")))
    colnames(lifespan) <- "lifespan"
    lifespan$ID <- stringr::str_remove(pattern="species.lifespan.sp", rownames(lifespan))
    rownames(lifespan) <- NULL
    name <- as.data.frame(as.matrix(get_par(conf, "species.name")))
    colnames(name) <- "Species"
    name$ID <- stringr::str_remove(pattern="species.name.sp", rownames(name))
    rownames(name) <- NULL
    lifespan <- full_join(name, lifespan)
    lifespan <- lifespan[which(lifespan$Species %in% species), ]
    lifespan$Species <- as.character(lifespan$Species)
    lifespan$lifespan <- as.numeric(lifespan$lifespan)-1
    df_trait_t <- full_join(df_trait_t, lifespan)
    
    df_trait_t <- df_trait_t[which(is.na(df_trait_t$lifespan) |
                                     as.numeric(df_trait_t$age) != as.numeric(df_trait_t$lifespan)),
                             colnames(df_trait_t) != "lifespan"]                 
    ### --------------------------------------------------------------------
    
    
    
    # clean dataset
    if (sum(is.nan(df_trait_t$value))>0) df_trait_t <- df_trait_t[-which(is.nan(df_trait_t$value)), ]
    if (sum(is.infinite(df_trait_t$value))>0) df_trait_t <- df_trait_t[-which(is.infinite(df_trait_t$value)), ]
    
    # exclude certain ages that are hiding main messages
    df_trait_t$age <- as.numeric(df_trait_t$age)
    # Atlantic mackerel: ages > 11 removed
    df_trait_t <- df_trait_t %>% filter(!(Species == "AtlanticMackerel" & age > 11))
    # sprat & horse mackerel: start at the first age where every θF curve has a point
    # this = no ? haddock & whiting: remove the first age (b'x close to 0 -> unstable ratio), start at the second one
    df_trait_t <- df_trait_t %>%
      group_by(Species, variable) %>%                                   # one panel = one Species x variable
      mutate(age_start = {
        ok <- !is.na(value) & value != 0                                # 0 values are removed in plot_trait_vs_fmax (age_dim)
        if (first(Species) %in% c("EuropeanSprat", "HorseMackerel")) {
          max(tapply(age[ok], as.character(Fmax[ok]), min))             # latest first age among θF curves
        #} else if (first(Species) %in% c("Haddock", "Whiting")) {
          #sort(unique(age[ok]))[2]                                      # second age
        } else -Inf                                                     # other species: unchanged
      }) %>%
      filter(is.na(age_start) | age >= age_start) %>%
      ungroup() %>%
      dplyr::select(-age_start)
    ### --------------------------------------------------------------------
    
    # conditions to stop
    if (nrow(df_trait_t)==0) next
    if (all(is.na(df_trait_t$value))) {
      Age_plot_list[[paste0(as.character(relative_to_ref), as.character(proportional))]][[paste(species, collapse="")]][[traits_selected]] <- NULL
      next}
    
    # improve aesthetics
    df_trait_t$variable <- improve_vriable_names(df_trait_t$variable, age=T)
    t <- improve_vriable_names(traits_selected, age=T)
    df_trait_t$variable <- factor(df_trait_t$variable, levels=t)
    df_trait_t$Species <- improve_sp_names(df_trait_t$Species)
    
    # Use the function to plot
    df_trait_t$Species = factor(levels = improve_sp_names(species[order(species)]), df_trait_t$Species)
    result <- plot_trait_vs_fmax(df_trait_sp = df_trait_t, trait_name = t, age_dim = T, spec=species_imp[order(species_imp)],
                                 Time1, Time2, relative_to_ref = relative_to_ref, proportional=proportional,
                                 position_legend = "right", facet = T, subtitle=F, legend_justif="right", 
                                 nrow=5, ncol=3, per_species_scale = T, delta_method = F, scale_z=F)
    slopes_list[[paste(traits_selected, collapse = "")]] <- result$data %>%
      dplyr::filter(!is.na(mean)) %>%                        # Fmax effondrés déjà en NA -> exclus
      dplyr::group_by(Species, variable) %>%
      dplyr::filter(dplyr::n() >= 3) %>%                     # pente sur au moins 3 valeurs de thetaF
      dplyr::summarise(
        slope    = coef(lm(mean ~ Fmax))[2],
        slope_se = summary(lm(mean ~ Fmax))$coefficients[2, 2],
        p_value  = summary(lm(mean ~ Fmax))$coefficients[2, 4],
        n_F      = dplyr::n(),
        .groups  = "drop")
    
    df_trait_t2 <- df_trait_t[which(df_trait_t$Species %in% sp_main_imp), ]
    df_trait_t2$Species = factor(levels = improve_sp_names(sp_main), df_trait_t2$Species)
    result2 <- plot_trait_vs_fmax(df_trait_sp = df_trait_t2, trait_name = t, age_dim = T, spec=sp_main_imp,
                                  Time1, Time2, relative_to_ref = relative_to_ref, proportional=proportional,
                                  position_legend = "right", facet = T, subtitle=F, legend_justif="right", 
                                  nrow=4, ncol=3, per_species_scale = T, delta_method = F, scale_z=F)
    
    # save in a list
    Age_plot_list[[paste0(as.character(relative_to_ref), as.character(proportional))]][[paste(species, collapse="")]][[paste(traits_selected, collapse="")]] <- result$plot
    Age_plot_list[[paste0(as.character(relative_to_ref), as.character(proportional))]][[paste(sp_main_imp, collapse="")]][[paste(traits_selected, collapse="")]] <- result2$plot
  } # trait selected in the list
} # relative ? 
save(Age_plot_list, file=paste0(pathoutputRDATA, "/Age_plots.RData"))
load(file=paste0(pathoutputRDATA, "/Age_plots.RData"))


for (relative_to_ref in c(F, TRUE, "proportional")){
  
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
  cairo_pdf(file = paste0(PATHOUT, "/Clean_for_paper/TOTb_rate_sc", plot_name_ch, ".pdf"),
            width = 10, height =15, onefile=T)
  print(Age_plot_list[[plot_name]]$AtlanticHerringAtlanticMackerelSandeelEuropeanSpratNorwayPoutEuropeanPlaiceCommonSoleSaitheAtlanticCodHaddockHorseMackerelWhitingCommonDabGreyGurnardHake$b_rate_sc)
  dev.off() 
  
  windowsFonts(myCustomWindowsFontName=windowsFont("DejaVu Sans"))
  cairo_pdf(file = paste0(PATHOUT, "/Clean_for_paper/b_rate_sc", plot_name_ch, ".pdf"),
            width = 10, height =12, onefile=T)
  print(Age_plot_list[[plot_name]]$`WhitingEuropean plaiceAtlantic mackerelHorse mackerel`$b_rate_sc)
  dev.off()   
}

