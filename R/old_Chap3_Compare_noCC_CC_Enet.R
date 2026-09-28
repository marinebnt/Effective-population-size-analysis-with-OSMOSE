# Initialisation --------------
library(MASS)
library(tidyverse)
library(caret)
library(mda)
library(klaR)
library(ggrepel)
library(broom)
library(pls)
library(factoextra)
library(plsVarSel)
library(ggpubr)
library(rstatix)
library(ggplot2)
library(dplyr)
library(hrbrthemes)
library(stargazer)
library(gsubfn)


# functions and environment
main = "New_gpu_conf.R" # config file
source("R/functions/source.R") # load config and key parameters and functions
source("R/functions/function_data_traits_summary.R") # datasets preparation & loading
source("R/functions/functions_plots.R") # plots 




## LOOP TO EXTRACT NET ENERGY DATA  of the two scenarios ---------
for (Scenario in c("Chap2", "Chap3")){
  cat("\nScenario: ", Scenario)
  # paths
  # path <- paste0(simu_repo, "/output/PLOTS/FINAL_Ne_Htz_mature_short/")
  PATHIN  = paste0("output/PLOTS/", toupper(Scenario), "/KEY_TRAITS_MATURE_removed0inNedata/")
  if (!dir.exists(PATHIN)) dir.create(PATHIN, recursive = TRUE)
  if (!dir.exists(paste0(PATHIN, "/Clean_for_paper"))) dir.create(paste0(PATHIN, "/Clean_for_paper"), recursive = TRUE)
  pathoutputRDATA = paste0(PATHIN, "RDATA_Enet")
  if(!dir.exists(pathoutputRDATA)) dir.create(pathoutputRDATA, recursive = T)
  base_path = paste0("output/DATASETS/", (Scenario), "/Trait_Enet")
  
  #### DATASETS load or format dataset ------------------------
  Time1=150
  Time2=160
  species_df=species
  sp_i = "Hake"
  dataset="data_Ne"
  FORCE=F #***************
  removeage0 = F
  
  
  ## loading the dataset of traits averaged over cohorts Time1 to Time2, and reorganised to remove age=0, add reassign natural mortality as add mortality, and estimate natural mortality
  assign(paste0("data_Enet", Scenario), load_trait_df(base_path = base_path, pathout=pathoutputRDATA, FORCE, removeage0))
}

## ORGANISATION NET ENERGY DATA BY SCENARIO ---------

# environmnent paths 
PATHOUT  = paste0("output/PLOTS/CHAP3/KEY_TRAITS_MATURE_removed0inNedata/")
if (!dir.exists(PATHOUT)) dir.create(PATHOUT, recursive = TRUE)
if (!dir.exists(paste0(PATHOUT, "/Clean_for_paper"))) dir.create(paste0(PATHOUT, "/Clean_for_paper"), recursive = TRUE)
pathoutputRDATA = paste0(PATHIN, "RDATA")

# add scenario variable
data_EnetChap3$Scenario = "With CC"
data_EnetChap2$Scenario = "No CC"
data_Enet <- rbind(data_EnetChap3, data_EnetChap2)

# mass-specific Enet
data_Enet <- pivot_wider(data_Enet, names_from = "variable", values_from = "value")
data_Enet$ratioEnet_wbeta <- data_Enet$mean_of_enetByAge_per_age/(data_Enet$mean_of_meanWeightByAge_per_age)^0.75
data_Enet <- pivot_longer(data_Enet, names_to = "variable", values_to = "value", cols= -c("Species", "Scenario", "replicate", "Fmax", "Fstart", "age", "cohort"))

# mean between **cohorts**
data_Enet_mage <- data_Enet %>%
  group_by( Species, Scenario, replicate, Fmax, Fstart, variable, age) %>%
  summarise(mean  = mean(value, na.rm=T)) %>%
  ungroup

# mean between all **ages** (option 1)
data_Enet_mage <- data_Enet_mage %>%
  group_by( Species, Scenario, replicate, Fmax, Fstart, variable) %>%
  summarise(mean  = mean(mean, na.rm=T)) %>%
  ungroup

# mean between **life stages** (option 2)
data_Enet_larv  = list()
data_Enet_juv   = list()
data_Enet_adult = list()
for (sp in species){
  load(file=paste0(pathoutputRDATA, "/Complete_INIT2", sp, ".RData")) # output is a list called 'complete'
  final_dataset_age = complete[[2]] 
  rm(complete)
  egg     <- final_dataset_age[final_dataset_age$variable == "eggNumber", ]
  mat_age <- min(as.numeric(egg$age[which(egg$value>0)]))
  
  data_Enet_larv[[sp]]  <- data_Enet[which(data_Enet$Species == sp & data_Enet$age == 0), ] # => based on get_par(conf, par = "species.larvae.growth.threshold")
  data_Enet_juv[[sp]]   <- data_Enet[which(data_Enet$Species == sp & data_Enet$age %in% c(1:c(mat_age-1))), ]%>%
    group_by(Species, Scenario, variable, Fmax, Fstart, replicate) %>%
    summarise(value  = mean(value, na.rm=T)) %>%
    ungroup
  data_Enet_adult[[sp]] <- data_Enet[which(data_Enet$Species == sp & data_Enet$age>=mat_age), ]%>%
    group_by(Species, Scenario, variable, Fmax, Fstart, replicate) %>%
    summarise(value  = mean(value, na.rm=T)) %>%
    ungroup
}
data_Enet_larv  <- do.call(rbind, data_Enet_larv)
data_Enet_juv   <- do.call(rbind, data_Enet_juv)
data_Enet_adult <- do.call(rbind, data_Enet_adult)

# mean between **replicates**
data_Enet_mage <- data_Enet_mage %>%
  group_by( Species, Scenario, variable, Fmax, Fstart ) %>%
  summarize(across(mean, list(sd = sd, mean = mean), na.rm=T)) %>% 
  ungroup
colnames(data_Enet_mage) <- stringr::str_remove(colnames(data_Enet_mage), pattern = "mean_")
data_Enet_larv <-data_Enet_larv %>%
  group_by(Species, Scenario, variable, Fmax, Fstart) %>%
  summarize(across(value, list(sd = sd, mean = mean), na.rm=T)) %>% 
  ungroup
colnames(data_Enet_larv) <- stringr::str_remove(colnames(data_Enet_larv), pattern = "value_")
data_Enet_juv <- data_Enet_juv %>%
  group_by(Species, Scenario, variable, Fmax, Fstart) %>%
  summarize(across(value, list(sd = sd, mean = mean), na.rm=T)) %>% 
  ungroup
colnames(data_Enet_juv) <- stringr::str_remove(colnames(data_Enet_juv), pattern = "value_")
data_Enet_adult <- data_Enet_adult %>%
  group_by(Species, Scenario, variable, Fmax, Fstart) %>%
  summarize(across(value, list(sd = sd, mean = mean), na.rm=T)) %>% 
  ungroup
colnames(data_Enet_adult) <- stringr::str_remove(colnames(data_Enet_adult), pattern = "value_")

data_Enet_mage_w  <- pivot_wider(data_Enet_mage, names_from = "variable", values_from = c("mean", "sd"))
data_Enet_mage_w  <- data_Enet_mage_w[which(data_Enet_mage_w$Fmax == 1), -which(colnames(data_Enet_mage_w) %in%c("Fmax", "Fstart"))]
data_Enet_larv_w  <- pivot_wider(data_Enet_larv, names_from = "variable", values_from = c("mean", "sd"))
data_Enet_larv_w  <- data_Enet_larv_w[which(data_Enet_larv_w$Fmax == 1), -which(colnames(data_Enet_larv_w) %in%c("Fmax", "Fstart"))]
data_Enet_juv_w   <- pivot_wider(data_Enet_juv, names_from = "variable", values_from = c("mean", "sd"))
data_Enet_juv_w   <- data_Enet_juv_w[which(data_Enet_juv_w$Fmax == 1), -which(colnames(data_Enet_juv_w) %in%c("Fmax", "Fstart"))]
data_Enet_adult_w <- pivot_wider(data_Enet_adult, names_from = "variable", values_from = c("mean", "sd"))
data_Enet_adult_w <- data_Enet_adult_w[which(data_Enet_adult_w$Fmax == 1), -which(colnames(data_Enet_adult_w) %in%c("Fmax", "Fstart"))]

list_Enet_data <- list("Average_all_ages"=data_Enet_mage_w, 
                       "Larvae"=data_Enet_larv_w, 
                       "Juvenile"=data_Enet_juv_w, 
                       "Adult"=data_Enet_adult_w)


## LOOP TO EXTRACT DATASETS effective population size & demographic parameters of the two scenarios  ---------------

Ne_NT_tot2 = NULL
for (Scenario in c("Chap2", "Chap3")){
  cat("\nScenario: ", Scenario)
  
  if (Scenario == "Chap2") main = "New_gpu_conf.R" #config file
  if (Scenario == "Chap3") main = "New_gpu_conf_RCP.R" #config file
  source(paste0(here::here(), "/R/functions/source.R")) #load config and key parameters and functions
  source(paste0(here::here(), "/R/functions/functions_plots.R")) #plots against Fmax
  
  # prepare environment
  Ne_NT2 = NULL 
  cat("=======> ", Scenario, "\n")
  pathoutputRDATA = paste0(PATHIN, "RDATA")
  if(!dir.exists(pathoutputRDATA)) dir.create(pathoutputRDATA, recursive = T)
  base_path = paste0("output/DATASETS/", Scenario, "/Trait_variables_by_cohort_by_age_withWbyAge/")
  
  ## loading the dataset of traits averaged over cohorts Time1 to Time2, and reorganised to remove age=0, 
  ## add reassign natural mortality as add mortality, and estimate natural mortality
  ## => not used in the end
  # data_complete_c <- load_trait_df(base_path = base_path, pathout=pathoutputRDATA, FORCE=F)
  
  for (sp in species){ #load and filter data 
    cat("\nSpecies: ", sp)
    
    load(file=paste0(pathoutputRDATA, "/Complete_INIT2", sp, ".RData")) # output is a list called 'complete'
    final_dataset_age = complete[[2]] 
    rm(complete)
    
    ## Extract only the dimensions I will be using
    final_dataset_age <- final_dataset_age[which(final_dataset_age$Fmax == 1), ]
    NeNT2 <- final_dataset_age[which(final_dataset_age$variable %in% c("Ne_N_alpha", "NT", "generation_length", "Vk_point", "Numerator_N_alpha")), 
                               -which(colnames(final_dataset_age) %in% c("Fstart", "age"))] %>% ungroup()    
    NeNT2[,"Scenario"] <- ifelse(Scenario == "Chap2", "No CC", "With CC") # *********
    
    ## Binding the Species together
    if (is.null(NeNT2)) {
      Ne_NT2 <- pivot_wider(NeNT2, values_from = "value", names_from = "variable") 
      Ne_NT2 |>
        dplyr::summarise(n = dplyr::n(), .by = c(Fmax, replicate, Species, Scenario)) |>
        dplyr::filter(n > 1L) 
    } else {
      Ne_NT2 <- rbind(Ne_NT2, pivot_wider(NeNT2, values_from = "value", names_from = "variable")) 
      dup <- Ne_NT2 |>
        dplyr::summarise(n = dplyr::n(), .by = c(Fmax, replicate, Species, Scenario)) |>
        dplyr::filter(n > 1L) 
      if(nrow(dup)>0) stop()
    }
  }#sp
  
  ## Binding the Scenarios together
  if (is.null(Ne_NT_tot2)) {Ne_NT_tot2 <- Ne_NT2 
  } else Ne_NT_tot2 <- rbind(Ne_NT_tot2, Ne_NT2) 
}#scenario

## ORGANISATION Ne DATA ---------

Ne_NT_tot2$Ne_NT  <- Ne_NT_tot2$Ne_N_alpha/Ne_NT_tot2$NT 
Ne_NT_tot2$'T/Vk' <- Ne_NT_tot2$generation_length/Ne_NT_tot2$Vk_point 

# data for Fmax == 1 & cleaning
Ne_NT_tot2 <- Ne_NT_tot2[which(Ne_NT_tot2$Fmax == 1), ]
Ne_NT_bx   <- pivot_longer(Ne_NT_tot2, cols=-c("Fmax","replicate", "Species", "Scenario"),
                         values_to = "value", names_to = "variable")
dataset_clean_Ne <- Ne_NT_bx[which(Ne_NT_bx$variable %in% c("Ne_N_alpha", "Ne_NT", "NT", 
                                           "Numerator_N_alpha", "Vk_point",
                                           "generation_length",
                                           "T/Vk")), ]


# ANALYSIS -------------


# test data before regression
ggplotdata <- dataset_clean_Ne
# ggplotdata$value <- log(ggplotdata$value)
ggplotdata$id <- 1:nrow(ggplotdata) 
ggplot(ggplotdata)+ 
  geom_density(aes(value))+ 
  facet_wrap(~variable, scales="free")
ggplot(ggplotdata, aes(sample=value))+ 
  stat_qq()+ 
  stat_qq_line()+ 
  facet_wrap(~variable, scales="free")



# 
dataset_clean_Ne <- pivot_wider(dataset_clean_Ne, names_from = "variable", values_from = "value")# summarize 
dataset_clean_Ne <- dataset_clean_Ne %>%
  group_by(Species, Scenario) %>%
  summarize(across(c(`Ne_N_alpha`, `Ne_NT`, NT,
                     `Numerator_N_alpha`, Vk_point, `generation_length`,
                     `T/Vk`),
                   list(sd = sd, mymean = mean))) %>%
  ungroup
colnames(dataset_clean_Ne) <- improve_vriable_names(colnames(dataset_clean_Ne))
colnames(dataset_clean_Ne) <- stringr::str_remove(colnames(dataset_clean_Ne), pattern = "_mymean")
colnames(dataset_clean_Ne) <- stringr::str_remove(colnames(dataset_clean_Ne), pattern = "mymean")




## COMPARE ENET effect on variables depending on life stage -------------------

PATHOUT  = paste0("output/PLOTS/CHAP3/KEY_TRAITS_MATURE_removed0inNedata/Enet/")
dir.create(PATHOUT)
list()-> tablediff
list()-> tabledifflog
list()-> table
list()-> tablelog
for (enet_dimension in names(list_Enet_data)){
  
  ## add the Enet data  ----------------
  data                  <- list_Enet_data[[enet_dimension]]
  dataset_clean_Ne_Enet <- full_join(data, dataset_clean_Ne)
  
  
 ### REGRESSION DATA_NE ~ DATA_ENET*Scenario ------

  ### compare CC and No CC scenarios with regression
  formula = c("`Nₑ/NT` ~ mean_ratioEnet_wbeta*Scenario",
              "`NT` ~ mean_ratioEnet_wbeta*Scenario",
              "`Nₑ` ~ mean_ratioEnet_wbeta*Scenario", 
              "`T` ~ mean_ratioEnet_wbeta*Scenario", 
              "`Vk.` ~ mean_ratioEnet_wbeta*Scenario",
              "`T/Vk.` ~ mean_ratioEnet_wbeta*Scenario",
              "`Nα` ~ mean_ratioEnet_wbeta*Scenario")
  nested_df <- tibble(formula = formula) %>%
    group_by(ID = formula) %>%
    group_modify(~ as_tibble(dataset_clean_Ne_Enet)) %>%
    nest()
  nested_df %>%
    mutate(estimates = data %>% map2(ID, ~ tidy(lm(.y, data = .x)))) %>%
    select(-data) %>%
    unnest()
  nested_df %>%
    mutate(model = ID %>% map2(data, lm),
           stats = model %>% map(augment)) %>%
    select(-data, -model) %>%
    unnest()
  star <- nested_df %>%
    mutate(model = ID %>% map2(data, ~ lm(.x, .y))) %>%
    pull(model) %>%
    stargazer(type = "text", digits = 3, header = F, digit.separator = '', align = T)
  replace_numbers = function(x, cutoff=4, digits=3, scipen=-7) {
    ifelse(nchar(x) < cutoff, x, prettyNum(as.numeric(x), digits=digits, scientific=scipen))
  }
  table[[enet_dimension]] <- gsubfn("([0-9.]+)", ~replace_numbers(x), star)
  
  # log-scale versions of the same formulas (kept, but see FIX-6 for
  # how to pick between raw and log rather than eyeballing both)
  formula = c("log(`Nₑ/NT`) ~ mean_ratioEnet_wbeta*Scenario",
              "log(NT) ~ mean_ratioEnet_wbeta*Scenario",
              "log(Nₑ) ~ mean_ratioEnet_wbeta*Scenario", 
              "log(T) ~ mean_ratioEnet_wbeta*Scenario", 
              "log(`Vk.`) ~ mean_ratioEnet_wbeta*Scenario",
              "log(`T/Vk.`) ~ mean_ratioEnet_wbeta*Scenario",
              "log(Nα) ~ mean_ratioEnet_wbeta*Scenario")
  nested_df <- tibble(formula = formula) %>%
    group_by(ID = formula) %>%
    group_modify(~ as_tibble(dataset_clean_Ne_Enet)) %>%
    nest()
  nested_df %>%
    mutate(estimates = data %>% map2(ID, ~ tidy(lm(.y, data = .x)))) %>%
    select(-data) %>%
    unnest()
  nested_df %>%
    mutate(model = ID %>% map2(data, lm),
           stats = model %>% map(augment)) %>%
    select(-data, -model) %>%
    unnest()
  star <- nested_df %>%
    mutate(model = ID %>% map2(data, ~ lm(.x, .y))) %>%
    pull(model) %>%
    stargazer(type = "text", digits = 3, header = F, digit.separator = '', align = T)
  tablelog[[enet_dimension]] <- gsubfn("([0-9.]+)", ~replace_numbers(x), star)
  
  # compare AIC (lower = better) per variable and log the
  # choice. Add this AFTER the two nested_df blocks above, once both sets
  # of models exist for this enet_dimension.
  model_choice <- tibble(
    variable = c("Ne_NT", "NT", "Ne", "T", "Vk", "T_Vk", "Nalpha"),
    aic_raw  = sapply(list(lm(dataset_clean_Ne_Enet$`Nₑ/NT`   ~ dataset_clean_Ne_Enet$mean_ratioEnet_wbeta*dataset_clean_Ne_Enet$Scenario),
                           lm(dataset_clean_Ne_Enet$`NT`      ~ dataset_clean_Ne_Enet$mean_ratioEnet_wbeta*dataset_clean_Ne_Enet$Scenario),
                           lm(dataset_clean_Ne_Enet$`Nₑ`      ~ dataset_clean_Ne_Enet$mean_ratioEnet_wbeta*dataset_clean_Ne_Enet$Scenario),
                           lm(dataset_clean_Ne_Enet$`T`       ~ dataset_clean_Ne_Enet$mean_ratioEnet_wbeta*dataset_clean_Ne_Enet$Scenario),
                           lm(dataset_clean_Ne_Enet$`Vk.`     ~ dataset_clean_Ne_Enet$mean_ratioEnet_wbeta*dataset_clean_Ne_Enet$Scenario),
                           lm(dataset_clean_Ne_Enet$`T/Vk.`   ~ dataset_clean_Ne_Enet$mean_ratioEnet_wbeta*dataset_clean_Ne_Enet$Scenario),
                           lm(dataset_clean_Ne_Enet$`Nα`      ~ dataset_clean_Ne_Enet$mean_ratioEnet_wbeta*dataset_clean_Ne_Enet$Scenario)),
                      AIC),
    aic_log  = sapply(list(lm(log(dataset_clean_Ne_Enet$`T`)     ~ dataset_clean_Ne_Enet$mean_ratioEnet_wbeta*dataset_clean_Ne_Enet$Scenario),
                           lm(log(dataset_clean_Ne_Enet$`Nₑ`)    ~ dataset_clean_Ne_Enet$mean_ratioEnet_wbeta*dataset_clean_Ne_Enet$Scenario),
                           lm(log(dataset_clean_Ne_Enet$`Vk.`)   ~ dataset_clean_Ne_Enet$mean_ratioEnet_wbeta*dataset_clean_Ne_Enet$Scenario),
                           lm(log(dataset_clean_Ne_Enet$`NT`)    ~ dataset_clean_Ne_Enet$mean_ratioEnet_wbeta*dataset_clean_Ne_Enet$Scenario),
                           lm(log(dataset_clean_Ne_Enet$`Nₑ/NT`) ~ dataset_clean_Ne_Enet$mean_ratioEnet_wbeta*dataset_clean_Ne_Enet$Scenario),
                           lm(log(dataset_clean_Ne_Enet$`T/Vk.`) ~ dataset_clean_Ne_Enet$mean_ratioEnet_wbeta*dataset_clean_Ne_Enet$Scenario),
                           lm(log(dataset_clean_Ne_Enet$`Nα`)    ~ dataset_clean_Ne_Enet$mean_ratioEnet_wbeta*dataset_clean_Ne_Enet$Scenario)),
                      AIC)
  ) %>%
    mutate(prefer = ifelse(aic_log < aic_raw, "log", "raw"),
           enet_dimension = enet_dimension)
  # AIC is only comparable within the same response scale (raw AIC vs raw
  # AIC, log AIC vs log AIC) if the response transform is a strictly
  # monotonic reparametrization -- for log-normal-shaped data this
  # comparison is standard practice, but confirm with a residual QQ plot
  # per model if you want to be rigorous rather than relying on AIC alone.
  print(model_choice)
  
  # PLOT summarize 
  p1 <- ggplot(dataset_clean_Ne_Enet) +
    aes(y = log(`Nₑ/NT`), size = `Vk.`, x = log(T), colour=mean_ratioEnet_wbeta) +
    geom_point(size=2) +
    geom_text_repel(max.overlaps=15, aes(label=Species), size=4)+
    scale_color_continuous(palette = "heat", trans  = 'reverse')+
    # scale_color_gradient() +
    labs(
      y = "log Nₑ/NT",
      size = "Vk.",
      color = "ratio Enet",
      x = "log Generation length"
    ) +
    theme_bw()+
    facet_wrap(~Scenario)
  p2 <- ggplot(dataset_clean_Ne_Enet) +
    aes(y = log(`Nₑ/NT`), size =T, x = log(`Vk.`), colour=mean_ratioEnet_wbeta) +
    geom_point(size=2) +
    geom_text_repel(max.overlaps=15, aes(label=Species), size=4)+
    scale_color_continuous(palette = "heat", trans  = 'reverse')+
    # scale_color_gradient() +
    labs(
      y = "log Nₑ/NT",
      size = "T",
      color = "ratio Enet",
      x = "log Vk."
    ) +
    theme_bw()+
    facet_wrap(~Scenario)
  p3 <- ggplot(dataset_clean_Ne_Enet) +
    aes(y = log(`Nₑ/NT`),x = log(Nα), colour=mean_ratioEnet_wbeta) +
    geom_point(size=2) +
    geom_text_repel(max.overlaps=15, aes(label=Species), size=4)+
    scale_color_continuous(palette = "heat", trans  = 'reverse')+
    # scale_color_gradient() +
    labs(
      y = "log Nₑ/NT",
      color = "ratio Enet",
      x = "log Nα", 
      title = enet_dimension
    ) +
    theme_bw()+
    facet_wrap(~Scenario)
  
  ggsave(p1, file=paste0(PATHOUT, "/", enet_dimension,"Temp_T.png"), 
         height = 4.29, width=10.3)
  ggsave(p2, file=paste0(PATHOUT, "/", enet_dimension,"Temp_Vk.png"), 
         height = 4.29, width=10.3)
  ggsave(p3, file=paste0(PATHOUT, "/", enet_dimension,"Temp_Nalpha.png"), 
         height = 4.29, width=10.3)
  

  ### REGRESSION DATA_NE ~ DATA_ENET*Scenario   ------
   ### focus on the DIFF between CC and No CC scenarios: using the delta method
  
  df_diff_Enet_scrio <- dataset_clean_Ne_Enet
  df_diff_Enet_scrio <- pivot_longer(df_diff_Enet_scrio, names_to = "variable", values_to = "value", cols = -c("Species", "Scenario"))
  df_diff_Enet_scrio <- pivot_wider(df_diff_Enet_scrio, names_from = "Scenario", values_from = "value")
  df_diff_Enet_scrio_sum <- df_diff_Enet_scrio %>%
    group_by(Species, variable) %>%                       # <- adapt to your grouping columns
    summarise(
      mu_x  = mean(`With CC`, na.rm = TRUE),
      var_x = var(`With CC`,  na.rm = TRUE),
      n_x   = sum(!is.na(`With CC`)),
      mu_y  = mean(`No CC`,   na.rm = TRUE),
      var_y = var(`No CC`,    na.rm = TRUE),
      n_y   = sum(!is.na(`No CC`)),
      .groups = "drop"
    ) %>%
    mutate(
      # Expectation (2nd-order Taylor), Cov = 0 because scenarios are independent
      diff    = sign(mu_y) * (mu_x / mu_y - 1 + mu_x * var_y / mu_y^3),
      # Between-replicate variability (delta method, consistent with your SDs)
      sd_diff = sqrt(var_x / mu_y^2 + mu_x^2 * var_y / mu_y^4),
      # Uncertainty on the mean response (use this instead if you want a SE / CI)
      se_diff = sqrt((var_x / n_x) / mu_y^2 + mu_x^2 * (var_y / n_y) / mu_y^4)
    )
  
  df_diff_Enet_scrio$diff <- (df_diff_Enet_scrio$`With CC` - df_diff_Enet_scrio$`No CC`)/abs(df_diff_Enet_scrio$`No CC`)
  dataset_clean_Ne_Enet_diff <- pivot_wider(df_diff_Enet_scrio[, -which(colnames(df_diff_Enet_scrio) %in% c("No CC", "With CC"))], names_from="variable", values_from="diff")
  
  
  # test data before regression
  ggplotdata <- df_diff_Enet_scrio
  # ggplotdata$value <- log(ggplotdata$value)
  ggplotdata$id <- 1:nrow(df_diff_Enet_scrio) 
  ggplot(df_diff_Enet_scrio)+ 
    geom_density(aes(diff))+ 
    facet_wrap(~variable, scales="free")
  ggplot(df_diff_Enet_scrio, aes(sample=diff))+ 
    stat_qq()+ 
    stat_qq_line()+ 
    facet_wrap(~variable, scales="free")
  
  
  # Created nested tibble
  formula = c("`Nₑ/NT` ~ mean_ratioEnet_wbeta",
              "`NT` ~ mean_ratioEnet_wbeta",
              "`Nₑ` ~ mean_ratioEnet_wbeta", 
              "`T` ~ mean_ratioEnet_wbeta", 
              "`Vk.` ~ mean_ratioEnet_wbeta",
              "`T/Vk.` ~ mean_ratioEnet_wbeta",
              "`Nα` ~ mean_ratioEnet_wbeta")
  
  nested_df <- tibble(formula = formula) %>%
    group_by(ID = formula) %>%
    group_modify(~ as_tibble(dataset_clean_Ne_Enet_diff)) %>%
    nest() 
  # Get model estimates
  nested_df %>%
    mutate(estimates = data %>% map2(ID, ~ tidy(lm(.y, data = .x)))) %>%
    select(-data) %>%
    unnest()
  # Get fitted values and residuals
  nested_df %>%
    mutate(model = ID %>% map2(data, lm),
           stats = model %>% map(augment)) %>%
    select(-data, -model) %>%
    unnest() 
  # View(nested_df)
  star <- nested_df %>%
    mutate(model = ID %>% map2(data, ~ lm(.x, .y))) %>%
    pull(model) %>%
    stargazer(type = "text", digits = 3, header = F, digit.separator = '', align = T)
  replace_numbers = function(x, cutoff=4, digits=3, scipen=-7) {
    ifelse(nchar(x) < cutoff, x, prettyNum(as.numeric(x), digits=digits, scientific=scipen))
  }
  tablediff[[enet_dimension]] <- gsubfn("([0-9.]+)", ~replace_numbers(x), star)
  
  
  
 #### LEAVE-ONE-OUT SENSITIVITY ANALYSIS ------------------------------------
  ## For each response variable regressed on mean_ratioEnet_wbeta (relative
  ## diff in net energy), refit the model after excluding each species in
  ## turn. Flags cases where a single species drives significance, R2, or
  ## even the sign of the relationship.
  
  response_vars <- c("Nₑ/NT", "NT", "Nₑ", "T", "Vk.", "T/Vk.", "Nα")
  
  loo_results <- list()
  for (resp in response_vars) {
    
    form <- as.formula(paste0("`", resp, "` ~ mean_ratioEnet_wbeta"))
    
    # full model (all 15 species)
    full_model   <- lm(form, data = dataset_clean_Ne_Enet_diff)
    full_summary <- summary(full_model)
    full_p       <- coef(full_summary)[2, 4]
    full_r2      <- full_summary$r.squared
    full_slope   <- coef(full_model)[2]
    
    # refit excluding each species, one at a time
    loo_df <- purrr::map_dfr(dataset_clean_Ne_Enet_diff$Species, function(sp) {
      data_loo   <- dataset_clean_Ne_Enet_diff[dataset_clean_Ne_Enet_diff$Species != sp, ]
      model_loo  <- lm(form, data = data_loo)
      s          <- summary(model_loo)
      tibble(
        enet_dimension = enet_dimension,
        variable       = resp,
        excluded       = sp,
        slope          = coef(model_loo)[2],
        r2             = s$r.squared,
        p_value        = coef(s)[2, 4],
        full_slope     = full_slope,
        full_r2        = full_r2,
        full_p         = full_p,
        sign_flip      = sign(coef(model_loo)[2]) != sign(full_slope),
        loses_signif   = full_p < 0.05 & coef(s)[2, 4] >= 0.05,
        gains_signif   = full_p >= 0.05 & coef(s)[2, 4] < 0.05
      )
    })
    
    loo_results[[resp]] <- loo_df
  }
  
  loo_all <- dplyr::bind_rows(loo_results)
  # save full leave-one-out table for this life stage
  write.csv(loo_all,
            file = paste0(PATHOUT, "/", enet_dimension, "_leave_one_out.csv"),
            row.names = FALSE)
  
  # print only the flagged cases (significance flip, sign flip, or large R2 swing)
  loo_flagged <- loo_all %>%
    dplyr::filter(loses_signif | gains_signif | sign_flip)
  
  if (nrow(loo_flagged) > 0) {
    cat("\n>>> LEAVE-ONE-OUT FLAGS —", enet_dimension, ":\n")
    print(loo_flagged[, c("variable", "excluded", "full_p", "p_value",
                          "loses_signif", "gains_signif", "sign_flip")])
  } else {
    cat("\n>>> No leave-one-out flags for", enet_dimension, "— all relationships stable.\n")
  }
  
 #### Lollipop - leave one out
  # visual: slope estimate per exclusion, one panel per response variable,
  # dashed line = full-sample slope, red = flagged exclusion
  loo_all <- loo_all %>%
    mutate(flag_type = case_when(
      sign_flip & loses_signif         ~ "Sign flip + loses significance",
      sign_flip & gains_signif         ~ "Sign flip + gains significance",
      sign_flip                        ~ "Sign flip",
      loses_signif                     ~ "Loses significance",
      gains_signif                     ~ "Gains significance",
      TRUE                             ~ "Stable"
    ))
  
  p_loo <- ggplot(loo_all, aes(x = excluded, y = slope)) +
    geom_hline(aes(yintercept = full_slope), color = "grey50", linetype = "dashed") +
    geom_point(aes(color = flag_type), size = 2) +
    scale_color_manual(values = c(
      "Stable"                          = "grey30",
      "Loses significance"              = "brown2",
      "Gains significance"              = "steelblue3",
      "Sign flip"                       = "darkorange2",
      "Sign flip + loses significance"  = "firebrick4",
      "Sign flip + gains significance"  = "purple4"
    )) +
    coord_flip() +
    facet_wrap(~variable, scales = "free_x", ncol = 3) +
    labs(x = "Species excluded", y = "Slope (leave-one-out)",
         color = "Flag",
         title = paste0("Leave-one-out sensitivity — ", enet_dimension)) +
    theme_bw()
  
  ggsave(p_loo,
         file = paste0(PATHOUT, "/", enet_dimension, "_leaveoneout_plot.png"),
         height = 9, width = 11)
  
  
  
  
 #### Lollipop - Enet
  # no regressions with log because diff => values <0
  df_diff_Enet_lolli <- df_diff_Enet_scrio[which(df_diff_Enet_scrio$variable %in% 
                                                   c("mean_ratioEnet_wbeta")),] #, "NT", "Nₑ", "Nₑ/NT", "T", "Vk.", "T/Vk.", "Nα")),]
  df_diff_Enet_lolli$variable <- str_replace(df_diff_Enet_lolli$variable, pattern = "mean_ratioEnet_wbeta", enet_dimension)
  # df_diff_Enet_lolli$variable <- factor(df_diff_Enet_lolli$variable, levels = c("Net energy", "NT", "Nₑ", "Nₑ/NT", "T", "Vk.", "T/Vk.", "Nα"))
  lollipop <- ggplot(df_diff_Enet_lolli, aes(x=Species, y=diff))  + #[which(df_diff_Enet_lolli$variable == "Net energy"),]
    geom_hline(yintercept = 0, color = "darkgrey")+
    geom_segment(aes(x=Species, xend=Species, y=0, yend=diff), color="grey") +
    geom_point(color=ifelse(round(df_diff_Enet_lolli$diff, digits=2)<0, #[which(df_diff_Enet_lolli$variable == "Net energy"),]
                            "brown2", "palegreen3"), size=4) +
    coord_flip() + 
    theme_bw() +
    theme(
      panel.grid.major.x = element_blank(),
      panel.border = element_blank(),
      axis.ticks.x = element_blank()
    ) +
    xlab("") +
    ylab("")+
    # labs(title = enet_dimension) +
    facet_wrap(~variable, scale="free_x", ncol=3)+
    scale_y_continuous() #n.breaks = 8
  ggsave(lollipop, file=paste0(PATHOUT, "/", enet_dimension,"lollipop.png"),
         height = 5, width=4)
  
  
 #### Geom_line - Enet
  theme_lollipop <- theme_minimal(base_size = 14) +
    theme(
      text = element_text(color = "grey20"),
      plot.title = element_text(size = 15, face = "plain", hjust = 0),
      axis.title = element_text(size = 12, color = "grey30"),
      axis.text = element_text(size = 11, color = "grey30"),
      panel.grid.major = element_line(color = "grey90", linewidth = 0.4),
      panel.grid.minor = element_blank(),
      axis.line = element_blank(),
      legend.position = "none",
      strip.background = element_rect(fill = "grey90", color = NA),
      strip.text = element_text(size = 12, color = "grey20"),
      panel.border = element_rect(colour = "black", fill=NA, linewidth=1)
    )
  # PLOT summarize 
  df_diff_Enet_scrio_w <- pivot_wider(df_diff_Enet_scrio[, -which(colnames(df_diff_Enet_scrio) %in% c("No CC", "With CC"))], 
              values_from = "diff", names_from = "variable")
  p1 <- ggplot(df_diff_Enet_scrio_w) +
    aes(y = `Nₑ/NT_sd`, size = `Vk.`, x = T, colour=log(mean_ratioEnet_wbeta)) +
    geom_point(size=2) +
    geom_text_repel(max.overlaps=15, aes(label=Species), size=4)+
    scale_color_continuous(palette = "heat", trans  = 'reverse')+
    # scale_color_gradient() +
    labs(
      y = "Nₑ/NT",
      size = "Vk.",
      color = "ratio Enet",
      x = "Generation length",
      title = "Comparing variables with and without climate change" 
    ) +
    theme_bw()
  p2 <- ggplot(df_diff_Enet_scrio_w) +
    aes(y = (`Nₑ/NT`), size =T, x = (`Vk.`), colour=log(mean_ratioEnet_wbeta)) +
    geom_point(size=2) +
    geom_text_repel(max.overlaps=15, aes(label=Species), size=4)+
    scale_color_continuous(palette = "heat", trans  = 'reverse')+
    # scale_color_gradient() +
    labs(
      y = "Nₑ/NT",
      size = "T",
      color = "ratio Enet",
      x = "Vk.",
      title = "Comparing variables with and without climate change"
    ) +
    theme_bw()
  p3 <- ggplot(df_diff_Enet_scrio_w) +
    aes(y = `Nₑ/NT`,x = Nα, colour=mean_ratioEnet_wbeta) +
    geom_point(size=2) +
    geom_text_repel(max.overlaps=15, aes(label=Species), size=4)+
    scale_color_continuous(palette = "heat", trans  = 'reverse')+
    # scale_color_gradient() +
    labs(
      y = "Nₑ/NT",
      color = "ratio Enet",
      x = "Nα",
      title = "Comparing variables with and without climate change"
    ) +
    theme_bw()
  p4 <- ggplot(df_diff_Enet_scrio_w) +
    aes(y = (`Nₑ/NT`), x = mean_ratioEnet_wbeta) +   # <- x et y inversés ici
    geom_smooth(color = "#e41a1c", fill = "#69b3a2", alpha = 0.3, se = TRUE, method = "lm") +
    geom_point(size = 3, color = "grey15") +
    geom_text_repel(max.overlaps=15, aes(label = Species), size = 4, color = "grey20") +
    labs(
      y = "Nₑ/NT (w CC - no CC) / (no CC) ",
      x = "Net energy (w CC - no CC) / (no CC) "
    ) +
    theme_lollipop
  p5 <- ggplot(df_diff_Enet_scrio_w) +
    aes(y = (`Vk.`), x = mean_ratioEnet_wbeta) + 
    geom_smooth(color = "#e41a1c", fill = "#69b3a2", alpha = 0.3, se = TRUE, method = "lm") +
    geom_point(size = 3, color = "grey15") +
    geom_text_repel(max.overlaps=15, aes(label = Species), size = 4, color = "grey20") +
    labs(
      y = "Vk. (w CC - no CC) / (no CC) ",
      x = "Net energy (w CC - no CC) / (no CC) "
    ) +
    theme_lollipop
  p6 <-   ggplot(df_diff_Enet_scrio_w) +
    aes(y = `T`, x = mean_ratioEnet_wbeta) + 
    geom_smooth(color = "#e41a1c", fill = "#69b3a2", alpha = 0.3, se = TRUE, method = "lm") +
    geom_point(size = 3, color = "grey15") +
    geom_text_repel(max.overlaps=15, aes(label = Species), size = 4, color = "grey20") +
    labs(
      y = "T (w CC - no CC) / (no CC) ",
      x = "Net energy (w CC - no CC) / (no CC) "
    ) +
    theme_lollipop
  p7 <- ggplot(df_diff_Enet_scrio_w) +
    aes(y = `Nα`, x = mean_ratioEnet_wbeta) + 
    geom_smooth(color = "#e41a1c", fill = "#69b3a2", alpha = 0.3, se = TRUE, method = "lm") +
    geom_point(size = 3, color = "grey15") +
    geom_text_repel(max.overlaps=15, aes(label = Species), size = 4, color = "grey20") +
    labs(
      y = "Nα (w CC - no CC) / (no CC) ",
      x = "Net energy (w CC - no CC) / (no CC) "
    ) +
    theme_lollipop
  p8 <-ggplot(df_diff_Enet_scrio_w) +
    aes(y = `Nₑ`, x = mean_ratioEnet_wbeta) +
  geom_smooth(color = "#e41a1c", fill = "#69b3a2", alpha = 0.3, se = TRUE, method = "lm") +
  geom_point(size = 3, color = "grey15") +
  geom_text_repel(max.overlaps=15, aes(label = Species), size = 4, color = "grey20") +
  labs(
    y = "Nₑ (w CC - no CC) / (no CC) ",
    x = "Net energy (w CC - no CC) / (no CC) ")+
  theme_lollipop
  p9 <- ggplot(df_diff_Enet_scrio_w) +
    aes(y = `NT`, x = mean_ratioEnet_wbeta) + 
  geom_smooth(color = "#e41a1c", fill = "#69b3a2", alpha = 0.3, se = TRUE, method = "lm") +
  geom_point(size = 3, color = "grey15") +
  geom_text_repel(max.overlaps=15, aes(label = Species), size = 4, color = "grey20") +
  labs(
    y = "NT (w CC - no CC) / (no CC) ",
    x = "Net energy (w CC - no CC) / (no CC) "
  ) +
  theme_lollipop
  p10 <- ggplot(df_diff_Enet_scrio_w) +
    aes(y = `T/Vk.`, x = mean_ratioEnet_wbeta) +  
  geom_smooth(color = "#e41a1c", fill = "#69b3a2", alpha = 0.3, se = TRUE, method = "lm") +
  geom_point(size = 3, color = "grey15") +
  geom_text_repel(max.overlaps=15, aes(label = Species), size = 4, color = "grey20") +
  labs(
    y = "T/Vk. (w CC - no CC) / (no CC) ",
    x = "Net energy (w CC - no CC) / (no CC) ",
  ) +
  theme_lollipop
  
  compo <- ggarrange(p9, p8, p4, p6, p5, p10, p7, ncol=3, nrow=3)
  ggsave(compo, file=paste0(PATHOUT, "/", enet_dimension,"Temp_Nalpha.png"),
         height = 13, width=13)
  

  compo <- ggarrange(p6,p5,p10,p4, ncol=2, nrow=2)
  ggsave(compo, file=paste0(PATHOUT, "/", enet_dimension,"Temp_Nalpha_2.png"),
         height = 8, width=8)

}
# table
# tablelog
tablediff # the one used in the supp mat
# tabledifflog



