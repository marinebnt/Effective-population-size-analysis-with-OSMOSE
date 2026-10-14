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
  
  data_Enet_larv[[sp]]  <- data_Enet[which(data_Enet$Species == sp & data_Enet$age == 0), ]%>% # => based on get_par(conf, par = "species.larvae.growth.threshold")
    group_by(Species, Scenario, variable, Fmax, Fstart, replicate) %>%
    summarise(value  = mean(value, na.rm=T)) %>%
    ungroup
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

data_Enet_mage_w  <- pivot_wider(data_Enet_mage, names_from = "variable", values_from = c("mean"))
data_Enet_mage_w  <- data_Enet_mage_w[which(data_Enet_mage_w$Fmax == 1), -which(colnames(data_Enet_mage_w) %in%c("Fmax", "Fstart"))]
data_Enet_larv_w  <- pivot_wider(data_Enet_larv, names_from = "variable", values_from = c("value"))
data_Enet_larv_w  <- data_Enet_larv_w[which(data_Enet_larv_w$Fmax == 1), -which(colnames(data_Enet_larv_w) %in%c("Fmax", "Fstart"))]
data_Enet_juv_w   <- pivot_wider(data_Enet_juv, names_from = "variable", values_from = c("value"))
data_Enet_juv_w   <- data_Enet_juv_w[which(data_Enet_juv_w$Fmax == 1), -which(colnames(data_Enet_juv_w) %in%c("Fmax", "Fstart"))]
data_Enet_adult_w <- pivot_wider(data_Enet_adult, names_from = "variable", values_from = c("value"))
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
  PATHIN  = paste0("output/PLOTS/", toupper(Scenario), "/KEY_TRAITS_MATURE_removed0inNedata/")
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
    final_dataset_age <- when_variables_missing_age_dim(dataset = final_dataset_age, add_missing_ages = F)
    final_dataset_age <- what_ages_have_null_dead_nb(final_dataset_age)
    final_dataset_age <- when_null_or_low_abund(final_dataset_age)

    
    ## Extract only the dimensions I will be using
    final_dataset_age <- final_dataset_age[which(final_dataset_age$Fmax == 1), ]
    NeNT2 <- final_dataset_age[which(final_dataset_age$variable %in% c("Ne_N_alpha", "NT", "generation_length", 
                                                                       "Vk_point", "Numerator_N_alpha", "survivorship_curve")),] 
    
    data_t <- estimatecolSum(NeNT2, var = paste0("colSum", "survivorship_curve"))
    data_t$variable <- c("L\u03C4")
    NeNT2 <- rbind(NeNT2[-which(NeNT2$variable == "survivorship_curve"), -which(colnames(NeNT2) == c("age"))], data_t)
    
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
Ne_NT_tot2 <- Ne_NT_tot2[which(Ne_NT_tot2$Fmax == 1), ]
Ne_NT_tot2 <- Ne_NT_tot2[, -which(colnames(Ne_NT_tot2) %in% c("Fmax", "age"))]
Ne_NT_tot2$Ne_NT  <- Ne_NT_tot2$Ne_N_alpha/Ne_NT_tot2$NT 
Ne_NT_tot2$'T/Vk' <- Ne_NT_tot2$generation_length/Ne_NT_tot2$Vk_point 

# data for Fmax == 1 & cleaning
Ne_NT_bx   <- pivot_longer(Ne_NT_tot2, cols=-c("replicate", "Species", "Scenario"),
                           values_to = "value", names_to = "variable")
dataset_clean_Ne <- Ne_NT_bx[which(Ne_NT_bx$variable %in% c("Ne_N_alpha", "Ne_NT", "NT", 
                                                            "Numerator_N_alpha", "Vk_point",
                                                            "generation_length",
                                                            "T/Vk", "L\u03C4")), ]


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


## COMPARE ENET effect on variables depending on life stage -------------------
PATHOUT  = paste0("output/PLOTS/CHAP3/KEY_TRAITS_MATURE_removed0inNedata/Enet/")
dir.create(PATHOUT)
list()-> tablediff
list()-> tabledifflog
list()-> table
list()-> tablelog
dataset_clean_Ne <- pivot_wider(dataset_clean_Ne, names_from = "variable", values_from = "value")# summarize 

for (enet_dimension in names(list_Enet_data)){
  
  ## add the Enet data  ----------------
  data                  <- list_Enet_data[[enet_dimension]]
  dataset_clean_Ne$replicate <- as.numeric(dataset_clean_Ne$replicate)
  data$replicate <- as.numeric(data$replicate)+1
  dataset_clean_Ne_Enet <- full_join(data, dataset_clean_Ne)
  colnames(dataset_clean_Ne_Enet) <- improve_vriable_names(colnames(dataset_clean_Ne_Enet))
  if ("age" %in% colnames(dataset_clean_Ne_Enet)) dataset_clean_Ne_Enet <- dataset_clean_Ne_Enet[, -which(colnames(dataset_clean_Ne_Enet) == "age")]
  if ("Individual weight" %in% colnames(dataset_clean_Ne_Enet)) dataset_clean_Ne_Enet <- dataset_clean_Ne_Enet[, -which(colnames(dataset_clean_Ne_Enet) == "Individual weight")]
  dataset_clean_Ne_Enet$Species <- improve_sp_names(dataset_clean_Ne_Enet$Species)
  dataset_clean_Ne_Enet$Species <- factor(dataset_clean_Ne_Enet$Species, levels= c("Common dab", "Sandeel", "European sprat", "Atlantic mackerel", "Norway pout", 
                                                                             "Saithe", "Haddock", "Whiting", "Grey gurnard", "Atlantic cod", "Hake",
                                                                             "Horse mackerel", "Atlantic herring", "Common sole", "European plaice"))
  
  ### REGRESSION DATA_NE ~ DATA_ENET*Scenario   ------
  ### focus on the DIFF between CC and No CC scenarios: using the delta method
  df_diff_Enet_scrio <- pivot_longer(dataset_clean_Ne_Enet, names_to = "variable", values_to = "value", cols = -c("Species", "Scenario", "replicate"))
  df_diff_Enet_scrio <- pivot_wider(df_diff_Enet_scrio, names_from = "Scenario", values_from = "value")
  # estimate difference and its uncertainty by using the delta method
  df_diff_Enet_scrio <- df_diff_Enet_scrio %>%
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
      se_diff = sqrt((var_x / n_x) / mu_y^2 + mu_x^2 * (var_y / n_y) / mu_y^4),
      # half-width of the IC95 (Student’s t-test, conservative lower limit)
      IC_95   = qt(0.975, df = pmin(n_x, n_y) - 1) * se_diff
    )
  dataset_clean_Ne_Enet_diff <- pivot_wider(df_diff_Enet_scrio[, c("Species", "variable", "diff", "se_diff")], names_from="variable", values_from=c("diff", "se_diff"))
  colnames(dataset_clean_Ne_Enet_diff) <- stringr::str_remove(colnames(dataset_clean_Ne_Enet_diff), "^diff_")
  colnames(dataset_clean_Ne_Enet_diff) <- stringr::str_replace(colnames(dataset_clean_Ne_Enet_diff), pattern = "se_diff_", replacement = "se_")
  
  # test data before regression
  ggplotdata <- dataset_clean_Ne_Enet_diff
  # ggplotdata$value <- log(ggplotdata$value)
  ggplotdata$id <- 1:nrow(dataset_clean_Ne_Enet_diff) 
  ggplot(dataset_clean_Ne_Enet_diff)+ 
    geom_density(aes(mean))+ 
    facet_wrap(~variable, scales="free")
  ggplot(dataset_clean_Ne_Enet_diff, aes(sample=mean))+ 
    stat_qq()+ 
    stat_qq_line()+ 
    facet_wrap(~variable, scales="free")
  
  
  # Created nested tibble
  formula = c("`Nₑ/NT` ~ ratioEnet_wbeta",
              "NT ~ ratioEnet_wbeta",
              "Nₑ ~ ratioEnet_wbeta", 
              "`T/Vk.` ~ ratioEnet_wbeta",
              "Lτ ~ ratioEnet_wbeta",
              "T ~ ratioEnet_wbeta", 
              "Vk. ~ ratioEnet_wbeta")
  
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
  models_diff <- nested_df %>%
    mutate(model = ID %>% map2(data, ~ lm(.x, .y)))
  p_lm <- setNames(sapply(models_diff$model, function(m) coef(summary(m))[2, 4]),
                   gsub("`| ~ ratioEnet_wbeta", "", models_diff$ID))   # noms : "Nₑ/NT", "NT", "Lτ", ...
  p_label <- function(v) sprintf("italic(p)~'%s'",
                                 ifelse(p_lm[[v]] < 0.001, "< 0.001", paste("=", sprintf("%.3f", p_lm[[v]]))))
  
  
  #### LEAVE-ONE-OUT SENSITIVITY ANALYSIS ------------------------------------
  ## For each response variable regressed on ratioEnet_wbeta (relative
  ## diff in net energy), refit the model after excluding each species in
  ## turn. Flags cases where a single species drives significance, R2, or
  ## even the sign of the relationship.
  
  response_vars <- c("Nₑ/NT", "NT", "Nₑ", "T/Vk.", "L\u03C4", "T", "Vk.")
  
  loo_results <- list()
  for (resp in response_vars) {
    
    form <- as.formula(paste0("`", resp, "` ~ ratioEnet_wbeta"))
    
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
  loo_all$variable <- factor(loo_all$variable, levels=c("Nₑ", "NT", "Nₑ/NT", "T/Vk.", "Lτ", "T", "Vk."))
  p_loo <- ggplot(loo_all, aes(x = excluded, y = slope)) +
    geom_hline(aes(yintercept = full_slope), color = "grey50", linetype = "dashed") +
    geom_point(aes(color = flag_type), size = 2.5) +
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
    labs(x = "Species", y = "Slope (leave-one-out)",
         color = "Flag",
         title = paste0("Leave-one-out sensitivity — ", enet_dimension)) +
    theme_bw()+
    theme(text=element_text(size=12),
          legend.position = "inside",
          legend.position.inside = c(0.83, 0.17),
          legend.justification = c(0.5, 0.5),
          legend.background = element_rect(colour = "grey70", fill = "white"))
  
  ggsave(p_loo,
         file = paste0(PATHOUT, enet_dimension, "_leaveoneout_plot.png"),
         height = 10, width = 8)
  
  
  
  
  #### Lollipop - Enet
  # no regressions with log because diff => values <0
  df_diff_Enet_lolli <- df_diff_Enet_scrio[which(df_diff_Enet_scrio$variable %in% 
                                                   c("ratioEnet_wbeta")),] #, "NT", "Nₑ", "Nₑ/NT", "T", "Vk.", "T/Vk.", "Nα")),]
  df_diff_Enet_lolli$variable <- str_replace(df_diff_Enet_lolli$variable, pattern = "ratioEnet_wbeta", enet_dimension)
  # df_diff_Enet_lolli$variable <- factor(df_diff_Enet_lolli$variable, levels = c("Net energy", "NT", "Nₑ", "Nₑ/NT", "T", "Vk.", "T/Vk.", "Nα"))
  lollipop <- ggplot(df_diff_Enet_lolli, aes(x=Species, y=diff))  + #[which(df_diff_Enet_lolli$variable == "Net energy"),]
    geom_hline(yintercept = 0, color = "darkgrey")+
    geom_segment(aes(x=Species, xend=Species, y=0, yend=diff), color="grey") +
    geom_errorbar(aes(ymin = diff - IC_95, ymax = diff + IC_95), width = 0.3, color = "grey30") +
    geom_point(color = ifelse(df_diff_Enet_lolli$diff - df_diff_Enet_lolli$IC_95 <= 0 &
                                df_diff_Enet_lolli$diff + df_diff_Enet_lolli$IC_95 >= 0, "grey60",
                              ifelse(df_diff_Enet_lolli$diff < 0, "brown2", "palegreen3")), size = 4) +
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
  # [MODIFIÉ] thème allégé : cadre plus fin, textes un peu plus petits
  theme_lollipop <- theme_minimal(base_size = 12) +
    theme(
      text = element_text(color = "grey20"),
      plot.title = element_text(size = 13, face = "plain", hjust = 0),
      axis.title = element_text(size = 11, color = "grey20"),
      axis.text = element_text(size = 10, color = "grey30"),
      panel.grid.major = element_line(color = "grey93", linewidth = 0.3),
      panel.grid.minor = element_blank(),
      axis.line = element_blank(),
      legend.position = "none",
      strip.background = element_rect(fill = "grey90", color = NA),
      strip.text = element_text(size = 12, color = "grey20"),
      panel.border = element_rect(colour = "grey20", fill = NA, linewidth = 0.6),
      plot.margin = margin(8, 10, 6, 6),
      plot.background = element_rect(fill = "white", colour = NA)   # [AJOUT] fond blanc (sinon transparent)
    )
  # [AJOUT] noms d'espèces lisibles : "AtlanticCod" -> "Atlantic cod"
  dataset_clean_Ne_Enet_diff$Species_lab <- gsub("([a-z])([A-Z])", "\\1 \\L\\2",
                                                 dataset_clean_Ne_Enet_diff$Species, perl = TRUE)
  # PLOT summarize 
  pblank <- ggplot() + theme_void()
  p1 <- ggplot(dataset_clean_Ne_Enet_diff) +
    aes(y = `Nₑ/NT`, size = `Vk.`, x = T, colour=log(ratioEnet_wbeta)) +
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
  p2 <- ggplot(dataset_clean_Ne_Enet_diff) +
    aes(y = (`Nₑ/NT`), size =T, x = (`Vk.`), colour=log(ratioEnet_wbeta)) +
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
  p3 <- ggplot(dataset_clean_Ne_Enet_diff) +
    aes(y = `Nₑ/NT`,x = Nα, colour=ratioEnet_wbeta) +
    geom_point(size=2) +
    geom_text_repel(max.overlaps=15, aes(label=Species), size=4)+
    scale_color_continuous(palette = "heat", trans  = 'reverse')+
    # scale_color_gradient() +
    labs(
      y = "Relative change in Nₑ/NT",
      color = "ratio Enet",
      x = "Relative change in Nα",
      title = "Comparing variables with and without climate change"
    ) +
    theme_bw()
  p4 <- ggplot(dataset_clean_Ne_Enet_diff) +
    aes(y = `Nₑ/NT`, x = ratioEnet_wbeta) +   # <- x et y inversés ici
    # [AJOUT] repères 0 : signe des deux variations lisible directement
    geom_hline(yintercept = 0, colour = "grey60", linetype = "dashed", linewidth = 0.4) +
    geom_vline(xintercept = 0, colour = "grey60", linetype = "dashed", linewidth = 0.4) +
    # [MODIFIÉ] ruban plus léger, droite plus fine
    geom_smooth(color = "#e41a1c", fill = "#69b3a2", alpha = 0.2, se = TRUE, method = "lm",
                linewidth = 0.9, formula = y ~ x) +
    # [MODIFIÉ] points à contour blanc (se détachent du ruban)
    geom_point(size = 2.6, shape = 21, fill = "grey15", colour = "white", stroke = 0.5) +
    # [MODIFIÉ] étiquettes plus petites, traits discrets, placement reproductible (seed)
    geom_text_repel(aes(label = Species_lab), size = 3.1, color = "grey25",
                    max.overlaps = Inf, box.padding = 0.4, point.padding = 0.25,
                    min.segment.length = 0, segment.color = "grey65", segment.size = 0.3, seed = 1) +
    # [AJOUT] p de la régression, coin inférieur droit
    annotate("text", x = Inf, y = -Inf, hjust = 1.05, vjust = -0.6, size = 3.3, colour = "grey30", parse = TRUE,
             label = p_label("Nₑ/NT")) +
    # [AJOUT] axes en % (plus de notation 3e-04) ; marge basse pour p
    #scale_x_continuous(labels = scales::label_percent(accuracy = 1)) +
    scale_y_continuous(#labels = scales::label_percent(accuracy = 0.1),
                       expand = expansion(mult = c(0.12, 0.05))) +
    labs(
      y = "Relative change in Nₑ/NT",
      x = "Relative change in net energy"
    ) +
    theme_lollipop
  p5 <- ggplot(dataset_clean_Ne_Enet_diff) +
    aes(y = `Vk.`, x = ratioEnet_wbeta) +
    # [AJOUT] repères 0 : signe des deux variations lisible directement
    geom_hline(yintercept = 0, colour = "grey60", linetype = "dashed", linewidth = 0.4) +
    geom_vline(xintercept = 0, colour = "grey60", linetype = "dashed", linewidth = 0.4) +
    # [MODIFIÉ] ruban plus léger, droite plus fine
    geom_smooth(color = "#e41a1c", fill = "#69b3a2", alpha = 0.2, se = TRUE, method = "lm",
                linewidth = 0.9, formula = y ~ x) +
    # [MODIFIÉ] points à contour blanc (se détachent du ruban)
    geom_point(size = 2.6, shape = 21, fill = "grey15", colour = "white", stroke = 0.5) +
    # [MODIFIÉ] étiquettes plus petites, traits discrets, placement reproductible (seed)
    geom_text_repel(aes(label = Species_lab), size = 3.1, color = "grey25",
                    max.overlaps = Inf, box.padding = 0.4, point.padding = 0.25,
                    min.segment.length = 0, segment.color = "grey65", segment.size = 0.3, seed = 1) +
    # [AJOUT] p de la régression, coin inférieur droit
    annotate("text", x = Inf, y = -Inf, hjust = 1.05, vjust = -0.6, size = 3.3, colour = "grey30", parse = TRUE,
             label = p_label("Vk.")) +
    # [AJOUT] axes en % (plus de notation 3e-04) ; marge basse pour p
    #scale_x_continuous(labels = scales::label_percent(accuracy = 1)) +
    scale_y_continuous(#labels = scales::label_percent(accuracy = 0.1),
                       expand = expansion(mult = c(0.12, 0.05))) +
    labs(
      y = "Relative change in Vk.",
      x = "Relative change in net energy"
    ) +
    theme_lollipop
  p6 <- ggplot(dataset_clean_Ne_Enet_diff) +
    aes(y = `T`, x = ratioEnet_wbeta) +
    # [AJOUT] repères 0 : signe des deux variations lisible directement
    geom_hline(yintercept = 0, colour = "grey60", linetype = "dashed", linewidth = 0.4) +
    geom_vline(xintercept = 0, colour = "grey60", linetype = "dashed", linewidth = 0.4) +
    # [MODIFIÉ] ruban plus léger, droite plus fine
    geom_smooth(color = "#e41a1c", fill = "#69b3a2", alpha = 0.2, se = TRUE, method = "lm",
                linewidth = 0.9, formula = y ~ x) +
    # [MODIFIÉ] points à contour blanc (se détachent du ruban)
    geom_point(size = 2.6, shape = 21, fill = "grey15", colour = "white", stroke = 0.5) +
    # [MODIFIÉ] étiquettes plus petites, traits discrets, placement reproductible (seed)
    geom_text_repel(aes(label = Species_lab), size = 3.1, color = "grey25",
                    max.overlaps = Inf, box.padding = 0.4, point.padding = 0.25,
                    min.segment.length = 0, segment.color = "grey65", segment.size = 0.3, seed = 1) +
    # [AJOUT] p de la régression, coin inférieur droit
    annotate("text", x = Inf, y = -Inf, hjust = 1.05, vjust = -0.6, size = 3.3, colour = "grey30", parse = TRUE,
             label = p_label("T")) +
    # [AJOUT] axes en % (plus de notation 3e-04) ; marge basse pour p
    #scale_x_continuous(labels = scales::label_percent(accuracy = 1)) +
    scale_y_continuous(#labels = scales::label_percent(accuracy = 0.1),
                       expand = expansion(mult = c(0.12, 0.05))) +
    labs(
      y = "Relative change in T",
      x = "Relative change in net energy"
    ) +
    theme_lollipop
  p7 <- ggplot(dataset_clean_Ne_Enet_diff) +
    aes(y = Lτ, x = ratioEnet_wbeta) +
    # [AJOUT] repères 0 : signe des deux variations lisible directement
    geom_hline(yintercept = 0, colour = "grey60", linetype = "dashed", linewidth = 0.4) +
    geom_vline(xintercept = 0, colour = "grey60", linetype = "dashed", linewidth = 0.4) +
    # [MODIFIÉ] ruban plus léger, droite plus fine
    geom_smooth(color = "#e41a1c", fill = "#69b3a2", alpha = 0.2, se = TRUE, method = "lm",
                linewidth = 0.9, formula = y ~ x) +
    # [MODIFIÉ] points à contour blanc (se détachent du ruban)
    geom_point(size = 2.6, shape = 21, fill = "grey15", colour = "white", stroke = 0.5) +
    # [MODIFIÉ] étiquettes plus petites, traits discrets, placement reproductible (seed)
    geom_text_repel(aes(label = Species_lab), size = 3.1, color = "grey25",
                    max.overlaps = Inf, box.padding = 0.4, point.padding = 0.25,
                    min.segment.length = 0, segment.color = "grey65", segment.size = 0.3, seed = 1) +
    # [AJOUT] p de la régression, coin inférieur droit
    annotate("text", x = Inf, y = -Inf, hjust = 1.05, vjust = -0.6, size = 3.3, colour = "grey30", parse = TRUE,
             label = p_label("Lτ")) +
    # [AJOUT] axes en % (plus de notation 3e-04) ; marge basse pour p
    #scale_x_continuous(labels = scales::label_percent(accuracy = 1)) +
    scale_y_continuous(#labels = scales::label_percent(accuracy = 0.1),
                       expand = expansion(mult = c(0.12, 0.05))) +
    labs(
      y = "Relative change in L\u03C4",
      x = "Relative change in net energy"
    ) +
    theme_lollipop
  p8 <- ggplot(dataset_clean_Ne_Enet_diff) +
    aes(y = `Nₑ`, x = ratioEnet_wbeta) +
    # [AJOUT] repères 0 : signe des deux variations lisible directement
    geom_hline(yintercept = 0, colour = "grey60", linetype = "dashed", linewidth = 0.4) +
    geom_vline(xintercept = 0, colour = "grey60", linetype = "dashed", linewidth = 0.4) +
    # [MODIFIÉ] ruban plus léger, droite plus fine
    geom_smooth(color = "#e41a1c", fill = "#69b3a2", alpha = 0.2, se = TRUE, method = "lm",
                linewidth = 0.9, formula = y ~ x) +
    # [MODIFIÉ] points à contour blanc (se détachent du ruban)
    geom_point(size = 2.6, shape = 21, fill = "grey15", colour = "white", stroke = 0.5) +
    # [MODIFIÉ] étiquettes plus petites, traits discrets, placement reproductible (seed)
    geom_text_repel(aes(label = Species_lab), size = 3.1, color = "grey25",
                    max.overlaps = Inf, box.padding = 0.4, point.padding = 0.25,
                    min.segment.length = 0, segment.color = "grey65", segment.size = 0.3, seed = 1) +
    # [AJOUT] p de la régression, coin inférieur droit
    annotate("text", x = Inf, y = -Inf, hjust = 1.05, vjust = -0.6, size = 3.3, colour = "grey30", parse = TRUE,
             label = p_label("Nₑ")) +
    # [AJOUT] axes en % (plus de notation 3e-04) ; marge basse pour p
    #scale_x_continuous(labels = scales::label_percent(accuracy = 1)) +
    scale_y_continuous(#labels = scales::label_percent(accuracy = 0.1),
                       expand = expansion(mult = c(0.12, 0.05))) +
    labs(
      y = "Relative change in Nₑ",
      x = "Relative change in net energy"
    ) +
    theme_lollipop
  p9 <- ggplot(dataset_clean_Ne_Enet_diff) +
    aes(y = `NT`, x = ratioEnet_wbeta) +
    # [AJOUT] repères 0 : signe des deux variations lisible directement
    geom_hline(yintercept = 0, colour = "grey60", linetype = "dashed", linewidth = 0.4) +
    geom_vline(xintercept = 0, colour = "grey60", linetype = "dashed", linewidth = 0.4) +
    # [MODIFIÉ] ruban plus léger, droite plus fine
    geom_smooth(color = "#e41a1c", fill = "#69b3a2", alpha = 0.2, se = TRUE, method = "lm",
                linewidth = 0.9, formula = y ~ x) +
    # [MODIFIÉ] points à contour blanc (se détachent du ruban)
    geom_point(size = 2.6, shape = 21, fill = "grey15", colour = "white", stroke = 0.5) +
    # [MODIFIÉ] étiquettes plus petites, traits discrets, placement reproductible (seed)
    geom_text_repel(aes(label = Species_lab), size = 3.1, color = "grey25",
                    max.overlaps = Inf, box.padding = 0.4, point.padding = 0.25,
                    min.segment.length = 0, segment.color = "grey65", segment.size = 0.3, seed = 1) +
    # [AJOUT] p de la régression, coin inférieur droit
    annotate("text", x = Inf, y = -Inf, hjust = 1.05, vjust = -0.6, size = 3.3, colour = "grey30", parse = TRUE,
             label = p_label("NT")) +
    # [AJOUT] axes en % (plus de notation 3e-04) ; marge basse pour p
    #scale_x_continuous(labels = scales::label_percent(accuracy = 1)) +
    scale_y_continuous(#labels = scales::label_percent(accuracy = 0.1),
                       expand = expansion(mult = c(0.12, 0.05))) +
    labs(
      y = "Relative change in NT",
      x = "Relative change in net energy"
    ) +
    theme_lollipop
  p10 <- ggplot(dataset_clean_Ne_Enet_diff) +
    aes(y = `T/Vk.`, x = ratioEnet_wbeta) +
    # [AJOUT] repères 0 : signe des deux variations lisible directement
    geom_hline(yintercept = 0, colour = "grey60", linetype = "dashed", linewidth = 0.4) +
    geom_vline(xintercept = 0, colour = "grey60", linetype = "dashed", linewidth = 0.4) +
    # [MODIFIÉ] ruban plus léger, droite plus fine
    geom_smooth(color = "#e41a1c", fill = "#69b3a2", alpha = 0.2, se = TRUE, method = "lm",
                linewidth = 0.9, formula = y ~ x) +
    # [MODIFIÉ] points à contour blanc (se détachent du ruban)
    geom_point(size = 2.6, shape = 21, fill = "grey15", colour = "white", stroke = 0.5) +
    # [MODIFIÉ] étiquettes plus petites, traits discrets, placement reproductible (seed)
    geom_text_repel(aes(label = Species_lab), size = 3.1, color = "grey25",
                    max.overlaps = Inf, box.padding = 0.4, point.padding = 0.25,
                    min.segment.length = 0, segment.color = "grey65", segment.size = 0.3, seed = 1) +
    # [AJOUT] p de la régression, coin inférieur droit
    annotate("text", x = Inf, y = -Inf, hjust = 1.05, vjust = -0.6, size = 3.3, colour = "grey30", parse = TRUE,
             label = p_label("T/Vk.")) +
    # [AJOUT] axes en % (plus de notation 3e-04) ; marge basse pour p
    #scale_x_continuous(labels = scales::label_percent(accuracy = 1)) +
    scale_y_continuous(#labels = scales::label_percent(accuracy = 0.1),
                       expand = expansion(mult = c(0.12, 0.05))) +
    labs(
      y = "Relative change in T/Vk.",
      x = "Relative change in net energy"
    ) +
    theme_lollipop
  
  # [MODIFIÉ] lettres de panneaux (a, b, ...) et alignement des axes
  compo <- ggarrange(p9, p8, p4, pblank, p10,p7,p6,p5, nrow=2, ncol=4,
                     labels = c("a", "b", "c", "", "d", "e", "f", "g"),
                     font.label = list(size = 13, face = "bold"), align = "hv")
  ggsave(compo, file=paste0(PATHOUT, "/", enet_dimension,"Temp_Nalpha.png"),
         height = 7, width=15, dpi = 300, device = ragg::agg_png, bg = "white")   # [MODIFIÉ] ragg : glyphe "ₑ"
  
  compo <- ggarrange(p4,p10,p6,p5, ncol=2, nrow=2,
                     labels = c("a", "b", "c", "d"),
                     font.label = list(size = 13, face = "bold"), align = "hv")
  ggsave(compo, file=paste0(PATHOUT, "/", enet_dimension,"Temp_Nalpha_2.png"),
         height = 7.5, width=8.5, dpi = 300, device = ragg::agg_png, bg = "white")
}
# table
# tablelog
tablediff # the one used in the supp mat
# tabledifflog



