library(MASS)
library(tidyverse)
library(caret)
library(mda)
library(klaR)
library(ggrepel)
library(broom)
library(pls)
# library(ggord)
library(factoextra)
library(plsVarSel)
library(ggpubr)
library(rstatix)

main = "New_gpu_conf.R" #config file
source(paste0(here::here(), "/R/functions/source.R")) #load config and key parameters and functions
source(paste0(here::here(), "/R/functions/functions_key_traits.R")) #plots against Fmax


# PREPARE DATASETS -------------

Ne_NT_tot = NULL
Ne_NT = NULL 
Ne_NT_tot2 = NULL
Ne_NT2 = NULL 


## LOOP TO EXTRACT DATA BY SCENARIO ---------
for (Scenario in c("Chap2", "Chap3")){
  cat("=======> ", Scenario, "\n")
  PATHOUT  = paste0("output/PLOTS/", toupper(Scenario), "/KEY_TRAITS_MATURE_removed0inNedata/")
  if (!dir.exists(PATHOUT)) dir.create(PATHOUT, recursive = TRUE)
  if (!dir.exists(paste0(PATHOUT, "/Clean_for_paper"))) dir.create(paste0(PATHOUT, "/Clean_for_paper"), recursive = TRUE)
  pathoutputRDATA = paste0(PATHOUT, "RDATA")
  if(!dir.exists(pathoutputRDATA)) dir.create(pathoutputRDATA, recursive = T)
  base_path = paste0("output/DATASETS/", Scenario, "/Trait_variables_by_cohort_by_age_withWbyAge/")
  
  ## loading the dataset of traits averaged over cohorts Time1 to Time2, and reorganised to remove age=0, add reassign natural mortality as add mortality, and estimate natural mortality
  data_complete_c <- load_trait_df(base_path = base_path, pathout=pathoutputRDATA, FORCE=F)
  
  for (sp in species){ #load and filter data 
    load(file=paste0(pathoutputRDATA, "/Complete_INIT2", sp, ".RData")) 
    final_dataset_age = complete[[2]] 
    rm(complete)
    NeNT <- final_dataset_age[which(final_dataset_age$variable %in% c("Deltasquare", 
                                                                      "Dead", "Dead_F", "Dead_Madd", "Dead_Mfor", "Dead_Mnat", "Dead_Mpred", 
                                                                      "Dead_Mstar", "Dead_Mtot")), 
                              -which(colnames(final_dataset_age) %in% c("Fstart"))]     
    NeNT2 <- final_dataset_age[which(final_dataset_age$variable %in% c("Ne_N_alpha", "NT", "generation_length", "Vk_point")), 
                              -which(colnames(final_dataset_age) %in% c("Fstart", "age"))] 
    rm(final_dataset_age) 
    if (is.null(Ne_NT)) {
      Ne_NT <- pivot_wider(NeNT, values_from = "value", names_from = "variable") 
      Ne_NT2 <- pivot_wider(NeNT2, values_from = "value", names_from = "variable") 
    } else {
      Ne_NT <- rbind(Ne_NT, pivot_wider(NeNT, values_from = "value", names_from = "variable")) 
      Ne_NT2 <- rbind(Ne_NT2, pivot_wider(NeNT2, values_from = "value", names_from = "variable")) 
    }
  }
  Ne_NT[,"Scenario"] <- ifelse(Scenario == "Chap2", "No CC", "With CC") # *********
  Ne_NT2[,"Scenario"] <- ifelse(Scenario == "Chap2", "No CC", "With CC") # *********
  
  if (is.null(Ne_NT_tot)) {Ne_NT_tot <- Ne_NT 
  } else Ne_NT_tot <- rbind(Ne_NT_tot, Ne_NT)   
  if (is.null(Ne_NT_tot2)) {Ne_NT_tot2 <- Ne_NT2 
  } else Ne_NT_tot2 <- rbind(Ne_NT_tot2, Ne_NT2) 
}

## ADDING THE NE/NT variable -----------
Ne_NT_tot2 <- Ne_NT_tot2[-which(Ne_NT_tot2$NT < 10), ]
Ne_NT_tot2$Ne_NT <- Ne_NT_tot2$Ne_N_alpha/Ne_NT_tot2$NT 
Ne_NT_tot2$'T/Vk' <- Ne_NT_tot2$generation_length/Ne_NT_tot2$Vk_point 

## ADDING THE INTEGRATED DeLTA² AND MORTALITIES --------
Ne_NT_tot <- pivot_longer(Ne_NT_tot, values_to = "value", names_to = "variable", cols = -c("age", "Fmax", "replicate", "Species", "Scenario"))
Ne_NT_tot$Fstart <- 1
Deltasquare <- estimatecolSum(trait_data_full_i = Ne_NT_tot, var = "Deltasquare")
Dead        <- estimatecolSum(trait_data_full_i = Ne_NT_tot, var = "Dead")
Dead_F      <- estimatecolSum(trait_data_full_i = Ne_NT_tot, var = "Dead_F")
Dead_Madd   <- estimatecolSum(trait_data_full_i = Ne_NT_tot, var = "Dead_Madd")
Dead_Mfor   <- estimatecolSum(trait_data_full_i = Ne_NT_tot, var = "Dead_Mfor")
Dead_Mnat   <- estimatecolSum(trait_data_full_i = Ne_NT_tot, var = "Dead_Mnat")
Dead_Mpred  <- estimatecolSum(trait_data_full_i = Ne_NT_tot, var = "Dead_Mpred")
Dead_Mstar  <- estimatecolSum(trait_data_full_i = Ne_NT_tot, var = "Dead_Mstar")
Dead_Mtot   <- estimatecolSum(trait_data_full_i = Ne_NT_tot, var = "Dead_Mtot")
Ne_NT_tot   <- pivot_wider(Ne_NT_tot, values_from = "value", names_from = "variable")

# PLOT NE/NT and NT -------------

Ne_NT_tot_slope <- Ne_NT_tot2 %>% 
  group_by(Species, replicate, Scenario) %>% 
  nest() %>% 
  mutate(model = map(data, ~lm(Ne_NT ~ 1 + Fmax, data = .x) %>% tidy)) %>% 
  unnest(model) %>% 
  filter(term == 'Fmax') %>% ungroup() 
Ne_NT_tot_slope <- Ne_NT_tot_slope[, c("Species", "estimate", "replicate", "Scenario")] 
Ne_NT_tot_slope$Species = factor(Ne_NT_tot_slope$Species, levels=species) # ordering the species as in all other plots
colnames(Ne_NT_tot_slope)[which(colnames(Ne_NT_tot_slope)=="estimate")] <- "slope"

Ne_NT_tot_Variability <- Ne_NT_tot2 %>% 
  group_by(Species, replicate, Scenario) %>% 
  summarize(Variability = var(Ne_NT, na.rm = T)) %>% #max(Ne_NT, na.rm = T) - min(Ne_NT, na.rm = T)) %>%
  ungroup() 
Ne_NT_tot_Variability <- Ne_NT_tot_Variability[, c("Species", "Variability", "replicate", "Scenario")] 
Ne_NT_tot_Variability$Species = factor(Ne_NT_tot_Variability$Species, levels=species) # ordering the species as in all other plots
colnames(Ne_NT_tot_Variability)[which(colnames(Ne_NT_tot_Variability)=="estimate")] <- "Ne_NT_Variability"

Ne_NT_tot_NT_slope <- Ne_NT_tot2 %>% 
  group_by(Species, replicate, Scenario) %>% 
  nest() %>% 
  mutate(model = map(data, ~lm(NT ~ 1 + Fmax, data = .x) %>% tidy)) %>% 
  unnest(model) %>% 
  filter(term == 'Fmax') %>% ungroup() 
Ne_NT_tot_NT_slope <- Ne_NT_tot_NT_slope[, c("Species", "estimate", "replicate", "Scenario")] 
Ne_NT_tot_NT_slope$Species = factor(Ne_NT_tot_NT_slope$Species, levels=species) # ordering the species as in all other plots
colnames(Ne_NT_tot_NT_slope)[which(colnames(Ne_NT_tot_NT_slope)=="estimate")] <- "slope"

Ne_NT_tot_NT <- Ne_NT_tot2 %>% 
  group_by(Species, replicate, Scenario) %>% 
  summarize(Mean_NT = mean(NT, na.rm = T)) %>% #max(Ne_NT, na.rm = T) - min(Ne_NT, na.rm = T)) %>% %>% 
  ungroup() 
Ne_NT_tot_NT <- Ne_NT_tot_NT[, c("Species", "Mean_NT", "replicate", "Scenario")] 
Ne_NT_tot_NT$Species = factor(Ne_NT_tot_NT$Species, levels=species) # ordering the species as in all other plots

Ne_NT_tot_Ne_slope <- Ne_NT_tot2 %>% 
  group_by(Species, replicate, Scenario) %>% 
  nest() %>% 
  mutate(model = map(data, ~lm(Ne_N_alpha ~ 1 + Fmax, data = .x) %>% tidy)) %>% 
  unnest(model) %>% 
  filter(term == 'Fmax') %>% ungroup() 
Ne_NT_tot_Ne_slope <- Ne_NT_tot_Ne_slope[, c("Species", "estimate", "replicate", "Scenario")] 
Ne_NT_tot_Ne_slope$Species = factor(Ne_NT_tot_Ne_slope$Species, levels=species) # ordering the species as in all other plots
colnames(Ne_NT_tot_Ne_slope)[which(colnames(Ne_NT_tot_Ne_slope)=="estimate")] <- "slope"

Ne_NT_tot_Ne <- Ne_NT_tot2 %>% 
  group_by(Species, replicate, Scenario) %>% 
  summarize(Mean_Ne = mean(Ne_N_alpha, na.rm = T)) %>% #max(Ne_NT, na.rm = T) - min(Ne_NT, na.rm = T)) %>% %>% 
  ungroup() 
Ne_NT_tot_Ne <- Ne_NT_tot_Ne[, c("Species", "Mean_Ne", "replicate", "Scenario")] 
Ne_NT_tot_Ne$Species = factor(Ne_NT_tot_Ne$Species, levels=species) # ordering the species as in all other plots



## PLOT : differentiate the slope of Ne/NT between scenarios with and without CC -----
stat.test <- Ne_NT_tot_slope %>%
  group_by(Species) %>%
  t_test(slope ~ Scenario) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance()
stat.test 
stat.test <- stat.test %>% add_xy_position(x = "Scenario")
stat.test$y.position <- stat.test$y.position*ifelse(stat.test$Species == "AtlanticCod", 1, 4)
ncolors <- nrow(Ne_NT_tot_slope)
my_palette <- grDevices::colorRampPalette(c("lightblue","black","red"))(ncolors)
my_df <- data.frame(value=Ne_NT_tot_slope$slope, color=my_palette[cut(Ne_NT_tot_slope$slope,ncolors)])
Ne_NT_tot_slope$col <- my_df$color
ggbarplot(
  Ne_NT_tot_slope, x = "Scenario", y = "slope", fill= "#00AFBB", #fill = "col", #
  facet.by = "Species"
  ) +
  theme(legend.position="none") +
  stat_pvalue_manual(stat.test) + 
  geom_abline(slope = 0, intercept = 0, col="grey50")

## PLOT : differentiate the Variability of values of Ne_NT between scenarios with and without CC ------
stat.test <- Ne_NT_tot_Variability %>%
  group_by(Species) %>%
  t_test(Variability ~ Scenario) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance()
stat.test <- stat.test %>% add_xy_position(x = "Scenario")
stat.test$y.position <- stat.test$y.position*ifelse(stat.test$Species == "AtlanticCod", 0.8, 1.1)
ggbarplot(
  Ne_NT_tot_Variability, x = "Scenario", y = "Variability", fill = "#00AFBB", 
  facet.by = "Species"
  ) +
  stat_pvalue_manual(stat.test) + 
  geom_abline(slope = 0, intercept = 0, col="grey50")

## PLOT : How different is the mean and slope of NT between the two scenarios ?  -------
stat.test <- Ne_NT_tot_NT_slope %>%
  group_by(Species) %>%
  t_test(slope ~ Scenario) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance()
stat.test 
stat.test <- stat.test %>% add_xy_position(x = "Scenario")
stat.test$y.position <- stat.test$y.position*ifelse(stat.test$Species == "AtlanticCod", 1, 4)
ncolors <- nrow(Ne_NT_tot_NT_slope)
my_palette <- grDevices::colorRampPalette(c("lightblue","black","red"))(ncolors)
my_df <- data.frame(value=Ne_NT_tot_NT_slope$slope, color=my_palette[cut(Ne_NT_tot_NT_slope$slope,ncolors)])
Ne_NT_tot_NT_slope$col <- my_df$color
ggboxplot(
  Ne_NT_tot_NT_slope, x = "Scenario", y = "slope", fill= "#00AFBB", #fill = "col", #
  facet.by = "Species"
) +
  theme(legend.position="none") +
  stat_pvalue_manual(stat.test) + 
  geom_abline(slope = 0, intercept = 0, col="grey50")

Ne_NT_tot_NT$log_NT <- log(Ne_NT_tot_NT$Mean_NT)
stat.test <- Ne_NT_tot_NT %>%
  group_by(Species) %>%
  t_test(log_NT ~ Scenario) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance()
stat.test <- stat.test %>% add_xy_position(x = "Scenario")
stat.test$y.position <- stat.test$y.position*20
ggbarplot(
  Ne_NT_tot_NT, x = "Scenario", y = "log_NT", fill = "#00AFBB", 
  facet.by = "Species"
) +
  stat_pvalue_manual(stat.test) + 
  geom_abline(slope = 0, intercept = 0, col="grey50")

## PLOT : How different is the mean and slope of Ne between the two scenarios ?  -------
stat.test <- Ne_NT_tot_Ne_slope %>%
  group_by(Species) %>%
  t_test(slope ~ Scenario) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance()
stat.test 
stat.test <- stat.test %>% add_xy_position(x = "Scenario")
stat.test$y.position <- stat.test$y.position*ifelse(
  stat.test$Species %in% c("Sandeel", "EuropeanSprat", "GreyGurnard"), 1, ifelse(  
  stat.test$Species %in% c("CommonSole", "EuropeanPlaice", "Saithe", "AtlanticCod", "Haddock", "Hake", "CommonDab"),c(0), c(-1)))
ggboxplot(
  Ne_NT_tot_Ne_slope, x = "Scenario", y = "slope", fill= "#00AFBB" #fill = "col", #
) +
  facet_wrap(~Species,scale="free")+
  theme(legend.position="none") +
  stat_pvalue_manual(stat.test) + 
  geom_abline(slope = 0, intercept = 0, col="grey50")+
  labs(title = "Boxplot of the SLOPE of Ne with Fmax, with the replicates, by species",
       x="")


Ne_NT_tot_Ne$log_Ne <- log(Ne_NT_tot_Ne$Mean_Ne)
stat.test <- Ne_NT_tot_Ne %>%
  group_by(Species) %>%
  t_test(log_Ne ~ Scenario) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance()
stat.test <- stat.test %>% add_xy_position(x = "Scenario")
stat.test$y.position <- stat.test$y.position*ifelse(stat.test$Species == "AtlanticCod", 0.96, 0.96)
ggboxplot(
  Ne_NT_tot_Ne, x = "Scenario", y = "log_Ne", fill = "#00AFBB"
) +
  facet_wrap(~Species,scale="free")+
  stat_pvalue_manual(stat.test) + 
  geom_abline(slope = 0, intercept = 0, col="grey50")+
  labs(title = "Boxplot of the MEAN of Ne with Fmax, with the replicates, by species",
       x="")


# PLOT Vk and T -------------


gene_length <- Ne_NT_tot2 %>% 
  group_by(Species, replicate, Scenario) %>% 
  nest() %>% 
  mutate(model = map(data, ~lm(generation_length ~ 1 + Fmax, data = .x) %>% tidy)) %>% 
  unnest(model) %>% 
  filter(term == 'Fmax') %>% ungroup() 
gene_length <- gene_length[, c("Species", "estimate", "replicate", "Scenario")] 
gene_length$Species = factor(gene_length$Species, levels=species) # ordering the species as in all other plots
colnames(gene_length)[which(colnames(gene_length)=="estimate")] <- "slope"

gene_length2 <- Ne_NT_tot2 %>% 
  group_by(Species, replicate, Scenario) %>% 
  summarize(mean_gene_length = mean(generation_length, na.rm = T)) %>% #max(Ne_NT, na.rm = T) - min(Ne_NT, na.rm = T)) %>% %>% 
  ungroup() 
gene_length2 <- gene_length2[, c("Species", "mean_gene_length", "replicate", "Scenario")] 
gene_length2$Species = factor(gene_length2$Species, levels=species) # ordering the species as in all other plots

Vkp <- Ne_NT_tot2 %>% 
  group_by(Species, replicate, Scenario) %>% 
  nest() %>% 
  mutate(model = map(data, ~lm(Vk_point ~ 1 + Fmax, data = .x) %>% tidy)) %>% 
  unnest(model) %>% 
  filter(term == 'Fmax') %>% ungroup() 
Vkp <- Vkp[, c("Species", "estimate", "replicate", "Scenario")] 
Vkp$Species = factor(Vkp$Species, levels=species) # ordering the species as in all other plots
colnames(Vkp)[which(colnames(Vkp)=="estimate")] <- "slope"

Vkp2 <- Ne_NT_tot2 %>% 
  group_by(Species, replicate, Scenario) %>% 
  summarize(mean_Vk_point = mean(Vk_point, na.rm = T)) %>% #max(Ne_NT, na.rm = T) - min(Ne_NT, na.rm = T)) %>% %>% 
  ungroup() 
Vkp2 <- Vkp2[, c("Species", "mean_Vk_point", "replicate", "Scenario")] 
Vkp2$Species = factor(Vkp2$Species, levels=species) # ordering the species as in all other plots


gonV <- Ne_NT_tot2 %>% 
  group_by(Species, replicate, Scenario) %>% 
  nest() %>% 
  mutate(model = map(data, ~lm(`T/Vk` ~ 1 + Fmax, data = .x) %>% tidy)) %>% 
  unnest(model) %>% 
  filter(term == 'Fmax') %>% ungroup() 
gonV <- gonV[, c("Species", "estimate", "replicate", "Scenario")] 
gonV$Species = factor(gonV$Species, levels=species) # ordering the species as in all other plots
colnames(gonV)[which(colnames(gonV)=="estimate")] <- "slope"

gonV2 <- Ne_NT_tot2 %>% 
  group_by(Species, replicate, Scenario) %>% 
  summarize(gonV = mean(`T/Vk`, na.rm = T)) %>% #max(Ne_NT, na.rm = T) - min(Ne_NT, na.rm = T)) %>% %>% 
  ungroup() 
gonV2 <- gonV2[, c("Species", "gonV", "replicate", "Scenario")] 
gonV2$Species = factor(gonV2$Species, levels=species) # ordering the species as in all other plots



## PLOT ---
stat.test <- gene_length %>%
  group_by(Species) %>%
  t_test(slope ~ Scenario) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance()
stat.test
stat.test <- stat.test %>% add_xy_position(x = "Scenario")
stat.test$y.position <- stat.test$y.position
ggboxplot(
  gene_length, x = "Scenario", y = "slope", fill = "#00AFBB",
  facet.by = "Species"
) +
  stat_pvalue_manual(stat.test) +
  geom_abline(slope = 0, intercept = 0, col="grey50")+
  labs(title="Comparing the SLOPE of the generation Length (T) with and without CC")

stat.test <- gene_length2 %>%
  group_by(Species) %>%
  t_test(mean_gene_length ~ Scenario) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance()
stat.test
stat.test <- stat.test %>% add_xy_position(x = "Scenario")
stat.test$y.position <- stat.test$y.position
ggboxplot(
  gene_length2, x = "Scenario", y = "mean_gene_length", fill = "#00AFBB",
  facet.by = "Species"
) +
  stat_pvalue_manual(stat.test) +
  geom_abline(slope = 0, intercept = 0, col="grey50")+
  labs(title="Comparing Generation Length (T) with and without CC")


stat.test <- Vkp %>%
  group_by(Species) %>%
  t_test(slope ~ Scenario) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance()
stat.test
stat.test <- stat.test %>% add_xy_position(x = "Scenario")
stat.test$y.position <- stat.test$y.position
ggboxplot(
  Vkp, x = "Scenario", y = "slope", fill = "#00AFBB",
  facet.by = "Species"
) +
  stat_pvalue_manual(stat.test) +
  geom_abline(slope = 0, intercept = 0, col="grey50")+
  labs(title="Comparing the SLOPE of Vk_point with and without CC")

stat.test <- Vkp2 %>%
  group_by(Species) %>%
  t_test(mean_Vk_point ~ Scenario) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance()
stat.test
stat.test <- stat.test %>% add_xy_position(x = "Scenario")
stat.test$y.position <- stat.test$y.position
ggboxplot(
  Vkp2, x = "Scenario", y = "mean_Vk_point", fill = "#00AFBB",
  facet.by = "Species"
) +
  stat_pvalue_manual(stat.test) +
  geom_abline(slope = 0, intercept = 0, col="grey50")+
  labs(title="Comparing Vk_point with and without CC")

stat.test <- gonV %>%
  group_by(Species) %>%
  t_test(slope ~ Scenario) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance()
stat.test
stat.test <- stat.test %>% add_xy_position(x = "Scenario")
stat.test$y.position <- stat.test$y.position
ggboxplot(
  gonV, x = "Scenario", y = "slope", fill = "#00AFBB",
  facet.by = "Species"
) +
  stat_pvalue_manual(stat.test) +
  geom_abline(slope = 0, intercept = 0, col="grey50")+
  labs(title="Comparing the SLOPE of T/Vk with and without CC")

stat.test <- gonV2 %>%
  group_by(Species) %>%
  t_test(gonV ~ Scenario) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance()
stat.test
stat.test <- stat.test %>% add_xy_position(x = "Scenario")
stat.test$y.position <- stat.test$y.position
ggbarplot(
  gonV2, x = "Scenario", y = "gonV", fill = "#00AFBB",
  facet.by = "Species"
) +
  stat_pvalue_manual(stat.test) +
  geom_abline(slope = 0, intercept = 0, col="grey50")+
  labs(title="Comparing T/Vk with and without CC")



# PLOT Delta² and Dead -------------

Deltasquare$log_value <- log(Deltasquare$value)
Deltasquare <- Deltasquare[-which(is.infinite(Deltasquare$log_value)), ]
Deltasquare <- Deltasquare %>% 
  group_by(Species, replicate, Scenario) %>% 
  summarize(log_Deltasquare = mean(log_value, na.rm = T)) %>% #max(Ne_NT, na.rm = T) - min(Ne_NT, na.rm = T)) %>% %>% 
  ungroup() 
Deltasquare <- Deltasquare[, c("Species", "log_Deltasquare", "replicate", "Scenario")] 
Deltasquare$Species = factor(Deltasquare$Species, levels=species) # ordering the species as in all other plots



##PLOT

stat.test <- Deltasquare %>%
  group_by(Species) %>%
  t_test(log_Deltasquare ~ Scenario) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance()
stat.test
stat.test <- stat.test %>% add_xy_position(x = "Scenario")
stat.test$y.position <- stat.test$y.position
ggboxplot(
  Deltasquare, x = "Scenario", y = "log_Deltasquare", fill = "#00AFBB",
  facet.by = "Species"
) +
  stat_pvalue_manual(stat.test) +
  geom_abline(slope = 0, intercept = 0, col="grey50")+
  labs(title="Comparing Deltasquare with and without CC")
