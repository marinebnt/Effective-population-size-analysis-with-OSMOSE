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


# STEPS ------------

#* We are willing to understand the consequences of the climate change (CC) on Ne/NT and the underlying processes, when increasing the fishing pressure. 
#* Initial dataset of the effective population size: by Species, Variable, Replicate, Fmax and Scenario.
#* Initial dataset of the rest of the functional traits : by Age, Species, Variable, Replicate, Fmax and Scenario.
#* Here is a description of the steps we used to produce the heatmaps : 
#*    1- If there is an age dimension : integrate between ages the values to have the cumulated variable between ages. 
#*       Dimensions left: Species, Variable, Replicate, Fmax and Scenario
#*    2- Estimate the slope as: lm(value~1+Fmax)
#*       Dimensions left: Species, Variable, Replicate and Scenario
#*    3- Estimate the significance of the difference of the slope between the no CC scenario and the CC scenario (t.test, slope~Scenario)
#*       No dimension removed and dimension of the significance of the slope difference : Species and Variable
#*    4- Difference between No CC and CC scenario ((abs(`With CC`) - abs(`No CC`))/abs(`No CC`) value)
#*       Dimensions left: Species, Variable and Replicate
#*    5- Mean of the responses between replicates
#*       Dimensions left: Species and Variable


# PREPARE DATASETS -------------

Ne_NT_tot = NULL
Ne_NT_tot2 = NULL

## LOOP TO EXTRACT DATA BY SCENARIO ---------
for (Scenario in c("Chap2", "Chap3")){
  Ne_NT = NULL 
  Ne_NT2 = NULL 
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
    NeNT <- final_dataset_age[which(final_dataset_age$variable %in% c("Deltasquare", "eggNumber_sc", "Mpred", "Mstar",
                                                                      "Dead", "Dead_F", "Dead_Madd", "Dead_Mfor", "Dead_Mnat", "Dead_Mpred", 
                                                                      "Dead_Mstar", "Dead_Mtot", "survivorship_curve", "SSDgx")), 
                              -which(colnames(final_dataset_age) %in% c("Fstart"))] %>% ungroup()     
    NeNT2 <- final_dataset_age[which(final_dataset_age$variable %in% c("Ne_N_alpha", "NT", "generation_length", "Vk_point", "Numerator_N_alpha")), 
                               -which(colnames(final_dataset_age) %in% c("Fstart", "age"))] %>% ungroup()    
    numerator_n_alpha <- NULL
    for (x in 1:30) {
      data_tempo <- final_dataset_age[which(final_dataset_age$variable=="Numerator_N_alpha"),]
      data_tempo$age <- x
      if (is.null(numerator_n_alpha)) numerator_n_alpha <- data_tempo
      else numerator_n_alpha <- rbind(numerator_n_alpha, data_tempo)
    }
    SSDgNalpha <- rbind(final_dataset_age[which(final_dataset_age$variable=="SSDgx"),], numerator_n_alpha)
    SSDgNalpha <- pivot_wider(SSDgNalpha, names_from = "variable", values_from = "value")
    SSDgNalpha$value <- SSDgNalpha$SSDgx/SSDgNalpha$Numerator_N_alpha
    SSDgNalpha$variable <- "SSDgx/Numerator_N_alpha"
    NeNT <- rbind(SSDgNalpha[, which(colnames(SSDgNalpha) %in% colnames(NeNT))], NeNT)
    
    rm(final_dataset_age) 
    NeNT[,"Scenario"] <- ifelse(Scenario == "Chap2", "No CC", "With CC") # *********
    NeNT2[,"Scenario"] <- ifelse(Scenario == "Chap2", "No CC", "With CC") # *********
    
    cat("\nSpecies: ", sp, " dim dataset: ", dim(NeNT2))
    
    if (is.null(Ne_NT)) {
      Ne_NT <- pivot_wider(NeNT, values_from = "value", names_from = "variable") 
      Ne_NT2 <- pivot_wider(NeNT2, values_from = "value", names_from = "variable") 
      Ne_NT2 |>
        dplyr::summarise(n = dplyr::n(), .by = c(Fmax, replicate, Species, Scenario)) |>
        dplyr::filter(n > 1L) 
    } else {
      Ne_NT <- rbind(Ne_NT, pivot_wider(NeNT, values_from = "value", names_from = "variable")) 
      Ne_NT2 <- rbind(Ne_NT2, pivot_wider(NeNT2, values_from = "value", names_from = "variable")) 
      dup <- Ne_NT2 |>
        dplyr::summarise(n = dplyr::n(), .by = c(Fmax, replicate, Species, Scenario)) |>
        dplyr::filter(n > 1L) 
      if(nrow(dup)>0) stop()
    }
  }#sp
  if (is.null(Ne_NT_tot)) {Ne_NT_tot <- Ne_NT 
  } else Ne_NT_tot <- rbind(Ne_NT_tot, Ne_NT)   
  if (is.null(Ne_NT_tot2)) {Ne_NT_tot2 <- Ne_NT2 
  } else Ne_NT_tot2 <- rbind(Ne_NT_tot2, Ne_NT2) 
}#scenario


## ADDING THE NE/NT variable -----------
Ne_NT_tot2 <- Ne_NT_tot2[-which(Ne_NT_tot2$NT < 10), ]
Ne_NT_tot2$Ne_NT <- Ne_NT_tot2$Ne_N_alpha/Ne_NT_tot2$NT 
Ne_NT_tot2$'T/Vk' <- Ne_NT_tot2$generation_length/Ne_NT_tot2$Vk_point 

# Ne/NT and related variales -----------------
## ESTIMATE SLOPE OF ALL VARIABLES -------------
Ne_NT_tot2 <- pivot_longer(Ne_NT_tot2, values_to = "value", names_to = "variable", cols = -c("Fmax", "replicate", "Species", "Scenario"))
#estimate the slope by variable/species/scenario/replicate


# center the values around the values of the Fmax=1
Ne_NT_tot2_centered <- Ne_NT_tot2 %>%
  group_by(Species, variable, Scenario, replicate) %>%
  mutate(center = value - value[which(Fmax == 1)]) %>%
  ungroup()
Ne_NT_tot2_centered$Fmax <- as.numeric(Ne_NT_tot2_centered$Fmax)
Ne_NT_tot2_centered$Scenario <- as.factor(Ne_NT_tot2_centered$Scenario)
# Ne_NT_tot2_centered <- Ne_NT_tot2_centered[, -which(colnames(Ne_NT_tot2_centered) == "value")]

# make slope estimate relative
# Ne_NT_tot2_centered <- pivot_wider(Ne_NT_tot2_centered, values_from = "center", names_from = "Scenario")
# Ne_NT_tot2_centered$`With CC` <- Ne_NT_tot2_centered$`With CC`/Ne_NT_tot2_centered$`No CC`
# Ne_NT_tot2_centered$`No CC` <- Ne_NT_tot2_centered$`No CC`/Ne_NT_tot2_centered$`No CC`
# Ne_NT_tot2_centered <- pivot_longer(Ne_NT_tot2_centered, cols=c("No CC", "With CC"), 
#                                     names_to = "Scenario", values_to = "center")

# diagnostics
ggplotdata <- Ne_NT_tot2_centered
ggplotdata$logslope <- log(ggplotdata$center-min(ggplotdata$center, na.rm = T)+0.01)
ggplotdata$id <- 1:nrow(ggplotdata) 
ggplot(ggplotdata)+ 
  geom_density(aes(center))+ 
  facet_wrap(~variable, scales="free")
ggplot(ggplotdata, aes(sample=center))+ 
  stat_qq()+ 
  stat_qq_line()+ 
  facet_wrap(~variable, scales="free")

# ggplotdata$logslope <- log(ggplotdata$value-min(ggplotdata$value, na.rm = T)+0.01)



# for (sp in species){
#   model_Ne_Nt <- lm(data=Ne_NT_tot2_centered[which(Ne_NT_tot2_centered$Species==sp),], center~Scenario+Fmax+Scenario:Fmax)
#   coef <- coef(model_Ne_Nt)
#   coef["ScenarioWith CC:Fmax"]
# }

# if(!all(table(Ne_NT_tot2_centered$Fmax, Ne_NT_tot2_centered$Scenario, Ne_NT_tot2_centered$Species)==140)) {
#   ID_pb   <- which(!table(Ne_NT_tot2_centered$Fmax, Ne_NT_tot2_centered$Scenario, Ne_NT_tot2_centered$Species)==140, arr.ind = T)
#   nb_pb   <- 140-table(Ne_NT_tot2_centered$Fmax, Ne_NT_tot2_centered$Scenario, Ne_NT_tot2_centered$Species)[ID_pb]
#   name_id <- dimnames(table(Ne_NT_tot2_centered$Fmax, Ne_NT_tot2_centered$Scenario, Ne_NT_tot2_centered$Species))
#   dim_1 <- name_id[[1]][ID_pb[1]]
#   dim_2 <- name_id[[2]][ID_pb[2]]
#   dim_3 <- name_id[[3]][ID_pb[3]]
#   list = list()
#   for (k in 1:28) list[[k]] = data.frame("replicate"=1, "Scenario"=dim_2, "Fmax"=dim_1,"Species"=dim_3, "center"=0, "value"=0)
#   df <- do.call(rbind, list)
#   df <- cbind(df, "variable"=rep(names(table(Ne_NT_tot2_centered$variable)), nb_pb/length(names(table(Ne_NT_tot2_centered$variable)))))
#   Ne_NT_tot2_centered <- rbind(Ne_NT_tot2_centered, df)
# }

Ne_NT_tot_slope <- Ne_NT_tot2_centered %>% 
  group_by(Species, variable) %>% 
  nest() %>% 
  mutate(model = map(data, ~lm(center ~ Fmax + Scenario + Fmax:Scenario, data = .x) %>% tidy)) %>% 
  unnest(model) %>% 
  # filter(term %in% c('Fmax:ScenarioWith CC', "(Intercept)")) %>%
  ungroup() 

# if slope Fmax decreasing : response of intercept and slope of noCC/wCC should be reversed
Ne_NT_tot_slope_abs <- pivot_wider(Ne_NT_tot_slope[, which(colnames(Ne_NT_tot_slope) %in% 
                                                                      c("Species", "variable", "estimate", "term"))], 
                                   values_from = "estimate", names_from = "term")
Ne_NT_tot_slope_abs$`(Intercept)` = Ne_NT_tot_slope_abs$`(Intercept)`*Ne_NT_tot_slope_abs$Fmax
Ne_NT_tot_slope_abs$`Fmax:ScenarioWith CC` = Ne_NT_tot_slope_abs$`ScenarioWith CC`*sign(Ne_NT_tot_slope_abs$Fmax)
Ne_NT_tot_slope_abs <- pivot_longer(Ne_NT_tot_slope_abs, cols=c("Fmax:ScenarioWith CC", "(Intercept)"), values_to = "estimate", names_to="term")
Ne_NT_tot_slope <- right_join(Ne_NT_tot_slope[, which(colnames(Ne_NT_tot_slope) %in% c("Species", "variable","p.value", "term"))], 
                             Ne_NT_tot_slope_abs[, which(colnames(Ne_NT_tot_slope_abs) %in% 
                                                           c("Species", "variable", "estimate", "term"))])

#improve format table
Ne_NT_tot_slope <- Ne_NT_tot_slope[, c("Species", "estimate", "variable", "p.value", "term")] 
Ne_NT_tot_slope[,"p.adj.signif"] <- ifelse(Ne_NT_tot_slope$p.value < 0.05, "*", "")
colnames(Ne_NT_tot_slope)[which(colnames(Ne_NT_tot_slope)=="estimate")] <- "slope_difference"


# add Ne/NT response type when NO CC
Ne_NT_tot_slope$Response <- ifelse(Ne_NT_tot_slope$Species %in% c("Hake", "Whiting", "Saithe", "Sandeel", "AtlanticHerring"), "Increasing N\u2091/NT", 
                                   ifelse(Ne_NT_tot_slope$Species %in% c("CommonDab", "GreyGurnard", "EuropeanPlaice"), "Decreasing N\u2091/NT", 
                                          ifelse(Ne_NT_tot_slope$Species %in% c("NorwayPout", "CommonSole", "AtlanticCod", "EuropeanSprat", "Haddock", "HorseMackerel", 
                                                                                "AtlanticMackerel"), "Constant N\u2091/NT", "out"))) 

## PREPARE DATASET FOR PLOT ------------
# scale dataset
Ne_NT_tot_heatmap <- Ne_NT_tot_slope[, -which(colnames(Ne_NT_tot_slope)=="p.value")]
Ne_NT_tot_heatmap <- as.data.frame(pivot_wider(Ne_NT_tot_heatmap, values_from = "slope_difference", names_from = "variable"))
Ne_NT_tot_heatmap[, -which(colnames(Ne_NT_tot_heatmap) %in% c("Species", "Response", "p.adj.signif", "term"))] <- 
  scale(Ne_NT_tot_heatmap[, -which(colnames(Ne_NT_tot_heatmap) %in% c("Species", "Response", "p.adj.signif", "term"))])
Ne_NT_tot_heatmap <- pivot_longer(Ne_NT_tot_heatmap, cols= -c("Species", "Response", "p.adj.signif", "term"), values_to = "slope_difference", names_to = "variable")
Ne_NT_tot_heatmap$variable <- improve_vriable_names(Ne_NT_tot_heatmap$variable)
if (any(is.na(Ne_NT_tot_heatmap$slope_difference)))  Ne_NT_tot_heatmap <- Ne_NT_tot_heatmap[-which(is.na(Ne_NT_tot_heatmap$slope_difference)), ]
# order variables in the plot
Ne_NT_tot_heatmap$variable <- factor(Ne_NT_tot_heatmap$variable, levels = 
                                       c("N\u2091/NT", "N\u2091", "NT", "T/Vk", "T", "Vk.", "Nα"))
sp_to_order <- unique(as.matrix(Ne_NT_tot_heatmap[which(Ne_NT_tot_heatmap$variable == "N\u2091/NT"), "Species"]))
levelz = sp_to_order[order(as.matrix(Ne_NT_tot_heatmap[which(Ne_NT_tot_heatmap$variable == "N\u2091/NT"), "slope_difference"]))]
Ne_NT_tot_heatmap$Species = factor(Ne_NT_tot_heatmap$Species, levels= levelz) # ordering the species as in all other plots

## PLOT -------------
pNe1 <- ggplot(Ne_NT_tot_heatmap[which(Ne_NT_tot_heatmap$term=="Fmax:ScenarioWith CC"), ], aes(x=variable, y=Species))+
  geom_tile(aes(fill=slope_difference))+
  theme(
    axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5)) + 
  ggplot2::scale_fill_gradient2(
    low = "blue", 
    high = "red", 
    # oob = scales::squish,
    midpoint = mean(Ne_NT_tot_heatmap$slope_difference[which(Ne_NT_tot_heatmap$term=="Fmax:ScenarioWith CC")], na.rm=T), 
    limits = c(min(Ne_NT_tot_heatmap$slope_difference[which(Ne_NT_tot_heatmap$term=="Fmax:ScenarioWith CC")], na.rm=T),
               max(Ne_NT_tot_heatmap$slope_difference[which(Ne_NT_tot_heatmap$term=="Fmax:ScenarioWith CC")], na.rm=T))
    # breaks = c("Weaker slope with CC (-1)"=-1.0, "Unchanged (0)"=0, "Stronger slope with CC (>2)"=2.0)
  ) + 
  theme_classic()+
  facet_grid(rows = vars(Response), scales = "free", space = "free") +
  geom_text(aes(label=p.adj.signif), color="black", size=5) + 
  labs(title = "Scaled relative difference between without and with CC\nin slope of the species' functional traits'\nresponse to ↑ Fmax ",
       fill = "Difference in slope\nof the variable's\nresponse to ↑ Fmax") + 
  ylab("") + xlab("")

ggsave(filename = paste0(PATHOUT, "/Heatmaps/effective_pop_size_slope.jpg"), plot = pNe1,
       width= 9, height = 6)

pNe2 <- ggplot(Ne_NT_tot_heatmap[which(Ne_NT_tot_heatmap$term=="(Intercept)"), ], aes(x=variable, y=Species))+
  geom_tile(aes(fill=slope_difference))+
  theme(
    axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5)) + 
  ggplot2::scale_fill_gradient2(
    low = "blue", 
    high = "red", 
    # oob = scales::squish,
    midpoint = mean(Ne_NT_tot_heatmap$slope_difference[which(Ne_NT_tot_heatmap$term=="(Intercept)")], na.rm=T), 
    limits = c(min(Ne_NT_tot_heatmap$slope_difference[which(Ne_NT_tot_heatmap$term=="(Intercept)")], na.rm=T),
               max(Ne_NT_tot_heatmap$slope_difference[which(Ne_NT_tot_heatmap$term=="(Intercept)")], na.rm=T))
    # breaks = c("Weaker slope with CC (-1)"=-1.0, "Unchanged (0)"=0, "Stronger slope with CC (>2)"=2.0)
  ) + 
  theme_classic()+
  facet_grid(rows = vars(Response), scales = "free", space = "free") +
  geom_text(aes(label=p.adj.signif), color="black", size=5) + 
  labs(title = "Scaled relative difference between without and with CC\nin intercept of the species' functional traits'\nresponse to ↑ Fmax ",
       fill = "Difference in intercept\nof the variable's\nresponse to ↑ Fmax") + 
  ylab("") + xlab("")

ggsave(filename = paste0(PATHOUT, "/Heatmaps/effective_pop_size_intercept.jpg"), plot = pNe2,
       width= 9, height = 6)










# Functional traits -----------------
## ESTIMATE SLOPE OF ALL VARIABLES -------------
# remove the age dimension by integrating the variables between ages
Ne_NT_tot <- pivot_longer(Ne_NT_tot, values_to = "value", names_to = "variable", cols = -c("age", "Fmax", "replicate", "Species", "Scenario"))
Ne_NT_tot$Fstart <- 1
eggNumber_sc <- estimatecolSum(trait_data_full_i = Ne_NT_tot, var = "eggNumber_sc")
Deltasquare <- estimatecolSum(trait_data_full_i = Ne_NT_tot, var = "Deltasquare")
Dead        <- estimatecolSum(trait_data_full_i = Ne_NT_tot, var = "Dead")
Dead_F      <- estimatecolSum(trait_data_full_i = Ne_NT_tot, var = "Dead_F")
Dead_Madd   <- estimatecolSum(trait_data_full_i = Ne_NT_tot, var = "Dead_Madd")
Dead_Mfor   <- estimatecolSum(trait_data_full_i = Ne_NT_tot, var = "Dead_Mfor")
Dead_Mnat   <- estimatecolSum(trait_data_full_i = Ne_NT_tot, var = "Dead_Mnat")
Dead_Mpred  <- estimatecolSum(trait_data_full_i = Ne_NT_tot, var = "Dead_Mpred")
Dead_Mstar  <- estimatecolSum(trait_data_full_i = Ne_NT_tot, var = "Dead_Mstar")
Dead_Mtot   <- estimatecolSum(trait_data_full_i = Ne_NT_tot, var = "Dead_Mtot")
Mstar   <- estimatecolSum(trait_data_full_i = Ne_NT_tot, var = "Mstar")
Mpred       <- estimatecolSum(trait_data_full_i = Ne_NT_tot, var = "Mpred")
SSDgx       <- estimatecolSum(trait_data_full_i = Ne_NT_tot, var = "SSDgx")
SSDgxNum       <- estimatecolSum(trait_data_full_i = Ne_NT_tot, var = "SSDgx/Numerator_N_alpha")
survivorship_curve<- estimatecolSum(trait_data_full_i = Ne_NT_tot, var = "survivorship_curve")
dt_tf <- rbind(Dead_Mtot, rbind(Dead, rbind(Dead_F, rbind(Dead_Madd, rbind(Dead_Mfor,
                                                                           rbind(Dead_Mnat, rbind(Dead_Mpred, rbind(Dead_Mstar, rbind(rbind(Mpred, Mstar), 
                                                                                                                                      rbind(SSDgx, rbind(SSDgxNum, rbind(Deltasquare, 
                                                                                                                                                                         rbind(eggNumber_sc, survivorship_curve)))))))))))))



# center the values around the values of the Fmax=1
dt_tf_centered <- dt_tf %>%
  group_by(Species, variable, Scenario, replicate) %>%
  mutate(center = value - value[which(Fmax == 1)]) %>%
  ungroup()
dt_tf_centered$Fmax <- as.numeric(dt_tf_centered$Fmax)
dt_tf_centered$Scenario <- as.factor(dt_tf_centered$Scenario)
# dt_tf_centered <- dt_tf_centered[, -which(colnames(dt_tf_centered) == "value")]

# make slope estimate relative
# dt_tf_centered <- pivot_wider(dt_tf_centered, values_from = "center", names_from = "Scenario")
# dt_tf_centered$`With CC` <- dt_tf_centered$`With CC`/dt_tf_centered$`No CC`
# dt_tf_centered$`No CC` <- dt_tf_centered$`No CC`/dt_tf_centered$`No CC`
# dt_tf_centered <- pivot_longer(dt_tf_centered, cols=c("No CC", "With CC"), 
#                                     names_to = "Scenario", values_to = "center")

# diagnostics
ggplotdata <- dt_tf_centered
ggplotdata$logslope <- log(ggplotdata$center-min(ggplotdata$center, na.rm = T)+0.01)
ggplotdata$id <- 1:nrow(ggplotdata) 
ggplot(ggplotdata)+ 
  geom_density(aes(center))+ 
  facet_wrap(~variable, scales="free")
ggplot(ggplotdata, aes(sample=center))+ 
  stat_qq()+ 
  stat_qq_line()+ 
  facet_wrap(~variable, scales="free")

# ggplotdata$logslope <- log(ggplotdata$value-min(ggplotdata$value, na.rm = T)+0.01)



# for (sp in species){
#   model_Ne_Nt <- lm(data=dt_tf_centered[which(dt_tf_centered$Species==sp),], center~Scenario+Fmax+Scenario:Fmax)
#   coef <- coef(model_Ne_Nt)
#   coef["ScenarioWith CC:Fmax"]
# }

# if(!all(table(dt_tf_centered$Fmax, dt_tf_centered$Scenario, dt_tf_centered$Species)==140)) {
#   ID_pb   <- which(!table(dt_tf_centered$Fmax, dt_tf_centered$Scenario, dt_tf_centered$Species)==140, arr.ind = T)
#   nb_pb   <- 140-table(dt_tf_centered$Fmax, dt_tf_centered$Scenario, dt_tf_centered$Species)[ID_pb]
#   name_id <- dimnames(table(dt_tf_centered$Fmax, dt_tf_centered$Scenario, dt_tf_centered$Species))
#   dim_1 <- name_id[[1]][ID_pb[1]]
#   dim_2 <- name_id[[2]][ID_pb[2]]
#   dim_3 <- name_id[[3]][ID_pb[3]]
#   list = list()
#   for (k in 1:28) list[[k]] = data.frame("replicate"=1, "Scenario"=dim_2, "Fmax"=dim_1,"Species"=dim_3, "center"=0, "value"=0)
#   df <- do.call(rbind, list)
#   df <- cbind(df, "variable"=rep(names(table(dt_tf_centered$variable)), nb_pb/length(names(table(dt_tf_centered$variable)))))
#   dt_tf_centered <- rbind(dt_tf_centered, df)
# }

dt_tf_slope <- dt_tf_centered %>% 
  group_by(Species, variable) %>% 
  nest() %>% 
  mutate(model = map(data, ~lm(center ~ Fmax + Scenario + Fmax:Scenario, data = .x) %>% tidy)) %>% 
  unnest(model) %>% 
  # filter(term %in% c('Fmax:ScenarioWith CC', "(Intercept)")) %>%
  ungroup() 

# if slope Fmax decreasing : response of intercept and slope of noCC/wCC should be reversed
dt_tf_slope_abs <- pivot_wider(dt_tf_slope[, which(colnames(dt_tf_slope) %in% 
                                                             c("Species", "variable", "estimate", "term"))], 
                                   values_from = "estimate", names_from = "term")
dt_tf_slope_abs$`(Intercept)` = dt_tf_slope_abs$`(Intercept)`*dt_tf_slope_abs$Fmax
dt_tf_slope_abs$`Fmax:ScenarioWith CC` = dt_tf_slope_abs$`ScenarioWith CC`*sign(dt_tf_slope_abs$Fmax)
dt_tf_slope_abs <- pivot_longer(dt_tf_slope_abs, cols=c("Fmax:ScenarioWith CC", "(Intercept)"), values_to = "estimate", names_to="term")
dt_tf_slope <- right_join(dt_tf_slope[, which(colnames(dt_tf_slope) %in% c("Species", "variable","p.value", "term"))], 
                              dt_tf_slope_abs[, which(colnames(dt_tf_slope_abs) %in% 
                                                            c("Species", "variable", "estimate", "term"))])

#improve format table
dt_tf_slope <- dt_tf_slope[, c("Species", "estimate", "variable", "p.value", "term")] 
dt_tf_slope[,"p.adj.signif"] <- ifelse(dt_tf_slope$p.value < 0.05, "*", "")
colnames(dt_tf_slope)[which(colnames(dt_tf_slope)=="estimate")] <- "slope_difference"


# add Ne/NT response type when NO CC
dt_tf_slope$Response <- ifelse(dt_tf_slope$Species %in% c("Hake", "Whiting", "Saithe", "Sandeel", "AtlanticHerring"), "Increasing N\u2091/NT", 
                                   ifelse(dt_tf_slope$Species %in% c("CommonDab", "GreyGurnard", "EuropeanPlaice"), "Decreasing N\u2091/NT", 
                                          ifelse(dt_tf_slope$Species %in% c("NorwayPout", "CommonSole", "AtlanticCod", "EuropeanSprat", "Haddock", "HorseMackerel", 
                                                                                "AtlanticMackerel"), "Constant N\u2091/NT", "out"))) 

## PREPARE DATASET FOR PLOT ------------
# scale dataset
dt_tf_heatmap <- dt_tf_slope[, -which(colnames(dt_tf_slope)=="p.value")]
dt_tf_heatmap <- as.data.frame(pivot_wider(dt_tf_heatmap, values_from = "slope_difference", names_from = "variable"))
dt_tf_heatmap[, -which(colnames(dt_tf_heatmap) %in% c("Species", "Response", "p.adj.signif", "term"))] <- 
  scale(dt_tf_heatmap[, -which(colnames(dt_tf_heatmap) %in% c("Species", "Response", "p.adj.signif", "term"))])
dt_tf_heatmap <- pivot_longer(dt_tf_heatmap, cols= -c("Species", "Response", "p.adj.signif", "term"), values_to = "slope_difference", names_to = "variable")
dt_tf_heatmap$variable <- improve_vriable_names(dt_tf_heatmap$variable)
if (any(is.na(dt_tf_heatmap$slope_difference)))  dt_tf_heatmap <- dt_tf_heatmap[-which(is.na(dt_tf_heatmap$slope_difference)), ]
# order variables in the plot
dt_tf_heatmap$variable <- factor(dt_tf_heatmap$variable, levels = 
                                       c("N\u2091/NT", "N\u2091", "NT", "T/Vk", "T", "Vk.", "Nα"))
sp_to_order <- unique(as.matrix(dt_tf_heatmap[which(dt_tf_heatmap$variable == "N\u2091/NT"), "Species"]))
levelz = sp_to_order[order(as.matrix(dt_tf_heatmap[which(dt_tf_heatmap$variable == "N\u2091/NT"), "slope_difference"]))]
dt_tf_heatmap$Species = factor(dt_tf_heatmap$Species, levels= levelz) # ordering the species as in all other plots


## PLOT -------------
ptf1 <- ggplot(dt_tf_heatmap[which(dt_tf_heatmap$term=="Fmax:ScenarioWith CC"), ], aes(x=variable, y=Species))+
  geom_tile(aes(fill=slope_difference))+
  theme(
    axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5)) + 
  ggplot2::scale_fill_gradient2(
    low = "blue", 
    high = "red", 
    # oob = scales::squish,
    midpoint = mean(Ne_NT_tot_heatmap$slope_difference[which(Ne_NT_tot_heatmap$term=="Fmax:ScenarioWith CC")], na.rm=T), 
    limits = c(min(Ne_NT_tot_heatmap$slope_difference[which(Ne_NT_tot_heatmap$term=="Fmax:ScenarioWith CC")], na.rm=T),
               max(Ne_NT_tot_heatmap$slope_difference[which(Ne_NT_tot_heatmap$term=="Fmax:ScenarioWith CC")], na.rm=T))
    # breaks = c("Weaker slope with CC (-1)"=-1.0, "Unchanged (0)"=0, "Stronger slope with CC (>2)"=2.0)
  ) + 
  theme_classic()+
  facet_grid(rows = vars(Response), scales = "free", space = "free") +
  geom_text(aes(label=p.adj.signif), color="black", size=5) + 
  labs(title = "Scaled relative difference between without and with CC\nin slope of the species' functional traits'\nresponse to ↑ Fmax ",
       fill = "Difference in slope\nof the variable's\nresponse to ↑ Fmax") + 
  ylab("") + xlab("")

ggsave(filename = paste0(PATHOUT, "/Heatmaps/functinal_traits_slope.jpg"), plot = ptf1,
       width= 9, height = 6)

ptf2 <- ggplot(dt_tf_heatmap[which(dt_tf_heatmap$term=="(Intercept)"), ], aes(x=variable, y=Species))+
  geom_tile(aes(fill=slope_difference))+
  theme(
    axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5)) + 
  ggplot2::scale_fill_gradient2(
    low = "blue", 
    high = "red", 
    # oob = scales::squish,
    midpoint = mean(Ne_NT_tot_heatmap$slope_difference[which(Ne_NT_tot_heatmap$term=="(Intercept)")], na.rm=T), 
    limits = c(min(Ne_NT_tot_heatmap$slope_difference[which(Ne_NT_tot_heatmap$term=="(Intercept)")], na.rm=T),
               max(Ne_NT_tot_heatmap$slope_difference[which(Ne_NT_tot_heatmap$term=="(Intercept)")], na.rm=T))
    # breaks = c("Weaker slope with CC (-1)"=-1.0, "Unchanged (0)"=0, "Stronger slope with CC (>2)"=2.0)
  ) + 
  theme_classic()+
  facet_grid(rows = vars(Response), scales = "free", space = "free") +
  geom_text(aes(label=p.adj.signif), color="black", size=5) + 
  labs(title = "Scaled relative difference between without and with CC\nin intercept of the species' functional traits'\nresponse to ↑ Fmax ",
       fill = "Difference in intercept\nof the variable's\nresponse to ↑ Fmax") + 
  ylab("") + xlab("")

ggsave(filename = paste0(PATHOUT, "/Heatmaps/functinal_traits_intercept.jpg"), plot = ptf2,
       width= 9, height = 6)

# 
# ## PLOT -------------
# ptf <- ggplot(dt_tf_heatmap, aes(x=variable, y=Species))+
#   geom_tile(aes(fill=slope_difference_scaled))+
#   facet_grid(rows = vars(Response), scales = "free", space = "free") +
#   ggplot2::scale_fill_gradient2(
#     low = "blue", 
#     high = "red", 
#     oob = scales::squish,
#     midpoint = mean(dt_tf_heatmap$slope_difference_scaled, na.rm=T),
#     limits = c(min(dt_tf_heatmap$slope_difference_scaled, na.rm=T), 2),
#     breaks = c("Weaker slope with CC (-1)"=-1.0, "Unchanged (0)"=0, "Stronger slope with CC (>2)"=2.0)
#   ) + 
#   theme_classic()+
#   theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5)) + 
#   geom_text(aes(label=p.adj.signif), color="black", size=5) + 
#   labs(title ="Scaled relative difference between without and with CC\nin slope of the species' functional traits'\nresponse to ↑ Fmax ",
#        fill = "Difference in slope\nof the variable's\nresponse to ↑ Fmax\nwith or without CC") + 
#   ylab("") + xlab("") 
# 
# ggsave(filename = paste0(PATHOUT, "/Heatmaps/functinal_traits.jpg"), plot = ptf,
#        width= 9, height = 6)
# 
