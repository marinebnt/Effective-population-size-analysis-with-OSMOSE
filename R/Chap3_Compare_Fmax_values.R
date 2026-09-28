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
#*    2- Estimate the significativity of each variable 'value~Scenario' by species. 
#*       No dimensions removed
#*    3- Calculate the relative difference in response between the scenarios : (abs(with CC) - abs(no CC))/abs(no CC), by variable by species. 
#*       Dimensions left : Species, Variable, Replicate, Fmax.
#*    4- Average the difference by Fmax and replicate. 
#*       Dimensions left : Species, Variable. 


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
## ESTIMATE DIFF Fmax OF ALL VARIABLES -------------
Ne_NT_tot2 <- pivot_longer(Ne_NT_tot2, values_to = "value", names_to = "variable", cols = -c("Fmax", "replicate", "Species", "Scenario"))
#add significance level of the difference of slope with and without CC, between replicates
ggplotdata <- Ne_NT_tot2
ggplotdata$logvalue <- log(ggplotdata$value+min(ggplotdata$value, na.rm=T)+0.01)
ggplotdata$id <- 1:nrow(ggplotdata) 
ggplot(ggplotdata)+ 
  geom_density(aes(logvalue))+ 
  facet_wrap(~variable, scales="free")
ggplot(ggplotdata, aes(sample=logvalue))+ 
  stat_qq()+ 
  stat_qq_line()+ 
  facet_wrap(~variable, scales="free")
stat.test <- ggplotdata %>%
  group_by(Species, variable) %>%
  do(estimate = wilcox.test(logvalue~Scenario, data=.)) %>% #wilcox.test
  reframe(Species=Species, variable=variable, estimate = estimate$p.value)%>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance()
stat.test$p.adj.signif=ifelse(stat.test$estimate<0.05, "*", "")
stat.test$p.adj.signif <- ifelse(stat.test$p.adj.signif == "ns", "", stat.test$p.adj.signif) 
#estimate the difference between scenarios variable/species/Fmax/replicate
Ne_NT_tot2 <- pivot_wider(Ne_NT_tot2, values_from = "value", names_from = "Scenario")
Ne_NT_tot_diff <- Ne_NT_tot2 %>% 
  group_by(Species, replicate, variable, Fmax) %>% 
  mutate(difference=(`With CC` - `No CC`)/abs(`No CC`)) %>%
  ungroup() 
Ne_NT_tot_diff$difference = ifelse(Ne_NT_tot_diff[, "With CC"] <0, -1, 1) * Ne_NT_tot_diff$difference
# average of the difference between replicates & Fmax
Ne_NT_tot_diff <- Ne_NT_tot_diff %>%
  group_by(Species, variable) %>%
  summarize(mean_diff=mean(difference, na.rm=T)) %>%
  ungroup()
# add Ne/NT response type when NO CC
Ne_NT_tot_diff$Response <- ifelse(Ne_NT_tot_diff$Species %in% c("Hake", "Whiting", "Saithe", "Sandeel", "AtlanticHerring"), "Increasing N\u2091/NT", 
                                   ifelse(Ne_NT_tot_diff$Species %in% c("CommonDab", "GreyGurnard", "EuropeanPlaice"), "Decreasing N\u2091/NT", 
                                          ifelse(Ne_NT_tot_diff$Species %in% c("NorwayPout", "CommonSole", "AtlanticCod", "EuropeanSprat", "Haddock", "HorseMackerel", 
                                                                                "AtlanticMackerel"), "Constant N\u2091/NT", "out"))) 
# add p-value
Ne_NT_tot_diff <- full_join(Ne_NT_tot_diff, stat.test[, which(colnames(stat.test) %in% c("Species", "variable", "p.adj.signif"))])


## PREPARE DATASET FOR PLOT ------------
# scale dataset
Ne_NT_tot_heatmap <- Ne_NT_tot_diff
Ne_NT_tot_heatmap <- as.data.frame(pivot_wider(Ne_NT_tot_heatmap, values_from = "mean_diff", names_from = "variable"))
Ne_NT_tot_heatmap[, -which(colnames(Ne_NT_tot_heatmap) %in% c("Species", "Response", "p.adj.signif"))] <- 
  scale(Ne_NT_tot_heatmap[, -which(colnames(Ne_NT_tot_heatmap) %in% c("Species", "Response", "p.adj.signif"))])
Ne_NT_tot_heatmap <- pivot_longer(Ne_NT_tot_heatmap, cols= -c("Species", "Response", "p.adj.signif"), values_to = "mean_diff", names_to = "variable")
Ne_NT_tot_heatmap$variable <- improve_vriable_names(Ne_NT_tot_heatmap$variable)
if (any(is.na(Ne_NT_tot_heatmap$mean_diff))) Ne_NT_tot_heatmap <- Ne_NT_tot_heatmap[-which(is.na(Ne_NT_tot_heatmap$mean_diff)), ]
# order variables in the plot
Ne_NT_tot_heatmap$variable <- factor(Ne_NT_tot_heatmap$variable, levels = 
                                       c("N\u2091/NT", "N\u2091", "NT", "T/Vk", "T", "Vk.", "Nα"))
sp_to_order <- as.matrix(Ne_NT_tot_heatmap[which(Ne_NT_tot_heatmap$variable == "N\u2091/NT"), "Species"])
levelz = sp_to_order[order(as.matrix(Ne_NT_tot_heatmap[which(Ne_NT_tot_heatmap$variable == "N\u2091/NT"), "mean_diff"]))]
Ne_NT_tot_heatmap$Species = factor(Ne_NT_tot_heatmap$Species, levels= levelz) # ordering the species as in all other plots

## PLOT -------------
pNe <- ggplot(Ne_NT_tot_heatmap, aes(x=variable, y=Species))+
  geom_tile(aes(fill=mean_diff))+
  theme(
    axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5)) + 
  ggplot2::scale_fill_gradient2(
    low = "blue", 
    high = "red", 
    midpoint = mean(Ne_NT_tot_heatmap$mean_diff), 
    limits = c(min(Ne_NT_tot_heatmap$mean_diff), max(Ne_NT_tot_heatmap$mean_diff)),
    breaks = c("Lower values with CC (-2)"=-2.0, "Unchanged (0)"=0, "Higher values with CC (3)"=3.0)
  ) + 
  theme_classic()+
  facet_grid(rows = vars(Response), scales = "free", space = "free") +
  geom_text(aes(label=p.adj.signif), color="black", size=5) + 
  labs(title = "Scaled relative difference between without and with CC\nin of the species' functional traits'\nresponse to ↑ Fmax ",
       fill = "Difference \nof the variable's\nresponse to ↑ Fmax") + 
  ylab("") + xlab("")

ggsave(filename = paste0(PATHOUT, "/Heatmaps/effective_pop_size_Fmaxdiff.jpg"), plot = pNe,
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
Mstar       <- estimatecolSum(trait_data_full_i = Ne_NT_tot, var = "Mstar")
Mpred       <- estimatecolSum(trait_data_full_i = Ne_NT_tot, var = "Mpred")
SSDgx       <- estimatecolSum(trait_data_full_i = Ne_NT_tot, var = "SSDgx")
SSDgxNum       <- estimatecolSum(trait_data_full_i = Ne_NT_tot, var = "SSDgx/Numerator_N_alpha")
survivorship_curve<- estimatecolSum(trait_data_full_i = Ne_NT_tot, var = "survivorship_curve")
dt_tf <- rbind(Dead_Mtot, rbind(Dead, rbind(Dead_F, rbind(Dead_Madd, rbind(Dead_Mfor,
         rbind(Dead_Mnat, rbind(Dead_Mpred, rbind(Dead_Mstar, rbind(rbind(Mpred, Mstar), 
         rbind(SSDgx, rbind(SSDgxNum, rbind(Deltasquare, 
         rbind(eggNumber_sc, survivorship_curve)))))))))))))


#add significance level of the difference of slope with and without CC, between replicates
ggplotdata <- dt_tf
ggplotdata$logvalue <- log(ggplotdata$value-min(ggplotdata$value, na.rm=T)+0.1)
ggplotdata$id <- 1:nrow(ggplotdata) 
ggplot(ggplotdata)+ 
  geom_density(aes(logvalue))+ 
  facet_wrap(~variable, scales="free")
ggplot(ggplotdata, aes(sample=logvalue))+ 
  stat_qq()+ 
  stat_qq_line()+ 
  facet_wrap(~variable, scales="free")
stat.test <- ggplotdata %>%
  group_by(Species, variable) %>%
  do(estimate = wilcox.test(logvalue~Scenario, data=.)) %>% # wilcox.test
  reframe(Species=Species, variable=variable, estimate = estimate$p.value)%>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance()
stat.test$p.adj.signif=ifelse(stat.test$estimate<0.05, "*", "")
stat.test$p.adj.signif <- ifelse(stat.test$p.adj.signif == "ns", "", stat.test$p.adj.signif) 
#estimate the difference between scenarios variable/species/Fmax/replicate
dt_tf <- pivot_wider(dt_tf, values_from = "value", names_from = "Scenario")
dt_tf_diff <- dt_tf %>% 
  group_by(Species, replicate, variable, Fmax) %>% 
  mutate(difference=(`With CC` - `No CC`)/abs(`No CC`)) %>%
  ungroup() 
dt_tf_diff$difference = ifelse(dt_tf_diff[, "With CC"] <0, -1, 1) * dt_tf_diff$difference
dt_tf_diff <- dt_tf_diff[-which(is.infinite(dt_tf_diff$difference)),]
# average of the difference between replicates & Fmax
dt_tf_diff <- dt_tf_diff %>%
  group_by(Species, variable) %>%
  summarize(mean_diff=mean(difference, na.rm=T)) %>%
  ungroup()
# add Ne/NT response type when NO CC
dt_tf_diff$Response <- ifelse(dt_tf_diff$Species %in% c("Hake", "Whiting", "Saithe", "Sandeel", "AtlanticHerring"), "Increasing N\u2091/NT", 
                                  ifelse(dt_tf_diff$Species %in% c("CommonDab", "GreyGurnard", "EuropeanPlaice"), "Decreasing N\u2091/NT", 
                                         ifelse(dt_tf_diff$Species %in% c("NorwayPout", "CommonSole", "AtlanticCod", "EuropeanSprat", "Haddock", "HorseMackerel", 
                                                                              "AtlanticMackerel"), "Constant N\u2091/NT", "out"))) 
# add p-value
dt_tf_diff <- full_join(dt_tf_diff, stat.test[, which(colnames(stat.test) %in% c("Species", "variable", "p.adj.signif"))])


# ## PREPARE DATASET FOR PLOT ------------
# # scale dataset
# dt_tf_heatmap <- dt_tf_diff
# dt_tf_heatmap <- as.data.frame(pivot_wider(dt_tf_heatmap, values_from = "mean_diff", names_from = "Species"))
# dt_tf_heatmap[, -which(colnames(dt_tf_heatmap) %in% c("variable", "Response", "p.adj.signif"))] <- 
#   scale(dt_tf_heatmap[, -which(colnames(dt_tf_heatmap) %in% c("variable", "Response", "p.adj.signif"))])
# dt_tf_heatmap <- pivot_longer(dt_tf_heatmap, cols= -c("variable", "Response", "p.adj.signif"), values_to = "mean_diff", names_to = "Species")
# dt_tf_heatmap$variable <- improve_vriable_names(dt_tf_heatmap$variable)
# if (any(is.na(dt_tf_heatmap$mean_diff))) dt_tf_heatmap <- dt_tf_heatmap[-which(is.na(dt_tf_heatmap$mean_diff)), ]
# # order variables in the plot
# # dt_tf_heatmap$variable <- factor(dt_tf_heatmap$variable, levels = 
# #                                        c("N\u2091/NT", "N\u2091", "NT", "Nα", "T/Vk", "T", "Vk."))
# sp_to_order <- as.matrix(Ne_NT_tot_heatmap[which(Ne_NT_tot_heatmap$variable == "N\u2091/NT"), "Species"])
# levelz = sp_to_order[order(as.matrix(dt_tf_heatmap[which(Ne_NT_tot_heatmap$variable == "N\u2091/NT"), "mean_diff"]))]
# dt_tf_heatmap$Species = factor(dt_tf_heatmap$Species, levels= levelz) # ordering the species as in all other plots
# dt_tf_heatmap[dt_tf_heatmap$Species == "AtlanticCod",]
# 
# 
# 
# ## PLOT -------------
# ptf <- ggplot(dt_tf_heatmap, aes(x=variable, y=Species))+
#   geom_tile(aes(fill=mean_diff))+
#   facet_grid(rows = vars(Response), scales = "free", space = "free") +
#   ggplot2::scale_fill_gradient2(
#     low = "blue", 
#     high = "red", 
#     midpoint = median(dt_tf_heatmap$mean_diff, na.rm=T), 
#     limits = c(min(dt_tf_heatmap$mean_diff, na.rm=T), max(dt_tf_heatmap$mean_diff, na.rm=T))
#   ) + 
#   theme_classic()+
#   theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5)) + 
#   geom_text(aes(label=p.adj.signif), color="black", size=5) + 
#   labs(title ="Scaled difference between without and with CC\nin of the species' functional traits'\nresponse to ↑ Fmax ",
#        fill = "Difference in slope\nof the variable's\nresponse to ↑ Fmax\nwith or without CC") + 
#   ylab("") + xlab("") 
# 
# ggsave(filename = paste0(PATHOUT, "/Heatmaps/functinal_traits_diffFmax_Scaled_by_species.jpg"), plot = ptf,
#        width= 9, height = 6)



## PREPARE DATASET FOR PLOT ------------
# scale dataset
dt_tf_heatmap <- dt_tf_diff
dt_tf_heatmap <- as.data.frame(pivot_wider(dt_tf_heatmap, values_from = "mean_diff", names_from = "variable"))
dt_tf_heatmap[, -which(colnames(dt_tf_heatmap) %in% c("Species", "Response", "p.adj.signif"))] <- 
   scale(dt_tf_heatmap[, -which(colnames(dt_tf_heatmap) %in% c("Species", "Response", "p.adj.signif"))], center = T, scale = T)
dt_tf_heatmap <- pivot_longer(dt_tf_heatmap, cols= -c("Species", "Response", "p.adj.signif"), values_to = "mean_diff", names_to = "variable")
dt_tf_heatmap$variable <- improve_vriable_names(dt_tf_heatmap$variable)
if (any(is.na(dt_tf_heatmap$mean_diff))) dt_tf_heatmap <- dt_tf_heatmap[-which(is.na(dt_tf_heatmap$mean_diff)), ]
# order variables in the plot
# dt_tf_heatmap$variable <- factor(dt_tf_heatmap$variable, levels = 
#                                        c("N\u2091/NT", "N\u2091", "NT", "Nα", "T/Vk", "T", "Vk."))
sp_to_order <- as.matrix(Ne_NT_tot_heatmap[which(Ne_NT_tot_heatmap$variable == "N\u2091/NT"), "Species"])
levelz = sp_to_order[order(as.matrix(dt_tf_heatmap[which(Ne_NT_tot_heatmap$variable == "N\u2091/NT"), "mean_diff"]))]
dt_tf_heatmap$Species = factor(dt_tf_heatmap$Species, levels= levelz) # ordering the species as in all other plots
dt_tf_heatmap[dt_tf_heatmap$Species == "AtlanticCod",]

# add the log because too much difference in the dataset
# dt_tf_heatmap$mean_diff <- sign(dt_tf_heatmap$mean_diff) * log(abs(dt_tf_heatmap$mean_diff))



## PLOT -------------
ptf <- ggplot(dt_tf_heatmap, aes(x=variable, y=Species))+
  geom_tile(aes(fill=mean_diff))+
  facet_grid(rows = vars(Response), scales = "free", space = "free") +
  ggplot2::scale_fill_gradient2(
    low = "blue", 
    high = "red", 
    oob = scales::squish,
    midpoint = mean(dt_tf_heatmap$mean_diff, na.rm=T),
    limits = c(-1, 1),
    breaks = c("Lower values with CC (<-1)"=-1, "Unchanged (0)"=0, "Higher values with CC (>1)"=1.0), 
  ) + 
  theme_classic()+
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5)) + 
  geom_text(aes(label=p.adj.signif), color="black", size=5) + 
  labs(title ="Scaled and LOG relative difference between without and with CC\nin of the species' functional traits'\nresponse to ↑ Fmax ",
       fill = "Difference in \nthe variable's\nresponse to ↑ Fmax\nwith or without CC") + 
  ylab("") + xlab("") 

ggsave(filename = paste0(PATHOUT, "/Heatmaps/functinal_traits_diffFmax_Scaled_by_variable.jpg"), plot = ptf,
       width= 9, height = 6)
