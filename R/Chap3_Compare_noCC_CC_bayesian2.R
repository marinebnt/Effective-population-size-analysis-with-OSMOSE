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
library(mvabund)
library(tidyverse) # glm
library(patchwork) # plots glm
library(brms)

main = "New_gpu_conf.R" #config file
source(paste0(here::here(), "/R/functions/source.R")) #load config and key parameters and functions
source(paste0(here::here(), "/R/functions/functions_key_traits.R")) #plots against Fmax


# STEPS ------------

#* 

# PREPARE DATASETS of the two scenarios  ---------------

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


## SLIGHT ADDITION ON NE/NT DATA -----------
Ne_NT_tot2 <- Ne_NT_tot2[-which(Ne_NT_tot2$NT < 10), ]
Ne_NT_tot2$Ne_NT <- Ne_NT_tot2$Ne_N_alpha/Ne_NT_tot2$NT 
Ne_NT_tot2$'T/Vk' <- Ne_NT_tot2$generation_length/Ne_NT_tot2$Vk_point 
Ne_NT_tot2 <- Ne_NT_tot2[which(Ne_NT_tot2$Fmax == 1), ]

## REMOVING AGE DIMENSION ON TF --------------
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
dt_tf <- as.data.frame(rbind(Dead_Mtot, rbind(Dead, rbind(Dead_F, rbind(Dead_Madd, rbind(Dead_Mfor,
                                                                                         rbind(Dead_Mnat, rbind(Dead_Mpred, rbind(Dead_Mstar, rbind(rbind(Mpred, Mstar), 
                                                                                                                                                    rbind(SSDgx, rbind(SSDgxNum, rbind(Deltasquare, 
                                                                                                                                                                                       rbind(eggNumber_sc, survivorship_curve))))))))))))))
dt_tf <- dt_tf[which(dt_tf$Fmax == 1), ]


## BOXPLOTS
# Ne_NT_bx <- pivot_longer(Ne_NT_tot2, cols=-c("Fmax","replicate", "Species", "Scenario"), 
#                          values_to = "value", names_to = "variable")
# dfbx <- rbind(Ne_NT_bx, dt_tf[, -which(colnames(dt_tf) == "Fstart")])
# 
# for (x in names(table(dfbx$variable))){
#   list_sp_x <- list()
#   for (sp in species){
#     dfbxi <- dfbx[which(dfbx$Species == sp & dfbx$variable == x), ]
#     list_sp_x[[sp]] <- ggplot(dfbxi) +
#       geom_boxplot(aes(y=value, x=Scenario))+
#       ylab(x)+
#       labs(title=sp)
#   }
#   xarr <- ggarrange(plotlist=list_sp_x, ncol=5, nrow=3)
#   ggsave(plot = xarr, filename=paste0(PATHOUT, "/Clean_for_paper/", x, ".png"), create.dir = T,
#           height = 9, width=12)
# }



## ADDING THE NE/NT variable -----------

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
Ne_NT_tot_diff_tempo <- Ne_NT_tot_diff
Ne_NT_tot_diff_tempo$difference = ifelse(Ne_NT_tot_diff[, "With CC"] <0, -1, 1) * Ne_NT_tot_diff_tempo$difference
Ne_NT_tot_diff_tempo$variable = paste0(Ne_NT_tot_diff_tempo$variable, "_abs")
Ne_NT_tot_diff <- rbind(Ne_NT_tot_diff, Ne_NT_tot_diff_tempo)
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
Ne_NT_tot_diff$Response <- ifelse(Ne_NT_tot_diff$Species %in% c("NorwayPout", "CommonSole", "AtlanticCod", "EuropeanSprat", "Hake", "Whiting", "Saithe", "Sandeel", "AtlanticHerring"), "Increasing N\u2091/NT", 
                                  ifelse(Ne_NT_tot_diff$Species %in% c("CommonDab", "GreyGurnard", "EuropeanPlaice", "AtlanticMackerel"), "Decreasing N\u2091/NT", 
                                         ifelse(Ne_NT_tot_diff$Species %in% c("Haddock", "HorseMackerel"), "Constant N\u2091/NT", "out"))) 
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
# Ne_NT_tot_heatmap$variable <- factor(Ne_NT_tot_heatmap$variable, levels = 
# c("N\u2091/NT", "N\u2091", "NT", "T/Vk", "T", "Vk.", "Nα"))
sp_to_order <- as.matrix(Ne_NT_tot_heatmap[which(Ne_NT_tot_heatmap$variable == "N\u2091/NT"), "Species"])
levelz = sp_to_order[order(as.matrix(Ne_NT_tot_heatmap[which(Ne_NT_tot_heatmap$variable == "N\u2091/NT"), "mean_diff"]))]
Ne_NT_tot_heatmap$Species = factor(Ne_NT_tot_heatmap$Species, levels= levelz) # ordering the species as in all other plots

## ADDING THE TF VARIABLE ----------
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
dt_tf_diff_tempo <- dt_tf_diff
dt_tf_diff_tempo$difference = ifelse(dt_tf_diff[, "With CC"] <0, -1, 1) * dt_tf_diff$difference
dt_tf_diff_tempo$variable = paste0(dt_tf_diff_tempo$variable, "_abs")
dt_tf_diff <- rbind(dt_tf_diff, dt_tf_diff_tempo)
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



# Prepare dataset of OSMOSE parameters -------- 

# | THV | OSMOSE parameter or output variable |
# |-------------------------|-----------------------------------------------|
# | TL | ratio.min, ratio.max |
# | a_mat/L_mat | m0, m1 |
# | a_max | lifespan |
# | fecundity | egg.weight, maturity.r |
# | Growth coefficient | predation.ingestion, Topt, Tmin, Tmax, c_m |
# | Mortality | variables : starvation, predation, foraging, additional mortality rate |
# |  | parameter : mortality.additional.larva.rate |
# | Selectivity | centered |

## 1- Prepare dataset with the parameters from the configuration process (optimal temperature) ----------

# Topt dataset
df_Topt <- read.csv(paste0(here::here(), "/DATA/FromOSMOSE-nsConfiguration_outputs_global_complet.csv"))[, c("Species", "Final_Topt", 
                                                                                                             "Final_Tmax", "Final_Tmin")]
# names   <- rfishbase::species(stringr::str_replace(df_Topt$Species, "_", " "))
# write.csv(names, file=paste0(here::here(), "/DATA/names_species_rfishbase.csv"))

# match up scientific names to species names in OSMOSE
names <- read.csv(file=paste0(here::here(), "/DATA/names_species_rfishbase.csv"))
names_commun <- names[, c("Species", "FBname")]
colnames(df_Topt)[which(colnames(df_Topt)=="Species")] <- "Spe_scien"
colnames(names_commun)[which(colnames(names_commun)=="Species")] <- "Spe_scien"
colnames(names_commun)[which(colnames(names_commun)=="FBname")] <- "Species"
names_commun$Spe_scien <- stringr::str_replace(names_commun$Spe_scien, " ", "_")
tempo <- stringr::str_split(names_commun$Species, " ")
tempo <- sapply(1:length(tempo), function(x) {
  if (length(tempo[[x]])==1) return(tempo[[x]])
  if (length(tempo[[x]])==3) tempo[[x]] <- c(paste0(toupper(substr(tempo[[x]][2], 1, 1)), substr(tempo[[x]][2], 2, nchar(tempo[[x]][2]))), tempo[[x]][3])
  new_name <- paste0(tempo[[x]][1], toupper(substr(tempo[[x]][2], 1, 1)), substr(tempo[[x]][2], 2, nchar(tempo[[x]][2])))
  if(new_name == c("EuropeanHake")) new_name <- "Hake"
  return(new_name)
})
names_commun$Species <- tempo
names_commun$Species <- stringr::str_replace(names_commun$Species, " ", "")
df_Topt <- full_join(df_Topt, names_commun)
df_Topt[which(df_Topt$Spe_scien %in% c("Ammodytes_spp", "Crangon_crangon")), "Species"] <- c("Sandeel", "Shrimp")


## 2- Prepare dataset with the configuration parameters --------

# extract ony the parameters I need
nspecies <- osmose::get_par(conf, "simulation.nspecies") -1 #because we removed the Shrimp
df_param <-as.data.frame(sapply(c("ratio.max", "ratio.min",
                                  "m0", "m1",
                                  "lifespan",                             
                                  "egg.weight", "maturity.r",
                                  "predation.ingestion",
                                  "maint.energy.c_m",                             
                                  "mortality.additional.larva.rate",  
                                  #"mortality.additional.rate", 
                                  "species.name"
), 
FUN = function(x) {
  a <- osmose::get_par(conf, x)    
  if (length(a)>nspecies) a <- a[c(1:nspecies)]
  b <- names(a)
  a <- sapply(1:length(a), FUN=function(x) a [[x]] <- a[[x]][length(a[[x]])])
  names(a) <- b
  a <- unlist(a)
}))
df_param$spi <- unlist(stringr::str_split(rownames(df_param), "max."))[seq(2, 2*nrow(df_param), by=2)]
colnames(df_param)[which(colnames(df_param) == "species.name")] <- "Species"


## 3- Prepare dataset with the simulations outputs ------

PATHOUT  = paste0(here::here(), "/output/PLOTS/CHAP2/KEY_TRAITS_MATURE_removed0inNedata/")
if (!dir.exists(PATHOUT)) dir.create(PATHOUT, recursive = TRUE)
if (!dir.exists(paste0(PATHOUT, "/Clean_for_paper"))) dir.create(paste0(PATHOUT, "/Clean_for_paper"), recursive = TRUE)
pathoutputRDATA = paste0(PATHOUT, "RDATA")
base_path =paste0(here::here(), "/output/DATASETS/Chap2/Trait_variables_by_cohort_by_age_withWbyAge/")
#loading the dataset of traits averaged over cohorts Time1 to Time2, and reorganised to remove age=0, add reassign natural mortality as add mortality, and estimate natural mortality
data_complete_c <- load_trait_df(base_path = base_path, pathout=pathoutputRDATA)

# check Madd > Mfor
data_complete_cw <- pivot_wider(data_complete_c, values_from = "value", names_from="variable") 
data_complete_cw$Fmax <- as.factor(data_complete_cw$Fmax)
data_complete_cw$Species <- as.factor(data_complete_cw$Species)
data_complete_cw$mean_of_Madd_per_age <- as.numeric(data_complete_cw$mean_of_Madd_per_age)
data_complete_cw$mean_of_Mfor_per_age <- as.numeric(data_complete_cw$mean_of_Mfor_per_age)
ggplot(data_complete_cw, aes(col = Fmax)) +
  geom_smooth(aes(x = age, y = mean_of_Madd_per_age, linetype = "Madd")) +
  geom_smooth(aes(x = age, y = mean_of_Mfor_per_age, linetype = "Mfor")) + 
  geom_smooth(aes(x = age, y = mean_of_Mpred_per_age, linetype = "Mpred")) + 
  scale_linetype_manual(values = c("Madd" = 1, "Mfor" = 2, "Mpred" = 3)) +
  facet_wrap(~Species)

#combine Mpred, Madd, Mfor and Mstarv into Mnat
df_M <- data_complete_c[data_complete_c$variable %in% 
                          c("mean_of_Madd_per_age", 
                            "mean_of_Mfor_per_age", 
                            "mean_of_Mpred_per_age", 
                            "mean_of_Mstar_per_age"), ]
df_M <- df_M %>%  #Mnat
  group_by(Species, age, replicate, Fmax, Fstart) %>%
  summarize(Mnat = sum(value))
df_M <- df_M %>%  #Mnat cumulated between ages
  group_by(Species, replicate, Fmax, Fstart) %>%
  summarize(Mnat = sum(Mnat))
df_M <- df_M %>%  #Mnat cumulated between ages averaged between Fmax, Fstart and replicate
  group_by(Species) %>%
  summarize(Mnat = mean(Mnat))


## 4- Prepare dataset of the selectivity curve ------------
path_csv <- "C:/Users/mbeneat/Documents/osmose/Chap2/osmose-ns/osmose-ns-v4.4.1_clim/fishing_rate"
df_select <- as.data.frame(matrix(-99, ncol=2, nrow=length(list.files(path_csv, pattern = "csv", full.names = T))))
colnames(df_select) <- c("selectq90", "spi")
k=0
for (csv in list.files(path_csv, pattern = "csv", full.names = T)){
  k=k+1
  f <- read.csv(csv)
  f_1 <- f[1, -1]
  f_1 <- f_1[-which(f_1 == 0)]
  q90 <- as.numeric(quantile(as.numeric(f_1), probs = 0.9))
  df_select$selectq90[k] = as.numeric(f_1[which(as.numeric(f_1)>=(q90))[1]])
  df_select$spi[k] <- stringr::str_remove(unlist(stringr::str_split(csv, pattern = "file."))[2], pattern = ".csv")
}


## 5- Combine datasets to shape a unique dataset : df_select + df_Topt + df_param + df_M ---------
df_t <- full_join(df_param, df_select)
df_t <- full_join(df_t, df_Topt, by="Species") 
df_t <- full_join(df_t, df_M, by="Species") 
df_t <- na.omit(df_t) 
rownames(df_t) <- df_t$Species 
# df_t <- df_t[, -which(colnames(df_t) %in% c("Spe_scien", "Species", "spi"))] 
df <- sapply(1:ncol(df_t), function(x) {
  if (is.numeric(df_t[1,x])) as.numeric(df_t[,x])
  else(df_t[,x])
}) 
colnames(df) <- colnames(df_t) 
rownames(df) <- rownames(df_t) 
df <- as.data.frame(df)



# MULTIVARITE GLM ----------------

## MERGE DATASETS + stat plots------
Ne <- Ne_NT_tot2
tf <- dt_tf
data <- full_join(Ne, tf)
data <- data[, -which(colnames(data)=="Fstart")]
data <- pivot_longer(data, values_to = "value", names_to = "Scenario", cols=c("No CC", "With CC"))

ggplotdata <- data
# ggplotdata$logvalue <- log(ggplotdata$value+min(ggplotdata$value, na.rm=T)+0.01)
ggplotdata$id <- 1:nrow(ggplotdata) 
ggplot(ggplotdata)+ 
  geom_density(aes(value))+ 
  facet_wrap(~variable, scales="free")
ggplot(ggplotdata, aes(sample=value))+ 
  stat_qq()+ 
  stat_qq_line()+ 
  facet_wrap(~variable, scales="free")
stat.test <- ggplotdata %>%
  group_by(Species, variable) %>%
  do(estimate = wilcox.test(value~Scenario, data=.)) %>% #wilcox.test
  reframe(Species=Species, variable=variable, estimate = estimate$p.value)%>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance()

data <- pivot_wider(data, values_from = "value", names_from = "variable")
data <- full_join(data, df)
data <- data[, -which(colnames(data) %in% c("Spe_scien", "spi"))] 
data[,-c(1:4)] <- sapply(5:ncol(data), FUN = function(x) as.numeric(as.matrix(data[,x])))



## MULTIVARIATE BAYESIAN MODEL -----------------------------

data$NT <- data$NT*10^-10 # *************
data$Ne_N_alpha <- data$NT*10^-10 # ***************

### VKpoint ---------

# DV = colnames(data)[5:10]
bform1 <- 
  bf(mvbind(Vk_point, eggNumber_sc, SSDgx, Deltasquare, Dead) ~ Species * Scenario)+
  brmsfamily("Gamma")
set_rescor(TRUE)

fiit1 <- brm(bform1, data = data, chains = 2, cores = 2, control = list(adapt_delta = 0.95))
save(fiit1, file=paste0(pathoutputRDATA, "/multivariate_bay_model_Vk_egg.RData"))
load(file=paste0(pathoutputRDATA, "/multivariate_bay_model_Vk_egg.RData"))
#fiit1 <- add_criterion(fiit1, "loo")
summary(fiit1)
# VarCorr(fiit1) # will not work with gamma model
# plot(fiit1) # super long
# fiit1$prior
# fiit1$fit

conditional_effects(fiit1, "Species:Scenario", resp = "Vkpoint")
conditional_effects(fiit1, "Species:Scenario", resp = "eggNumbersc")
conditional_effects(fiit1, "Scenario", resp = "NeNT")


pp_check(fiit1, resp = c("Vkpoint"))
bayes_R2(fiit1)
posterior_summary(fiit1)
as_draws_df(fiit1) |> 
  glimpse()

### generation length ---------

# DV = colnames(data)[5:10]
bform1 <- 
  bf(mvbind(generation_length, survivorship_curve, Mpred, NT) ~ Species * Scenario)+
  brmsfamily("Gamma")
set_rescor(TRUE)

# fiit2 <- brm(bform1, data = data, chains = 2, cores = 2, control = list(adapt_delta = 0.95))
# save(fiit2, file=paste0(pathoutputRDATA, "/multivariate_bay_model_T_mortality.RData"))
load(file=paste0(pathoutputRDATA, "/multivariate_bay_model_T_mortality.RData"))
#fiit1 <- add_criterion(fiit1, "loo")
summary(fiit2)
# VarCorr(fiit1) # will not work with gamma model
# plot(fiit1) # super long
# fiit1$prior
# fiit1$fit

conditional_effects(fiit2, "Species:Scenario", resp = "generationlength")
conditional_effects(fiit2, "Species:Scenario", resp = "Mpred")
conditional_effects(fiit2, "Scenario", resp = "NeNT")


pp_check(fiit2, resp = c("Vkpoint"))
bayes_R2(fiit2)
posterior_summary(fiit2)
as_draws_df(fiit2) |> 
  glimpse()


## Projection # https://mjskay.github.io/tidybayes/reference/add_predicted_draws.html
library(tidybayes)
library(boxplotdbl)
NewData<-expand.grid(Species=unique(data$Species),
                     Scenario=unique(data$Scenario))


## PLOTTING # https://mjskay.github.io/tidybayes/articles/tidybayes.html
# m the model 

data_pred <- NewData %>%
  add_predicted_draws(  ndraws = 100, fiit1)
data_pred <- pivot_wider(data_pred, names_from = ".category", values_from = ".prediction")

par(mfrow=c(1,2))
boxplotdou(Vkpoint~Species, 
           data_pred[which(data_pred$Scenario=="No CC"),], 
           log(SSDgx)~Species, 
           data_pred[which(data_pred$Scenario=="No CC"),], 
           main="No Climate Change", name.on.axis=FALSE)
boxplotdou(Vkpoint~Species, 
           data_pred[which(data_pred$Scenario=="With CC"),], 
           log(SSDgx)~Species, 
           data_pred[which(data_pred$Scenario=="With CC"),], 
           main = "With Climate Change", name.on.axis=FALSE)

par(mfrow=c(1,2))
boxplotdou(Vkpoint~Species, 
           data_pred[which(data_pred$Scenario=="No CC"),], 
           log(Dead)~Species, 
           data_pred[which(data_pred$Scenario=="No CC"),], 
           main="No Climate Change", name.on.axis=FALSE)
boxplotdou(Vkpoint~Species, 
           data_pred[which(data_pred$Scenario=="With CC"),], 
           log(Dead)~Species, 
           data_pred[which(data_pred$Scenario=="With CC"),], 
           main = "With Climate Change", name.on.axis=FALSE)


par(mfrow=c(1,2))
boxplotdou(log(Deltasquare)~Species, 
           data_pred[which(data_pred$Scenario=="No CC"),], 
           Vkpoint~Species, 
           data_pred[which(data_pred$Scenario=="No CC"),], 
           main="No Climate Change", name.on.axis=FALSE)
boxplotdou(log(Deltasquare)~Species, 
           data_pred[which(data_pred$Scenario=="With CC"),], 
           Vkpoint~Species, 
           data_pred[which(data_pred$Scenario=="With CC"),], 
           main = "With Climate Change", name.on.axis=FALSE)



ggplot(data=data_pred, aes(y = Vkpoint, x = Scenario)) +
  stat_halfeye(.width = c(.95, .66), position = position_nudge(y = -0.3)) 

ggplot(data=data_pred, aes(y = Vkpoint, x = Scenario)) +
  stat_halfeye(.width = c(.95, .66), position = position_nudge(y = -0.3)) +
  geom_point() 
# stat_lineribbon() + #aes(width = c(.99, .95, .8, .5)), alpha = 0.25) +
# geom_point() +
# scale_fill_brewer(palette = "Greys")

 

## GGPLOT with crossbar

NewData<-expand.grid(Species=unique(data$Species),
                     Scenario=unique(data$Scenario))

data_pred <- NewData %>%
  add_predicted_draws(ndraws = 10, fiit1)
data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE))
    c(Q9 = as.numeric(quantile(x[[col]], 0.9)),
      Q1 = as.numeric(quantile(x[[col]], 0.1)))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}
data_pred <- full_join(data_pred, data_summary(data_pred, ".prediction", c("Species", "Scenario", ".category")))
data_pred <- pivot_wider(data_pred, names_from = ".category", values_from = c(".prediction", "Q9", "Q1"))
data_pred <- full_join(data_pred, data[, 
                                       c("Species", "Scenario", "Final_Topt", "Final_Tmin", "Final_Tmax")])


data_pred$Species_n <- NA
data_pred$Species_n[c(
  which(data_pred$Species %in% unique(data_pred$Species)[1])[1], 
  which(data_pred$Species %in% unique(data_pred$Species)[2])[1],
  which(data_pred$Species %in% unique(data_pred$Species)[3])[1],
  which(data_pred$Species %in% unique(data_pred$Species)[4])[1],
  which(data_pred$Species %in% unique(data_pred$Species)[5])[1],
  which(data_pred$Species %in% unique(data_pred$Species)[6])[1],
  which(data_pred$Species %in% unique(data_pred$Species)[7])[1],
  which(data_pred$Species %in% unique(data_pred$Species)[8])[1],
  which(data_pred$Species %in% unique(data_pred$Species)[9])[1],
  which(data_pred$Species %in% unique(data_pred$Species)[10])[1],
  which(data_pred$Species %in% unique(data_pred$Species)[11])[1],
  which(data_pred$Species %in% unique(data_pred$Species)[12])[1],
  which(data_pred$Species %in% unique(data_pred$Species)[13])[1],
  which(data_pred$Species %in% unique(data_pred$Species)[14])[1],
  which(data_pred$Species %in% unique(data_pred$Species)[15])[1]
)] <- data_pred$Species[c(
  which(data_pred$Species %in% unique(data_pred$Species)[1])[1], 
  which(data_pred$Species %in% unique(data_pred$Species)[2])[1],
  which(data_pred$Species %in% unique(data_pred$Species)[3])[1],
  which(data_pred$Species %in% unique(data_pred$Species)[4])[1],
  which(data_pred$Species %in% unique(data_pred$Species)[5])[1],
  which(data_pred$Species %in% unique(data_pred$Species)[6])[1],
  which(data_pred$Species %in% unique(data_pred$Species)[7])[1],
  which(data_pred$Species %in% unique(data_pred$Species)[8])[1],
  which(data_pred$Species %in% unique(data_pred$Species)[9])[1],
  which(data_pred$Species %in% unique(data_pred$Species)[10])[1],
  which(data_pred$Species %in% unique(data_pred$Species)[11])[1],
  which(data_pred$Species %in% unique(data_pred$Species)[12])[1],
  which(data_pred$Species %in% unique(data_pred$Species)[13])[1],
  which(data_pred$Species %in% unique(data_pred$Species)[14])[1],
  which(data_pred$Species %in% unique(data_pred$Species)[15])[1]
)]

data_pred$log.prediction_eggNumbersc<- log(data_pred$.prediction_eggNumbersc )
data_pred$log.prediction_Deltasquare<- log(data_pred$.prediction_Deltasquare )
data_pred$log.prediction_Dead<- log(data_pred$.prediction_Dead)
data_pred$logQ1_eggNumbersc<- log(data_pred$Q1_eggNumbersc )
data_pred$logQ9_Deltasquare<- log(data_pred$Q9_Deltasquare )
data_pred$logQ9_eggNumbersc<- log(data_pred$Q9_eggNumbersc )
data_pred$logQ9_Dead<- log(data_pred$Q9_Dead )
data_pred$logQ1_Deltasquare<- log(data_pred$Q1_Deltasquare )
data_pred$logQ1_Dead<- log(data_pred$Q1_Dead )

ggplot(data_pred, aes_string(y=paste0(".prediction_", "Vkpoint"), 
                             x=paste0("log.prediction_", "Deltasquare"), 
                            group="Species", col="Final_Tmin")) +
 geom_point() + 
 facet_wrap(~Scenario) + 
 geom_crossbar(aes_string(ymin=paste0("Q1_", "Vkpoint"), 
                          ymax=paste0("Q9_", "Vkpoint"), 
                          xmin=paste0("logQ1_", "Deltasquare"), 
                          xmax=paste0("logQ9_", "Deltasquare")), 
               position = 'dodge', 
               width=0.1, 
               alpha=0.9, 
               size=1) + 
 geom_text(aes(label=Species_n), 
           position=position_dodge(width=0.9), 
           vjust=-0.25) +
 scale_colour_continuous(palette="YlOrRd")
  

ggplot(data_pred, aes_string(y=paste0(".prediction_", "Vkpoint"), 
                             x=paste0("log.prediction_", "Dead"), 
                             group="Species", col="Final_Tmin")) +
  geom_point() + 
  facet_wrap(~Scenario) + 
  geom_crossbar(aes_string(ymin=paste0("Q1_", "Vkpoint"), 
                           ymax=paste0("Q9_", "Vkpoint"), 
                           xmin=paste0("logQ1_", "Dead"), 
                           xmax=paste0("logQ9_", "Dead")), 
                position = 'dodge', 
                width=0.1, 
                alpha=0.9, 
                size=1) + 
  geom_text(aes(label=Species_n), 
            position=position_dodge(width=0.9), 
            vjust=-0.25) +
  scale_colour_continuous(palette="YlOrRd")



data_pred <- NewData %>%
 add_predicted_draws(  ndraws = 100, fiit2)
data_pred <- pivot_wider(data_pred, names_from = ".category", values_from = ".prediction")

par(mfrow=c(1,2))
boxplotdou(survivorshipcurve ~Species, 
          data_pred[which(data_pred$Scenario=="No CC"),], 
          generationlength ~Species, 
          data_pred[which(data_pred$Scenario=="No CC"),], 
          main="No Climate Change", name.on.axis=FALSE)
boxplotdou(survivorshipcurve~Species, 
          data_pred[which(data_pred$Scenario=="With CC"),], 
          generationlength~Species, 
          data_pred[which(data_pred$Scenario=="With CC"),], 
          main = "With Climate Change", name.on.axis=FALSE)
par(mfrow=c(1,2))
boxplotdou(Mpred~Species, 
          data_pred[which(data_pred$Scenario=="No CC"),], 
          generationlength~Species, 
          data_pred[which(data_pred$Scenario=="No CC"),], 
          main="No Climate Change", name.on.axis=FALSE)
boxplotdou(Mpred~Species, 
          data_pred[which(data_pred$Scenario=="With CC"),], 
          generationlength~Species, 
          data_pred[which(data_pred$Scenario=="With CC"),], 
          main = "With Climate Change", name.on.axis=FALSE)

par(mfrow=c(1,2))
boxplotdou(NT~Species, 
          data_pred[which(data_pred$Scenario=="No CC"),], 
          generationlength~Species, 
          data_pred[which(data_pred$Scenario=="No CC"),], 
          main="No Climate Change", name.on.axis=FALSE)
boxplotdou(NT~Species, 
          data_pred[which(data_pred$Scenario=="With CC"),], 
          generationlength~Species, 
          data_pred[which(data_pred$Scenario=="With CC"),], 
          main = "With Climate Change", name.on.axis=FALSE)


ggplot(data=data_pred, aes(y = Vkpoint, x = Scenario)) +
 stat_halfeye(.width = c(.95, .66), position = position_nudge(y = -0.3)) 

ggplot(data=data_pred, aes(y = Vkpoint, x = Scenario)) +
 stat_halfeye(.width = c(.95, .66), position = position_nudge(y = -0.3))
#+
geom_point() # stat_lineribbon() + #aes(width = c(.99, .95, .8, .5)), alpha = 0.25) +
# geom_point() +
# scale_fill_brewer(palette = "Greys")


# OTHER THINGS I TRIED --------------------
### TRY AGAIN WITH BAYESIAN  + STUDENT BUT NOT CONVINCING ----------
# 
# bform2 <- 
#   bf(mvbind(Numerator_N_alpha, generation_length, Ne_N_alpha, Vk_point, NT, Ne_NT) ~ Species * Scenario)+
#   brmsfamily("student") + 
#   set_rescor(TRUE)
# 
# datasc<- data
# datasc[, -c(1:4)] <- scale(datasc[, -c(1:4)])
# mean(datasc$Mpred)
# var(datasc$Mpred)
# fit2 <- brm(bform2, data = datasc, chains = 2, cores = 2, control = list(adapt_delta = 0.95))
# 
# require(stats)
# x <- matrix(1:10, ncol = 2)
# (centered.x <- scale(x, scale = FALSE))
# cov(centered.scaled.x <- scale(x)) # all 1
# 
# # Messages d'avis :
# # 1: There were 2000 transitions after warmup that exceeded the maximum treedepth. Increase max_treedepth above 10. See
# # https://mc-stan.org/misc/warnings.html#maximum-treedepth-exceeded 
# # 2: There were 1 chains where the estimated Bayesian Fraction of Missing Information was low. See
# # https://mc-stan.org/misc/warnings.html#bfmi-low 
# # 3: Examine the pairs() plot to diagnose sampling problems
# #  
# # 4: The largest R-hat is NA, indicating chains have not mixed.
# # Running the chains for more iterations may help. See
# # https://mc-stan.org/misc/warnings.html#r-hat 
# # 5: Bulk Effective Samples Size (ESS) is too low, indicating posterior means and medians may be unreliable.
# # Running the chains for more iterations may help. See
# # https://mc-stan.org/misc/warnings.html#bulk-ess 
# # 6: Tail Effective Samples Size (ESS) is too low, indicating posterior variances and tail quantiles may be unreliable.
# # Running the chains for more iterations may help. See
# # https://mc-stan.org/misc/warnings.html#tail-ess
# 
# fit2 <- add_criterion(fit2, "loo")
# summary(fit2)
# VarCorr(fit2)
# # plot(fit2) # super long
# # fiit1$prior
# # fiit1$fit
# 
# conditional_effects(fit2, "Species:Scenario", resp = "NeNT")
# conditional_effects(fit2, "Species:Scenario", resp = "Vkpoint")
# conditional_effects(fit2, "Scenario", resp = "NeNT")
# 
# pp_check(fit2, resp = c("Vkpoint"))
# pp_check(fit2, resp = c("NumeratorNalpha"))
# pp_check(fit2, resp = c("NeNalpha"))
# pp_check(fit2, resp = c("NT"))
# pp_check(fit2, resp = c("NeNT"))
# bayes_R2(fit2)
# posterior_summary(fit2)
# as_draws_df(fit2) |> 
#   glimpse()
# 
### Example if two distributions among the response variables ----------
# bf_normal <- bf(Numerator_N_apha ~ sex + (1|p|fosternest) + (1|q|dam))+ gaussian()
# bf_back   <- bf( ~ hatchdate + (1|p|fosternest) + (1|q|dam))+ gamma()
# fit2 <- brm(bf_tarsus + bf_back + set_rescor(TRUE), 
#             data = BTdata, chains = 2, cores = 2) 
# 
# compare two fits
# loo(fiit1, fit2)
# 
### MULTIVARIATE GLM -----------------------------
# NOT CONCLUSIVE - DATA is GAMMA while manyglm with gamma is not working
#
# TRY WITH GAUSSIAN AND INVERSE GAUSSIAN :failing, not successful
# IV = c("Species", "Scenario")
# DV = colnames(data)[5:9]
# data_spp <- mvabund(data[, DV])
# # data_spp <- (scale(data_spp))
# 
# glm1<-manylm(data_spp~Species*Scenario,data=data, family = "gaussian")
# glm1<-manylm(data_spp~Species*Scenario,data=data, family = "inverse.gaussian")
# # summary(glm1)
# 
# Fitted=data.frame(fitted(glm1)) %>% 
#   rownames_to_column("ID") %>% 
#   pivot_longer(-ID,names_to = "Spec",values_to = "Fitted")
# 
# Residuals=data.frame(residuals(glm1))%>% 
#   rownames_to_column("ID") %>% 
#   pivot_longer(-ID,names_to = "Spec",values_to = "Residuals")
# 
# ModelOutputs<-Fitted %>% 
#   left_join(Residuals,by=c("ID","Spec"))
# p3<-ggplot(ModelOutputs)+
#   geom_point(aes(x=Fitted,y=Residuals,colour=Spec))+
#   theme_classic()+
#   scale_color_viridis_d()+
#   labs(y="Residuals",x="Fitted Values")
# p4<-ggplot(ModelOutputs) +
#   stat_qq(aes(sample=Residuals,colour=Spec))+
#   stat_qq_line(aes(sample=Residuals,colour=Spec))+
#   theme_classic()+
#   scale_color_viridis_d()+
#   facet_wrap(~Spec)+
#   theme(legend.position = "none")+
#   labs(y="Sample Quartiles",x="Theoretical Quartiles")
# p3+p4+plot_layout(guides = "collect")
# 
# 
# # these gamma models are all failing : issue in manyglm function for GAMMA distribution
# # glm2 <- manyglm(formula=data_spp~Species*Scenario,data=data, family= "gamma")
# # glm2 <- manyglm(data_spp~Species*Scenario,data=data, family= Gamma(link="log"))
# # glm2 <- manyany("glm", data_spp, data = data, data_spp ~ Species+Scenario, 
# #                 family = tweedie(var.power = 2, link.power = 0), var.power = 2)
# 
# Fitted=data.frame(fitted(glm2)) %>% 
#   rownames_to_column("ID") %>% 
#   pivot_longer(-ID,names_to = "Spec",values_to = "Fitted")
# Residuals=data.frame(residuals(glm2))%>% 
#   rownames_to_column("ID") %>% 
#   pivot_longer(-ID,names_to = "Spec",values_to = "Residuals")
# ModelOutputs<-Fitted %>% 
#   left_join(Residuals,by=c("ID","Spec"))
# p3<-ggplot(ModelOutputs)+
#   geom_point(aes(x=Fitted,y=Residuals,colour=Spec))+
#   theme_classic()+
#   scale_color_viridis_d()+
#   labs(y="Residuals",x="Fitted Values")
# p4<-ggplot(ModelOutputs) +
#   stat_qq(aes(sample=Residuals,colour=Spec))+
#   stat_qq_line(aes(sample=Residuals,colour=Spec))+
#   theme_classic()+
#   scale_color_viridis_d()+
#   facet_wrap(~Spec)+
#   theme(legend.position = "none")+
#   labs(y="Sample Quartiles",x="Theoretical Quartiles")
# p3+p4+plot_layout(guides = "collect")

