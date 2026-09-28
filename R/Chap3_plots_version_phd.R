library(boot)
library(tidyverse)
library(caret)
library(mda)
library(klaR)
library(ggrepel)
library(pls)
library(factoextra)
library(plsVarSel)
library(ggpubr)
library(rstatix)
library(hrbrthemes)
library(stargazer)
library(gsubfn)
require(grid)

main = "New_gpu_conf.R" #config file
source(paste0(here::here(), "/R/functions/source.R")) #load config and key parameters and functions
source(paste0(here::here(), "/R/functions/functions_plots.R")) #plots against Fmax

#* 
#* 
#* warning : all what was done with `DOcrit` is wrong, because what I thought was the species-specific critical dissolved oxygen parameter 
#*           had nothing to do with it. The calculation of O2_dist is also wrong.    
#*           
#*            

# FUNCTIONS -----------------

#### DATASETS clean ------------------------
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
    dataset <- pivot_longer(dataset_t, cols = -c("age", "replicate", "Fmax", "Fstart", 'Species'), names_to = "variable", values_to = "value")
  }
  data_t  <- dataset[-which(dataset$variable %in% var_less_rep), ]
  dataset <- rbind(data_t, na.omit(dataset[which(dataset$variable %in% var_less_rep), ]))
  if (any(as.matrix(dataset_t[,grep(colnames(dataset_t), pattern="abundanceByAge")])<100, na.rm=T)){
    dataset_t[which(as.matrix(dataset_t[,grep(colnames(dataset_t), pattern="abundanceByAge")])<100), -which(colnames(dataset_t) %in%
                                                                                                              c("age", "replicate", "Fmax", "Fstart", 'Species'))] <- NA
    dataset <- pivot_longer(dataset_t, cols = -c("age", "replicate", "Fmax", "Fstart", 'Species'), names_to = "variable", values_to = "value")
  }
  data_t  <- dataset[-which(dataset$variable %in% var_less_rep), ]
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

# PREPARE DATASETS of the two scenarios  ---------------

Ne_NT_tot = NULL
Ne_NT_tot2 = NULL

## LOOP TO EXTRACT DATA BY SCENARIO ---------
for (Scenario in c("Chap2", "Chap3")){
  cat("\nScenario: ", Scenario)
  
  # prepare environment
  Ne_NT = NULL 
  Ne_NT2 = NULL 
  cat("=======> ", Scenario, "\n")
  PATHOUT  = paste0("output/PLOTS/", toupper(Scenario), "/KEY_TRAITS_MATURE_removed0inNedata/")
  if (!dir.exists(PATHOUT)) dir.create(PATHOUT, recursive = TRUE)
  if (!dir.exists(paste0(PATHOUT, "/Clean_for_paper"))) dir.create(paste0(PATHOUT, "/Clean_for_paper"), recursive = TRUE)
  pathoutputRDATA = paste0(PATHOUT, "RDATA")
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
    
    final_dataset_age <- when_variables_missing_age_dim(final_dataset_age, add_missing_ages = F)
    final_dataset_age <- when_null_or_low_abund(final_dataset_age)
    final_dataset_age <- what_ages_have_null_dead_nb(final_dataset_age)
    
    ## Extract only the dimensions I will be using 
    NeNT <- final_dataset_age[which(final_dataset_age$variable %in% c("Deltasquare", "Mpred", "Mstar", "eggNumber_sc", 
                                                                      "Dead", "Dead_F", "Dead_Madd", "Dead_Mfor", "Dead_Mnat", "Dead_Mpred", 
                                                                      "Dead_Mstar", "Dead_Mtot", "survivorship_curve", "SSDgx")), 
                              -which(colnames(final_dataset_age) %in% c("Fstart"))] %>% ungroup()     
    NeNT2 <- final_dataset_age[which(final_dataset_age$variable %in% c("Ne_N_alpha", "NT", "generation_length", "Vk_point", "Numerator_N_alpha")), 
                               -which(colnames(final_dataset_age) %in% c("Fstart", "age"))] %>% ungroup()    
    
    ## Increase the dimension of the population size at maturation age (Numerator_N_alpha) (add age dimension)
    ## to be able to manipulate it with other variables => estimate SSDgx/Nalpha
    numerator_n_alpha <- NULL
    for (x in min(as.numeric(as.matrix(NeNT$age))):30) {
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
    
    ## Cleaning environment
    rm(final_dataset_age) 
    NeNT[,"Scenario"] <- ifelse(Scenario == "Chap2", "No CC", "With CC") # *********
    NeNT2[,"Scenario"] <- ifelse(Scenario == "Chap2", "No CC", "With CC") # *********
    
    ## Binding the Species together
    if (is.null(Ne_NT)) {
      Ne_NT <- pivot_wider(NeNT, values_from = "value", names_from = "variable") 
      Ne_NT2 <- pivot_wider(NeNT2, values_from = "value", names_from = "variable") 
    } else {
      Ne_NT <- rbind(Ne_NT, pivot_wider(NeNT, values_from = "value", names_from = "variable")) 
      Ne_NT2 <- rbind(Ne_NT2, pivot_wider(NeNT2, values_from = "value", names_from = "variable")) 
    }
  }#sp
  
  ## Binding the Scenarios together
  if (is.null(Ne_NT_tot)) {Ne_NT_tot <- Ne_NT 
  } else Ne_NT_tot <- rbind(Ne_NT_tot, Ne_NT)   
  if (is.null(Ne_NT_tot2)) {Ne_NT_tot2 <- Ne_NT2 
  } else Ne_NT_tot2 <- rbind(Ne_NT_tot2, Ne_NT2) 
}#scenario

## SLIGHT ADDITION ON NE/NT DATA -----------
if(any(which(Ne_NT_tot2$NT < 10))) Ne_NT_tot2 <- Ne_NT_tot2[-which(Ne_NT_tot2$NT < 10), ]
Ne_NT_tot2$Ne_NT <- Ne_NT_tot2$Ne_N_alpha/Ne_NT_tot2$NT 
Ne_NT_tot2$'T/Vk' <- Ne_NT_tot2$generation_length/Ne_NT_tot2$Vk_point 

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

# KEEPING ONLY THE FMAX=1 DIMENSION ----------------------
Ne_NT_tot2 <- Ne_NT_tot2[which(Ne_NT_tot2$Fmax == 1), ]
dt_tf <- dt_tf[which(dt_tf$Fmax == 1), ]



# TABLE of the relative change of abundance with CC and no CC-----------------

#comparing with and without climate change biomass
biomassplot = NULL
rel_biomassplot = NULL
list_relative_biomass <- list()
NTplot = list()
biomassplot_curve = list()
abundanceplot_curve = list()
for (nm in "data_Ne"){
  for (sp in species){
    
    d_i_i <- Ne_NT_tot2
    id_biom <- c(which(colnames(d_i_i) %in% c("NT", "Species", "Fmax", "Scenario", "replicate")))
    d_i_i <- d_i_i[, id_biom]
    d_i_i <- d_i_i[which(d_i_i$Species == sp), ]
    d_i_i <- d_i_i[which(d_i_i$Fmax == 1), ]
    d_i_i$Iszero <- ifelse(d_i_i$NT < 100, "Abundance <100", "Abundance >100")
    
    NTplot[[sp]] <- ggplot(d_i_i, aes(y=NT, x=Scenario)) +
      geom_point(aes(colour = Iszero))+
      labs(title = "Abundance per fishing pressure", subtitle = sp)
    
    ref_data_percent <- d_i_i %>%
      filter(Scenario == "No CC") %>%
      dplyr::select(replicate, NT) %>%
      dplyr::rename(ref_value_percent_abund = NT)
    # Subtract reference from all means
    rel_d_i_i_perc <- d_i_i %>%
      left_join(ref_data_percent, by = c("replicate"), relationship = "many-to-one") %>%
      mutate(relative_abundance_percent = 100*NT/ref_value_percent_abund-100) %>%
      dplyr::select(-ref_value_percent_abund)
    list_relative_biomass[[sp]] <- rel_d_i_i_perc
    rel_biomassplot[[sp]] <- ggplot(rel_d_i_i_perc, aes(y=relative_abundance_percent, x=Scenario)) +
      geom_point(aes(colour = Iszero))+
      labs(title = "Relative biomass per fishing pressure", subtitle = sp)
  }
}

## Estimate the % in decline (in biomass and abundance, at first age at reproduction) from Fmax=1.5 to Fmax = 0.5
df_relative_N <- do.call(rbind, list_relative_biomass)
summary_df_relative_abunance <- df_relative_N %>%
  group_by(Scenario, Species)%>%
  summarise(across(relative_abundance_percent, list(mean=mean, sd=sd), na.rm=TRUE))
# summarize_each(percent_difference_withFmax0.5 = mean(relative_abundance_percent, na.rm=T))
summary_df_relative_abunance <- summary_df_relative_abunance %>%
  dplyr::mutate(percent_difference_withScenario=round(relative_abundance_percent_mean, 2)) %>%
  dplyr::mutate(percent_difference_withScenario_sd=round(relative_abundance_percent_sd, 2)) %>%
  select(Scenario, Species, percent_difference_withScenario, percent_difference_withScenario_sd)
summary_df_relative_abunance <- pivot_wider(summary_df_relative_abunance, 
                                            values_from = c("percent_difference_withScenario", "percent_difference_withScenario_sd"), 
                                            names_from = "Species")
write.csv(summary_df_relative_abunance, file=paste0(PATHOUT, "/Clean_for_paper/Abund_declinetrends.csv"))


biomassplot = NULL
rel_biomassplot = NULL
list_relative_biomass <- list()
NTplot = list()
biomassplot_curve = list()
abundanceplot_curve = list()
for (nm in "data_Ne"){
  for (sp in species){
    
    d_i_i <- Ne_NT_tot2
    id_biom <- c(which(colnames(d_i_i) %in% c("Ne_N_alpha", "Species", "Fmax", "Scenario", "replicate")))
    d_i_i <- d_i_i[, id_biom]
    d_i_i <- d_i_i[which(d_i_i$Species == sp), ]
    d_i_i <- d_i_i[which(d_i_i$Fmax == 1), ]
    d_i_i$Iszero <- ifelse(d_i_i$Ne_N_alpha < 100, "Abundance <100", "Abundance >100")
    
    NTplot[[sp]] <- ggplot(d_i_i, aes(y=Ne_N_alpha, x=Scenario)) +
      geom_point(aes(colour = Iszero))+
      labs(title = "Abundance per fishing pressure", subtitle = sp)
    
    ref_data_percent <- d_i_i %>%
      filter(Scenario == "No CC") %>%
      dplyr::select(replicate, Ne_N_alpha) %>%
      dplyr::rename(ref_value_percent_abund = Ne_N_alpha)
    # Subtract reference from all means
    rel_d_i_i_perc <- d_i_i %>%
      left_join(ref_data_percent, by = c("replicate"), relationship = "many-to-one") %>%
      mutate(relative_abundance_percent = 100*Ne_N_alpha/ref_value_percent_abund-100) %>%
      dplyr::select(-ref_value_percent_abund)
    list_relative_biomass[[sp]] <- rel_d_i_i_perc
    rel_biomassplot[[sp]] <- ggplot(rel_d_i_i_perc, aes(y=relative_abundance_percent, x=Scenario)) +
      geom_point(aes(colour = Iszero))+
      labs(title = "Relative biomass per fishing pressure", subtitle = sp)
  }
}

## Estimate the % in decline (in biomass and abundance, at first age at reproduction) from Fmax=1.5 to Fmax = 0.5
df_relative_N <- do.call(rbind, list_relative_biomass)
summary_df_relative_abunance <- df_relative_N %>%
  group_by(Scenario, Species)%>%
  summarise(across(relative_abundance_percent, list(mean=mean, sd=sd), na.rm=TRUE))
# summarize_each(percent_difference_withFmax0.5 = mean(relative_abundance_percent, na.rm=T))
summary_df_relative_abunance <- summary_df_relative_abunance %>%
  dplyr::mutate(percent_difference_withScenario=round(relative_abundance_percent_mean, 2)) %>%
  dplyr::mutate(percent_difference_withScenario_sd=round(relative_abundance_percent_sd, 2)) %>%
  select(Scenario, Species, percent_difference_withScenario, percent_difference_withScenario_sd)
summary_df_relative_abunance <- pivot_wider(summary_df_relative_abunance, 
                                            values_from = c("percent_difference_withScenario", "percent_difference_withScenario_sd"), 
                                            names_from = "Species")
write.csv(summary_df_relative_abunance, file=paste0(PATHOUT, "/Clean_for_paper/Ne_declinetrends.csv"))

# PLOTS -----------------


## BOXPLOT ---------------

Ne_NT_bx <- pivot_longer(Ne_NT_tot2, cols=-c("Fmax","replicate", "Species", "Scenario"),
                         values_to = "value", names_to = "variable")
dfbx <- rbind(Ne_NT_bx, dt_tf[, -which(colnames(dt_tf) == "Fstart")])

for (x in names(table(dfbx$variable))){
  list_sp_x <- list()
  for (sp in species){
    dfbxi <- dfbx[which(dfbx$Species == sp & dfbx$variable == x), ]
    list_sp_x[[sp]] <- ggplot(dfbxi) +
      geom_boxplot(aes(y=value, x=Scenario))+
      ylab(x)+
      labs(title=sp)
  }
  xarr <- ggarrange(plotlist=list_sp_x, ncol=5, nrow=3)
  ggsave(plot = xarr, filename=paste0(PATHOUT, "/BOXPLOT/", stringr::str_replace(x, "/", "_"), ".png"), create.dir = T,
         height = 9, width=12)
}


## Lollipop plot of the difference -----------------

# dataset for lollipop plot
dfbx <- pivot_wider(dfbx, names_from = "Scenario", values_from = "value")

data_barplot <- dfbx
data_barplot <- data_barplot[which(data_barplot$variable %in% c("Ne_NT", "NT", "Ne_N_alpha")),]
data_barplot$Species <- as.factor(data_barplot$Species)
data_barplot$variable <- improve_vriable_names(data_barplot$variable)
data_barplot2 <- data_barplot

# estimate the difference
data_barplot2$diff = (data_barplot2$`With CC` - data_barplot2$`No CC`)
data_barplot2$diff_rel = (data_barplot2$`With CC` - data_barplot2$`No CC`)/abs(data_barplot2$`No CC`)
# reduce the scale of the largest variables
data_barplot2[which(data_barplot2$variable %in% c("NT", "Nₑ")),which(colnames(data_barplot2)=="diff")] <-
  c(10^(-10)) * (data_barplot2[which(data_barplot2$variable %in% c("NT", "Nₑ")),which(colnames(data_barplot2)=="diff")])
# average of the difference by variable 
data_barplot2 <- data_barplot2 %>% 
  group_by(Species, variable) %>% 
  summarize(across(c(`diff`, `diff_rel`), 
                   list(sd = ~sd(., na.rm = T), mymean = ~mean(., na.rm = T)))) %>% 
  ungroup
# order dataset
data_barplot2$variable <- factor(data_barplot2$variable, levels = c("NT", "Nₑ", "Nₑ/NT"))
data_barplot2$Species <- factor(data_barplot2$Species, 
                                levels = as.matrix(data_barplot2[which(data_barplot2$variable=="Nₑ/NT"), "Species"])[order(as.matrix(data_barplot2[which(data_barplot2$variable=="Nₑ/NT"), "diff_mymean"]))])

# plot
lollipop <- ggplot(data_barplot2, aes(x=Species, y=diff_rel_mymean))  + 
  geom_hline(yintercept = 0, color = "darkgrey")+
  geom_segment(aes(x=Species, xend=Species, y=0, yend=diff_rel_mymean), color="grey") +
  geom_point(color=ifelse(round(data_barplot2$diff_rel_mymean, digits=2)<0,  ## *************
                          "brown2", ifelse(round(data_barplot2$diff_rel_mymean, digits=2)==0, "black", "#69b3a2")), size=4) +
  coord_flip() + 
  theme_bw() +
  theme(
    panel.grid.major.x = element_blank(),
    panel.border = element_blank(),
    axis.ticks.x = element_blank()
  ) +
  xlab("") +
  ylab("")+
  facet_wrap(~variable, scale="free_x", ncol=5)+
  scale_y_continuous(n.breaks = 8)

ggsave(lollipop, file=paste0(PATHOUT, "/LOLLI/lollipopNe_rel.png"), 
       height = 4.29, width=10.3, create.dir = T)

data_barplot2[,"variable"] <- as.character(as.matrix(data_barplot2[,"variable"]))
data_barplot2[which(data_barplot2$variable %in% c("NT", "Nₑ")),"variable"] <-
  paste0(as.matrix(data_barplot2[which(data_barplot2$variable %in% c("NT", "Nₑ")),"variable"]), " (e-10)")
data_barplot2$variable <- factor(data_barplot2$variable, levels = c("NT (e-10)", "Nₑ (e-10)", "Nₑ/NT"))
lollipop <- ggplot(data_barplot2, aes(x=Species, y=diff_mymean))  + 
  geom_hline(yintercept = 0, color = "darkgrey")+
  geom_segment(aes(x=Species, xend=Species, y=0, yend=diff_mymean), color="grey") +
  geom_point(color=ifelse(round(data_barplot2$diff_mymean, digits=2)<0,  ## *************
                          "brown2", ifelse(round(data_barplot2$diff_mymean, digits=2)==0, "black", "#69b3a2")), size=4) +
  coord_flip() + 
  theme_bw() +
  theme(
    panel.grid.major.x = element_blank(),
    panel.border = element_blank(),
    axis.ticks.x = element_blank()
  ) +
  xlab("") +
  ylab("")+
  facet_wrap(~variable, scale="free_x", ncol=5)+
  scale_y_continuous(n.breaks = 8)

ggsave(lollipop, file=paste0(PATHOUT, "/LOLLI/lollipopNe.png"), 
       height = 4.29, width=10.3, create.dir = T)


# dataset for lollipop plot : second set of variables
data_barplot <- dfbx
data_barplot <- data_barplot[which(data_barplot$variable %in% c("Numerator_N_alpha", "generation_length", "Vk_point", "T/Vk")),]
data_barplot$Species <- as.factor(data_barplot$Species)
data_barplot$variable <- improve_vriable_names(data_barplot$variable)
data_barplot3 <- data_barplot

# estimate the difference
data_barplot3$diff = (data_barplot3$`With CC` - data_barplot3$`No CC`)
data_barplot3$diff_rel = (data_barplot3$`With CC` - data_barplot3$`No CC`)/abs(data_barplot3$`No CC`)
# reduce the scale of the largest variables
data_barplot3[which(data_barplot3$variable %in% c("Nα")),which(colnames(data_barplot3)=="diff")] <-
  c(10^(-10)) * (data_barplot3[which(data_barplot3$variable %in% c("Nα")),which(colnames(data_barplot3)=="diff")])
# average of the difference by variable 
data_barplot3 <- data_barplot3 %>% 
  group_by(Species, variable) %>% 
  summarize(across(c(`diff`, `diff_rel`), 
                   list(sd = ~sd(., na.rm = T), mymean = ~mean(., na.rm = T)))) %>% 
  ungroup
# order dataset
data_barplot3$variable <- factor(data_barplot3$variable, levels = c("T", "Vk.", "T/Vk.", "Nα"))
data_barplot3$Species <- factor(data_barplot3$Species, 
                                levels = as.matrix(data_barplot2[which(data_barplot2$variable=="Nₑ/NT"), "Species"])[order(as.matrix(data_barplot2[which(data_barplot2$variable=="Nₑ/NT"), "diff_mymean"]))])

# plot
lollipop <- ggplot(data_barplot3, aes(x=Species, y=diff_rel_mymean))  + 
  geom_hline(yintercept = 0, color = "darkgrey")+
  geom_segment(aes(x=Species, xend=Species, y=0, yend=diff_rel_mymean), color="grey") +
  geom_point(color=ifelse(round(data_barplot3$diff_rel_mymean, digits=2)<0,  ## *************
                          "brown2", ifelse(round(data_barplot3$diff_rel_mymean, digits=2)==0, "black", "#69b3a2")), size=4) +
  coord_flip() + 
  theme_bw() +
  theme(
    panel.grid.major.x = element_blank(),
    panel.border = element_blank(),
    axis.ticks.x = element_blank()
  ) +
  xlab("") +
  ylab("")+
  facet_wrap(~variable, scale="free_x", ncol=5)+
  scale_y_continuous(n.breaks = 8)

ggsave(lollipop, file=paste0(PATHOUT, "/LOLLI/lollipopVk_rel.png"), 
       height = 4.29, width=10.3, create.dir = T)

data_barplot3[,"variable"] <- as.character(as.matrix(data_barplot3[,"variable"]))
data_barplot3[which(data_barplot3$variable %in% c("Nα")),"variable"] <-
  paste0(as.matrix(data_barplot3[which(data_barplot3$variable %in% c("Nα")),"variable"]), " (e-10)")
data_barplot3$variable <- factor(data_barplot3$variable, levels = c("T", "Vk.", "T/Vk.", "Nα (e-10)"))
lollipop <- ggplot(data_barplot3, aes(x=Species, y=diff_mymean))  + 
  geom_hline(yintercept = 0, color = "darkgrey")+
  geom_segment(aes(x=Species, xend=Species, y=0, yend=diff_mymean), color="grey") +
  geom_point(color=ifelse(round(data_barplot3$diff_mymean, digits=2)<0,  ## *************
                          "brown2", ifelse(round(data_barplot3$diff_mymean, digits=2)==0, "black", "#69b3a2")), size=4) +
  coord_flip() + 
  theme_bw() +
  theme(
    panel.grid.major.x = element_blank(),
    panel.border = element_blank(),
    axis.ticks.x = element_blank()
  ) +
  xlab("") +
  ylab("")+
  facet_wrap(~variable, scale="free_x", ncol=5)+
  scale_y_continuous(n.breaks = 8)

ggsave(lollipop, file=paste0(PATHOUT, "/LOLLI/lollipopVk.png"), 
       height = 4.29, width=10.3, create.dir = T)


## Regression and associated plot -------------

dfbx2 <- rbind(Ne_NT_bx, dt_tf[, -which(colnames(dt_tf) == "Fstart")])
dfbx2 <- dfbx2[which(dfbx2$variable %in% c("Ne_N_alpha", "Ne_NT", "NT", 
                                           "Numerator_N_alpha", "Vk_point",
                                           "generation_length",
                                           "T/Vk")), ]
dfbx2 <- pivot_wider(dfbx2, values_from = "value", names_from = "variable")

## add the temperature and oxygen data
# environmental temperature and oxygen concentration : estimated before (see "Compare_Tmin_CC_NoCC.R")
matrix_mean_T_sp  <- get(load(file=paste0(pathoutputRDATA, "/Tmean_by_species.RData")))
matrix_mean_T_sp$env <- "mean_temp"
matrix_mean_O2_sp <- get(load(file=paste0(pathoutputRDATA, "/Oxygen_by_species.RData")))
matrix_mean_O2_sp$env <- "mean_O2"
matrix_mean_T_sp<- rbind(matrix_mean_T_sp, matrix_mean_O2_sp)
# Topt dataset (estimated by Alaia while creating config) & DO crit (from configuration)
df_Topt   <- read.csv(paste0(here::here(), "/DATA/FromOSMOSE-nsConfiguration_outputs_global_complet.csv"))[, c("Species", "Final_Topt", "Final_Tmax", "Final_Tmin")]
df_Oxygen <- as.data.frame(unlist(get_par(conf=conf, c("oxygen.c2"))))
colnames(df_Oxygen)  <- "oxygen.c2" #Je pensais que ce param correspondait à DO crit mais en fait pas du tout, c'est pour scaler je crois
df_Oxygen[, "sp_id"] <- stringr::str_remove(rownames(df_Oxygen), pattern="species.oxygen.c2.sp")
name_sp   <- as.data.frame(unlist(get_par(conf=conf, c("species.name"))))
colnames(name_sp)  <- "Species" 
name_sp[, "sp_id"]   <- stringr::str_remove(rownames(name_sp), pattern="species.name.sp")
df_Oxygen <- left_join(df_Oxygen, name_sp)
df_Oxygen <- df_Oxygen[, which(colnames(df_Oxygen) %in% c("oxygen.c2", "Species"))]
rownames(df_Oxygen) <- NULL
# match up scientific names to species names in OSMOSE (I have saved the data, but can be rerun)
# names   <- rfishbase::species(stringr::str_replace(df_Topt$Species, "_", " "))
# write.csv(names, file=paste0(here::here(), "/DATA/names_species_rfishbase.csv"))
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
# attach temperature and oxygen data 
df_Topt <- full_join(df_Topt, df_Oxygen)
# add the mean T° encountered by each species and calculate : 
# the distance of env T° to Tmean, divided by the range of T° => (Tmean-Topt)/(Tmax-Tmin) 
df_Topt <- full_join(df_Topt, matrix_mean_T_sp)
df_Topt <- pivot_longer(df_Topt, cols = c("No CC", "With CC"), values_to = "mean_Temp_env", names_to = "Scenario")
# estimate the temperature distance : how far is the species from its optimal temperature ? and weight that by the temperature range of the sp
df_Topt$Topt_dist <- (df_Topt$mean_Temp_env - df_Topt$Final_Topt)/(df_Topt$Final_Tmax-df_Topt$Final_Tmin)
df_Topt$Tmin_dist <- (df_Topt$mean_Temp_env - df_Topt$Final_Tmin)/(df_Topt$Final_Tmax-df_Topt$Final_Tmin)
df_Topt$Tmax_dist <- (df_Topt$mean_Temp_env - df_Topt$Final_Tmax)/(df_Topt$Final_Tmax-df_Topt$Final_Tmin)
df_Topt$DOcrit    <- (df_Topt$mean_Temp_env - df_Topt$oxygen.c2)

# add Temperature distance to the dataset
dfbx2 <- full_join(df_Topt, dfbx2)
dfbx2 <- dfbx2[, -which(colnames(dfbx2) == c("Spe_scien"))]
cor(df_Topt[,c("Final_Tmax", "Final_Topt", "Final_Tmin")])

# add dimensions to the dataset: ratio, and difference of With CC // Without CC (absolute, relative and scaled)
dfbx2 <- pivot_longer(dfbx2, cols=-c("Species", "Scenario", "replicate", "env"), names_to = "variable", values_to = "value")
dfbx2 <- pivot_wider(dfbx2, names_from = "Scenario", values_from = "value")
dfbx2$ratio <- dfbx2$`With CC`/c(dfbx2$`No CC` + 0.0001)
dfbx2$diff_abs <- dfbx2$`With CC` - dfbx2$`No CC`+0.0001
dfbx2$diff_rel <- (dfbx2$`With CC` - dfbx2$`No CC`+0.0001)/abs(dfbx2$`No CC` +0.0001)
dfbx2 <- dfbx2 %>%
  group_by(variable, env) %>%
  mutate(diff_abs_scaled = scale(diff_abs, )) %>%
  ungroup()
dfbx2 <- dfbx2 %>%
  group_by(variable, env) %>%
  mutate(diff_rel_scaled = scale(diff_rel)) %>%
  ungroup()
dfbx2$variable[-which(dfbx2$variable == "Tmax_dist")] <- improve_vriable_names(dfbx2$variable[-which(dfbx2$variable == "Tmax_dist")])

# test data before regression (should be tested for all the dimensions (absolute, relative, scaled))
ggplotdata <- dfbx2
ggplotdata$value <- log(ggplotdata$diff_abs_scaled)
ggplotdata$id <- 1:nrow(ggplotdata) 
ggplot(ggplotdata)+ 
  geom_density(aes(value))+ 
  facet_wrap(~variable+env, scales="free")
ggplot(ggplotdata, aes(sample=value))+ 
  stat_qq()+ 
  stat_qq_line()+ 
  facet_wrap(~variable+env, scales="free")

# estimate regressions parameters + make it pretty with stargazer
table_regressions <- list()
table_intercept_full_dataset <- data.frame(matrix(NA, nrow=4, ncol=7))
table_slope_full_dataset     <- data.frame(matrix(NA, nrow=4, ncol=7))
colnames(table_slope_full_dataset) <- c("Nₑ/NT", "NT", "Nₑ", "T", "Vk.", "T/Vk.", "Nα")
colnames(table_intercept_full_dataset) <- c("Nₑ/NT", "NT", "Nₑ", "T", "Vk.", "T/Vk.", "Nα")
dim_focus <- c("diff_abs", "diff_rel", "diff_abs_scaled", "diff_rel_scaled")
rownames(table_slope_full_dataset) <- dim_focus
rownames(table_intercept_full_dataset) <- dim_focus
iteraction_list <- list()
iteraction_list_plot <- list()
for (temp_dist in c("Topt_dist", "Tmin_dist", "Tmax_dist")){
  for (diff_focus in dim_focus){ # loop to compare the regressions between the different ways to average calculate the difference
    
    data_reg <- dfbx2[which(dfbx2$env == "mean_temp"), ]
    if(temp_dist == "DOcrit")    data_reg <- dfbx2[which(dfbx2$env == "mean_O2"), ]
    data_reg[which(data_reg$variable %in% c(temp_dist)), diff_focus] <- data_reg$ratio[which(data_reg$variable %in% c(temp_dist))]
    data_reg <- pivot_wider(
      data_reg[which(data_reg$variable %in% c("Nₑ/NT", "NT", "Nₑ",  "T", "Vk.", "T/Vk.", "Nα", temp_dist)), which(colnames(data_reg) %in% 
                                            c("Species", "replicate", "variable", "value", diff_focus))], 
      values_from = diff_focus, 
      names_from = "variable")
    colnames(data_reg)[-which(colnames(data_reg) %in% c(temp_dist, "replicate", "Species"))] <- 
      paste0(diff_focus, colnames(data_reg)[-which(colnames(data_reg)  %in% c(temp_dist, "replicate", "Species"))])
    
    # Created nested tibble
    nested_df <- tibble(formula = c(paste0("`", diff_focus, "Nₑ/NT` ~", temp_dist),
                                    paste0("`", diff_focus, "NT` ~", temp_dist),
                                    paste0("`", diff_focus, "Nₑ` ~", temp_dist), 
                                    paste0("`", diff_focus, "T` ~", temp_dist), 
                                    paste0("`", diff_focus, "Vk.` ~", temp_dist),
                                    paste0("`", diff_focus, "T/Vk.` ~", temp_dist),
                                    paste0("`", diff_focus, "Nα` ~", temp_dist))) %>%
      group_by(ID = formula) %>%
      group_modify(~ as_tibble(data_reg)) %>%
      nest() 
                  
    # Get model estimates
    estimates <- nested_df %>%
      mutate(estimates = data %>% map2(ID, ~ tidy(lm(.y, data = .x)))) %>%
      select(-data) %>%
      unnest()
    
    id_intercept <- which(estimates$term=="(Intercept)")
    id_slope <- which(estimates$term==temp_dist)
    ID_intercept <- as.matrix(estimates[id_intercept, "ID"])
    ID_slope <- as.matrix(estimates[id_slope, "ID"])
    order_estimates <- c(grep(pattern=paste0(diff_focus, colnames(table_slope_full_dataset)[1], "` ~"), ID_intercept),
                         grep(pattern=paste0(diff_focus, colnames(table_slope_full_dataset)[2], "` ~"), ID_intercept),
                         grep(pattern=paste0(diff_focus, colnames(table_slope_full_dataset)[3], "` ~"), ID_intercept),
                         grep(pattern=paste0(diff_focus, colnames(table_slope_full_dataset)[4], "` ~"), ID_intercept),
                         grep(pattern=paste0(diff_focus, colnames(table_slope_full_dataset)[5], "` ~"), ID_intercept),
                         grep(pattern=paste0(diff_focus, colnames(table_slope_full_dataset)[6], "` ~"), ID_intercept),
                         grep(pattern=paste0(diff_focus, colnames(table_slope_full_dataset)[7], "` ~"), ID_intercept))
    table_intercept_full_dataset[diff_focus, ] <- as.matrix(estimates[id_intercept, "estimate"])[order_estimates]
    table_slope_full_dataset[diff_focus, ] <- as.matrix(estimates[id_slope, "estimate"])[order_estimates]
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
    table_regressions[[temp_dist]][[diff_focus]] <- gsubfn("([0-9.]+)", ~replace_numbers(x), star)
  }
}

# plot the biplots and the regression 
for (temp_dist in c("Topt_dist", "Tmin_dist", "Tmax_dist")){
  for (diff_focus in dim_focus){ # loop to compare the regressions between the different ways to average calculate the difference
    
    data_reg <- dfbx2[which(dfbx2$env == "mean_temp"), ]
    if(temp_dist == "DOcrit")    data_reg <- dfbx2[which(dfbx2$env == "mean_O2"), ]
    data_reg[which(data_reg$variable %in% c(temp_dist)), diff_focus] <- data_reg$ratio[which(data_reg$variable %in% c(temp_dist))]
    data_reg <- pivot_wider(
      data_reg[which(data_reg$variable %in% c("Nₑ/NT", "NT", "Nₑ",  "T", "Vk.", "T/Vk.", "Nα", temp_dist)), which(colnames(data_reg) %in% 
                                                                                                                      c("Species", "replicate", "variable", "value", diff_focus))], 
      values_from = diff_focus, 
      names_from = "variable")
    colnames(data_reg)[-which(colnames(data_reg) %in% c(temp_dist, "Species"))] <- 
      paste0(diff_focus, colnames(data_reg)[-which(colnames(data_reg)  %in% c(temp_dist, "Species"))])
    
    
    
    # data biplot
    data_biplo <- dfbx2[which(dfbx2$env == "mean_temp"), ]
    if(temp_dist == "DOcrit")    data_biplo <- dfbx2[which(dfbx2$env == "mean_O2"), ]
    data_biplo[which(data_biplo$variable %in% c(temp_dist)), diff_focus] <- data_biplo$ratio[which(data_biplo$variable %in% c(temp_dist))] # replace the diff of the temperature by the ratio : that's the indicator we want to focus on 
    data_biplo <- pivot_wider(
      data_biplo[which(data_biplo$variable %in% c("Nₑ/NT", "NT", "Nₑ",  "T", "Vk.", "T/Vk.", "Nα", temp_dist)), which(colnames(data_biplo) %in% 
                                                                                                                      c("Species", "replicate", "variable", "value", diff_focus))], 
      values_from = diff_focus, 
      names_from = "variable")
    
    
    # transform dataset for 
    data_biplo <- data_biplo %>% 
      group_by(Species) %>% 
      summarize(across(c(`Nₑ/NT`,`NT`, `Nₑ`, `T`, `Vk.`, `T/Vk.`, `Nα`, temp_dist), 
                       list(sd = sd, mymean = mean))) %>% 
      ungroup
    colnames(data_biplo) <- stringr::str_remove(colnames(data_biplo), pattern = "_mymean")
    data_biplo <- na.omit(data_biplo)
    
    
    # PLOT 
    p1 <- ggplot(data_biplo) +
      aes_string(y = "`Nₑ/NT`", size = "`Vk.`", x = "T", colour=temp_dist) +
      geom_point() +
      geom_text_repel(aes(label=Species), size=4)+
      scale_color_continuous(palette = "heat", trans  = 'reverse')+
      # scale_color_gradient() +
      labs(
        y = paste0(diff_focus, " Nₑ/NT"),
        size = "Vk.",
        color = paste0("ratio of ", temp_dist, " distance to\nmean env T°"),
        x = paste0(diff_focus, " Generation length")
      ) +
      theme_light()
    p2 <- ggplot(data_biplo) +
      aes_string(y = "`Nₑ/NT`", x = "Vk.", size = "T", colour=temp_dist) +
      geom_point() +
      geom_text_repel(aes(label=Species), size=4)+
      scale_color_continuous(palette = "heat", trans  = 'reverse')+
      # scale_color_gradient() +
      labs(
        y = paste0(diff_focus, " Nₑ/NT"),
        size = paste0(diff_focus, "T"),
        color = paste0("ratio of ", temp_dist, " distance to\nmean env T°"),
        x = paste0(diff_focus, " Vk.")
      ) +
      theme_light()
    p3 <- ggplot(data_biplo) +
      aes_string(y = "`Nₑ/NT`", x = "Nα", colour=temp_dist) +
      geom_point() +
      geom_text_repel(aes(label=Species), size=4)+
      scale_color_continuous(palette = "heat", trans  = 'reverse')+
      # scale_color_gradient() +
      labs(
        y = paste0(diff_focus, " Nₑ/NT"),
        color = paste0("ratio of ", temp_dist, " distance to\nmean env T°"),
        x = paste0(diff_focus, " Nα")
      ) +
      theme_light()
    p4 <- ggplot(data_biplo) +
      aes_string(y = "`Nₑ/NT`", x = "Nₑ", colour=temp_dist) +
      geom_point() +
      geom_text_repel(aes(label=Species), size=4)+
      scale_color_continuous(palette = "heat", trans  = 'reverse')+
      # scale_color_gradient() +
      labs(
        y = paste0(diff_focus, " Nₑ"),
        color = paste0("ratio of ", temp_dist, " distance to\nmean env T°"),
        x = paste0(diff_focus, " Nα")
      ) +
      theme_light()
    
    ggsave(p1, file=paste0(PATHOUT, "/BIPLOT/", temp_dist, "/Temp_T", diff_focus, temp_dist, ".png"), 
           height = 4.29, width=10.3, create.dir = T)
    ggsave(p2, file=paste0(PATHOUT, "/BIPLOT/", temp_dist, "/Temp_Vk", diff_focus, temp_dist, ".png"), 
           height = 4.29, width=10.3)
    ggsave(p3, file=paste0(PATHOUT, "/BIPLOT/", temp_dist, "/Temp_Nalpha", diff_focus, temp_dist, ".png"),
           height = 4.29, width=10.3)  
    ggsave(p4, file=paste0(PATHOUT, "/BIPLOT/", temp_dist, "/Temp_Ne", diff_focus, temp_dist, ".png"),
           height = 4.29, width=10.3)
    
    list_sp_x <- list()
    for (x in colnames(data_biplo)[-c(1, seq(2, c(ncol(data_biplo)), by=2))]){
      if(x==temp_dist) next
      data_biploi <- data_biplo[,c(which(colnames(data_biplo) %in% c(x, "Species")), ncol(data_biplo))]
      list_sp_x[[paste0(x, "1")]] <- ggplot(data_biploi, aes_string(y= paste0("`", x,"`"), x=temp_dist)) +
        geom_smooth(data=data_reg, aes_string(y= paste0("`",diff_focus, x,"`"), x=temp_dist), color="red", fill="#69b3a2", se=TRUE, method = "lm")+
        geom_point(size=2) +
        geom_text_repel(aes(label=Species), size=3.5) +
        theme_light()+
        ylab(paste0(x, " ", stringr::str_replace(stringr::str_replace(
          stringr::str_replace(stringr::str_replace(diff_focus, "_", " "), "_", " "), 
          pattern = "diff", replacement = "(With CC - Without CC)\n"),
          "abs", "absolute and")))+
        xlab(paste0("ratio of the distance to\n", temp_dist, "with CC and without CC"))
    }
    xarr <- ggarrange(plotlist= list_sp_x[c(1:3)], ncol=3)
    ggsave(xarr, file=paste0(PATHOUT, "/BIPLOT/", temp_dist, "/", stringr::str_replace(x, "/", "_"), diff_focus, temp_dist, ".png"), height = 5, width=12)  
    
    xarr <- ggarrange(plotlist= list_sp_x[-c(1:3)], ncol=4)
    ggsave(xarr, file=paste0(PATHOUT, "/BIPLOT/", temp_dist, "/", stringr::str_replace(x, "/", "_"), diff_focus, temp_dist, "2.png"), height = 5, width=16)
    
  }
}

# plot the biplots and the regression of the absolute value
for (temp_dist in c("Topt_dist", "Tmin_dist", "Tmax_dist")){
  for (diff_focus in dim_focus){ # loop to compare the regressions between the different ways to average calculate the difference
    
    data_reg <- dfbx2[which(dfbx2$env == "mean_temp"), ]
    if(temp_dist == "DOcrit")    data_reg <- dfbx2[which(dfbx2$env == "mean_O2"), ]
    data_reg[which(data_reg$variable %in% c(temp_dist)), diff_focus] <- data_reg$ratio[which(data_reg$variable %in% c(temp_dist))]
    data_reg <- pivot_wider(
      data_reg[which(data_reg$variable %in% c("Nₑ/NT", "NT", "Nₑ",  "T", "Vk.", "T/Vk.", "Nα", temp_dist)), which(colnames(data_reg) %in% 
                                                                                                                    c("Species", "replicate", "variable", "value", diff_focus))], 
      values_from = diff_focus, 
      names_from = "variable")
    colnames(data_reg)[-which(colnames(data_reg) %in% c(temp_dist, "Species"))] <- 
      paste0(diff_focus, colnames(data_reg)[-which(colnames(data_reg)  %in% c(temp_dist, "Species"))])
    
    
    #* side test : lm(Var ~abs(O2)*abs(Temp))
    #* + plot it 
    #* => not convinced graphically & DOcrit varies very little
    # data_reg_t <- dfbx2[which(dfbx2$env == "mean_O2"), ]
    # data_reg_t[which(data_reg_t$variable %in% c("DOcrit")), diff_focus] <- data_reg_t$ratio[which(data_reg_t$variable %in% c("DOcrit"))]
    # data_reg_t <- pivot_wider(
    #   data_reg_t[which(data_reg_t$variable %in% c("DOcrit")), which(colnames(data_reg_t) %in% 
    #                                                                   c("Species", "replicate", "variable", "value", "ratio"))], 
    #   values_from = "ratio", 
    #   names_from = "variable")
    # data_reg <- full_join(data_reg, data_reg_t)
    # data_reg[,paste0(temp_dist, "_abs")] <- as.matrix(abs(data_reg[,temp_dist]))
    # data_reg$DOcrit_abs <- as.matrix(abs(data_reg$DOcrit))
    # for (dim in c("Nₑ/NT", "NT", "Nₑ",  "T", "Vk.", "T/Vk.", "Nα")){
    #   a <- lm(as.matrix(data_reg[, paste0(diff_focus, dim)])~as.matrix(data_reg[,temp_dist])*data_reg$DOcrit_abs)
    #   iteraction_list[[dim]][[temp_dist]][[diff_focus]] <- summary(a)
    #   iteraction_list_plot[[dim]][[temp_dist]][[diff_focus]] <- ggplot(data_reg, aes(label=Species)) +
    #     aes_string(y = "DOcrit_abs", x = paste0(temp_dist, "_abs"), colour=paste0("`", diff_focus, dim, "`")) +
    #     geom_point(size=2) +
    #     geom_text(size=4)+
    #     scale_color_continuous(palette = "heat", trans  = 'reverse')+
    #     # scale_color_gradient() +
    #     labs(
    #       y = "DOcrit_abs",
    #       color = dim,
    #       x = paste0(temp_dist, "_abs")
    #     ) +
    #     theme_light()
    # }
    # next()
    
    
    # data biplot
    data_biplo <- dfbx2[which(dfbx2$env == "mean_temp"), ]
    if(temp_dist == "DOcrit")    data_biplo <- dfbx2[which(dfbx2$env == "mean_O2"), ]
    data_biplo[which(data_biplo$variable %in% c(temp_dist)), diff_focus] <- data_biplo$ratio[which(data_biplo$variable %in% c(temp_dist))] # replace the diff of the temperature by the ratio : that's the indicator we want to focus on 
    data_biplo <- pivot_wider(
      data_biplo[which(data_biplo$variable %in% c("Nₑ/NT", "NT", "Nₑ",  "T", "Vk.", "T/Vk.", "Nα", temp_dist)), which(colnames(data_biplo) %in% 
                                                                                                                        c("Species", "replicate", "variable", "value", diff_focus))], 
      values_from = diff_focus, 
      names_from = "variable")
    
    
    # transform dataset for 
    data_biplo <- data_biplo %>% 
      group_by(Species) %>% 
      summarize(across(c(`Nₑ/NT`,`NT`, `Nₑ`, `T`, `Vk.`, `T/Vk.`, `Nα`, temp_dist), 
                       list(sd = sd, mymean = mean))) %>% 
      ungroup
    colnames(data_biplo) <- stringr::str_remove(colnames(data_biplo), pattern = "_mymean")
    data_biplo <- na.omit(data_biplo)
    
    
    # PLOT 
    p1 <- ggplot(data_biplo) +
      aes_string(y = "`Nₑ/NT`", size = "`Vk.`", x = "T", colour=paste0("abs(",temp_dist, ")")) +
      geom_point() +
      geom_text_repel(aes(label=Species), size=4)+
      scale_color_continuous(palette = "heat", trans  = 'reverse')+
      # scale_color_gradient() +
      labs(
        y = paste0(diff_focus, " Nₑ/NT"),
        size = "Vk.",
        color = paste0("ratio of ", temp_dist, " distance to\nmean env T°"),
        x = paste0(diff_focus, " Generation length")
      ) +
      theme_light()
    p2 <- ggplot(data_biplo) +
      aes_string(y = "`Nₑ/NT`", x = "Vk.", size = "T", colour=paste0("abs(",temp_dist, ")")) +
      geom_point() +
      geom_text_repel(aes(label=Species), size=4)+
      scale_color_continuous(palette = "heat", trans  = 'reverse')+
      # scale_color_gradient() +
      labs(
        y = paste0(diff_focus, " Nₑ/NT"),
        size = paste0(diff_focus, "T"),
        color =paste0("ratio of ", temp_dist, " distance to\nmean env T°"),
        x = paste0(diff_focus, " Vk.")
      ) +
      theme_light()
    p3 <- ggplot(data_biplo) +
      aes_string(y = "`Nₑ/NT`", x = "Nα", colour=paste0("abs(",temp_dist, ")")) +
      geom_point() +
      geom_text_repel(aes(label=Species), size=4)+
      scale_color_continuous(palette = "heat", trans  = 'reverse')+
      # scale_color_gradient() +
      labs(
        y = paste0(diff_focus, " Nₑ/NT"),
        color = paste0("ratio of ", temp_dist, " distance to\nmean env T°"),
        x = paste0(diff_focus, " Nα")
      ) +
      theme_light()
    p4 <- ggplot(data_biplo) +
      aes_string(y = "`Nₑ/NT`", x = "Nₑ", colour=paste0("abs(",temp_dist, ")")) +
      geom_point() +
      geom_text_repel(aes(label=Species), size=4)+
      scale_color_continuous(palette = "heat", trans  = 'reverse')+
      # scale_color_gradient() +
      labs(
        y = paste0(diff_focus, " Nₑ"),
        color = paste0("ratio of ", temp_dist, " distance to\nmean env T°"),
        x = paste0(diff_focus, " Nα")
      ) +
      theme_light()
    
    ggsave(p1, file=paste0(PATHOUT, "/BIPLOT_ABS/", temp_dist, "/Temp_T", diff_focus, temp_dist, ".png"), 
           height = 4.29, width=10.3, create.dir = T)
    ggsave(p2, file=paste0(PATHOUT, "/BIPLOT_ABS/", temp_dist, "/Temp_Vk", diff_focus, temp_dist, ".png"), 
           height = 4.29, width=10.3)
    ggsave(p3, file=paste0(PATHOUT, "/BIPLOT_ABS/", temp_dist, "/Temp_Nalpha", diff_focus, temp_dist, ".png"),
           height = 4.29, width=10.3)  
    ggsave(p4, file=paste0(PATHOUT, "/BIPLOT_ABS/", temp_dist, "/Temp_Ne", diff_focus, temp_dist, ".png"),
           height = 4.29, width=10.3)
    
    list_sp_x <- list()
    for (x in colnames(data_biplo)[-c(1, seq(2, c(ncol(data_biplo)), by=2))]){
      if(x==temp_dist) next
      data_biploi <- data_biplo[,c(which(colnames(data_biplo) %in% c(x, "Species")), ncol(data_biplo))]
      list_sp_x[[paste0(x, "1")]] <- ggplot(data_biploi, aes_string(y= paste0("`", x,"`"), x=paste0("abs(",temp_dist, ")"))) +
        geom_smooth(data=data_reg, aes_string(y= paste0("`",diff_focus, x,"`"), x=paste0("abs(",temp_dist, ")")), color="red", fill="#69b3a2", se=TRUE, method = "lm")+
        geom_point(size=2) +
        geom_text_repel(aes(label=Species), size=3.5) +
        theme_light()+
        ylab(paste0(x, " ", stringr::str_replace(stringr::str_replace(
          stringr::str_replace(stringr::str_replace(diff_focus, "_", " "), "_", " "), 
          pattern = "diff", replacement = "(With CC - Without CC)\n"),
          "abs", "absolute and")))+
        xlab(paste0("ratio of the distance to\n", temp_dist, " with CC and without CC"))
    }
    xarr <- ggarrange(plotlist= list_sp_x[c(1:3)], ncol=3)
    ggsave(xarr, file=paste0(PATHOUT, "/BIPLOT_ABS/", temp_dist, "/", stringr::str_replace(x, "/", "_"), diff_focus, ".png"), height = 5, width=12, create.dir = T)  
    
    xarr <- ggarrange(plotlist= list_sp_x[-c(1:3)], ncol=4)
    ggsave(xarr, file=paste0(PATHOUT, "/BIPLOT_ABS/", temp_dist, "/", stringr::str_replace(x, "/", "_"), diff_focus, "2.png"), height = 5, width=16)
    
  }
}



# Complementary analysis : suggested by BRUNO --------------

### Description of the idea ---------------
#* Hello Marine, le fait que Atlantic herring et Atlantic mackerel drive les relations ne me semble pas 
#* un réel problème. Ce sont effectivement des "outliers" statistiquement mais biologiquement  
#* représentatifs du gradient de position par rapport à l'optimum thermique. Effectivement, on s'attend 
#* donc à ce que la relation soit drivée par les extrêmes du gradient Temp-dist et que pour des valeurs 
#* plus modérées de Temp_dist on ait du bruit autour de la relation. SI on avait pu échantillonner les 
#* espèces, on aurait justement choisi d'essayer de reprrésenter tout le gradient en prenant des espèces 
#* distribuées de manière homogène le long du gradient. 
#* Je pense donc qu'on peut partir là-dessus moyennant 
#* - ce type d'explication dans le texte 
#* - une analyse de sensibilité : le modèle avec les 15 espèces comme tu le fais, le modèle sans herring, un modèle sans mackerel, un modèle sans les 2 ou bien une analyse "leave-one-out" avec le package boot
#* - calculer ou plotter des statistiques d'influence : Cook's distance et Leverage value par exemple avec la fonction influence
#* - Plot du package cartenter une régression robuste avec la fonction rlm du package MASS et comparer les pentes avec ton lm 
#* => rlm based on weighting the residuals, maybe not the best https://stats.oarc.ucla.edu/r/dae/robust-regression/
#* "When fitting a least squares regression, we might find some outliers or high leverage data points. We have decided that these data points are not data entry errors, 
#* neither they are from a different population than most of our data. So we have no compelling reason to exclude them from the analysis. 
#* Robust regression might be a good strategy since it is a compromise between excluding these points entirely from the analysis and 
#* including all the data points and treating all them equally in OLS regression. The idea of robust regression is to weigh the 
#* observations differently based on how well behaved these observations are. Roughly speaking, it is a form of weighted and reweighted least squares regression"
#* 
#* TRY THESE ANALYSIS TO WEIGHT THE RESULTS I JUST HAVE PRODUCED TO RELATE dist_temp to the variables in the model
#* 
#* BEDE's suggestions : 
#* with easystats, performance package : gives me more model diagnostics, to have a better idea of the dataset quality

## Sensitivity analysis -----------

# regressions

for (temp_dist in c("Topt_dist", "Tmin_dist", "Tmax_dist")){

  table_regressionssp <- list()
  table_intercept_full_dataset_sp <- list()
  table_slope_full_dataset_sp <- list()
  for (sp in species){
    
    table_intercept_tempo <- data.frame(matrix(NA, nrow=4, ncol=7))
    table_slope_tempo     <- data.frame(matrix(NA, nrow=4, ncol=7))
    colnames(table_slope_tempo) <- c("Nₑ/NT", "NT", "Nₑ", "T", "Vk.", "T/Vk.", "Nα")
    colnames(table_intercept_tempo) <- c("Nₑ/NT", "NT", "Nₑ", "T", "Vk.", "T/Vk.", "Nα")
    dim_focus <- c("diff_abs", "diff_rel", "diff_abs_scaled", "diff_rel_scaled")
    rownames(table_slope_tempo) <- dim_focus
    rownames(table_intercept_tempo) <- dim_focus
    
    data_regsp <- dfbx2
    data_regsp <- data_regsp[data_regsp$Species != sp, ]
    
    for (diff_focus in dim_focus){ # loop to compare the regressions between the different ways to average calculate the difference
    
      data_reg <- dfbx2[which(dfbx2$env == "mean_temp"), ]
      if(temp_dist == "DOcrit")    data_reg <- dfbx2[which(dfbx2$env == "mean_O2"), ]
      data_reg[which(data_reg$variable %in% c(temp_dist)), diff_focus] <- data_reg$ratio[which(data_reg$variable %in% c(temp_dist))]
      data_reg <- pivot_wider(
      data_reg[which(data_reg$variable %in% c("Nₑ/NT", "NT", "Nₑ",  "T", "Vk.", "T/Vk.", "Nα", temp_dist)), which(colnames(data_reg) %in% 
                                                                                                                      c("Species", "replicate", "variable", "value", diff_focus))], 
      values_from = diff_focus, 
      names_from = "variable")
      colnames(data_reg)[-which(colnames(data_reg) %in% c(temp_dist, "Species"))] <- 
      paste0(diff_focus, colnames(data_reg)[-which(colnames(data_reg)  %in% c(temp_dist, "Species"))])
      
      # Created nested tibble
      nested_df <- tibble(formula = c(paste0("`", diff_focus, "Nₑ/NT` ~", temp_dist),
                                      paste0("`", diff_focus, "NT` ~", temp_dist),
                                      paste0("`", diff_focus, "Nₑ` ~", temp_dist), 
                                      paste0("`", diff_focus, "T` ~ ", temp_dist), 
                                      paste0("`", diff_focus, "Vk.` ~", temp_dist),
                                      paste0("`", diff_focus, "T/Vk.` ~", temp_dist),
                                      paste0(diff_focus, "Nα ~", temp_dist))) %>%
        group_by(ID = formula) %>%
        group_modify(~ as_tibble(data_reg)) %>%
        nest() 
      
      # Get model estimates
      estimates <- nested_df %>%
        mutate(estimates = data %>% map2(ID, ~ tidy(lm(.y, data = .x)))) %>%
        select(-data) %>%
        unnest()
      
      id_intercept <- which(estimates$term=="(Intercept)")
      id_slope <- which(estimates$term==temp_dist)
      ID_intercept <- as.matrix(estimates[id_intercept, "ID"])
      ID_slope <- as.matrix(estimates[id_slope, "ID"])
      order_estimates <- c(grep(pattern=paste0(diff_focus, colnames(table_slope_full_dataset)[1], "` ~"), ID_intercept),
                           grep(pattern=paste0(diff_focus, colnames(table_slope_full_dataset)[2], "` ~"), ID_intercept),
                           grep(pattern=paste0(diff_focus, colnames(table_slope_full_dataset)[3], "` ~"), ID_intercept),
                           grep(pattern=paste0(diff_focus, colnames(table_slope_full_dataset)[4], "` ~"), ID_intercept),
                           grep(pattern=paste0(diff_focus, colnames(table_slope_full_dataset)[5], "` ~"), ID_intercept),
                           grep(pattern=paste0(diff_focus, colnames(table_slope_full_dataset)[6], "` ~"), ID_intercept),
                           grep(pattern=paste0(diff_focus, colnames(table_slope_full_dataset)[7], " ~"), ID_intercept))
      table_intercept_tempo[diff_focus, ] <- as.matrix(estimates[id_intercept, "estimate"])[order_estimates]
      table_slope_tempo[diff_focus, ] <- as.matrix(estimates[id_slope, "estimate"])[order_estimates]
      
      
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
      table_regressionssp[[sp]][[diff_focus]] <- gsubfn("([0-9.]+)", ~replace_numbers(x), star)
      
    }
    
    table_intercept_full_dataset_sp[[sp]] <- table_intercept_tempo
    table_slope_full_dataset_sp[[sp]] <- table_slope_tempo
  }
  
  
  # plot the difference in slope and intercept : 1 plot per variable
  # like a lollipop plot, one row by species 
  # I cannot do stat to compare the significance of the difference because the datasets are not the same between the different iterations (leave one out)
  
  
  list() -> diff_slope_intercept
  for (foc in dim_focus){
    list() ->foc_dataset
    for (dim in c("Nₑ/NT", "NT", "Nₑ", "T", "Vk.", "T/Vk.", "Nα")){
      tempo_dim <- as.data.frame(matrix(ncol=2, nrow=length(species)))
      colnames(tempo_dim) <- c("(Intercept)", "slope")
      rownames(tempo_dim) <- species
      for (sp in species){
        tempo_dim[sp, "(Intercept)"] = table_intercept_full_dataset[foc, dim] - table_intercept_full_dataset_sp[[sp]][foc, dim]
        tempo_dim[sp, "slope"]       = table_slope_full_dataset[foc, dim] - table_slope_full_dataset_sp[[sp]][foc, dim]
        tempo_dim[sp, "sign_slope"]  = ifelse(sign(table_slope_full_dataset[foc, dim]) == sign(table_slope_full_dataset_sp[[sp]][foc, dim]), "slope (same sign)", "slope (sign changes)")
        tempo_dim[sp, "sign_Intercept"]= ifelse(sign(table_intercept_full_dataset[foc, dim]) == sign(table_intercept_full_dataset_sp[[sp]][foc, dim]), "intercept (same sign)", "intercept (sign changes)")
      }
      foc_dataset[[dim]] = tempo_dim
    }
    diff_slope_intercept[[foc]] = foc_dataset
  }

  ##plot lollipop difference in intercept, slope of the Var~T° regression with and without CC
  for (foc in dim_focus){
    lolliplot =NULL
    diff_i <- diff_slope_intercept[[foc]]
    n_diff_i <- names(diff_i)
    diff_ii <-  do.call(rbind, lapply(1:length(diff_i), FUN = function(x){
      y<- as.data.frame(diff_i[[x]])
      y$Species <- rownames(y)
      y$variable <- n_diff_i[x]
      rownames(y) <- NULL
      return(y)
    }))
    if (!grepl(pattern="scale", foc)){
      diff_ii[which(diff_ii$variable %in% c("NT", "Nₑ", "Nα")), which(colnames(diff_ii) %in% c("slope", "(Intercept)"))] <-
        c(10^(-10)) * (diff_ii[which(diff_ii$variable %in% c("NT", "Nₑ", "Nα")), which(colnames(diff_ii) %in% c("slope", "(Intercept)"))]) # intercept is rescaled to have all variables in the same range more or less
      diff_ii[which(diff_ii$variable %in% c("NT", "Nₑ", "Nα")),"variable"] <-
        paste0(as.matrix(diff_ii[which(diff_ii$variable %in% c("NT", "Nₑ", "Nα")),"variable"]), " (e-10)")
    }
    lolliplot[["Ne"]] <- ggplot(diff_ii[which(diff_ii$variable %in% c("NT", "Nₑ", "Nα", "Nₑ/NT")), ]) +
      geom_hline(yintercept = 0, color = "darkgrey")+
      geom_point(aes(x=Species, y=`(Intercept)`, fill=sign_Intercept, pch=sign_Intercept), col="black", size=3, alpha=0.5) +
      geom_point(aes(x=Species, y=slope, fill=sign_slope, pch=sign_slope), col="black", size=3, alpha=1) +
      coord_flip()+
      theme_bw() +
      xlab("") +
      ylab("Leave-one-out slope and intercept difference")+
      facet_wrap(~variable, ncol=4)+
      scale_fill_manual(values = c("lightgreen", "darkgreen", "lightblue", "darkblue"),
                         guide  = guide_legend(), name="") +
      scale_shape_manual(values = c(21, 21, 22, 22),
                         guide  = guide_legend(), name="") +
      theme(legend.position = "right",
            panel.border    = element_blank())
    lolliplot[["Vk"]] <- ggplot(diff_ii[-which(diff_ii$variable %in% c("NT", "Nₑ", "Nα", "Nₑ/NT")), ]) +
      geom_hline(yintercept = 0, color = "darkgrey")+
      geom_point(aes(x=Species, y=`(Intercept)`, fill=sign_Intercept, pch=sign_Intercept), col="black", size=3, alpha=0.5) +
      geom_point(aes(x=Species, y=slope, fill=sign_slope, pch=sign_slope), col="black", size=3, alpha=1) +
      coord_flip()+
      theme_bw() +
      xlab("") +
      ylab("Leave-one-out slope and intercept difference")+
      facet_wrap(~variable)+
      scale_fill_manual(values = c("lightgreen", "darkgreen", "lightblue", "darkblue"),
                         guide  = guide_legend(), name="") +
      scale_shape_manual(values = c(21,21,22,22),
                         guide  = guide_legend(), name="") +
      theme(legend.position = "right",
            panel.border    = element_blank())
    ggsave(lolliplot, file=paste0(PATHOUT, "LOLLI/", temp_dist, "/Var_Temp_", foc, ".pdf"), create.dir = T)
  }
  
  
  # biplots of the slope and intercept change with the leave-on-out method
  for (diff_focus in dim_focus){ # loop to compare the biplots between the different ways to average calculate the difference
    
    # select the column of interest 
    data_biplo <- dfbx2[which(dfbx2$env == "mean_temp"), ]
    if(temp_dist == "DOcrit")    data_biplo <- dfbx2[which(dfbx2$env == "mean_O2"), ]
    data_biplo[which(data_biplo$variable %in% c(temp_dist)), diff_focus] <- data_biplo$ratio[which(data_biplo$variable %in% c(temp_dist))] # replace the diff of the temperature by the ratio : that's the indicator we want to focus on 
    data_biplo <- pivot_wider(
      data_biplo[which(data_biplo$variable %in% c("Nₑ/NT", "NT", "Nₑ",  "T", "Vk.", "T/Vk.", "Nα", temp_dist)), which(colnames(data_biplo) %in% 
                                                                                                                          c("Species", "replicate", "variable", "value", diff_focus))], 
      values_from = diff_focus, 
      names_from = "variable")
    
    
    # transform dataset for 
    data_biplo <- data_biplo %>% 
      group_by(Species) %>% 
      summarize(across(c(`Nₑ/NT`,`NT`, `Nₑ`, `T`, `Vk.`, `T/Vk.`, `Nα`, temp_dist), 
                       list(sd = sd, mymean = mean))) %>% 
      ungroup
    colnames(data_biplo) <- stringr::str_remove(colnames(data_biplo), pattern = "_mymean")
    data_biplo <- na.omit(data_biplo)
    
    for (sp in species){
      
      data_biplosp <- data_biplo[which(data_biplo$Species != sp), ]
    
      list_sp_x <- list()
      for (x in colnames(data_biplosp)[-c(1, seq(2, c(ncol(data_biplosp)), by=2))]){
        if(x==temp_dist) next
        data_biploi <- data_biplosp[,c(which(colnames(data_biplosp) %in% c(x, "Species")), ncol(data_biplosp))]
        list_sp_x[[x]] <- ggplot(data_biploi, aes_string(y= paste0("`", x,"`"), x=temp_dist)) +
          geom_smooth(col="grey", alpha=0.5, method = "lm")+
          geom_point(size=2) +
          geom_text_repel(aes(label=Species), size=3.5) +
          theme_light()+
          ylab(paste0(diff_focus, " ", x)) +
          xlab(paste0("Temperature distance to ", temp_dist, " ratio"))
      }
      xarr <- ggarrange(plotlist=list_sp_x)
      ggsave(plot = xarr, filename=paste0(PATHOUT, "/BIPLOT_cross_valid/", temp_dist, "/", sp, "_", stringr::str_replace(x, "/", "_"), diff_focus, ".png"), create.dir = T,
             height = 9, width=12)
    }
  }
}


## Influence statistics -----------

## Robust regression : weighting the regression elements to their weight : less weighted when they influence a lot the regression ----------
# I doubt this is what I am loooking for : the analysis is based on the residuals, and the residuals = predicted - point value
# atalntic H and atalntic M both have residuals =0, while it is what we want to use to weight them lower 
# 
# lm_test <- lm(data=data_reg, formula=`diff_rel_scaledNₑ/NT`~Temp_dist)
# opar <- par(mfrow = c(2,2), oma = c(0, 0, 1.1, 0))
# plot(lm_test, las = 1)
# par(opar)
# data_reg[c(208, 117, 110), ]

















# # more complex regressions
# table_regressions2 <- list()
# for (diff_focus in c("diff_abs", "diff_rel", "diff_abs_scale", "diff_rel_scale")){ # loop to compare the regressions between the different ways to average calculate the difference
#   
#   data_reg <- dfbx2
#   data_reg[which(data_reg$variable %in% c(temp_dist)), diff_focus] <- data_reg$ratio[which(data_reg$variable %in% c(temp_dist))]
#   data_reg <- pivot_wider(
#     data_reg[which(data_reg$variable %in% c("Nₑ/NT", "NT", "Nₑ",  "T", "Vk.", "T/Vk.", "Nα", temp_dist)), which(colnames(data_reg) %in% 
#                                                                                                                     c("Species", "replicate", "variable", "value", diff_focus))], 
#     values_from = diff_focus, 
#     names_from = "variable")
#   colnames(data_reg)[-which(colnames(data_reg) %in% c(temp_dist, "Species"))] <- 
#     paste0(diff_focus, colnames(data_reg)[-which(colnames(data_reg)  %in% c(temp_dist, "Species"))])
#   
#   
#   # Created nested tibble
#   nested_df <- tibble(formula = c(paste0("`", diff_focus, "Nₑ/NT` ~ `", diff_focus, "Nₑ`*Temp_dist"),
#                                   paste0("`", diff_focus, "Nₑ/NT` ~ `", diff_focus, "T`*Temp_dist"),
#                                   paste0("`", diff_focus, "Nₑ/NT` ~ `", diff_focus, "NT`*Temp_dist"), 
#                                   paste0("`", diff_focus, "Nₑ/NT` ~ `", diff_focus, "T/Vk.`*Temp_dist"), 
#                                   paste0("`", diff_focus, "Nₑ/NT` ~ `", diff_focus, "Vk.`*Temp_dist"),
#                                   paste0("`", diff_focus, "Nₑ/NT` ~ `", diff_focus, "Nα`*Temp_dist"))) %>%
#     group_by(ID = formula) %>%
#     group_modify(~ as_tibble(data_reg)) %>%
#     nest() 
#   
#   # Get model estimates
#   nested_df %>%
#     mutate(estimates = data %>% map2(ID, ~ tidy(lm(.y, data = .x)))) %>%
#     select(-data) %>%
#     unnest()
#   # Get fitted values and residuals
#   nested_df %>%
#     mutate(model = ID %>% map2(data, lm),
#            stats = model %>% map(augment)) %>%
#     select(-data, -model) %>%
#     unnest() 
#   # View(nested_df)
#   star <- nested_df %>%
#     mutate(model = ID %>% map2(data, ~ lm(.x, .y))) %>%
#     pull(model) %>%
#     stargazer(type = "text", digits = 3, header = F, digit.separator = '', align = T)
#   replace_numbers = function(x, cutoff=4, digits=3, scipen=-7) {
#     ifelse(nchar(x) < cutoff, x, prettyNum(as.numeric(x), digits=digits, scientific=scipen))
#   }
#   table_regressions2[[diff_focus]] <- gsubfn("([0-9.]+)", ~replace_numbers(x), star)
# }

