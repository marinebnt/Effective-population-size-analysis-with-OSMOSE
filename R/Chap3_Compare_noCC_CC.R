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
# Library
library(ggplot2)
library(dplyr)
library(hrbrthemes)


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
    
    ## Extract only the dimensions I will be using 
    NeNT <- final_dataset_age[which(final_dataset_age$variable %in% c("Deltasquare", "eggNumber_sc", "Mpred", "Mstar",
                                                                      "Dead", "Dead_F", "Dead_Madd", "Dead_Mfor", "Dead_Mnat", "Dead_Mpred", 
                                                                      "Dead_Mstar", "Dead_Mtot", "survivorship_curve", "SSDgx")), 
                              -which(colnames(final_dataset_age) %in% c("Fstart"))] %>% ungroup()     
    NeNT2 <- final_dataset_age[which(final_dataset_age$variable %in% c("Ne_N_alpha", "NT", "generation_length", "Vk_point", "Numerator_N_alpha")), 
                               -which(colnames(final_dataset_age) %in% c("Fstart", "age"))] %>% ungroup()    
    
    ## Increase the dimension of the population size at maturation age (Numerator_N_alpha) (add age dimension)
    ## to be able to manipulate it with other variables => estimate SSDgx/Nalpha
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
    
    ## Cleaning environment
    rm(final_dataset_age) 
    NeNT[,"Scenario"] <- ifelse(Scenario == "Chap2", "No CC", "With CC") # *********
    NeNT2[,"Scenario"] <- ifelse(Scenario == "Chap2", "No CC", "With CC") # *********
    
    ## Binding the Species together
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
  
  ## Binding the Scenarios together
  if (is.null(Ne_NT_tot)) {Ne_NT_tot <- Ne_NT 
  } else Ne_NT_tot <- rbind(Ne_NT_tot, Ne_NT)   
  if (is.null(Ne_NT_tot2)) {Ne_NT_tot2 <- Ne_NT2 
  } else Ne_NT_tot2 <- rbind(Ne_NT_tot2, Ne_NT2) 
}#scenario

## SLIGHT ADDITION ON NE/NT DATA -----------
Ne_NT_tot2 <- Ne_NT_tot2[-which(Ne_NT_tot2$NT < 10), ]
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

## KEEPING ONLY THE FMAX=1 DIMENSION ----------------------
Ne_NT_tot2 <- Ne_NT_tot2[which(Ne_NT_tot2$Fmax == 1), ]
dt_tf <- dt_tf[which(dt_tf$Fmax == 1), ]

## BOXPLOTS ------------------
Ne_NT_bx <- pivot_longer(Ne_NT_tot2, cols=-c("Fmax","replicate", "Species", "Scenario"),
                         values_to = "value", names_to = "variable")
dfbx <- rbind(Ne_NT_bx, dt_tf[, -which(colnames(dt_tf) == "Fstart")])

for (x in names(table(dfbx$variable))[14]){
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



# dataset for barplot & lollipop plot
dfbx <- pivot_wider(dfbx, names_from = "Scenario", values_from = "value")

# Reorder data using average? Learn more about reordering in chart #267
data <- dfbx %>% 
  group_by(Species, variable) %>% 
  summarize(across(c(`No CC`, `With CC`), 
            list(sd = sd, mymean = mean))) %>% 
  ungroup

data_barplot <- dfbx
data_barplot <- data_barplot[which(data_barplot$variable %in% c("Ne_NT", "NT", "Ne_N_alpha")),]
# plot horizontal barplot
# Create data
# Barplot
data_barplot$Species <- as.factor(data_barplot$Species)
data_barplot$variable <- improve_vriable_names(data_barplot$variable)
# Stacked + percent
# 
# ggplot(data_barplot, aes(x=Species, y=value, group=interaction(Scenario,variable))) + 
#   geom_col(position="fill", stat="identity", aes(fill=as.factor(Scenario))) + 
#   coord_flip() + 
#   scale_fill_manual(name = "", values = c("No CC"="#f79326", "With CC"="gray"), labels = c("No CC", "With CC")) +
#   facet_grid(~variable, scale="free")


# lollipop plot mais avec juste la différence
data_barplot2 <- data_barplot #pivot_wider(data_barplot, names_from = "Scenario", values_from = "value")
data_barplot2$diff = (data_barplot2$`With CC` - data_barplot2$`No CC`)
data_barplot2[which(data_barplot2$variable %in% c("NT", "Nₑ")),grep("diff", colnames(data_barplot2))] <-
  c(10^(-10)) * (data_barplot2[which(data_barplot2$variable %in% c("NT", "Nₑ")),grep("diff", colnames(data_barplot2))])
data_barplot2[which(data_barplot2$variable %in% c("NT", "Nₑ")),"variable"] <-
  paste0(as.matrix(data_barplot2[which(data_barplot2$variable %in% c("NT", "Nₑ")),"variable"]), " (e-10)")

data_barplot2 <- data_barplot2 %>% 
  group_by(Species, variable) %>% 
  summarize(across(c(`diff`), 
                   list(sd = sd, mymean = mean))) %>% 
  ungroup

data_barplot2$variable <- factor(data_barplot2$variable, levels = c("Nₑ/NT", "NT (e-10)", "Nₑ (e-10)"))
data_barplot2$Species <- factor(data_barplot2$Species, 
                                 levels = as.matrix(data_barplot2[which(data_barplot2$variable=="Nₑ/NT"), "Species"])[order(as.matrix(data_barplot2[which(data_barplot2$variable=="Nₑ/NT"), "diff_mymean"]))])
lollipop <- ggplot(data_barplot2, aes(x=Species, y=diff_mymean))  + 
  geom_hline(yintercept = 0, color = "darkgrey")+
  geom_segment(aes(x=Species, xend=Species, y=0, yend=diff_mymean), color="grey") +
  geom_point(color=ifelse(round(data_barplot2$diff_mymean, digits=2)<0, 
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
  facet_wrap(~variable, scale="free_x", ncol=5)+
  scale_y_continuous(n.breaks = 8)

ggsave(lollipop, file=paste0(PATHOUT, "/lollipopNe.png"), 
       height = 4.29, width=10.3)


data_barplot <- dfbx
data_barplot <- data_barplot[which(data_barplot$variable %in% c("Numerator_N_alpha", "generation_length", "Vk_point", "T/Vk")),]
# plot horizontal barplot
# Create data
# Barplot
data_barplot$Species <- as.factor(data_barplot$Species)
data_barplot$variable <- improve_vriable_names(data_barplot$variable)
# Stacked + percent
# lollipop plot mais avec juste la différence
data_barplot3 <- data_barplot #pivot_wider(data_barplot, names_from = "Scenario", values_from = "value")
data_barplot3$diff = (data_barplot3$`With CC` - data_barplot3$`No CC`)
data_barplot3[which(data_barplot3$variable %in% c("Nα")),grep("diff", colnames(data_barplot3))] <-
  c(10^(-10)) * (data_barplot3[which(data_barplot3$variable %in% c("Nα")),grep("diff", colnames(data_barplot3))])
data_barplot3[which(data_barplot3$variable %in% c("Nα")),"variable"] <-
  paste0(as.matrix(data_barplot3[which(data_barplot3$variable %in% c("Nα")),"variable"]), " (e-10)")

data_barplot3 <- data_barplot3 %>% 
  group_by(Species, variable) %>% 
  summarize(across(c(`diff`), 
                   list(sd = sd, mymean = mean))) %>% 
  ungroup

data_barplot3$variable <- factor(data_barplot3$variable, levels = c("T", "Vk.", "T/Vk", "Nα (e-10)"))
data_barplot3$Species <- factor(data_barplot3$Species, 
                                levels = as.matrix(data_barplot2[which(data_barplot2$variable=="Nₑ/NT"), "Species"])[
                                  order(as.matrix(data_barplot2[which(data_barplot2$variable=="Nₑ/NT"), "diff_mymean"]))])
lollipop <- ggplot(data_barplot3, aes(x=Species, y=diff_mymean))  + 
  geom_hline(yintercept = 0, color = "darkgrey")+
  geom_segment(aes(x=Species, xend=Species, y=0, yend=diff_mymean), color="grey") +
  geom_point(color=ifelse(round(data_barplot3$diff_mymean, digits=2)<0, 
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
  facet_wrap(~variable, scale="free_x", ncol=5)+
  scale_y_continuous(n.breaks = 8)

ggsave(lollipop, file=paste0(PATHOUT, "/lollipopVk.png"), 
       height = 4.29, width=10.3)


# data for ggplots -------------

dfbx2 <- rbind(Ne_NT_bx, dt_tf[, -which(colnames(dt_tf) == "Fstart")])
dfbx2 <- dfbx2[which(dfbx2$variable %in% c("Ne_N_alpha", "Ne_NT", "NT", 
                                           "Numerator_N_alpha", "Vk_point",
                                           "generation_length",
                                           "T/Vk")), ]
# test data before regression 
ggplotdata <- dfbx2
ggplotdata$value <- log(ggplotdata$value)
ggplotdata$id <- 1:nrow(ggplotdata) 
ggplot(ggplotdata)+ 
  geom_density(aes(value))+ 
  facet_wrap(~variable, scales="free")
ggplot(ggplotdata, aes(sample=value))+ 
  stat_qq()+ 
  stat_qq_line()+ 
  facet_wrap(~variable, scales="free")
# transform dataset for 
dfbx2 <- pivot_wider(dfbx2, names_from = "variable", values_from = "value")# summarize 
dfbx2 <- dfbx2 %>% 
  group_by(Species, Scenario) %>% 
  summarize(across(c(`Ne_N_alpha`, `Ne_NT`, NT, 
                     `Numerator_N_alpha`, Vk_point, `generation_length`,
                     `T/Vk`), 
                   list(sd = sd, mymean = mean))) %>% 
  ungroup
colnames(dfbx2) <- improve_vriable_names(colnames(dfbx2))
colnames(dfbx2) <- stringr::str_remove(colnames(dfbx2), pattern = "_mymean")

## add the temperature data  ----------------
## 1- Prepare dataset with the parameters from the configuration process (optimal temperature) ----------
# Topt dataset & Mean T from North sea by z layer
df_Topt <- read.csv(paste0(here::here(), "/DATA/FromOSMOSE-nsConfiguration_outputs_global_complet.csv"))[, c("Species", "Final_Topt",                                                                                                            "Final_Tmax", "Final_Tmin")]
matrix_mean_T_sp <- get(load(file=paste0(pathoutputRDATA, "/Tmean_by_species.RData")))
# match up scientific names to species names in OSMOSE
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
## add the mean T° encountered by each species and calculate : 
## the distance of env T° to Tmean, divided by the range of T° => (Tmean-Topt)/(Tmax-Tmin) 
df_Topt <- full_join(df_Topt, matrix_mean_T_sp)
df_Topt <- pivot_longer(df_Topt, cols = c("No CC", "With CC"), values_to = "mean_Temp_env", names_to = "Scenario")
df_Topt$Temp_dist <- (df_Topt$Final_Topt - df_Topt$mean_Temp_env)/(df_Topt$Final_Tmax-df_Topt$Final_Tmin)
dfbi2 <- full_join(df_Topt, dfbx2)



d <- pivot_wider(dfbi2[, which(colnames(dfbi2) %in% c("Species", "Scenario", "Temp_dist"))], 
                 names_from = "Scenario", values_from = "Temp_dist")
d$`With CC`/d$`No CC`

cor(df_Topt[,c("Final_Tmax", "Final_Topt", "Final_Tmin")])


# simple regressions
lm_gene_length_temp <- lm(log(dfbx2$T) ~ dfbx2$Temp_dist*dfbx2$Scenario)
summary(lm_gene_length_temp)
lm_Ne_temp <- lm(log(dfbx2$Nₑ) ~ dfbx2$Temp_dist*dfbx2$Scenario)
summary(lm_Ne_temp)
lm_Vkpoint_temp <- lm(log(dfbx2$Vk.) ~ dfbx2$Temp_dist*dfbx2$Scenario)
summary(lm_Vkpoint_temp)
lm_NT_temp <- lm(log(dfbx2$NT) ~ dfbx2$Temp_dist*dfbx2$Scenario)
summary(lm_NT_temp)
lm_NeNT_temp <- lm(log(dfbx2$`Nₑ/NT`) ~ dfbx2$Temp_dist*dfbx2$Scenario)
summary(lm_NeNT_temp)
lm_Nalpha_temp <- lm(log(dfbx2$Nα) ~ dfbx2$Temp_dist*dfbx2$Scenario)
summary(lm_Nalpha_temp)


library(broom)
# Created nested tibble
nested_df <- tibble(formula = c("log(Nₑ/NT) ~ Temp_dist*Scenario",
                                "log(NT) ~ Temp_dist*Scenario",
                                "log(Nₑ) ~ Temp_dist*Scenario", 
                                "log(T) ~ Temp_dist*Scenario", 
                                "log(Vk.) ~Temp_dist*Scenario",
                                "log(T/Vk.) ~Temp_dist*Scenario",
                                "log(Nα) ~Temp_dist*Scenario")) %>%
  group_by(ID = formula) %>%
  group_modify(~ as_tibble(dfbx2)) %>%
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
library(stargazer)
library(gsubfn)
star <- nested_df %>%
  mutate(model = ID %>% map2(data, ~ lm(.x, .y))) %>%
  pull(model) %>%
  stargazer(type = "text", digits = 3, header = F, digit.separator = '', align = T)
replace_numbers = function(x, cutoff=4, digits=3, scipen=-7) {
  ifelse(nchar(x) < cutoff, x, prettyNum(as.numeric(x), digits=digits, scientific=scipen))
}
table <- gsubfn("([0-9.]+)", ~replace_numbers(x), star)



# PLOT summarize 
p1 <- ggplot(dfbx2) +
  aes(y = log(`Nₑ/NT`), size = `Vk.`, x = log(T), colour=Temp_dist) +
  geom_point() +
  geom_text_repel(aes(label=Species), size=4)+
  scale_color_continuous(palette = "heat", trans  = 'reverse')+
  # scale_color_gradient() +
  labs(
    y = "log Nₑ/NT",
    size = "Vk.",
    color = "Topt distance to\nmean env T°",
    x = "log Generation length"
  ) +
  theme_bw()+
  facet_wrap(~Scenario)
p2 <- ggplot(dfbx2) +
  aes(y = log(`Nₑ/NT`), size =T, x = log(`Vk.`), colour=Temp_dist) +
  geom_point() +
  geom_text_repel(aes(label=Species), size=4)+
  scale_color_continuous(palette = "heat", trans  = 'reverse')+
  # scale_color_gradient() +
  labs(
    y = "log Nₑ/NT",
    size = "T",
    color = "Topt distance to\nmean env T°",
    x = "log Vk."
  ) +
  theme_bw()+
  facet_wrap(~Scenario)
p3 <- ggplot(dfbx2) +
  aes(y = log(`Nₑ/NT`),x = log(Nα), colour=Temp_dist) +
  geom_point() +
  geom_text_repel(aes(label=Species), size=4)+
  scale_color_continuous(palette = "heat", trans  = 'reverse')+
  # scale_color_gradient() +
  labs(
    y = "log Nₑ/NT",
    color = "Topt distance to\nmean env T°",
    x = "log Nα"
  ) +
  theme_bw()+
  facet_wrap(~Scenario)

ggsave(p1, file=paste0(PATHOUT, "/Temp_T.png"), 
       height = 4.29, width=10.3)
ggsave(p2, file=paste0(PATHOUT, "/Temp_Vk.png"), 
       height = 4.29, width=10.3)
ggsave(p3, file=paste0(PATHOUT, "/Temp_Nalpha.png"), 
       height = 4.29, width=10.3)



# plot redundancy analysis # https://pmassicotte.github.io/stats-denmark-2019/07_rda.html#/global-significance
library(vegan)
dfrda <- rbind(Ne_NT_bx, dt_tf[, -which(colnames(dt_tf) == "Fstart")])
dfrda <- dfrda %>%ungroup()
dfrda <- dfrda %>%
  group_by(Species, Scenario, variable) %>%
  mutate(value_sc = scale(value)) %>%
  ungroup()
dfrda <- pivot_wider(dfrda[, -which(colnames(dfrda) %in% c("value"))], names_from = "variable", values_from="value_sc")
dfrda <- full_join(dfrda, df_Topt)
dfrda <- dfrda[-which(dfrda$Species == "Shrimp"),]
dfrda$Scenario <- stringr::str_remove(dfrda$Scenario, " ")
dfrda <- na.omit(dfrda) # shrimp removed
dfrda_Y <- dfrda[, which(colnames(dfrda) %in% c("Ne/NT", "NT", "Ne_N_alpha", "Vk_point", "generation_length"))]
rda_1 <- rda(dfrda_Y~Scenario+Temp_dist+Mpred+Mstar+eggNumber_sc, data=dfrda)

summary(rda_1)
scl <- 3 ## scaling = 3
colvec <- c("red2", "green4", "mediumblue")
ordiplot(rda_1)
points(rda_1, "species", pch=17, col="black")

eco <- names(table(dfrda$Scenario))
bg <- c("#ff7f00","#1f78b4","#ffff33","#a6cee3","#33a02c","#e31a1c", "pink", 
        'purple', "darkblue", "orange4", "grey", "black", "red4", "lightgreen", "green4") # 6 nice colors for our ecotypes

plot(rda_1, type="n", scaling=3, choices=c(1,2))
points(rda_1, display="species", pch=20, cex=0.7, col="gray32", scaling=3, choices=c(1,2))
points(rda_1, display="sites", pch=21, cex=1.3, col="gray32", scaling=3, bg=bg[c(1,2)], choices=c(1,2))
text(rda_1, scaling=3, display="bp", col="#0868ac", cex=1, choices=c(1,2))
legend("topleft", legend=eco, bty="n", col="gray32", pch=21, cex=1, pt.bg=bg)

plot(rda_1, type="n", scaling=3, choices=c(1,2))
points(rda_1, display="species", pch=21, cex=0.7, col="gray32", scaling=3, choices=c(1,2))
points(rda_1, display="sites", pch=21, cex=1.3, col="gray32", scaling=3,  bg=bg,choices=c(1,2))
text(rda_1, scaling=3, display="bp", col="#0868ac", cex=1, choices=c(1,2))
text(rda_1, scaling=3, display = "sites", col="#0868ac", cex=1, choices=c(1,2))
legend("topleft", legend=eco, bty="n", col="gray32", pch=21, cex=1, pt.bg=bg)

library(ggord)
ggord(rda_1, dfrda$Species) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())


RsquareAdj(rda_1) # total variance explained by RDA : extremely low
# permutation test 
anova.cca(rda_1, permutations = 999)
anova.cca(rda_1, by="axis")
anova.cca(rda_1, by="term")
# multicolinearity, if >2 : considered high
sqrt(vif.cca(rda_1))



# Plot
## data_lollipop
data_lollipop <- data[which(data$variable %in% c("Ne_NT", "NT", "Ne_N_alpha")),]
data_lollipop$groups_ranges <- ifelse(data_lollipop$Species %in% c("AtlanticHerring", "Sandeel", "EuropeanSprat"), 
                                      "", " ")
data_lollipop[which(data_lollipop$variable %in% c("NT", "Ne_N_alpha")),grep("CC", colnames(data_lollipop))] <-
  c(10^(-10)) * (data_lollipop[which(data_lollipop$variable %in% c("NT", "Ne_N_alpha")),grep("CC", colnames(data_lollipop))])
data_lollipop[which(data_lollipop$variable %in% c("NT", "Ne_N_alpha")),"variable"] <-
  paste0("(10^-11) ", as.matrix(data_lollipop[which(data_lollipop$variable %in% c("NT", "Ne_N_alpha")),"variable"]))
data_lollipop$variable <- improve_vriable_names(data_lollipop$variable)
data_lollipop$group_for_plot <- paste0(data_lollipop$Species, data_lollipop$groups_ranges)
## plot lollipop 
a <- ggplot(data_lollipop[which(data_lollipop$groups_ranges == ""),]) +
  geom_segment( aes(x=Species, xend=Species, y=`No CC_mymean`, yend=`With CC_mymean`), 
                color="grey") +
  geom_point( aes(x=Species, y=`No CC_mymean`), color=rgb(0.2,0.7,0.1,0.5), size=3 ) +
  geom_point( aes(x=Species, y=`With CC_mymean`), color=rgb(0.7,0.2,0.1,0.5), size=3 ) +
  coord_flip()+
  theme_ipsum() +
  geom_errorbar(aes(x=Species, y=`With CC_mymean`,
                    ymin = `With CC_mymean`-`With CC_sd`, 
                    ymax = `With CC_mymean`+`With CC_sd`),
                color=rgb(0.7,0.2,0.1,0.5))+
  geom_errorbar(aes(x=Species, y=`No CC_mymean`, 
                    ymin = `No CC_mymean`-`No CC_sd`, 
                    ymax = `No CC_mymean`+`No CC_sd`), 
                color=rgb(0.2,0.7,0.1,0.5))+
  xlab("") +
  ylab("")+
  # facet_wrap(~variable, scale="free")+
  facet_grid(~variable, scales = "free")
b <- ggplot(data_lollipop[which(data_lollipop$groups_ranges == " "),]) +
  geom_segment( aes(x=Species, xend=Species, y=`No CC_mymean`, yend=`With CC_mymean`), 
                color="grey") +
  geom_point( aes(x=Species, y=`No CC_mymean`), color=rgb(0.2,0.7,0.1,0.5), size=3 ) +
  geom_point( aes(x=Species, y=`With CC_mymean`), color=rgb(0.7,0.2,0.1,0.5), size=3 ) +
  coord_flip()+
  theme_ipsum() +
  geom_errorbar(aes(x=Species, y=`With CC_mymean`,
                    ymin = `With CC_mymean`-`With CC_sd`, 
                    ymax = `With CC_mymean`+`With CC_sd`),
                color=rgb(0.7,0.2,0.1,0.5))+
  geom_errorbar(aes(x=Species, y=`No CC_mymean`, 
                    ymin = `No CC_mymean`-`No CC_sd`, 
                    ymax = `No CC_mymean`+`No CC_sd`), 
                color=rgb(0.2,0.7,0.1,0.5))+
  xlab("") +
  ylab("")+
  # facet_wrap(~variable, scale="free")+
  facet_grid(~variable, scales = "free")

ggarrange(a,b, ncol=1)

# plot horizontal barplot
# Create data
# Barplot : marche pas en fill, erreurs dans quelles sont les min and max
data_lollipop$Species <- as.factor(data_lollipop$Species)
ggplot(data_lollipop, aes(x=Species, y=`No CC_mymean`)) + 
  geom_bar(alpha=0.5, col="darkred") +
  # geom_bar(aes(x=Species, y=`With CC_mymean`), alpha=0.5, col="lightblue") +
  coord_flip()

ggplot(data_lollipop, aes(x=Species, y=`No CC_mymean`)) + 
  geom_bar(stat = "identity", width=0.2) 

# Stacked + percent
ggplot(data, aes(fill=condition, y=value, x=specie)) + 
  geom_bar(position="fill", stat="identity")






## NE/NT variables READY FOR PLOT : CALCULATE DIFFERENCE & P-VALUE -----------

Ne_NT_tot2 <- pivot_longer(Ne_NT_tot2, values_to = "value", names_to = "variable", cols = -c("Fmax", "replicate", "Species", "Scenario"))

## add significance level of the difference with and without CC, between replicates
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

## Estimate the difference between scenarios variable/species/Fmax/replicate
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

## HEATMAP DATASET ------------
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


## TF VARIABLE READY FOR PLOT : CALCULATE DIFFERENCE & P-VALUE ----------

## Add significance level of the difference of slope with and without CC, between replicates
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

## Estimate the difference between scenarios variable/species/Fmax/replicate
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

## Average of the difference between replicates & Fmax
dt_tf_diff <- dt_tf_diff %>%
  group_by(Species, variable) %>%
  summarize(mean_diff=mean(difference, na.rm=T)) %>%
  ungroup()
## Add Ne/NT response type when NO CC
dt_tf_diff$Response <- ifelse(dt_tf_diff$Species %in% c("Hake", "Whiting", "Saithe", "Sandeel", "AtlanticHerring"), "Increasing N\u2091/NT", 
                              ifelse(dt_tf_diff$Species %in% c("CommonDab", "GreyGurnard", "EuropeanPlaice"), "Decreasing N\u2091/NT", 
                                     ifelse(dt_tf_diff$Species %in% c("NorwayPout", "CommonSole", "AtlanticCod", "EuropeanSprat", "Haddock", "HorseMackerel", 
                                                                      "AtlanticMackerel"), "Constant N\u2091/NT", "out"))) 
## Add p-value
dt_tf_diff <- full_join(dt_tf_diff, stat.test[, which(colnames(stat.test) %in% c("Species", "variable", "p.adj.signif"))])

## HEATMAP DATASET ------------
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

# Topt dataset & Mean T from North sea by z layer
df_Topt <- read.csv(paste0(here::here(), "/DATA/FromOSMOSE-nsConfiguration_outputs_global_complet.csv"))[, c("Species", "Final_Topt", 
                                                                                                             "Final_Tmax", "Final_Tmin")]
matrix_mean_T_sp <- get(load(file=paste0(pathoutputRDATA, "/Tmean_by_species.RData")))

# match up scientific names to species names in OSMOSE
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

## add the mean T° encountered by each species and calculate : 
## the distance of env T° to Tmean, divided by the range of T° => (Tmean-Topt)/(Tmax-Tmin) 
df_Topt <- full_join(df_Topt, matrix_mean_T_sp)
df_Topt$Temp_dist_NoCC <- (df_Topt$`No CC` - df_Topt$Final_Topt)/(df_Topt$Final_Tmax-df_Topt$Final_Tmin)
df_Topt$Temp_dist_WithCC <- (df_Topt$`With CC` - df_Topt$Final_Topt)/(df_Topt$Final_Tmax-df_Topt$Final_Tmin)


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
# ggplot(data_complete_cw, aes(col = Fmax)) +
#   geom_smooth(aes(x = age, y = mean_of_Madd_per_age, linetype = "Madd")) +
#   geom_smooth(aes(x = age, y = mean_of_Mfor_per_age, linetype = "Mfor")) + 
#   geom_smooth(aes(x = age, y = mean_of_Mpred_per_age, linetype = "Mpred")) + 
#   scale_linetype_manual(values = c("Madd" = 1, "Mfor" = 2, "Mpred" = 3)) +
#   facet_wrap(~Species)

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





# join datasets
# Ne <- Ne_NT_tot_heatmap[which(Ne_NT_tot_heatmap$term == "(Intercept)"), ]
Ne <- Ne_NT_tot_heatmap[, c("Species", "variable", "mean_diff", "Response")]
tf <- dt_tf_heatmap[, c("Species", "variable", "mean_diff")]
Ne <- pivot_wider(Ne, values_from = "mean_diff", names_from = "variable")
tf <- pivot_wider(tf, values_from = "mean_diff", names_from = "variable")

data <- full_join(Ne, df)
data <- full_join(data, tf)
data <- data[, -which(colnames(data) %in% c("Spe_scien", "spi"))] 
data[,-c(1,2)] <- sapply(3:ncol(data), FUN = function(x) as.numeric(as.matrix(data[,x])))


# PLOT summarize

ggplot(data) +
  aes(y = `Nₑ/NT`, size = `Vk.`, colour = Temp_dist_WithCC, x = T) +
  geom_point() +
  geom_text(aes(label=Species), size=4)+
  scale_color_continuous(palette = "heat", trans  = 'reverse')+
  # scale_color_gradient() +
  labs(
    y = "Nₑ/NT",
    size = "Vk.",
    color = "Topt",
    x = "Generation length"
  ) +
  theme_bw()






# PCA ---------------

data_sp <- data
rown <- as.matrix(data[, "Species"])
data <- as.matrix(data[,-c(1, which(is.na(data), arr.ind = T)[,2])])
rownames(data) <- rown

# res.pca <- prcomp(data[,-c(1, which(colnames(data)=="Nₑ/NT"))], 
                         # scale = F)
data_pca <- data_sp[,c(grep(pattern="Final", colnames(data_sp)), 
                       which(colnames(data_sp) %in% c("Starvation mortality rate", "Predation mortality rate")))]
nm <- colnames(data_pca)
data_pca <- as.matrix(as.data.frame(sapply(1:ncol(data_pca), FUN=function(x) log(abs(min(data_pca[, x])+data_pca[, x]+0.01)))))
colnames(data_pca) <- paste0("log(", nm, ")")
rownames(data_pca) <- data_sp$Species
if(any(is.na(data_pca))) data_pca <- data_pca[, -unique(which(is.na(data_pca), arr.ind = T)[,2])]
res.pca <- prcomp(data_pca, scale = F)
data_sp[, "Diff"] <- ifelse(data_sp[, "Nₑ/NT"]>0, "Increase", "Decrease")
fviz_pca_biplot(res.pca, axes=c(1,2), habillage=data_sp$Diff, label.ind =data_sp$Species,
                addEllipses=F, ellipse.level=0.95)+
  theme_minimal()

data_pca <- data_sp[,c(grep(pattern="Final", colnames(data_sp)))]
data_pca <- as.matrix(data_pca)
rownames(data_pca) <- data_sp$Species
if(any(is.na(data_pca))) data_pca <- data_pca[, -unique(which(is.na(data_pca), arr.ind = T)[,2])]
res.pca <- prcomp(data_pca, scale = F)
data_sp[, "Diff"] <- ifelse(data_sp[, "Nₑ/NT"]>0, "Increase", "Decrease")
fviz_pca_biplot(res.pca, axes=c(1,2), habillage=data_sp$Diff, label.ind =data_sp$Species,
                addEllipses=F, ellipse.level=0.95)+
  theme_minimal()


data_pca <- data_sp[,c(grep(pattern="Final", colnames(data_sp)), 
                       which(colnames(data_sp) %in% c("lifespan", "Mnat", "Eggs (scaled for constant Nα)", "m0")))]
data_pca <- as.matrix(data_pca)
rownames(data_pca) <- data_sp$Species
if(any(is.na(data_pca))) data_pca <- data_pca[, -unique(which(is.na(data_pca), arr.ind = T)[,2])]
res.pca <- prcomp(data_pca, scale = F)
data_sp[, "Diff"] <- ifelse(data_sp[, "Nₑ/NT"]>0, "Increase", "Decrease")
fviz_pca_biplot(res.pca, axes=c(1,2), habillage=data_sp$Diff, label.ind =data_sp$Species,
                addEllipses=F, ellipse.level=0.95)+
  theme_minimal()


data_pca <- data_sp[,c(grep(pattern="Final", colnames(data_sp)), 
                       which(colnames(data_sp) %in% c("Eggs (scaled for constant Nα)", )))]
data_pca <- as.matrix(data_pca)
rownames(data_pca) <- data_sp$Species
if(any(is.na(data_pca))) data_pca <- data_pca[, -unique(which(is.na(data_pca), arr.ind = T)[,2])]
res.pca <- prcomp(data_pca, scale = F)
data_sp[, "Diff"] <- ifelse(data_sp[, "Nₑ/NT"]>0, "Increase", "Decrease")
fviz_pca_biplot(res.pca, axes=c(1,2), habillage=data_sp$Diff, label.ind =data_sp$Species,
                addEllipses=F, ellipse.level=0.95)+
  theme_minimal()

tempo <- as.data.frame(sapply(3:ncol(data_sp[, -1]), FUN=function(x) log(abs(min(data_sp[, x])+data_sp[, x]+0.01))))
colnames(tempo) <- paste0("log(", colnames(tempo), ")")
data_sp <- cbind(data_sp, tempo)
lm<-lm(data=(data_sp), formula = `log(Nₑ.NT)` ~ `log(Vk.)`+`log(lx)`)
summary(lm)

pca_var <- as.data.frame(res.pca$x)
pca_var$Species <- rownames(pca_var)

data$Species <- rown
df_Ne_NT <- full_join(pca_var, (data_sp))
df_Ne_NT <- df_Ne_NT[, -which(colnames(df_Ne_NT) == c("Species"))]

df_sp     <- data_sp
df_sp$Species -> rownames(df) 
df_Ne_NT2 <- full_join(df_sp, Ne_NT, by = "Species") 
df_Ne_NT2 <- df_Ne_NT2[, -grep(colnames(df_Ne_NT2), pattern = "_abs")]

lm(data = df_Ne_NT, formula = `Nₑ/NT` ~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + PC11 + PC12 + PC13)

fit <- plsr(`Nₑ/NT` ~ ., data = df_Ne_NT, scale = TRUE, validation = "CV")
coef(fit, ncomp = 5)
VIP(fit, opt.comp = 5)


# GLM --------------
library(glm)
library(gtsummary)
data_glm <- pivot_wider(Ne_NT_tot2, values_from="value", names_from="variable")
data_glm[, -c(1:4)] <- scale(data_glm[, -c(1:4)])
pretty_coefficients(deposit_model, type_iii = 'Wald')
data_glm$Species <- as.factor(data_glm$Species)
data_glm$Scenario <- as.factor(data_glm$Scenario)
glm_s <- glm(data=data_glm, Ne_N_alpha ~ 0 + Species + Scenario + Species:Scenario)
summary(glm_s)
glm_s %>%
  tbl_regression()


# BIPLOT -----------------


# join datasets
# Ne <- Ne_NT_tot_heatmap[which(Ne_NT_tot_heatmap$term == "(Intercept)"), ]


Ne <- Ne_NT_tot_heatmap[, c("Species", "variable", "mean_diff", "Response")]
tf <- dt_tf_heatmap[, c("Species", "variable", "mean_diff")]
Ne <- pivot_wider(Ne, values_from = "mean_diff", names_from = "variable")
tf <- pivot_wider(tf, values_from = "mean_diff", names_from = "variable")
data <- full_join(Ne, df)
data <- full_join(data, tf)
data <- data[, -which(colnames(data) %in% c("Spe_scien", "spi"))] 
data[,-c(1,2)] <- sapply(3:ncol(data), FUN = function(x) as.numeric(as.matrix(data[,x])))

ggplot(data, aes(y=`Nₑ/NT`, x=Final_Topt, group=Response))+
  geom_point(aes(col=(Response)), size=3)+
  theme_bw()+
  geom_smooth(method ="lm")+
  geom_text_repel(aes(label = Species))
ggplot(data, aes(y=`Nₑ/NT`, x=Final_Tmax))+
  geom_point(aes(col=(Response)), size=3)+
  theme_bw()+
  geom_smooth(method ="lm")+
  geom_text_repel(aes(label = Species))
ggplot(data, aes(y=`Nₑ/NT`, x=Final_Tmin, group=Response))+
  geom_point(aes(col=(Response)), size=3)+
  theme_bw()+
  geom_smooth(method ="lm")+
  geom_text_repel(aes(label = Species))

Ne_NT_tot_diff <- Ne_NT_tot2 %>% 
  group_by(Species, replicate, variable, Fmax) %>% 
  mutate(difference=(`With CC` - `No CC`)/abs(`No CC`)) %>%
  ungroup() 
dt_tf_diff <- dt_tf %>% 
  group_by(Species, replicate, variable, Fmax) %>% 
  mutate(difference=(`With CC` - `No CC`)/abs(`No CC`)) %>%
  ungroup() 

Ne <- Ne_NT_tot_diff[, c("Species", "variable", "difference", "replicate")]
tf <- dt_tf_diff[, c("Species", "variable", "difference", "replicate")]
Ne <- pivot_wider(Ne, values_from = "difference", names_from = "variable")
tf <- pivot_wider(tf, values_from = "difference", names_from = "variable")

data <- full_join(Ne, df)
data <- full_join(data, tf)
data <- data[, -which(colnames(data) %in% c("Spe_scien", "spi"))] 
data[,-c(1,2)] <- sapply(3:ncol(data), FUN = function(x) as.numeric(as.matrix(data[,x])))

ggplot(data, aes(y=Ne_NT, x=Final_Topt))+
  geom_point()+
  geom_smooth(method ="lm")+
  geom_text_repel(aes(label = data$Species))




cor(data[,-c(1, 2)])

data_log <- as.matrix(as.data.frame(sapply(3:ncol(data), FUN=function(x) log(abs(min(data[, x])+data[, x]+0.01)))))
colnames(data_log) <- paste0("log(", colnames(data_log), ")")
ggplot(data_log, aes(y=`log(T)`, x=`log(Final_Topt)`, group=data$Response))+
  geom_point(aes(col=data$Response))+
  geom_smooth(method ="lm")+
  geom_text_repel(aes(label = data$Species))


ggplot(data, aes(y=`Nₑ/NT`, x=Final_Topt, group=Response))+
  geom_point(aes(col=(Response)))+
  geom_smooth(method ="lm")+
  geom_text_repel(aes(label = Species))
ggplot(data, aes(y=`Nₑ/NT`, x=Final_Tmin))+
  geom_point(aes(col=(`Nₑ/NT`>0)))+
  geom_smooth(method ="lm")+
  geom_text_repel(aes(label = Species))
ggplot(data, aes(y=`Nₑ/NT`, x=Final_Tmax))+
  geom_point()+
  geom_smooth(method ="lm")+
  geom_text_repel(aes(label = Species))
ggplot(data, aes(y=`Nₑ/NT`, x=log(ratio.min)))+
  geom_point()+
  geom_smooth(method ="lm")+
  geom_text_repel(aes(label = Species))
ggplot(data, aes(y=`Nₑ/NT`, x=log(maint.energy.c_m)))+
  geom_point()+
  geom_smooth(method ="loess")+
  geom_text_repel(aes(label = Species))
ggplot(data, aes(y=`Nₑ/NT`, x=`Δ²`))+
  geom_point()+
  geom_smooth(method ="lm")+
  geom_text_repel(aes(label = Species))
ggplot(data, aes(y=`Nₑ/NT`, x=`SSDɢ (Δ²Dead)/Nα`))+
  geom_point()+
  geom_smooth(method ="lm")+
  geom_text_repel(aes(label = Species))
ggplot(data, aes(y=`SSDɢ (Δ²Dead)/Nα`, x=Final_Topt))+
  geom_point()+
  geom_smooth(method ="lm")+
  geom_text_repel(aes(label = Species))
ggplot(data, aes(y=`SSDɢ (Δ²Dead)`, x=Final_Topt, group=Response))+
  geom_point(aes(col=Response))+
  geom_smooth(method ="lm")+
  geom_text_repel(aes(label = Species))+
  theme_bw()
ggplot(data, aes(y=`T`, x=Final_Topt))+
  geom_point(aes(col=Response))+
  geom_smooth(method ="lm")+
  geom_text_repel(aes(label = Species))+
  theme_bw()
ggplot(data, aes(y=`Vk.`, x=Final_Topt))+
  geom_point(aes(col=Response))+
  geom_smooth(method ="lm")+
  geom_text_repel(aes(label = Species))+
  theme_bw()
ggplot(data, aes(y=`SSDɢ (Δ²Dead)`, x=Final_Topt, group=Response))+
  geom_point(aes(col=Response))+
  geom_smooth(method ="lm")+
  geom_text_repel(aes(label = Species))+
  theme_bw()

ggplot(data, aes(y=`Nₑ/NT`, x=log(selectq90)))+
  geom_point()+
  geom_smooth(method ="lm")+
  geom_text_repel(aes(label = Species))
ggplot(data, aes(y=`Nₑ/NT`, x=Vk.))+
  geom_point()+
  geom_smooth(method ="lm")+
  geom_text_repel(aes(label = Species))
ggplot(data, aes(y=`Nₑ/NT`, x=T))+
  geom_point()+
  geom_smooth(method ="lm")+
  geom_text_repel(aes(label = Species))
ggplot(data, aes(y=`Nₑ/NT`, x=Nₑ))+
  geom_point()+
  geom_smooth(method ="lm")+
  geom_text_repel(aes(label = Species))
ggplot(data, aes(y=`Nₑ/NT`, x=`SSDɢ (Δ²Dead)/Nα`))+
  geom_point()+
  geom_smooth(method ="lm")+
  geom_text_repel(aes(label = Species))


ggplot(data, aes(y=`Nₑ/NT_abs`, x=Final_Topt))+
  geom_point()+
  geom_smooth(method ="loess")+
  geom_text_repel(aes(label = Species))
ggplot(data, aes(y=`NT_abs`, x=Final_Topt))+
  geom_point()+
  geom_smooth(method ="loess")+
  geom_text_repel(aes(label = Species))
ggplot(data, aes(y=`Nₑ/NT_abs`, x=log(ratio.min)))+
  geom_point()+
  geom_smooth(method ="lm")+
  geom_text_repel(aes(label = Species))
ggplot(data, aes(y=`Nₑ/NT_abs`, x=log(maint.energy.c_m)))+
  geom_point()+
  geom_smooth(method ="loess")+
  geom_text_repel(aes(label = Species))
ggplot(data, aes(y=`Nₑ/NT_abs`, x=`Δ²`))+
  geom_point()+
  geom_smooth(method ="lm")+
  geom_text_repel(aes(label = Species))
ggplot(data, aes(y=`Nₑ/NT_abs`, x=`SSDɢ (Δ²Dead)/Nα_abs`))+
  geom_point()+
  geom_smooth(method ="lm")+
  geom_text_repel(aes(label = Species))
ggplot(data, aes(y=`Nₑ/NT_abs`, x=log(selectq90)))+
  geom_point()+
  geom_smooth(method ="lm")+
  geom_text_repel(aes(label = Species))

ggplot(data)+
  geom_point(aes(y=`Nₑ/NT`, x=m0))+
  geom_smooth(method ="lm", col="red", aes(y=`Nₑ/NT`, x=m0))+  
  geom_point(aes(y=`T/Vk`, x=m0))+
  geom_smooth(method ="lm", col="orange", aes(y=`T/Vk`, x=m0))+  
  geom_point(aes(y=T, x=m0))+
  geom_smooth(method ="lm", col="blue", aes(y=T, x=m0))+
  geom_text_repel(aes(label = Species, y=`Nₑ/NT`, x=m0), col="red")+
  geom_text_repel(aes(label = Species, y=`T/Vk`, x=m0), col="orange")+
  geom_text_repel(aes(label = Species, y=T, x=m0), col="blue")

ggplot(data, aes(x=Final_Topt))+
  geom_point(aes(y=`Nₑ/NT`))+
  geom_smooth(method ="lm", col="red", aes(y=`Nₑ/NT`))+  
  geom_point(aes(y=`T/Vk`))+
  geom_smooth(method ="lm", col="orange", aes(y=`T/Vk`))+  
  geom_point(aes(y=T))+
  geom_smooth(method ="lm", col="blue", aes(y=T))+
  geom_text_repel(aes(label = Species, y=`Nₑ/NT`), col="red")+
  geom_text_repel(aes(label = Species, y=`T/Vk`), col="orange")+
  geom_text_repel(aes(label = Species, y=T), col="blue")

ggplot(data, aes(x=lifespan))+
  geom_point(aes(y=`Nₑ/NT`))+
  geom_smooth(method ="lm", col="red", aes(y=`Nₑ/NT`))+  
  geom_point(aes(y=`T/Vk`))+
  geom_smooth(method ="lm", col="orange", aes(y=`T/Vk`))+  
  geom_point(aes(y=T))+
  geom_smooth(method ="lm", col="blue", aes(y=T))+
  geom_text_repel(aes(label = Species, y=`Nₑ/NT`), col="red")+
  geom_text_repel(aes(label = Species, y=`T/Vk`), col="orange")+
  geom_text_repel(aes(label = Species, y=T), col="blue")

ggplot(data, aes(x=predation.ingestion))+
  geom_point(aes(y=`Nₑ/NT`))+
  geom_smooth(method ="lm", col="red", aes(y=`Nₑ/NT`))+  
  geom_point(aes(y=`T/Vk`))+
  geom_smooth(method ="lm", col="orange", aes(y=`T/Vk`))+  
  geom_point(aes(y=T))+
  geom_smooth(method ="lm", col="blue", aes(y=T))+
  geom_text_repel(aes(label = Species, y=`Nₑ/NT`), col="red")+
  geom_text_repel(aes(label = Species, y=`T/Vk`), col="orange")+
  geom_text_repel(aes(label = Species, y=T), col="blue")




ggplot(data, aes(y=`Nₑ/NT`, x=lifespan))+
  geom_point()+
  geom_smooth(method ="lm")+
  geom_text_repel(aes(label = Species))

