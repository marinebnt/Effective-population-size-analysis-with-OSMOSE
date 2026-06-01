
main = "New_gpu_conf.R"
source("osmose-ns-v4.4.1_clim/R/functions/function_plot_trait_age_or_time.R")
source("osmose-ns-v4.4.1_clim/R/functions/source.R")
source("osmose-ns-v4.4.1_clim/R/functions/function_data_traits_summary.R")

base_dir <- "osmose-ns-v4.4.1_clim/output/PLOTS/All_clean_genetic_drift_indicators/INIT/"
folders  <- list.dirs(base_dir, recursive = F)
repo        = "Simu_LONG"
multip_name = "y"
multipliers = c(1,0.5,0.75,1.25,1.5)
ideal_order_factor <- paste0("Fmax=", rep(c(0.5,0.75,1,1.25,1.5), 5),
                             " Fstart=", c(rep(0.5, 5),rep(0.75, 5),rep(1, 5),rep(1.25, 5),rep(1.5, 5)))



trait_dim= c("eggNumber", "abundanceByAge", "biomassByAge", "mortalityByAge", "schoolNumber", "Mpred", "Mstar", "Mnat", "F", "Z", "Mfor", "Dis")


# path_in <- c("C:/Users/mbeneat/Documents/osmose/Chap2/osmose-ns/osmose-ns-v4.4.1_clim/output/test")
# path_in <- c("/home1/datawork/mbeneat/SIMU_CHAP2_FIN")
# list_d <- list.dirs(path = path_in, recursive = F, full.names = T)

F_vec <- list()





# BY YEAR -------------------------

pathouty <- paste0("osmose-ns-v4.4.1_clim/output/PLOTS/Trait_variables_by_year_by_age2/otherplots150160")
if (!dir.exists(pathouty)) dir.create(pathouty, recursive = T)
# load dataset
# side note : when I produce the RData with trait ranges, with files called "PLOT_trait_replicate"-something, I log(fecundity + 10^-8)
data_complete_sp = NULL
data_complete_y = NULL
base_path = "osmose-ns-v4.4.1_clim/output/PLOTS/Trait_variables_by_year_by_age150-160/"
folder_name <- basename(base_path)
suffix <- gsub("Trait_variables_by_", "", folder_name)
suffix <- gsub("_", "_", suffix)  # keep underscores
save_dir <- file.path(dirname(dirname(base_path)), "RDATA")
save_file <- file.path(save_dir, paste0("summary_trait_", suffix, ".RData"))
if (!dir.exists(save_dir)){dir.create(save_dir, recursive = T)}
if (!file.exists(save_file)){
  data_complete_y <- data_traits_summary(base_path = base_path)
} else{ 
  data_complete_y = get(load(file=save_file))
}

# dimensions dataset required : [time, Sp, value, rep, Fish_press, TRAIT]
head(data_complete_y)
unique(data_complete_y$TRAIT)

#adding survival dimension
tempo <- data_complete_y[data_complete_y$TRAIT == "mortality",]
tempo <- pivot_wider(tempo, names_from = "TRAIT", values_from = "value")
tempo$survival <- exp(-tempo$mortality)
tempo$age <- as.numeric(tempo$age)
tempo <- tempo %>%
  arrange(Fish_press, rep, age) %>%
  group_by(Fish_press, rep) %>%
  mutate(
    cumulative_survival = cumsum(survival),
    survival_threshold = cumsum(survival >= 1) > 0,
    cumulative_survival = ifelse(survival_threshold, 
                                 first(cumulative_survival[survival >= 1]), 
                                 cumulative_survival),
    cumulative_survival = cumulative_survival / max(cumulative_survival, na.rm = TRUE)
  )
tempo <- pivot_longer(tempo, cols = "cumulative_survival", names_to = "TRAIT", values_to = "value")
data_complete_y <- rbind(data_complete_y, tempo[, which(colnames(tempo) %in% colnames(data_complete_y))])



# reference data 
data_complete_y <- data_complete_y[-which(is.infinite(data_complete_y$value)),]
F1_data_complete_y <- data_complete_y[which(data_complete_y$Fish_press == c("Fmax=1 Fstart=1")),] %>%
  group_by(TRAIT, Species)%>%
  summarise(F1 = mean(value, na.rm = TRUE))
max_data_complete_y <- data_complete_y %>%
  group_by(TRAIT, Species, Fish_press)%>%
  summarise(max = max(value, na.rm = TRUE))
F1_max <- max_data_complete_y %>%
  left_join(F1_data_complete_y, by = c("Species", "TRAIT"), relationship = "many-to-one")
F1_max$perc_diff <- ((F1_max$max-F1_max$F1)*100)/F1_max$F1
summary(F1_max[which(F1_max$TRAIT == "biomassByAge"),"perc_diff"])
summary(F1_max[which(F1_max$TRAIT == "sizeMature"),"perc_diff"])
summary(F1_max[which(F1_max$TRAIT == "sizeMature"),"perc_diff"])






##################################################################################################################
#### First version of the estimator of TL with mortality by species : plotting each species mortality in // ######
####
library(RColorBrewer)
n <- 15
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
pie(rep(1,n), col=sample(col_vector, n))

# Plot the mortality by species over time
mort = NULL
for (Fi in unique(data_complete_y$Fish_press)){
  data_F <- data_complete_y[which(data_complete_y$Fish_press == Fi),]
  data_F <- data_F[-which(data_F$age == 0),]
  data_F$age   <- as.numeric(data_F$age)
  data_F$value <- as.numeric(data_F$value)
  data_F = data_F[which(data_F$TRAIT == "MpredbyY"),]
  data_F <- data_F%>%
    group_by(Species, age, TRAIT)%>%
    summarise(mean=mean(value, na.rm=T),
           sd=sd(value, na.rm=T))
  # data_F <- data_F%>%
  #   group_by(age, TRAIT)%>%
  #   mutate(rela = 100*(mean - mean(mean, na.rm=T))/mean(mean, na.rm=T))
  # data_F[which(data_F$mean == 0),"mean"] <- NA
  mort[[Fi]] <- ggplot(data_F,
                 aes(x=age, y=mean, colour=Species))+
                geom_line(lwd=0.55)+
                theme_bw()+
                scale_color_manual(values = sample(col_vector, n))+
                ylab(label = paste("predation mortality rate", collapse="\n"))+
                xlab(label = "Age")+
                labs(title = paste0("Relative predation mortality by species by age for F = ", Fi))
}

pdf(file = paste0("osmose-ns-v4.4.1_clim/output/PLOTS/mortalitybyspecies.pdf"), width = 10, height = 8)
print(mort)
dev.off()

# other version : try to collect only the last mortlity rate before death
last_Mpred <- data.frame(matrix(NA,
                                ncol=length(unique(data_complete_y$Fish_press)),
                                nrow = length(species)))
colnames(last_Mpred) <- unique(data_complete_y$Fish_press)
rownames(last_Mpred) <- species
for (Fi in unique(data_complete_y$Fish_press)){
  Fi_chr <- as.character(Fi)
  data_F <- data_complete_y[which(data_complete_y$Fish_press == Fi),]
  data_F <- data_F[-which(data_F$age == 0),]
  data_F$age   <- as.numeric(data_F$age)
  data_F$value <- as.numeric(data_F$value)
  data_F = data_F[which(data_F$TRAIT == "MpredbyY"),]
  data_F <- data_F%>%
    group_by(Species, age, TRAIT)%>%
    summarise(mean=mean(value, na.rm=T),
              sd=sd(value, na.rm=T))
  last_Mpred[, Fi_chr] <- unlist(sapply(species, FUN = function(x) {
    data_F_sp <- data_F[which(data_F$Species == x),]
    data_F_sp <- data_F_sp[-which(data_F_sp$mean == 0),]
    return(data_F_sp[which.max(data_F_sp$age),"mean"])
  }))
}

####
##################################