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
pathout=pathoutputRDATA, Time2, Time2)

## loading the dataset of traits averaged over cohorts Time1 to Time2, and reorganised to remove age=0, 
# add reassign natural mortality as add mortality, and estimate natural mortality
data_variables <- load_trait_df(base_path = base_path, pathout=pathoutputRDATA, FORCE)

CLEAN = T
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

# _  ------------------
# 
# 
######## PLOT -0 selectivity curves --------------
mat_sp_i <- as.data.frame(as.matrix(get_par(conf, "species.name")))
mat_sp_i2 <- as.data.frame(as.matrix(get_par(conf, "inf")))
rownames(mat_sp_i) <- stringr::str_remove(rownames(mat_sp_i), pattern="species.name.sp")
rownames(mat_sp_i2) <- stringr::str_remove(rownames(mat_sp_i2), pattern="species.linf.sp")
surv_curve = list()
for (sp in species[order(species)]){
  i <- rownames(mat_sp_i)[which(mat_sp_i == sp)]
  df_select <- NULL
  data_surv_curve <- do.call(what = cbind, lapply(c(0.5, 0.75, 1, 1.25, 1.5), FUN = function(x) {
    file_i <- paste0("new_config_Fmax_Fstart/F_Fstart1_Fmax", x, "_sp", i, ".csv")
    df_select_i <- read.csv(file_i, sep = ",", dec = ".", check.names = FALSE)
    df_select_ii <- as.data.frame(t(df_select_i[1,]))
    colnames(df_select_ii) <- x
    if (is.null(df_select)) {df_select <- df_select_ii
    } else df_select <- cbind(df_select_ii, df_select)
  }))
  df_select <- as.data.frame(data_surv_curve)
  for (i in c(1:ncol(data_surv_curve))) {df_select[,i] <-  data_surv_curve[,i]/max(data_surv_curve[,c("1")])}
  data_surv_curve <- df_select
  data_surv_curve[which(apply(df_select,2,is.nan), arr.ind = T)] <- 0
  data_surv_curve$length <- as.numeric(rownames(data_surv_curve))
  if (max(data_surv_curve$length)<mat_sp_i2[as.numeric(which(rownames(mat_sp_i2) == i)),1]){
    dlength <- data_surv_curve$length[2] - data_surv_curve$length[1]
    vec <- seq(max(data_surv_curve$length)+dlength,as.numeric(mat_sp_i2[as.numeric(which(rownames(mat_sp_i2) == i)),1]), by=dlength)
    mat_supp <- as.data.frame(do.call(rbind, lapply(vec, FUN = function(i) as.matrix(data_surv_curve[nrow(data_surv_curve),-which(colnames(data_surv_curve) == c("length"))]))))
    mat_supp$length <- as.matrix(vec)
    data_surv_curve <- rbind(data_surv_curve, mat_supp)
  }
  data_surv_curve <- pivot_longer(data_surv_curve, names_to = "θF", cols = -c("length"), values_to = "F_sc")
  data_surv_curve$Species <- improve_sp_names(sp)
  surv_curve[[sp]] <- ggplot(data_surv_curve, aes(x=length, y=F_sc, group=θF)) +
    geom_line(aes(linetype = θF), linewidth=1.3)+
    labs(title="", y="Fishing mortality rate, F", x="Body length", linetype="θF")+ #, y="SSD\u0262"
    theme(strip.background = element_rect(fill="white"),
          panel.background = element_rect(fill = "white", colour = "grey50"),
          panel.text = element_text(size = 12),
          plot.title = element_text(hjust = 0.5, size=12),
          axis.text.x = element_text(angle = 60, hjust = 1),
          axis.text = element_text(colour = 1, size = 12),
          axis.title=element_text(size=14,face="bold"),
          panel.border = element_rect(colour = "black", fill=NA),
          strip.text = element_text(size = 12),
          legend.key.size = unit(1.1, "cm"),
          legend.text=element_text(size=9))+
    facet_wrap(~Species) +                # nom de l'espèce en bandeau
    labs(title = NULL, subtitle = NULL) +
    theme(strip.background = element_blank(),
          strip.text = element_text(size = 12),
          axis.title = element_blank())   # titres d'axes communs, ajoutés plus bas

}


names(surv_curve) <- species_imp
lab_left  <- surv_curve[[1]]$labels$y
lab_right <- tryCatch(panels[[1]]$scales$get_scales("y")$secondary.axis$name,
                      error = function(e) "")
if (is.null(lab_right) || inherits(lab_right, "waiver")) lab_right <- ""
txt <- function(label, angle) ggplot() +
  annotate("text", x = 0, y = 0, label = label, angle = angle, size = 4.2) + theme_void()
x_lab <- ggplot() + annotate("text", x = 0, y = 0, label = paste("Length (cm)"), parse = TRUE,
                             size = 4.2) + theme_void()
per_page <- 15
pages <- split(seq_along(surv_curve), ceiling(seq_along(surv_curve) / per_page))
p_list <- lapply(pages, function(idx) {
  grid <- wrap_plots(surv_curve[idx], ncol = 3, nrow = 5)
  (((txt(lab_left, 90) | grid | txt(lab_right, -90)) +
      plot_layout(widths = c(0.035, 1, 0.035))) /
      x_lab + plot_layout(heights = c(1, 0.025))) +
    plot_layout(guides = "collect") &
    theme(legend.position = "top", legend.direction = "horizontal",
          legend.title = element_blank(),
          legend.box.background = element_rect(colour = "black"))
})

windowsFonts(myCustomWindowsFontName=windowsFont("DejaVu Sans"))
cairo_pdf(file = paste0(PATHOUT, "/Clean_for_paper/Surv_Ne_NT_with_F.pdf"),
          width = 10, height =13, onefile=T)
print(p_list)
dev.off()


#### PLOT 0- Biomass vs Fmax & keep it on the side -------------------
variables_to_process <- c(NT = "Abund", Ne_N_alpha = "Ne")

#comparing with and without climate change biomass
biomassplot = NULL
rel_biomassplot = list()
NTplot = list()
biomassplot_curve = list()
abundanceplot_curve = list()
for (var in names(variables_to_process)){
  list_relative_biomass <- list()
  list_long_for_delta <- list()
  for (sp in species){

    d <- data_parameters

    id_biom <- c(grep(pattern=var, d$variable))
    d_i_i <- d[id_biom, ]
    d_i_i <- d_i_i[which(d_i_i$Species == sp), ]
    list_long_for_delta[[sp]] <- d_i_i[which(d_i_i$variable == var), ]
    d_i_i <- pivot_wider(d_i_i, values_from = "value", names_from = "variable")
    d_i_i$Iszero <- ifelse(d_i_i[[var]] < 100, "Abundance <100", "Abundance >100")

    NTplot[[var]][[sp]] <- ggplot(d_i_i, aes(y=.data[[var]], x=Fmax)) +
      geom_point(aes(colour = Iszero))+
      labs(title = "Abundance per fishing pressure", subtitle = sp, y = var)

    ref_data_percent <- d_i_i %>%
      filter(Fmax == 1, Fstart == 1) %>%
      dplyr::select(replicate, all_of(var)) %>%
      dplyr::rename(ref_value_percent_abund = all_of(var))
    # Subtract reference from all means
    rel_d_i_i_perc <- d_i_i %>%
      left_join(ref_data_percent) %>%
      mutate(relative_abundance_percent = 100*.data[[var]]/ref_value_percent_abund-100)
    list_relative_biomass[[sp]] <- rel_d_i_i_perc
    rel_biomassplot[[var]][[sp]] <- ggplot(rel_d_i_i_perc, aes(y=relative_abundance_percent, x=Fmax)) +
      geom_point(aes(colour = Iszero))+
      labs(title = "Relative biomass per fishing pressure", subtitle = sp)
  }

  ## Estimate the % in decline (in biomass and abundance, at first age at reproduction) from Fmax=1.5 to Fmax = 0.5
  df_relative_N <- do.call(rbind, list_relative_biomass)
  summary_df_relative_abunance <- df_relative_N %>%
    group_by(Fmax, Species)%>%
    summarise(across(relative_abundance_percent, list(mean=mean, sd=sd), na.rm=TRUE),
              .groups = "drop")
  df_long_for_delta <- do.call(rbind, list_long_for_delta)
  delta_res <- value_to_sd_esperance_delta_method(df_long_for_delta, var_x = "Fmax") %>%
    filter(Fstart == 1) %>%   # same reference as the empirical estimate (Fmax == 1, Fstart == 1)
    mutate(relative_abundance_percent_delta_mean = 100*mean,   # fraction -> %
           relative_abundance_percent_delta_sd   = 100*sd) %>%
    dplyr::select(Fmax, Species, relative_abundance_percent_delta_mean, relative_abundance_percent_delta_sd)
  summary_df_relative_abunance <- summary_df_relative_abunance %>%
    left_join(delta_res, by = c("Fmax", "Species"))
  summary_df_relative_abunance$Species <- improve_sp_names(summary_df_relative_abunance$Species)
  summary_df_relative_abunance <- summary_df_relative_abunance %>%
    dplyr::mutate(percent_difference_withFmax_delta=round(relative_abundance_percent_delta_mean, 2)) %>%   # >>> DELTA
    dplyr::mutate(percent_difference_withFmax_delta_sd=round(relative_abundance_percent_delta_sd, 2)) %>%  # >>> DELTA
    select(Fmax, Species, #percent_difference_withFmax, percent_difference_withFmax_sd,
           percent_difference_withFmax_delta, percent_difference_withFmax_delta_sd)
  summary_df_relative_abunance <- pivot_wider(summary_df_relative_abunance,
                                              values_from = c(#"percent_difference_withFmax", "percent_difference_withFmax_sd",
                                                              "percent_difference_withFmax_delta", "percent_difference_withFmax_delta_sd"),
                                              names_from = "Species")
  colnames(summary_df_relative_abunance) <- stringr::str_remove(colnames(summary_df_relative_abunance), pattern = "percent_difference_withFmax_delta_")
  write.csv(summary_df_relative_abunance,
            file=paste0(PATHOUT, "/Clean_for_paper/", variables_to_process[[var]], "_declinetrends.csv"))

}


#### PLOT 1- NE/N VS FMAX + Biomass VS FMAX ------------------
if (CLEAN){
Ne_plot_list = list()
traits_to_plot <- list(c("generation_length", "Vk_point", "T/Vk."),
                       c("generation_length", "Vk_point+2"),
                       c("log(T/(Vk.+2))", "log(Numerator_N_alpha\nDivNT)", "log(Ne_N_alpha\nDivNT)"),
                       c("Ne_N_alpha", "NT", "Ne_N_alpha\nDivNT"), 
                       c("colSumsurvivorship_curve", "T/Vk."),
                       c("colSumSSDgx", "Numerator_N_alpha", "colSumSSDgx\nDivNumerator_N_alpha"), 
                       c("Numerator_N_alpha", "colSumDeltasquare", "colSumSSDgx")
                       )

# select data_parameters
data_plot_1 <- as_long_if_needed(data_parameters, columns_keep_pivot = 
                          colnames(data_parameters)[-which(colnames(data_parameters) %in% c("replicate", "Fmax", "Fstart", "Species", "age"))])

# clean dataset
data_plot_1 <- data_plot_1 %>% filter(variable %in% c(unlist(traits_to_plot), "SSDgx", "survivorship_curve", "Deltasquare", "Dead", "Mtot"))

# add the Ne/NT
df_trait_t    <- estimateDiv(trait_data_full_i = data_plot_1, var=paste0("Ne_N_alpha\nDivNT"))
data_plot_1       <- rbind(data_plot_1, df_trait_t)

# add the Ne/Nalpha
df_trait_t    <- estimateDiv(trait_data_full_i = data_plot_1, var=paste0("Ne_N_alpha\nDivNumerator_N_alpha"))
data_plot_1       <- rbind(data_plot_1, df_trait_t)

# add the Numerator_N_alpha/NT
df_trait_t    <- estimateDiv(trait_data_full_i = data_plot_1, var=paste0("Numerator_N_alpha\nDivNT"))
data_plot_1       <- rbind(data_plot_1, df_trait_t)

# add the SSDgx/Nalpha
df_trait_t    <- estimateDiv(trait_data_full_i = data_plot_1, var=paste0("SSDgx\nDivNumerator_N_alpha"))
data_plot_1       <- rbind(data_plot_1, df_trait_t)

# Vk.+2
df_trait_t <- data_plot_1[which(data_plot_1$variable == "Vk_point"), ]
df_trait_t$variable <- "Vk_point+2"
df_trait_t$value <- df_trait_t$value+2
data_plot_1 <- rbind(data_plot_1, df_trait_t)

# log(Nalpha/NT)
df_trait_t <- data_plot_1[which(data_plot_1$variable == "Numerator_N_alpha\nDivNT"), ]
df_trait_t$variable <- "log(Numerator_N_alpha\nDivNT)" 
df_trait_t$value <- log(df_trait_t$value)
data_plot_1 <- rbind(data_plot_1, df_trait_t)

# log(Ne/NT)
df_trait_t <- data_plot_1[which(data_plot_1$variable == "Ne_N_alpha\nDivNT"), ]
df_trait_t$variable <- "log(Ne_N_alpha\nDivNT)" 
df_trait_t$value <- log(df_trait_t$value)
data_plot_1 <- rbind(data_plot_1, df_trait_t)

# add the SSDgx/Nalpha
df_trait_t    <- estimateDiv(trait_data_full_i = data_plot_1, var=paste0("generation_length\nDivVk_point"))
df_trait_t$variable <-c("T/Vk.")
data_plot_1       <- rbind(data_plot_1, df_trait_t)

# add the SSDgx/Nalpha
df_trait_t    <- estimateDiv(trait_data_full_i = data_plot_1, var=paste0("generation_length\nDivVk_point+2"))
df_trait_t$variable <-c("T/(Vk.+2)")
data_plot_1       <- rbind(data_plot_1, df_trait_t)

# log(T/(Vk+2))
df_trait_t <- data_plot_1[which(data_plot_1$variable == "T/(Vk.+2)"), ]
df_trait_t$variable <- "log(T/(Vk.+2))" 
df_trait_t$value <- log(df_trait_t$value)
data_plot_1 <- rbind(data_plot_1, df_trait_t)

# colSum survivorship curve
df_trait_t    <- estimatecolSum(trait_data_full_i = data_plot_1, var=paste0("colSumsurvivorship_curve"))
min_age_mat   <- (data_plot_1[sapply(species, function(x) {which(data_plot_1$Species==x)[1]}), c("min_age", "Species")])
df_trait_t    <- full_join(df_trait_t, min_age_mat)
df_trait_t$age <- df_trait_t$min_age
data_plot_1   <- rbind(data_plot_1, df_trait_t)

# colSum SSDgx
df_trait_t    <- estimatecolSum(trait_data_full_i = data_plot_1, var=paste0("colSumSSDgx"))
min_age_mat   <- (data_plot_1[sapply(species, function(x) {which(data_plot_1$Species==x)[1]}), c("min_age", "Species")])
df_trait_t    <- full_join(df_trait_t, min_age_mat)
df_trait_t$age <- df_trait_t$min_age
data_plot_1   <- rbind(data_plot_1, df_trait_t)

# colSum SSDgx/Napha
df_trait_t    <- estimatecolSum(trait_data_full_i = data_plot_1, var=paste0("colSumSSDgx\nDivNumerator_N_alpha"))
min_age_mat   <- (data_plot_1[sapply(species, function(x) {which(data_plot_1$Species==x)[1]}), c("min_age", "Species")])
df_trait_t    <- full_join(df_trait_t, min_age_mat)
df_trait_t$age<- df_trait_t$min_age
data_plot_1   <- rbind(data_plot_1, df_trait_t)

# colSum Dead
df_trait_t    <- estimatecolSum(trait_data_full_i = data_plot_1, var=paste0("colSumDead"))
min_age_mat   <- (data_plot_1[sapply(species, function(x) {which(data_plot_1$Species==x)[1]}), c("min_age", "Species")])
df_trait_t    <- full_join(df_trait_t, min_age_mat)
df_trait_t$age<- df_trait_t$min_age
data_plot_1   <- rbind(data_plot_1, df_trait_t)

# colSum Mortality rate
df_trait_t    <- estimatecolSum(trait_data_full_i = data_plot_1, var=paste0("colSumMtot"))
min_age_mat   <- (data_plot_1[sapply(species, function(x) {which(data_plot_1$Species==x)[1]}), c("min_age", "Species")])
df_trait_t    <- full_join(df_trait_t, min_age_mat)
df_trait_t$age<- df_trait_t$min_age
data_plot_1   <- rbind(data_plot_1, df_trait_t)

# colSum Deltasquare
df_trait_t    <- estimatecolSum(trait_data_full_i = data_plot_1, var=paste0("colSumDeltasquare"))
min_age_mat   <- (data_plot_1[sapply(species, function(x) {which(data_plot_1$Species==x)[1]}), c("min_age", "Species")])
df_trait_t    <- full_join(df_trait_t, min_age_mat)
df_trait_t$age<- df_trait_t$min_age
data_plot_1   <- rbind(data_plot_1, df_trait_t)

data_plot_1 <- data_plot_1[which(data_plot_1$variable %in% c(unlist(traits_to_plot), names(table(data_plot_1$variable))[grep(names(table(data_plot_1$variable)), pattern="colSum")])), ]
data_plot_1 <- data_plot_1[which(data_plot_1$age==data_plot_1$min_age),]
# condition to stop
if (nrow(data_plot_1)<1) next

d_residuals <- pivot_wider(data_plot_1, values_from = "value", names_from = "variable")
d_residuals$`log(Ne_N_alpha\nDivNT)`-log(4)-d_residuals$`Ne_N_alpha\nDivNT`-d_residuals$`log(T/(Vk.+2))`


slopes_list <- list()
#loop to plot relative and absolute data
for (relative_to_ref in c(TRUE)){

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
    df_trait_t <- data_plot_1[which(data_plot_1$variable %in% traits_selected),]

    # clean dataset
    if (sum(is.nan(df_trait_t$value))>0) df_trait_t <- df_trait_t[-which(is.nan(df_trait_t$value)), ]
    if (sum(is.infinite(df_trait_t$value))>0) df_trait_t <- df_trait_t[-which(is.infinite(df_trait_t$value)), ]

    # conditions to stop
    if (nrow(df_trait_t)==0) next
    if (all(is.na(df_trait_t$value))) {
      Ne_plot_list[[paste0(as.character(relative_to_ref), as.character(proportional))]][[paste(species, collapse="")]][[traits_selected]] <- NULL
      next}
    
    # improve aesthetics
    df_trait_t$variable <- improve_vriable_names(vector = df_trait_t$variable)
    t <- improve_vriable_names(traits_selected)
    df_trait_t$variable <- factor(df_trait_t$variable, levels=t)
    df_trait_t$Species <- improve_sp_names(df_trait_t$Species)
    
    # scale_z if >2 var
    scale_z_i <- ifelse(length(t)>2, T, F)
    
    # Use the function to plot
    if (any(grep("log", t))) {
      proportional = T
      relative_to_ref = F
      cat("\n=================> LOG <=============================\n")
    } else {
      proportional = F
      relative_to_ref = T
    }
    sp_main_imp_i <- sp_main_imp
    if (any(grepl(traits_selected, pattern = "SSD") | grepl(traits_selected, pattern = "Delta"))) sp_main_imp_i <- stringr::str_replace(sp_main_imp_i, pattern = "Horse mackerel", "Grey gurnard")
    
    df_trait_t$Species = factor(levels = improve_sp_names(species[order(species)]), df_trait_t$Species)
    result <- plot_trait_vs_fmax(df_trait_sp = df_trait_t, trait_name = t, age_dim = F, spec=species_imp[order(species_imp)],
                                 Time1, Time2, relative_to_ref = relative_to_ref, proportional=proportional,
                                 position_legend = "right", per_species_scale = TRUE, facet = T, subtitle=F, legend_justif="right", 
                                 nrow=5, ncol=3, scale_z = scale_z_i)
    slopes_list[[paste(traits_selected, collapse = "")]] <- result$data %>%
      dplyr::filter(!is.na(mean)) %>%                        # Fmax effondrés déjà en NA -> exclus
      dplyr::group_by(Species, variable) %>%
      dplyr::filter(dplyr::n() >= 3) %>%                     # pente sur au moins 3 valeurs de thetaF
      dplyr::summarise(
        slope    = coef(lm(mean ~ Fmax))[2],
        slope_se = summary(lm(mean ~ Fmax))$coefficients[2, 2],
        p_value  = summary(lm(mean ~ Fmax))$coefficients[2, 4],
        n_F      = dplyr::n(),
        n_F_ci_excl0 = sum(Fmax != 1 & abs(mean) > IC_95, na.rm = TRUE),   # nb de thetaF dont l'IC exclut 0
        n_F_tested   = sum(Fmax != 1),
        .groups  = "drop") %>%
      dplyr::mutate(response = ifelse(n_F_ci_excl0 >= 2, "responsive", "constant"))

    
    
    
    df_trait_t2 <- df_trait_t[which(df_trait_t$Species %in% sp_main_imp_i), ]
    df_trait_t2$Species = factor(levels = sp_main_imp_i, df_trait_t2$Species)
    result2 <- plot_trait_vs_fmax(df_trait_sp = df_trait_t2, trait_name = t, age_dim = F, spec=sp_main_imp_i,
                                 Time1, Time2, relative_to_ref = relative_to_ref, proportional=proportional,
                                 position_legend = "right", facet = T, subtitle=F, 
                                 nrow=2, ncol=2, per_species_scale = TRUE, scale_z = scale_z_i, delta_method = T)
    
    # save in a list
      Ne_plot_list[[paste0(as.character(relative_to_ref), as.character(proportional))]][[paste(species, collapse="")]][[paste(traits_selected, collapse="")]] <- result$plot
      Ne_plot_list[[paste0(as.character(relative_to_ref), as.character(proportional))]][[paste(sp_main_imp, collapse="")]][[paste(traits_selected, collapse="")]] <- result2$plot
  } # trait selected in the list
} # relative ? 
save(Ne_plot_list, file=paste0(pathoutputRDATA, "/Ne_plots.RData"))
load(file=paste0(pathoutputRDATA, "/Ne_plots.RData"))

d_residuals <- pivot_wider(data_plot_1, values_from = "value", names_from = "variable")
residuals <- d_residuals$`log(Ne_N_alpha\nDivNT)`-log(4)-d_residuals$`log(Numerator_N_alpha\nDivNT)`-d_residuals$`log(T/(Vk.+2))`

for (relative_to_ref in c(TRUE)){
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
  cairo_pdf(file = paste0(PATHOUT, "/Clean_for_paper/Ne_NT_NeonNT", plot_name_ch, ".pdf"),
            width = 7, height =6, onefile=T)
  print(Ne_plot_list[[plot_name]]$`WhitingEuropean plaiceAtlantic mackerelHorse mackerel`$`Ne_N_alphaNTNe_N_alpha\nDivNT`)
  dev.off() 
  
  windowsFonts(myCustomWindowsFontName=windowsFont("DejaVu Sans"))
  cairo_pdf(file = paste0(PATHOUT, "/Clean_for_paper/TOTNe_NT_NeonNT", plot_name_ch, ".pdf"),
            width = 10, height =12, onefile=T)
  print(Ne_plot_list[[plot_name]]$AtlanticHerringAtlanticMackerelSandeelEuropeanSpratNorwayPoutEuropeanPlaiceCommonSoleSaitheAtlanticCodHaddockHorseMackerelWhitingCommonDabGreyGurnardHake$`Ne_N_alphaNTNe_N_alpha\nDivNT`)
  dev.off()   
  
  windowsFonts(myCustomWindowsFontName=windowsFont("DejaVu Sans"))
  cairo_pdf(file = paste0(PATHOUT, "/Clean_for_paper/SSDgx_Numerator_Nalpha", plot_name_ch, ".pdf"),
            width = 7, height =6, onefile=T)
  print(Ne_plot_list[[plot_name]]$`WhitingEuropean plaiceAtlantic mackerelHorse mackerel`$`colSumSSDgxNumerator_N_alphacolSumSSDgx\nDivNumerator_N_alpha`)
  dev.off() 
  
  windowsFonts(myCustomWindowsFontName=windowsFont("DejaVu Sans"))
  cairo_pdf(file = paste0(PATHOUT, "/Clean_for_paper/TOTSSDgx_Numerator_Nalpha", plot_name_ch, ".pdf"),
            width = 10, height =12, onefile=T)
  print(Ne_plot_list[[plot_name]]$AtlanticHerringAtlanticMackerelSandeelEuropeanSpratNorwayPoutEuropeanPlaiceCommonSoleSaitheAtlanticCodHaddockHorseMackerelWhitingCommonDabGreyGurnardHake$`colSumSSDgxNumerator_N_alphacolSumSSDgx\nDivNumerator_N_alpha`)
  dev.off()  
  
  windowsFonts(myCustomWindowsFontName=windowsFont("DejaVu Sans"))
  cairo_pdf(file = paste0(PATHOUT, "/Clean_for_paper/Vk_T_TonVk", plot_name_ch, ".pdf"),
            width = 7, height =6, onefile=T)
  print(Ne_plot_list[[plot_name]]$`WhitingEuropean plaiceAtlantic mackerelHorse mackerel`$`generation_lengthVk_pointT/Vk.`)
  dev.off()
  
  windowsFonts(myCustomWindowsFontName=windowsFont("DejaVu Sans"))
  cairo_pdf(file = paste0(PATHOUT, "/Clean_for_paper/TOTVk_T_TonVk", plot_name_ch, ".pdf"),
            width = 10, height =12, onefile=T)
  print(Ne_plot_list[[plot_name]]$AtlanticHerringAtlanticMackerelSandeelEuropeanSpratNorwayPoutEuropeanPlaiceCommonSoleSaitheAtlanticCodHaddockHorseMackerelWhitingCommonDabGreyGurnardHake$`generation_lengthVk_pointT/Vk.`)
  dev.off()  
  
  windowsFonts(myCustomWindowsFontName=windowsFont("DejaVu Sans"))
  cairo_pdf(file = paste0(PATHOUT, "/Clean_for_paper/TVk+2", plot_name_ch, ".pdf"),
            width = 7, height =6, onefile=T)
  print(Ne_plot_list[[plot_name]]$`WhitingEuropean plaiceAtlantic mackerelHorse mackerel`$`generation_lengthVk_point+2`)
  dev.off()
  
  windowsFonts(myCustomWindowsFontName=windowsFont("DejaVu Sans"))
  cairo_pdf(file = paste0(PATHOUT, "/Clean_for_paper/TOTTVk+2", plot_name_ch, ".pdf"),
            width = 10, height =12, onefile=T)
  print(Ne_plot_list[[plot_name]]$AtlanticHerringAtlanticMackerelSandeelEuropeanSpratNorwayPoutEuropeanPlaiceCommonSoleSaitheAtlanticCodHaddockHorseMackerelWhitingCommonDabGreyGurnardHake$`generation_lengthVk_point+2`)
  dev.off()  
  
  windowsFonts(myCustomWindowsFontName=windowsFont("DejaVu Sans"))
  cairo_pdf(file = paste0(PATHOUT, "/Clean_for_paper/colSumsurvivorship_curve", plot_name_ch, ".pdf"),
            width = 7, height =6, onefile=T)
  print(Ne_plot_list[[plot_name]]$`WhitingEuropean plaiceAtlantic mackerelHorse mackerel`$`colSumsurvivorship_curveT/Vk.`)
  dev.off()
  
  windowsFonts(myCustomWindowsFontName=windowsFont("DejaVu Sans"))
  cairo_pdf(file = paste0(PATHOUT, "/Clean_for_paper/TOTcolSumsurvivorship_curve", plot_name_ch, ".pdf"),
            width = 10, height =12, onefile=T)
  print(Ne_plot_list[[plot_name]]$AtlanticHerringAtlanticMackerelSandeelEuropeanSpratNorwayPoutEuropeanPlaiceCommonSoleSaitheAtlanticCodHaddockHorseMackerelWhitingCommonDabGreyGurnardHake$`colSumsurvivorship_curveT/Vk.`)
  dev.off()    
  
  windowsFonts(myCustomWindowsFontName=windowsFont("DejaVu Sans"))
  cairo_pdf(file = paste0(PATHOUT, "/Clean_for_paper/colSumDeadDeltaM", plot_name_ch, ".pdf"),
            width = 7, height =6, onefile=T)
  print(Ne_plot_list[[plot_name]]$`WhitingEuropean plaiceAtlantic mackerelHorse mackerel`$Numerator_N_alphacolSumDeltasquarecolSumSSDgx)
  dev.off()
  
  windowsFonts(myCustomWindowsFontName=windowsFont("DejaVu Sans"))
  cairo_pdf(file = paste0(PATHOUT, "/Clean_for_paper/TOTcolSumDeadDeltaM", plot_name_ch, ".pdf"),
            width = 10, height =12, onefile=T)
  print(Ne_plot_list[[plot_name]]$AtlanticHerringAtlanticMackerelSandeelEuropeanSpratNorwayPoutEuropeanPlaiceCommonSoleSaitheAtlanticCodHaddockHorseMackerelWhitingCommonDabGreyGurnardHake$Numerator_N_alphacolSumDeltasquarecolSumSSDgx)
  dev.off()  
}
for (relative_to_ref in c("proportional")){
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
  cairo_pdf(file = paste0(PATHOUT, "/Clean_for_paper/logNalpha_NTlogT_Vk2", plot_name_ch, ".pdf"),
            width = 7, height =6, onefile=T)
  print(Ne_plot_list[[plot_name]]$`WhitingEuropean plaiceAtlantic mackerelHorse mackerel`$`log(T/(Vk.+2))log(Numerator_N_alpha\nDivNT)log(Ne_N_alpha\nDivNT`)
  dev.off()
  
  windowsFonts(myCustomWindowsFontName=windowsFont("DejaVu Sans"))
  cairo_pdf(file = paste0(PATHOUT, "/Clean_for_paper/TOTlogNalpha_NTlogT_Vk2", plot_name_ch, ".pdf"),
            width = 10, height =12, onefile=T)
  print(Ne_plot_list[[plot_name]]$AtlanticHerringAtlanticMackerelSandeelEuropeanSpratNorwayPoutEuropeanPlaiceCommonSoleSaitheAtlanticCodHaddockHorseMackerelWhitingCommonDabGreyGurnardHake$`log(T/(Vk.+2))log(Numerator_N_alpha\nDivNT)log(Ne_N_alpha\nDivNT`)
  dev.off()
}

# from the plot loop, we have extracted the slopes, and now plot it
slopes_all <- dplyr::bind_rows(slopes_list) %>%
  dplyr::distinct(Species, variable, .keep_all = TRUE) %>%   # une variable peut apparaître dans plusieurs jeux de traits
  dplyr::group_by(variable) %>%
  dplyr::mutate(slope_std = pmax(pmin(slope / quantile(abs(slope), 0.8, na.rm = TRUE), 1), -1)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(signif = cut(p_value, c(-Inf, 0.01, 0.05, 0.1, Inf), c("***", "**", "*", "")))
write.csv(slopes_all, file = paste0(PATHOUT, "/slopes_relative_response_vs_thetaF1.csv"), row.names = FALSE)


slopes_order_chosen_var <- slopes_all %>%
  filter(variable %in% c("Nₑ", "NT", "Nₑ/NT", "T", "Vk.", "T/Vk.", "Lτ", 
                       "SSDɢτ", "SSDɢτ/Nα", "Nα"))
slopes_order_chosen_var$slope <- round(slopes_order_chosen_var$slope, digits = 2)
slopes_order_chosen_var$variable <- factor(slopes_order_chosen_var$variable, 
                levels=c("NT", "Nₑ", "Nₑ/NT", "T", "Vk.", "T/Vk.", "Lτ",
                         "SSDɢτ", "Nα", "SSDɢτ/Nα"))
slopes_order_chosen_var$Species <- factor(slopes_order_chosen_var$Species, 
                                  levels=as.matrix(slopes_order_chosen_var[which(
                                    slopes_order_chosen_var[,"variable"]=="Nₑ/NT"),"Species"])[order(as.matrix(
                                    slopes_order_chosen_var[which(slopes_order_chosen_var$variable=="Nₑ/NT"),"slope"]))])

p_heat <- ggplot(slopes_order_chosen_var, aes(x = variable, y = Species, fill = slope_std)) +
  geom_tile(colour = "white", linewidth = 0.6) +
  geom_text(aes(label = paste0(signif(slope, 2), ifelse(response=="responsive", "***", ""))), size = 3, colour = "grey15", fontface="bold") +
  scale_fill_gradient2(low = "#2166AC", mid = "grey97", high = "brown", midpoint = 0,
                       limits = c(-1, 1), name = "Slope\n(scaled within\nvariable,\ncapped)") +
  scale_y_discrete(limits = rev) +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 11) +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        plot.background = element_rect(fill = "white", colour = NA))
ggsave(p_heat, file = paste0(PATHOUT, "/heatmap_slopes_thetaF1.png"),
       width = 7, height = 6, dpi = 300, device = ragg::agg_png, bg = "white")
}

#### PLOT AGE DIMENSION ------------------

## loading the dataset of traits averaged over cohorts Time1 to Time2, and reorganised to remove age=0, add reassign natural mortality as add mortality, and estimate natural mortality
data_variables <- load_trait_df(base_path = base_path, pathout=pathoutputRDATA, FORCE)

CLEAN = F
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
    data_t <- what_ages_have_null_dead_nb(data_t)
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
traits_to_plot <- list(c("Dead", "Deltasquare", "Mtot"), 
                       c("Dead", "Dead\nDivNumerator_N_alpha", "kx_bar"), 
                        c("b_rate_sc")
)

# select data_parameters
data_plot_1b <- as_long_if_needed(data_parameters, columns_keep_pivot = 
                                   colnames(data_parameters)[-which(colnames(data_parameters) %in% 
                                                                      c("replicate", "Fmax", "Fstart", "Species", "age"))])
tab_rep <- table(data_plot_1b$Fmax, data_plot_1b$age, data_plot_1b$Species)
tab_rep_all <-NULL
if (!CLEAN){
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
    data_plot_1b[which((data_plot_1b$Species == s & 
                          data_plot_1b$Fmax %in% unique(tab_remove$Fmax)) & 
                         data_plot_1b$age %in% unique(tab_remove$age)), ]
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
      data_sp_fmax <- abs(as.matrix(tempo_abund_mean[which(tempo_abund_mean$Species == sp), which(colnames(tempo_abund_mean) == fmax)]))
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
}

df_trait_t <- estimateX(trait_data_full_i = data_plot_1b, var = paste0("b_rate", "x", "survivorship_curve"))
data_plot_1b <- rbind(data_plot_1b, df_trait_t)

df_trait_t <- estimateDiv(trait_data_full_i = data_plot_1b, var = paste0("Dead", "\nDiv", "Numerator_N_alpha"))
data_plot_1b <- rbind(data_plot_1b, df_trait_t)
# clean dataset
data_plot_1b <- data_plot_1b %>% filter(variable %in% c(unlist(traits_to_plot), "SSDgx", "survivorship_curve", "Deltasquare", "Dead", "Mtot"))

# condition to stop
if (nrow(data_plot_1b)<1) next

#loop to plot relative and absolute data
slopes_list <- list()
for (relative_to_ref in c(T, "proportional")){
  
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
    
    # clean dataset
    if (sum(is.nan(df_trait_t$value))>0) df_trait_t <- df_trait_t[-which(is.nan(df_trait_t$value)), ]
    if (sum(is.infinite(df_trait_t$value))>0) df_trait_t <- df_trait_t[-which(is.infinite(df_trait_t$value)), ]
    
    # conditions to stop
    if (nrow(df_trait_t)==0) next
    if (all(is.na(df_trait_t$value))) {
      Ne_plot_list[[paste0(as.character(relative_to_ref), as.character(proportional))]][[sp]][[traits_selected]] <- NULL
      next}
    if (all(is.na(df_trait_t$value))) next
    
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


for (relative_to_ref in c(TRUE, "proportional")){
  
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
  
  if (!CLEAN){
    windowsFonts(myCustomWindowsFontName=windowsFont("DejaVu Sans"))
    cairo_pdf(file = paste0(PATHOUT, "/Clean_for_paper/TOTAge_plot", plot_name_ch, ".pdf"),
              width = 10, height =15, onefile=T)
    print(Age_plot_list[[plot_name]]$AtlanticHerringAtlanticMackerelSandeelEuropeanSpratNorwayPoutEuropeanPlaiceCommonSoleSaitheAtlanticCodHaddockHorseMackerelWhitingCommonDabGreyGurnardHake$DeadDeltasquareMtot)
    dev.off() 
    
    windowsFonts(myCustomWindowsFontName=windowsFont("DejaVu Sans"))
    cairo_pdf(file = paste0(PATHOUT, "/Clean_for_paper/Age_plot", plot_name_ch, ".pdf"),
              width = 10, height =12, onefile=T)
    print(Age_plot_list[[plot_name]]$`WhitingEuropean plaiceAtlantic mackerelHorse mackerel`$DeadDeltasquareMtot)
    dev.off()   
    
    windowsFonts(myCustomWindowsFontName=windowsFont("DejaVu Sans"))
    cairo_pdf(file = paste0(PATHOUT, "/Clean_for_paper/TOTAge2_plot", plot_name_ch, ".pdf"),
              width = 10, height =12, onefile=T)
    print(Age_plot_list[[plot_name]]$AtlanticHerringAtlanticMackerelSandeelEuropeanSpratNorwayPoutEuropeanPlaiceCommonSoleSaitheAtlanticCodHaddockHorseMackerelWhitingCommonDabGreyGurnardHake$`DeadDead\nDivNumerator_N_alphakx_bar`)
    dev.off()   
    
    windowsFonts(myCustomWindowsFontName=windowsFont("DejaVu Sans"))
    cairo_pdf(file = paste0(PATHOUT, "/Clean_for_paper/Age2_plot", plot_name_ch, ".pdf"),
              width = 10, height =12, onefile=T)
    print(Age_plot_list[[plot_name]]$`WhitingEuropean plaiceAtlantic mackerelHorse mackerel`$`DeadDead\nDivNumerator_N_alphakx_bar`)
    dev.off()   
  }
  
  
  if (!CLEAN){
    windowsFonts(myCustomWindowsFontName=windowsFont("DejaVu Sans"))
    cairo_pdf(file = paste0(PATHOUT, "/Clean_for_paper/TOTbx_sc", plot_name_ch, ".pdf"),
              width = 10, height =12, onefile=T)
    print(Age_plot_list[[plot_name]]$AtlanticHerringAtlanticMackerelSandeelEuropeanSpratNorwayPoutEuropeanPlaiceCommonSoleSaitheAtlanticCodHaddockHorseMackerelWhitingCommonDabGreyGurnardHake$`b_rate_sc`)
    dev.off()   
  }
}
