
# the idea is to store separately all the key parameters per species into plots
# and to produce one pdf sheet per species with all the key variables to have it all summarized at one glance




#### Initialisation --------------
library(legendry)
library(ggpubr)
library(ggrepel)
library(ggforce)
library(ggnewscale)
library(patchwork)

# functions and environment
main = "New_gpu_conf.R" # config file
source("R/functions/source.R") # load config and key parameters and functions
source("R/functions/function_data_traits_summary.R") # datasets preparation & loading
source("R/functions/functions_plots.R") # plots 

# paths
# path <- paste0(simu_repo, "/output/PLOTS/FINAL_Ne_Htz_mature_short/")
# PATHOUT  = "output/PLOTS/CHAP2/KEY_TRAITS_MATURE_removed0_harmo_T_Vk_Nalpha/"
PATHOUT  = "output/PLOTS/CHAP2/KEY_TRAITS_MATURE_removed0inNedata/"#KEY_TRAITS_MATURE_150-150/"
if (!dir.exists(PATHOUT)) dir.create(PATHOUT, recursive = TRUE)
if (!dir.exists(paste0(PATHOUT, "/Clean_for_paper"))) dir.create(paste0(PATHOUT, "/Clean_for_paper"), recursive = TRUE)
pathoutputRDATA = paste0(PATHOUT, "RDATA")
if(!dir.exists(pathoutputRDATA)) dir.create(pathoutputRDATA, recursive = T)
base_path = paste0("output/DATASETS/Chap2/Trait_variables_by_cohort_by_age2")
base_path = paste0("output/DATASETS/Chap2/Trait_variables_by_cohort_by_age_withWbyAge/")
# folders_long <- list.dirs(paste0(path_long, "INIT"), recursive = F)

folders_traits_c <- paste0("output/DATASETS/Chap2/Trait_variables_by_cohort/RData/")
folders_traits_c  <- list.dirs(folders_traits_c, recursive = F)
folders_traits_c_a <- paste0("output/DATASETS/Chap2/Trait_variables_by_cohort_by_age/RData/")
folders_traits_c_a <- paste0("output/DATASETS/Chap2/Trait_variables_by_cohort_by_age_withWbyAge/RData/")
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

## % activate if needed : will produce .RData files in the folder paste0(simu_repo, "/output/PLOTS/FINAL_Ne_Htz_mature/")
## % 1 file per species, with the dataset at the scale of the cohorts, starting at age = alpha (first mature age)
dataset_Ne_mean_cohorts(species, FORCE, pathdata=paste0("output/DATASETS/Chap2/FINAL_Ne_Htz_mature/"),
                        pathout=pathoutputRDATA, Time2, Time2)

## loading the dataset of traits averaged over cohorts Time1 to Time2, and reorganised to remove age=0, add reassign natural mortality as add mortality, and estimate natural mortality
data_variables <- load_trait_df(base_path = base_path, pathout=pathoutputRDATA, FORCE)


#### DATASETS clean ------------------------
# CLEAN DATASET OF TRAITS VARIABLES
data_variables <- when_variables_missing_age_dim(data_variables)
data_variables <- when_null_or_low_abund(dataset=data_variables)
data_variables <- data_variables[which(data_variables$Fstart == 1),]

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
  data_t <- what_ages_have_null_dead_nb(data_t)
  data_t <- when_null_or_low_abund(data_t)
  data_t$Fmax <- as.factor(data_t$Fmax) 
  data_t$min_age <- min_age[[sp]]
  
  ## dataset
  data <-  data_t %>% filter(Species == sp, Fstart == 1)
  if (is.null(data_parameters)){ data_parameters <- data 
  } else data_parameters <- rbind(data_parameters, data)
}


# df_lxbx : one row per Species x age x Fmax x variable ("lx", "bx"), value = mean over replicates
# [option] add lx*bx (age-specific reproductive output, the quantity that weights T)
df_lxbx <- data_parameters[which(data_parameters$variable %in% c("survival_rate", "b_rate_sc"#, "survivorship_curve"
                                                                 )),]
df_lxbx$variable <- improve_vriable_names(df_lxbx$variable, age = T)

df_plot <- df_lxbx %>%
  tidyr::pivot_wider(names_from = variable, values_from = value) %>%
  tidyr::pivot_longer(c(sx, #lx, 
                        bx), names_to = "variable", values_to = "value") %>%
  group_by(Species, age, Fmax, variable) %>%
  summarise(n     = sum(!is.na(value)),
            sd    = sd(value, na.rm = TRUE),
            value = mean(value, na.rm = TRUE),
            .groups = "drop") %>%
  mutate(lo = value - sd, hi = value + sd,
         lo = ifelse(variable == "lx", pmax(lo, 0), lo)) %>%  # lx cannot be < 0
  mutate(variable = factor(variable, levels = c("sx", #"lx", 
                                                "bx"),
                           labels = c("Survival~rate~(italic(s)[x])", 
                                      # "Survivorship~(italic(l)[x])",
                                      "Fecundity~(italic(b)[x])")),
         Species = gsub("([a-z])([A-Z])", "\\1 \\L\\2", Species, perl = TRUE),
         Fmax = factor(Fmax))

# diverging palette centred on the reference thetaF = 1 (neutral grey)
col_F <- c("0.5" = "#2166AC", "0.75" = "#67A9CF", "1" = "#4D4D4D",
           "1.25" = "#EF8A62", "1.5" = "#B2182B")

y_cap   <- 40
var_cap <- "Fecundity~(italic(b)[x])"   # exact label of the row to cap (see levels(df_plot$variable))
p_sp <- function(sp, first) {
  cat(sp)
  d <- df_plot[df_plot$Species == sp, ]
  rows_cap <- d$variable == var_cap
  if (max(d$hi[rows_cap], na.rm = TRUE) > y_cap) {
    y_cap_i <- min(ifelse(is.infinite(max(d$value[rows_cap], na.rm=T)), 
                        Inf, max(d$value[rows_cap]))+1, y_cap)
    cat(y_cap_i, "\n")
    cat(d$value[rows_cap])
    d$hi[rows_cap] <- pmin(d$hi[rows_cap], y_cap_i)        # ribbon stops at the cap
    d$lo[rows_cap] <- pmax(d$lo[rows_cap], -1)
    d$value[rows_cap & d$value > y_cap_i] <- NA             # lines/points above the cap are not drawn
  }

  ggplot(d, aes(x = age, y = value, colour = Fmax, group = Fmax)) +
    # [AJOUT] between-replicate SD as a light ribbon
    geom_ribbon(aes(ymin = lo, ymax = hi, fill = Fmax), alpha = 0.1, colour = NA) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 1.3) +
    facet_grid(variable ~ Species, scales = "free_y", switch = "y",   # y free within EACH species
               labeller = labeller(variable = label_parsed)) +
    scale_colour_manual(values = col_F, name = expression(theta[F])) +
    scale_fill_manual(values = col_F, name = expression(theta[F])) +   # [AJOUT] same legend as the lines
    scale_y_continuous(labels = scales::label_number(scale_cut = scales::cut_short_scale()),
                       expand = expansion(mult = c(0.02, 0.06))) +
    scale_x_continuous(breaks = scales::breaks_pretty(n = 5)) +
    labs(x = "Age (years)", y = NULL) +
    theme_minimal(base_size = 11) +
    theme(panel.border = element_rect(colour = "grey30", fill = NA, linewidth = 0.5),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_line(colour = "grey92", linewidth = 0.3),
          strip.placement = "outside",
          strip.text.y.left = element_text(angle = 90, size = 10),
          strip.text.x = element_text(face = "bold", size = 11),
          axis.text = element_text(colour = "grey30"),
          legend.position = "top",
          plot.background = element_rect(fill = "white", colour = NA)) +
    guides(colour = guide_legend(nrow = 1, override.aes = list(linewidth = 1.2))) +
    # row labels only on the first column
    (if (!first) theme(strip.text.y.left = element_blank()) else theme())
}

# one column per species, each with its own y scales (facet_grid alone shares y along a row)
df_plot <- na.omit(df_plot)
if (any(df_plot$n)<3) df_plot <- df_plot[-which(df_plot$n<=5),]
sp_list <- unique(df_plot$Species)
sp_list <- sp_list[which(sp_list %in% c("Atlantic mackerel", "European plaice"))]
p_lxbx  <- wrap_plots(lapply(seq_along(sp_list), function(i) p_sp(sp = sp_list[i], i == 1)), nrow = 1) +
  plot_layout(guides = "collect", axis_titles = "collect") &
  theme(legend.position = "top")

ggsave("output/for_manuuscript_lx_bx.png", p_lxbx, width = 7.5, height = 7, dpi = 300, device = ragg::agg_png, bg = "white")

# 
# 
# p_sp <- function(sp, first) ggplot(df_plot[df_plot$Species == sp, ], aes(x = age, y = value, colour = Fmax, group = Fmax)) +
#   geom_line(linewidth = 0.8) +
#   geom_point(size = 1.3) +
#   facet_grid(variable ~ Species, scales = "free_y", switch = "y",   # y free within EACH species
#              labeller = labeller(variable = label_parsed)) +
#   scale_colour_manual(values = col_F, name = expression(theta[F])) +
#   scale_y_continuous(labels = scales::label_number(scale_cut = scales::cut_short_scale()),
#                      expand = expansion(mult = c(0.02, 0.06))) +
#   scale_x_continuous(breaks = scales::breaks_pretty(n = 5)) +
#   labs(x = "Age (years)", y = NULL) +
#   theme_minimal(base_size = 11) +
#   theme(panel.border = element_rect(colour = "grey30", fill = NA, linewidth = 0.5),
#         panel.grid.minor = element_blank(),
#         panel.grid.major = element_line(colour = "grey92", linewidth = 0.3),
#         strip.placement = "outside",
#         strip.text.y.left = element_text(angle = 90, size = 10),
#         strip.text.x = element_text(face = "bold", size = 11),
#         axis.text = element_text(colour = "grey30"),
#         legend.position = "top",
#         plot.background = element_rect(fill = "white", colour = NA)) +
#   guides(colour = guide_legend(nrow = 1, override.aes = list(linewidth = 1.2))) +
#   # row labels only on the first column
#   (if (!first) theme(strip.text.y.left = element_blank()) else theme())
# 
# # one column per species, each with its own y scales (facet_grid alone shares y along a row)
# sp_list <- unique(df_plot$Species)
# p_lxbx  <- wrap_plots(lapply(seq_along(sp_list), function(i) p_sp(sp_list[i], i == 1)), nrow = 1) +
#   plot_layout(guides = "collect", axis_titles = "collect") &
#   theme(legend.position = "top")
# 
# ggsave("lx_bx.png", p_lxbx, width = 7.5, height = 7, dpi = 300, device = ragg::agg_png, bg = "white")