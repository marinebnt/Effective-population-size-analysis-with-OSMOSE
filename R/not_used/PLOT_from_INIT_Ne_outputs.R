library(ggplot2)
library(gridExtra)   # pour arranger plusieurs plots
library(patchwork)   # alternative moderne (recommandée)
library(stringr)
library(dplyr)
library(GGally)
library(tidyverse)


base_dir <- "osmose-ns-v4.4.1_clim/output/PLOTS/All_genetic_drift_indicators/INIT/"
folders  <- list.dirs(base_dir, recursive = F)



trait = "longevity"


for (fol in folders[1]){
  df_fol = get(load(list.files(fol)))
  rapply(df_fol, function(df_fol) {rbind()})
  rapply(df_fol, c(1,2), sd)
}