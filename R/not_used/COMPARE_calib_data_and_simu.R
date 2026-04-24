

main="New_gpu_conf.R"
source(file.path("osmose-ns-v4.4.1_clim/R/functions/source.R"))
repo        = "Simu_LONG"
multip_name <- c("m")
multipliers = c("0.5", "1", "1.5")
multipliers = 1

# prepare environment for this particular plot
rep=1
incr_multiplier <- 0
cat("load simulation output")
simulation = paste0(repo, "/output_", multip_name, multipliers)
Ne_version = paste0("plot_biomass", multip_name, multipliers)

ns1 = NULL
outputDir_Simu  = file.path(configDir, "output", simulation)
outputDir_Ne =paste0(configDir, "/output/", "Ne","/PLOT_REPLICATES/")
ns1  = read_osmose(path = outputDir_Simu)
conf = read_osmose(path = paste0(outputDir_Simu, "/ns-configuration.osm"))
nrep = conf$simulation.nsimulation
# multipliers = "output_etaequalsone_bioen_neutral"
trait_dim <- c("biomass", "landings", "meanSize")

fi              <-list.files("../osmose-ns_calibration_data/OSMOSE_NS_calibration/data/", full.names = T)
numbers         <- gregexpr("[0-9]+", fi)
result          <- regmatches(fi, numbers) 
numeric_result  <- as.numeric(unlist(result))
max_age         <- numeric_result[which.max(numeric_result)]

list_comparison <- list()
for (rep in 1:nrep){
  cat("\nREP ", rep, "\n")
  mean_fr <- as.data.frame(matrix(-99, ncol=length(species), nrow=c(max_age+length(trait_dim)-1)))
  colnames(mean_fr) <- species
  rownames(mean_fr) <- c(trait_dim[1], "yield", paste0(trait_dim[3], "ByAge", 0:c(max_age-1)))  
  sd_fr <- as.data.frame(matrix(-99, ncol=length(species), nrow=c(max_age+length(trait_dim)-1)))
  colnames(sd_fr) <- species
  rownames(sd_fr) <- c(trait_dim[1], "yield", paste0(trait_dim[3], "ByAge", 0:c(max_age-1)))
  for (sp in species){
    for (trait in trait_dim){
      cat("\n\t--->\t", sp, "\t", trait, '\n')
      f     <- fi[grep(sp, fi)][grep(trait, fi[grep(sp, fi)])]
      if (length(f)==0){next}
      if (trait == "landings") trait = "yield"
      
      # extract the corresponding data from the simulation output
      if (length(f)>1){
        trait <- paste0(trait, "ByAge")
        data  <- ns1[[trait]][[sp]][,,rep]
        for (fa in f){
          # calibration data
          fr    <- read.csv(fa)
          
          # extract that age data from simulation
          age     <- as.numeric(paste(str_extract_all(string = fa, pattern = '\\(\\d{2}\\)|\\d')[[1]], collapse=""))-1
          trait_i <- paste0(trait, age)
          data_a  <- data[, which(colnames(data) == age)]
          
          # compare simulation and calibration data
          data_trunc <- as.data.frame(data_a)
          data_trunc$time <- as.numeric(rownames(data_trunc))
          if (any(colnames(data_trunc) == "Time")) {
            data_trunc$time <- data_trunc$Time
            data_trunc      <- data_trunc[,-which(colnames(data_trunc) == c("Time"))]}
          # keep only the data corresponding to the calibration data
          data_trunc      <- data_trunc[which(data_trunc$time %in% time), ]
          
          # Compare calibration goal and simulated
          mean_fr[trait_i, sp] <- mean(data_trunc[, -which(colnames(data_trunc) =="time")] - fr[, -which(colnames(fr) =="Time")], na.rm = T)
          sd_fr[trait_i, sp]   <- sd(data_trunc[, -which(colnames(data_trunc) =="time")] - fr[, -which(colnames(fr) =="Time")], na.rm = T)
        }
      }else {
        fr    <- read.csv(f)
        time  <- fr$Time
        data <- ns1[[trait]][,sp, rep]
        if(is.null(data)){
          data     <- list()
          files    <- list.files(outputDir_Simu, full.names = T, recursive = T)
          files    <- files[c(grep(".csv", files))]
          files_rep0    <- files[c(grep("Simu0", files))]
          ID_files <- grep(trait, files_rep0)
          dim_data_ID_files <- unlist(sapply(X = files_rep0[ID_files], 
                                             FUN=function(x) {dim(read.csv(x, skip = c(1)))}))
          ID_file_selected  <- which.min(dim_data_ID_files[1,]) # chose the file the smallest : exclude the data by age or by size ...
          if (nrep>1){
            str_a <- stringr::str_remove(pattern = "0.csv", string = files_rep0[ID_files[ID_file_selected]])
            str_b <- stringr::str_remove(pattern = "0.nc", string = files_rep0[ID_files[ID_file_selected]])
            ID_file_selected_Simu <- grep(c(str_a, str_b)[which.min(c(nchar(str_a), nchar(str_b)))], files)
          } else {ID_file_selected_Simu <- which(files == files_rep0[ID_files[ID_file_selected]])}
          for (ID_i in 1:length(ID_file_selected_Simu)){
            if (length(ID_i)>1){cat("error dim selected outputs")
              stop()}
            ID <- ID_file_selected_Simu[ID_i]
            if(grepl("ByAge", files[ID]) || grepl("BySize", files[ID]) ){
              data_A <- as.data.frame(read.csv(files[ID], skip=c(1)))
              data_A <- data_A[, -c(which(colnames(data_A)=="Age"), which(colnames(data_A)=="Size"))]
              avg_by_species_age <- data_A %>%
                group_by(Time) %>%
                summarise(across(where(is.numeric), \(x) mean(x, na.rm = TRUE)), .groups = "drop")
              data <- avg_by_species_age[,sp]
              if (grepl("yieldN", files[ID])){
                data[1:c(nrow(data)-1), sp] <- data[1:c(nrow(data)-1), sp] - data[2:nrow(data), sp]
              }
            }else {data <- as.data.frame(read.csv(files[ID], skip=c(1)))}
          }
        }
        data_trunc <- as.data.frame(data)
        data_trunc$time <- as.numeric(rownames(data_trunc))
        if (any(colnames(data_trunc) == "Time")) {
          data_trunc$time <- data_trunc$Time
          data_trunc      <- data_trunc[,-which(colnames(data_trunc) == c("Time"))]}
        # keep only the data corresponding to the calibration data
        data_trunc         <- data_trunc[which(data_trunc$time %in% time), ]
        mean_fr[trait, sp] <- mean(data_trunc[,1] - fr[, -which(colnames(fr) =="Time")], na.rm = T)
        sd_fr[trait, sp]   <- sd(data_trunc[,1] - fr[, -which(colnames(fr) =="Time")], na.rm = T)
      }
    }
  } 
  list_comparison[[rep]] <- list(MEAN = mean_fr, SD = sd_fr) 
}


save(file="COMPARE_RESTART_REP_WITH_CALIB.RData", list_comparison)








######### PLOT ##################

# --- Load the dataset ---
sim_list <- get(load(c("osmose-ns-v4.4.1_clim/output/Ne/PLOT_REPLICATES/COMPARE_RESTART_REP_WITH_CALIB.RData")))

library(tidyverse)        # dplyr, tidyr, ggplot2, purrr
library(tibble)           # rownames_to_column

# --- convert one replicate (handles traits in rownames, species in column names) ---
replicate_to_long <- function(rep_item, rep_id) {
  mean_df <- rep_item$MEAN %>% 
    rownames_to_column(var = "trait") %>%      # move rownames -> trait column
    pivot_longer(cols = -trait, names_to = "species", values_to = "mean")
  
  sd_df <- rep_item$SD %>% 
    rownames_to_column(var = "trait") %>%
    pivot_longer(cols = -trait, names_to = "species", values_to = "sd")
  
  dat <- left_join(mean_df, sd_df, by = c("trait", "species"))
  dat <- dat %>% mutate(replicate = rep_id,
                        mean = ifelse(mean==-99, 0, abs(as.numeric(mean))),
                        sd   = ifelse(mean==-99, 0, as.numeric(mean)),
                        ymin = ifelse(mean == -99, 0, abs(mean) - sd*mean),
                        ymax = ifelse(mean == -99, 0, abs(mean) + sd*mean))
  dat
}

# --- build full tidy table from the list ---
df_all <- map2_dfr(sim_list, seq_along(sim_list), replicate_to_long)
df_all[which(df_all$species == "Shrimp" & df_all$trait == "biomass"), c("ymin")] <- df_all[which(df_all$species == "Shrimp" & df_all$trait == "biomass"), "ymin"]/500
df_all[which(df_all$species == "Shrimp" & df_all$trait == "biomass"), c("ymax")] <- df_all[which(df_all$species == "Shrimp" & df_all$trait == "biomass"), "ymax"]/500
# quick check
glimpse(df_all)
# columns: trait, species, mean, sd, replicate, ymin, ymax

# ============================
# Plot option A: one panel per replicate (facets = replicates),
# x-axis = trait, y = mean ± sd (species shown via color)
# ============================
pA <- ggplot(df_all, aes(x = trait, y = mean, ymin = ymin, ymax = ymax, color = species)) +
  geom_pointrange(position = position_dodge(width = 0.7), size = 0.4) +
  facet_wrap(~ replicate, ncol = 4) +
  coord_flip() +
  # scale_y_log10() +  # ← log10 scale for y-axis
  theme_minimal() +
  theme(axis.text.y = element_text(size = 7),
        legend.position = "right") +
  labs(x = "Trait", y = "Mean ± SD", title = "Per-replicate: trait means ± SD (colored by species)")

pdf(file = paste0("COMPARE_rel_mean_by_rep(all).pdf"))
pA
dev.off()


# ============================
# Plot option B: one panel per species (facets = species),
# x-axis = trait, inside each panel show the 20 replicate intervals (x grouped by trait, color by replicate)
# ============================
# to show replicates distinctively, make replicate a factor
df_all <- df_all %>% mutate(replicate_f = factor(replicate))
# because we are using the log : avoid the zeros
offset <- min(df_all$ymin[df_all$ymin > 0]) / 10  # small fraction of smallest positive
df_all$mean_offset <- df_all$mean + offset
df_all$ymin_offset <- df_all$ymin + offset
df_all$ymax_offset <- df_all$ymax + offset


nrep = length(sim_list)
list_pB <- list()
df_all_pB <- df_all[-which(df_all$species == "Shrimp" & df_all$trait == "biomass"),]
df_all_pB <- df_all_pB[-which(df_all_pB$species == "GreyGurnard" & df_all_pB$trait %in% c("yield", "biomass")),]
y_min <- min(df_all_pB$ymin, na.rm = TRUE)
y_max <- max(df_all_pB$ymax, na.rm = TRUE)

df_all[which(df_all$species == "Shrimp" & df_all$trait == "biomass"),]
df_all[which(df_all$species == "GreyGurnard" & df_all$trait %in% c("yield", "biomass")),]

pA <- ggplot(df_all_pB, aes(x = trait, y = mean, ymin = ymin, ymax = ymax, color = species)) +
  geom_pointrange(position = position_dodge(width = 0.7), size = 0.4) +
  facet_wrap(~ replicate, ncol = 4) +
  coord_flip() +
  # scale_y_log10() +  # ← log10 scale for y-axis
  theme_minimal() +
  theme(axis.text.y = element_text(size = 7),
        legend.position = "right") +
  labs(x = "Trait", y = "Mean ± SD", title = "Per-replicate: trait means ± SD (colored by species)")


for (rep in 1:c(nrep)){
  pB <- ggplot(df_all_pB[which(df_all_pB$replicate == rep), ], aes(x = trait, y = mean, ymin = ymin, ymax = ymax)) +
    geom_pointrange(position = position_dodge(width = 0.8), size = 0.4) +
    # scale_y_log10(limits=c(floor(y_min), cieling(y_max))) +  # ← log10 scale for y-axis
    ylim(c(y_min, y_max))+
    facet_wrap(~ species, ncol = 4, scales = "fixed")+
    ylab(c(floor(y_min), ceiling(y_max)))+
    coord_flip() +
    theme_minimal() +
    theme(axis.text.y = element_text(size = 7),
          legend.position = "bottom") +
    labs(x = "Trait", y = "Mean ± SD", color = "Replicate",
         title = paste0("Replicate  ", rep,", per-species: relative trait means ± SD"))
  list_pB[[rep]] <- pB
}

pdf(file = paste0("COMPARE_rel_mean_by_rep(no_biomassofGG&Shr_noyieldGG).pdf"))
for (l in list_pB){
  print(l)
}
dev.off()

pdf(file = paste0("COMPARE_rel_mean_by_sp(no_biomassofGG&Shr_noyieldGG).pdf"))
pA
dev.off()

# Print the plot you prefer:
print(pA)   # option A
# or
print(pB)   # option B




# load(c("osmose-ns-v4.4.1_clim/output/Ne/PLOT_REPLICATES/COMPARE_RESTART_REP_WITH_CALIB.RData"))
# length(list_comparison)
# 
# library(dplyr)
# library(tidyr)
# library(purrr)
# library(ggplot2)
# 
# sim_list <- list_comparison
# 
# 
# # PLOT ALL THE ERRORS (mean and SD) = differences between simulated and calibration data
# df_all <- map2_dfr(sim_list, seq_along(sim_list), function(x, i) {
#   mean_df <- x$MEAN %>%
#     mutate(var_type = "mean", replicate = i)
#   sd_df <- x$SD %>%
#     mutate(var_type = "sd", replicate = i)
#   bind_rows(mean_df, sd_df)
# })
# df_long <- df_all %>%
#   pivot_longer(
#     cols = -c(var_type, replicate),
#     names_to = "species",
#     values_to = "value"
#   )
# 
# ggplot(df_long, aes(x = species, y = value, fill = var_type)) +
#   geom_boxplot(position = position_dodge(width = 0.8)) +
#   facet_wrap(~ replicate, ncol = 4) +
#   theme_bw() +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
#   labs(title = "Trait values (mean & SD) across 20 replicates",
#        y = "Value", x = "Species")
# 
# 
# ggplot(df_long, aes(x = species, y = value, fill = var_type)) +
#   geom_point(position = position_dodge(width = 0.8), pch=21) +
#   facet_wrap(~ replicate, ncol = 4) +
#   scale_y_log10() +  # ← log10 scale for y-axis
#   theme_bw() +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
#   labs(title = "Trait values (mean & SD) across 20 replicates",
#        y = "Value (log scale)", x = "Species")
# 
# 
# 
# 
# 
# 
# 
# # PLOT THE ERRORS > 20% (mean and SD) = differences between simulated and calibration data
# df_all <- 
#   map2_dfr( # returns data frames by binding all the replicate + mean and sd errors together (880 = 22(traits)*2(mean + sd)*20(replicates))
#   sim_list,
#   .y = seq_along(sim_list),
#   ~ bind_rows(
#     mutate(.x$MEAN, replicate = .y, var_type = "mean"),
#     mutate(.x$SD, replicate = .y, var_type = "sd")
#   ),
#   .id = "replicate"
# ) %>%
#   pivot_longer( # from a column per species, to one column with all the species
#     cols = -c(replicate, var_type),
#     names_to = "species",
#     values_to = "value"
#   )
# 
# df_diff_filtered <- df_all %>% # select only the errors > 10%
#   filter(var_type == "sd", value > 0.1)
# 
# df_plot <- df_all %>% # remove the errors >10%
#   semi_join(df_diff_filtered, by = c("species", "var_type", "value"))
# 
# ggplot(df_plot, aes(x = species, y = value, fill = var_type)) +
#   geom_boxplot(position = position_dodge(width = 0.8)) +
#   facet_wrap(~ replicate, ncol = 4) +
#   theme_bw() +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
#   labs(title = "Trait values (mean & SD) across 20 replicates",
#        y = "Value", x = "Species")
# 
# 
# ggplot(df_plot, aes(x = species, y = value, fill = var_type)) +
#   geom_boxplot(position = position_dodge(width = 0.8)) +
#   facet_wrap(~ replicate, ncol = 4) +
#   scale_y_log10() +  # ← log10 scale for y-axis
#   theme_bw() +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
#   labs(title = "Trait values (mean & SD) across 20 replicates",
#        y = "Value (log scale)", x = "Species")








