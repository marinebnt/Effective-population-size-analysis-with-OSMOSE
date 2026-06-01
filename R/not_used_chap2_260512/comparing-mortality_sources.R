# the idea is to store separately all the key parameters per species into plots
# and to produce one pdf sheet per species with all the key variables to have it all summarized at one glance




#### Initialisation --------------
library(legendry)
library(ggpubr)

# functions and environment 
main = "New_gpu_conf.R" # config file
source("osmose-ns-v4.4.1_clim/R/functions/source.R") # load config and key parameters and functions
source("osmose-ns-v4.4.1_clim/R/functions/functions_key_traits.R") # plots against Fmax
source("osmose-ns-v4.4.1_clim/R/functions/auxiliary_for_lists.R") # from list to df 
source("osmose-ns-v4.4.1_clim/R/functions/functions_output_SimuRData_into_df_for_plots.R") # to transform dataset from RData to long df
source("osmose-ns-v4.4.1_clim/R/functions/function_data_traits_summary.R")

# paths
PATHOUT  = "osmose-ns-v4.4.1_clim/output/PLOTS/KEY_TRAITS/"
if (!dir.exists(PATHOUT)) dir.create(PATHOUT, recursive = TRUE)
base_dir <- "osmose-ns-v4.4.1_clim/output/PLOTS/FINAL_Ne_Htz2/INIT"
folders  <- list.dirs(base_dir, recursive = F) # files INIT2 with raw Ne datasets

folders_traits_y <- "osmose-ns-v4.4.1_clim/output/PLOTS/Trait_variables_by_year_BYSP//RData/"
folders_traits_y <- list.dirs(folders_traits_y, recursive = F)
folders_traits_c <- "osmose-ns-v4.4.1_clim/output/PLOTS/Trait_variables_by_cohort/RData/"
folders_traits_c  <- list.dirs(folders_traits_c, recursive = F)
folders_traits_y_a <- "osmose-ns-v4.4.1_clim/output/PLOTS/Trait_variables_by_year_by_age2/RData/"
folders_traits_y_a <- list.dirs(folders_traits_y_a, recursive = F)
folders_traits_c_a <- "osmose-ns-v4.4.1_clim/output/PLOTS/Trait_variables_by_cohort_by_age2/RData/"
folders_traits_c_a  <- list.dirs(folders_traits_c_a, recursive = F)

# load or format dataset
pathoutputRDATA = paste0(configDir, "/output/PLOTS/KEY_TRAITS/RDATA")
if(!dir.exists(pathoutputRDATA)) dir.create(pathoutputRDATA, recursive = T)
Time1=150
Time2=160
species_df=species





pathout <- paste0("osmose-ns-v4.4.1_clim/output/PLOTS/Trait_variables_by_cohort_by_age2/otherplots")
if (!dir.exists(pathout)) dir.create(pathout, recursive = T)

# load dataset
data_complete_sp = NULL
data_complete_c = NULL
base_path = "osmose-ns-v4.4.1_clim/output/PLOTS/Trait_variables_by_cohort_by_age2"
folder_name <- basename(base_path)
suffix <- gsub("Trait_variables_by_", "", folder_name)
suffix <- gsub("_", "_", suffix)  # keep underscores
save_dir <- file.path(dirname(dirname(base_path)), "RDATA")
save_file <- file.path(save_dir, paste0("summary_trait_", suffix, ".RData"))
if (!dir.exists(save_dir)){dir.create(save_dir, recursive = T)}
FORCE=F
if (!file.exists(save_file) || FORCE==T){
  data_complete_c <- data_traits_summary(base_path)
} else {
  data_complete_c = get(load(file=save_file))
}

# dimensions dataset required : [time, Sp, value, rep, Fish_press, TRAIT]
head(data_complete_c)
data_complete_c <- data_complete_c[, which(colnames(data_complete_c) %in% c("age", "time", "Species", "value", "TRAIT", "replicate", "Fish_press"))]
data_complete_c[grep(data_complete_c$TRAIT, pattern="fec"), "TRAIT"] <- "mean_of_eggNumber_per_age"
unique(data_complete_c$TRAIT)

# adding natural mortality dimension
data_complete_c[which((data_complete_c$TRAIT) == "mean_of_Mnat_per_age"), "TRAIT"] <- "mean_of_Madd_per_age"
data_complete_c <- pivot_wider(data = data_complete_c, names_from = "TRAIT", values_from = "value")
data_complete_c <- data_complete_c[,-grep(colnames(data_complete_c), pattern="mean_of_Dis_per_age")]
data_complete_c <- data_complete_c[,-grep(colnames(data_complete_c), pattern="mean_of_Z_per_age")]
data_complete_c$mean_of_Mnat_per_age <- data_complete_c$mean_of_Madd_per_age + data_complete_c$mean_of_Mpred_per_age + data_complete_c$mean_of_Mstar_per_age
# data_complete_c <- pivot_longer(data = data_complete_c, names_to = "TRAIT", values_to = "value", cols=-colnames(data_complete_c)[1:4])

# changing colnames to match expectations
colnames(data_complete_c)[which(colnames(data_complete_c) == "TRAIT")] <- "variable"
data_complete_c$Fmax   <- as.numeric(unlist(stringr::str_remove(unlist(stringr::str_split(data_complete_c$Fish_press, " Fstart"))[seq(1, 2*length(data_complete_c$Fish_press), by=2)], "Fmax=")))
data_complete_c$Fstart <- as.numeric(unlist(stringr::str_split(data_complete_c$Fish_press, " Fstart="))[seq(2, 2*length(data_complete_c$Fish_press), by=2)])
data_complete_c <- (data_complete_c[,-which(colnames(data_complete_c) == "Fish_press")])

colnames(data_complete_c)


Time1=150
Time2=160

PATHOUT  = "osmose-ns-v4.4.1_clim/output/PLOTS/KEY_TRAITS/"
path <- c("C:/Users/mbeneat/Documents/osmose/Chap2/osmose-ns/osmose-ns-v4.4.1_clim/output/PLOTS/FINAL_Ne_Htz2")


### produce dataset reading multiple .RData to make 1 dataframe, stored as a .RData file ----------
sp_i = "Hake"
FORCE=F








## PLOT DELTA²xDEAD_i/N and compare the mortality source of variation in the results -----------------

combined_plot =list()
Time1=150
Time2=160
traits_sp_to_plot = list(c("Deltasquare","Dead"),c("Deltasquare\nDivNumerator_N_1","Dead\nDivNumerator_N_1"),
                         c("DeltasquarexDead", "DeltasquarexDead\nDivNumerator_N_1"), 
                         c("survival_rate"), c("kx_point", "kx_bar"))
PATHOUT  = "osmose-ns-v4.4.1_clim/output/PLOTS/KEY_TRAITS/"
path <- c("C:/Users/mbeneat/Documents/osmose/Chap2/osmose-ns/osmose-ns-v4.4.1_clim/output/PLOTS/FINAL_Ne_Htz2")
plot_Deathsdelta_age = list()
plot_Deathsdelta_age2=list()
inset.tot.F.Nat = list()
for (sp in species){
  cat(sp)

  load(file=paste0(path, "/Complete_INIT2", sp, ".RData"))
  final_dataset_age = complete[[2]]
  final_dataset_age <- final_dataset_age %>%
    dplyr::mutate(dplyr::across(c(cohort, replicate), ~ stringr::str_remove(.x, "^(cohort_|rep_)")))
  final_dataset_age <- final_dataset_age[which(final_dataset_age$cohort %in% Time1:Time2),]
  rm(complete)

  # remove ages that are not included
  final_dataset_age$age <- as.numeric(final_dataset_age$age)
  ages_to_remove <- table(final_dataset_age[which(final_dataset_age$variable == "Dead" & final_dataset_age$value ==0),"age"])
  ages_to_remove <- ages_to_remove[which(ages_to_remove == max(ages_to_remove))]
  final_dataset_age <- final_dataset_age[-which(final_dataset_age$age %in% as.numeric(names(ages_to_remove))), ]

  df_trait_sp <- final_dataset_age %>% filter(Species == sp, Fstart == 1)

  # avoid the estimation bias of very populations of very low size
  if (any(((df_trait_sp[which(df_trait_sp$variable == "abundanceByAge"), "value"])<10))){
    df_trait_sp <- pivot_wider(df_trait_sp, names_from = "variable", values_from = "value")
    df_trait_sp <- df_trait_sp[-which((df_trait_sp$abundanceByAge)<10), ]
    df_trait_sp <- pivot_longer(df_trait_sp, cols = -c("age", "replicate", "Fmax", "Fstart", 'Species', "cohort"), names_to = "variable", values_to = "value")
  }

  # remove cohort dimension
  df_trait_sp <- df_trait_sp %>%
    mutate(age = as.numeric(age))%>%
    mutate(replicate = as.numeric(replicate)-1)%>% # to make the datasets match : unify the replicate range
    mutate(Fmax = as.numeric(Fmax))%>%
    mutate(Fstart = as.numeric(Fstart))
  weight_df <- build_weight_df(df_trait_sp, weight_var = "abundanceByAge", key_cols = c("age", "Fmax", "Fstart", "Species", "replicate"))
  df_trait_sp <- weighted_reduce(
    df = df_trait_sp,
    weight_df = weight_df,
    join_cols = c("Species", "replicate", "Fmax", "Fstart", "age"),
    group_cols = c("Species", "replicate", "Fstart", "Fmax", "variable", "age"),
    drop_age_zero = TRUE
  )

  # combine all traits: Mortality + Ne traits : both need to be pivot_wider to be combined
  df_trait_sp <- pivot_wider(df_trait_sp, names_from = "variable", values_from = "value")
  data_complete_c_i <- data_complete_c %>%
    filter(Species == sp)
  data_complete_c_i <- pivot_wider(data_complete_c_i, names_from = "variable", values_from = "value")
  df_trait_sp <- full_join(data_complete_c_i,df_trait_sp, by=c("age", "Fmax", "Fstart", "Species", "replicate"))
  if (any(df_trait_sp$age == 0)) df_trait_sp <- df_trait_sp[-which(df_trait_sp$age == 0), ]
  df_trait_sp <- pivot_longer(data = df_trait_sp, cols = -which(colnames(df_trait_sp) %in% c("age", "Fmax", "Fstart", "Species", "replicate")),
                              names_to = "variable", values_to = "value")
  df_trait_sp <- na.omit(df_trait_sp)

  # add deaths estimates per mortaltiy source
  df_trait_t    <- estimateSum(trait_data_full_i = df_trait_sp, var=paste0("mean_of_F_per_age", "\nsum", "mean_of_Mnat_per_age"))
  df_trait_t$variable    <- "sum_mortality_rates"
  df_trait_sp   <- rbind(df_trait_sp, df_trait_t)
  var_i <- unique(df_trait_sp$variable)
  for (v in var_i){
    if (any(grepl(pattern="mean_of_M", v) | grepl(pattern="mean_of_F", v) | grepl(pattern="sum_mortality_rates", v))) {
      df_trait_t    <- estimateMnumber(trait_data_full_i = df_trait_sp, var=paste0(v))
      df_trait_sp   <- rbind(df_trait_sp, df_trait_t)
    }
  }  
  for (v in var_i){
    if (any(grepl(pattern="mean_of_M", v) | grepl(pattern="mean_of_F", v) | grepl(pattern="sum_mortality_rates", v))) {
      df_trait_t    <- estimateMnumber(trait_data_full_i = df_trait_sp, var=paste0(v), var2="survivorship_curve")
      df_trait_sp   <- rbind(df_trait_sp, df_trait_t)
    }
  }
  df_trait_t    <- estimatecolSum(trait_data_full_i = df_trait_sp, var=paste0("colSumall_deaths"))
  df_trait_sp   <- rbind(df_trait_sp, df_trait_t)
  var_i <- unique(df_trait_sp$variable)
  for (v in var_i){
    if (any(grepl(pattern="death", v) & c(grepl(pattern="M", v)) | (grepl(pattern="F_death", v)) | (v=="all_deaths")| (v=="all_deaths_rate"))) {
      df_trait_t    <- estimateDiv(trait_data_full_i = df_trait_sp, var=paste0(v, "\nDiv", "colSumall_deaths"))
      df_trait_t2   <- estimatecolSum(trait_data_full_i = df_trait_t, var=paste0("colSum", (df_trait_t$variable)[1]))
      df_trait_sp   <- rbind(df_trait_sp, df_trait_t)
      df_trait_sp   <- rbind(df_trait_sp, df_trait_t2)
    }
  }

  #add delta square
  df_trait_t    <- estimateSquare(trait_data_full_i = df_trait_sp, var=paste0("Deltasquare"))
  df_trait_sp   <- rbind(df_trait_sp, df_trait_t)

  #add delta square
  for (v in var_i){
    if (any(grepl(pattern="deaths", v) | (v=="Dead"))) {
      df_trait_t    <- estimateX(trait_data_full_i = df_trait_sp, var=paste0("Deltasquarex", v))
      df_trait_sp   <- rbind(df_trait_sp, df_trait_t)
    }
  }

  #add sum by age of SSDgx and SSDix
  df_trait_t    <- estimatecolSum(trait_data_full_i = df_trait_sp, var=paste0("colSumSSDix"))
  df_trait_sp   <- rbind(df_trait_sp, df_trait_t)
  df_trait_t    <- estimatecolSum(trait_data_full_i = df_trait_sp, var=paste0("colSumSSDgx"))
  df_trait_sp   <- rbind(df_trait_sp, df_trait_t)
  df_trait_t    <- estimatecolSum(trait_data_full_i = df_trait_sp, var=paste0("colSumDeltasquare"))
  df_trait_sp   <- rbind(df_trait_sp, df_trait_t)

  var_i <- unique(df_trait_sp$variable)

  # keep abundance by age variable on the side
  WEIGHT = df_trait_sp[which(df_trait_sp$variable=="abundanceByAge"),]

  # compare all traits to their value /N1 to understand where is the increase coming from
  df_trait_t    <- estimatecolSum(trait_data_full_i = df_trait_sp, var=paste0("colSumDeltasquarexDead"))
  df_trait_sp   <- rbind(df_trait_sp, df_trait_t)

  for (v in var_i){
    if (v != "Numerator_N_1") {
      df_trait_t    <- estimateDiv(trait_data_full_i = df_trait_sp, var=paste0(v, "\nDivNumerator_N_1")) # div by population size
      df_trait_sp   <- rbind(df_trait_sp, df_trait_t)
    }
  }

  var_i <- unique(df_trait_sp$variable)

  for (relative_to_ref in c(FALSE)){
    proportional = FALSE
    if (relative_to_ref == "proportional") {
      relative_to_ref = FALSE
      proportional = TRUE
    }
    relative_to_ref = as.logical(relative_to_ref)
    plot_name_ch=c()
    plot_name <- paste0(as.character(relative_to_ref), as.character(proportional))
    if (plot_name=="TRUEFALSE") plot_name_ch <- c("Relative")
    if (plot_name=="FALSEFALSE") plot_name_ch <- c("Absolute")
    if (plot_name=="FALSETRUE") plot_name_ch <- c("Proportional")
    # for the plot
    plot_biom = ifelse(relative_to_ref, "rel_biomassplot", "biomassplot")
    plot_biom = ifelse(relative_to_ref, "biomassplot", "biomassplot")

    ref_data <- df_trait_sp %>%
      filter(Fstart == 1 & Fmax==1) %>%
      dplyr::select(variable, value, replicate, age) %>%
      dplyr::rename(ref_value = value)
    df_trait_sp_a <- df_trait_sp

    df_trait_sp_a <- clean_values(df_trait_sp_a)

    # conditions to stop
    if (nrow(df_trait_sp_a)==0) next
    if (all(is.na(df_trait_sp_a$value) | df_trait_sp_a$value == 0)) {
      next
    }

    # produce plot
    v <- names(table(df_trait_sp_a$variable))[c(which(grepl(names(table(df_trait_sp_a$variable)), pattern="Deltasquare") &
                                             c(grepl(names(table(df_trait_sp_a$variable)), pattern="M") | grepl(names(table(df_trait_sp$variable)), pattern="F")) 
                                             & grepl(names(table(df_trait_sp_a$variable)), pattern="Div") & grepl(names(table(df_trait_sp_a$variable)), pattern="deaths_rate")
                                             & !grepl(names(table(df_trait_sp_a$variable)), pattern="colSum")),
                                             which(names(table(df_trait_sp_a$variable)) == "Deltasquarexall_deaths_rate\nDivNumerator_N_1"))]
    v2 <- names(table(df_trait_sp_a$variable))[c(which(grepl(names(table(df_trait_sp_a$variable)), pattern="Deltasquare") & 
                                                        c(grepl(names(table(df_trait_sp_a$variable)), pattern="M") | grepl(names(table(df_trait_sp$variable)), pattern="F"))
                                                      & grepl(names(table(df_trait_sp_a$variable)), pattern="Div") & !grepl(names(table(df_trait_sp_a$variable)), pattern="deaths_rate")
                                                      & !grepl(names(table(df_trait_sp_a$variable)), pattern="colSum")),
                                                which(names(table(df_trait_sp_a$variable)) == "Deltasquarexall_deaths\nDivNumerator_N_1"))]

    df_trait_sp_plot <- df_trait_sp_a[which(df_trait_sp_a$variable %in% v),]
    ssdg <-  add_SSDg(df=df_trait_sp_plot, weight=WEIGHT)
    df_trait_sp_plot <- left_join(df_trait_sp_plot, ssdg, by=c("Fmax", "Fstart", "variable", "Species"), relationship = "many-to-many")
    id_FNatTot <- which(df_trait_sp_plot$variable %in% c("Deltasquarexall_deaths_rate\nDivNumerator_N_1", 
                                                         "DeltasquarexF_deaths_rate\nDivNumerator_N_1", 
                                                         "DeltasquarexMnat_deaths_rate\nDivNumerator_N_1"))
    df_trait_sp_plot$variable <- improve_vriable_names(df_trait_sp_plot$variable)
    t <- improve_vriable_names(v)
    result <- plot_trait_vs_age_simple(df_trait_sp=df_trait_sp_plot, sp_name=sp,
                                       var_names=t, Time1=150, Time2=160, relative_to_ref = relative_to_ref,
                                       scale_factor = 1, weight=WEIGHT, nrow=3, ncol=3,
                                       proportional = proportional, position_legend="auto", facet=T,
                                       chosen_title = paste0(plot_name_ch, " variance of reproductive success per mortality rates sources (lx instead of abundance),\n",
                                                              "averaged over years ",  Time1, " and ", Time2))
    if ("SSDg"%in%colnames(df_trait_sp_plot)){
      df_trait_sp_plot$Fmax <- as.factor(df_trait_sp_plot$Fmax)
      inset.plot <- ggplot(df_trait_sp_plot, group = (variable)) +
        geom_boxplot(aes(y=SSDg, x=Fmax)) +
        facet_wrap(~variable, scales="free") +
        theme(strip.background = element_rect(fill="white"), 
              panel.background = element_rect(fill = "white", colour = "grey50"),
              plot.title = element_text(hjust = 0.5), 
              axis.text.x = element_text(angle = 60, hjust = 1),
              axis.text = element_text(colour = 1, size = 9),
              panel.border = element_rect(colour = "black", fill=NA))
      
      inset.tot.F.Nat[[plot_name]][[sp]][[paste(t, collapse="_")]] <- ggplot(df_trait_sp_plot[id_FNatTot,], group = (variable)) +
        geom_boxplot(aes(y=SSDg, x=Fmax)) +
        facet_wrap(~variable, scales="free", ncol=1) +
        labs(title=sp)+ #, y="SSD\u0262"
        theme(strip.background = element_rect(fill="white"), 
              panel.background = element_rect(fill = "white", colour = "grey50"),
              plot.title = element_text(hjust = 0.5), 
              axis.text.x = element_text(angle = 60, hjust = 1),
              axis.text = element_text(colour = 1, size = 9),
              panel.border = element_rect(colour = "black", fill=NA))
      
      result <- result + theme(legend.position = c(0.53, 0))
      result <- ggdraw() +
        draw_plot(result) +
        draw_plot(inset.plot, x = 0.65, y = 0.0, width = 0.3, height = 0.32)
    }
    plot_Deathsdelta_age[[plot_name]][[sp]][[paste(t, collapse="_")]]<- result
    
    
    df_trait_sp_plot <- df_trait_sp_a[which(df_trait_sp_a$variable %in% v2),]
    ssdg <-  add_SSDg(df=df_trait_sp_plot, weight=WEIGHT)
    df_trait_sp_plot <- left_join(df_trait_sp_plot, ssdg, by=c("Fmax", "Fstart", "variable", "Species"), relationship = "many-to-many")
    id_FNatTot <- which(df_trait_sp_plot$variable %in% c("Deltasquarexall_deaths\nDivNumerator_N_1", 
                                                         "DeltasquarexF_deaths\nDivNumerator_N_1", 
                                                         "DeltasquarexMnat_deaths\nDivNumerator_N_1"))
    df_trait_sp_plot$variable <- improve_vriable_names(df_trait_sp_plot$variable)
    t <- improve_vriable_names(v2)
    result <- plot_trait_vs_age_simple(df_trait_sp=df_trait_sp_plot, sp_name=sp,
                                       var_names=t, Time1=150, Time2=160, relative_to_ref = relative_to_ref,
                                       scale_factor = 1, weight=WEIGHT, nrow=3, ncol=3,
                                       proportional = proportional, position_legend="auto", facet=T,
                                       chosen_title = paste0(plot_name_ch, " variance of reproductive success per mortality sources,\n",
                                                             "averaged over years ",  Time1, " and ", Time2))
    if ("SSDg"%in%colnames(df_trait_sp_plot)){
      df_trait_sp_plot$Fmax <- as.factor(df_trait_sp_plot$Fmax)
      inset.plot <- ggplot(df_trait_sp_plot, group = (variable)) +
        geom_boxplot(aes(y=SSDg, x=Fmax)) +
        facet_wrap(~variable, scales="free") +
        theme(strip.background = element_rect(fill="white"), 
              panel.background = element_rect(fill = "white", colour = "grey50"),
              plot.title = element_text(hjust = 0.5), 
              axis.text.x = element_text(angle = 60, hjust = 1),
              axis.text = element_text(colour = 1, size = 9),
              panel.border = element_rect(colour = "black", fill=NA))
      
      inset.tot.F.Nat[[plot_name]][[sp]][[paste(t, collapse="_")]] <- ggplot(df_trait_sp_plot[id_FNatTot,], group = (variable)) +
        geom_boxplot(aes(y=SSDg, x=Fmax)) +
        facet_wrap(~variable, scales="free", ncol=1) +
        labs(title=sp)+ #, y="SSD\u0262"
        theme(strip.background = element_rect(fill="white"), 
              panel.background = element_rect(fill = "white", colour = "grey50"),
              plot.title = element_text(hjust = 0.5), 
              axis.text.x = element_text(angle = 60, hjust = 1),
              axis.text = element_text(colour = 1, size = 9),
              panel.border = element_rect(colour = "black", fill=NA))
      
      result <- result + theme(legend.position = c(0.53, 0))
      result <- ggdraw() +
        draw_plot(result) +
        draw_plot(inset.plot, x = 0.65, y = 0.0, width = 0.3, height = 0.32)
    }
    plot_Deathsdelta_age2[[plot_name]][[sp]][[paste(t, collapse="_")]]<- result
  }#relative/proportional/none
}#sp

save(plot_Deathsdelta_age, file=paste0(pathoutputRDATA, "/plot_Deathsdelta_age.RData"))
save(plot_Deathsdelta_age2, file=paste0(pathoutputRDATA, "/plot_Deathsdelta_age2.RData"))
deaths_deltasquare_byage <- get(load(file=paste0(pathoutputRDATA, "/plot_Deathsdelta_age.RData")))
deaths_deltasquare_byage2 <- get(load(file=paste0(pathoutputRDATA, "/plot_Deathsdelta_age2.RData")))

save(inset.tot.F.Nat, file=paste0(pathoutputRDATA, "/inset.tot.F.Nat.RData"))
load(file=paste0(pathoutputRDATA, "/inset.tot.F.Nat.RData"))



for (relative_to_ref in c(FALSE)){ #, FALSE, "proportional"
  proportional = FALSE
  if (relative_to_ref == "proportional") {
    relative_to_ref = FALSE
    proportional = TRUE
  }
  relative_to_ref = as.logical(relative_to_ref)
  plot_name_ch=c()
  plot_name <- paste0(as.character(relative_to_ref), as.character(proportional))
  if (plot_name=="TRUEFALSE") plot_name_ch <- c("relative")
  if (plot_name=="FALSEFALSE") plot_name_ch <- c("")
  if (plot_name=="FALSETRUE") plot_name_ch <- c("proportional")

  windowsFonts(myCustomWindowsFontName=windowsFont("DejaVu Sans"))
  cairo_pdf( onefile=T, file = paste0(PATHOUT, "/Clean_for_paper/DeathsDeltasquareonN_2_LX_", plot_name_ch, ".pdf"), width = 15, height = 10)
    p <-(plot_Deathsdelta_age[[plot_name]])
    print(p)
  dev.off()
}



for (relative_to_ref in c(FALSE)){ #, FALSE, "proportional"
  proportional = FALSE
  if (relative_to_ref == "proportional") {
    relative_to_ref = FALSE
    proportional = TRUE
  }
  relative_to_ref = as.logical(relative_to_ref)
  plot_name_ch=c()
  plot_name <- paste0(as.character(relative_to_ref), as.character(proportional))
  if (plot_name=="TRUEFALSE") plot_name_ch <- c("relative")
  if (plot_name=="FALSEFALSE") plot_name_ch <- c("")
  if (plot_name=="FALSETRUE") plot_name_ch <- c("proportional")
  
  windowsFonts(myCustomWindowsFontName=windowsFont("DejaVu Sans"))
  cairo_pdf( onefile=T, file = paste0(PATHOUT, "/Clean_for_paper/DeathsDeltasquareonN_2", plot_name_ch, ".pdf"), width = 15, height = 10)
  p <-(plot_Deathsdelta_age2[[plot_name]])
  print(p)
  dev.off()
}





###### PLOT DELTA²xDEAD/N raw, and comparing the SSDgs -------------------



# # age delta²/N and Death by age
combined_plot =list()
Time1=150
Time2=160
traits_sp_to_plot = list(
                         c("DeltasquarexDead", "DeltasquarexDead\nDivNumerator_N_1", "Dead\nDivNumerator_N_1"))
PATHOUT  = "osmose-ns-v4.4.1_clim/output/PLOTS/KEY_TRAITS/"
path <- c("C:/Users/mbeneat/Documents/osmose/Chap2/osmose-ns/osmose-ns-v4.4.1_clim/output/PLOTS/FINAL_Ne_Htz2")
plots_by_age_delta = list()
for (sp in species){
  cat(sp)
  
  load(file=paste0(path, "/Complete_INIT2", sp, ".RData"))
  final_dataset_age = complete[[2]]
  final_dataset_age <- final_dataset_age %>%
    dplyr::mutate(dplyr::across(c(cohort, replicate), ~ stringr::str_remove(.x, "^(cohort_|rep_)")))
  final_dataset_age <- final_dataset_age[which(final_dataset_age$cohort %in% Time1:Time2),]
  rm(complete)
  
  # remove ages that are not included
  final_dataset_age$age <- as.numeric(final_dataset_age$age)
  ages_to_remove <- table(final_dataset_age[which(final_dataset_age$variable == "Dead" & final_dataset_age$value ==0),"age"])
  ages_to_remove <- ages_to_remove[which(ages_to_remove == max(ages_to_remove))]
  final_dataset_age <- final_dataset_age[-which(final_dataset_age$age %in% as.numeric(names(ages_to_remove))), ]
  
  df_trait_sp <- final_dataset_age %>% filter(Species == sp, Fstart == 1)
  
  # avoid the estimation bias of very populations of very low size
  if (any(((df_trait_sp[which(df_trait_sp$variable == "abundanceByAge"), "value"])<10))){
    df_trait_sp <- pivot_wider(df_trait_sp, names_from = "variable", values_from = "value")
    df_trait_sp <- df_trait_sp[-which((df_trait_sp$abundanceByAge)<10), ]
    df_trait_sp <- pivot_longer(df_trait_sp, cols = -c("age", "replicate", "Fmax", "Fstart", 'Species', "cohort"), names_to = "variable", values_to = "value")
  }
  
  #add delta square
  df_trait_t    <- estimateSquare(trait_data_full_i = df_trait_sp, var=paste0("Deltasquare"))
  df_trait_sp   <- rbind(df_trait_sp, df_trait_t)
  
  #add delta square X Dead
  df_trait_t    <- estimateX(trait_data_full_i = df_trait_sp, var=paste0("DeltasquarexDead"))
  df_trait_sp   <- rbind(df_trait_sp, df_trait_t)
  
  # delta squareXdead / lx
  df_trait_t    <- estimateDiv(trait_data_full_i = df_trait_sp, var=paste0("DeltasquarexDead\nDivsurvivorship_curve"))
  df_trait_sp   <- rbind(df_trait_sp, df_trait_t)
  
  #add sum by age of SSDgx and SSDix
  df_trait_t    <- estimatecolSum(trait_data_full_i = df_trait_sp, var=paste0("colSumSSDix"))
  df_trait_sp   <- rbind(df_trait_sp, df_trait_t)
  df_trait_t    <- estimatecolSum(trait_data_full_i = df_trait_sp, var=paste0("colSumSSDgx"))
  df_trait_sp   <- rbind(df_trait_sp, df_trait_t)
  df_trait_t    <- estimatecolSum(trait_data_full_i = df_trait_sp, var=paste0("colSumDeltasquare"))
  df_trait_sp   <- rbind(df_trait_sp, df_trait_t)
  
  var_i <- unique(df_trait_sp$variable)
  
  # keep abundance by age variable on the side
  WEIGHT = df_trait_sp[grep(df_trait_sp$variable, pattern="abundanceByAge"),]
  
  # compare all traits to their value /N1 to understand where is the increase coming from
  df_trait_t    <- estimatecolSum(trait_data_full_i = df_trait_sp, var=paste0("colSumDeltasquarexDead"))
  df_trait_sp   <- rbind(df_trait_sp, df_trait_t)
  
  for (v in var_i){
    if (v != "Numerator_N_1") {
      df_trait_t    <- estimateDiv(trait_data_full_i = df_trait_sp, var=paste0(v, "\nDivNumerator_N_1")) # div by population size
      df_trait_sp   <- rbind(df_trait_sp, df_trait_t)
    }
    
  }
  
  ssdg <-  add_SSDg(df=df_trait_sp, weight=WEIGHT)
  df_trait_sp <- left_join(df_trait_sp, ssdg, by=c("Fmax", "Fstart", "variable", "Species"), relationship = "many-to-many")
  
  
  var_i <- unique(df_trait_sp$variable)
  
  for (relative_to_ref in c(FALSE)){
    proportional = FALSE
    if (relative_to_ref == "proportional") {
      relative_to_ref = FALSE
      proportional = TRUE
    }
    relative_to_ref = as.logical(relative_to_ref)
    # for the plot
    plot_biom = ifelse(relative_to_ref, "rel_biomassplot", "biomassplot")
    plot_biom = ifelse(relative_to_ref, "biomassplot", "biomassplot")
    
    ref_data <- df_trait_sp %>%
      filter(Fstart == 1 & Fmax==1) %>%
      dplyr::select(variable, value, replicate, age) %>%
      dplyr::rename(ref_value = value)
    df_trait_sp_a <- df_trait_sp
    
    df_trait_sp_a <- clean_values(df_trait_sp_a)
    
    # conditions to stop
    if (nrow(df_trait_sp_a)==0) next
    if (all(is.na(df_trait_sp_a$value) | df_trait_sp_a$value == 0)) {
      next
    }
    
    
    # produce plot
    var_i <- unlist(traits_sp_to_plot)
    df_trait_sp_trait <- df_trait_sp_a[which(df_trait_sp_a$variable %in% var_i),]
    ref_data_split <- ref_data[which(ref_data$variable %in% var_i),]
    if (all(is.na(df_trait_sp_trait$value))) next
    for (t in traits_sp_to_plot){
      df_trait_sp_plot <- df_trait_sp_trait[which(df_trait_sp_trait$variable == t),]
      df_trait_sp_plot$variable <- improve_vriable_names(df_trait_sp_plot$variable)
      t <- improve_vriable_names(t)
      df_trait_sp_plot$variable <- as.factor(df_trait_sp_plot$variable)
      special_layout = NULL
      if ("Δ²Dead" %in% t) {
        special_layout <- c("aa
                             a#")
        factor_order_forced = factor(df_trait_sp_plot$variable,levels=c(levels(df_trait_sp_plot$variable), c('', ' ', '  ')[1:((2*2)-length(levels(df_trait_sp_plot$variable)))]))
        df_trait_sp_plot$variable <- factor_order_forced
      }
      result <- plot_trait_vs_age_simple(df_trait_sp=df_trait_sp_plot, sp_name = sp, var_names= t, #age_dim = T,
                                         Time1=150, Time2=160, relative_to_ref = relative_to_ref,
                                         scale_factor = 1, 
                                         ncol=2, nrow=2, 
                                         weight=WEIGHT, , facet=T, scale_facet="free",
                                         proportional = proportional, position_legend="auto")
      if ("SSDg"%in%colnames(df_trait_sp_plot)){
        df_trait_sp_plot$Fmax <- as.factor(df_trait_sp_plot$Fmax)
        inset.plot <- ggplot(df_trait_sp_plot, group = (variable)) +
          geom_boxplot(aes(y=SSDg, x=Fmax)) +
          facet_wrap(~variable, scales="free") +
          theme(strip.background = element_rect(fill="white"), 
                          panel.background = element_rect(fill = "white", colour = "grey50"),
                          plot.title = element_text(hjust = 0.5), 
                          axis.text.x = element_text(angle = 60, hjust = 1),
                          axis.text = element_text(colour = 1, size = 9),
                          panel.border = element_rect(colour = "black", fill=NA))+
          ylab("")
        result <- result + theme(legend.position = c(0.18, 0.78))
        result <- remove_facets(plot=result, layout = special_layout)
        result <- ggdraw() +
          draw_plot(result) +
          draw_plot(inset.plot, x = 0.53, y = 0.0, width = 0.4, height = 0.4)
      }
      plots_by_age_delta[[paste0(as.character(relative_to_ref), as.character(proportional))]][[sp]][[paste(t, collapse="_")]]<- result
    }
  }#relative/proportional/none
}#sp

save(plots_by_age_delta, file=paste0(pathoutputRDATA, "/plots_by_age_delta.RData"))
plots_by_age_delta <- get(load(file=paste0(pathoutputRDATA, "/plots_by_age_delta.RData")))


for (relative_to_ref in c(FALSE)){ #, FALSE, "proportional"
  proportional = FALSE
  if (relative_to_ref == "proportional") {
    relative_to_ref = FALSE
    proportional = TRUE
  }
  relative_to_ref = as.logical(relative_to_ref)
  plot_name_ch=c()
  plot_name <- paste0(as.character(relative_to_ref), as.character(proportional))
  if (plot_name=="TRUEFALSE") plot_name_ch <- c("relative")
  if (plot_name=="FALSEFALSE") plot_name_ch <- c("absolute")
  if (plot_name=="FALSETRUE") plot_name_ch <- c("proportional")
  
  windowsFonts(myCustomWindowsFontName=windowsFont("DejaVu Sans"))
  cairo_pdf( onefile=T, file = paste0(PATHOUT, "/Clean_for_paper/DeathsDeltasquareCOMPAREwithN", plot_name_ch, ".pdf"), width = 15, height = 10)
  p <-(plots_by_age_delta[[plot_name]])
  print(p)
  dev.off()
}

