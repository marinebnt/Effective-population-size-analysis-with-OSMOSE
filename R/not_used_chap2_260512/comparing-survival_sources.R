# the idea is to store separately all the key parameters per species into plots
# and to produce one pdf sheet per species with all the key variables to have it all summarized at one glance




#### Initialisation --------------
library(legendry)
library(ggpubr)

# functions and environment 
main = "New_gpu_conf.R" # config file
source("R/functions/source.R") # load config and key parameters and functions
source("R/functions/functions_plots.R") # plots against Fmax
source("R/functions/auxiliary_for_lists.R") # from list to df 
# source("R/functions/functions_output_SimuRData_into_df_for_plots.R") # to transform dataset from RData to long df
source("R/functions/function_data_traits_summary.R")

# paths
PATHOUT  = "output/PLOTS/KEY_TRAITS/"
if (!dir.exists(PATHOUT)) dir.create(PATHOUT, recursive = TRUE)
base_dir <- "output/PLOTS/FINAL_Ne_Htz2/INIT"
folders  <- list.dirs(base_dir, recursive = F) # files INIT2 with raw Ne datasets

folders_traits_y <- "output/PLOTS/Trait_variables_by_year_BYSP//RData/"
folders_traits_y <- list.dirs(folders_traits_y, recursive = F)
folders_traits_c <- "output/PLOTS/Trait_variables_by_cohort/RData/"
folders_traits_c  <- list.dirs(folders_traits_c, recursive = F)
folders_traits_y_a <- "output/PLOTS/Trait_variables_by_year_by_age2/RData/"
folders_traits_y_a <- list.dirs(folders_traits_y_a, recursive = F)
folders_traits_c_a <- "output/PLOTS/Trait_variables_by_cohort_by_age2/RData/"
folders_traits_c_a  <- list.dirs(folders_traits_c_a, recursive = F)

# load or format dataset
pathoutputRDATA = paste0(configDir, "/output/PLOTS/KEY_TRAITS/RDATA")
if(!dir.exists(pathoutputRDATA)) dir.create(pathoutputRDATA, recursive = T)
Time1=150
Time2=160
species_df=species





pathout <- paste0("output/PLOTS/Trait_variables_by_cohort_by_age2/otherplots")
if (!dir.exists(pathout)) dir.create(pathout, recursive = T)

# load dataset
data_complete_sp = NULL
data_complete_c = NULL
base_path = "output/PLOTS/Trait_variables_by_cohort_by_age2"
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

PATHOUT  = "output/PLOTS/KEY_TRAITS/"
path <- c("C:/Users/mbeneat/Documents/osmose/Chap2/osmose-ns/output/PLOTS/FINAL_Ne_Htz2")


### produce dataset reading multiple .RData to make 1 dataframe, stored as a .RData file ----------
sp_i = "Hake"
FORCE=F

combined_plot =list()
Time1=150
Time2=160
PATHOUT  = "output/PLOTS/KEY_TRAITS/"
path <- c("C:/Users/mbeneat/Documents/osmose/Chap2/osmose-ns/output/PLOTS/FINAL_Ne_Htz2")
plot_Abund_age = list()
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
  
  df_trait_sp <- full_join(data_complete_c_i,df_trait_sp, by=c("age", "Fmax", "Fstart", "Species", "replicate"))
  if (any(df_trait_sp$age == 0)) df_trait_sp <- df_trait_sp[-which(df_trait_sp$age == 0), ]
  df_trait_sp <- na.omit(df_trait_sp)
  df_trait_sp <- pivot_longer(data = df_trait_sp, cols = -which(colnames(df_trait_sp) %in% c("age", "Fmax", "Fstart", "Species", "replicate")), 
                              names_to = "variable", values_to = "value")
  
  #add lxbx
  df_trait_t    <- estimateX(trait_data_full_i=df_trait_sp, var=paste0("survivorship_curvexb_rate"))
  df_trait_sp   <- rbind(df_trait_sp, df_trait_t)
  df_trait_t    <- estimateX(trait_data_full_i=df_trait_sp, var=paste0("survivorship_curvexb_rate_sc"))
  df_trait_sp   <- rbind(df_trait_sp, df_trait_t)
  
  var_i <- unique(df_trait_sp$variable)

  # colSumDeltasquarexDead
  df_trait_t    <- estimatecolSum(trait_data_full_i = df_trait_sp, var=paste0("colSumDeltasquarexDead"))
  df_trait_sp   <- rbind(df_trait_sp, df_trait_t)  
  
  # Bx = bxNx
  df_trait_t    <- estimateX(trait_data_full_i = df_trait_sp, var=paste0("b_ratexabundanceByAge"))
  df_trait_t$variable <- "B_age"
  df_trait_sp   <- rbind(df_trait_sp, df_trait_t)  
  
  # Bx = agebxNx
  df_trait_t    <- df_trait_sp[which(df_trait_sp$variable == "B_age"), ]
  df_trait_t$value    <- df_trait_t$value[which(df_trait_t$variable == "B_age")]*df_trait_sp$age[which(df_trait_t$variable == "B_age")]
  df_trait_t$variable <- "ageB_age (generation length)"
  df_trait_sp   <- rbind(df_trait_sp, df_trait_t)  
  
  # Generation length
  df_trait_t    <- estimateDiv(trait_data_full_i = df_trait_sp, var=paste0("ageB_age (generation length)\nDivNumerator_N_1"))
  df_trait_sp   <- rbind(df_trait_sp, df_trait_t) 
  df_trait_t    <- estimatecolSum(trait_data_full_i = df_trait_sp, var=paste0("colSumxBx\nDivNumerator_N_1"))
  df_trait_t$variable <- "Generation length"
  df_trait_sp   <- rbind(df_trait_sp, df_trait_t)  
  
  # Nx = lxNx
  df_trait_t    <- estimateX(trait_data_full_i = df_trait_sp, var=paste0("survivorship_curvexNumerator_N_1"))
  df_trait_sp   <- rbind(df_trait_sp, df_trait_t)
  
  # delta square
  df_trait_t    <- estimateSquare(trait_data_full_i = df_trait_sp, var=paste0("Deltasquare"))
  df_trait_sp   <- rbind(df_trait_sp, df_trait_t)
  
  var_i <- unique(df_trait_sp$variable)
  # second kind of traits vector construction
  for (v in var_i){
    if (v %in% c("k_point", "kx_bar", "Dead", "Deltasquare")){
      df_trait_t    <- estimateDiv(trait_data_full_i = df_trait_sp, var=paste0(v, "\nDivNumerator_N_1")) # div by population size
      df_trait_sp   <- rbind(df_trait_sp, df_trait_t)
    }
  }
  
  
  # keep abundance by age variable on the side
  WEIGHT = df_trait_sp[which(df_trait_sp$variable=="abundanceByAge"),]
  
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
    v <- list( 
      unique(names(table(df_trait_sp_a$variable)[c(
        which(names(table(df_trait_sp_a$variable))=="b_rate"), 
        which(names(table(df_trait_sp_a$variable))=="B_age"), 
        # which(names(table(df_trait_sp_a$variable))=="ageB_age (generation length)"), 
        grep(names(table(df_trait_sp_a$variable)), pattern="survivorship"))])),
      
      unique(names(table(df_trait_sp_a$variable)[c(which(names(table(df_trait_sp_a$variable)) == "k_point"), 
                                                   which(names(table(df_trait_sp_a$variable)) == "kx_bar"), 
                                                   which(names(table(df_trait_sp_a$variable)) == "k_point\nDivNumerator_N_1"), 
                                                   which(names(table(df_trait_sp_a$variable)) == "kx_bar\nDivNumerator_N_1"))])),    
      
      unique(names(table(df_trait_sp_a$variable)[c(which(names(table(df_trait_sp_a$variable)) == "Deltasquare"), 
                                                   which(names(table(df_trait_sp_a$variable)) == "Deltasquare\nDivNumerator_N_1"), 
                                                   which(names(table(df_trait_sp_a$variable)) == "Dead"), 
                                                   which(names(table(df_trait_sp_a$variable)) == "Dead\nDivNumerator_N_1"))]))
    )
    # v <- v[-which(v == "survivorship_curvexNumerator_N_1")]
    
    
    if (all(is.na(df_trait_sp_a$value))) next
    for (v_ii in 1:length(v)){
      v_i = v[[v_ii]]
      df_trait_sp_plot <- df_trait_sp_a[which(df_trait_sp_a$variable %in% v_i),]
      t <- improve_vriable_names(v_i)
      df_trait_sp_plot$variable <- improve_vriable_names(df_trait_sp_plot$variable)
      df_trait_sp_plot$variable <- factor(as.factor(df_trait_sp_plot$variable),levels=c(t))
      df_trait_sp_plot <- df_trait_sp_plot %>%
        group_by(variable, Fmax, Species, Fstart) %>% 
        mutate(q99 = quantile(value, 0.99)) %>%
        filter(value <= q99) %>%
        select(-c(q99))
      ncol=ifelse(length(t)>3, 3, 2)
      result <- plot_trait_vs_fmax(df_trait_sp=df_trait_sp_plot, spec=sp, age_dim = T, sp_dim = F,  
                                   trait_name=t, Time1=150, Time2=160, relative_to_ref = relative_to_ref,
                                   # scale_factor = 1, 
                                   weight=WEIGHT, nrow=ceiling(length(t)/ncol), ncol=ncol, 
                                   proportional = proportional, position_legend="auto", facet=T, scale_z=F,
                                   chosen_title = paste0(plot_name_ch, " abundance and related variables per age, per Fmax,\n", 
                                                         "averaged over years ",  Time1, " and ", Time2))
      plot_Abund_age[[plot_name]][[sp]][[paste(t, collapse="_")]]<- result
    }#plot
  }#relative/proportional/none
}#sp

save(plot_Abund_age, file=paste0(pathoutputRDATA, "/plot_Abund_age.RData"))
plot_Abund_age <- get(load(file=paste0(pathoutputRDATA, "/plot_Abund_age.RData")))

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
  for (trai in names(plot_Abund_age[[plot_name]][["Hake"]])){
    windowsFonts(myCustomWindowsFontName=windowsFont("DejaVu Sans"))
    cairo_pdf(file = paste0(PATHOUT, "/Clean_for_paper/Abund_age_lxbx_", plot_name_ch, substr(trai, start = 1, stop=5), ".pdf"), width = 15, height = 10, onefile = T)
    for (sp in species)  {p <-(plot_Abund_age[[plot_name]][[sp]][[trai]])
    print(p)}
    dev.off()
  }
}


