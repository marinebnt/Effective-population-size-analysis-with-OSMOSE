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
if (!dir.exists(paste0(PATHOUT, "/Clean_for_paper_rm_term_age"))) dir.create(paste0(PATHOUT, "/Clean_for_paper_rm_term_age"), recursive = TRUE)
pathoutputRDATA = paste0(PATHOUT, "RDATA")
if(!dir.exists(pathoutputRDATA)) dir.create(pathoutputRDATA, recursive = T)
base_path = paste0("output/DATASETS/Ne_on_NT_with_F/Trait_variables_by_cohort_by_age2")
base_path = paste0("output/DATASETS/Ne_on_NT_with_F/Other_variables/")

folders_traits_c <- paste0("output/DATASETS/Ne_on_NT_with_F/Trait_variables_by_cohort/RData/")
folders_traits_c  <- list.dirs(folders_traits_c, recursive = F)
folders_traits_c_a <- paste0("output/DATASETS/Ne_on_NT_with_F/Trait_variables_by_cohort_by_age/RData/")
folders_traits_c_a <- paste0("output/DATASETS/Ne_on_NT_with_F/Other_variables/RData/")
folders_traits_c_a  <- list.dirs(folders_traits_c_a, recursive = F)


pathdata=paste0("output/DATASETS/Ne_on_NT_with_F/Ne_NT_calculated_variables/") # output from simulations, only selected cohorts
base_dir <- paste0(pathdata, "INIT")
folders  <- list.dirs(paste0(here::here(), "/", base_dir), recursive = F) # files INIT2 with raw Ne datasets
full_dataset = list()
full_dataset_age = list()
dataset="data_Ne"


param_b_l = NULL
for (f_i in folders){ # 374,362
  
  if (!grepl(pattern = "Fstart1_", f_i)) next
  
  Fmax   <- as.numeric(stringr::str_split(stringr::str_split(f_i, pattern="Fmax")[[1]][2], pattern="_Fstart")[[1]][1])
  Fstart <- as.numeric(stringr::str_split(stringr::str_split(f_i, pattern="Fstart")[[1]][2], pattern="_sp")[[1]][1])
  
  if (Fstart != 1) next
  
  cat("Loading   ",  f_i, " Fmax= ", Fmax, "\n")
  
  load(list.files(f_i, full.names = T))
  dataNe = INIT2[[dataset]]
  rm(INIT2)
  sp = names(dataNe)[1]
  
  for (var in c("survival_rate", "survivorship_curve", "b_rate", "Numerator_N_alpha", "Deltasquare", "Vk_point", "Dead")){
    param_b_l_t = NULL 
    ## retrieve variables needed to estimate Vk_point for cohorts 150-160 and all replicates/ages/Fmax (Fstart = 1)
    dataNe[[sp]][[var]]       -> data_t
    
    if (class(data_t) == "list"){
      tempo_l <- lapply(names(data_t), FUN = function(l) {
        tempo <- as.data.frame(data_t[[l]])
        tempo$replicate <- stringr::str_remove(l, pattern = "rep_")
        tempo$age <- rownames(tempo)
        tempo <- pivot_longer(tempo, names_to = "cohort", values_to = "value", cols=dimnames(data_t$rep_1)[[2]])
        tempo$cohort <- stringr::str_remove(tempo$cohort, pattern = "cohort_")
        tempo <- tempo %>% filter(cohort %in% c(150:160))
        return(tempo)
      }
      )
      param_b_l_t <- do.call(rbind, tempo_l)
    } else {
      data_t[, "cohort"] <- rownames(data_t)
      param_b_l_t <- pivot_longer(data_t, cols = paste0("rep_", 1:20), values_to = "value", names_to = "replicate") 
      param_b_l_t$cohort    <- stringr::str_remove(param_b_l_t$cohort, pattern = "cohort_")
      param_b_l_t <- param_b_l_t %>% filter(cohort %in% c(150:160))
      param_b_l_t$replicate <- stringr::str_remove(param_b_l_t$replicate, pattern = "rep_")
      ages <- rownames(dataNe[[sp]][["survival_rate"]]$rep_1)
      param_b_l_tt <- param_b_l_t
      param_b_l_t[,"age"] <- as.numeric(ages[1])
      for (a in ages[-1]){
        param_b_l_tt[,"age"] <- as.numeric(a)
        param_b_l_t <- rbind(param_b_l_t, param_b_l_tt)
      }
    }
    param_b_l_t[, "Fmax"]     <- Fmax
    param_b_l_t[, "Species"]  <- sp
    param_b_l_t[, "variable"] <- var
    
    if (is.null(param_b_l_t)) {param_b_l <- param_b_l_t
    } else {param_b_l <- rbind(param_b_l_t, param_b_l)}
  }
}

p_t <- pivot_wider(param_b_l, values_from = "value", names_from = "variable")
p_t[which(p_t$Dead == 0), c("Vk_point", "Numerator_N_alpha")] <- NA      
p_t <- pivot_longer(p_t, cols = -c("replicate", "age",   "cohort",  "Fmax", "Species"), values_to = "value", names_to = "variable")
param_b_l <- p_t

### *** remove age corresponding to the max age *** ----------------
lifespan <- as.data.frame(as.matrix(get_par(conf, "lifespan")))
colnames(lifespan) <- "lifespan"
lifespan$ID <- stringr::str_remove(pattern="species.lifespan.sp", rownames(lifespan))
rownames(lifespan) <- NULL
name <- as.data.frame(as.matrix(get_par(conf, "species.name")))
colnames(name) <- "Species"
name$ID <- stringr::str_remove(pattern="species.name.sp", rownames(name))
rownames(name) <- NULL
lifespan <- full_join(name, lifespan)
lifespan <- lifespan[which(lifespan$Species %in% species), ]
lifespan$Species <- as.character(lifespan$Species)
lifespan$lifespan <- as.numeric(lifespan$lifespan)-1
param_b_l <- full_join(param_b_l, lifespan)

param_b_l <- param_b_l[which(is.na(param_b_l$lifespan) | as.numeric(param_b_l$age) != param_b_l$lifespan),
                       !(colnames(param_b_l) %in% c("lifespan", "ID"))]            # [CORRIGÉ] no longer empties the table if nothing matches + ID column removed
### --------------------------------------------------------------------


# Counterfactual decomposition --------------

### prepare dataset : re-estimate Vk using bx, bx_ref, kx or lx_ref ----
data_counter <- NULL
param_b_l_calculations   <- param_b_l
param_b_l_calculations_t <- param_b_l_calculations %>% filter(variable == "survivorship_curve")
param_b_l_calculations_t$value[is.na(param_b_l_calculations_t$value)] <- 0   
param_b_l_calculations_t$age <- as.numeric(param_b_l_calculations_t$age) - 1
param_b_l_calculations_t$variable <- "survivorship_curve+1"
param_b_l_calculations <- rbind(param_b_l_calculations_t, param_b_l_calculations)

# loop to estimate Vk_point with different sets of bx and lx as input
for (bx_lx in list(c("b_rate", "survivorship_curve"), c("b_rate", "survivorship_curveref"),
                   c("b_rateref", "survivorship_curve"), c("b_rateref", "survivorship_curveref"))){
  
  bx_n <- bx_lx[1]
  lx_n <- bx_lx[2]
  lx_n_1 <- paste0(lx_n, "+1")
  if (grepl("ref", lx_n)) lx_n_1 <- stringr::str_replace(lx_n, pattern = "ref", replacement = "+1ref")
  
  store_Fm = NULL
  
  for (Fm in c(1.5, 1.25, 1, 0.75, 0.5)){
    ## section that depends on lx and bx
    cat(Fm)
    
    data_Fm <- param_b_l_calculations %>% filter(Fmax == Fm)
    data_Fref <- param_b_l_calculations %>% filter(Fmax == 1)
    
    # reference data: merge it with the dataset
    data_Fref$Fmax <- Fm
    data_Fref$variable <- paste0(data_Fref$variable, "ref")
    data_Fm <- rbind(data_Fref, data_Fm)
    data_Fm <- pivot_wider(data_Fm, values_from = "value", names_from = "variable")
    data_Fm$age <- as.numeric(data_Fm$age)
    
    # lx+1 is missing the last age value because I removed it by removing terminal age of lx
    # => fill in the gaps with 0
    data_Fm$`survivorship_curve+1`[which(is.na(data_Fm$`survivorship_curve+1`))] <- 0
    
    # remove missing values (doc, section 7.2)
    data_Fm <- data_Fm %>% filter(!is.na(.data[[lx_n]]))            # âges de la source L uniquement
    data_Fm[[lx_n_1]][is.na(data_Fm[[lx_n_1]])] <- 0                # l_(omega+1) = 0
    # add N_alpha only if there is a value for the other variables
    data_Fm <- data_Fm %>%
      group_by(replicate, cohort, Fmax, Species) %>%
      mutate(Numerator_N_alpha = dplyr::coalesce(dplyr::first(na.omit(Numerator_N_alpha)),
                                                 dplyr::first(na.omit(Numerator_N_alpharef)))) %>%
      ungroup()
    
    #scaling
    data_Fm <- data_Fm %>%
      group_by(replicate, cohort, Fmax, Species) %>%
      arrange(as.numeric(age), .by_group = TRUE) %>%           # ordre des âges indispensable avant cumsum
      mutate(rescaling_f = sum(.data[[bx_n]] * .data[[lx_n]])) %>%   # somme sur les âges (doc, étape 3)
      ungroup()
    data_Fm$birthrate_scaled <- data_Fm[[bx_n]] / data_Fm$rescaling_f
    
    # kx
    data_Fm <- data_Fm %>%
      group_by(replicate, cohort, Fmax, Species)%>%
      mutate(kx_bar = cumsum(birthrate_scaled))%>%
      ungroup()    
    data_Fm$Vx     <- data_Fm$kx_bar      
    
    # Mean lifetime repro success overall
    data_Fm$w_x <- data_Fm[[lx_n]] - data_Fm[[lx_n_1]]
    data_Fm[, "dead*kx_bar/pop_size"]    <- data_Fm$kx_bar*data_Fm$w_x
    data_Fm$`dead*kx_bar/pop_size`[is.nan(data_Fm$`dead*kx_bar/pop_size`)] <- NA
    data_Fm <- data_Fm %>%
      group_by(replicate, cohort, Fmax, Species)%>%
      mutate(k_point = sum(`dead*kx_bar/pop_size`))%>%   # NA si la cohorte est invalide
      ungroup() 
    
    # generation length  
    data_Fm$abundanceByAge    <- data_Fm$Numerator_N_alpha*data_Fm[[lx_n]]
    data_Fm$Dead     <- data_Fm$Numerator_N_alpha*data_Fm$w_x            
    data_Fm$eggNb    <- data_Fm$abundanceByAge*data_Fm$birthrate_scaled 
    data_Fm$age <- as.numeric(data_Fm$age)
    data_Fm[, "egg*ages/pop_size"] <- data_Fm$eggNb*data_Fm$age/data_Fm$Numerator_N_alpha
    data_Fm <- data_Fm %>%
      group_by(replicate, cohort, Fmax, Species)%>%
      mutate(generation_length = sum(`egg*ages/pop_size`))%>%
      ungroup() 
    
    # Mean square deviation lifetime repro success overall 
    data_Fm$Delta    <- data_Fm$kx_bar - data_Fm$k_point
    data_Fm$Deltasquare <- data_Fm$Delta^2
    data_Fm$SSDix    <- data_Fm$Vx * data_Fm$Numerator_N_alpha * data_Fm$w_x          
    data_Fm$SSDgx    <- data_Fm$Numerator_N_alpha * data_Fm$w_x * data_Fm$Delta^2
    data_Fm$SSD      <- data_Fm$SSDix+data_Fm$SSDgx
    data_Fm$SSD[is.nan(data_Fm$SSD)] <- NA
    data_Fm <- data_Fm %>%
      group_by(replicate, cohort, Fmax, Species)%>%
      mutate(Vk_point = sum(SSD/Numerator_N_alpha))%>%     
      ungroup() 
    
    rownames(data_Fm) <- NULL 
    rownames(store_Fm) <- NULL 
    
    if (is.null(store_Fm)) {store_Fm <- as.data.frame(data_Fm)
    } else {store_Fm <- rbind(as.data.frame(data_Fm), store_Fm)}
  }    
  data_counter[[paste0(bx_lx, collapse="__")]] <- store_Fm[, c("Vk_point", "generation_length", "Dead", "Deltasquare",
                                                               "Numerator_N_alpha", "survivorship_curve",
                                                              "replicate", "cohort", "Fmax", "Species", "age")]
}

# ***  TEST : compare Vk_point (i) estmated from loop of Ne estimation and (ii) recalculated with the bx and lx permutation test for Counterfactual decomposition  ***
old_Vk <- pivot_wider(param_b_l[which(param_b_l$variable == c("Vk_point")),], values_from = "value", names_from = "variable")
new_Vk <- as.data.frame(data_counter$b_rate__survivorship_curve)[,which(colnames(data_counter$b_rate__survivorship_curve) %in% colnames(old_Vk))]
old_Vk$age <- as.numeric(old_Vk$age)
new_Vk$age <- as.numeric(new_Vk$age)
colnames(new_Vk) <- stringr::str_replace(colnames(new_Vk), pattern = "Vk_point", replacement = "recalcu_Vk_point")
compare_Vk <- full_join(old_Vk, new_Vk)

plot(compare_Vk$Vk_point, compare_Vk$recalcu_Vk_point) # 1:1 line well represented. Only diff in new dataset : I have put Vk_point=0 when bx = 0

### prepare dataset : reshape Vk estimated from bx, bx_ref, kx or lx_ref into dataframe ----
data_counter_df <- lapply(names(data_counter), FUN = function(l) {
  tempo <- as.data.frame(data_counter[[l]][, c("Vk_point", "replicate", "cohort", "Fmax", "Species", "age")])
  tempo$bxlx <- l
  return(tempo)
})
data_counter_df <- do.call(rbind, data_counter_df)
data_counter_df$Vk_point <- as.numeric(data_counter_df$Vk_point)
data_counter_df$Vk_point[is.nan(data_counter_df$Vk_point)] <- NA
# Vk_point is the same for all ages within a cohort: one row is retained per cohort
data_counter_df <- data_counter_df[, colnames(data_counter_df) != "age"] %>% distinct()
data_counter_df <- pivot_wider(data_counter_df, values_from = "Vk_point", names_from = "bxlx")

### prepare dataset : decomposition of Vk (doc, Eq. 8-10) ----
data_counter_df$Elx = 0.5*(data_counter_df$b_rateref__survivorship_curve-
                             data_counter_df$b_rateref__survivorship_curveref)+
  0.5*(data_counter_df$b_rate__survivorship_curve-
         data_counter_df$b_rate__survivorship_curveref)
data_counter_df$Ebx = 0.5*(data_counter_df$b_rate__survivorship_curveref-
                             data_counter_df$b_rateref__survivorship_curveref)+
  0.5*(data_counter_df$b_rate__survivorship_curve-
         data_counter_df$b_rateref__survivorship_curve)
data_counter_df$Ix  = data_counter_df$b_rate__survivorship_curve - data_counter_df$b_rateref__survivorship_curve -
  data_counter_df$b_rate__survivorship_curveref + data_counter_df$b_rateref__survivorship_curveref

# TEST (doc 7.4): El + Eb = DeltaVk for each cohort -> must display TRUE
cat("\nControle El + Eb = DeltaVk :",
    isTRUE(all.equal(data_counter_df$Elx + data_counter_df$Ebx,
                     data_counter_df$b_rate__survivorship_curve - data_counter_df$b_rateref__survivorship_curveref)), "\n")

### prepare dataset : average between cohorts ----
# average between cohorts (no sd estimated here because we are at the cohort level)
df_summarised <- data_counter_df %>%
  group_by(replicate, Fmax, Species) %>%
  summarise(
    El=mean(Elx, na.rm = TRUE), 
    Eb=mean(Ebx, na.rm = TRUE), 
    I =mean(Ix,  na.rm = TRUE),                           # [AJOUT]
    V1=mean(b_rateref__survivorship_curveref, na.rm=T),
    Vk_point=mean(b_rate__survivorship_curve, na.rm=T)) %>%
  ungroup()

### prepare dataset : relative to Vk_point(thetaF=1) ---- 
# COMMON denominator: the average of V(1) across the replicates (doc, Eqs. 6–7),
df_summarised <- df_summarised %>%
  group_by(Fmax, Species) %>%
  mutate(V1bar = mean(V1, na.rm = TRUE)) %>%
  ungroup()
df_summarised[, "Eb'"] <- df_summarised$Eb/df_summarised$V1bar
df_summarised[, "El'"] <- df_summarised$El/df_summarised$V1bar
df_summarised[, "I'"]  <- df_summarised$I /df_summarised$V1bar
df_summarised[, "Vk_point'"] <- (df_summarised$Vk_point - df_summarised$V1) / df_summarised$V1bar

# TEST (doc 7.4): Eb' + El' = Vk_point' -> must display TRUE
cat("Controle Eb' + El' = Vk_point' :",
    isTRUE(all.equal(df_summarised$`Eb'` + df_summarised$`El'`, df_summarised$`Vk_point'`)), "\n")

### FIGURE ---- 
df_lollipop <- pivot_longer(df_summarised[, which(colnames(df_summarised) %in% c("replicate", "Fmax", 
                                                                                 "Species", "Eb'",  "El'", "Vk_point'"))], 
                            cols=c("Eb'",  "El'", "Vk_point'"), 
                            values_to = "value", names_to = "variable")
df_lollipop$Fstart <- 1
df_lollipop$Fmax <- as.numeric(df_lollipop$Fmax)
plots_counter<-function_lollipop_decompositions(data = df_lollipop[-which(df_lollipop$Species %in% c("AtlanticMackerel")), ], #, "Sandeel")), ], 
                                                #sp_log = c("Atlantic mackerel", "Sandeel"),
                                                variables_to_bars = c("Eb'",  "El'"), 
                                                variables_to_point = c("Vk_point'"), 
                                                colors = c("#5EB1BF", "#EF7B45"), 
                                                val_names = c("Fecundity relative effect", "Mortality relative effect"))

ggsave(plots_counter$lollipop,     file = paste0(PATHOUT, "/Clean_for_paper_rm_term_age/counter_main.png"), width = 7, height = 7)
ggsave(plots_counter$lollipop_rel, file = paste0(PATHOUT, "/Clean_for_paper_rm_term_age/counter_supp.png"), width = 12, height = 7)




# Kitagawa decomposition (figure B) --------------

# def variables to work on
weight       = c("lx-lx+1")                       # <=> Dx/Nalpha
weight_verif = c("Dead\nDivNumerator_N_alpha")    # <=> Dx/Nalpha
variable     = c("Deltasquare")
keys_coh     = c("Species", "Fmax", "Fstart", "replicate", "cohort")  

# select data from data_Ne
data_kitagawa <- pivot_longer(data_counter$b_rate__survivorship_curve, cols=-c("replicate", "cohort", "Species", "age", "Fmax"),
                              values_to = "value", names_to="variable")
data_kitagawa$Fstart <- 1
data_kitagawa$age    <- as.numeric(data_kitagawa$age)

# estimate lx-lx+1
data_lx_lx1 <- data_kitagawa %>%
  filter(variable == "survivorship_curve", !is.na(value)) %>%
  group_by(replicate, Species, Fmax, Fstart, cohort) %>%  
  arrange(age, .by_group = TRUE) %>%
  mutate(value = value - dplyr::lead(value, default = 0)) %>%
  ungroup()
data_lx_lx1$variable <- c("weight")

## *** TEST => Compare lx-lx+1 and Dx/Nalpha : should be equal ****
data_DxNalpha_t <- data_kitagawa %>%
  filter(variable %in% c("Dead", "Numerator_N_alpha")) %>%
  pivot_wider(names_from = "variable", values_from = "value") %>%
  mutate(DxNalpha = Dead / Numerator_N_alpha)
test_w <- inner_join(data_lx_lx1, data_DxNalpha_t[, c("replicate", "Species", "Fmax", "Fstart", "age", "DxNalpha", "cohort")], 
                     by = c("replicate", "Species", "Fmax", "Fstart", "age", "cohort"))
cat("\nTest lx-lx+1 = Dx/Nalpha : ecart max =", max(abs(test_w$value - test_w$DxNalpha), na.rm = TRUE), "\n")
plot(test_w$value, test_w$DxNalpha)

# select variables & other fixes dataset
data_kitagawa <- rbind(data_lx_lx1,
                       data_kitagawa)   # 3 variables : weight, Vk_point and variable

# calculation of Px and Qx
data_kitagawa_loop <- data_kitagawa
data_kitagawa_loop <- pivot_wider(data_kitagawa_loop, names_from = c("Fmax", "variable"), values_from = c("value"))
data_kitagawa_plot_tt <- NULL
for (Fmaxi in c(0.5, 0.75, 1, 1.25, 1.5)){
  for (sp in unique(data_kitagawa$Species)){
    data_kitagawa_plot_t <- NULL
    # select sp and Fmax of the dataset
    data_kitagawa_loopi <- data_kitagawa_loop[which(data_kitagawa_loop$Species == sp), ]
    wF <- data_kitagawa_loopi[[paste0(Fmaxi, "_weight")]];  w1 <- data_kitagawa_loopi[[paste0(1, "_weight")]]
    dF <- data_kitagawa_loopi[[paste0(Fmaxi, "_", variable)]]; d1 <- data_kitagawa_loopi[[paste0(1, "_", variable)]]
    # cleaning
    dF[is.na(dF)] <- d1[is.na(dF)]
    d1[is.na(d1)] <- dF[is.na(d1)]
    wF[is.na(wF)] <- 0
    w1[is.na(w1)] <- 0
    # calculations of Px and Qx
    data_kitagawa_plot_t <- data.frame("Px" = (wF - w1) * (dF + d1) / 2)
    data_kitagawa_plot_t[, "Qx"] <- (dF - d1) * (wF + w1) / 2
    # storage with the dataset of the right shape
    data_kitagawa_plot_t$Species   <- sp
    data_kitagawa_plot_t$Fmax      <- Fmaxi
    data_kitagawa_plot_t$replicate <- data_kitagawa_loopi$replicate
    data_kitagawa_plot_t$cohort    <- data_kitagawa_loopi$cohort          # [AJOUT]
    data_kitagawa_plot_t$age       <- data_kitagawa_loopi$age
    data_kitagawa_plot_t$Fstart    <- data_kitagawa_loopi$Fstart
    if (is.null(data_kitagawa_plot_tt)) {data_kitagawa_plot_tt <- data_kitagawa_plot_t
    } else data_kitagawa_plot_tt <- rbind(data_kitagawa_plot_tt, data_kitagawa_plot_t)
  }
}
data_kitagawa_plot <- data_kitagawa_plot_tt[!(is.na(data_kitagawa_plot_tt$Px) & is.na(data_kitagawa_plot_tt$Qx)), ]
# organizing dataset to have the right variables
data_kitagawa_plot      <- pivot_longer(data_kitagawa_plot, cols=c("Px", "Qx"), values_to = c("value"), names_to = c("variable"))
# calculations of Sum(Px), Sum(Qx) : estimatecolSum somme par cohorte (colonne cohort présente)
data_kitagawa_sums      <- rbind(estimatecolSum(trait_data_full_i = data_kitagawa_plot, var = paste0("colSumPx")),
                                 estimatecolSum(trait_data_full_i = data_kitagawa_plot, var = paste0("colSumQx")))
data_kitagawa_sums      <- pivot_wider(data_kitagawa_sums, names_from = "variable", values_from = "value")

# DeltaVk_point = Vk_point(thetaF) - Vk_point(1) observed in data_Ne, per cohort
data_vk <- data_kitagawa %>% filter(variable == "Vk_point") %>% distinct(across(all_of(keys_coh)), value)
data_vk <- left_join(data_vk, data_vk %>% filter(Fmax == 1) %>% select(-Fmax) %>% rename(Vk_point_ref = value),
                     by = c("Species", "Fstart", "replicate", "cohort")) %>%
  rename(Vk_point = value)
data_kitagawa_sums <- left_join(data_kitagawa_sums, data_vk, by = keys_coh) %>%
  mutate(DeltaVk_point = Vk_point - Vk_point_ref)

# make sure missing delta Vk are also missing in colsumPx and colsumQx
data_kitagawa_sums[which(is.na(data_kitagawa_sums$DeltaVk_point)), c("colSumPx", "colSumQx")] <- NA          # [AJOUT] as in F_decomposition_age_related_figures2.R
data_kitagawa_sums <- data_kitagawa_sums[which(!is.na(data_kitagawa_sums$colSumPx + data_kitagawa_sums$colSumQx)), ]   # [AJOUT]

## *** TEST => Compare calculated and observed DeltaVk_point : should be equal (doc 6.5)
cat("Test Kitagawa colSumPx + colSumQx = DeltaVk_point : ecart max =",
    max(abs(data_kitagawa_sums$colSumPx + data_kitagawa_sums$colSumQx - data_kitagawa_sums$DeltaVk_point), na.rm = TRUE), "\n")
# not perfect, but looks good on figure : 
plot(data_kitagawa_sums$colSumPx + data_kitagawa_sums$colSumQx, data_kitagawa_sums$DeltaVk_point)

# relative data values
# average across the cohorts per replicate, then divided by the COMMON average of Vk_point(1)
# on the replicates (doc, Eqs. 6–7), rather than the Vk_point(1) for each replicate
data_kitagawa_rel <- data_kitagawa_sums %>%
  group_by(Species, Fmax, Fstart, replicate) %>%
  summarise(across(c(colSumPx, colSumQx, DeltaVk_point, Vk_point_ref), ~ mean(.x, na.rm = TRUE)), .groups = "drop") %>%
  group_by(Species, Fmax, Fstart) %>%
  mutate(V1bar = mean(Vk_point_ref, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(colSumPx_rel = colSumPx / V1bar, colSumQx_rel = colSumQx / V1bar, DeltaVk_point_rel = DeltaVk_point / V1bar)

# lollipop ------------
data_lolli_plot <- pivot_longer(data_kitagawa_rel[, c("Species", "Fmax", "Fstart", "replicate",
                                                      "colSumPx_rel", "colSumQx_rel", "DeltaVk_point_rel")],
                                cols = c("colSumPx_rel", "colSumQx_rel", "DeltaVk_point_rel"),
                                values_to = "value", names_to = "variable")
plots <- function_lollipop_decompositions(data = data_lolli_plot, 
                                          sp_log = c("Atlantic mackerel"), 
                                          variables_to_point = paste0("DeltaVk_point_rel"),
                                          variables_to_bars = c("P'", "Q'"), 
                                          colors = c("palegreen3", "brown2"), 
                                          val_names = c("Dx/N\u03B1 relative effect", "\u0394² relative effect"))

# SAVE
ggsave(plots$lollipop,     file=paste0(PATHOUT, "/Clean_for_paper_rm_term_age/Kitagawa_main.png"), width = 7, height=7)
ggsave(plots$lollipop_rel, file=paste0(PATHOUT, "/Clean_for_paper_rm_term_age/Kitagawa_supp.png"), width=10.5, height=7)





#### PLOT - T/Vk VS FMAX + Biomass VS FMAX ------------------
Ne_plot_list = list()
traits_to_plot <- list(c("generation_length", "Vk_point", "T/Vk."))
dataset="data_parameters"
sp_main = (c("Whiting", "EuropeanPlaice", "AtlanticMackerel", "HorseMackerel"))
sp_main_imp = improve_sp_names(sp_main)
species_imp = improve_sp_names(species)

# select data_parameters
data_plot_1 <- data_counter$b_rate__survivorship_curve 
data_plot_1 <- as_long_if_needed(data_plot_1, columns_keep_pivot = 
                                   colnames(data_plot_1)[-which(colnames(data_plot_1) %in% 
                                      c("replicate", "cohort", "Fmax", "Fstart", "Species", "age"))])

# clean dataset
data_plot_1 <- data_plot_1 %>% filter(variable %in% c(unlist(traits_to_plot), "b_rate", "SSDgx", "survivorship_curve", "Deltasquare", "Dead", "Mtot"))
data_plot_1$Fstart <- 1

data_plot_1 <- data_plot_1[which(data_plot_1$variable %in% 
                                   c(unlist(traits_to_plot), names(table(data_plot_1$variable))[grep(names(table(data_plot_1$variable)), 
                                                                                                     pattern="colSum")])), ]
# remove the data points that only have age that is different
# (as most variables have the same value for all ages)
data_plot_1 <- data_plot_1[, colnames(data_plot_1) != "age"] %>% distinct()

# average acros cohorts
data_plot_1 <- data_plot_1 %>%
  group_by(replicate, Fmax, Fstart, Species, variable) %>%
  summarise(value=mean(value, na.rm=T))%>%
  ungroup()

# add the "T/Vk."                                                                  # [MODIFIÉ] moved after cohort averaging, as in F_decomposition_age_related_figures2.R
df_trait_t    <- estimateDiv(trait_data_full_i = data_plot_1, var=paste0("generation_length\nDivVk_point"))
df_trait_t$variable <-c("T/Vk.")
data_plot_1       <- rbind(data_plot_1, df_trait_t)

# condition to stop
if (nrow(data_plot_1)<1) next

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
  for (traits_selected in traits_to_plot[1]) {
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
    scale_z_i <- ifelse(any(grepl(traits_selected, pattern="b_rate")), T, F)
    
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

for (relative_to_ref in c(TRUE)){
  if (!dir.exists(paste0(PATHOUT, "/Clean_for_paper_rm_term_age"))) dir.create(paste0(PATHOUT, "/Clean_for_paper_rm_term_age"))
  proportional = FALSE
  if (relative_to_ref == "proportional") {
    relative_to_ref = FALSE
    proportional = TRUE
  }
  relative_to_ref <- as.logical(relative_to_ref)
  
  plot_name <- paste0(as.character(relative_to_ref), as.character(proportional))
  if (plot_name=="TRUEFALSE") plot_name_ch <- c("relative")
  if (plot_name=="FALSEFALSE") plot_name_ch <- c("")
  if (plot_name=="FALSETRUE") plot_name_ch <- c("proportional")
  
  plot_biom = ifelse(relative_to_ref, "rel_biomassplot", "biomassplot")
  plot_biom = ifelse(relative_to_ref, "biomassplot", "biomassplot")
  suffix_relative <- plot_name_ch
  
  windowsFonts(myCustomWindowsFontName=windowsFont("DejaVu Sans"))
  cairo_pdf(file = paste0(PATHOUT, "/Clean_for_paper_rm_term_age/recalcuVk_T_TonVk", plot_name_ch, ".pdf"),
            width = 7, height =6, onefile=T)
  print(Ne_plot_list[[plot_name]]$`WhitingEuropean plaiceAtlantic mackerelHorse mackerel`$`generation_lengthVk_pointT/Vk.`)
  dev.off()
  
  windowsFonts(myCustomWindowsFontName=windowsFont("DejaVu Sans"))
  cairo_pdf(file = paste0(PATHOUT, "/Clean_for_paper_rm_term_age/recalcuTOTVk_T_TonVk", plot_name_ch, ".pdf"),
            width = 10, height =12, onefile=T)
  print(Ne_plot_list[[plot_name]]$AtlanticHerringAtlanticMackerelSandeelEuropeanSpratNorwayPoutEuropeanPlaiceCommonSoleSaitheAtlanticCodHaddockHorseMackerelWhitingCommonDabGreyGurnardHake$`generation_lengthVk_pointT/Vk.`)
  dev.off()

}

