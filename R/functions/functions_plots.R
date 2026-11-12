
####---------------------------------------------------------#
### Functions to plot the traits variation with Fmax or age
####---------------------------------------------------------#


# TOOLS functions ----------------

# function of the partial derivative of the relative difference of response to fishing mortality
# g(x, r) = (x - r)/r, x = mean at Fmax (mean_Fmax_i), r = mean at Fmax = 1 (mean_Fmax_ref)
# 'wrt' : derivative with respect to x ("x") or at the reference r ("ref").
# Both derivatives are evaluated on the same line (same x), for each Fmax.
relative_diff_derivated_f <- function(Fmax, mean_Fmax_ref, mean_Fmax_i, wrt = c("x", "ref")){
  wrt <- match.arg(wrt)
  # Fmax == 1 : x et r sont les mêmes réplicats, g = 0 exactement -> pas d'incertitude
  if (Fmax == 1) return(0 * mean_Fmax_ref)
  if (wrt == "x")   return(1/mean_Fmax_ref)                 # dg/dx = 1/r
  if (wrt == "ref") return(mean_Fmax_i/mean_Fmax_ref^2)     # |dg/dr| = x/r^2 
}

# function of the second partial derivative of the relative difference of response to fishing mortality
relative_diff_second_derivate_f <- function(Fmax, mean_Fmax_ref, mean_Fmax_i, wrt = c("x", "ref")){
  wrt <- match.arg(wrt)
  if (Fmax == 1) return(0 * mean_Fmax_ref)
  if (wrt == "x")   return(0 * mean_Fmax_ref)               # d2g/dx2 = 0
  if (wrt == "ref") return(2*mean_Fmax_i/mean_Fmax_ref^3)   # d2g/dr2 = 2x/r^3
}

# Delta method to estimate the standard deviation of the relative difference of response to fishing mortality when independent replicates
# + Esperance to estimate the mean of the relative difference of response to fishing mortality when independent replicates
value_to_sd_esperance_delta_method_relative <- function(dataset_w_value,  # migariane.github.io/DeltaMethodEpiTutorial.nb.html
                                               var_x, h = 1e-6){
  dataset_mean_sd <- NULL
  for (Fmax_i in names(table(dataset_w_value$Fmax))){
    dataset_Fmax <- dataset_w_value[which(dataset_w_value$Fmax==Fmax_i),]
    
    # mean and sd estimators + partial derivative of the
    if (length(var_x)==1){ # depending on the number of variables
      dataset_Fmax <- dataset_Fmax %>%
        group_by(get(var_x), Species, Fstart, variable) %>%
        summarise(mean = mean(value, na.rm = TRUE),
                  nbrep = length(value),
                  var   = var(value, na.rm = TRUE)/nbrep) %>% # we want the the variance of the mean : the uncertainty of `value`
        ungroup()
    }
    if (length(var_x) == 2){ # depending on the number of variables
      dataset_Fmax <- dataset_Fmax %>%
        group_by(get(var_x[1]),get(var_x[2]), Species, Fstart, variable) %>%
        summarise(mean = mean(value, na.rm = TRUE),
                  nbrep = length(value),
                  var   = var(value, na.rm = TRUE)/nbrep) %>% # we want the the variance of the mean : the undertainty of `value`
        ungroup()
    }
    # dataset with mean and sd by fishing pressure
    if (is.null(dataset_mean_sd)) {dataset_mean_sd <- dataset_Fmax
    } else dataset_mean_sd <- rbind(dataset_Fmax, dataset_mean_sd)
  }
  
  # add the case where Fmax = 1 as a reference value of mean and sd
  if (length(var_x)==1){
    colnames(dataset_mean_sd)[which(colnames(dataset_mean_sd) == "get(var_x)")] <- var_x
    data_ref_tempo <- dataset_mean_sd[which(dataset_mean_sd$Fmax == 1), ]
    colnames(data_ref_tempo)[which(colnames(data_ref_tempo) == "mean")] <- c("mean_ref")
    colnames(data_ref_tempo)[which(colnames(data_ref_tempo) == "var")] <- c("var_ref")
    dataset_mean_sd <- full_join(dataset_mean_sd, data_ref_tempo[, c("Fstart", "variable", "Species", "mean_ref", "var_ref")])
  }
  if (length(var_x) == 2){
    colnames(dataset_mean_sd)[which(colnames(dataset_mean_sd) == "get(var_x[1])")] <- var_x[1]
    colnames(dataset_mean_sd)[which(colnames(dataset_mean_sd) == "get(var_x[2])")] <- var_x[2]
    data_ref_tempo <- dataset_mean_sd[which(dataset_mean_sd$Fmax == 1), ]
    colnames(data_ref_tempo)[which(colnames(data_ref_tempo) == "mean")] <- c("mean_ref")
    colnames(data_ref_tempo)[which(colnames(data_ref_tempo) == "var")] <- c("var_ref")
    dataset_mean_sd <- full_join(dataset_mean_sd, data_ref_tempo[, c(var_x[2], "Fstart", "variable", "Species", "mean_ref", "var_ref")])
  }
  
  # dgdx and dg²dx + add the case where Fmax = 1 as a reference value
  dataset_dgdx <- NULL
  for (Fmax_i in names(table(dataset_w_value$Fmax))) {
    data_ref_tempo <- dataset_mean_sd[which(dataset_mean_sd$Fmax == Fmax_i), ]
    data_ref_tempo$dgdx            = relative_diff_derivated_f(Fmax=Fmax_i, wrt = "x",
                                                               mean_Fmax_ref = data_ref_tempo$mean_ref, 
                                                               mean_Fmax_i = data_ref_tempo$mean)
    data_ref_tempo$dgdx_ref        = relative_diff_derivated_f(Fmax=Fmax_i, wrt = "ref",
                                                               mean_Fmax_ref = data_ref_tempo$mean_ref, 
                                                               mean_Fmax_i = data_ref_tempo$mean)
    data_ref_tempo$dgdx_square     = relative_diff_second_derivate_f(Fmax=Fmax_i, wrt = "x",
                                                                     mean_Fmax_ref = data_ref_tempo$mean_ref, 
                                                                     mean_Fmax_i = data_ref_tempo$mean)
    data_ref_tempo$dgdx_square_ref = relative_diff_second_derivate_f(Fmax=Fmax_i, wrt = "ref",
                                                                     mean_Fmax_ref = data_ref_tempo$mean_ref, 
                                                                     mean_Fmax_i = data_ref_tempo$mean)
    if (is.null(dataset_dgdx)) {dataset_dgdx <- data_ref_tempo
    } else dataset_dgdx <- rbind(data_ref_tempo, dataset_dgdx)
  }
  dataset_mean_sd <- dataset_dgdx
  
  # from there : estimate the partial derivative and the following estimate of the standard deviation
  dataset_mean_sd$sd  = sqrt(dataset_mean_sd$dgdx_ref^2*dataset_mean_sd$var_ref+dataset_mean_sd$dgdx^2*dataset_mean_sd$var)
  dataset_mean_sd$mean =  (dataset_mean_sd$mean-dataset_mean_sd$mean_ref)/dataset_mean_sd$mean_ref + 
    0.5*(dataset_mean_sd$dgdx_square*dataset_mean_sd$var+
           2*0+ # independent variables
           dataset_mean_sd$dgdx_square_ref*dataset_mean_sd$var_ref)
  
  # Fmax = 1 : g = 0
  is_ref <- as.numeric(as.character(dataset_mean_sd$Fmax)) == 1
  dataset_mean_sd$sd[is_ref]   <- 0
  dataset_mean_sd$mean[is_ref] <- 0
  
  return(dataset_mean_sd[, c(colnames(dataset_mean_sd)[c(which(colnames(dataset_mean_sd) %in% colnames(dataset_w_value)))], "sd", "mean", "nbrep")])
}


# Ensure each dataset is in long format once, before entering the nested plotting loops.
as_long_if_needed <- function(df, columns_keep_pivot=colnames(df)[-c(1:4)]) {
  if ("variable" %in% colnames(df)) return(df)
  pivot_longer(df, cols = columns_keep_pivot, names_to = "variable", values_to = "value")
}


# estiamte ssdg as the sum of the delta*Dead/N per age
add_SSDg <- function(df, relative=F, mean_sd=T){   
  df_s <- df %>%
    group_by(variable, Fmax, Fstart, Species, replicate) %>%
    summarize(SSDg = sum(value, na.rm = T))%>%
    ungroup()
  
  if (relative==T) {
    relatived <- df_s %>% filter(Fmax==1) 
    relatived <- relatived[, -which(colnames(relatived) == c("Fmax"))]
    colnames(relatived)[which(colnames(relatived) == "SSDg")] <- "relative"
    dfi <- full_join(df_s, relatived)
    dfi$SSDg <- (dfi$SSDg - dfi$relative)/abs(dfi$relative)
    df_s <- dfi
  }  
  
  if (mean_sd == T){
    df_m <- df_s %>%
      group_by(Fmax, Fstart, variable, Species)%>%
      summarise(across(SSDg, list(sd = sd, mean = mean)))
    
    for(sp in unique(df_m$Species)) {
      Fmax_coll <- Fmax_remove_collapses(sp)
      if (!all(is.na(Fmax_coll))) {
        if (any("Species" %in% colnames(df_m))) {
          df_temp <- df_m[which(df_m$Species == sp), ]
          df_temp[which(df_temp$Fmax %in% Fmax_coll), "SSDg_mean"] <- NA
          df_temp[which(df_temp$Fmax %in% Fmax_coll), "SSDg_sd"] <- NA
          df_m[which(df_m$Species == sp), ] <- df_temp
        } else {
          df_m[which(df_m$Fmax %in% Fmax_coll), "SSDg_mean"] <- NA
          df_m[which(df_m$Fmax %in% Fmax_coll), "SSDg_sd"] <- NA
        }
      }
    } 
  } else {
    for(sp in unique(df_s$Species)) {
      Fmax_coll <- Fmax_remove_collapses(sp)
      if (!all(is.na(Fmax_coll))) {
        if (any("Species" %in% colnames(df_s))) {
          df_temp <- df_s[which(df_s$Species == sp), ]
          df_temp[which(df_temp$Fmax %in% Fmax_coll), "SSDg"] <- NA
          df_s[which(df_s$Species == sp), ] <- df_temp
        } else {
          df_s[which(df_s$Fmax %in% Fmax_coll), "SSDg"] <- NA
        }
      }
    }
    df_m <- df_s
  }
  return(df_m)
}

# shorter and more explicit variable names
improve_sp_names <- function(vector){
  vector <- stringr::str_replace(vector, pattern="EuropeanPlaice", replacement = "European plaice")
  vector <- stringr::str_replace(vector, pattern="AtlanticMackerel", replacement = "Atlantic mackerel")
  vector <- stringr::str_replace(vector, pattern="AtlanticHerring", replacement = "Atlantic herring")
  vector <- stringr::str_replace(vector, pattern="AtlanticCod", replacement = "Atlantic cod")
  vector <- stringr::str_replace(vector, pattern="EuropeanSole", replacement = "European sole")
  vector <- stringr::str_replace(vector, pattern="HorseMackerel", replacement = "Horse mackerel")
  vector <- stringr::str_replace(vector, pattern="NorwayPout", replacement = "Norway pout")
  vector <- stringr::str_replace(vector, pattern="CommonDab", replacement = "Common dab")
  vector <- stringr::str_replace(vector, pattern="CommonSole", replacement = "Common sole")
  vector <- stringr::str_replace(vector, pattern="EuropeanSprat", replacement = "European sprat")
  vector <- stringr::str_replace(vector, pattern="GreyGurnard", replacement = "Grey gurnard")
}

# shorter and more explicit variable names
improve_parm_names <- function(vector){
  vector <- stringr::str_replace(vector, pattern="Fmax", replacement = "\u03B8F")
}

improve_vriable_names <- function(vector, age=F){
  vector <- stringr::str_replace(vector, pattern="mean_of_abundanceMature_per_age", replacement = "Nmat_\u03C4")
  vector <- stringr::str_replace(vector, pattern="Mtot", replacement = "M\u03C4")
  vector <- stringr::str_replace(vector, pattern="abundanceByAge", replacement = "Abundance")
  vector <- stringr::str_replace(vector, pattern="mean_of_meanWeightByAge_per_age", replacement = "Individual weight")
  vector <- stringr::str_replace(vector, pattern="SSDgx_Mpred", replacement = "Predation SSD\u0262 (\u0394²Dead_Pred)")
  vector <- stringr::str_replace(vector, pattern="SSDgx_Mstar", replacement = "Starvation SSD\u0262 (\u0394²Dead_Starv)")
  vector <- stringr::str_replace(vector, pattern="SSDgx_F", replacement = "Fisheries SSD\u0262 (\u0394²Dead_F)")
  vector <- stringr::str_replace(vector, pattern="SSDgx_Mnat", replacement = "Natural SSD\u0262 (\u0394²Dead_Natural)")
  vector <- stringr::str_replace(vector, pattern="Numerator_N_alpha", replacement = "N\u03B1")
  vector <- stringr::str_remove(vector, pattern="\nx")      
  vector <- stringr::str_remove(vector, pattern="x")      
  vector <- stringr::str_replace(vector, pattern="b_rate_sc", replacement = "bx'")  
  vector <- stringr::str_replace(vector, pattern="_rel", replacement = "'")  
  vector <- stringr::str_replace(vector, pattern="Delta", replacement = "\u0394")
  vector <- stringr::str_replace(vector, pattern="square", replacement = "²\u03C4")
  vector <- stringr::str_replace(vector, pattern="\nDiv", replacement = "/")
  vector <- stringr::str_remove(vector, pattern="colSum")      
  vector <- stringr::str_remove(vector, pattern="mean_of_")      
  vector <- stringr::str_remove(vector, pattern="_per_age")      
  vector <- stringr::str_replace(vector, pattern="sum_mortality_rates", replacement = "M")
  vector <- stringr::str_replace(vector, pattern="Mstar", replacement = "Starvation mortality rate")
  vector <- stringr::str_replace(vector, pattern="Mstar", replacement = "Starvation mortality rate")
  vector <- stringr::str_replace(vector, pattern="Madd", replacement = "Madditional")
  vector <- stringr::str_replace(vector, pattern="generation_length", replacement = "T")
  vector <- stringr::str_replace(vector, pattern="all_deaths", replacement = "Deathtot")
  vector <- stringr::str_replace(vector, pattern="b_rate", replacement = "bx")
  vector <- stringr::str_replace(vector, pattern="survivorship_curve", replacement = "L\u03C4")
  vector <- stringr::str_replace(vector, pattern="abundanceByAge", replacement = "N_age")
  vector <- stringr::str_replace(vector, pattern="mean_", replacement = "")
  vector <- stringr::str_replace(vector, pattern="_point", replacement = ".")
  vector <- stringr::str_replace(vector, pattern="Ne_N_alpha", replacement = "N\u2091")
  vector <- stringr::str_replace(vector, pattern="_point", replacement = "\u2022")
  vector <- stringr::str_replace(vector, pattern="F_deaths", replacement = "Fisheries deaths")
  vector <- stringr::str_replace(vector, pattern="Mpred_deaths", replacement = "Predation deaths")
  vector <- stringr::str_replace(vector, pattern="Mnat_deaths", replacement = "Natural deaths")
  vector <- stringr::str_replace(vector, pattern="Mstar_deaths", replacement = "Starvation deaths")
  vector <- stringr::str_replace(vector, pattern="Mstarvation_deaths", replacement = "Starvation deaths")
  vector <- stringr::str_replace(vector, pattern="tot_deaths", replacement = "All deaths")
  vector <- stringr::str_replace(vector, pattern="\nsum", replacement = " - ")
  vector <- stringr::str_replace(vector, pattern="eggNumber_sc", replacement = "Eggs (scaled for constant N\u03B1)")  
  vector <- stringr::str_replace(vector, pattern="F", replacement = "Fishing\nmortality rate\u03C4")
  vector <- stringr::str_replace(vector, pattern="Mpred", replacement = "Predation\nmortality rate\u03C4")
  vector <- stringr::str_replace(vector, pattern="Mnat", replacement = "Natural\nmortality rate\u03C4")
  vector <- stringr::str_replace(vector, pattern="Ne_NT", replacement = "N\u2091/NT")
  vector <- stringr::str_replace(vector, pattern="T/Vk.", replacement = "T/Vk")
  vector <- stringr::str_replace(vector, pattern="T/Vk", replacement = "T/Vk.")
  vector <- stringr::str_replace(vector, pattern="eggNumber", replacement = "egg number_\u03C4")
  vector <- stringr::str_replace(vector, pattern="SSDg/Nα", replacement = "SSD\u0262\u03C4/N\u03B1") # ((\u0394²Dead/N\u03B1)\u03C4)
  vector <- stringr::str_replace(vector, pattern="SSDg", replacement = "SSD\u0262\u03C4") #((\u0394\u03C4²Dead)\u03C4)
  vector <- stringr::str_replace(vector, pattern="T_Vk", replacement = "T/Vk.")
  vector <- stringr::str_replace(vector, pattern="Ne_NT", replacement = "Ne/NT")
  vector <- stringr::str_replace(vector, pattern="survival_rate", replacement = "sx")
  vector <- stringr::str_replace(vector, pattern="Dead", replacement = "D\u03C4")
  vector <- stringr::str_replace(vector, pattern="k_bar", replacement = "k\u03C4")
  # vector <- expression(vector)
  if(age)   vector <- stringr::str_replace(vector, pattern="\u03C4", replacement = "x")
  if(age)   vector <- stringr::str_replace(vector, pattern="\u03C4", replacement = "x")
  if(age)   vector <- stringr::str_replace(vector, pattern="L", replacement = "l")
  return(vector)
}

Fmax_remove_collapses <- function(sp){
  Fmax_coll = NA
  if (sp %in% c("Atlantic cod", "Grey gurnard", "Hake", "AtlanticCod", "GreyGurnard")) Fmax_coll <- c(1.25, 1.5)  
  if (sp %in% c("Norway pout", "Sandeel", "Whiting", "NorwayPout"))                    Fmax_coll <- c(1.5)        
  return(Fmax_coll)
}


# BIPLOTS -----------------------------


# Plot traits vs Fmax with optional relative-to-reference scaling.
# works from raw trait values (value) and does the aggregation itself (mean/sd by Fmax, Fstart, variable). 
# It also auto-computes a secondary-axis scaling and returns both the plot and the scale. 
# It builds a subtitle with per‑group slopes and facet labels for Fstart.
plot_trait_vs_fmax <- function(df_trait_sp, trait_name, Time1, Time2, 
                               relative_to_ref = FALSE, proportional=F,
                               spec, age_dim=FALSE, sp_dim=FALSE, 
                               reference_input=NULL, FORCElm=F, title = "regular", 
                               scale_facet="free",
                               chosen_title=NULL, subtitle=NULL,legend_justif="inside",
                               position_legend="auto", ncol=2, nrow=2, facet=T, scale_z = T, 
                               per_species_scale = F, delta_method=T) {
  
  
  # output name
  plot_name_ch=c()
  plot_name <- paste0(as.character(relative_to_ref), as.character(proportional))
  if (plot_name=="TRUEFALSE") plot_name_ch <- c("Relative ")
  if (plot_name=="FALSEFALSE") plot_name_ch <- c("")
  if (plot_name=="FALSETRUE") plot_name_ch <- c("Proportional ")
  
  # which variables ploted ? 
  var_x <- c("Fmax")
  if (age_dim) var_x <- c("Fmax", "age")
  if (sp_dim) var_x <- c(var_x, "Species")
  PARAM <- c(var_x, ifelse(any(grepl("cohort", colnames(df_trait_sp))), c("cohort"), c(NA)))
  if(any(is.na(PARAM))) PARAM <- PARAM[-which(is.na(PARAM))]
  if("age" %in% colnames(df_trait_sp)) PARAM <- c(PARAM, "age")
  if (length(table(df_trait_sp$Species))>1) PARAM <- c(PARAM, "Species")
  PARAM <- unique(PARAM)
  
  # Filter to target species and cast key numeric columns and rows
  df_trait_sp[, "Fmax"]  <- factor(as.numeric(as.matrix(df_trait_sp[, "Fmax"])), levels = c(0.5, 0.75, 1, 1.25, 1.5))
  df_trait_sp$Species  <- improve_sp_names(df_trait_sp$Species)
  df_trait_sp_F        <- df_trait_sp[which(df_trait_sp$Species %in% spec),]
  if ("Fstart" %in% colnames(df_trait_sp)) df_trait_sp_F$Fstart <- as.numeric(df_trait_sp_F$Fstart)
  if (any(is.infinite(df_trait_sp_F$value))) df_trait_sp_F <- df_trait_sp_F[-which(is.infinite(df_trait_sp_F$value)), ]
  if (age_dim)  {
    df_trait_sp_F$age   <- as.numeric(df_trait_sp_F$age)
    if (!any(grepl(trait_name, pattern="Mstar"), grepl(trait_name, pattern="Vsquare"),  grepl(trait_name, pattern="VsquarexDead")) & any(which(df_trait_sp_F$value == 0))) {
      df_trait_sp_F <- df_trait_sp_F[-which(df_trait_sp_F$value == 0), ]}
    if (any(is.na(df_trait_sp_F$value))) df_trait_sp_F <- df_trait_sp_F[-which(is.na(df_trait_sp_F$value)), ]
  }
  
  #---------------------------------------------------------------------------#
  #  PANELS BY SPECIES ----------
  # sec_axis() applies a single transformation to all facets: with
  # facet_wrap, the scale factor was calculated across all species combined
  # Solution: the function calls itself for each species,
  # so a scale factor is calculated for that species alone, and the panels are then
  # assembled using patchwork. A species thus has the same axes in all the
  # figures.
  # ---------------------------------------------------------------------------#
  nvar <- length(table(df_trait_sp_F$variable))
  if (facet & per_species_scale &  #!age_dim &&
      length(unique(df_trait_sp_F$Species)) > 1) {
    data_sp <- list() # store relative values non rescaled per sp
    panels<-lapply(spec, function(s){
      cat(s)
      #  # Assemble the plot with loess trends and a secondary axis.
      fill_vector=c("palegreen3", "#4e4b5c", "brown2", "cornflowerblue", "lightyellow", "grey20", "orange", "yellow2", "yellow4", "grey30")
      shape_vector=c(15:20)  

      res <- plot_trait_vs_fmax(df_trait_sp = df_trait_sp_F[df_trait_sp_F$Species == s, ],
                                trait_name = trait_name, Time1 = Time1, Time2 = Time2,
                                relative_to_ref = relative_to_ref, proportional = proportional,
                                spec = s, age_dim = age_dim, sp_dim = sp_dim,
                                reference_input = reference_input, FORCElm = FORCElm,
                                scale_facet = scale_facet, subtitle = FALSE,
                                position_legend = position_legend, facet = ifelse(nvar>1 & age_dim, T, F),
                                scale_z = scale_z, per_species_scale = FALSE, 
                                delta_method=delta_method, 
                                ncol=ifelse(nvar>1, nvar, ncol), nrow=ifelse(nvar>1 & age_dim, 1, nrow))
      if (is.null(res$plot)) return(NULL)
      if (nvar>1 & age_dim){
        if (nlevels(df_trait_sp_F$Fmax)>1){
          all_var = c("0.5", "0.75", "1", "1.25", "1.5")
        } else {
          all_var  <- if (is.factor(df_trait_sp_F$variable)) levels(droplevels(df_trait_sp_F$variable)) else unique(as.character(df_trait_sp_F$variable))
        }
        cols_all <- setNames(fill_vector[seq_along(all_var)], all_var)
        shp_all  <- setNames(shape_vector[seq_along(all_var)], all_var)
        p_base <- if (inherits(res$plot, "ggplot")) res$plot else res$plot[[1]]   
        p <- p_base +
          facet_wrap(~Species*variable, scale="free") +                # nom de l'espèce en bandeau
          facet_wrap(~Species*variable, scale="free") +                # nom de l'espèce en bandeau
          labs(title = NULL, subtitle = NULL) +
          theme(strip.background = element_blank(),
                strip.text = element_text(size = 10),
                axis.title = element_blank())+   # titres d'axes communs, ajoutés plus bas
          # scale_colour_manual(name = "", values = cols_all, limits = all_var, drop = FALSE) +
          # scale_shape_manual(name = "", values = shp_all, limits = all_var, drop = FALSE)
          scale_colour_manual(name = NULL,
                              values = colorRampPalette(c("palegreen3", "brown2"))(length(levels(df_trait_sp_F$Fmax))), 
                              limits = levels(df_trait_sp_F$Fmax)) +
          scale_shape_manual(name = NULL, values = shp_all, limits = levels(df_trait_sp_F$Fmax), drop = FALSE)
      }      
      if (nvar==1 || !age_dim){
        if (nvar == 1){
          all_var = c("0.5", "0.75", "1", "1.25", "1.5")
          lvar <- length(all_var)
          cols_all <- setNames(fill_vector[seq_along(all_var)], all_var)
          shp_all  <- setNames(shape_vector[seq_along(all_var)], all_var)
          p <- res$plot+
            scale_colour_manual(name = NULL,
                                values = colorRampPalette(c("palegreen3", "brown2"))(lvar), 
                                limits = as.numeric(all_var)) +
            scale_shape_manual(name = NULL, 
                               values = shp_all, 
                               limits = as.numeric(all_var), 
                               drop = FALSE)
        } else {
          all_var = trait_name
          cols_all <- setNames(fill_vector[seq_along(all_var)], all_var)
          shp_all  <- setNames(shape_vector[seq_along(all_var)], all_var)
          p <- res$plot+
            scale_colour_manual(name = NULL,
                                values = cols_all) +
            scale_shape_manual(name = NULL, values = shp_all, 
                               drop = FALSE)
        }
        p <- p +
          facet_wrap(~Species, scale="free") +                # nom de l'espèce en bandeau
          labs(title = NULL, subtitle = NULL) +
          theme(strip.background = element_blank(),
                strip.text = element_text(size = 10),
                axis.title = element_blank())#+
        # geom_point(show.legend=T)+
        # geom_line(show.legend=T)
      }
      p
    })
    data_sp <- lapply(spec, function(s) {
      cat(s)
      # browser()
      res <- plot_trait_vs_fmax(df_trait_sp = df_trait_sp_F[df_trait_sp_F$Species == s, ],
                                trait_name = trait_name, Time1 = Time1, Time2 = Time2,
                                relative_to_ref = relative_to_ref, proportional = proportional,
                                spec = s, age_dim = age_dim, sp_dim = sp_dim,
                                reference_input = reference_input, FORCElm = FORCElm,
                                scale_facet = scale_facet, subtitle = FALSE,
                                position_legend = position_legend, facet = ifelse(nvar>1, T, F),
                                scale_z = scale_z, per_species_scale = FALSE, 
                                delta_method=delta_method, 
                                ncol=ifelse(nvar>1, nvar, ncol), nrow=ifelse(nvar>1, 1, nrow))
      if (is.null(res$plot)) return(NULL)
      data_sp[[s]] <- res$data
      return(res$data)
    })
    
    if (age_dim) ncol = floor(ncol/nvar)
    names(panels) <- spec
    panels <- panels[!sapply(panels, is.null)]
    # titres d'axes communs (on les récupère sur un panneau)
    lab_left  <- panels[[1]]$labels$y
    lab_right <- tryCatch(panels[[1]]$scales$get_scales("y")$secondary.axis$name,
                          error = function(e) "")
    if (is.null(lab_right) || inherits(lab_right, "waiver")) lab_right <- ""
    txt <- function(label, angle) ggplot() +
      annotate("text", x = 0, y = 0, label = label, angle = angle, size = 4.2) + theme_void()
    x_lab <- ggplot() + annotate("text", x = 0, y = 0, label = ifelse(age_dim == T, "age", "theta[F]"), parse = TRUE,
                                 size = 4.2) + theme_void()
    per_page <- ncol * nrow
    pages <- split(seq_along(panels), ceiling(seq_along(panels) / per_page))
    p_list <- lapply(pages, function(idx) {
      grid <- wrap_plots(panels[idx], ncol = ncol, nrow = nrow)
      (((txt(lab_left, 90) | grid | txt(lab_right, -90)) +
          plot_layout(widths = c(0.035, 1, 0.035))) /
          x_lab + plot_layout(heights = c(1, 0.025))) +
        plot_layout(guides = "collect") &
        theme(legend.position = "top", legend.direction = "horizontal",
              legend.title = element_blank(),
              legend.box.background = element_rect(colour = "black"))
    })
    return(list(plot = unname(p_list), scale = NA,
                data = dplyr::bind_rows(data_sp))) 
  }
  
  #---------------------------------------------------------------------------#
  #  PANELS in FACETS -----------
  # all dataset treated at once
  #---------------------------------------------------------------------------#
  
  
  ## Variables transformations :  proportional or relative ? ------
  # proportional: 
  # => case 1: X/Y, the value proportional to the reference (Fmax = 1)
  # => case 2: log(X/Y), the logarithm of the value proportional to the reference (Fmax = 1)
  # Replicates are independent across fishing scenarios: no pairing.
  # We compare the MEANS per scenario (mean, SE = sd/sqrt(n)) with those for Fmax = 1,
  # and the SEs are combined as for two independent samples.
  if (proportional == T){
    group_cols <- unique(c(var_x, "Species", "Fstart", "variable"))
    df_mean_se <- df_trait_sp_F %>%
      group_by(across(all_of(group_cols))) %>%
      summarise(m  = mean(value, na.rm = TRUE),
                nbrep = length(value),
                se = sd(value, na.rm = TRUE)/sqrt(length(value)),
                .groups = "drop")
    # reference: same species, Fstart, variable (and age if age_dim), at Fmax = 1
    ref_data_percent <- df_mean_se %>%
      filter(as.numeric(as.character(Fmax)) == 1) %>%
      dplyr::select(-Fmax, -nbrep) %>%        
      dplyr::rename(m_ref = m, se_ref = se)
    if (delta_method) {
      cat("delta method")
      if (any(grepl("^log", df_mean_se$variable))){
        df_trait_sp_F_mean_sd <- df_mean_se %>%
          left_join(ref_data_percent, by = setdiff(group_cols, "Fmax"), relationship = "many-to-one") %>%
          mutate(is_log = grepl("^log", variable),
                 mean = ifelse(is_log, m - m_ref,
                               log(m / m_ref) + 0.5*(-se^2/m^2 + se_ref^2/m_ref^2)),
                 sd   = ifelse(is_log, sqrt(se^2 + se_ref^2),
                               sqrt(se^2/m^2 + se_ref^2/m_ref^2))) %>%
          dplyr::select(all_of(group_cols), mean, sd, nbrep)   
      } else{
        df_trait_sp_F_mean_sd <- df_mean_se %>%
          left_join(ref_data_percent, by = setdiff(group_cols, "Fmax"), relationship = "many-to-one") %>%
          mutate(
            mean = (m / m_ref) * (1 + (se_ref^2 / m_ref^2)),
            sd   = sqrt(se^2 / m_ref^2 + (m^2 * se_ref^2) / m_ref^4) 
          ) %>%
          dplyr::select(all_of(group_cols), mean, sd, nbrep)   
      }
    } else {
      cat("no delta method")
      ref_data_percent <- df_trait_sp_F %>%
        filter(as.numeric(as.character(Fmax)) == 1) %>%
        dplyr::select(-Fmax) %>%
        dplyr::rename(ref = value) %>%
        ungroup()
      df_trait_sp_F_mean_sd <- df_trait_sp_F %>%
        left_join(ref_data_percent, by = c(setdiff(group_cols, "Fmax"), "replicate", "min_age"), relationship = "many-to-one")%>%
        mutate(proportional = value/ref) %>%
        group_by(across(all_of(group_cols))) %>%
        summarise(mean  = mean(proportional, na.rm = TRUE),
                  sd    = sd(proportional, na.rm = TRUE),
                  nbrep = length(value),
                  .groups = "drop")
    }
    
    is_ref <- as.numeric(as.character(df_trait_sp_F_mean_sd$Fmax)) == 1
    df_trait_sp_F_mean_sd$mean[is_ref] <- ifelse(grepl("^log", df_trait_sp_F_mean_sd$variable[is_ref]), 0, 1)
    df_trait_sp_F_mean_sd$sd[is_ref]   <- 0
    if (any(is.nan(df_trait_sp_F_mean_sd$mean)))
      warning("proportional : moyenne <= 0 dans une variable non loguée -> NaN")
  }
  
  if (relative_to_ref){
    df_trait_sp_F_mean_sd <- value_to_sd_esperance_delta_method_relative(dataset_w_value = df_trait_sp_F, var_x = var_x)
  }
  
  if (!proportional & !relative_to_ref){
    if (length(var_x)==1){ 
      df_trait_sp_F_mean_sd <- df_trait_sp_F %>%
        group_by(get(var_x), Species, Fstart, variable) %>%
        summarise(mean = mean(value, na.rm = TRUE),
                  nbrep = length(value),
                  sd  = sd(value, na.rm = TRUE)) %>%
        ungroup()
    }
    if (length(var_x) == 2){
      df_trait_sp_F_mean_sd <- df_trait_sp_F %>%
        group_by(get(var_x[1]),get(var_x[2]), Species, Fstart, variable) %>%
        summarise(mean = mean(value, na.rm = TRUE),
                  nbrep = length(value),
                  sd   = sd(value, na.rm = TRUE)) %>%
        ungroup()
    }
  }
  
  # update colnames
  if (any(colnames(df_trait_sp_F_mean_sd) == "get(var_x)")) colnames(df_trait_sp_F_mean_sd)[which(colnames(df_trait_sp_F_mean_sd) == "get(var_x)")] <- var_x
  if (any(colnames(df_trait_sp_F_mean_sd) == "get(var_x[1])")) colnames(df_trait_sp_F_mean_sd)[which(colnames(df_trait_sp_F_mean_sd) == "get(var_x[1])")] <- var_x[1]
  if (any(colnames(df_trait_sp_F_mean_sd) == "get(var_x[2])")) colnames(df_trait_sp_F_mean_sd)[which(colnames(df_trait_sp_F_mean_sd) == "get(var_x[2])")] <- var_x[2]
  
  # standard error is better than standard deviation for the IC
  # L'écart type nous permet de représenter la stochasticité intrinsèque du modèle, 
  # ce qui en soi est une info importante. Il n'indique cependant pas la précision de la moyenne, 
  # qui est donnée plutôt par l'erreur standard (SE = SD/sqrt(n) avec n le nombre de points, soit le nombre de répliquas pour nous), 
  # voire la confiance qu'on a dans cette moyenne (IC95 = t_(19,0.975) SE = 2.093 SD/sqrt(20) = 0.468 SD). 
  # L'utilisation de l'IC95 serait plus rigoureuse statistiquement pour la comparaison à la valeur 0 que 
  # l'écart type et en plus moins conservative (puisqu'environ moitié de la SD) de sorte qu'on récupèrerait vraisemblablement des espèces supplémentaires avec une réponse (herring, common dab, sandeel au moins). 
  sd_is_se <- relative_to_ref | (proportional & delta_method)                      # [CORRIGÉ] sd is already a SE -> no further division by sqrt(n)
  if (!any(grepl("se", colnames(df_trait_sp_F_mean_sd)))){
    grouping = c(var_x, "Species", "Fstart", "variable")
    df_trait_sp_F_mean_sd <- df_trait_sp_F_mean_sd %>%
      group_by(across(all_of(grouping))) %>%
      mutate(IC_95 = qt(0.975, df = nbrep - 1) * sd / ifelse(sd_is_se, 1, sqrt(nbrep))) %>%  
      ungroup()
  }
  
 # Prepare data for aesthetic plot --------
  # Normalize types for plotting and factor ordering.
  df_trait_sp_F_mean_sd$Fstart      = as.factor(df_trait_sp_F_mean_sd$Fstart)
  df_trait_sp_F_mean_sd[,var_x[1]]  = as.numeric(as.matrix(df_trait_sp_F_mean_sd[,var_x[1]]))
  if (length(var_x)>=2) df_trait_sp_F_mean_sd[,var_x[2]]  = as.numeric(as.matrix(df_trait_sp_F_mean_sd[,var_x[2]]))
  
  # remove the Fmax for which species abundance collapsed
  for(sp in spec) {                                                     
    Fmax_coll <- Fmax_remove_collapses(sp)
    if (!all(is.na(Fmax_coll))) {
      id_coll <- which(as.numeric(as.character(df_trait_sp_F_mean_sd$Fmax)) %in% Fmax_coll)   
      if ("Species" %in% colnames(df_trait_sp_F_mean_sd)) {
        id_coll <- intersect(id_coll, which(df_trait_sp_F_mean_sd$Species == sp))
      }
      df_trait_sp_F_mean_sd[id_coll, "mean"]  <- NA
      df_trait_sp_F_mean_sd[id_coll, "IC_95"] <- NA
      if ("sd" %in% colnames(df_trait_sp_F_mean_sd)) df_trait_sp_F_mean_sd[id_coll, "sd"] <- NA   
    }
  }
  
  # data for output
  data_unscaled <- df_trait_sp_F_mean_sd
  
  # what variable have scale y axis and y and z axis?
  scale = 1
  scale_factor <- 1
  shift <- 0
  if (nrow(df_trait_sp_F_mean_sd)==0) return(list(c(), c()))
  mean_var <- aggregate(df_trait_sp_F_mean_sd[, "mean"], list(df_trait_sp_F_mean_sd$variable), mean, na.rm=T)
  max_var  <- aggregate(df_trait_sp_F_mean_sd[, "mean"], list(df_trait_sp_F_mean_sd$variable), max, na.rm=T)
  
  # adjust y and z axis labels
  var_scale_regular <- NULL
  var_scale_regular = mean_var[which((abs(mean_var[, 2]-max(mean_var[, 2])) > (0.8 * (max(mean_var[, 2])))) & !(mean_var[,2] == max(mean_var[, 2]))), 1]
  if (all(trait_name == c("T", "Vk.", "T/Vk."))) {
    var_scale_regular = ifelse(any(var_scale_regular %in% c("T", "Vk.")), 1, 2)
    var_scale_regular = c("T", "Vk.")
  }
  if (all(trait_name == c("Nₑ", "NT", "Nₑ/NT")))  {
    var_scale_regular = ifelse(any(var_scale_regular %in% c("Nₑ", "NT")), 1, 2)
    var_scale_regular = c("Nₑ", "NT")
  }  
  if (all(trait_name == c( "Nₑ","NT","Nₑ/NT", "Nα")))  {
    var_scale_regular = ifelse(any(var_scale_regular %in% c("Nₑ", "NT", "Nα")), 1, 2)
    var_scale_regular = c("Nₑ", "NT", "Nα")
  }  
  if (all(trait_name == c("SSDɢ,τ", "Nα", "(SSDɢ/Nα)τ")))  {
    var_scale_regular = ifelse(any(var_scale_regular %in% c("SSDɢ,τ", "Nα", "(SSDɢ/Nα)τ")), 1, 2)
    var_scale_regular = c("SSDɢ,τ", "Nα")
  }  
  if (all(trait_name == c("Nα", "Δ²τ", "SSDɢ,τ")))  {
    var_scale_regular = c("Nα", "Δ²τ")
  }  
  if (all(trait_name == c("log(T/(Vk.+2))", "log(Nα/NT)", "log(Nₑ/NT)")))  {
    var_scale_regular = c("log(T/(Vk.+2))", "log(Nα/NT)")
  }
  if (all(trait_name == c("Deadx", "Δ²x", "Mortality ratex")))  {
    var_scale_regular = c("Deadx", "Δ²x", "Mortality ratex")
  }
  
 ## Rescale variables depending if y or z axis variables ------
  # Default functions
  scale_function <- function(x, scale, shift){ 
    return((x)*1)
  }
  inv_scale_function <- scale_function
  # Process
  if (length(var_scale_regular) == length(trait_name)) var_scale_regular <- NULL
  id_keep    <- seq_len(nrow(df_trait_sp_F_mean_sd))
  max_first  <- max(df_trait_sp_F_mean_sd$mean[id_keep], na.rm = TRUE)   
  min_first  <- min(df_trait_sp_F_mean_sd$mean[id_keep], na.rm = TRUE)   
  id_keep    <- seq_len(nrow(df_trait_sp_F_mean_sd))[which(df_trait_sp_F_mean_sd$variable %in% var_scale_regular)]
  do_scale   <- length(var_scale_regular) > 0 & length(id_keep) > 0 &
              ((length(df_trait_sp_F_mean_sd$mean[-id_keep]))!=0 ) & scale_z == T
  if (do_scale) { #https://finchstudio.io/blog/ggplot-dual-y-axes/
    # Function to scale secondary axis
    scale_function <- function(x, scale, shift){ return ((x-shift)/scale) }
    # Function to scale secondary variable values
    inv_scale_function <- function(x, scale, shift){ return ((x)*scale + shift) }
    max_first  <- max(df_trait_sp_F_mean_sd$mean[ id_keep], na.rm = TRUE) 
    min_first  <- min(df_trait_sp_F_mean_sd$mean[ id_keep], na.rm = TRUE) 
    min_second <- min(df_trait_sp_F_mean_sd$mean[-id_keep], na.rm = TRUE) 
    max_second <- max(df_trait_sp_F_mean_sd$mean[-id_keep], na.rm = TRUE) 
    if(max_second==0) max_second <- max_first
    shift=0
    scale_factor <- diff(range(df_trait_sp_F_mean_sd$mean[id_keep], na.rm = T)) / diff(range(df_trait_sp_F_mean_sd$mean[-id_keep], na.rm = T))
    if (is.nan(scale_factor)) scale_factor <- diff(range(df_trait_sp_F_mean_sd$mean[id_keep][-which(is.nan(df_trait_sp_F_mean_sd$mean[id_keep]))], na.rm = T)) /
      diff(range(df_trait_sp_F_mean_sd$mean[-id_keep][-which(is.nan(df_trait_sp_F_mean_sd$mean[-id_keep]))], na.rm = T))
    if (is.infinite(scale_factor) || is.nan(scale_factor) || scale_factor==0) scale_factor <- 1
    df_trait_sp_F_mean_sd$mean[-id_keep] <-
      inv_scale_function(df_trait_sp_F_mean_sd$mean[-id_keep], scale_factor, shift)
    df_trait_sp_F_mean_sd$IC_95[-id_keep] <-
      inv_scale_function(df_trait_sp_F_mean_sd$IC_95[-id_keep], scale_factor, shift)
  }
  if (!(all(seq_len(nrow(df_trait_sp_F_mean_sd)) %in% id_keep))) id_keep <- seq_len(nrow(df_trait_sp_F_mean_sd))[-id_keep]
  max_second <- max(c(df_trait_sp_F_mean_sd$mean[id_keep]), na.rm = TRUE) 
  min_second <- min(c(df_trait_sp_F_mean_sd$mean[id_keep]), na.rm = TRUE) 
  if (any(is.infinite(max_first)))  max_first  <- max(df_trait_sp_F_mean_sd$mean[-id_keep][-which(is.infinite(df_trait_sp_F_mean_sd$mean[-id_keep]))], na.rm = TRUE)   
  if (any(is.infinite(max_second))) max_second <- max(df_trait_sp_F_mean_sd$mean[ id_keep][-which(is.infinite(df_trait_sp_F_mean_sd$mean[id_keep]))], na.rm = TRUE)   
  if (any(is.infinite(min_first)))  min_first  <- min(df_trait_sp_F_mean_sd$mean[-id_keep][-which(is.infinite(df_trait_sp_F_mean_sd$mean[-id_keep]))], na.rm = TRUE)   
  if (any(is.infinite(min_second))) min_second <- min(df_trait_sp_F_mean_sd$mean[ id_keep][-which(is.infinite(df_trait_sp_F_mean_sd$mean[id_keep]))], na.rm = TRUE)   
  
  # Decide which variable name goes on the secondary axis (if there is one = TRUE)
  supp_name = ""
  if (length(unique(df_trait_sp_F$Fmax))==1) supp_name = paste0("Fmax = ", unique(df_trait_sp_F$Fmax))  
  if (length(unique(df_trait_sp_F$age))==1)  supp_name = paste0("age = ", unique(df_trait_sp_F$age))
  z_axis <- do_scale
  if(z_axis)  {z_axis <- paste0(plot_name_ch, "change in ", paste(trait_name[-which(trait_name %in% var_scale_regular)], collapse = ", "))
  } else       z_axis <- ""
  if (nchar(as.character(z_axis[1]))>1) {
    name_y <- paste0(plot_name_ch, "change in ", paste(var_scale_regular, collapse=", "))
  }else  name_y <-  paste0(plot_name_ch, "change in ", paste(mean_var[,1] , collapse=", "))
  limits = c(min(c(df_trait_sp_F_mean_sd$mean), na.rm = T), 
             max(c(df_trait_sp_F_mean_sd$mean), na.rm = T))
  
 ## Build slope summary text for the subtitle ---------
  title = paste0(plot_name_ch, " ", paste(trait_name, collapse=", "), " variation over y 150 to 160")
  if (!is.null(chosen_title)) title = paste0(plot_name_ch, " ", chosen_title)
  method = ifelse(max(nchar(trait_name[1]))> 30 || FORCElm==T, "lm", "loess")
  labeller <- as_labeller(c(
    `0.5` = "Fmax 0.5",
    `0.75` = "Fmax 0.75",
    `1` = "Fmax 1",
    `1.25` = "Fmax 1.25",
    `1.5` = "Fmax 1.5"
  ))
  # Assemble the plot with loess trends and a secondary axis.
  df_trait_sp_F_mean_sd$Fmax = as.factor(df_trait_sp_F_mean_sd$Fmax)
  if ("age" %in% colnames(df_trait_sp_F_mean_sd)) df_trait_sp_F_mean_sd$age = as.numeric(df_trait_sp_F_mean_sd$age)
  df_trait_sp_F_mean_sd$variable = as.factor(df_trait_sp_F_mean_sd$variable)

 # Assemble the plot -----------
  df_trait_sp_F_mean_sd$variable <- as.factor(df_trait_sp_F_mean_sd$variable)   
  df_trait_sp_F_mean_sd$Species  <- factor(df_trait_sp_F_mean_sd$Species, levels=spec)
  lev_var  <- levels(df_trait_sp_F_mean_sd$variable)
  lev_Fmax <- levels(as.factor(df_trait_sp_F_mean_sd$Fmax))                     
  
  # prepare colors and shapes and fills
  fill_vector  <- c("palegreen3", "#4e4b5c", "brown2", "cornflowerblue", "lightyellow", "grey20", "orange", "yellow2", "yellow4", "grey30")
  shape_vector <- c(15:20)
  line_vector  <- c(1:9)
  col_vector   <- rep(fill_vector, 2)
  if (length(lev_var) > 6) {                                                     
    shape_vector <- rep(15, length(fill_vector))
    fill_vector  <- rep("white", length(fill_vector))
    line_vector  <- rep(1, length(fill_vector))
  }
  # named vectors: one value per level
  fill_vectorfmax <- setNames(fill_vector[seq_along(lev_Fmax)], lev_Fmax)
  fill_vector     <- setNames(fill_vector[seq_along(lev_var)],  lev_var)
  shape_vector    <- setNames(shape_vector[seq_along(lev_var)], lev_var)
  line_vector     <- setNames(line_vector[seq_along(lev_var)],  lev_var)
  # one value per row (vectorised lookup instead of tempo/tempo2 + loops)      
  v <- as.character(df_trait_sp_F_mean_sd$variable)
  fill_vectori     <- unname(fill_vector[v])
  fill_vectorifmax <- unname(fill_vectorfmax[as.character(df_trait_sp_F_mean_sd$Fmax)])
  shape_vectori    <- as.character(shape_vector[v])
  line_vectori     <- as.character(line_vector[v])
  # shapes on "variable" if both vary, otherwise on whichever has the most levels
  n_var  <- length(table(df_trait_sp_F_mean_sd$variable))
  n_Fmax <- length(table(df_trait_sp_F_mean_sd$Fmax))
  shape_var <- if ((n_var > 1 & n_Fmax > 1) | n_var >= n_Fmax) "variable" else "Fmax"   
  df_trait_sp_F_mean_sd$variable <- factor(df_trait_sp_F_mean_sd$variable)

  if (length(var_x) == 1) { #only one dimension. If facet = by variable
    
    # remove the variables that have only one age class or NAs in the ssdg sums for instance
    if (any(sapply(df_trait_sp_F_mean_sd$variable, FUN = function(x) (all(is.na(df_trait_sp_F_mean_sd$mean[which(df_trait_sp_F_mean_sd$variable == x)])))))) {
      v_null <- df_trait_sp_F_mean_sd$variable[which(sapply(df_trait_sp_F_mean_sd$variable, FUN = function(x) {
        cat(x)
        (all(is.na(df_trait_sp_F_mean_sd$mean[which(df_trait_sp_F_mean_sd$variable == x)])))
      }))]
      add_null <- which(is.na(df_trait_sp_F_mean_sd$Fmax))
      df_trait_sp_F_mean_sd <- df_trait_sp_F_mean_sd[-c(which(df_trait_sp_F_mean_sd$variable %in% v_null), add_null),]
    }

    p <- ggplot(df_trait_sp_F_mean_sd, aes(x = .data[[var_x]], y = mean,
                group = interaction(variable, Fmax))) +
      geom_errorbar(aes(ymin=mean-IC_95, ymax=mean+IC_95, group = .data[[shape_var]], color= .data[[shape_var]]
                        # alpha=as.numeric(.data[[shape_var]])/10,
                        ), alpha=0.5, #as.numeric(line_vectori)*0.3,
                    show.legend = FALSE, position = position_dodge(0.3), 
                    width = 0.2)+
      scale_color_manual(name=NULL, values=as.character(fill_vector))+
      geom_point(data=df_trait_sp_F_mean_sd, aes(shape = .data[[shape_var]], color = .data[[shape_var]]), 
                size=2.1, position = position_dodge(0.3),  
                stroke=1.2, show.legend = T)+ #, , stroke=1.5
      scale_shape_manual(name=NULL, values=as.numeric(shape_vectori)) +
      # scale_fill_manual(name="",values=as.character(fill_vectori)) +
    # geom_smooth(method = method, se=TRUE, col = "black",
    #   aes(group = interaction(variable, Fstart))
    # ) + 
    # geom_point(size = 1.5, alpha=0.5) +
    theme(panel.background = element_rect(fill = "white", colour = "grey50"), 
          strip.background = element_rect(fill="white"),
          strip.text = element_text(size = 9.5),
          plot.title = element_text(hjust = 0.5), 
          axis.text.x = element_text(angle = 60, hjust = 1),
          axis.text = element_text(colour = 1, size = 10),
          panel.border = element_rect(colour = "black", fill=NA),
          legend.background = element_blank(),
          legend.justification = c("top"),
          legend.title = element_blank(), 
          legend.text=element_text(size=9), 
          legend.box.background = element_rect(colour = "black"), 
          legend.key = element_rect(colour = NA))+
    # facet_wrap(~Fstart, labeller = Fstart_names) +
    labs(
      # title = title,
      x = improve_parm_names(c(var_x)),
      y =  name_y, 
    ) 
    
    if (!(z_axis == "")){
      p <- p+
        scale_y_continuous(
          # limits = limits,
          sec.axis = sec_axis(~scale_function(., scale_factor, shift), #. * scale, # (. - a)/b,  #sec_axis(~. * scale, 
                              name = z_axis
          )
        )
    }
    p <- p + scale_x_discrete(limits = c("0.5", "0.75", "1", "1.25", "1.5"), drop = FALSE)
    
    if (length(table(df_trait_sp_F_mean_sd$Species))>1 & facet){
      p_list=list()
      for (i in 1:ceiling(c(length(spec)/(ncol*nrow)))){
        pi <- p+
          facet_wrap_paginate(~Species,nrow = nrow, ncol=ncol,drop=FALSE, scales = scale_facet, page=i)+
          theme(legend.position="top")
        p_list[[i]] <- pi
      }
      p <- p_list
      if (legend_justif =="inside") p <- p+ theme(legend.box = "horizontal")
    }
    if (length(table(df_trait_sp_F_mean_sd$Species))==1 &
        length(table(df_trait_sp_F_mean_sd$variable))>1 & facet==T){
      p <- p+
          facet_wrap_paginate(~variable,nrow = nrow, ncol=ncol,drop=FALSE, scales = scale_facet)+
          theme(legend.position.inside = c(1,0.3), 
                legend.direction = "vertical")
      if (legend_justif =="inside") p <- p+ theme(legend.box = "horizontal")
      
    }
    
    if (relative_to_ref | proportional){
      if (is.list(p) && !ggplot2::is.ggplot(p)){
        for (i in 1:ceiling(c(length(spec)/(ncol*nrow)))){
          p[[i]] <- p[[i]] + 
            theme(legend.spacing.y = unit(0, "mm")) +
            geom_hline(col="grey80", linetype="dotted", lwd=0.8, alpha=0.6, yintercept=0)+
            geom_hline(col="grey80", linetype="dotted", lwd=0.8, alpha=0.6, 
                       yintercept=inv_scale_function(0, scale = scale, shift = shift))
        }
      } else p <- p + 
        theme(legend.position = "inside",
              legend.position.inside = c(1,1),
              legend.spacing.y = unit(0, "mm")) +
        geom_hline(col="grey80", linetype="dotted", lwd=0.8, alpha=0.6, yintercept=0)+
        geom_hline(col="grey80", linetype="dotted", lwd=0.8, alpha=0.6, yintercept=inv_scale_function(0, scale = scale, shift = shift))
    }
    if (is.list(p) && !ggplot2::is.ggplot(p)) {
      for (i in 1:ceiling(c(length(spec)/(ncol*nrow)))){
        p[[i]] <- p[[i]] +
          new_scale_color()+
          geom_line(aes(group = .data[[shape_var]], 
                        col=as.factor(as.numeric(.data[[shape_var]])/10)), 
                    show.legend = T, size=0.9, position = position_dodge(0.3))+
          scale_color_manual(name="",values=as.character(fill_vector[-which(is.na(names(fill_vector)))]))
      }
    } else p <- p +
      geom_line(aes(group = .data[[shape_var]], color = .data[[shape_var]]),
                show.legend = FALSE, linewidth = 0.9, position = position_dodge(0.3))
  } else {
    df_trait_sp_F_mean_sd$variable <- factor(df_trait_sp_F_mean_sd$variable)
    df_trait_sp_F$Species <- as.factor(df_trait_sp_F$Species)
    p <- ggplot(df_trait_sp_F_mean_sd, aes(x = .data[[var_x[2]]], y = mean, #mean,
                                   shape = .data[[var_x[1]]], colour= .data[[var_x[1]]], group = interaction(variable, .data[[var_x[1]]]))) +
      geom_point(aes(y=mean, x=.data[[var_x[2]]]), size=2, show.legend = T)+
      geom_line(linewidth = 1, show.legend = T)+
      scale_colour_manual(name = NULL,
                          values = colorRampPalette(c("palegreen3", "brown2"))(length(levels(df_trait_sp_F_mean_sd[[var_x[1]]]))), 
                          limits = levels(df_trait_sp_F_mean_sd[[var_x[1]]])) +
      theme_bw() +
      labs(shape = c(""), colour=improve_parm_names("Fmax"))+
      labs(
        x = improve_parm_names(c(var_x[[2]])),
        y =  name_y, 
      ) +
      theme(panel.background = element_rect(fill = "white", colour = "grey50"), 
            strip.background = element_rect(fill="white"),
            strip.text = element_text(size = 9.5),
            plot.title = element_text(hjust = 0.5), 
            axis.text.x = element_text(angle = 60, hjust = 1),
            axis.text = element_text(colour = 1, size = 10),
            panel.border = element_rect(colour = "black", fill=NA),
            legend.background = element_blank(),
            legend.justification = c("top"),
            legend.text=element_text(size=9), 
            legend.box.background = element_rect(colour = "black"), 
            legend.key = element_rect(colour = NA))
    
    if (!(z_axis == "")){
      p <- p+
        scale_y_continuous(
          sec.axis = sec_axis(~scale_function(., scale_factor, shift), #. * scale, # (. - a)/b,  #sec_axis(~. * scale, 
                              name = z_axis
          )
        )
    }
    
    if (proportional) {
      p <- p +       geom_hline(col="grey", linetype="dotted", lwd=0.8, alpha=0.4, yintercept=1)
      data_tempo <- na.omit(df_trait_sp_F_mean_sd$mean)
      if (any(is.infinite(data_tempo))) data_tempo <- data_tempo[-which(is.infinite(data_tempo))]
      if (any(data_tempo>10)) data_tempo <- data_tempo[-which(data_tempo>10)]
    }
    
    if (relative_to_ref) p <- p +    geom_hline(col="grey", linetype="dotted", lwd=0.8, alpha=0.4, yintercept=0)
    
    if (facet & legend_justif=="inside") {
      p <- p+
        theme(legend.position = legend_justif,
              legend.spacing.y = unit(0, "mm")) 
    }
    class_p <- class(p)
    if (length(table(df_trait_sp_F_mean_sd$Species))>1 & facet==T & !(c("list") %in% class_p) ){
      p_list=list()
      for (i in 1:ceiling(c((length(t)*length(spec))/(ncol*nrow)))){
        cat(i)
        p_list[[i]] <- p+
          facet_wrap_paginate(~Species*variable,nrow = nrow, ncol=ncol,drop=FALSE, scales = scale_facet, page=i)+
          theme(legend.position = "top", 
                legend.direction = "horizontal")
      }
      p <- p_list
    }
    
    if (length(table(df_trait_sp_F_mean_sd$Species))==1 &
      length(table(df_trait_sp_F_mean_sd$variable))>1 & facet==T & (c("list") %in% class_p)){
      p_list=list()
      for (i in 1:ceiling(c((length(spec))/(ncol*nrow)))){
        p_list[[i]] <-  p[[i]]+
          facet_wrap_paginate(~variable,nrow = nrow, ncol=ncol,drop=FALSE, scales = scale_facet, page=i)+
          theme(legend.position.inside = c(1,0.3), 
                legend.direction = "vertical", )
      }
      p <- p_list
      if (legend_justif =="inside") p <- p+ theme(legend.box = "horizontal")
    }
    if (length(table(df_trait_sp_F_mean_sd$Species))==1 &
        length(table(df_trait_sp_F_mean_sd$variable))>1 & facet==T & !(c("list") %in% class_p)){
      p_list=list()
      for (i in 1:ceiling(c((length(spec))/(ncol*nrow)))){
        p_list[[i]] <-  p+
          facet_wrap_paginate(~variable,nrow = nrow, ncol=ncol,drop=FALSE, scales = scale_facet, page=i)+
          theme(legend.position = "top", 
                legend.direction = "horizontal")
      }
      p <- p_list
    }
  }
  
  # Return both the ggplot object and the applied scaling factor.
  return(list(plot = p, scale = scale_factor, # [CORRIGÉ] renvoyait 'scale' (toujours 1)
              data = data_unscaled))            # [AJOUT] valeurs tracées, à l'échelle réelle
}


function_lollipop_decompositions <- function(data, sp_log = NA, 
                                             variables_to_bars, 
                                             variables_to_point, 
                                             colors, 
                                             val_names, 
                                             var_x = c("Fmax", "Fstart")){ 
# data in input : relative values, organised by replicates, variable, Species, Fstart and Fmax
  
  # init parameters
  var_to_point <- improve_vriable_names(variables_to_point)
  group_cols   <- c(var_x, "Species", "variable")
  if (all(var_x == c("Fmax", "Fstart"))) fishing_lab <- labeller(Fmax = c(`θF = 0.5`  = "Low fishing", `θF = 1.5` = "High fishing"))
  if (all(var_x == c("Scenario")))       fishing_lab <- labeller(Scenario = c(`CC`  = "Climate change"))
  
  # add log values to dataset
  data_lolli_plot_log          <- data
  data_lolli_plot_log$variable <- paste0(data_lolli_plot_log$variable, "_log") 
  data_lolli_plot_log$value    <- sign(data_lolli_plot_log$value)*log(abs(data_lolli_plot_log$value)+1)
  data_lolli_plot              <- rbind(data, data_lolli_plot_log)

  # average the replicates
  data_lolli_plot <- data_lolli_plot %>%
    group_by(across(all_of(group_cols))) %>%
    summarise(m  = mean(value, na.rm = TRUE),
              nrep = sum(!is.na(value)),
              se = sd(value, na.rm = TRUE)/sqrt(nrep),
              .groups = "drop")
  data_lolli_plot <- data_lolli_plot %>%
    mutate(IC_95 = qt(0.975, df = nrep-1) * se)
  data_lolli_plot$IC_text <- NA
  if (all(variables_to_bars == c("Ebx'",  "Elx'")) & any(na.omit(data_lolli_plot$IC_95)>0.55)){
    data_lolli_plot$IC_text[which((data_lolli_plot$IC_95>0.55))] <- data_lolli_plot$IC_95[which((data_lolli_plot$IC_95>0.55))]
    data_lolli_plot$IC_95[which((data_lolli_plot$IC_95>0.55))] <- 0.45
  }
	
  # organise dataset for figure
  data_lolli_plot$Species <- improve_sp_names(data_lolli_plot$Species)
  
  # for species collapsing with Fmax, replace their theta_F value by the theta_F' values without collapse + add a tag to the species name
  if (all(var_x == c("Fmax", "Fstart"))){
    for (sp in improve_sp_names(species)) {
      if (!any(sp %in% improve_sp_names(data$Species))) next
      cat(sp, dim(data_lolli_plot), "\n")
      Fmax_collapse <- Fmax_remove_collapses(sp)[1]
      if (!is.na(Fmax_collapse)){                                                       # [MODIFIÉ] no longer depends on rows existing at the collapsed θF
        # remove the collapsed θF, if there are any rows left for them
        id_rm <- which((as.numeric(data_lolli_plot$Fmax) >= Fmax_collapse) & (data_lolli_plot$Species == sp))
        if (length(id_rm) > 0) data_lolli_plot <- data_lolli_plot[-id_rm, ]              # [MODIFIÉ] guard: df[-integer(0), ] would empty the table
        # the last θF before collapse is used as the high fishing scenario
        data_lolli_plot[which((data_lolli_plot$Fmax == c(Fmax_collapse-0.25)) & (data_lolli_plot$Species == sp)), "Fmax"] <- c(1.5)
        # symbol indicating the collapse (always added)
        if (Fmax_collapse == 1.25) data_lolli_plot$Species[which(data_lolli_plot$Species == sp)] <-
          paste0(data_lolli_plot$Species[which(data_lolli_plot$Species == sp)], " ‡")
        if (Fmax_collapse == 1.5)  data_lolli_plot$Species[which(data_lolli_plot$Species == sp)] <-
          paste0(data_lolli_plot$Species[which(data_lolli_plot$Species == sp)], " †")
      }
    }
    # group species depending on their Vk. response to F
    sp_no_tag <- stringr::str_remove(data_lolli_plot$Species, " [‡†]$")                     
    data_lolli_plot$group <- ifelse(sp_no_tag %in% improve_sp_names(c("EuropeanPlaice", "CommonDab", "GreyGurnard")),
                                    "Vk_point ↑",
                                    ifelse(sp_no_tag %in% improve_sp_names(c("AtlanticHerring", "CommonSole",
                                                                             "Sandeel", "Saithe", "Hake", "Whiting",
                                                                             "EuropeanSprat", "NorwayPout")),
                                           "Vk_point ↓", "Vk_point constant"))
    # filter dataset 
    data_lolli_plot         <- data_lolli_plot %>% filter(Fmax %in% c(1.5, 0.5))
    # order species according to the relative response of Vk. at theta_F=0.5
    data_lolli_plot$Species <- factor(data_lolli_plot$Species, 
                                      levels=data_lolli_plot$Species[which(data_lolli_plot$Fmax == 0.5 & data_lolli_plot$variable == variables_to_point)][
                                        order(abs(data_lolli_plot$m[which(data_lolli_plot$Fmax == 0.5 & data_lolli_plot$variable == variables_to_point)]))])
    data_lolli_plot$Fmax    <- paste0("\u03B8F = ", data_lolli_plot$Fmax)
    ymin=0.55
    ymax=0.6
  }
  if (all(var_x == c("Scenario"))) {
    # group species depending on their Vk. response to CC
    data_lolli_plot$group <- c("No group at the moment")
    # order species 
    data_lolli_plot$Species <- factor(data_lolli_plot$Species, 
                                      levels=data_lolli_plot$Species[which(data_lolli_plot$Scenario == "CC" & data_lolli_plot$variable == variables_to_point)][
                                        order(abs(data_lolli_plot$m[which(data_lolli_plot$Scenario == "CC" & data_lolli_plot$variable == variables_to_point)]))])
    # filter dataset 
    data_lolli_plot         <- data_lolli_plot %>% filter(Scenario %in% c("CC"))
    ymin=0.75
    ymax=0.80
  }
  
  # pretty dataset
  data_lolli_plot$group   <- improve_vriable_names(data_lolli_plot$group)
  data_lolli_plot$variable<- improve_vriable_names(data_lolli_plot$variable)
  var_to_point = improve_vriable_names(var_to_point)
  variables_to_bars = improve_vriable_names(variables_to_bars)
  
  # add a variable to say if it is Q or P that is explaining deltaV 
  # => if Q and P change in the same direction, need to show the two colors on the same side of the lollipop
  data_lolli_plot_tt <- pivot_wider(data_lolli_plot[, c("Species", var_x, "group", "m", "variable")], 
                                    names_from = "variable", values_from = "m")
  data_lolli_plot_tt$sign <- ifelse(sign(data_lolli_plot_tt[, var_to_point]) == sign(data_lolli_plot_tt[,variables_to_bars[1]]), variables_to_bars[1], variables_to_bars[2])
  data_lolli_plot_tt$change_position_square <- 0
  if(any(sign(data_lolli_plot_tt[, variables_to_bars[1]]) == sign(data_lolli_plot_tt[,variables_to_bars[2]]))) {
    ID_same_s <- which(sign(data_lolli_plot_tt[, variables_to_bars[1]]) == sign(data_lolli_plot_tt[,variables_to_bars[2]]))
    data_lolli_plot_tt$sign[ID_same_s] <- ifelse(abs(data_lolli_plot_tt[ID_same_s, variables_to_bars[1]]) > 
                                                   abs(data_lolli_plot_tt[ID_same_s,variables_to_bars[2]]), variables_to_bars[1], variables_to_bars[2])
    data_lolli_plot_tt$change_position_square[ID_same_s] <- 0.05
  }
  
  # plot => MAIN TEXT
  position_label <- mean(abs(as.matrix(data_lolli_plot[which(data_lolli_plot$variable==var_to_point), "m"])))
  position_label <- c(-position_label-0.2, position_label+0.2)
  data_t <- data_lolli_plot[c(grep(data_lolli_plot$variable, pattern=stringr::str_remove(variables_to_bars[1], pattern = "'")), 
                              grep(data_lolli_plot$variable, pattern=stringr::str_remove(variables_to_bars[2], pattern = "'"))), ]
  data_t       <- full_join(data_t, data_lolli_plot_tt[, c("Species", var_x, "group", "sign", "change_position_square")], )
  data_ti<- data_t
  
  lollipop <-
    ggplot(data_ti[which(data_ti$variable %in% variables_to_bars), ], aes(x=Species, col=variable))  + 
    geom_hline(yintercept = 0, color =  "#4e4b5c", lty = "dotted")+
    geom_vline(xintercept = -Inf, colour = "grey75", linewidth = 0.4) +
    geom_linerange(data = data_lolli_plot[which(data_lolli_plot$variable==var_to_point), ],
                   aes(x=Species, ymin=m-IC_95, ymax=m+IC_95), col="darkblue", linewidth=0.6,
                   inherit.aes = FALSE) +
    geom_point(data = data_lolli_plot[which(data_lolli_plot$variable==var_to_point), ],
                   aes(x=Species, y=m, 
                       shape = paste0("relative ",                                   
                                      stringr::str_remove(stringr::str_remove(var_to_point, pattern="'"), pattern="\u0394"))), 
               size=3.5, col="darkblue")+
    geom_linerange(data= data_ti[which(data_ti$variable %in% variables_to_bars[1]), ],
                  aes(x=Species, 
                      ymin=c(sign(m)*(ymin+ifelse(variables_to_bars[1]==sign, change_position_square, 0))), 
                      ymax=c(sign(m)*(ymax+ifelse(variables_to_bars[1]==sign, change_position_square, 0))),
                      colour = variable), linewidth=10, 
                  alpha=ifelse(data_ti[which(data_ti$variable %in% 
                                               variables_to_bars[1]), ]$sign==variables_to_bars[1], 0.75, 0.25)) +  
    geom_linerange(data= data_ti[which(data_ti$variable %in% c(variables_to_bars[2])), ],
                  aes(x=Species, 
                      ymin=c(sign(m)*(ymin+ifelse(variables_to_bars[2]==sign, change_position_square, 0))), 
                      ymax=c(sign(m)*(ymax+ifelse(variables_to_bars[2]==sign, change_position_square, 0))),
                      colour = variable), linewidth=10,
                  alpha=ifelse(data_ti[which(data_ti$variable %in% 
                                               c(variables_to_bars[2])), ]$sign==variables_to_bars[2], 0.75, 0.25)) +
    scale_colour_manual(name = NULL, breaks = variables_to_bars, 
                       labels = val_names, values = colors)+
    scale_shape_manual("", breaks = paste0("relative ",                                   
                                           stringr::str_remove(stringr::str_remove(var_to_point, pattern="'"), pattern="\u0394")), 
                       values = 18)+
    scale_y_continuous(breaks = scales::breaks_pretty(n = 4), #breaks = c(-0.6, -0.3, 0, 0.3, 0.6),
                      sec.axis = dup_axis(breaks = position_label, 
                                          labels = c("Negative", "Positive")), name = NULL) +
    guides(alpha="none", color = guide_legend(override.aes = list(linewidth = 6))) +
    theme_bw() +
    theme(
     panel.grid.major.x = element_blank(),
     panel.border = element_blank(),
     strip.background.x = element_rect(fill = "grey95", color="white", linetype = "solid", linewidth = 1),
     strip.background.y = element_rect(fill = "grey95", color="white"),
     strip.text.x = element_text(size = 9.5),
     strip.text = element_text(face = "bold"),
     strip.placement = "outside",                                  
     axis.ticks.x.top = element_blank(),                        
     axis.text.x.top = element_text(face = "italic", colour = "grey30", size = 9),  
     axis.line.x.bottom = element_line(colour = "grey30", linewidth = 0.4),        
     axis.ticks.x.bottom = element_line(colour = "grey30", linewidth = 0.4),       
     legend.title = element_text(size = 9.5),
     panel.grid.minor = element_blank(),
     panel.grid.major = element_blank(),
     plot.caption = element_text(colour = "grey60", hjust = 1, size=10),
     panel.spacing.x=unit(0.6, "lines"),panel.spacing.y=unit(0.2, "lines"),
     legend.position = "top"
    ) +
    xlab("") +
    ylab("")+
    coord_flip() +
    facet_grid(group~.data[[var_x[1]]], scales="free", space="free_y", labeller = fishing_lab)+ 
    labs(caption = "†: θF=1.25 used as the high fishing scenario\n‡: θF=1 used as the high fishing scenario")
   
  
  # plot => SUPP TEXT (with P' and Q' values)
  
  # dataset
  data_lolli_plot$lab     <- ifelse(is.na(data_lolli_plot$m), round(data_lolli_plot$m, 2), 
                                    signif(data_lolli_plot$m, 2)) # set asides values before log
  if (!any(is.na(sp_log))){
    #replace species in sp_log by their log values of P, P' and Q, Q'
    data_lolli_plot[which(data_lolli_plot$Species %in% sp_log & 
                            (data_lolli_plot$variable == variables_to_bars[1] | 
                               data_lolli_plot$variable ==variables_to_bars[2])), c("m", "IC_95", "Species")] <- 
      data_lolli_plot[which(data_lolli_plot$Species %in% sp_log & 
                              (data_lolli_plot$variable == paste0(variables_to_bars[1], "_log") | 
                                 data_lolli_plot$variable == paste0(variables_to_bars[2], "_log"))), c("m", "IC_95", "Species")]
    # reorder species levels wih the new sp names
    for (s in sp_log) {
      s_log <- paste0(s, " (log)")
      if (s=="Atlantic mackerel") data_lolli_plot$Species <- recode_factor(data_lolli_plot$Species,  "Atlantic mackerel"="Atlantic mackerel (log)")
    }
  }
  data_t <- data_lolli_plot[c(grep(data_lolli_plot$variable, pattern=stringr::str_remove(variables_to_bars[1], pattern = "'")), 
                              grep(data_lolli_plot$variable, pattern=stringr::str_remove(variables_to_bars[2], pattern = "'"))), ]
  d_PQ   <- data_t[which(data_t$variable %in% variables_to_bars), ]
  
  # invisible point beyond each label so that the axis includes the text
  d_PQ   <- data_t[which(data_t$variable %in% variables_to_bars), ]
  d_PQ <- d_PQ %>%
    mutate(y_lab = m + sign(m) * dplyr::coalesce(IC_95, 0)) %>%       # label position (same as geom_text)
    group_by(.data[[var_x[1]]]) %>%                                                # one value axis per θF panel
    mutate(y_blank = y_lab + sign(m) * 0.35 * max(abs(y_lab), na.rm = TRUE)) %>%
    ungroup()
  
  dodge  <- position_dodge(width = 0.7)   
  lollipop_rel <-
    ggplot(d_PQ, aes(x=Species, y=m, col=variable,
                     group=factor(variable, levels = variables_to_bars)))  +  
    geom_hline(yintercept = 0, color =  "#4e4b5c", lty = "dotted")+
    geom_vline(xintercept = -Inf, colour = "grey75", linewidth = 0.4) +   
    geom_linerange(data = data_lolli_plot[which(data_lolli_plot$variable==var_to_point), ],
                   aes(x=Species, ymin=m-IC_95, ymax=m+IC_95), col="darkblue", linewidth=0.6,
                   inherit.aes = FALSE) +
    geom_linerange(aes(ymin=0, ymax=m), linewidth=4, alpha=0.35, position = dodge) +   
    geom_linerange(aes(ymin = m-IC_95, ymax = m+IC_95), color = "#4e4b5c", linewidth = 0.4, position = dodge) +
    geom_point(data = data_lolli_plot[which(data_lolli_plot$variable==var_to_point), ],
               aes(x=Species, y=m, 
                   shape= paste0("relative ", 
                                 stringr::str_remove(stringr::str_remove(var_to_point, pattern="'"), pattern="\u0394"))), 
                   size=3.5, col="darkblue", inherit.aes = FALSE)+
    geom_text(aes(y = y_lab, label = lab,                                   
                  hjust = ifelse(m >= 0, -0.15, 1.15)),
              col = "darkblue", size = 2.6, position = dodge, show.legend = FALSE) +
    geom_blank(aes(y = y_blank)) +                                          # extends the axis to include the labels
    scale_colour_manual("", breaks = variables_to_bars, values = colors, labels = val_names)+
    scale_shape_manual("", breaks  = paste0("relative ",                                   
                                           stringr::str_remove(stringr::str_remove(var_to_point, pattern="'"), pattern="\u0394")), 
                       values = 18)+
    scale_y_continuous(breaks = scales::breaks_pretty(n = 4)#,          
                      #expand = expansion(mult = 0.3)
    ) +                 
    theme_bw() +
    theme(
       panel.grid.major.x = element_blank(),
       panel.border = element_blank(),
       strip.background.x = element_rect(fill = "grey95", color="white", linetype = "solid", linewidth = 1),
       strip.background.y = element_rect(fill = "grey95", color="white"),
       strip.text.x = element_text(size = 9.5),
       strip.text = element_text(face = "bold"),
       axis.line.x.bottom = element_line(colour = "grey30", linewidth = 0.4),   
       axis.ticks.x.bottom = element_line(colour = "grey30", linewidth = 0.4),  
       legend.title = element_blank(),
       panel.grid.minor = element_blank(),
       panel.grid.major = element_blank(),
       plot.caption = element_text(colour = "grey60", hjust = 1, size=10),
       panel.spacing.x = unit(1.2, "lines"),
       legend.position = "top"                                     
    ) +
    xlab("") +
    ylab("")+
    coord_flip() +
    facet_grid(group~.data[[var_x[1]]], scales="free", space="free_y", labeller = fishing_lab)+ 
    labs(caption = "†: θF=1.25 used as the high fishing scenario\n‡: θF=1 used as the high fishing scenario")
  
  return(list("lollipop"=lollipop, "lollipop_rel"=lollipop_rel))
}



# OPERATIONS on DATASETS -------------------


# estimate the T/Vk ratio 
estimateT_Vk <- function(trait_data_full_i){
  vec_i <- c(grep(trait_data_full_i$variable, pattern="generation"), 
            grep(trait_data_full_i$variable, pattern="Vk"))
  df_trait_t <- trait_data_full_i[vec_i,]
  df_trait_t <- pivot_wider(df_trait_t, names_from = "variable", values_from = "value") 
  df_trait_t[, "T/Vk."] <- df_trait_t[, grep(pattern = "generation", colnames(df_trait_t))]/df_trait_t[, grep(pattern = "V", colnames(df_trait_t))]
  df_trait_t[which(is.na(as.matrix(df_trait_t[, "T/Vk."]))), "T/Vk."] <- 0
  add = "T/Vk."
  if ("cohort" %in% colnames(trait_data_full_i)) add <- c(add, "cohort")
  if ("age" %in% colnames(trait_data_full_i)) add <- c(add, "age")
  if ("min_age" %in% colnames(trait_data_full_i)) add <- c(add, "min_age")
  df_trait_t <- pivot_longer(df_trait_t[, c("Species", "replicate", "Fstart", "Fmax", add)], cols = "T/Vk.", values_to = "value", names_to = "variable")
  return(df_trait_t)
}

# estimate the Z 
estimateZ <- function(trait_data_full_i){
  vec_i <- c(grep(trait_data_full_i$variable, pattern="Mpred"),
            grep(trait_data_full_i$variable, pattern="Mnat"),
            grep(trait_data_full_i$variable, pattern="Mstar"), 
            grep(trait_data_full_i$variable, pattern="F"))
  df_trait_t <- trait_data_full_i[vec_i,]
  df_trait_t <- pivot_wider(df_trait_t, names_from = "variable", values_from = "value") 
  df_trait_t[, "Z"] <- rowSums(df_trait_t[, grep(pattern = "mean", colnames(df_trait_t))], na.rm = TRUE)
  df_trait_t <- pivot_longer(df_trait_t[, c("Species", "replicate", "Fstart", "Fmax", "age", "Z")], cols = "Z", values_to = "value", names_to = "variable")
  return(df_trait_t)
}

# estimate the Mnat
estimateMnat <- function(trait_data_full_i){
  vec_i <- c(grep(trait_data_full_i$variable, pattern="Mpred"),
             grep(trait_data_full_i$variable, pattern="Mnat"),
             grep(trait_data_full_i$variable, pattern="Mstar"))
  df_trait_t <- trait_data_full_i[vec_i,]
  df_trait_t <- pivot_wider(df_trait_t, names_from = "variable", values_from = "value") 
  df_trait_t[, "Mnat"] <- rowSums(df_trait_t[, grep(pattern = "mean", colnames(df_trait_t))], na.rm = TRUE)
  df_trait_t <- pivot_longer(df_trait_t[, c("Species", "replicate", "Fstart", "Fmax", "age", "Mnat")], cols = "Mnat", values_to = "value", names_to = "variable")
  return(df_trait_t)
}

# estimate the survival curve
estimateSurv <- function(trait_data_full_i){
  var_i <- c(var, var2, var3)
  post = "_deaths"
  if (!grepl("mean_of", var)) var <- paste0("tot", post)
  if (grepl("mean_of", var)) var <- paste0(stringr::str_remove(stringr::str_remove(var, pattern = "mean_of_"), pattern= "_per_age"), post)
  vec_i <- as.vector(unlist(sapply(unique(var_i), FUN = function(x) (which(trait_data_full_i$variable== x)))))
  df_trait_t <- trait_data_full_i[vec_i,]
  if(length(table(df_trait_t$variable)) != 3)  { 
    if (var == "total_M") cat("normal : ")
    if (any(grepl(pattern="square", df_trait_t$variable))) df_trait_t <- df_trait_t[-grep(pattern="square", df_trait_t$variable), ]
    if (length(table(df_trait_t$variable)) != 2) {
      vec_i <- as.vector(unlist(sapply(var_i, FUN = function(x) (which(trait_data_full_i$variable == x)))))
      df_trait_t <- trait_data_full_i[vec_i,]
    }
    if(length(table(df_trait_t$variable)) != 2) cat("too many or too little var")
  }
  df_trait_t   <- (pivot_wider(df_trait_t, names_from = "variable", values_from = "value"))
  df_trait_t_i <- as.matrix(df_trait_t)
  df_trait_t[,var] <- (as.numeric(df_trait_t_i[, var_i[1]])/as.numeric(df_trait_t_i[, "total_M"]))*as.numeric(df_trait_t_i[, "Dead"])
  df_trait_t <- pivot_longer(df_trait_t[, c("Species", "replicate", "Fstart", "Fmax", "age", var)], cols = var, values_to = "value", names_to = "variable")
  
  
  return(df_trait_t)
}

# estimate the individual weight
estimateWei <- function(trait_data_full_i){
  vec_i <- c(grep(trait_data_full_i$variable, pattern="abundanceByAge"), 
             grep(trait_data_full_i$variable, pattern="biomass"))
  df_trait_t <- trait_data_full_i[vec_i,]
  
  df_trait_t <- pivot_wider(df_trait_t, names_from = "variable", values_from = "value") 
  df_trait_t[, "weight"] <- df_trait_t[, grep(pattern = "biomass", colnames(df_trait_t))]/df_trait_t[, grep(pattern = "abund", colnames(df_trait_t))]
  df_trait_t <- pivot_longer(df_trait_t[, c("Species", "replicate", "Fstart", "Fmax", "age", "weight")], cols = "weight", values_to = "value", names_to = "variable")
  
  return(df_trait_t)
}

# Unified function for division and multiplication operations
estimate_BinaryOp <- function(trait_data_full_i, var, operator = "/", dim=c("Fmax", "Fstart")) {
  # Parse variable names based on operator
  if (operator == "/") {
    var_i <- unlist(stringr::str_split(pattern = "\nDiv", var))
    if (length(var_i) > 2) var_i <- c(paste(var_i[1:2], collapse = "\nDiv"), var_i[3])
  } else if (operator == "*") {
    if (grepl("\n", var)) var_i <- unlist(stringr::str_split(pattern = "\nx", var))
    else var_i <- unlist(stringr::str_split(pattern = "x", var))
    if (length(var_i) > 2) return(NULL)
  }
  
  # Extract relevant rows
  vec_i <- as.vector(unlist(sapply(var_i, FUN = function(x) grep(trait_data_full_i$variable, pattern = x))))
  df_trait_t <- trait_data_full_i[vec_i, ]
  if (length(table(df_trait_t$variable))>2){
      vec_i <- as.vector(unlist(sapply(var_i, FUN = function(x) which(trait_data_full_i$variable == x))))
      df_trait_t <- trait_data_full_i[vec_i, ]
  }
  
  # Fallback to exact matching if pattern matching doesn't yield 2 variables
  if (length(table(df_trait_t$variable)) != 2) {
    vec_i <- as.vector(unlist(sapply(var_i, FUN = function(x) which(trait_data_full_i$variable == x))))
    df_trait_t <- trait_data_full_i[vec_i, ]
    
    # For multiplication, exclude 'square' variables if encountered
    if (operator == "*" && length(table(df_trait_t$variable)) != 2) {
      if (any(grepl(pattern = "square", df_trait_t$variable))) 
        df_trait_t <- df_trait_t[-grep(pattern = "square", df_trait_t$variable), ]
    }
    
    if (length(table(df_trait_t$variable)) != 2) 
      cat("too many or too little var\n")
  }
  
  vari=NULL
  # do we have the thirty age classes as expected ? If not, append matrix 30 times to uniformise dataset structure
  if (length(table(df_trait_t$age[which(df_trait_t$variable == unique(df_trait_t$variable)[1])])) != 
      length(table(df_trait_t$age[which(df_trait_t$variable == unique(df_trait_t$variable)[2])]))){

    if (length(table(df_trait_t$age[which(df_trait_t$variable == unique(df_trait_t$variable)[1])]))%in%c(1,0)) vari=1
    if (length(table(df_trait_t$age[which(df_trait_t$variable == unique(df_trait_t$variable)[2])]))%in%c(1,0)) vari=2
    if (is.null(vari)) 
      warning("Incompatible age dimensions of ", paste(unique(df_trait_t$variable), collapse=" and "))
    if (!is.null(vari)){
            df_var_i <- df_trait_t[which(df_trait_t$variable == unique(df_trait_t$variable)[vari]), ]
      df_else <- df_trait_t[which(df_trait_t$variable != unique(df_trait_t$variable)[vari]), ]
      df_var = NULL
      for (age in unique(df_else$age)) {
        df_var_i$age <- age
        df_var <- rbind(df_var, df_var_i)
      }
      df_trait_t = rbind(df_var, df_else)
    }
  }
  
  # Compute on long format directly to avoid repeated pivot_wider/pivot_longer cycles.
  id_cols <- c("Species", "replicate", dim)
  if ("age" %in% colnames(trait_data_full_i)) id_cols <- c(id_cols, "age")
  if ("cohort" %in% colnames(trait_data_full_i)) id_cols <- c(id_cols, "cohort")
  if ("min_age" %in% colnames(trait_data_full_i)) id_cols <- c(id_cols, "min_age")
  
  v1_name <- if (any(df_trait_t$variable == var_i[1])) {
    var_i[1]
  } else {
    unique(df_trait_t$variable[grepl(var_i[1], df_trait_t$variable)])[1]
  }
  v2_name <- if (any(df_trait_t$variable == var_i[2])) {
    var_i[2]
  } else {
    unique(df_trait_t$variable[grepl(var_i[2], df_trait_t$variable)])[1]
  }

  d1 <- df_trait_t[df_trait_t$variable == v1_name, c(id_cols, "value")]
  colnames(d1)[colnames(d1) == "value"] <- "value_1"
  d2 <- df_trait_t[df_trait_t$variable == v2_name, c(id_cols, "value")]
  colnames(d2)[colnames(d2) == "value"] <- "value_2"

  df_trait_t <- na.omit(df_trait_t)
  df_trait_t <- dplyr::inner_join(d1, d2, by = id_cols)
  if (operator == "/") {
    df_trait_t$value <- df_trait_t$value_1 / df_trait_t$value_2
  } else if (operator == "*") {
    df_trait_t$value <- df_trait_t$value_1 * df_trait_t$value_2
  }
  df_trait_t$variable <- var
  df_trait_t <- df_trait_t[, c(id_cols, "variable", "value")]
  
  return(df_trait_t)
}

# Wrapper: estimate division (fecundity/survival ratio)
estimateDiv <- function(trait_data_full_i, var, dim=c("Fmax", "Fstart")) {
  estimate_BinaryOp(trait_data_full_i, var, operator = "/", dim)
}

# Wrapper: estimate multiplication
estimateX <- function(trait_data_full_i, var) {
  estimate_BinaryOp(trait_data_full_i, var, operator = "*")
}

# estimate the number of dead individuals, according to the mortality rate  
estimateMnumber <- function(trait_data_full_i, var, var2=c("total_M"), var3=c("Dead")){
  var_i <- c(var, var2, var3)
  post = "_deaths"
  if (!grepl("mean_of", var)) var <- paste0("tot", post)
  if (grepl("mean_of", var)) var <- paste0(stringr::str_remove(stringr::str_remove(var, pattern = "mean_of_"), pattern= "_per_age"), post)
  vec_i <- as.vector(unlist(sapply(unique(var_i), FUN = function(x) (which(trait_data_full_i$variable== x)))))
  df_trait_t <- trait_data_full_i[vec_i,]
  if(length(table(df_trait_t$variable)) != 3)  { 
    if (var == "total_M") cat("normal : ")
    if (any(grepl(pattern="square", df_trait_t$variable))) df_trait_t <- df_trait_t[-grep(pattern="square", df_trait_t$variable), ]
    if (length(table(df_trait_t$variable)) != 2) {
      vec_i <- as.vector(unlist(sapply(var_i, FUN = function(x) (which(trait_data_full_i$variable == x)))))
      df_trait_t <- trait_data_full_i[vec_i,]
    }
    if(length(table(df_trait_t$variable)) != 2) cat("too many or too little var")
  }
  df_trait_t   <- (pivot_wider(df_trait_t, names_from = "variable", values_from = "value"))
  df_trait_t_i <- as.matrix(df_trait_t)
  df_trait_t[,var] <- (as.numeric(df_trait_t_i[, var_i[1]])/as.numeric(df_trait_t_i[, "total_M"]))*as.numeric(df_trait_t_i[, "Dead"])
  if ("min_age" %in% colnames(trait_data_full_i)) {vart <- c(var, "min_age")
  } else vart <- var
  df_trait_t <- pivot_longer(df_trait_t[, c("Species", "replicate", "Fstart", "Fmax", "age", vart)], cols = var, values_to = "value", names_to = "variable")
  
  return(df_trait_t)
}

# estimate the survival curve
estimateSurvnumber <- function(trait_data_full_i, var, var2=c("total_M"), var3=c("Dead"),  var4=c("mean_of_abundanceByAge_per_age")){
  var_i <- c(var, var2, var3, var4)
  post = "_surv"
  if (!grepl("mean_of", var)) var <- paste0("tot", post)
  if (grepl("mean_of", var)) var <- paste0(stringr::str_remove(stringr::str_remove(var, pattern = "mean_of_"), pattern= "_per_age"), post)
  vec_i <- as.vector(unlist(sapply(unique(var_i), FUN = function(x) (which(trait_data_full_i$variable== x)))))
  df_trait_t <- trait_data_full_i[vec_i,]
  if(length(table(df_trait_t$variable)) != 4 & !(var == "tot_surv" & length(table(df_trait_t$variable)) == 3))  { 
    if (var == "total_M") cat("normal : ")
    if (any(grepl(pattern="square", df_trait_t$variable))) df_trait_t <- df_trait_t[-grep(pattern="square", df_trait_t$variable), ]
    if (length(table(df_trait_t$variable)) != 2) {
      vec_i <- as.vector(unlist(sapply(var_i, FUN = function(x) (which(trait_data_full_i$variable == x)))))
      df_trait_t <- trait_data_full_i[vec_i,]
    }
    if(length(table(df_trait_t$variable)) != 2) cat("too many or too little var")
  }
  df_trait_t   <- (pivot_wider(df_trait_t, names_from = "variable", values_from = "value"))
  df_trait_t_i <- as.matrix(df_trait_t)
  df_trait_t[,var] <- as.numeric(df_trait_t_i[, grep("abund", colnames(df_trait_t_i))]) - (as.numeric(df_trait_t_i[, var_i[1]])/as.numeric(df_trait_t_i[, "total_M"]))*as.numeric(df_trait_t_i[, "Dead"])
  df_trait_t <- pivot_longer(df_trait_t[, c("Species", "replicate", "Fstart", "Fmax", "age", var)], cols = var, values_to = "value", names_to = "variable")
  
  
  return(df_trait_t)
}

# estimate the sum  
estimateSum <- function(trait_data_full_i, var, soustraction=F){
  var_i <- unlist(stringr::str_split(pattern = "sum", var))
  if (grepl("\n", var)) var_i <- unlist(stringr::str_split(pattern = "\nsum", var))
  vec_i <- as.vector(unlist(sapply(var_i, FUN = function(x) (which(trait_data_full_i$variable== x)))))
  df_trait_t <- trait_data_full_i[vec_i,]
  
  df_trait_t <- pivot_wider(df_trait_t, names_from = "variable", values_from = "value") 
  if (soustraction == T) df_trait_t[, which(colnames(df_trait_t)== var_i[2])] <- -df_trait_t[, which(colnames(df_trait_t)== var_i[2])]
  df_trait_t[, var] <- df_trait_t[, which(colnames(df_trait_t)== var_i[1])]+df_trait_t[, which(colnames(df_trait_t)== var_i[2])]
  if ("min_age" %in% colnames(trait_data_full_i)) var<-c(var, "min_age") 
  df_trait_t <- pivot_longer(df_trait_t[, c("Species", "replicate", "Fstart", "Fmax", "age", var)], cols = var, values_to = "value", names_to = "variable")
  
  return(df_trait_t)
}

# estimate the dead matrix by age
estimateDead <- function(trait_data_full_i){
  vec_i <- c(grep(trait_data_full_i$variable, pattern="abundanceByAge"))
  df_trait_t <- trait_data_full_i[vec_i,]

  # df_trait_t <- pivot_wider(df_trait_t, names_from = "variable", values_from = "value")
  df_trait_t$age <- as.numeric(as.matrix(df_trait_t$age))
  df_trait_t<- df_trait_t%>%
    group_by(Species, Fmax, Fstart, replicate)%>%
    mutate(Dead = c(value[age %in% c(1:(max(age)-1))] - value[age %in% c(2:max(age))], 0))
  df_trait_t <- pivot_longer(df_trait_t[, c("Species", "replicate", "Fstart", "Fmax", "age", "Dead")], cols = "Dead", values_to = "value", names_to = "variable")

  return(df_trait_t)
}

# estimate the cumsum
estimateCumSum <- function(trait_data_full_i, var){
  var <- stringr::str_remove(pattern = "square", var)[1]
  var_i <- stringr::str_remove(pattern = "cum", var)[1]
  vec_i <- grep(trait_data_full_i$variable, pattern= var_i)
  df_trait_t <- trait_data_full_i[vec_i,]
  if (any(unique(df_trait_t$variable) %in% traits_to_plot[[traits_selected]])) df_trait_t = df_trait_t[-grep(df_trait_t$variable, pattern=traits_to_plot[[traits_selected]]), ]
  if (length(unique(df_trait_t$variable))>1)  df_trait_t = df_trait_t[which(df_trait_t$variable == unique(df_trait_t$variable)[1])]
  
  df_trait_t$age <- as.numeric(as.matrix(df_trait_t$age))
  df_trait_t <- df_trait_t[order(df_trait_t$age), ]
  df_trait_t <- df_trait_t%>%
    group_by(replicate, Species, Fmax, Fstart) %>%
    mutate(var = cumsum(value, na.rm=T))
  colnames(df_trait_t)[which(colnames(df_trait_t) == "var")] <- var[1]
  df_trait_t <- pivot_longer(df_trait_t[, c("Species", "replicate", "Fstart", "Fmax", "age", var[1])], cols = var, values_to = "value", names_to = "variable")

  return(df_trait_t)
}

# estimate the colSum
estimatecolSum <- function(trait_data_full_i, var, dim=c("Fstart", "Fmax")){
  add = var[1]
  grouping = c(dim, "replicate", "Species")
  if (any(grepl("cohort", colnames(trait_data_full_i)))) {
    add <- c(add, "cohort")
    grouping <- c(grouping, "cohort")}
  if (any(grepl("Scenario", colnames(trait_data_full_i)))) {grouping <- c(grouping, "Scenario")
  add <- c(add, "Scenario")
  }
  var_i <- stringr::str_remove(pattern = "colSum", var)[1]
  vec_i <- grep(trait_data_full_i$variable, pattern= var_i)
  df_trait_t <- trait_data_full_i[vec_i,]
  if (length(unique(df_trait_t$variable))>1)  df_trait_t = df_trait_t[which(df_trait_t$variable == var_i),]
  
  df_trait_t$age <- as.numeric(as.matrix(df_trait_t$age))
  df_trait_t <- df_trait_t[order(df_trait_t$age), ]
  df_trait_t <- df_trait_t%>%
    group_by(across(all_of(grouping))) %>%
    summarize(var = sum(value, na.rm = T)) %>%
    ungroup()
  colnames(df_trait_t)[which(colnames(df_trait_t) == "var")] <- var[1]
  df_trait_t <- pivot_longer(df_trait_t[, c("Species", "replicate", dim, add)], cols = var, values_to = "value", names_to = "variable")
  
  return(df_trait_t)
}

# estimate age 1
estimateage1 <- function(trait_data_full_i, var){
  var_i <- stringr::str_remove(pattern = "age1", var)[1]
  vec_i <- grep(trait_data_full_i$variable, pattern= var_i)
  df_trait_t <- trait_data_full_i[vec_i,]
  if (length(unique(df_trait_t$variable))>1)  {
    vec_i <- which(trait_data_full_i$variable == var_i)
    df_trait_t <- trait_data_full_i[vec_i,]
  }
  
  df_trait_t$age <- as.numeric(as.matrix(df_trait_t$age))
  df_trait_t <- df_trait_t[order(df_trait_t$age), ]
  df_trait_t <- df_trait_t%>%
    group_by(replicate, Species, Fmax, Fstart) %>%
    mutate(var = (value[which(age==1)]))
  colnames(df_trait_t)[which(colnames(df_trait_t) == "var")] <- var[1]
  df_trait_t <- pivot_longer(df_trait_t[, c("Species", "replicate", "Fstart", "Fmax", "age", var[1])], cols = var, values_to = "value", names_to = "variable")
  
  return(df_trait_t)
}

# estimate the square
estimateSquare <- function(trait_data_full_i, var){
  var_i <- stringr::str_remove(pattern = "square", var)[1]
  vec_i <- which(trait_data_full_i$variable == var_i)
  if (!(length(vec_i)>0))   vec_i <- grep(trait_data_full_i$variable, pattern=var_i)
  df_trait_t <- trait_data_full_i[vec_i,]

  df_trait_t <- pivot_wider(df_trait_t, names_from = "variable", values_from = "value") 
  df_trait_t[, var] <- df_trait_t[, ncol(df_trait_t)]^2
  colnames(df_trait_t)[which(colnames(df_trait_t) == "var")] <- var[1]
  add <- "age"
  if (any(grepl("cohort", colnames(trait_data_full_i)))) add <- c(add, "cohort")
  df_trait_t <- pivot_longer(df_trait_t[, c("Species", "replicate", "Fstart", "Fmax", add, colnames(df_trait_t)[ncol(df_trait_t)])], 
                             cols = colnames(df_trait_t)[ncol(df_trait_t)], values_to = "value", names_to = "variable")

  return(df_trait_t)
}
