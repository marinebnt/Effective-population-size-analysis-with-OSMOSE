
###################################
### Functions to plot the traits variation with Fmax or age
##################################


# TOOLS functions ----------------


# Ensure each dataset is in long format once, before entering the nested plotting loops.
as_long_if_needed <- function(df, columns_keep_pivot=colnames(df)[-c(1:4)]) {
  if ("variable" %in% colnames(df)) return(df)
  pivot_longer(df, cols = columns_keep_pivot, names_to = "variable", values_to = "value")
}

# remove facets
remove_facets <- function(plot, layout) {
  layout <- strsplit(layout, split = '\n')[[1]]
  layout <- lapply(layout, trimws)
  layout <- matrix(unlist(sapply(layout, strsplit, "")),
                   nrow = length(layout), byrow = T)
  layout <- which(layout == "#", arr.ind = TRUE)
  prm <- apply(layout,1,\(x) {
    c(glue::glue("panel-{x[1]}-{x[2]}"),
      glue::glue("axis-r-{x[1]}-{x[2]}"),
      glue::glue("axis-l-{x[1]}-{x[2]}"),
      glue::glue("axis-b-{x[1]}-{x[2]}"),
      glue::glue("axis-b-{x[2]}-{x[1]}"),
      glue::glue("strip-t-{x[2]}-{x[1]}"), 
      "xlab-b")
  })  
  # https://stackoverflow.com/a/30372692/1296582
  g <- ggplot2::ggplotGrob(plot)
  rm_grobs     <- g$layout$name %in% prm
  g$grobs[rm_grobs] <- NULL
  g$layout <- g$layout[!rm_grobs, ]
  ggpubr::as_ggplot(g)
}

# no nan, no na, no infinite
clean_values <- function(df, value_col = "value") {
  vals <- df[[value_col]]
  df[!(is.na(vals) | is.nan(vals) | is.infinite(vals)), , drop = FALSE]
}


# estiamte ssdg as the sum of the delta*Dead/N per age
add_SSDg <- function(df, relative=F){   
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
  
  return(df_m)
}


# weighted mean process
# build_weight_df <- function(df, weight_var = "abundanceByAge",
#                             key_cols = c("Species", "replicate", "Fstart", "Fmax", "cohort", "age")) {
#   out <- df[df$variable == weight_var, c(key_cols, "value"), drop = FALSE]
#   colnames(out)[colnames(out) == "value"] <- "weight"
#   out
# }

# weighted_reduce <- function(df, weight_df, join_cols, group_cols, drop_age_zero = T) {
#   data_with_w <- df %>%
#     dplyr::left_join(weight_df, relationship = "many-to-many", by=join_cols)
#   data_with_w <- data_with_w[data_with_w$weight != 0, , drop = FALSE]
#   data_with_w$age <- as.numeric(data_with_w$age)
#   if (drop_age_zero) data_with_w <- data_with_w[data_with_w$age != 0, , drop = FALSE]
#   data_with_w %>%
#     dplyr::group_by(dplyr::across(dplyr::all_of(group_cols))) %>%
#     dplyr::summarise(value = stats::weighted.mean(value, w = weight, na.rm = TRUE), .groups = "drop")
# }


# optimize legend position to avoid overlapping with the text or curves
estimate_inside_legend_position <- function(df, x_col, y_col,
                                            legend_labels = NULL,
                                            y_axis_name = "",
                                            z_axis_name = "") {
  x_col <- as.character(x_col)[1]
  y_col <- as.character(y_col)[1]

  if (is.null(df) || nrow(df) == 0 || !all(x_col %in% colnames(df)) || !all(y_col %in% colnames(df))) {
    return(c(0.85, 0.85))
  }

  to_numeric_axis <- function(v) {
    if (is.factor(v) || is.character(v)) {
      return(as.numeric(as.factor(v)))
    }
    as.numeric(v)
  }

  x <- to_numeric_axis(df[[x_col]])
  y <- as.numeric(df[[y_col]])
  keep <- is.finite(x) & is.finite(y)
  x <- x[keep]
  y <- y[keep]

  if (length(x) < 3 || length(unique(x)) < 2 || length(unique(y)) < 2) {
    return(c(0.85, 0.85))
  }

  xr <- range(x, na.rm = TRUE)
  yr <- range(y, na.rm = TRUE)
  x <- (x - xr[1]) / (xr[2] - xr[1])
  y <- (y - yr[1]) / (yr[2] - yr[1])

  n_items <- max(1, length(unique(as.character(legend_labels))))
  label_len <- max(nchar(as.character(legend_labels)), na.rm = TRUE)
  if (!is.finite(label_len)) label_len <- 10

  y_name_len <- max(nchar(unlist(strsplit(as.character(y_axis_name), "\\n", fixed = FALSE))), na.rm = TRUE)
  z_name_len <- max(nchar(unlist(strsplit(as.character(z_axis_name), "\\n", fixed = FALSE))), na.rm = TRUE)
  if (!is.finite(y_name_len)) y_name_len <- 0
  if (!is.finite(z_name_len)) z_name_len <- 0

  legend_w <- min(0.40, max(0.18, 0.15 + 0.007 * label_len))
  legend_h <- min(0.52, max(0.11, 0.09 + 0.032 * n_items))

  left_pad <- min(0.20, max(0.03, 0.02 + 0.004 * y_name_len))
  right_pad <- min(0.20, max(0.03, 0.02 + 0.004 * z_name_len))
  top_pad <- 0.04
  bottom_pad <- 0.04

  x_left <- min(0.95, max(0.05, left_pad + legend_w / 2))
  x_right <- min(0.95, max(0.05, 1 - right_pad - legend_w / 2))
  y_top <- min(0.95, max(0.05, 1 - top_pad - legend_h / 2))
  y_bottom <- min(0.95, max(0.05, bottom_pad + legend_h / 2))
  y_middle <- 0.5

  x_mid <- (x_left + x_right) / 2
  y_mid_top <- (y_top + y_middle) / 2
  y_mid_bottom <- (y_bottom + y_middle) / 2

  candidates <- rbind(
    c(x_left, y_top),
    c(x_right, y_top),
    c(x_left, y_bottom),
    c(x_right, y_bottom),
    c(x_left, y_middle),
    c(x_right, y_middle),
    c(x_mid, y_top),
    c(x_mid, y_bottom),
    c(x_mid, y_middle),
    c(x_left, y_mid_top),
    c(x_right, y_mid_top),
    c(x_left, y_mid_bottom),
    c(x_right, y_mid_bottom)
  )

  # Build a simple line proxy: ordered x-y points for each legend group, interpolated
  # along x to discourage placing the legend over curve trajectories.
  approx_group_lines <- function(df_local, x_name, y_name, group_vec) {
    gx <- to_numeric_axis(df_local[[x_name]])
    gy <- as.numeric(df_local[[y_name]])
    gg <- as.character(group_vec)
    keep2 <- is.finite(gx) & is.finite(gy) & !is.na(gg)
    gx <- gx[keep2]
    gy <- gy[keep2]
    gg <- gg[keep2]
    if (length(gx) < 3) return(list())
    gx <- (gx - xr[1]) / (xr[2] - xr[1])
    gy <- (gy - yr[1]) / (yr[2] - yr[1])

    groups_u <- unique(gg)
    out <- vector("list", length(groups_u))
    j <- 1
    for (g in groups_u) {
      idx <- which(gg == g)
      if (length(idx) < 2) next
      xx <- gx[idx]
      yy <- gy[idx]
      ord <- order(xx)
      xx <- xx[ord]
      yy <- yy[ord]
      dup <- duplicated(xx)
      if (any(dup)) {
        yy <- tapply(yy, xx, mean, na.rm = TRUE)
        xx <- as.numeric(names(yy))
        yy <- as.numeric(yy)
      }
      if (length(xx) < 2) next
      xg <- seq(max(0, min(xx)), min(1, max(xx)), length.out = 80)
      yg <- approx(x = xx, y = yy, xout = xg, method = "linear", rule = 2)$y
      out[[j]] <- cbind(x = xg, y = yg)
      j <- j + 1
    }
    out[seq_len(j - 1)]
  }

  group_proxy <- if (!is.null(legend_labels) && length(legend_labels) == nrow(df)) {
    legend_labels
  } else if ("variable" %in% colnames(df)) {
    df$variable
  } else {
    rep("all", nrow(df))
  }
  line_proxies <- approx_group_lines(df, x_col, y_col, group_proxy)

  score_candidate <- function(cx, cy) {
    in_box <- (abs(x - cx) <= legend_w / 2) & (abs(y - cy) <= legend_h / 2)
    in_halo <- (abs(x - cx) <= (legend_w / 2 + 0.03)) & (abs(y - cy) <= (legend_h / 2 + 0.03))

    # Distance penalty to curve proxies: close curves to legend center are costly.
    min_d2 <- Inf
    if (length(line_proxies) > 0) {
      for (lp in line_proxies) {
        d2 <- (lp[, "x"] - cx)^2 + (lp[, "y"] - cy)^2
        d2m <- min(d2, na.rm = TRUE)
        if (is.finite(d2m) && d2m < min_d2) min_d2 <- d2m
      }
    }
    dist_pen <- if (is.finite(min_d2)) {
      1 / (sqrt(min_d2) + 0.02)
    } else {
      0
    }

    5 * sum(in_box) + 2 * sum(in_halo) + 0.8 * dist_pen
  }

  scores <- apply(candidates, 1, function(v) score_candidate(v[1], v[2]))
  best <- which.min(scores)
  best_pos <- candidates[best, ]

  c(min(0.98, max(0.02, best_pos[1])), min(0.98, max(0.02, best_pos[2])))
}


resolve_legend_position <- function(position_legend, df, x_col, y_col,
                                    legend_labels = NULL,
                                    y_axis_name = "",
                                    z_axis_name = "") {
  if (is.numeric(position_legend) && length(position_legend) == 2) {
    return(c(min(0.98, max(0.02, position_legend[1])), min(0.98, max(0.02, position_legend[2]))))
  }

  if (is.character(position_legend) && position_legend %in% c("bottom", "top", "left", "right"))
    return(c(1, 0.5))
  
  if (is.character(position_legend) && length(position_legend) == 1) {
    if (tolower(position_legend) == "auto") {
      return(estimate_inside_legend_position(
        df = df,
        x_col = x_col,
        y_col = y_col,
        legend_labels = legend_labels,
        y_axis_name = y_axis_name,
        z_axis_name = z_axis_name
      ))
    }
  }

  estimate_inside_legend_position(
    df = df,
    x_col = x_col,
    y_col = y_col,
    legend_labels = legend_labels,
    y_axis_name = y_axis_name,
    z_axis_name = z_axis_name
  )
  
}


# Build a per-group slope + p-value summary string for use in plot subtitles.
format_slope_subtitle <- function(df, x_col = "Fmax", y_col = "mean", group_col = "variable") {
  # Drop rows with missing x/y before fitting models.
  df <- df[!is.na(df[[x_col]]) & !is.na(df[[y_col]]), ]
  # If too few rows, return a single NA label.
  if (nrow(df) < 2) return("Slope (lm mean ~ Fmax): NA")
  # Fit a simple linear model per group and format slope + p-value.
  groups <- unique(df[[group_col]])
  slope_texts <- character(0)
  for (g in groups) {
    dfg <- df[df[[group_col]] == g, ]
    if (any(is.na(dfg[[y_col]]))) dfg <- dfg[-which(is.na(dfg[[y_col]])),]
    if (any(is.infinite(dfg[[y_col]]))) dfg <- dfg[-which(is.infinite(dfg[[y_col]])),]
    # Require variation in x to avoid singular fits.
    if (nrow(dfg) < 2 || length(unique(dfg[[x_col]])) < 2 || all(is.na(dfg[[y_col]])) || all(is.infinite(dfg[[y_col]]))) {
      slope <- NA
      p_val <- NA
    } else {
      slope <- coef(lm(dfg[[y_col]] ~ dfg[[x_col]]))[2]
      coef_p_val <- lm(dfg[[y_col]] ~ dfg[[x_col]])
      coef_p_val <- as.matrix(coef_p_val$coefficients)
      if (nrow(coef_p_val)>1 & ncol(coef_p_val)>1) {p_val <- coef_p_val[2,4]
      } else {p_val <- NA}
    }
    # Append a formatted label for this group.
    slope_texts <- c(slope_texts, paste0(g, "=", round(slope, 2), "(p-val=", round(p_val, 2), ")\n"))
  }
  # One subtitle string for ggplot.
  paste0("Slope (lm mean ~ Fmax): ", paste(slope_texts, collapse = "; "))
}


# shorter and more explicit variable names
improve_parm_names <- function(vector){
  vector <- stringr::str_replace(vector, pattern="Fmax", replacement = "\u03B8F")
}

improve_vriable_names <- function(vector){
  vector <- stringr::str_replace(vector, pattern="b_rate_sc", replacement = "Birth rate (scaled)")  
  vector <- stringr::str_replace(vector, pattern="mean_of_abundanceMature_per_age", replacement = "Nmat_\u03C4")
  vector <- stringr::str_replace(vector, pattern="Mtot", replacement = "Mortality rate")
  vector <- stringr::str_replace(vector, pattern="abundanceByAge", replacement = "Abundance")
  vector <- stringr::str_replace(vector, pattern="mean_of_meanWeightByAge_per_age", replacement = "Individual weight")
  vector <- stringr::str_replace(vector, pattern="SSDgx_Mpred", replacement = "Predation SSD\u0262 (\u0394²Dead_Pred)")
  vector <- stringr::str_replace(vector, pattern="SSDgx_Mstar", replacement = "Starvation SSD\u0262 (\u0394²Dead_Starv)")
  vector <- stringr::str_replace(vector, pattern="SSDgx_F", replacement = "Fisheries SSD\u0262 (\u0394²Dead_F)")
  vector <- stringr::str_replace(vector, pattern="SSDgx_Mnat", replacement = "Natural SSD\u0262 (\u0394²Dead_Natural)")
  vector <- stringr::str_replace(vector, pattern="Numerator_N_alpha", replacement = "N\u03B1")
  vector <- stringr::str_remove(vector, pattern="x")      
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
  vector <- stringr::str_replace(vector, pattern="survivorship_curve", replacement = "lx")
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
  vector <- stringr::str_replace(vector, pattern="Dead", replacement = "Dead\u03C4")
  vector <- stringr::str_replace(vector, pattern="SSDg/Nα", replacement = "(SSD\u0262/N\u03B1)\u03C4 ((\u0394²Dead/N\u03B1)\u03C4)")
  vector <- stringr::str_replace(vector, pattern="SSDg", replacement = "SSD\u0262,\u03C4 ((\u0394\u03C4²Dead)\u03C4)")
  # vector <- expression(vector)
}

Fmax_remove_collapses<- function(sp){
  Fmax_coll=NA
  if (sp %in% c("AtlanticCod", "GreyGurnard", "Hake")) Fmax_coll <- c(1.25,1.5)
  if (sp %in% c("NorwayPout", "Sandeel", "Whiting")) Fmax_coll <- c(1.5)
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
                               reference_input=NULL, FORCElm=F, title = "regular", scale_facet="free",
                               # weight, 
                               chosen_title=NULL, subtitle=NULL,legend_justif="inside",
                               position_legend="auto", ncol=2, nrow=2, facet=T, scale_z = T) {
  
  plot_name_ch=c()
  plot_name <- paste0(as.character(relative_to_ref), as.character(proportional))
  if (plot_name=="TRUEFALSE") plot_name_ch <- c("relative")
  if (plot_name=="FALSEFALSE") plot_name_ch <- c("")
  if (plot_name=="FALSETRUE") plot_name_ch <- c("proportional")
  
  if (age_dim==T & any(sapply(df_trait_sp$variable, FUN = function(x) (all(is.na(df_trait_sp$age[which(df_trait_sp$variable == x)])))))) {
    v_null <- df_trait_sp$variable[which(sapply(df_trait_sp$variable, FUN = function(x) {
      (all(is.na(df_trait_sp$age[which(df_trait_sp$variable == x)])))
    }))]
    df_trait_sp <- df_trait_sp[-which(df_trait_sp$variable == v_null),]
  }
  
  var_x <- c("Fmax")
  if (age_dim) var_x <- c("age", "Fmax")
  if (sp_dim) var_x <- c(var_x, "Species")
  PARAM <- c(var_x, ifelse(any(grepl("cohort", colnames(df_trait_sp))), c("cohort"), c(NA)))
  if(any(is.na(PARAM))) PARAM <- PARAM[-which(is.na(PARAM))]
  if("age" %in% colnames(df_trait_sp)) PARAM <- c(PARAM, "age")
  PARAM <- unique(PARAM)
  
  # Filter to target species and cast key numeric columns.
  df_trait_sp_F        <- df_trait_sp[which(df_trait_sp$Species %in% spec),]
  if ("Fstart" %in% colnames(df_trait_sp)) df_trait_sp_F$Fstart <- as.numeric(df_trait_sp_F$Fstart)
  if (any(is.infinite(df_trait_sp_F$value))) df_trait_sp_F <- df_trait_sp_F[-which(is.infinite(df_trait_sp_F$value)), ]

  supp_name = ""
  if (length(unique(df_trait_sp_F$Fmax))==1) supp_name = paste0("Fmax = ", unique(df_trait_sp_F$Fmax))  
  if (length(unique(df_trait_sp_F$age))==1) supp_name = paste0("age = ", unique(df_trait_sp_F$age))
  if (age_dim)  {
    df_trait_sp_F$age   <- as.numeric(df_trait_sp_F$age)
    if (!any(grepl(trait_name, pattern="Mstar"), grepl(trait_name, pattern="Vsquare"),  grepl(trait_name, pattern="VsquarexDead")) & any(which(df_trait_sp_F$value == 0))) {
      df_trait_sp_F <- df_trait_sp_F[-which(df_trait_sp_F$value == 0), ]}
    if (any(is.na(df_trait_sp_F$value))) df_trait_sp_F <- df_trait_sp_F[-which(is.na(df_trait_sp_F$value)), ]
  }
  
  
  # Extract reference (Fmax=1, Fstart=1) if relative mode is on.
  if (is.null(reference_input)){
    if(!(c("cohort") %in% PARAM) & age_dim) ref_data <- df_trait_sp_F %>%
      filter(get(var_x[2]) == 1 & Fstart == 1) %>%
      dplyr::select(variable, value, replicate, age, Species) %>%
      dplyr::rename(ref_value = value)
    if(!age_dim) ref_data <- df_trait_sp_F %>%
      filter(Fmax == 1 & Fstart == 1) %>%
      dplyr::select(variable, value, replicate, Species) %>%
      dplyr::rename(ref_value = value)    
    if(all(c("age", "cohort") %in% PARAM)) ref_data <- df_trait_sp_F %>%
      filter(Fmax == 1 & Fstart == 1) %>%
      dplyr::select(variable, value, replicate, age, cohort, Species) %>%
      dplyr::rename(ref_value = value)
    if(any(colnames(ref_data) %in% c("Fmax", "Fstart"))) ref_data <- ref_data[, -which(colnames(ref_data) %in% c("Fmax", "Fstart"))]
  } else {
    ref_data <- reference_input
  }
  if ((relative_to_ref | proportional) & (nrow(ref_data) == 0)) cat("no data in ref_data")
  if ((relative_to_ref | proportional) & nrow(ref_data)>0) {
    if (any(ref_data$ref_value == 0, na.omit=T)) {
      ref_data$ref_value      <- ref_data$ref_value + 0.000001
      df_trait_sp_F$value <- df_trait_sp_F$value + 0.000001
    }
    df_trait_sp_F$replicate <- as.numeric(df_trait_sp_F$replicate)
    ref_data$replicate      <- as.numeric(ref_data$replicate)
    if (all(var_x %in% colnames(ref_data) & nrow(ref_data)>1)) {
      ref_data[, var_x[1]]       <- as.numeric(as.matrix(ref_data[, var_x[1]]))
      df_trait_sp_F <- df_trait_sp_F %>%
        left_join(ref_data, by = c("replicate", "variable","Species", PARAM[-which(PARAM=="Fmax")]), relationship = "many-to-many") ##############
    }
    if (all(nrow(ref_data)>1))  {
      # Convert to relative change vs reference per replicate/variable.
      df_trait_sp_F <- df_trait_sp_F %>%
        full_join(ref_data, relationship = "many-to-many") 
    }
    if (nrow(ref_data)>1 & relative_to_ref){
      if (length(var_x)==1) {
        df_trait_sp_F <- df_trait_sp_F %>%
        group_by(replicate, variable, .data[[var_x]]) %>%
        mutate(value =(value - ref_value)/abs(ref_value)) %>%
        dplyr::select(-ref_value)
      } else {
        df_trait_sp_F <- df_trait_sp_F %>%
          group_by(replicate, variable, .data[[var_x[1]]], .data[[var_x[2]]]) %>%
          mutate(value =(value - ref_value)/ref_value) %>%
          dplyr::select(-ref_value)
      } 
      if(any(is.nan(df_trait_sp_F$value))) df_trait_sp_F[which(is.nan(df_trait_sp_F$value)), "value"] <- 0
      if(any((na.omit(df_trait_sp_F$value)>10))) df_trait_sp_F[which(df_trait_sp_F$value>10), "value"] <- 10
    }
      if (nrow(ref_data)>1 & proportional){
        if (length(var_x)==1) {df_trait_sp_F <- df_trait_sp_F %>%
          group_by(replicate, variable, .data[[var_x]]) %>%
          mutate(value =(value)/ref_value) %>%
          dplyr::select(-ref_value)
        } else {
          df_trait_sp_F <- df_trait_sp_F %>%
            group_by(replicate, variable, .data[[var_x[1]]], .data[[var_x[2]]]) %>%
            mutate(value =(value)/ref_value) %>%
            dplyr::select(-ref_value)
        } 
    }
  }
  
  # arithmetic mean mean
  if (length(var_x)==1){
    df_trait_sp_F_mean_sd <- df_trait_sp_F %>%
      group_by(get(var_x), Fstart, variable) %>%
      # summarise(mean = weighted.mean(w=weight, x=value, na.rm = TRUE),      
      summarise(mean = mean(value, na.rm = TRUE),
                sd   = sd(value, na.rm = TRUE))
  } else {df_trait_sp_F_mean_sd <- df_trait_sp_F %>%
    group_by(get(var_x[1]),get(var_x[2]), Fstart, variable) %>%
    # summarise(mean = weighted.mean(w=weight, x=value, na.rm = TRUE),
    summarise(mean = mean(value, na.rm = TRUE),
              sd   = sd(value, na.rm = TRUE))}
  
    
  # update colnames
  if (any(colnames(df_trait_sp_F_mean_sd) == "get(var_x)")) colnames(df_trait_sp_F_mean_sd)[which(colnames(df_trait_sp_F_mean_sd) == "get(var_x)")] <- var_x
  if (any(colnames(df_trait_sp_F_mean_sd) == "get(var_x[1])")) colnames(df_trait_sp_F_mean_sd)[which(colnames(df_trait_sp_F_mean_sd) == "get(var_x[1])")] <- var_x[1]
  if (any(colnames(df_trait_sp_F_mean_sd) == "get(var_x[2])")) colnames(df_trait_sp_F_mean_sd)[which(colnames(df_trait_sp_F_mean_sd) == "get(var_x[2])")] <- var_x[2]
  
  
  # Normalize types for plotting and factor ordering.
  df_trait_sp_F_mean_sd$variable = as.factor(df_trait_sp_F_mean_sd$variable)
  df_trait_sp_F_mean_sd$Fstart   = as.factor(df_trait_sp_F_mean_sd$Fstart)
  df_trait_sp_F_mean_sd$mean     = as.numeric(df_trait_sp_F_mean_sd$mean)
  df_trait_sp_F_mean_sd[,var_x[1]]  = as.numeric(as.matrix(df_trait_sp_F_mean_sd[,var_x[1]]))
  if (length(var_x)>2) df_trait_sp_F_mean_sd[,var_x[2]]  = as.numeric(as.matrix(df_trait_sp_F_mean_sd[,var_x[2]]))
  df_trait_sp_F_mean_sd$mean     = as.numeric(df_trait_sp_F_mean_sd$mean)
  
  # what variable have scale y axis and which ones have scale secondary axis?
  scale = 1
  scale_factor <- 1
  shift <- 0
  var_scale_regular = NA
  
  if (nrow(df_trait_sp_F_mean_sd)==0) return(list(c(), c()))
  
  mean_var <- aggregate(df_trait_sp_F_mean_sd[, "mean"], list(df_trait_sp_F_mean_sd$variable), mean, na.rm=T)
  max_var  <- aggregate(df_trait_sp_F_mean_sd[, "mean"], list(df_trait_sp_F_mean_sd$variable), max, na.rm=T)
  var_scale_regular <- NULL
  var_scale_regular = mean_var[which((abs(mean_var[, 2]-max(mean_var[, 2])) > (0.8 * (max(mean_var[, 2])))) & !(mean_var[,2] == max(mean_var[, 2]))), 1]
  range=0
  if (length(var_scale_regular) == length(trait_name)) var_scale_regular <- NULL
  id_keep <- seq_len(nrow(df_trait_sp_F_mean_sd))
  max_first  <- max(df_trait_sp_F_mean_sd$mean[id_keep], na.rm = TRUE)   # Specify max of first y axis
  min_first  <- min(df_trait_sp_F_mean_sd$mean[id_keep], na.rm = TRUE)   # Specify min of first y axis
  scale_function <- function(x, scale, shift){
    return((x)*1)
  }
  inv_scale_function <-scale_function
  id_keep <- seq_len(nrow(df_trait_sp_F_mean_sd))[which(df_trait_sp_F_mean_sd$variable %in% var_scale_regular)]
  if (length(var_scale_regular) > 0 & ((length(df_trait_sp_F_mean_sd$mean[-id_keep]))!=0 ) & scale_z == T) { #https://finchstudio.io/blog/ggplot-dual-y-axes/
    max_first  <- max(df_trait_sp_F_mean_sd$mean[id_keep], na.rm = TRUE)   # Specify max of first y axis
    min_first  <- min(df_trait_sp_F_mean_sd$mean[id_keep], na.rm = TRUE)   # Specify min of first y axis
    min_second <- min(df_trait_sp_F_mean_sd$mean[-id_keep], na.rm = TRUE) # Specify min of second y axis
    max_second <- max(df_trait_sp_F_mean_sd$mean[-id_keep], na.rm = TRUE) # Specify max of second y axis
    if(max_second==0) max_second <- max_first
    # shift <- min(df_trait_sp_F_mean_sd$mean[id_keep], na.rm = T) - min(df_trait_sp_F_mean_sd$mean[-id_keep], na.rm = T) * max_first/max_second
    shift=0
    scale_factor <- diff(range(df_trait_sp_F_mean_sd$mean[id_keep], na.rm = T)) / diff(range(df_trait_sp_F_mean_sd$mean[-id_keep], na.rm = T))
    if (is.nan(scale_factor)) scale_factor <- diff(range(df_trait_sp_F_mean_sd$mean[id_keep][-which(is.nan(df_trait_sp_F_mean_sd$mean[id_keep]))], na.rm = T)) /
      diff(range(df_trait_sp_F_mean_sd$mean[-id_keep][-which(is.nan(df_trait_sp_F_mean_sd$mean[-id_keep]))], na.rm = T))
    if (is.infinite(scale_factor) || is.nan(scale_factor) || scale_factor==0) scale_factor <- 1
    # Function to scale secondary axis
    scale_function <- function(x, scale, shift){
      return ((x-shift)/scale)
    }
    # Function to scale secondary variable values
    inv_scale_function <- function(x, scale, shift){
      return ((x)*scale + shift)
    }
    df_trait_sp_F_mean_sd$mean[-id_keep] <-
      # inv_scale_function(df_trait_sp_F$mean[-id_keep], scale_factor, shift)
      inv_scale_function(df_trait_sp_F_mean_sd$mean[-id_keep], scale_factor, shift)
    df_trait_sp_F_mean_sd$sd[-id_keep] <-
      # inv_scale_function(df_trait_sp_F$mean[-id_keep], scale_factor, shift)
      inv_scale_function(df_trait_sp_F_mean_sd$sd[-id_keep], scale_factor, shift)
  }
  if (!(all(seq_len(nrow(df_trait_sp_F_mean_sd)) %in% id_keep))) id_keep <- seq_len(nrow(df_trait_sp_F_mean_sd))[-id_keep]
  max_second <- max(c(df_trait_sp_F_mean_sd$mean[id_keep]), na.rm = TRUE) # Specify max of second y axis
  min_second <- min(c(df_trait_sp_F_mean_sd$mean[id_keep]), na.rm = TRUE) # Specify min of second y axis
  if (any(is.infinite(max_first)))  max_first  <- max(df_trait_sp_F_mean_sd$mean[-id_keep][-which(is.infinite(df_trait_sp_F_mean_sd$mean[-id_keep]))], na.rm = TRUE)   # Specify max of first y axis
  if (any(is.infinite(max_second)))  max_second  <- max(df_trait_sp_F_mean_sd$mean[id_keep][-which(is.infinite(df_trait_sp_F_mean_sd$mean[id_keep]))], na.rm = TRUE)   # Specify max of first y axis
  if (any(is.infinite(min_first)))  min_first  <- min(df_trait_sp_F_mean_sd$mean[-id_keep][-which(is.infinite(df_trait_sp_F_mean_sd$mean[-id_keep]))], na.rm = TRUE)   # Specify max of first y axis
  if (any(is.infinite(min_second)))  min_second  <- min(df_trait_sp_F_mean_sd$mean[id_keep][-which(is.infinite(df_trait_sp_F_mean_sd$mean[id_keep]))], na.rm = TRUE)   # Specify max of first y axis
  
  # Decide which variable goes on the secondary axis (if there is one = TRUE)
  z_axis <- ifelse(length(var_scale_regular) > 0, TRUE, FALSE)
  if(z_axis)  {z_axis <- paste(trait_name[-which(trait_name %in% var_scale_regular)], collapse = "\n")
  } else z_axis <- ""
  if (nchar(as.character(z_axis[1]))>1) {
    name_y <- paste(var_scale_regular, collapse="\n")
  }else  name_y <- paste(mean_var[,1] , collapse="\n")
  limits = c(min(c(df_trait_sp_F_mean_sd$mean), na.rm = T), 
             max(c(df_trait_sp_F_mean_sd$mean), na.rm = T))
  # Build slope summary text for the subtitle.
  slope_subtitle <- format_slope_subtitle(df = df_trait_sp_F, x_col = var_x[1], y_col = "value", group_col = "variable")
  if (is.null(subtitle)) subtitle <- ifelse(nchar(paste(spec, slope_subtitle, sep = " | ", supp_name))>100, FALSE, TRUE)
  if (subtitle) {
    subtitle = paste(spec, slope_subtitle, sep = " | ", supp_name)
  } else { subtitle = spec}
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

 #  # Assemble the plot with loess trends and a secondary axis.
  fill_vector=c("white", "black", "grey", "grey30", "lightyellow", "grey20", "orange", "yellow2", "yellow4")
  shape_vector=c(21:31, 12)  
  line_vector=c(1:9)
  if(length(fill_vector[-which(is.na(names(fill_vector)))])>6) {
    shape_vector <- rep(21, length(fill_vector))
    fill_vector <- rep("white", length(fill_vector))
    line_vector <- rep(1, length(fill_vector))
  }
  fill_vectorfmax <- fill_vector
  fill_vector <- setNames(fill_vector, levels(df_trait_sp_F_mean_sd$variable))
  fill_vectorfmax <- setNames(fill_vectorfmax, levels(df_trait_sp_F_mean_sd$Fmax))
  shape_vector <- setNames(shape_vector, levels(df_trait_sp_F_mean_sd$variable))
  line_vector <- setNames(line_vector, levels(df_trait_sp_F_mean_sd$variable))
  fill_vectori <- as.character(df_trait_sp_F_mean_sd$variable)
  fill_vectorifmax <- as.character(df_trait_sp_F_mean_sd$Fmax)
  shape_vectori <- as.character(df_trait_sp_F_mean_sd$variable)
  line_vectori <- as.character(df_trait_sp_F_mean_sd$variable)
  tempo <- sapply(levels(df_trait_sp_F_mean_sd$variable), FUN = function(x) which(fill_vectori==x))
  tempo2 <- sapply(levels(df_trait_sp_F_mean_sd$Fmax), FUN = function(x) which(fill_vectori==x))
  if (is.list(tempo)) {
    for (i in names(tempo))  fill_vectori[fill_vectori==i] <- fill_vector[i]
    for (i in names(tempo2)) fill_vectorifmax[fill_vectorifmax==i] <- fill_vectorfmax[i]
    for (i in names(tempo))  shape_vectori[shape_vectori==i] <- shape_vector[i]
    for (i in names(tempo))  line_vectori[line_vectori==i] <- line_vector[i]
  } else{
    for (i in colnames(tempo))  fill_vectori[tempo[,i]] <- fill_vector[i]
    for (i in names(tempo2)) fill_vectorifmax[fill_vectorifmax==i] <- fill_vectorfmax[i]
    for (i in colnames(tempo))  shape_vectori[tempo[,i]] <- shape_vector[i]
    for (i in colnames(tempo))  line_vectori[tempo[,i]] <- line_vector[i]
  }

  shape_var <- c(length(table(df_trait_sp_F_mean_sd$variable))>1, length(table(df_trait_sp_F_mean_sd$Fmax))>1)
  if (all(shape_var)>0) {shape_var <- "variable" 
  } else {
    shape_var <- c("variable", "Fmax")[which.max( c(length(table(df_trait_sp_F_mean_sd$variable)), length(table(df_trait_sp_F_mean_sd$Fmax))))]}
  df_trait_sp_F_mean_sd$variable <- as.factor(df_trait_sp_F_mean_sd$variable)
  df_trait_sp_F_mean_sd$variable <- factor(df_trait_sp_F_mean_sd$variable)
  
  for(sp in spec) {
    Fmax_coll <- Fmax_remove_collapses(sp)
    if (!all(is.na(Fmax_coll))) {
      if (any("Species" %in% colnames(df_trait_sp_F_mean_sd))) {
        df_temp <- df_trait_sp_F_mean_sd[which(df_trait_sp_F_mean_sd$Species == sp), ]
        df_temp[which(df_temp$Fmax %in% Fmax_coll), "mean"] <- NA
        df_temp[which(df_temp$Fmax %in% Fmax_coll), "sd"] <- NA
        df_trait_sp_F_mean_sd[which(df_trait_sp_F_mean_sd$Species == sp), ] <- df_temp
      } else {
        df_trait_sp_F_mean_sd[which(df_trait_sp_F_mean_sd$Fmax %in% Fmax_coll), "mean"] <- NA
        df_trait_sp_F_mean_sd[which(df_trait_sp_F_mean_sd$Fmax %in% Fmax_coll), "sd"] <- NA
      }
    }
  }
  
  if (length(var_x) == 1) {
    # remove the variables that have only one age class or NAs in the ssdg sums for instance
    if (any(sapply(df_trait_sp_F_mean_sd$variable, FUN = function(x) (all(is.na(df_trait_sp_F_mean_sd$mean[which(df_trait_sp_F_mean_sd$variable == x)])))))) {
      v_null <- df_trait_sp_F_mean_sd$variable[which(sapply(df_trait_sp_F_mean_sd$variable, FUN = function(x) {
        cat(x)
        (all(is.na(df_trait_sp_F_mean_sd$mean[which(df_trait_sp_F_mean_sd$variable == x)])))
      }))]
      add_null <- ifelse(any(is.na(df_trait_sp_F_mean_sd$Fmax)), which(is.na(df_trait_sp_F_mean_sd$Fmax)), c())
      df_trait_sp_F_mean_sd <- df_trait_sp_F_mean_sd[-c(which(df_trait_sp_F_mean_sd$variable == v_null), add_null),]
    }
    legend_pos <- resolve_legend_position(
      position_legend = position_legend,
      df = df_trait_sp_F_mean_sd,
      x_col = var_x,
      y_col = "mean",
      legend_labels = df_trait_sp_F_mean_sd[[shape_var]],
      y_axis_name = name_y,
      z_axis_name = z_axis
    )
    p <- ggplot(df_trait_sp_F_mean_sd, aes(x = .data[[var_x]], y = mean,
                group = interaction(variable, Fmax))) +
      geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd, group = .data[[shape_var]], 
                        alpha=as.numeric(.data[[shape_var]])/10), 
                    show.legend = FALSE, 
                    position = position_dodge(0.1), width = 0.2)+
      geom_point(data=df_trait_sp_F_mean_sd, aes(shape = .data[[shape_var]], fill = .data[[shape_var]]), 
                size=2.5, position = position_dodge(0.1),  
                stroke=1.2)+ #, , stroke=1.5
      scale_shape_manual(name="",values=as.numeric(shape_vectori)) +
      scale_fill_manual(name="",values=as.character(fill_vectori)) +
    # geom_smooth(method = method, se=TRUE, col = "black",
    #   aes(group = interaction(variable, Fstart))
    # ) + 
    # geom_point(size = 1.5, alpha=0.5) +
    theme(panel.background = element_rect(fill = "white", colour = "grey50"), 
          plot.title = element_text(hjust = 0.5), 
          axis.text.x = element_text(angle = 60, hjust = 1),
          axis.text = element_text(colour = 1, size = 10),
          panel.border = element_rect(colour = "black", fill=NA),
          legend.background = element_blank(),
          legend.justification = c("right", "top"),
          legend.title = element_blank(), 
          legend.text=element_text(size=9), 
          legend.box.background = element_rect(colour = "black"), 
          legend.key = element_rect(colour = NA))+
    # facet_wrap(~Fstart, labeller = Fstart_names) +
    labs(
      title = title,
      x = improve_parm_names(var_x),
      y =  name_y, 
      subtitle = spec
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
    
    if (all(legend_pos !=c(1,0.5))) {
      p <- p+
        theme(legend.position = legend_justif,
              legend.position.inside = legend_pos,
              legend.spacing.y = unit(0, "mm")) 
    }
    
    if (length(df_trait_sp_F_mean_sd$variable)>1 & facet==T){
      p <- p+
        facet_wrap(~variable,nrow = nrow, ncol=ncol,drop=FALSE, scales = scale_facet)+
        theme(legend.position.inside = c(1,0.3), 
              legend.direction = "vertical") + 
        ylab(label = "")
      
      if (all(legend_pos !=c(1,0.5))) p <- p+ theme(legend.box = "horizontal")
    }
    
    if (relative_to_ref == T | proportional == T){
      p <- p + 
        theme(legend.position = "inside",
              legend.position.inside = c(1,1),
              legend.spacing.y = unit(0, "mm")) +
        geom_hline(col="grey80", linetype="dotted", lwd=0.8, alpha=0.6, yintercept=0)+
        geom_hline(col="grey80", linetype="dotted", lwd=0.8, alpha=0.6, yintercept=inv_scale_function(0, scale = scale, shift = shift))
    }
    p <- p +
      geom_line(aes(group = .data[[shape_var]], alpha=as.numeric(.data[[shape_var]])/10), show.legend = FALSE)
  } else {
    df_trait_sp_F_mean_sd$variable <- factor(df_trait_sp_F_mean_sd$variable)
    df_trait_sp_F$Species <- as.factor(df_trait_sp_F$Species)
    legend_pos <- resolve_legend_position(
      position_legend = position_legend,
      df = df_trait_sp_F_mean_sd,
      x_col = improve_parm_names(var_x[1]),
      y_col = "mean",
      legend_labels = df_trait_sp_F_mean_sd$Species,
      y_axis_name = name_y,
      z_axis_name = z_axis
    )
    if (facet==T & legend_justif=="inside") legend_pos <- c(0.9, 0.4)
    p <- ggplot(df_trait_sp_F_mean_sd, aes(x = .data[[var_x[1]]], y = mean, #mean,
                                   shape = .data[[var_x[2]]], colour=variable, group = interaction(variable, .data[[var_x[2]]]))) +
      # geom_smooth(method = method, se=TRUE) +
      # geom_pointrange(aes(ymin=mean-sd, ymax=mean+sd), size=0.5)+
      geom_point(aes(y=mean, x=.data[[var_x[1]]]), size=2)+
      geom_line(linewidth=0.7)+
      theme_bw() +
      geom_hline(col="grey", linetype="dotted", lwd=0.8, alpha=0.4, yintercept=0)+
      geom_hline(col="grey", linetype="dotted", lwd=0.8, alpha=0.4, yintercept=inv_scale_function(0, scale = scale, shift = shift))+
      # facet_wrap(~Fstart, labeller = Fstart_names) +
      labs(
        title = title,
        x = improve_parm_names(var_x),
        y =  name_y, 
        subtitle = subtitle,
        shape=improve_parm_names()
      ) +
      theme(panel.background = element_rect(fill = "white", colour = "grey50"), 
            plot.title = element_text(hjust = 0.5), 
            axis.text.x = element_text(angle = 60, hjust = 1),
            axis.text = element_text(colour = 1, size = 10),
            panel.border = element_rect(colour = "black", fill=NA),
            legend.background = element_blank(),
            legend.justification = c("right", "top"),
            legend.title = element_text(size=13),
            legend.text=element_text(size=12), 
            legend.box.background = element_rect(colour = "black"), 
            legend.key = element_rect(colour = NA))
    
    if (length(levels(df_trait_sp_F_mean_sd$variable))==1) {
      p <- p + 
        guides(colour=FALSE) + 
        geom_point(size=3, colour="black") +
        geom_line(colour="black")
    }
    if (!(z_axis == "")){
      p <- p+
        scale_y_continuous(
          # limits = limits,
          sec.axis = sec_axis(~scale_function(., scale_factor, shift), #. * scale, # (. - a)/b,  #sec_axis(~. * scale, 
                              name = z_axis
          )
        )
    }
    
    if (all(legend_pos !=c(1,0.5)) || (facet==T & legend_justif=="inside")) {
      p <- p+
        theme(legend.position = legend_justif,
              legend.position.inside = legend_pos,
              legend.spacing.y = unit(0, "mm")) 
    }
    
    if (length(df_trait_sp_F_mean_sd$variable)>1 & facet==T){
      p <- p+
        facet_wrap(~variable,nrow = nrow, ncol=ncol,drop=FALSE, scales = scale_facet)+
          theme(legend.position.inside = c(1,0.3), 
              legend.direction = "vertical")
      
      if (all(legend_pos !=c(1,0.5)) & legend_justif =="inside") p <- p+ theme(legend.box = "horizontal")
    }
    # if (relative_to_ref == T | proportional == T){
    #   p <- p + 
    #     theme(legend.position = legend_justif,
    #           legend.position.inside = c(1,1),
    #           legend.spacing.y = unit(0, "mm")) +
    #     geom_hline(col="grey80", linetype="dotted", lwd=0.8, alpha=0.6, yintercept=0)+
    #     geom_hline(col="grey80", linetype="dotted", lwd=0.8, alpha=0.6, yintercept=inv_scale_function(0, scale = scale, shift = shift))
    # }
  }
  
  # Return both the ggplot object and the applied scaling factor.
  return(list(plot = p, scale = scale))
}







plot_trait_facet_fmax <- function(df_trait_sp, trait_name, Time1, Time2, relative_to_ref = FALSE, 
                               spec, age_dim=FALSE, reference_input=NULL, FORCElm=F, title = "regular", legend=NULL, #weight,
                               position_legend="auto", facet=F) {
  
  var_x <- ifelse(age_dim, "age", "Fmax")
  PARAM <- c(var_x, ifelse(any(grepl("cohort", colnames(df_trait_sp))), c("cohort"), c(NA)))
  if(any(is.na(PARAM))) PARAM <- PARAM[1]
  
  # Filter to target species and cast key numeric columns.
  df_trait_sp_F        <- df_trait_sp[which(df_trait_sp$Species %in% spec),]
  if ("Fstart" %in% colnames(df_trait_sp)) df_trait_sp_F$Fstart <- as.numeric(df_trait_sp_F$Fstart)
  df_trait_sp_F$Fmax   <- as.numeric(df_trait_sp_F$Fmax)
  if (any(is.infinite(df_trait_sp_F$value))) df_trait_sp_F <- df_trait_sp_F[-which(is.infinite(df_trait_sp_F$value)), ]
  
  supp_name = ""
  if (length(unique(df_trait_sp_F$Fmax))==1) supp_name = paste0("Fmax = ", unique(df_trait_sp_F$Fmax))  
  if (length(unique(df_trait_sp_F$age))==1) supp_name = paste0("age = ", unique(df_trait_sp_F$age))
  if (age_dim)  {
    df_trait_sp_F$age   <- as.numeric(df_trait_sp_F$age)
    if (!any(grepl(trait_name, pattern="Mstar"), grepl(trait_name, pattern="Vsquare"),  grepl(trait_name, pattern="VsquarexDead")) & any(which(df_trait_sp_F$value == 0))) {
      df_trait_sp_F <- df_trait_sp_F[-which(df_trait_sp_F$value == 0), ]}
    if (any(is.na(df_trait_sp_F$value))) df_trait_sp_F <- df_trait_sp_F[-which(is.na(df_trait_sp_F$value)), ]
  }
  
  # Extract reference (Fmax=1, Fstart=1) if relative mode is on.
  if (is.null(reference_input)){
    ref_data <- df_trait_sp_F %>%
      filter(Fstart == 1 & Fmax == 1) %>%
      dplyr::select(variable, value, replicate) %>%
      dplyr::rename(ref_value = value)
  } else {
    ref_data <- reference_input
  }
  if (relative_to_ref & (nrow(ref_data) == 0)) cat("no data in ref_data")
  if (relative_to_ref & nrow(ref_data)>0) {
    ref_data <- na.omit(ref_data)
    if (any(ref_data$ref_value == 0)) {
      ref_data$ref_value      <- ref_data$ref_value + 0.000001
      df_trait_sp_F$value <- df_trait_sp_F$value + 0.000001
    }
    df_trait_sp_F$replicate <- as.numeric(df_trait_sp_F$replicate)
    ref_data$replicate      <- as.numeric(ref_data$replicate)
    if (var_x %in% colnames(ref_data) & nrow(ref_data)>1) {
      ref_data[, var_x]       <- as.numeric(as.matrix(ref_data[, var_x]))
      df_trait_sp_F <- df_trait_sp_F %>%
        left_join(ref_data, by = c("replicate", "variable", PARAM), relationship = "many-to-many") ##############
    }
    if (nrow(ref_data)>1 & !(var_x %in% colnames(ref_data)))  {
      # Convert to relative change vs reference per replicate/variable.
      add = "variable"
      if ("age" %in% colnames(ref_data)) add = c("variable", "age")
      df_trait_sp_F <- df_trait_sp_F %>%
        full_join(ref_data, by = c("Species", "replicate", add)) 
    }
    if (nrow(ref_data)>1){
      df_trait_sp_F <- df_trait_sp_F %>%
        group_by(replicate, variable, .data[[var_x]]) %>%
        mutate(value =(value - ref_value)/ref_value) %>%
        dplyr::select(-ref_value)
    }
    if(any(is.nan(df_trait_sp_F$value))) df_trait_sp_F[which(is.nan(df_trait_sp_F$value)), "value"] <- 0
    if(any(na.omit(df_trait_sp_F$value>10))) df_trait_sp_F[which(df_trait_sp_F$value>10), "value"] <- 10
  }
  ref_data$Fstart <- 1
  
  # Aggregate across replicates for each Fmax/Fstart/variable. 
  # colnames(weight)[which(colnames(weight) == "value")] <- "weight"
  # weight <- weight[,-which(colnames(weight) == "variable")]
  # id_merge <- which((colnames(df_trait_sp_F) %in% colnames(WEIGHT)) & (colnames(df_trait_sp_F) != "variable"))
  # df_trait_sp_F$age <- as.numeric(df_trait_sp_F$age)
  # df_trait_sp_F <- df_trait_sp_F %>%
  #   left_join(weight, by = colnames(df_trait_sp_F)[id_merge], relationship = "many-to-many") ##############
  if (length(var_x)==1) df_trait_sp_F_mean_sd <- df_trait_sp_F %>%
    group_by(get(var_x), Fstart, variable, Fmax) %>%
    # summarise(mean = weighted.mean(w=weight, value, na.rm = TRUE),
    summarise(mean = mean(value, na.rm = TRUE),
              sd   = sd(value, na.rm = TRUE))#%>%
    # dplyr::select(-weight)  
  if (length(var_x)!=1) df_trait_sp_F_mean_sd <- df_trait_sp_F %>%
    group_by(get(var_x[1]), get(var_x[2]), Fstart, variable, Fmax) %>%
    # summarise(mean = weighted.mean(w=weight, value, na.rm = TRUE),
    summarise(mean = mean(value, na.rm = TRUE),
              sd   = sd(value, na.rm = TRUE))#%>%
    # dplyr::select(-weight)
  if (any(colnames(df_trait_sp_F_mean_sd) == "get(var_x)")) colnames(df_trait_sp_F_mean_sd)[which(colnames(df_trait_sp_F_mean_sd) == "get(var_x)")] <- var_x
  
  # Normalize types for plotting and factor ordering.
  df_trait_sp_F_mean_sd$variable = as.factor(df_trait_sp_F_mean_sd$variable)
  df_trait_sp_F_mean_sd$Fstart   = as.factor(df_trait_sp_F_mean_sd$Fstart)
  df_trait_sp_F_mean_sd$mean     = as.numeric(df_trait_sp_F_mean_sd$mean)
  df_trait_sp_F_mean_sd[,var_x]  = as.numeric(as.matrix(df_trait_sp_F_mean_sd[,var_x]))
  df_trait_sp_F_mean_sd$mean     = as.numeric(df_trait_sp_F_mean_sd$mean)
  
  # Build slope summary text for the subtitle.
  slope_subtitle <- format_slope_subtitle(df = df_trait_sp_F_mean_sd, x_col = var_x, y_col = "mean", group_col = "variable")
  
  # what variable have scale y axis and which ones have scale secondary axis?
  scale = 1
  scale_factor <- 1
  shift <- 0
  var_scale_regular = NA
  
  if (nrow(df_trait_sp_F_mean_sd)==0) return(list(c(), c()))
  
  mean_var <- aggregate(df_trait_sp_F_mean_sd[, "mean"], list(df_trait_sp_F_mean_sd$variable), mean, na.rm=T)
  max_var  <- aggregate(df_trait_sp_F_mean_sd[, "mean"], list(df_trait_sp_F_mean_sd$variable), max, na.rm=T)
  var_scale_regular <- NULL
  var_scale_regular = mean_var[which(mean_var[, 2] > (0.8 * (max(mean_var[, 2])))), 1]
  if (length(var_scale_regular) == length(trait_name)) var_scale_regular <- NULL
  id_keep <- seq_len(nrow(df_trait_sp_F_mean_sd))
  max_first  <- max(df_trait_sp_F_mean_sd$mean[id_keep], na.rm = TRUE)   # Specify max of first y axis
  min_first  <- min(df_trait_sp_F_mean_sd$mean[id_keep], na.rm = TRUE)   # Specify min of first y axis
  scale_function <- function(x, scale, shift){
    return((x)*1)
  }
  inv_scale_function <-scale_function
  if (length(var_scale_regular) > 0) { #https://finchstudio.io/blog/ggplot-dual-y-axes/
    id_keep <- seq_len(nrow(df_trait_sp_F_mean_sd))[which(df_trait_sp_F_mean_sd$variable %in% var_scale_regular)]
    max_first  <- max(df_trait_sp_F_mean_sd$mean[id_keep], na.rm = TRUE)   # Specify max of first y axis
    min_first  <- min(df_trait_sp_F_mean_sd$mean[id_keep], na.rm = TRUE)   # Specify min of first y axis
    min_second <- min(df_trait_sp_F_mean_sd$mean[-id_keep], na.rm = TRUE) # Specify min of second y axis
    max_second <- max(df_trait_sp_F_mean_sd$mean[-id_keep], na.rm = TRUE) # Specify max of second y axis
    if(max_second==0) max_second <- max_first
    shift <- min(df_trait_sp_F_mean_sd$mean[id_keep], na.rm = T) - min(df_trait_sp_F_mean_sd$mean[-id_keep], na.rm = T) * max_first/max_second
    scale_factor <- diff(range(df_trait_sp_F_mean_sd$mean[id_keep], na.rm = T)) / diff(range(df_trait_sp_F_mean_sd$mean[-id_keep], na.rm = T))
     if (is.nan(scale_factor)) scale_factor <- diff(range(df_trait_sp_F_mean_sd$mean[id_keep][-which(is.nan(df_trait_sp_F_mean_sd$mean[id_keep]))], na.rm = T)) /
      diff(range(df_trait_sp_F_mean_sd$mean[-id_keep][-which(is.nan(df_trait_sp_F_mean_sd$mean[-id_keep]))], na.rm = T))
    if (is.infinite(scale_factor) || is.nan(scale_factor) || scale_factor==0) scale_factor <- 1
     # Function to scale secondary axis
    scale_function <- function(x, scale, shift){
      return ((x-shift)/scale)
  }
    # Function to scale secondary variable values
    inv_scale_function <- function(x, scale, shift){
      return ((x)*scale + shift)
    }
    df_trait_sp_F$value[-id_keep] <-
      # inv_scale_function(df_trait_sp_F$mean[-id_keep], scale_factor, shift)
      inv_scale_function(df_trait_sp_F_mean_sd$mean[-id_keep], scale_factor, shift)
  }
  if (!(all(seq_len(nrow(df_trait_sp_F_mean_sd)) %in% id_keep))) id_keep <- seq_len(nrow(df_trait_sp_F_mean_sd))[-id_keep]
  max_second <- max(c(df_trait_sp_F_mean_sd$mean[id_keep]), na.rm = TRUE) # Specify max of second y axis
  min_second <- min(c(df_trait_sp_F_mean_sd$mean[id_keep]), na.rm = TRUE) # Specify min of second y axis
  if (any(is.infinite(max_first)))  max_first  <- max(df_trait_sp_F_mean_sd$mean[-id_keep][-which(is.infinite(df_trait_sp_F_mean_sd$mean[-id_keep]))], na.rm = TRUE)   # Specify max of first y axis
  if (any(is.infinite(max_second)))  max_second  <- max(df_trait_sp_F_mean_sd$mean[id_keep][-which(is.infinite(df_trait_sp_F_mean_sd$mean[id_keep]))], na.rm = TRUE)   # Specify max of first y axis
  if (any(is.infinite(min_first)))  min_first  <- min(df_trait_sp_F_mean_sd$mean[-id_keep][-which(is.infinite(df_trait_sp_F_mean_sd$mean[-id_keep]))], na.rm = TRUE)   # Specify max of first y axis
  if (any(is.infinite(min_second)))  min_second  <- min(df_trait_sp_F_mean_sd$mean[id_keep][-which(is.infinite(df_trait_sp_F_mean_sd$mean[id_keep]))], na.rm = TRUE)   # Specify max of first y axis

  # Decide which variable goes on the secondary axis (if there is one = TRUE)
  z_axis <- ifelse(length(var_scale_regular) > 0, TRUE, FALSE)
  if(z_axis)  {z_axis <- paste0(var_names[-which(trait_name %in% var_scale_regular)], collapse=" ")
  } else z_axis <- ""
  if (nchar(as.character(z_axis[1]))>1) {
    name_y <- paste0(var_scale_regular, collapse=" ")
  }else  name_y <- paste(mean_var[,1] , collapse="\n")
  limits = c(min(c(df_trait_sp_F_mean_sd$mean), na.rm = T), 
             max(c(df_trait_sp_F_mean_sd$mean), na.rm = T))
  subtitle <- ifelse(nchar(paste(spec, slope_subtitle, sep = " | ", supp_name))>100, FALSE, TRUE)
  if (subtitle) {
    subtitle = paste(spec, slope_subtitle, sep = " | ", supp_name)
  } else { subtitle = ""}
  title = paste0(paste(trait_name, collapse=" "), " variation per age over y 150 to 160")
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
  df_trait_sp_F_mean_sd$age = as.numeric(df_trait_sp_F_mean_sd$age)
  df_trait_sp_F_mean_sd$variable = as.factor(df_trait_sp_F_mean_sd$variable)
  legend_pos <- resolve_legend_position(
    position_legend = position_legend,
    df = df_trait_sp_F_mean_sd,
    x_col = var_x,
    y_col = "mean",
    legend_labels = df_trait_sp_F_mean_sd$variable,
    y_axis_name = name_y,
    z_axis_name = z_axis
  )
  p <- ggplot(df_trait_sp_F_mean_sd, aes(x = .data[[var_x]], y = mean, #mean,
                                 colour = variable, group = interaction(Fmax,variable))) +
    # geom_smooth(method = method, se=TRUE, 
    #             aes(group = interaction(variable, Fstart))
    # ) +
    geom_line()+
    geom_pointrange(aes(ymin=mean-sd, ymax=mean+sd), size=0.5)+
    # geom_point(size = 1.5, alpha=0.5) +
    theme_bw() +
    geom_hline(col="grey50", linetype="dotted", lwd=0.8, alpha=0.6, yintercept=0)+
    geom_hline(col="grey50", linetype="dotted", lwd=0.8, alpha=0.6, yintercept=inv_scale_function(0, scale = scale, shift = shift))+
    # facet_wrap(~Fstart, labeller = Fstart_names) +
    labs(
      title = title,
      x = var_x,
      y =  name_y, 
      subtitle = spec
    ) +
    theme(axis.text.x = element_text(angle = 60, hjust = 1),
          legend.position = "inside",
          legend.position.inside = legend_pos,
          legend.justification = c("right", legend_justif)) +
    scale_y_continuous(
      # limits = limits,
      sec.axis = sec_axis(~scale_function(., scale_factor, shift), #. * scale, # (. - a)/b,  #sec_axis(~. * scale, 
                          name = z_axis
      )
    )+
    facet_wrap(vars(Fmax), ncol=ncol, nrow=nrow, scales="free_y", labeller=labeller)
  
  if (legend=="remove") p <- p+theme(legend.position = "none")
  if (facet==T){
    p <- p+
      facet_wrap(~variable, scales="free")+
      theme(legend.position.inside = c(1,1))
  }
  
  # Return both the ggplot object and the applied scaling factor.
  return(list(plot = p, scale = scale))
}








# Plot traits with a secondary axis scaling and optional relative reference.
# assumes the data are already aggregated (mean, sd). It doesn’t compute scaling (you pass scale_factor) 
# and returns only the plot. 
# It uses ribbons for SD and a simpler labeling scheme, and it doesn’t add slope text.
plot_trait_vs_fmax_simple <- function(df_trait_sp, sp_name, var_names, Time1, Time2, relative_to_ref = FALSE, scale_factor = 1) {
  # Extract reference (Fmax=1, Fstart=1) if relative mode is on.
  if (relative_to_ref) {
    ref_data <- df_trait_sp %>%
      filter(Fmax == 1 & Fstart == "1") %>%
      dplyr::select(variable, mean) %>%
      dplyr::rename(ref_mean = mean)
    
    # Convert to relative change vs reference per variable.
    df_trait_sp <- df_trait_sp %>%
      left_join(ref_data[, -which(colnames(ref_data) %in% c("Fmax", "Fstart"))], by = "variable", relationship = "many-to-one") %>%
      mutate(mean =(mean - ref_mean)/ref_mean) %>%
      dplyr::select(-ref_mean)
  }

  legend_pos <- resolve_legend_position(
    position_legend = "auto",
    df = df_trait_sp,
    x_col = "Fmax",
    y_col = "mean",
    legend_labels = df_trait_sp$variable,
    y_axis_name = ifelse(relative_to_ref, paste(var_names, "- RELATIVE", collapse = "\n"), paste(var_names, collapse = "\n")),
    z_axis_name = "T/Vk."
  )

  # Build the plot with ribbons for variability.
  p <- ggplot(df_trait_sp, aes(x = Fmax, group = variable, colour=variable)) +
    geom_ribbon(aes(ymin = mean - sd, ymax = mean + sd),
      alpha = 0.5, fill = "grey70"
    ) +
    geom_smooth(method = "loess", se = TRUE,
                aes(group = interaction(variable, Fstart), y=mean)
    ) +
    geom_point(aes(pch = as.factor(variable), y = mean)) +
    theme_bw() +
    labs(
      title = paste0(paste(var_names, collapse = " and "), " mean\nover ", Time1, " and ", Time2),
      x = "Fmax",
      y = ifelse(relative_to_ref, paste(var_names, "- RELATIVE", collapse = "\n"), paste(var_names, collapse = "\n")),
      subtitle = sp_name
    ) +
    theme(axis.text.x = element_text(angle = 60, hjust = 1),
          legend.position = "inside",
          legend.position.inside = legend_pos,
          legend.justification = c("right", "top")) +
    scale_y_continuous(
      sec.axis = sec_axis(~. * scale_factor, name = "T/Vk.")
    )
  
  # Return the ggplot object for downstream customization.
  return(p)
}


plot_trait_vs_age_simple <- function(df_trait_sp, 
                                     sp_name, 
                                     var_names, 
                                     Time1, 
                                     Time2, 
                                     relative_to_ref = FALSE, 
                                     proportional= FALSE, 
                                     scale_factor = 1, 
                                     # weight,
                                     position_legend="auto", 
                                     facet=F, scale_facet=NULL,
                                     chosen_title = NA, 
                                     z_axis=F, nrow=2, ncol=3) {
  
  # Extract reference (Fmax=1, Fstart=1) if relative mode is on
  if ("cohort" %in% colnames(df_trait_sp)){
    ref_data <- df_trait_sp %>%
      filter(Fmax == 1 & Fstart == "1") %>%
      dplyr::select(value,variable, age, replicate, Species, cohort) %>%
      dplyr::rename(ref_value = value)
  } else {
    ref_data <- df_trait_sp %>%
      filter(Fmax == 1 & Fstart == "1") %>%
      dplyr::select(value,variable, age, replicate, Species) %>%
      dplyr::rename(ref_value = value)
  }

  if (relative_to_ref) {
    # Convert to relative change vs reference per variable.
    df_trait_sp <- df_trait_sp %>%
      left_join(ref_data, relationship = "many-to-many") %>%
      mutate(value =(value - ref_value)/ref_value) %>%
      dplyr::select(-ref_value)
  }  
  if (proportional) {
    # Convert to relative change vs reference per variable.
    df_trait_sp <- df_trait_sp %>%
      left_join(ref_data, relationship = "many-to-many") %>%
      mutate(value =(value)/ref_value) %>%
      dplyr::select(-ref_value)
  }
  
  # colnames(weight)[which(colnames(weight)=="value")] <- "weight"
  # if ("variable" %in% colnames(weight)) weight <- weight[, -which(colnames(weight)=="variable")]
  # df_trait_sp <- df_trait_sp %>%
    # left_join(weight, relationship = "many-to-one")
  
  
  df_trait_sp_summ <- df_trait_sp %>%
    group_by(Fmax, Fstart, variable, age) %>%
    # summarise(mean = weighted.mean(w=weight, x=value, na.rm = TRUE),
    summarise(mean = mean(value, na.rm = TRUE),
              sd   = sd(value, na.rm = TRUE))
  # if ("SSDg" %in% colnames(df_trait_sp)) {
  #   df_trait_sp_summ <- full_join(df_trait_sp_summ, 
  #                       df_trait_sp[, which(colnames(df_trait_sp) %in% c(colnames(df_trait_sp_summ), "SSDg"))])
  # }
  df_trait_sp_summ$age <- as.numeric(df_trait_sp_summ$age)
  df_trait_sp_summ$Fmax <- as.factor(df_trait_sp_summ$Fmax)
  df_trait_sp_summ <- clean_values(df_trait_sp_summ, value_col = "mean")
  
  # # what variable have scale y axis and which ones have scale secondary axis?
  # scale = 1
  # scale_factor = 1
  # shift = 0
  # var_scale_regular = NA
  # if (nrow(df_trait_sp_summ)==0) return(list(c(), c()))
  # mean_var <- aggregate(df_trait_sp_summ[, "mean"], list(df_trait_sp_summ$variable), mean, na.rm=T)
  # max_var  <- aggregate(df_trait_sp_summ[, "mean"], list(df_trait_sp_summ$variable), max, na.rm=T)
  # var_scale_regular <- NULL
  # var_scale_regular = mean_var[which(mean_var[, 2] > (0.8 * (max(mean_var[, 2])))), 1]
  # if (length(var_scale_regular) == length(var_names)) var_scale_regular <- NULL
  # id_keep <- seq_len(nrow(df_trait_sp_summ))
  # max_first  <- max(df_trait_sp_summ$mean[id_keep], na.rm = TRUE)   # Specify max of first y axis
  # min_first  <- min(df_trait_sp_summ$mean[id_keep], na.rm = TRUE)   # Specify min of first y axis
  # scale_function <- function(x, scale, shift){
  #   return((x)*1)
  # }
  # inv_scale_function <-scale_function
  # if (length(var_scale_regular) > 0 & z_axis==T) { #https://finchstudio.io/blog/ggplot-dual-y-axes/
  #   id_keep <- seq_len(nrow(df_trait_sp_summ))[which(df_trait_sp_summ$variable %in% var_scale_regular)]
  #   max_first  <- max(df_trait_sp_summ$mean[id_keep], na.rm = TRUE)   # Specify max of first y axis
  #   min_first  <- min(df_trait_sp_summ$mean[id_keep], na.rm = TRUE)   # Specify min of first y axis
  #   min_second <- min(df_trait_sp_summ$mean[-id_keep], na.rm = TRUE) # Specify min of second y axis
  #   max_second <- max(df_trait_sp_summ$mean[-id_keep], na.rm = TRUE) # Specify max of second y axis
  #   if(max_second==0) max_second <- max_first
  #   shift <- min(df_trait_sp_summ$mean[id_keep], na.rm = T) - min(df_trait_sp_summ$mean[-id_keep], na.rm = T) * max_first/max_second
  #   scale_factor <- diff(range(df_trait_sp_summ$mean[id_keep], na.rm = T)) / diff(range(df_trait_sp_summ$mean[-id_keep], na.rm = T))
  #   if (is.nan(scale_factor)) scale_factor <- diff(range(df_trait_sp_summ$mean[id_keep][-which(is.nan(df_trait_sp_summ$mean[id_keep]))], na.rm = T)) /
  #     diff(range(df_trait_sp_summ$mean[-id_keep][-which(is.nan(df_trait_sp_summ$mean[-id_keep]))], na.rm = T))
  #   if (is.infinite(scale_factor) || is.nan(scale_factor) || scale_factor==0) scale_factor <- 1
  #   # Function to scale secondary axis
  #   scale_function <- function(x, scale, shift){
  #     return ((x-shift)/scale)
  #   }
  #   # Function to scale secondary variable values
  #   inv_scale_function <- function(x, scale, shift){
  #     return ((x)*scale + shift)
  #   }
  #   df_trait_sp_summ$mean[-id_keep] <-
  #     # inv_scale_function(df_trait_sp_summ_F$mean[-id_keep], scale_factor, shift)
  #     inv_scale_function(df_trait_sp_summ$mean[-id_keep], scale_factor, shift)
  # }
  # if (!(all(seq_len(nrow(df_trait_sp_summ)) %in% id_keep))) id_keep <- seq_len(nrow(df_trait_sp_summ))[-id_keep]
  # max_second <- max(c(df_trait_sp_summ$mean[id_keep]), na.rm = TRUE) # Specify max of second y axis
  # min_second <- min(c(df_trait_sp_summ$mean[id_keep]), na.rm = TRUE) # Specify min of second y axis
  # if (any(is.infinite(max_first)))  max_first  <- max(df_trait_sp_summ$mean[-id_keep][-which(is.infinite(df_trait_sp_summ$mean[-id_keep]))], na.rm = TRUE)   # Specify max of first y axis
  # if (any(is.infinite(max_second)))  max_second  <- max(df_trait_sp_summ$mean[id_keep][-which(is.infinite(df_trait_sp_summ$mean[id_keep]))], na.rm = TRUE)   # Specify max of first y axis
  # if (any(is.infinite(min_first)))  min_first  <- min(df_trait_sp_summ$mean[-id_keep][-which(is.infinite(df_trait_sp_summ$mean[-id_keep]))], na.rm = TRUE)   # Specify max of first y axis
  # if (any(is.infinite(min_second)))  min_second  <- min(df_trait_sp_summ$mean[id_keep][-which(is.infinite(df_trait_sp_summ$mean[id_keep]))], na.rm = TRUE)   # Specify max of first y axis
  #   
  # what variable have scale y axis and which ones have scale secondary axis?
  scale = 1
  scale_factor <- 1
  shift <- 0
  var_scale_regular = NA
  
  if (nrow(df_trait_sp_summ)==0) return(list(c(), c()))
  
  mean_var <- aggregate(df_trait_sp_summ[, "mean"], list(df_trait_sp_summ$variable), mean, na.rm=T)
  max_var  <- aggregate(df_trait_sp_summ[, "mean"], list(df_trait_sp_summ$variable), max, na.rm=T)
  var_scale_regular <- NULL
  var_scale_regular = mean_var[which((abs(mean_var[, 2]-max(mean_var[, 2])) > (0.8 * (max(mean_var[, 2])))) & !(mean_var[,2] == max(mean_var[, 2]))), 1]
  range=0
  if (length(var_scale_regular) == length(var_names)) var_scale_regular <- NULL
  id_keep <- seq_len(nrow(df_trait_sp_summ))
  max_first  <- max(df_trait_sp_summ$mean[id_keep], na.rm = TRUE)   # Specify max of first y axis
  min_first  <- min(df_trait_sp_summ$mean[id_keep], na.rm = TRUE)   # Specify min of first y axis
  scale_function <- function(x, scale, shift){
    return((x)*1)
  }
  inv_scale_function <-scale_function
  id_keep <- seq_len(nrow(df_trait_sp_summ))[which(df_trait_sp_summ$variable %in% var_scale_regular)]
  if (length(var_scale_regular) > 0 & ((length(df_trait_sp_summ$mean[-id_keep]))!=0 )) { #https://finchstudio.io/blog/ggplot-dual-y-axes/
    max_first  <- max(df_trait_sp_summ$mean[id_keep], na.rm = TRUE)   # Specify max of first y axis
    min_first  <- min(df_trait_sp_summ$mean[id_keep], na.rm = TRUE)   # Specify min of first y axis
    min_second <- min(df_trait_sp_summ$mean[-id_keep], na.rm = TRUE) # Specify min of second y axis
    max_second <- max(df_trait_sp_summ$mean[-id_keep], na.rm = TRUE) # Specify max of second y axis
    if(max_second==0) max_second <- max_first
    shift <- min(df_trait_sp_summ$mean[id_keep], na.rm = T) - min(df_trait_sp_summ$mean[-id_keep], na.rm = T) * max_first/max_second
    scale_factor <- diff(range(df_trait_sp_summ$mean[id_keep], na.rm = T)) / diff(range(df_trait_sp_summ$mean[-id_keep], na.rm = T))
    if (is.nan(scale_factor)) scale_factor <- diff(range(df_trait_sp_summ$mean[id_keep][-which(is.nan(df_trait_sp_summ$mean[id_keep]))], na.rm = T)) /
      diff(range(df_trait_sp_summ$mean[-id_keep][-which(is.nan(df_trait_sp_summ$mean[-id_keep]))], na.rm = T))
    if (is.infinite(scale_factor) || is.nan(scale_factor) || scale_factor==0) scale_factor <- 1
    # Function to scale secondary axis
    scale_function <- function(x, scale, shift){
      return ((x-shift)/scale)
    }
    # Function to scale secondary variable values
    inv_scale_function <- function(x, scale, shift){
      return ((x)*scale + shift)
    }
    df_trait_sp_summ$mean[-id_keep] <-
      # inv_scale_function(df_trait_sp_F$mean[-id_keep], scale_factor, shift)
      inv_scale_function(df_trait_sp_summ$mean[-id_keep], scale_factor, shift)
    df_trait_sp_summ$sd[-id_keep] <-
      # inv_scale_function(df_trait_sp_F$mean[-id_keep], scale_factor, shift)
      inv_scale_function(df_trait_sp_summ$sd[-id_keep], scale_factor, shift)
  }
  if (!(all(seq_len(nrow(df_trait_sp_summ)) %in% id_keep))) id_keep <- seq_len(nrow(df_trait_sp_summ))[-id_keep]
  max_second <- max(c(df_trait_sp_summ$mean[id_keep]), na.rm = TRUE) # Specify max of second y axis
  min_second <- min(c(df_trait_sp_summ$mean[id_keep]), na.rm = TRUE) # Specify min of second y axis
  if (any(is.infinite(max_first)))  max_first  <- max(df_trait_sp_summ$mean[-id_keep][-which(is.infinite(df_trait_sp_summ$mean[-id_keep]))], na.rm = TRUE)   # Specify max of first y axis
  if (any(is.infinite(max_second)))  max_second  <- max(df_trait_sp_summ$mean[id_keep][-which(is.infinite(df_trait_sp_summ$mean[id_keep]))], na.rm = TRUE)   # Specify max of first y axis
  if (any(is.infinite(min_first)))  min_first  <- min(df_trait_sp_summ$mean[-id_keep][-which(is.infinite(df_trait_sp_summ$mean[-id_keep]))], na.rm = TRUE)   # Specify max of first y axis
  if (any(is.infinite(min_second)))  min_second  <- min(df_trait_sp_summ$mean[id_keep][-which(is.infinite(df_trait_sp_summ$mean[id_keep]))], na.rm = TRUE)   # Specify max of first y axis
  
  
  
  
  z_axis <- ifelse(length(var_scale_regular) > 0, TRUE, FALSE)
  if(z_axis)  {z_axis <- paste(var_names[-which(var_names %in% var_scale_regular)], collapse = "\n")
  } else z_axis <- ""
  if (nchar(as.character(z_axis[1]))>1) {
    name_y <- paste(var_scale_regular, collapse="\n")
  }else  name_y <- paste(mean_var[,1] , collapse="\n")
  legend_labels <- unique(c(as.character(df_trait_sp_summ$variable), as.character(df_trait_sp_summ$Fmax)))
  legend_pos <- resolve_legend_position(
    position_legend = position_legend,
    df = df_trait_sp_summ,
    x_col = "age",
    y_col = "mean",
    legend_labels = legend_labels,
    y_axis_name = name_y,
    z_axis_name = z_axis
  )

  if (is.na(chosen_title)) title = paste0(paste(var_names, collapse = " and "), " mean\nover ", Time1, " and ", Time2)
  if (!is.na(chosen_title))title = chosen_title
  # Build the plot with ribbons for variability.
  if (length(table(df_trait_sp_summ$variable))==1){
    p <- ggplot(df_trait_sp_summ, aes(x = age, group = interaction(variable, Fmax))) +
      # geom_ribbon(aes(ymin = mean - sd, ymax = mean + sd),
      #             alpha = 0.5, fill = "grey70"
      # ) +
      # geom_line(method = "loess", se = TRUE,
      #             aes(group = interaction(variable, Fmax), y=mean)
      # ) +
      # geom_point    
      geom_line(aes(y=mean), lwd=0.5) +
      geom_point(aes(pch = Fmax, y = mean)) +
      theme_bw() +
      labs(
        title = title,
        x = "age",
        y = ifelse(relative_to_ref, paste(name_y, "- RELATIVE", collapse = "\n"), 
                   ifelse(proportional, paste(name_y, " - PROPORTIONAL", collapse = "\n"), paste(name_y, collapse = "\n"))),
        subtitle = sp_name
      ) +
      theme(axis.text.x = element_text(angle = 60, hjust = 1),
        legend.position = "inside",
        legend.position.inside = legend_pos,
            legend.justification = c("right", "top"), 
        legend.title = element_text( size=8), legend.text=element_text(size=8), 
        panel.background = element_rect(fill="white"), 
        panel.border = element_rect(fill="black")) +
      scale_y_continuous(
        # limits = limits,
        sec.axis = sec_axis(~scale_function(., scale_factor, shift), #. * scale, # (. - a)/b,  #sec_axis(~. * scale, 
                            name = ifelse(relative_to_ref, paste(z_axis, "- RELATIVE", collapse = "\n"), 
                                          ifelse(proportional, paste(z_axis, " - PROPORTIONAL", collapse = "\n"), paste(z_axis, collapse = "\n")))
                            
        )
      )
  } else {
    p <- ggplot(df_trait_sp_summ, aes(x = age, group = interaction(variable, Fmax))) +
      # geom_ribbon(aes(ymin = mean - sd, ymax = mean + sd),
      #             alpha = 0.5, fill = "grey70"
      # ) +
      # geom_line(method = "loess", se = TRUE,
      #             aes(group = interaction(variable, Fmax), y=mean)
      # ) +
      # geom_point 
      theme_bw() +
      theme(axis.text.x = element_text(angle = 60, hjust = 1),
        legend.position = "inside",
        legend.position.inside = legend_pos,
            legend.justification = c("right", "top"), 
        legend.title = element_text( size=8), legend.text=element_text(size=8)) +
      geom_line(aes(y=mean, colour = variable), lwd=0.5) +
      geom_point(aes(pch = Fmax, y = mean)) +
      labs(
        title = title,
        x = "age",
        y = ifelse(relative_to_ref, paste(name_y, "- RELATIVE", collapse = "\n"), 
                   ifelse(proportional, paste(name_y, " - PROPORTIONAL", collapse = "\n"), paste(name_y, collapse = "\n"))),
        subtitle = sp_name
      ) +
      scale_y_continuous(
        # limits = limits,
        sec.axis = sec_axis(~scale_function(., scale_factor, shift), #. * scale, # (. - a)/b,  #sec_axis(~. * scale, 
                            name = ifelse(relative_to_ref, paste(z_axis, "- RELATIVE", collapse = "\n"), 
                                          ifelse(proportional, paste(z_axis, " - PROPORTIONAL", collapse = "\n"), paste(z_axis, collapse = "\n")))
        )
      )
  }
  if (facet==T & !("SSDg"%in%colnames(df_trait_sp))){
    p <- p+
      facet_wrap(~variable, nrow = nrow, ncol=ncol,drop=FALSE, scale=scale_facet)+
      theme(legend.position = c(1, 0),
        legend.justification = c(1, 0), 
        legend.direction = "vertical", 
        legend.box = "horizontal")
  }
  if ("SSDg"%in%colnames(df_trait_sp)){
    # df_trait_sp$Fmax <- as.factor(df_trait_sp$Fmax)
    # inset.plot <- ggplot(df_trait_sp, group = (variable)) +
    #   # theme(legend.position = "none") +
    #   geom_boxplot(aes(y=SSDg, x=Fmax)) +
    #   facet_wrap(~variable) +
    #   theme_bw()
    p <- p+
      facet_wrap(~variable, nrow = nrow, ncol=ncol, drop=FALSE, strip.position = "top", scale=scale_facet) +
      scale_x_continuous(sec.axis = dup_axis())+
      theme(legend.justification = c(1, 0),
            legend.direction = "vertical",
            legend.box = "horizontal", 
            strip.placement = "outside")
    # p <- ggdraw() +
    #   draw_plot(p) +
    #   draw_plot(inset.plot, x = 0.63, y = 0.0, width = .3, height = .32)
  }
  
  
  # Return the ggplot object for downstream customization.
  return(p)
}

# estimate the T/Vk ratio 
estimateT_Vk <- function(trait_data_full_i){
  vec_i <- c(grep(trait_data_full_i$variable, pattern="generation"), 
            grep(trait_data_full_i$variable, pattern="Vk"))
  df_trait_t <- trait_data_full_i[vec_i,]

  # df_trait_t_rep <- df_trait_t %>%
  #   dplyr::group_by(replicate, Fstart, Fmax, Species, variable) %>%
  #   dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
  #   dplyr::filter(n > 1L)

  df_trait_t <- pivot_wider(df_trait_t, names_from = "variable", values_from = "value") 
  df_trait_t[, "T/Vk."] <- df_trait_t[, grep(pattern = "generation", colnames(df_trait_t))]/df_trait_t[, grep(pattern = "V", colnames(df_trait_t))]
  df_trait_t[which(is.na(as.matrix(df_trait_t[, "T/Vk."]))), "T/Vk."] <- 0
  add = "T/Vk."
  if ("cohort" %in% colnames(trait_data_full_i)) add <- c(add, "cohort")
  if ("age" %in% colnames(trait_data_full_i)) add <- c(add, "age")
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

  # df_trait_t_rep <- df_trait_t %>%
  #   dplyr::group_by(replicate, Fstart, Fmax, Species, variable) %>%
  #   dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
  #   dplyr::filter(n > 1L)

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
  
  # df_trait_t_rep <- df_trait_t %>%
  #   dplyr::group_by(replicate, Fstart, Fmax, Species, variable) %>%
  #   dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
  #   dplyr::filter(n > 1L)
  
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
estimate_BinaryOp <- function(trait_data_full_i, var, operator = "/") {
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
  id_cols <- c("Species", "replicate", "Fstart", "Fmax")
  if ("age" %in% colnames(trait_data_full_i)) id_cols <- c(id_cols, "age")
  if ("cohort" %in% colnames(trait_data_full_i)) id_cols <- c(id_cols, "cohort")

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
estimateDiv <- function(trait_data_full_i, var) {
  estimate_BinaryOp(trait_data_full_i, var, operator = "/")
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
  df_trait_t <- pivot_longer(df_trait_t[, c("Species", "replicate", "Fstart", "Fmax", "age", var)], cols = var, values_to = "value", names_to = "variable")
  
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
  var_i <- unlist(stringr::str_split(pattern = "\nsum", var))
  vec_i <- as.vector(unlist(sapply(var_i, FUN = function(x) (which(trait_data_full_i$variable== x)))))
  df_trait_t <- trait_data_full_i[vec_i,]
  
  df_trait_t <- pivot_wider(df_trait_t, names_from = "variable", values_from = "value") 
  if (soustraction == T) df_trait_t[, which(colnames(df_trait_t)== var_i[2])] <- -df_trait_t[, which(colnames(df_trait_t)== var_i[2])]
  df_trait_t[, var] <- df_trait_t[, which(colnames(df_trait_t)== var_i[1])]+df_trait_t[, which(colnames(df_trait_t)== var_i[2])]
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
estimatecolSum <- function(trait_data_full_i, var){
  var_i <- stringr::str_remove(pattern = "colSum", var)[1]
  vec_i <- grep(trait_data_full_i$variable, pattern= var_i)
  df_trait_t <- trait_data_full_i[vec_i,]
  if (length(unique(df_trait_t$variable))>1)  df_trait_t = df_trait_t[which(df_trait_t$variable == var_i),]
  
  df_trait_t$age <- as.numeric(as.matrix(df_trait_t$age))
  df_trait_t <- df_trait_t[order(df_trait_t$age), ]
  add = var[1]
  if (any(grepl("cohort", colnames(trait_data_full_i)))) {
    add <- c(add, "cohort")
    df_trait_t <- df_trait_t%>%
      group_by(replicate, Species, Fstart, Fmax, cohort) %>%
      summarize(var = sum(value, na.rm = T)) %>%
      ungroup()
  } 
  if (any(grepl("Scenario", colnames(trait_data_full_i)))) {
    add <- c(add, "Scenario")
    df_trait_t <- df_trait_t%>%
      group_by(replicate, Species, Fstart, Fmax, Scenario) %>%
      summarize(var = sum(value, na.rm = T)) %>%
      ungroup()
    
  }
  if (!any(grepl("cohort", colnames(trait_data_full_i))) & !any(grepl("Scenario", colnames(trait_data_full_i)))){
    df_trait_t <- df_trait_t%>%
      group_by(replicate, Species, Fstart, Fmax) %>%
      summarize(var = sum(value, na.rm = T))
  }
  colnames(df_trait_t)[which(colnames(df_trait_t) == "var")] <- var[1]
  data_end <- df_trait_t
  # data_end <- NULL
  # for (age in as.numeric(names(table(trait_data_full_i$age)))) {
  #   df_trait_t$age <- age
  #   data_end <- rbind(data_end, df_trait_t)
  # } "age",
  df_trait_t <- pivot_longer(data_end[, c("Species", "replicate", "Fstart", "Fmax", add)], cols = var, values_to = "value", names_to = "variable")
  
  return(df_trait_t)
}
# estimate age 1
estimateage1 <- function(trait_data_full_i, var){
  var_i <- stringr::str_remove(pattern = "age1", var)[1]
  vec_i <- grep(trait_data_full_i$variable, pattern= var_i)
  df_trait_t <- trait_data_full_i[vec_i,]
  # if (length(unique(df_trait_t$variable))>1)  {
  #   df_trait_t = df_trait_t[which(df_trait_t$variable == unique(df_trait_t$variable)[1])]
  # }  
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








#### FOR PCA --------------------------------

# pca axis built from the complete dataset with all species, plot species one by one
plot_PCA <- function (data, DIFF = FALSE, Traits_raw = FALSE, Traits_Ne = FALSE, pathoutputplots){
  if (!any(grepl(names(data), pattern="Ne"))) {
    tempo <- c(names(data[[1]])[1], names(data[[1]])[2], names(data[[1]])[3])
    data <- list(do.call(rbind, list(data[[1]][[names(data[[1]])[1]]], data[[2]][[names(data[[1]])[1]]])), 
                                  do.call(rbind, list(data[[1]][[names(data[[1]])[2]]], data[[2]][[names(data[[1]])[2]]])),
                                  do.call(rbind, list(data[[1]][[names(data[[1]])[3]]], data[[2]][[names(data[[1]])[3]]])))
    names(data) <- tempo

  }
  
  for (j in names(data)){
    
    # store dataset and add dimensions if needed
    P <- list()
    dat_merged_w_sc <- data[[j]]
    if (j=="data_Nb" & !any(grepl(colnames(dat_merged_w_sc), pattern="genera"))) dat_merged_w_sc <- rbind(dat_merged_w_sc,  data[["tempoNe"]][-grep(pattern= "tempo", data[["tempoNe"]]$variable),])
    if (j=="tempoNe" & !any(grepl(colnames(dat_merged_w_sc), pattern="V"))) dat_merged_w_sc <- rbind(dat_merged_w_sc,  data[["data_Nb"]][-grep(pattern= "Nb", data[["data_Nb"]]$variable),])
    if ("variable" %in% colnames(dat_merged_w_sc)) dat_merged_w_sc <- pivot_wider(dat_merged_w_sc, names_from = "variable", values_from = "value")
    
    # ACP
    dat_merged_pca <- dat_merged_w_sc
    dat_merged_pca$Fmax <- as.numeric(dat_merged_pca$Fmax)
    dat_merged_pca$Fstart <- as.numeric(dat_merged_pca$Fstart)
    dat_merged_pca$replicate <- as.numeric(dat_merged_pca$replicate)
    
    if (Traits_raw == TRUE) {
      ID_keep_dat <- c(grep(colnames(dat_merged_pca), pattern="cohort"), 
                       grep(colnames(dat_merged_pca), pattern="year"), 
                       which(colnames(dat_merged_pca) %in% c("replicate", "Fstart", "Fmax", "Species")), 
                       which(colnames(dat_merged_pca) %in% c("mean_Ne_N_alpha", "mean_ratio_Ne_N_alpha", 
                                                             "mean_Nb", "mean_ratio_Nb", "htz_slope", "mean_Ne_temporal", "mean_ratio_Ne_temporal","Nb/NT", "Nb",
                                                             "Ne_temporal/NT", "Ne_temporal", "Ne_N_alpha", "Ne_N_alpha/Numerator_N_alpha", 
                                                             "diff_mean_Ne_N_alpha", "diff_mean_ratio_Ne_N_alpha", 
                                                             "diff_mean_Nb", "diff_mean_ratio_Nb", "diff_htz_slope", "diff_mean_Ne_temporal", "diff_mean_ratio_Ne_temporal"))
      )
      dat_merged_pca <- dat_merged_pca[, ID_keep_dat]
    } else if (Traits_Ne == TRUE) {
      ID_drop_dat <- c(grep(colnames(dat_merged_pca), pattern="cohort"),
                       grep(colnames(dat_merged_pca), pattern="year"))
      if (length(ID_drop_dat) > 0) {
        dat_merged_pca <- dat_merged_pca[, -ID_drop_dat, drop = FALSE]
      }
    }
    
    
    colnm_div_genet <- which(colnames(dat_merged_pca) %in% c("mean_Ne_N_alpha", "mean_ratio_Ne_N_alpha",
                                                             "mean_Nb", "mean_ratio_Nb", "htz_slope", "mean_Ne_temporal", "mean_ratio_Ne_temporal","Nb/NT", "Nb",
                                                             "Ne_temporal/NT", "Ne_temporal", "Ne_N_alpha", "Ne_N_alpha/Numerator_N_alpha", 
                                                             "diff_mean_Ne_N_alpha", "diff_mean_ratio_Ne_N_alpha", 
                                                             "diff_mean_Nb", "diff_mean_ratio_Nb", "diff_htz_slope", "diff_mean_Ne_temporal", "diff_mean_ratio_Ne_temporal"))
    for (div_genet in colnm_div_genet){
      
      
      data_PCA <- dat_merged_pca
      data_PCA <- data_PCA[,-colnm_div_genet]
      assign(colnames(dat_merged_pca)[div_genet], dat_merged_pca[div_genet])
      
      pca1 <- PCA(data_PCA |> dplyr::select(where(is.numeric)),
                  scale.unit = TRUE,
                  graph = FALSE)
      
      for (spi in unique(data_PCA$Species)){
        ind <- as.data.frame(pca1$ind$coord)
        ind$Species <- dat_merged_pca$Species
        ind[colnames(dat_merged_pca)[div_genet]] <- get(colnames(dat_merged_pca)[div_genet])
        ind[colnames(dat_merged_pca)[div_genet]][which(is.na(ind[colnames(dat_merged_pca)[div_genet]])),] <- 0
        if(any(ind[colnames(dat_merged_pca)[div_genet]]==0)) ind <- ind[-which(ind[colnames(dat_merged_pca)[div_genet]]==0),]
        ind <- ind[which(ind$Species %in% spi),]
        
        Var <- as.data.frame(pca1$var$coord)
        scale <- mean(unlist(as.vector(Var[,1:2])), na.rm=T)/quantile(unlist(as.vector(ind[,1:2])), na.rm=T, 0.70)
        if (scale>1) scale = 1/(2*scale)
        if (abs(scale) < 0.05) scale = 3*scale
        Var <- Var/scale
        Var$var <- rownames(Var)
        
        P[[spi]] <- ggplot() +  
          stat_summary_hex(fun = mean, bins = 25, data=ind, 
                           aes(Dim.1, Dim.2, z=.data[[colnames(dat_merged_pca)[div_genet]]])) +
          geom_point(data=ind, 
                     aes(Dim.1, Dim.2, color = .data[[colnames(dat_merged_pca)[div_genet]]]), alpha=0.5) +
          
          geom_segment(data = Var, aes(x=0, y=0, xend = Dim.1, yend = Dim.2), arrow = arrow(length = unit(0.1, "inches"))
          )+  
          geom_text(data = Var, aes(x = Dim.1, y = Dim.2, label = var), size=3)+
          scale_color_viridis_c() +
          theme_minimal() +
          geom_abline(slope = 0, intercept = 0, linetype = "dashed", color = "grey") +
          geom_hline(yintercept = 0, linetype = "dashed", color = "grey") +
          labs(title = paste0("PCA of traits vs ", colnames(dat_merged_pca)[div_genet]),
               subtitle = spi,
               x = paste0("PC1, ", round(pca1$eig[1,2], 2), "%"),
               y = paste0("PC2, ", round(pca1$eig[2,2], 2), "%"))
        
      }
      if(!dir.exists(pathoutputplots)) dir.create(pathoutputplots, recursive = T)
      var_names <- stringr::str_replace(colnames(dat_merged_pca)[div_genet], pattern = "/", "ratio")
      name_pdf = paste0(pathoutputplots, "/PCA", 
                        ifelse(grepl(var_names, pattern="htz"), paste0(j, var_names), var_names),
                        ifelse(Traits_raw==TRUE, "_raw_traits", ""),ifelse(Traits_Ne==TRUE, "_Ne_traits", ""),ifelse(DIFF==TRUE, "_DIFF", ""), 
                        "_one_genet_drift_PCA_htz_vs_traits.pdf")
      pdf(name_pdf, width = 10, height = 8)
      print(P)
      dev.off()
    }
  }
}
# pca axis built from the 1 species dataset at a time, plot species one by one
plot_PCA_sep_species <- function (data, DIFF = FALSE, Traits_raw = FALSE, Traits_Ne = FALSE, pathoutputplots){
  if (!any(grepl(names(data), pattern="Ne"))) {
    tempo <- c(names(data[[1]])[1], names(data[[1]])[2], names(data[[1]])[3])
    data <- list(do.call(rbind, list(data[[1]][[names(data[[1]])[1]]], data[[2]][[names(data[[1]])[1]]])), 
                 do.call(rbind, list(data[[1]][[names(data[[1]])[2]]], data[[2]][[names(data[[1]])[2]]])),
                 do.call(rbind, list(data[[1]][[names(data[[1]])[3]]], data[[2]][[names(data[[1]])[3]]])))
    names(data) <- tempo
    
  }
  for (j in names(data)){
    P <- list()
    dat_merged_w_sc <- data[[j]]
    if (j=="data_Nb" & !any(grepl(colnames(dat_merged_w_sc), pattern="genera"))) dat_merged_w_sc <- rbind(dat_merged_w_sc,  data[["tempoNe"]][-grep(pattern= "tempo", data[["tempoNe"]]$variable),])
    if (j=="tempoNe" & !any(grepl(colnames(dat_merged_w_sc), pattern="V"))) dat_merged_w_sc <- rbind(dat_merged_w_sc,  data[["data_Nb"]][-grep(pattern= "Nb", data[["data_Nb"]]$variable),])
    if ("variable" %in% colnames(dat_merged_w_sc)) dat_merged_w_sc <- pivot_wider(dat_merged_w_sc, names_from = "variable", values_from = "value")
    # ACP
    dat_merged_pca <- dat_merged_w_sc
    dat_merged_pca$Fmax <- as.numeric(dat_merged_pca$Fmax)
    dat_merged_pca$Fstart <- as.numeric(dat_merged_pca$Fstart)
    dat_merged_pca$replicate <- as.numeric(dat_merged_pca$replicate)
    
    if (Traits_raw == TRUE) {
      ID_keep_dat <- c(grep(colnames(dat_merged_pca), pattern="cohort"), 
                       grep(colnames(dat_merged_pca), pattern="year"), 
                       which(colnames(dat_merged_pca) %in% c("replicate", "Fstart", "Fmax", "Species")), 
                       which(colnames(dat_merged_pca) %in% c("mean_Ne_N_alpha", "mean_ratio_Ne_N_alpha", 
                                                             "mean_Nb", "mean_ratio_Nb", "htz_slope", "mean_Ne_temporal", "mean_ratio_Ne_temporal","Nb/NT", "Nb",
                                                             "Ne_temporal/NT", "Ne_temporal", "Ne_N_alpha", "Ne_N_alpha/Numerator_N_alpha", 
                                                             "diff_mean_Ne_N_alpha", "diff_mean_ratio_Ne_N_alpha", 
                                                             "diff_mean_Nb", "diff_mean_ratio_Nb", "diff_htz_slope", "diff_mean_Ne_temporal", "diff_mean_ratio_Ne_temporal"))
      )
      dat_merged_pca <- dat_merged_pca[, ID_keep_dat]
    } else if (Traits_Ne == TRUE) {
      ID_drop_dat <- c(grep(colnames(dat_merged_pca), pattern="cohort"),
                       grep(colnames(dat_merged_pca), pattern="year"))
      if (length(ID_drop_dat) > 0) {
        dat_merged_pca <- dat_merged_pca[, -ID_drop_dat, drop = FALSE]
      }
    }
    
    
    colnm_div_genet <- which(colnames(dat_merged_pca) %in% c("mean_Ne_N_alpha", "mean_ratio_Ne_N_alpha",
                                                             "mean_Nb", "mean_ratio_Nb", "htz_slope", "mean_Ne_temporal", "mean_ratio_Ne_temporal","Nb/NT", "Nb",
                                                             "Ne_temporal/NT", "Ne_temporal", "Ne_N_alpha", "Ne_N_alpha/Numerator_N_alpha", 
                                                             "diff_mean_Ne_N_alpha", "diff_mean_ratio_Ne_N_alpha", 
                                                             "diff_mean_Nb", "diff_mean_ratio_Nb", "diff_htz_slope", "diff_mean_Ne_temporal", "diff_mean_ratio_Ne_temporal"))
    for (div_genet in colnm_div_genet){
      
      
      data_PCA <- dat_merged_pca
      data_PCA <- data_PCA[,-colnm_div_genet]
      assign(colnames(dat_merged_pca)[div_genet], dat_merged_pca[div_genet])
      
      
      for (spi in unique(data_PCA$Species)){
        
        data_PCA_sp <- data_PCA[which(data_PCA$Species == spi), ]
        pca1 <- PCA(data_PCA_sp |> dplyr::select(where(is.numeric)),
                    scale.unit = TRUE,
                    graph = FALSE)
        
        ind <- as.data.frame(pca1$ind$coord)
        ind$Species <- data_PCA_sp$Species
        ind[colnames(dat_merged_pca)[div_genet]] <- get(colnames(dat_merged_pca)[div_genet])[which(dat_merged_pca$Species == spi),1]
        if (any((is.na(ind[colnames(dat_merged_pca)[div_genet]])))) ind[colnames(dat_merged_pca)[div_genet]][which(is.na(ind[colnames(dat_merged_pca)[div_genet]])),] <- 0
        if (any((ind[colnames(dat_merged_pca)[div_genet]]==0)))     ind <- ind[-which(ind[colnames(dat_merged_pca)[div_genet]]==0),]
        
        Var <- as.data.frame(pca1$var$coord)
        scale <- mean(unlist(as.vector(Var[,1:2])), na.rm=T)/quantile(unlist(as.vector(ind[,1:2])), na.rm=T, 0.70)
        if (scale>1) scale = 1/(2*scale)
        if (abs(scale) < 0.05) scale = 3*scale
        Var <- Var/scale
        Var$var <- rownames(Var)
        
        P[[spi]] <- ggplot() +  
          stat_summary_hex(fun = mean, bins = 25, data=ind, 
                           aes(Dim.1, Dim.2, z=.data[[colnames(dat_merged_pca)[div_genet]]])) +
          geom_point(data=ind, 
                     aes(Dim.1, Dim.2, color = .data[[colnames(dat_merged_pca)[div_genet]]]), alpha=0.5) + #
          
          geom_segment(data = Var, aes(x=0, y=0, xend = Dim.1, yend = Dim.2), arrow = arrow(length = unit(0.1, "inches"))
          )+  
          geom_text(data = Var, aes(x = Dim.1, y = Dim.2, label = var), size=3)+
          scale_color_viridis_c() +
          theme_minimal() +
          geom_abline(slope = 0, intercept = 0, linetype = "dashed", color = "grey") +
          geom_hline(yintercept = 0, linetype = "dashed", color = "grey") +
          labs(title = paste0("PCA of traits vs ", colnames(dat_merged_pca)[div_genet]),
               subtitle = spi,
               x = paste0("PC1, ", round(pca1$eig[1,2], 2), "%"),
               y = paste0("PC2, ", round(pca1$eig[2,2], 2), "%"))
        
      }
      if(!dir.exists(pathoutputplots)) dir.create(pathoutputplots, recursive = T)
      var_names <- stringr::str_replace(colnames(dat_merged_pca)[div_genet], pattern = "/", "ratio")
      name_pdf = paste0(pathoutputplots, "/PCA_SEP_SP_", 
                        ifelse(grepl(var_names, pattern="htz"), paste0(j, var_names), var_names),
                        ifelse(Traits_raw==TRUE, "_raw_traits", ""),ifelse(Traits_Ne==TRUE, "_Ne_traits", ""),ifelse(DIFF==TRUE, "_DIFF", ""), 
                        "_one_genet_drift_PCA_htz_vs_traits.pdf")
      pdf(name_pdf, width = 10, height = 8)
      print(P)
      dev.off()
    }
  }
}
# ALL TRAITS, all dimensions, all species
plot_PCA_ALL <- function (data, DIFF = FALSE, Traits_raw = FALSE, Traits_Ne = FALSE, pathoutputplots){
  if (!any(grepl(names(data), pattern="Ne"))) {
    tempo <- c(names(data[[1]])[1], names(data[[1]])[2], names(data[[1]])[3])
    data <- list(do.call(rbind, list(data[[1]][[names(data[[1]])[1]]], data[[2]][[names(data[[1]])[1]]])), 
                 do.call(rbind, list(data[[1]][[names(data[[1]])[2]]], data[[2]][[names(data[[1]])[2]]])),
                 do.call(rbind, list(data[[1]][[names(data[[1]])[3]]], data[[2]][[names(data[[1]])[3]]])))
    names(data) <- tempo
    
  }
  for (j in names(data)){
    P <- list()
    dat_merged_w_sc <- data[[j]]
    if (j=="data_Nb" & !any(grepl(colnames(dat_merged_w_sc), pattern="genera"))) dat_merged_w_sc <- rbind(dat_merged_w_sc,  data[["tempoNe"]][-grep(pattern= "tempo", data[["tempoNe"]]$variable),])
    if (j=="tempoNe" & !any(grepl(colnames(dat_merged_w_sc), pattern="V"))) dat_merged_w_sc <- rbind(dat_merged_w_sc,  data[["data_Nb"]][-grep(pattern= "Nb", data[["data_Nb"]]$variable),])
    if ("variable" %in% colnames(dat_merged_w_sc)) dat_merged_w_sc <- pivot_wider(dat_merged_w_sc, names_from = "variable", values_from = "value")
    # ACP
    dat_merged_pca <- dat_merged_w_sc
    dat_merged_pca$Fmax <- as.numeric(dat_merged_pca$Fmax)
    dat_merged_pca$Fstart <- as.numeric(dat_merged_pca$Fstart)
    dat_merged_pca$replicate <- as.numeric(dat_merged_pca$replicate)
    
    if (Traits_raw == TRUE) {
      ID_keep_dat <- c(grep(colnames(dat_merged_pca), pattern="cohort"),
                       grep(colnames(dat_merged_pca), pattern="year"),
                       which(colnames(dat_merged_pca) %in% c("replicate", "Fstart", "Fmax", "Species")),
                       which(colnames(dat_merged_pca) %in% c("mean_Ne_N_alpha", "mean_ratio_Ne_N_alpha",
                                                             "mean_Nb", "mean_ratio_Nb", "htz_slope", "mean_Ne_temporal", "mean_ratio_Ne_temporal","Nb/NT", "Nb",
                                                             "Ne_temporal/NT", "Ne_temporal", "Ne_N_alpha", "Ne_N_alpha/Numerator_N_alpha",  
                                                             "diff_mean_Ne_N_alpha", "diff_mean_ratio_Ne_N_alpha", 
                                                             "diff_mean_Nb", "diff_mean_ratio_Nb", "diff_htz_slope", "diff_mean_Ne_temporal", "diff_mean_ratio_Ne_temporal"))
      )
      dat_merged_pca <- dat_merged_pca[, ID_keep_dat]
    } else if (Traits_Ne == TRUE) {
      ID_drop_dat <- c(grep(colnames(dat_merged_pca), pattern="cohort"),
                       grep(colnames(dat_merged_pca), pattern="year"))
      if (length(ID_drop_dat) > 0) {
        dat_merged_pca <- dat_merged_pca[, -ID_drop_dat, drop = FALSE]
      }
    }
    
    
    colnm_div_genet <- which(colnames(dat_merged_pca) %in% c("mean_Ne_N_alpha", "mean_ratio_Ne_N_alpha",
                                                             "mean_Nb", "mean_ratio_Nb", "htz_slope", "mean_Ne_temporal", "mean_ratio_Ne_temporal","Nb/NT", "Nb",
                                                             "Ne_temporal/NT", "Ne_temporal", "Ne_N_alpha", "Ne_N_alpha/Numerator_N_alpha",
                                                             "diff_mean_Ne_N_alpha", "diff_mean_ratio_Ne_N_alpha",
                                                             "diff_mean_Nb", "diff_mean_ratio_Nb", "diff_htz_slope", "diff_mean_Ne_temporal", "diff_mean_ratio_Ne_temporal"))
    for (div_genet in colnm_div_genet){
      
      
      data_PCA <- dat_merged_pca
      data_PCA <- data_PCA[,-colnm_div_genet]
      assign(colnames(dat_merged_pca)[div_genet], dat_merged_pca[div_genet])
      
      
      # for (spi in unique(data_PCA$Species)){
      
      data_PCA_sp <- data_PCA#[which(data_PCA$Species == spi), ]
      pca1 <- PCA(data_PCA_sp |> dplyr::select(where(is.numeric)),
                  scale.unit = TRUE,
                  graph = FALSE)
      
      ind <- as.data.frame(pca1$ind$coord)
      ind$Species <- data_PCA_sp$Species
      ind[colnames(dat_merged_pca)[div_genet]] <- get(colnames(dat_merged_pca)[div_genet])#[which(dat_merged_pca$Species == spi),1]
      if (any((is.na(ind[colnames(dat_merged_pca)[div_genet]])))) ind[colnames(dat_merged_pca)[div_genet]][which(is.na(ind[colnames(dat_merged_pca)[div_genet]])),] <- 0
      if (any((ind[colnames(dat_merged_pca)[div_genet]]==0)))     ind <- ind[-which(ind[colnames(dat_merged_pca)[div_genet]]==0),]
      
      Var <- as.data.frame(pca1$var$coord)
      scale <- mean(unlist(as.vector(Var[,1:2])), na.rm=T)/quantile(unlist(as.vector(ind[,1:2])), na.rm=T, 0.70)
      if (scale>1) scale = 1/(2*scale)
      Var <- Var/scale
      Var$var <- rownames(Var)
      
      #[[spi]]
      P<- ggplot() +
        stat_summary_hex(fun = mean, bins = 25, data=ind,
                         aes(Dim.1, Dim.2, z=.data[[colnames(dat_merged_pca)[div_genet]]])) +
        geom_point(data=ind,
                   aes(Dim.1, Dim.2, color = Species), alpha=0.5) + #.data[[colnames(dat_merged_pca)[div_genet]]]
        
        geom_segment(data = Var, aes(x=0, y=0, xend = Dim.1, yend = Dim.2), arrow = arrow(length = unit(0.1, "inches"))
        )+
        geom_text(data = Var, aes(x = Dim.1, y = Dim.2, label = var), size=3)+
        scale_color_viridis_d() +
        theme_minimal() +
        geom_abline(slope = 0, intercept = 0, linetype = "dashed", color = "grey") +
        geom_hline(yintercept = 0, linetype = "dashed", color = "grey") +
        labs(title = paste0("PCA of traits vs ", colnames(dat_merged_pca)[div_genet]),
             #subtitle = spi,
             x = paste0("PC1, ", round(pca1$eig[1,2], 2), "%"),
             y = paste0("PC2, ", round(pca1$eig[2,2], 2), "%"))
      
      #}
      if(!dir.exists(pathoutputplots)) dir.create(pathoutputplots, recursive = T)
      var_names <- stringr::str_replace(colnames(dat_merged_pca)[div_genet], pattern = "/", "ratio")
      name_pdf = paste0(pathoutputplots, "/PCA_ALL_", 
                        ifelse(grepl(var_names, pattern="htz"), paste0(j, var_names), var_names),
                        ifelse(Traits_raw==TRUE, "_raw_traits", ""),ifelse(Traits_Ne==TRUE, "_Ne_traits", ""),ifelse(DIFF==TRUE, "_DIFF", ""),
                        # "_one_genet_drift_PCA_htz_vs_traits.pdf"
                        "ALL.pdf")
      pdf(name_pdf, width = 10, height = 8)
      print(P)
      dev.off()
    }
  }
}




#### FOR MORTALITY as TL--------------------------------

source(paste0(here::here(), "/R/functions/function_data_traits_summary.R"))

prepare_only_simutraits_dataset <- function(path_simutraits, survival =F, natural_mort =F, F1=FALSE, FORCE=F){
  data_complete_y = NULL
  
  # build dataset by combining multiple RData
  folder_name <- basename(path_simutraits)
  suffix <- gsub("Trait_variables_by_", "", folder_name)
  suffix <- gsub("_", "_", suffix)  # keep underscores
  save_dir  <- file.path(dirname(dirname(path_simutraits)), "RDATA")
  save_file <- file.path(save_dir, paste0("summary_trait_", suffix, ".RData"))
  if (!dir.exists(save_dir)){dir.create(save_dir, recursive = T)}
  if (!file.exists(save_file) | FORCE == TRUE){
    data_complete_y <- data_traits_summary(base_path = path_simutraits)
  } else{ 
    data_complete_y = get(load(file=save_file))
  }
  
  # dimensions dataset required : [time, Sp, value, rep, Fish_press, TRAIT]
  head(data_complete_y)
  unique(data_complete_y$TRAIT)
  
  # if needed : add certain outputs : 
  if (survival == TRUE){
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
  }
  
  if (natural_mort == TRUE){
    # adding natural mortality dimension
    data_complete_y[grep((data_complete_y$TRAIT), pattern= "Mnat"), "TRAIT"] <- "Madd"
    data_complete_y <- pivot_wider(data = data_complete_y, names_from = "TRAIT", values_from = "value")
    if(any(grep(colnames(data_complete_y), pattern="Dis"))) data_complete_y <- data_complete_y[,-grep(colnames(data_complete_y), pattern="Dis")]
    if(any(grep(colnames(data_complete_y), pattern="Z"))) data_complete_y <- data_complete_y[,-grep(colnames(data_complete_y), pattern="Z")]
    data_complete_y[,"Mnat"] <- data_complete_y[, grep(pattern="Madd", colnames(data_complete_y))] + 
                               data_complete_y[, grep(pattern="Mpred", colnames(data_complete_y))] + 
                               data_complete_y[, grep(pattern="Mstar", colnames(data_complete_y))]
    data_complete_y <- pivot_longer(data = data_complete_y, names_to = "TRAIT", values_to = "value", cols=-colnames(data_complete_y)[1:4])
  }
  
  if (F1 == TRUE){
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
  }
  
  return(data_complete_y)
}

