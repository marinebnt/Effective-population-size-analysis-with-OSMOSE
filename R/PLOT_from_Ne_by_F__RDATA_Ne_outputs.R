library(ggplot2)
library(gridExtra)   # pour arranger plusieurs plots
library(patchwork)   # alternative moderne (recommandée)
library(stringr)
library(dplyr)
library(GGally)


base_dir <- "osmose-ns-v4.4.1_clim/output/PLOTS/All_clean_genetic_drift_indicators/RDATA/"
folders  <- list.dirs(base_dir, recursive = FALSE)


plot_explain_F <- function(df_sp, y_var, species_name, ref_sp_name) {
  
  ggplot(df_sp, aes(x = Fmax, y = mean_global)) +
    geom_point(size = 3) +
    geom_errorbar(
      aes(
        ymin = mean_global - sd_global,
        ymax = mean_global + sd_global
      ),
      width = 0.05
    ) +
    facet_wrap(~ Fstart) +
    theme_bw()+
    labs(
      title = paste0(y_vars[1], " mean over 150 and 175, \n, ", species_name, " sp with variable F is ", ref_sp_name),
      x = expression(F[max]~"(fishing intensity)"),
      y = y_vars[1],
      subtitle = "Facets = fishing onset (Fstart)"
    )
}






for (folder in folders[1]) {
  
  rdata_files <- list.files(folder, pattern = "\\.RData$", full.names = TRUE)
  
  for (file in rdata_files[1]) {
    
    # Charger le fichier
    plot_obj <- get(load(file))
    
    if (length(plot_obj)==0){
      cat(file, " is empty\n")
      next
    }
    
    # Si c'est une liste → on la convertit en data.frame
    if (is.list(plot_obj)) {
      df <- do.call(
        rbind,
        Map(
          function(x, nm) {
            x$species <- nm
            x
          },
          plot_obj,
          names(plot_obj)
        )
      )
    } else {
      df <- plot_obj
    }

    # Variables numériques à tracer (hors métadonnées)
    y_vars <- setdiff(
      names(df),
      c("species", "years", "Fish_press", "rep", "Fstart", "Fmax")
    )
    ref_sp_name <- stringr::str_split(stringr::str_split(file, "_sp")[[1]][2], "/")[[1]][1]
    
    # prepare df to fit the requisits for the plot 
    df <- df %>%
      mutate(
        Fmax = as.numeric(str_extract(Fish_press, "(?<=Fmax=)[0-9.]+")),
        Fstart    = as.numeric(str_extract(Fish_press, "(?<=Fstart=)[0-9.]+"))
      )
    
    species <- unique(df$species)
    
    mean_var <- sym(y_vars[1])
    sd_var <- sym(y_vars[2])
    summary_df <- df %>%
      filter(years >= 150, years <= 175) %>%
      group_by(species, Fmax, Fstart) %>%
      summarise(
        mean_global = mean(!!mean_var),
        sd_global   = sqrt(!!mean_var) + mean((!!sd_var)^2),
        .groups = "drop"
      )
        
    # Nom du PDF
    pdf_name <- paste0(
      tools::file_path_sans_ext(basename(file)),
      ".pdf"
    )

    # enregistrer pdf
    pdf_path <- file.path(folder, pdf_name)
    pdf(pdf_path, width = 12, height = 8)
    
    for (sp in species) {
      df_sp <- summary_df[summary_df$species == sp, ]
      
      p <- plot_explain_F(df_sp, y_var, sp, ref_sp_name)
      print(p)
    }
    
    dev.off()
  }

}
