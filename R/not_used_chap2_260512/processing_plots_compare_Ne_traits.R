library(ggplot2)
library(gridExtra)   # pour arranger plusieurs plots
library(patchwork)   # alternative moderne (recommandée)
library(stringr)
library(dplyr)
library(GGally)
library(tidyverse)

main = "New_gpu_conf.R"
source("osmose-ns-v4.4.1_clim/R/functions/source.R")
base_dir <- "osmose-ns-v4.4.1_clim/output/PLOTS/All_genetic_drift_indicators/INIT/"
folders  <- list.dirs(base_dir, recursive = F)
INIT2 <- get(load(list.files(folders[1])))

folders_traits <- "osmose-ns-v4.4.1_clim/output/PLOTS/Trait_variables/RData/"
folders_traits  <- list.dirs(folders_traits, recursive = F)

# Auxiliar functions ------------------------------------------------------

myf = function(x) do.call(cbind, x)

myf2 = function(x) {
  for(i in seq_along(x)) {
    for(j in seq_along(x[[i]])) {
      check2 = any(sapply(x[[i]][[j]], is_tibble))
      if(check2) {
        for(k in seq_along(x[[i]][[j]])) {
          # print(x[[i]][[j]][[k]])
          x[[i]][[j]][[k]] = as.numeric(x[[i]][[j]][[k]][[2]])
          # print(x[[i]][[j]][[k]])
        }
      }
      check = any(sapply(x[[i]][[j]], is.matrix))
      if(isTRUE(check)) {
        print(sprintf("Skipping %s (%d) %s (%d)", names(x)[[i]], i, names(x[[i]])[j], j))
        next
      }
      tmp = try(myf(x[[i]][[j]]), silent = TRUE)
      if(!inherits(tmp, "try-error")) {
        print(sprintf("Merging %s (%d) %s (%d)", names(x)[[i]], i, names(x[[i]])[j], j))
        x[[i]][[j]] = tmp 
      }
    }
  }
  return(x)
}

reshape_list = function(x) {
  xxx = unlist(x, recursive = FALSE)
  nm = lapply(lapply(strsplit(names(xxx), split = ".", fixed = TRUE), FUN="[", i=c(2,1)), FUN=paste, collapse=".")
  names(xxx) = nm
  # xxx = xxx[sort(names(xxx))]
  
  sp = sapply(strsplit(names(xxx), split=".", fixed=TRUE), FUN = "[", i=1)
  var = sapply(strsplit(names(xxx), split=".", fixed=TRUE), FUN = "[", i=2)
  
  out = setNames(vector(length(unique(sp)), mode="list"), nm=unique(sp))
  for(i in seq_along(sp)) {
    out[[sp[i]]][[var[i]]] = xxx[[i]]
  }
  return(out)
}

reshape_data = function(x, varName="var_name", as.variable=FALSE) {
  if(!is.list(x)) return(x)
  nm = try(strsplit(names(x), split=".", fixed=TRUE), silent=TRUE)
  if(inherits(nm, "try-error")) return(x)
  ll = sapply(nm, length)
  exit = TRUE
  if(all(ll==1)) {
    dat = expand.grid(cohort=colnames(x[[1]]), replicate=unlist(nm))
    exit=FALSE
  }
  if(all(ll==2)) {
    dat = data.frame(cohort=sapply(nm, FUN = "[", i=2), 
                     replicate=sapply(nm, FUN = "[", i=1))
    exit=FALSE
  }
  if(isTRUE(exit)) return(x)
  if(isTRUE(as.variable)) {
    dat$variable = varName
    varName = "value"
  }
  dat[[varName]] = unlist(x)
  dat$cohort = as.numeric(gsub(dat$cohort, pattern="cohort_", replacement=""))
  dat$replicate = as.numeric(gsub(dat$replicate, pattern="rep_", replacement=""))
  return(dat)
}

reshape_all = function(x) {
  for(i in seq_along(x)) {
    for(j in seq_along(x[[i]])) {
      if(!is.list(x[[i]][[j]])) next
      x[[i]][[j]] = lapply(x[[i]][[j]], FUN = reshape_data)
      }
  }
  return(x)
}


# Processing one list -----------------------------------------------------

INIT2b = INIT2
INIT2b[2:3] = lapply(INIT2b[2:3], FUN = myf2)
INIT2b[1:3] = lapply(INIT2b[1:3], FUN = reshape_list)
INIT2b[[5]] = lapply(INIT2b[[5]], FUN = as.list)
INIT2b[[5]] = reshape_list(INIT2b[[5]])
INIT2b[[5]] = lapply(INIT2b[[5]], FUN = myf)
INIT2b[[5]]$age = NULL
INIT2b[[1]] = NULL
INIT2b = reshape_list(INIT2b)

longeH = reshape_data(INIT2b$AtlanticHerring$data_Ne$longevity, "longevity", as.variable = TRUE)
Ne_N1H = reshape_data(INIT2b$AtlanticHerring$data_Ne$Ne_N_1, "Ne_N1", as.variable = TRUE)
dat = rbind(longeH, Ne_N1H)
dat$fishing = 0

# Loop for Fishing Pressure -----------------------------------------------

folders
for (SP in species){
  dat_tot = NULL
  dat_plot = NULL
  p = list()
  ID_fold <- grep(SP, folders)
  
  trait_file <- list.files(folders_traits[grep(folders_traits, pattern=SP)], full.names = T)
  trait_var <- get(load(trait_file))
  
  mortalityH <- trait_var[which(trait_var$TRAIT == "mortality"), ]
  mortalityH <- mortalityH[which(mortalityH$Species == SP), ]
  mortalityH$Fmax <- unlist(stringr::str_split(unlist(stringr::str_split(mortalityH$Fish_press, "Fmax="))[seq(by=2, 2, 2*length(mortalityH$Fish_press))], " "))[seq(by=2, 1, 2*length(mortalityH$Fish_press))]
  mortalityH$Fstart <- unlist(stringr::str_split(mortalityH$Fish_press, "Fstart="))[seq(by=2, 2, 2*length(mortalityH$Fish_press))]
  mortalityH <- mortalityH[, which(colnames(mortalityH) %in% c("time", "mean", "Fmax", "Fstart"))]
  colnames(mortalityH) <- c("cohort", "value", "Fmax", "Fstart")
  mortalityH$variable = "mean_natural_mortality"
  mortalityH_l <-bind_rows(replicate(20,mortalityH, simplify=FALSE)) 
  mortalityH_l$replicate <- as.vector(sapply(1:20, FUN=function(x) rep(x, times=nrow(mortalityH))))

  
  for (f_sp in ID_fold){
    files <- list.files(folders[f_sp], full.names = T)
    data_sp <- get(load(files))
    Fmax <- unlist(stringr::str_split(unlist(stringr::str_split(folders[f_sp], "Fmax"))[2], "_"))[1]
    Fstart <- unlist(stringr::str_split(unlist(stringr::str_split(folders[f_sp], "Fstart"))[2], "_"))[1]
  
    
    INIT2b = data_sp
    INIT2b[2:3] = lapply(INIT2b[2:3], FUN = myf2)
    INIT2b[1:3] = lapply(INIT2b[1:3], FUN = reshape_list)
    INIT2b[[5]] = lapply(INIT2b[[5]], FUN = as.list)
    INIT2b[[5]] = reshape_list(INIT2b[[5]])
    INIT2b[[5]] = lapply(INIT2b[[5]], FUN = myf)
    INIT2b[[5]]$age = NULL
    INIT2b[[1]] = NULL
    INIT2b = reshape_list(INIT2b)

    
    longeH = reshape_data(INIT2b[[SP]]$data_Ne$longevity, "longevity", as.variable = TRUE)
    Ne_N1H = reshape_data(INIT2b[[SP]]$data_Ne$Ne_N_1, "Ne_N1", as.variable = TRUE)
    Ne_NalphaH = reshape_data(INIT2b[[SP]]$data_Ne$Ne_N_alpha, "Ne_Nalpha", as.variable = TRUE)
    Vk_pointH = reshape_data(INIT2b[[SP]]$data_Ne$Vk_point, "Vk_point", as.variable = TRUE)
  
    
    Nb_H_l = expand.grid(cohort=(1:nrow(INIT2b[[SP]]$data_Nb$Nb))+80,replicate=as.numeric(gsub(colnames(INIT2b[[SP]]$data_Nb$Nb), pattern="rep_", replacement="")))
    Nb_H_l$value = as.vector(INIT2b[[SP]]$data_Nb$Nb)
    Nb_H_l$variable = "Nb"
    Ne_t = myf(INIT2b[[SP]]$tempoNe)
    Ne_t_l = expand.grid(cohort=(1:nrow(Ne_t))+80,replicate=as.numeric(gsub(colnames(Ne_t), pattern="rep_", replacement="")))
    Ne_t_l$value = as.vector(Ne_t)
    Ne_t_l$variable = "Ne_t"
    htzH_l = expand.grid(cohort=(1:nrow(INIT2b[[SP]]$htz_list))+80,replicate=as.numeric(gsub(colnames(INIT2b[[SP]]$htz_list), pattern="rep_", replacement="")))
    htzH_l$value = as.vector(INIT2b[[SP]]$htz_list)
    htzH_l$variable = "htz"

    
    dat = rbind(longeH, Ne_N1H, Ne_NalphaH, htzH_l, Ne_t_l, Nb_H_l, Vk_pointH)
    dat$Fmax = Fmax
    dat$Fstart = Fstart
    # dat$Species = SP
    
    if(f_sp != ID_fold[1]) dat_tot = rbind(dat_tot, dat)
    if(f_sp == ID_fold[1]) dat_tot = dat
  }
  
  dat_tot <- rbind(dat_tot, mortalityH_l)
  
  dat_plot <- dat_tot
  
  Y_TRAIT = "htz"#c("Ne_N1")
  
  # dat_plot <- dat_tot[which(dat_tot$Species==SP),]
  dat_plot$y <- c(dat_plot$variable %in% Y_TRAIT)
  for (trait in names(table(dat_plot$variable[-which(dat_plot$y)]))){
    dat_plot_t <- dat_plot[which((dat_plot$cohort %in% c((min(dat_plot$cohort)+20):(max(dat_plot$cohort)-20)))),]
    dat_plot_t <- dat_plot_t[which((dat_plot_t$variable %in% c(trait)) | dat_plot_t$y),]
    dat_plot_t <- dat_plot_t[, -which(colnames(dat_plot_t) == "y")]
    dat_plot_w <- dat_plot_t %>%
      pivot_wider(
        names_from  = variable,
        values_from = value
      )
    ID_y <- which(colnames(dat_plot_w) == Y_TRAIT)
    ID_x = ifelse(ID_y==5, c(6),c(5))
    y_var = colnames(dat_plot_w)[ID_y]
    if (Y_TRAIT %in% c("htz")) {
      dat_plot_w[, Y_TRAIT] =  log(dat_plot_w[,Y_TRAIT])
      y_var <- paste0("log ",y_var)
    }
    if (all(is.na(dat_plot_w[, trait]))) next
    p[[SP]][[trait]] <- ggplot(dat_plot_w, aes_string(x=colnames(dat_plot_w)[ID_x], y=colnames(dat_plot_w)[ID_y]))+
      geom_point()+
      geom_smooth(method = loess, formula = y ~ x, se = FALSE)+
      facet_grid(Fstart~Fmax,labeller = labeller(.rows = label_both, .cols = label_both))+
      labs(title = paste0(SP, " ", y_var, " against ", colnames(dat_plot_w)[ID_x]))+
      ylab(label = y_var)+
      geom_abline(intercept = 0, slope = 1, col="red")
  }
  pathout = paste0("osmose-ns-v4.4.1_clim/output/PLOTS/Var_VS_Var/log",Y_TRAIT)
  if(!dir.exists(pathout)) dir.create(pathout)
  for (n in names(p)){
    pdf(file = paste0(pathout, "/plot_", n, ".pdf"), width = 10, height = 8)
    par(mfrow=c(1,1))
    print(p[[SP]])
    dev.off()
  }
}



# Plots -------------------------------------------------------------------
