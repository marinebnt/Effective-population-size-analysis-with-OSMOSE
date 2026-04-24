library(GGally)
library(PerformanceAnalytics)
library(lme4)
library("buildmer")
library(FactoMineR)
library(factoextra)


extract_INIT2_data_Ne_weighted_by_age <- function(pathoutputRDATA,
                                                  folders_in,
                                                  species = NULL,
                                                  abundance_var = "abundanceByAge",
                                                  Time1, Time2,
                                                  keep_species_column = TRUE, FORCE =F) {
  file_output <- paste0(pathoutputRDATA, "/data_Ne_weighted_by_age_", Time1, "_", Time2, ".RData")                                                      
  if (file.exists(file_output) & !FORCE) {
    data_out_all <- get(load(file = file_output))
  } else {
    data_out_all = NULL
    for (SP in species){      cat(SP, "\n")

      # initialisation
      p = list()
      ID_fold <- grep(SP, folders_in)
      k=0
      
      # Loop for folders
      for (f_sp in ID_fold){

        # initialisation
        k=k+1
        files  <- list.files(folders_in[f_sp], full.names = T)
        data_sp<- get(load(files))
        Fmax   <- unlist(stringr::str_split(unlist(stringr::str_split(folders_in[f_sp], "Fmax"))[2], "_"))[1]
        Fstart <- unlist(stringr::str_split(unlist(stringr::str_split(folders_in[f_sp], "Fstart"))[2], "_"))[1]
        
        # load file and match the dimensions of var = f(rep, Fish_pressure, time/cohort)
        INIT2  = get(load(files))
        # INIT2b = lapply(INIT2, FUN = myf2)
        # INIT2b = reshape_list(INIT2b)

        # by sp, by genetic drift dimension
        for (species_plotted in SP){
          for (genetic_drift_i in names(INIT2)[1]){
            data_ne = INIT2[[genetic_drift_i]][[species_plotted]]

            matrix_to_long <- function(mat, replicate_id, variable_name, species_name) {
              if (is.null(mat)) return(NULL)
              mat_df <- as.data.frame(mat)
              if (nrow(mat_df) == 0 || ncol(mat_df) == 0) return(NULL)

              age_vals <- suppressWarnings(as.numeric(rownames(mat_df)))
              if (all(is.na(age_vals))) age_vals <- seq_len(nrow(mat_df))
              mat_df$age <- age_vals

              long_df <- tidyr::pivot_longer(mat_df, cols = -age, names_to = "cohort", values_to = "value")
              long_df$cohort <- suppressWarnings(as.numeric(stringr::str_remove(long_df$cohort, "cohort_")))
              long_df$replicate <- as.numeric(replicate_id)
              long_df$variable <- variable_name
              long_df$Species <- species_name
              long_df
            }

            variable_to_long <- function(obj, variable_name, species_name) {
              out <- NULL

              if (is.list(obj) && !is.data.frame(obj)) {
                rep_names <- names(obj)
                for (k in seq_along(obj)) {
                  rep_id <- if (!is.null(rep_names) && nzchar(rep_names[k])) {
                    suppressWarnings(as.numeric(stringr::str_extract(rep_names[k], "\\d+")))
                  } else {
                    k - 1
                  }
                  if (is.na(rep_id)) rep_id <- k - 1
                  one <- matrix_to_long(obj[[k]], rep_id, variable_name, species_name)
                  if (!is.null(one)) out <- if (is.null(out)) one else rbind(out, one)
                }
                return(out)
              }

              if (is.matrix(obj) || is.data.frame(obj)) {
                return(matrix_to_long(obj, 0, variable_name, species_name))
              }

              NULL
            }

            raw_long <- NULL
            for (var in names(data_ne)) {
              if (!is.list(data_ne[[var]])) next
              var_long <- variable_to_long(data_ne[[var]], var, species_plotted)
              if (!is.null(var_long)) raw_long <- if (is.null(raw_long)) var_long else rbind(raw_long, var_long)
            }
            if (is.null(raw_long) || nrow(raw_long) == 0) return(NULL)

            raw_long <- raw_long[!is.na(raw_long$cohort) & raw_long$cohort %in% c(Time1:Time2), ]
            if (nrow(raw_long) == 0) return(NULL)

            weight_df <- build_weight_df(raw_long, abundance_var)

            data_with_w <- raw_long %>%
              dplyr::left_join(weight_df, by = c("Species", "replicate", "age", "cohort"))

            data_out = NULL
            data_out <- data_with_w %>%
              dplyr::group_by(Species, variable, replicate, age) %>%
              dplyr::summarise(
                value = {
                  if (dplyr::first(variable) == abundance_var) {
                    mean(value, na.rm = TRUE)
                  } else if (all(is.na(weight_value)) || sum(weight_value, na.rm = TRUE) == 0) {
                    mean(value, na.rm = TRUE)
                  } else {
                    mean(value, na.rm = TRUE)
                    # stats::weighted.mean(value, w = weight_value, na.rm = TRUE)
                  }
                },
                .groups = "drop"
              )
            data_out <- clean_values(data_out)

            data_out$Fmax <- Fmax
            data_out$Fstart <- Fstart

            if (!keep_species_column) {
              data_out <- data_out[, c("variable", "value", "replicate", "age", "Fmax", "Fstart")]
            }
            
            if (!is.null(data_out)) data_out_all <- if (is.null(data_out_all)) data_out else rbind(data_out_all, data_out)
            
          }
        }
      }
    }
    save(data_out_all, file = file_output)
  }
  return(data_out_all)
}


data_trait_and_genetic_into_df <- function(folders_in, Time1, Time2, pathoutput, FORCE = F, 
                                           species, RATIO = FALSE, keep_cohort_or_age_dimension = FALSE){
  data_full_merged = NULL
  RAT_name = ""
  if (RATIO == TRUE){
    RAT_name = "ratio"
  }
  if (!file.exists(paste0(pathoutput, "/genet_drift_data", RAT_name, "_",Time1,"_", Time2, ".RData")) & 
      !file.exists(paste0(pathoutput, "/trait_data", RAT_name, "_",Time1,"_", Time2, ".RData")) | 
      FORCE == T){
    genet_drift_data = list()
    trait_data = list()
    for (SP in species){      cat(SP, "\n")

      # initialisation
      dat_tot = NULL
      dat_plot = NULL
      p = list()
      ID_fold <- grep(SP, folders_in)
      k=0
      
      # Loop for folders
      for (f_sp in ID_fold){

        # initialisation
        k=k+1
        files  <- list.files(folders_in[f_sp], full.names = T)
        data_sp<- get(load(files))
        Fmax   <- unlist(stringr::str_split(unlist(stringr::str_split(folders_in[f_sp], "Fmax"))[2], "_"))[1]
        Fstart <- unlist(stringr::str_split(unlist(stringr::str_split(folders_in[f_sp], "Fstart"))[2], "_"))[1]
        
        # load file and match the dimensions of var = f(rep, Fish_pressure, time/cohort)
        INIT2  = get(load(files))
        INIT2b = lapply(INIT2, FUN = myf2)
        INIT2b = reshape_list(INIT2b)
        
        # by sp, by genetic drift dimension
        for (species_plotted in SP){
          for (genetic_drift_i in names(INIT2b[[species_plotted]])){
            
            if(genetic_drift_i != "htz_list"){
              # add variable : convert abundance by age into a NT variable (starting at age alpha, NT=Nalpha*sum(l_x))
              matrix_NT <- as.data.frame(sapply(1:length(INIT2b[[species_plotted]]$data_Ne$abundanceByAge),
                                                FUN = function(x) {
                                                  product <- INIT2b[[species_plotted]]$data_Ne$abundanceByAge[[x]] * INIT2b[[species_plotted]]$data_Ne$survivorship_curve[[x]]
                                                  colsums <- colSums(product)
                                                }))
              colnames(matrix_NT) <- paste0(1:ncol(matrix_NT))
              matrix_NT$cohort    <- stringr::str_remove(rownames(matrix_NT), pattern="cohort_")
              matrix_NT = pivot_longer(data = matrix_NT, cols = -c("cohort"), names_to = "replicate", values_to = "value")
              matrix_NT$variable = "NT"
              INIT2b[[species_plotted]][[genetic_drift_i]][["NT"]]<- matrix_NT
            }
            
            # initialisation
            names_var_genet<- names(INIT2b[[species_plotted]][[genetic_drift_i]])

            # differenciate which variables are genetic drift indicators and which are traits
            ID_Ne <- c(grep(pattern = "Ne", names_var_genet),  grep(pattern = "Nb", names_var_genet))
            ID_df <- as.numeric(which(lapply(INIT2b[[species_plotted]][[genetic_drift_i]], function(x) is.data.frame(x))==TRUE))
            
            # from HTZ time series to HTZ_slope
            if (grepl(pattern = "htz", genetic_drift_i)){

              # load data & clean it 
              data_t <- NULL
              data_t <- reshape_data(INIT2b[[species_plotted]][[genetic_drift_i]], varName = genetic_drift_i, as.variable = TRUE)
              if (any(grepl(pattern="cohort", colnames(data_t)))) { 
                  data_t$cohort <- stringr::str_replace(string = data_t$cohort, pattern = "cohort_", "")
                  data_t <- subset(data_t, cohort >= Time1 & cohort <= Time2) 
              }
              data_t <- na.omit(data_t)
              if (length(which(table(data_t$replicate)==1))>1) data_t <- data_t[-which(data_t$replicate %in% names(table(data_t$replicate))[which(table(data_t$replicate)==1)]),]
              
              # if no htz data
              if (nrow(data_t)==0) {cat("EMPTY data_t", folders[f_sp], "\n") 
                 next}

              # estimate the slope
              if (!keep_cohort_or_age_dimension){
                data_t <- data_t %>% # SLOPE HTZ
                  group_by(replicate) %>%
                  do({
                    model = lm(value ~ cohort, data = .)
                    data.frame(value = coef(model)[2])
                  }) %>%
                  ungroup()
              }
              data_t$Fstart  <- Fstart
              data_t$Fmax    <- Fmax
              data_t$Species <- species_plotted
              data_t$variable<- "htz_slope"
              
              if (is.null(trait_data[[genetic_drift_i]])){
                trait_data[[genetic_drift_i]] = data_t
              } else {
                trait_data[[genetic_drift_i]] = rbind(trait_data[[genetic_drift_i]], data_t)
              }

            } else {

              for (i in names_var_genet[ID_df]){
                if (grepl(pattern = "b_rate", i) | grepl(pattern = "eggNumber_sc", i) | grepl(pattern = "birthrate_scaled", i)) next

                # load data & clean it
                data_t_i <- INIT2b[[species_plotted]][[genetic_drift_i]][[i]]
                
                if (any(grepl(pattern="cohort", colnames(data_t_i)))) {
                  data_t_i$cohort <- stringr::str_replace(string = data_t_i$cohort, pattern = "cohort_", "")
                  data_t_i <- subset(data_t_i, cohort >= Time1 & cohort <= Time2) 
                }
                if (any(grepl(pattern="time", colnames(data_t_i))))   data_t_i <- subset(data_t_i, time >= Time1 & time <= Time2) 
                
                # calculate Ne/N if required
                if (RATIO == TRUE &  (i %in% names(INIT2b[[species_plotted]][[genetic_drift_i]])[ID_Ne])){
                  abund = NULL
                  if (grepl(i, pattern = "Ne_N_1") | grepl(i, pattern = "Ne_N_alpha")) abund <- INIT2b[[species_plotted]][[genetic_drift_i]][[grep(names_var_genet, pattern = "Numerator_N_alpha")]]
                  if (grepl(i, pattern = "Nb") | grepl(i, pattern = "Ne_temporal"))    abund <- INIT2b[[species_plotted]][["data_Nb"]][[grep(names(INIT2b[[species_plotted]][["data_Nb"]]), pattern = "NT")]]
                  if (!is.null(abund)){
                    abund$cohort <- stringr::str_remove(abund$cohort, pattern = "cohort_")
                    abund <- pivot_wider(data=abund, names_from  = variable, values_from = value)
                    tempo_data_t <- pivot_wider(data=data_t_i, names_from  = variable, values_from = value)
                    abund <- suppressMessages(na.omit(full_join(abund, tempo_data_t))) # join both dataset to assure replicates and cohorts are matching
                    ratio_var_name <- paste0(colnames(abund)[4], "/", colnames(abund)[3])
                    abund[,ratio_var_name] <- abund[, 4]/abund[,3]
                    abund <- abund[, which(colnames(abund) %in% c("cohort", "replicate", ratio_var_name))]
                    abund <- pivot_longer(data=abund, names_to = "variable", cols = ratio_var_name, values_to = "value")
                    data_t_i <- rbind(data_t_i, abund)
                  }
                }
                data_t <- data_t_i
                # mean over time or cohort
                if (!keep_cohort_or_age_dimension & data_t_i$variable[1] != "NT"){
                  data_t_i <- rbind(data_t_i, (INIT2b[[species_plotted]][[genetic_drift_i]][["NT"]]))
                  data_t_i <- pivot_wider(data_t_i, names_from = "variable", values_from = "value")
                  
                  data_t_i <- data_t_i %>%
                    group_by(replicate) %>%
                    # summarise(across(where(is.numeric), mean, na.rm = TRUE)) %>%
                    summarise(weighted.mean(w= NT, 
                                            x= .data[[data_t$variable[1]]], na.rm = TRUE)) %>%
                    ungroup()
                  colnames(data_t_i)<- c("replicate", "value")
                  data_t_i$variable <- data_t$variable[1] 
                  if (!(i %in% names(INIT2b[[species_plotted]][[genetic_drift_i]])[ID_Ne])) data_t_i$variable<- paste0("mean_", data_t_i$variable) 
                  data_t_i -> data_t
                }
                if (!keep_cohort_or_age_dimension & data_t_i$variable[1] == "NT"){
                  data_t <- data_t_i %>%
                    group_by(replicate, variable) %>%
                    summarise(across(where(is.numeric), mean, na.rm = TRUE)) %>%
                    ungroup()                 
                  if (!(i %in% names(INIT2b[[species_plotted]][[genetic_drift_i]])[ID_Ne])) data_t$variable<- paste0("mean_", data_t$variable) 
                }
                  
                data_t$Fstart  <- Fstart
                data_t$Fmax    <- Fmax
                data_t$Species <- species_plotted

                # store data
                if (i %in% names(INIT2b[[species_plotted]][[genetic_drift_i]])[ID_Ne]){
                  if (is.null(genet_drift_data[[genetic_drift_i]])){
                    genet_drift_data[[genetic_drift_i]] = data_t
                  } else {
                    genet_drift_data[[genetic_drift_i]] = rbind(genet_drift_data[[genetic_drift_i]], data_t)
                  }
                }else{      
                  if (is.null(trait_data[[genetic_drift_i]])){
                    trait_data[[genetic_drift_i]] = data_t
                  } else {
                    trait_data[[genetic_drift_i]] = rbind(trait_data[[genetic_drift_i]], data_t)
                  }
                } 
              } #each sub element of the list
            } # if/else

            # finalise data_t

          } # each element of the list
        } # each species of the list
      } # each file for this specific species
    }#species

    # save data
    save(genet_drift_data, file=paste0(pathoutput, "/genet_drift_data", RAT_name, "_",Time1,"_", Time2, ".RData"))
    save(trait_data, file=paste0(pathoutput, "/trait_data_",Time1,"_", Time2, ".RData"))
  } else {
    genet_drift_data <- get(load(file=paste0(pathoutput, "/genet_drift_data", RAT_name, "_" ,Time1,"_", Time2, ".RData")))
    trait_data <- get(load(file=paste0(pathoutput, "/trait_data_",Time1,"_", Time2, ".RData")))
  }
  return(list("genet_drift_data"=genet_drift_data, "trait_data"=trait_data))
}



dataset_attach_all_traits_genetic_drift <- function(pathoutputRDATA, s_pecies = species, 
                              list_outputs, FORCE = FALSE, trait_cv, keep_cohort_or_age_dimension = FALSE){

  expand_with_dimension <- function(df, dim_name, dim_values) {
    if (dim_name %in% colnames(df)) return(df)
    if (length(dim_values) == 0) {
      df[[dim_name]] <- NA_real_
      return(df)
    }
    df_expanded <- df[rep(seq_len(nrow(df)), each = length(dim_values)), , drop = FALSE]
    df_expanded[[dim_name]] <- rep(dim_values, times = nrow(df))
    df_expanded
  }

  # initialisation
  age_nm = ""
  cohort_or_age = "cohort"
  trait_by_age = list()
  if (any(grepl(names(trait_cv), pattern="_a"))) {
    age_nm = "_age"
    cohort_or_age = "age"
  }
  f = names(list_outputs)[1]
  if (keep_cohort_or_age_dimension) f <- paste0(f, "WITHAGE")
  # load or create dataset merging traits and genetic drift data
  if (!file.exists(paste0(pathoutputRDATA, "/PLOT_",f,"_", age_nm,"_only_trait_correlations_chart.RData")) | FORCE == T){
    data_full_merged = NULL
    for (f in names(list_outputs)){

      # load elements of the list 
      genet_drift_data_full       <- list_outputs[[f]][["genet_drift_data"]] # genetic drift indicators
      trait_data_full             <- list_outputs[[f]][["trait_data"]]       # trait data involved in the estimation of the genetic drift indicators

      for (SP in s_pecies){
        cat(SP)
        for (species_plotted in SP){
          cat(species_plotted, "\n")
          for (i in names(genet_drift_data_full)[1]){
            cat("   ", i, "\n")

          # selecting only the genetic drift data of the species of interest
            genet_drift_di = genet_drift_data_full[[i]]
            genet_drift_di <- genet_drift_di[which(genet_drift_di$Species == species_plotted), ]
            if(nrow(genet_drift_di) == 0) next

          # selecting only the trait data involved in the estimation of the genetic drift indicators
            trait_data_full_i <- trait_data_full[[i]]
            trait_data_full_i <- trait_data_full_i[which(trait_data_full_i$Species == species_plotted), ]

          # loading trait data NOT involved in the estimation of the genetic drift indicators
            # choosing trait data suited to compare genetic drift indicators with : if Ne_theoretical : by cohort, else by year
            folders_traits_y = trait_cv[[grep(names(trait_cv), pattern="y")]]
            folders_traits_c = trait_cv[[grep(names(trait_cv), pattern="c")]]
            ID_fold_y <- grep(folders_traits_y, pattern=species_plotted)
            ID_fold_c <- grep(folders_traits_c, pattern=species_plotted)
            if (any(grepl(i, pattern="Ne_N_1"),grepl(i, pattern="Ne_N_alpha"),grepl(i, pattern="data_Ne"))) {
              trait_file_i <- list.files(folders_traits_c[ID_fold_c], full.names = T)
              var_dim = "cohort_"
            } else {
              trait_file_i <- list.files(folders_traits_y[ID_fold_y], full.names = T)
              var_dim = "year_"
            }
            trait_df <- as.tibble(get(load(trait_file_i)) )
            
            # harmonize the replicates order
            if(all(range(trait_df$replicate) != range(as.numeric(trait_data_full_i$replicate)))){
              diff <- (min(trait_df$replicate) - min(as.numeric(trait_data_full_i$replicate)))
              trait_df$replicate <- as.numeric(trait_df$replicate) - diff
            }
            
            # cleaning trait data
            # uniformize column names
            colnames(trait_df)[which(colnames(trait_df) == "TRAIT")] <- "variable"
            trait_df$variable <- paste0(var_dim, trait_df$variable)                    # specify if traits are cohort- or  year-based
            if (any(trait_df$time==0)) trait_df <- trait_df[-which(trait_df$time==0),] # remove time 0 if any
            if (length(trait_df$rep)>0) {
              trait_df$replicate <- trait_df$rep
              trait_df$rep <- NULL
            }
            # adapt the dataset to the time or cohort or age dimension and store it differently if required
            if (length(trait_df$time)>0) {
              trait_df[,cohort_or_age] <- trait_df$time
              trait_df$time <- NULL
            }
            if (cohort_or_age != "age") trait_df <- trait_df[as.numeric(as.matrix(trait_df[,cohort_or_age])) %in% c(Time1:Time2),] 
            if (cohort_or_age == "age") trait_df <- trait_df[as.numeric(as.matrix(trait_df[,cohort_or_age])) > 0,]

            # add Fmax and Fstart column from Fish_press column
            trait_df$Fmax   <-  unlist(stringr::str_split(unlist(stringr::str_split(trait_df$Fish_press, pattern="Fmax="))[seq(by=2, 2, 2*nrow(trait_df))], pattern = " "))[seq(by=2, 1, 2*nrow(trait_df))]
            trait_df$Fstart <-  unlist(stringr::str_split(trait_df$Fish_press, pattern="Fstart="))[seq(by=2, 2, 2*nrow(trait_df))]
            trait_df        <- trait_df[, -which(colnames(trait_df) %in% c("Fish_press"))]
            # select only the trait data of the species of interest
            trait_df <- trait_df[which(trait_df$Species == species_plotted), ]
            if (length(trait_df$var_F_for_sp) > 0) trait_df <- trait_df[which(trait_df$var_F_for_sp == species_plotted), ]

            # if keeping the age dimension 
            if(keep_cohort_or_age_dimension) {
              if (is.null(trait_by_age[[i]])){
                trait_by_age[[i]] <- trait_df
              } else {
                trait_by_age[[i]] <- rbind(trait_by_age[[i]], trait_df)
              }
            }

            dim_values <- NULL
            if (keep_cohort_or_age_dimension && cohort_or_age %in% colnames(trait_df)) {
              dim_values <- sort(unique(trait_df[[cohort_or_age]]))
            }

            # the trait data NOT involved in the estimation of the genetic drift indicators has to be averaged over the period of interest
            # weighted by the abundance if the trait is a rate (except for birthrate_scaled and eggNumber_sc which are already scaled by abundance) to be comparable to the genetic drift indicators which are mean over the period of interest.
            weight_var <- character(0)
            weight_candidates <- c("abundanceBy")
            unique_vars <- unique(trait_df$variable)
            weight_var <- unique_vars[grep(pattern = weight_candidates, unique_vars)]
            weight_df <- trait_df[which(trait_df$variable %in% weight_var), ]
            weight_df <- pivot_wider(weight_df, names_from = variable, values_from = value)
            colnames(weight_df)[ncol(weight_df)] <- "weight_value"
            trait_df <- trait_df %>% 
                      full_join(weight_df)
            exceptions <- c("schoolNumber", "biomass", "abundance", "eggNumber")
            id_to_not_wm_all <- NULL
            id_to_not_wm <- NULL
            id_to_not_wm_all <- unlist(sapply(exceptions, function(x) {
              if (any(grepl(pattern = x, trait_df$variable)))  id_to_not_wm <- c(id_to_not_wm, grep(pattern = x, trait_df$variable))}))
            trait_df_to_wm  <- trait_df[-id_to_not_wm_all, ]
            trait_df_not_wm <- trait_df[id_to_not_wm_all, ]
            trait_df <- trait_df %>%
                    group_by(Species, variable, replicate, Fmax, Fstart) %>%
                    summarise(value = weighted.mean(value, w = weight_value, na.rm = TRUE))

          # Merging trait_df with already averaged trait data involved in the estimation of the genetic drift indicators
            # selecting the trait data involved in the estimation of the genetic drift indicators
            trait_df <- trait_df[, which(colnames(trait_df) %in% colnames(trait_data_full_i))] # be assured of their compatibility
            trait_data_full_i$replicate <- as.numeric(trait_data_full_i$replicate)
            trait_df <- rbind(trait_df, trait_data_full_i)
            # add htz slope if missing
            if (!any(grepl(unique(trait_df$variable), pattern="htz"))) trait_df <- rbind(trait_df, trait_data_full[["htz_list"]][which(trait_data_full[["htz_list"]]$Species == species_plotted), which(colnames(trait_data_full[["htz_list"]]) %in% colnames(trait_df))])
            # add NT if missing and needed
            trait_data_full[["data_Nb"]]$replicate <- as.numeric(trait_data_full[["data_Nb"]]$replicate)
            if (i == "tempoNe") { t_NT <- trait_data_full[["data_Nb"]][grep(trait_data_full$data_Nb$variable, pattern="NT"), ]
                trait_df <- rbind(trait_df,     t_NT[which(t_NT$Species == species_plotted),])
            } # because Nb and Netempo same reference : yearly
            
          # make sure the variables have the replicates in the same range of values
            trait_df$replicate <- as.numeric(trait_df$replicate)
            genet_drift_di$replicate <- as.numeric(genet_drift_di$replicate)
            if (any((range(trait_df$replicate)[2] - range(trait_df$replicate)[1])!=20, 
            (range(trait_df$replicate)[2] - range(trait_df$replicate)[1])!=20)){
              id_min <- which(c((range(trait_df$replicate)[2] - range(trait_df$replicate)[1]), 
            (range(trait_df$replicate)[2] - range(trait_df$replicate)[1]))!=19)
              if (1 %in% id_min)    {
                var_0 <- names(table(trait_df$variable[which(trait_df$replicate == range(trait_df$replicate))]))
                trait_df$replicate[which(trait_df$variable == var_0)] <- trait_df$replicate[which(trait_df$variable == var_0)] -1
              }
              if (2 %in% id_min){
                  var_0 <- names(table(genet_drift_di$variable[which(genet_drift_di$replicate == range(genet_drift_di$replicate))]))
                  trait_df$replicate[which(genet_drift_di$variable == var_0)] <- genet_drift_di$replicate[which(genet_drift_di$variable == var_0)] -1
              } 
            }

          # attach together trait and genetic drift data
            trait_df <- trait_df[, which(colnames(trait_df) %in% colnames(genet_drift_di))] # be assured of their compatibility
            dat_merged <- rbind(genet_drift_di, trait_df)
            dat_merged_w <- pivot_wider(dat_merged, names_from = variable, values_from = value)
    
            # CONTROL : are there replicates ? 
            qq <- dat_merged_w %>%
              dplyr::group_by(replicate, Fmax, Fstart, Species) %>%
              dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
              dplyr::filter(n > 1L) 
            # if (nrow(qq)!=0) {cat("duplicates")
            #   next}
            
            # skip missing columns 
            if (nrow(dat_merged_w) ==0) {
              cat("No merged data for ", species_plotted, " with Fmax=", Fmax, " and Fstart=", Fstart, "\n")
              next}
            
            # store final data
            if (is.null(data_full_merged[[i]])){
              data_full_merged[[i]] <- dat_merged_w
            } else {
              data_full_merged[[i]] <- rbind(data_full_merged[[i]], dat_merged_w)
            }
            rep <- data_full_merged[[i]] %>%
              dplyr::group_by(replicate, Fstart, Fmax, Species) %>%
              dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
              dplyr::filter(n > 1L)
          # if(nrow(rep)>1) stop(cat("duplicates : quit"))
            
          } 
        }# files
      }#species
      if (keep_cohort_or_age_dimension)  data_full_merged<- list("data_full_merged"=data_full_merged, "data_full_merged_age"=trait_by_age)
      if (keep_cohort_or_age_dimension) f <- paste0(f, "WITHAGE")
      save(data_full_merged, file=paste0(pathoutputRDATA, "/PLOT_",f,"_", age_nm,"_only_trait_correlations_chart.RData"))
    } 
  }else {
    f = names(list_outputs)[1]
    if (keep_cohort_or_age_dimension) f <- paste0(f, "WITHAGE")
    data_full_merged = get(load(paste0(pathoutputRDATA, "/PLOT_",f,"_", age_nm,"_only_trait_correlations_chart.RData")))
  }
  # return(list(data_full_merged, trait_by_age))
}
  
  
scale_dataset_list <- function(pathoutputRDATA, data_mean, FORCE=FALSE){
    data_mean_sc <- list()
    for (j in names(data_mean)){    
      data_mean_sc[[j]] <- data_mean[[j]]
      data_mean_sc[[j]][, -which(colnames(data_mean_sc[[j]]) %in%
                                   c("Species", "replicate", "Fmax", "Fstart"))] <- 
        apply(data_mean_sc[[j]][, -which(colnames(data_mean_sc[[j]]) %in%
                   c("Species", "replicate", "Fmax", "Fstart"))], 2, function(x) scale(x, ))
      data_mean_sc[[j]][,as.numeric(names(table(which(is.na(data_mean_sc[[j]]), arr.ind=TRUE)[,2])))] <-
        lapply(data_mean_sc[[j]][,as.numeric(names(table(which(is.na(data_mean_sc[[j]]), arr.ind=TRUE)[,2])))],
               function(x) if(all(is.na(x))) NULL else x)
      data_mean_sc[[j]] <- na.omit(data_mean_sc[[j]])
    }
  return(data_mean_sc)
} 
  
  
relative_dataset <- function(data_mean_sc, pathoutputRDATA, FORCE=FALSE, RATIO=F){

  # initialisation
  rt_name = ""
  if (RATIO==T) rt_name = "ratio_"

  if (!file.exists(paste0(pathoutputRDATA, "/LinearModel_",rt_name, "Ne_traits.RData")) | FORCE == TRUE){
    
    # initialisation
    p <- list()
    diff_data = list()
    
    for (j in names(data_mean_sc)){

      # initialisation
      my_data_all = NULL
      dat_merged_w_sc <- data_mean_sc[[j]]
      dat_merged_w_sc$F_p <- paste0("Fmax=", dat_merged_w_sc$Fmax, "_Fstart=", dat_merged_w_sc$Fstart)
      my_data = NULL
      
      for (spi in unique(dat_merged_w_sc$Species)){
        
        # subset species
        dat_merged_w_sc_sp <- dat_merged_w_sc[which(dat_merged_w_sc$Species == spi), ]
        
        # store aside the reference data (Fmax=1, Fstart=1)
        ref <- dat_merged_w_sc_sp[which(dat_merged_w_sc_sp$F_p == paste0("Fmax=", 1, "_Fstart=", 1)), ]
        ref <-  ref[,-which(colnames(ref) %in% c('Fmax', "Fstart", "F_p"))]
        # change colnames coherently
        colnames(ref)[-which(colnames(ref) %in% c("replicate", 'Fmax', "Fstart", "Species", "F_p"))] <- paste0("ref_", colnames(ref[, -which(colnames(ref) %in% c("replicate", 'Fmax', "Fstart", "Species", "F_p"))]))
        
        # loop over the fishing pressures, to each time substract the reference data
        for (F_i in unique(dat_merged_w_sc$F_p)){
          # subset fishing pressure
          dat_merged_w_sc_sp_F <- dat_merged_w_sc_sp[which(dat_merged_w_sc_sp$F_p == F_i), ]
          # merge with reference
          dat_merged_w_sc_sp_F <- left_join(dat_merged_w_sc_sp_F, ref, by=c("Species", "replicate"))
          # substract reference
          for (nm_ref in colnames(dat_merged_w_sc_sp_F)[-grep(colnames(dat_merged_w_sc_sp_F), pattern="ref_")]) {
            if (nm_ref %in% c("replicate", 'Fmax', "Fstart", "Species", "F_p")) next
            dat_merged_w_sc_sp_F[, paste0("diff_", nm_ref)] <- dat_merged_w_sc_sp_F[, paste0("ref_", nm_ref)] - dat_merged_w_sc_sp_F[, nm_ref]
          }
          # reshape data
          my_data <- dat_merged_w_sc_sp_F[, grep("diff_", x = colnames(dat_merged_w_sc_sp_F))]
          my_data <- cbind(dat_merged_w_sc_sp_F[, which(colnames(dat_merged_w_sc_sp_F)%in%c("replicate", 'Fmax', "Fstart", "Species", "F_p"))], my_data)
          my_data <- pivot_longer(my_data, cols= colnames(my_data)[-which(colnames(my_data)%in%c("replicate", 'Fmax', "Fstart", "Species", "F_p"))], 
                                  names_to = "variable", values_to = "value")
          # store data per F pressures 
          if (length(my_data_all)==0) {my_data_all <- my_data
          } else my_data_all <- rbind(my_data_all, my_data) 

        }
      } # species
      # store all the species and F pressures data of the j dimension 
      diff_data[[j]] <- my_data_all
    } # dimensions dataset

    save(diff_data, file=paste0(pathoutputRDATA, "/LinearModel_",rt_name, "Ne_traits.RData"))
  }else {
    diff_data = get(load(paste0(pathoutputRDATA, "/LinearModel_",rt_name, "Ne_traits.RData")))
  }
  return(diff_data)
}
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
