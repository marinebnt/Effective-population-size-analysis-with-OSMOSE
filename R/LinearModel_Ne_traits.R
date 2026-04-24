


main = "New_gpu_conf.R"
source("osmose-ns-v4.4.1_clim/R/functions/source.R")
base_dir <- "osmose-ns-v4.4.1_clim/output/PLOTS/All_clean_genetic_drift_indicators/INIT/"
folders  <- list.dirs(base_dir, recursive = F)
INIT2 <- get(load(list.files(folders[1], full.names = T)))

folders_traits_y <- "osmose-ns-v4.4.1_clim/output/PLOTS/Trait_variables_by_year/RData/"
folders_traits_y <- list.dirs(folders_traits_y, recursive = F)

folders_traits_c <- "osmose-ns-v4.4.1_clim/output/PLOTS/Trait_variables_by_cohort/RData/"
folders_traits_c  <- list.dirs(folders_traits_c, recursive = F)

folders_traits_y_a <- "osmose-ns-v4.4.1_clim/output/PLOTS/Trait_variables_by_year_by_age/RData/"
folder_traits_y_a <- list.dirs(folders_traits_y_a, recursive = T)

folders_traits_c_a <- "osmose-ns-v4.4.1_clim/output/PLOTS/Trait_variables_by_cohort_by_age/RData/"
folder_traits_c_a <- list.dirs(folders_traits_c_a, recursive = F)

# Auxiliar functions ------------------------------------------------------

# load
source("osmose-ns-v4.4.1_clim/R/functions/auxiliary_for_lists.R")
source("osmose-ns-v4.4.1_clim/R/functions/functions_output_SimuRData_into_df_for_plots.R")






# Loop for Fishing Pressure -----------------------------------------------


pathout <- paste0(configDir, "/output/PLOTS/Htz_vs_Ne_variables")

# from RData output of the loops producing the plots of Trait by time and Genetic drift indicators, to df
# produce a dataframe with the right format : Ne variables and htz variables together
Time1=150
Time2=160
species_df=species
list_data_traits <- data_trait_and_genetic_into_df(Time1=150, Time2=160, pathoutput=paste0(configDir, "/output/PLOTS/RDATA"), 
                                                   FORCE = F, species = species_df)

# produce a dataframe with the right format : Ne variables and htz variables together
species_df=species
list_data_traits_ratio <- data_trait_and_genetic_into_df(Time1=150, Time2=160, pathoutput=paste0(configDir, "/output/PLOTS/RDATA"), 
                                                                FORCE = F, species = species_df, RATIO=T)




# Bind together the trait data coming from similar dimensions : cohort or year data together (by genetic drift dimension)
pathoutputRDATA = paste0(configDir, "/output/PLOTS/RDATA")
data_full_merged <- dataset_attach_all_traits_genetic_drift(pathoutputRDATA, list_outputs=list("traits"=list_data_traits), 
                                            trait_cv = list("folders_traits_c"=folders_traits_c, "folders_traits_y"=folders_traits_y))
dataRATIO_full_merged <- dataset_attach_all_traits_genetic_drift(pathoutputRDATA, list_outputs=list("ratio"=list_data_traits_ratio), 
                                                            trait_cv = list("folders_traits_c"=folders_traits_c, "folders_traits_y"=folders_traits_y))


data_full_merged_age <- dataset_attach_all_traits_genetic_drift(pathoutputRDATA, list_outputs = list("traits"=list_data_traits), s_pecies = species,
                                                                trait_cv = list("folders_traits_c_a"=folder_traits_c_a, "folders_traits_y_a"=folder_traits_c_a),
                                                                keep_cohort_or_age_dimension = T, FORCE=F)
# 
# # Dataset ready for modeling --------------------------------------------
# # scaling as a standardized loop
data_mean_sc = scale_dataset_list(pathoutputRDATA, data_mean=data_full_merged[[1]], FORCE=FALSE)
dataRATIO_mean_sc = scale_dataset_list(pathoutputRDATA, data_mean=dataRATIO_full_merged[[1]], FORCE=F)
# 
# # new dataset of the difference between Fishing pressure 1,1 and Else--------------------
diff_data = relative_dataset(data_mean_sc, pathoutputRDATA, FORCE=FALSE)
diffRATIO_data = relative_dataset(dataRATIO_mean_sc, pathoutputRDATA, FORCE=F, RATIO=T)
# 
# 
# # ANALYSIS -----------------------
# 
# ### DIFF dataset -----------------
# 
# 
# a <- pivot_wider(diff_data$Ne_N_alpha, names_from = "variable", values_from = "value")
# model <- lmer(data=a, diff_htz_slope~+(diff_frac|F_p))
# summary(model)
# sjPlot:: tab_model(model)
# 
# 
# 
# 
# 
# 
# 
# 
# ### REGULAR dataset ---------
# 
# colnames(data_mean_sc$Ne_temporal)
# colnames(data_mean_sc$Nb)
# 
# 
# 
# model3a <- lmer(data=data_mean_sc$Ne_N_1,  Ne_N_1~mean_of_F_per_age*mean_of_Mpred_per_age+(mean_of_F_per_age*mean_of_Mpred_per_age|Species))
# model3b <- lmer(data=data_mean_sc$Ne_N_1,  Ne_N_1~mean_of_F_per_age*mean_of_Mpred_per_age+(Fmax+Fstart|Species)+(1|Fmax)+(1|Fstart))
# model3c <- lmer(data=data_mean_sc$Ne_N_1,  Ne_N_1~mean_of_F_per_age*mean_of_Mpred_per_age+(Species|Fmax)+(Species|Fstart))
# 
# # model3e11 <- lmer(data=data_mean_sc$Ne_N_1,  htz_slope~Ne_N_1+(Species|Fstart))
# # sjPlot:: tab_model(model3e11)
# data_mean_sc$Ne_N_1$F_p <- paste0("Fmax_", data_mean_sc$Ne_N_1$Fmax, "_Fstart_", data_mean_sc$Ne_N_1$Fstart)
# model3e12 <- lmer(data=data_mean_sc$Ne_N_1,  htz_slope~0+Ne_N_1+(Species|F_p))
# model3e15 <- lmer(data=data_mean_sc$Ne_N_1,  htz_slope~0+Ne_N_1+mean_of_F_per_age+(Species|F_p))
# model3e16 <- lmer(data=data_mean_sc$Ne_N_1,  htz_slope~0+Ne_N_1+mean_of_Mpred_per_age+(Species|F_p))
# sjPlot:: tab_model(model3e16)
# model3e14 <- lmer(data=data_mean_sc$Ne_N_1,  htz_slope~0+Ne_N_1+(Ne_N_1+Species|F_p))
# sjPlot:: tab_model(model3e14)
# # model3e13 <- lmer(data=data_mean_sc$Ne_N_1,  htz_slope~Ne_N_1+(Species|Fstart)+(1|Species))
# # sjPlot:: tab_model(model3e13)
# AIC(model3e13, model3e12, model3e15, model3e16)
# model3e2 <- lmer(data=data_mean_sc$Ne_N_alpha,  htz_slope~Ne_N_alpha+(Species|Fmax))
# model3e3 <- lmer(data=data_mean_sc$Nb,  htz_slope~Nb+(Species|Fmax))
# model3e4 <- lmer(data=data_mean_sc$Ne_N_1,  htz_slope~Ne_N_1+(Species|Fmax))
# 
# AIC(model3a, model3b)
# summary(model3b)
# coef(model3c)
# 
# sjPlot::plot_model(model3e12,
#                    # axis.labels=c("Urchin", "Depth", "Fish"),
#                    show.values=TRUE, show.p=TRUE,
#                    title="Effect of Herbivores on Coral Cover")
# sjPlot:: tab_model(model3a)
# 
# 
# #1
# 
# df_neN1 <- as.data.frame(effects::effect(term= "Ne_N_1", mod= model3e12, ))
# 
# plot <- ggplot() +
#   geom_point(data=(data_mean_sc$Ne_N_1), aes(Fmax, htz_slope)) +
#   geom_point(data=df_neN1, aes(x=Ne_N_1, y=fit), color="blue") +
#   geom_line(data=df_neN1, aes(x=Ne_N_1, y=fit), color="blue") +
#   geom_ribbon(data= df_neN1, aes(x=Ne_N_1, ymin=lower, ymax=upper), alpha= 0.3, fill="blue") +
#   labs(x="Fmax", y="Ne_N_1")
# plot
# 
# 
# 
# 
# atlcod   <- data_mean_sc$Ne_N_1[which(data_mean_sc$Ne_N_1$Species == "AtlanticCod"),]
# atlcodm1 <- lmer(data=data_mean_sc$Ne_N_1,  htz_slope~0+Ne_N_1+mean_of_Mpred_per_age+(1 + mean_of_Mpred_per_age | Fmax))
# sjPlot:: tab_model(atlcodm1)
# 
# 
# atlcod   <- data_mean_sc$Ne_N_alpha[which(data_mean_sc$Ne_N_1$Species == "AtlanticCod"),]
# atlcodm1 <- lmer(data=data_mean_sc$Ne_N_alpha,  htz_slope~0+Ne_N_alpha+mean_of_Mpred_per_age+(1 + mean_of_Mpred_per_age | Fmax))
# sjPlot:: tab_model(atlcodm1)
# 
# atlcod   <- data_mean_sc$Ne_N_alpha[which(data_mean_sc$Ne_N_1$Species == "AtlanticCod"),]
# atlcodm1 <- lm(data=data_mean_sc$Ne_N_alpha,  htz_slope~0+Ne_N_alpha+Species+Fmax)
# sjPlot:: tab_model(atlcodm1)
# 
# 
# atlcod2   <- atlcod[which(atlcod$Fmax == "0.5"),]
# atlcodm1 <- lmer(data=atlcod,  htz_slope~Ne_N_alpha+mean_of_Mpred_per_age+(1+mean_of_Mpred_per_age|Fstart))
# sjPlot:: tab_model(atlcodm1)
# 
# 
# 
# sole   <- data_mean_sc$Ne_N_alpha[which(data_mean_sc$Ne_N_alpha$Species == "CommonSole"),]
# sole$Fp <- paste0(sole$Fmax, "_", sole$Fstart)
# data_mean_sc$Ne_N_alpha$Fp <- paste0(data_mean_sc$Ne_N_alpha$Fmax, "_", data_mean_sc$Ne_N_alpha$Fstart)
# solem1 <- lmer(data=data_mean_sc$Ne_N_alpha,  htz_slope~Ne_N_alpha+mean_of_Mpred_per_age+(mean_of_Mpred_per_age | Fp))
# sjPlot:: tab_model(solem1)
# sole2   <- sole[which(sole$Fmax == "0.5"),]
# solem1 <- lmer(data=sole,  htz_slope~Ne_N_alpha+(1|Fstart))
# sjPlot:: tab_model(atlcodm1)
# 
# 
# 
# data_mean_sc$Ne_N_alpha
# predm1 <- lmer(data=data_mean_sc$Ne_N_alpha,  mean_of_Mpred_per_age~1+Ne_N_alpha+htz_slope+(1|Species))
# 
# sjPlot:: tab_model(predm1)
# coef(predm1)
# 
# 
# htzm1 <- lmer(data=data_mean_sc$Ne_N_1,  htz_slope~1+Ne_N_1+mean_of_Mpred_per_age+generation_length+
#                 (1+Ne_N_1+mean_of_Mpred_per_age+generation_length|Species))
# coef(htzm1)
# sjPlot:: tab_model(htzm1)
# 
# 
# 

for (i in names(diffRATIO_data) ){
  diffRATIO_data[[i]] <- pivot_wider(diffRATIO_data[[i]], names_from = "variable", values_from = "value")
  colnames(diffRATIO_data[[i]]) <- stringr::str_replace(colnames(diffRATIO_data[[i]]), pattern = "diff_", replacement = "")
}


# as expected : Ne ~ negatively with htz slope, except for Ne tempo, it has to be Ne/NT
# significative for Netempo and Ne_N1 and Ne_Nalpha
  lm(data=(diffRATIO_data$Nb),  htz_slope~1+mean_Nb)->nb
  summary(nb)
  sjPlot:: tab_model(nb)
  lm(data=(diffRATIO_data$Nb),  htz_slope~1+frac)->nb
  summary(nb)
  sjPlot:: tab_model(nb)
lm(data=(diffRATIO_data$Ne_N_1),  htz_slope~1+mean_Ne_N_1)-> nen1
sjPlot:: tab_model(nen1)
diffRATIO_data$Ne_N_0$Fp <- paste0(diffRATIO_data$Ne_N_1$Fmax, "_", diffRATIO_data$Ne_N_1$Fstart)
# diffRATIO_data$Ne_N_1$frac <- diffRATIO_data$Ne_N_1$mean_Ne_N_1/diffRATIO_data$Ne_N_1$mean_Numerator_N_1
nen1_2 <- lm(data=(diffRATIO_data$Ne_N_1),  htz_slope~1+frac)
coef(nen1_2)
sjPlot:: tab_model(nen1_2)
  diffRATIO_data$Ne_temporal$Fp <- paste0(diffRATIO_data$Ne_temporal$Fmax, "_", diffRATIO_data$Ne_temporal$Fstart)
  # diffRATIO_data$Ne_temporal$frac <- diffRATIO_data$Ne_temporal$mean_Ne_temporal/diffRATIO_data$Ne_temporal$mean_NT
  netempo <- lm(data=(diffRATIO_data$Ne_temporal),  htz_slope~1+frac)
  netempo2 <- lm(data=(diffRATIO_data$Ne_temporal),  htz_slope~1+mean_Ne_temporal)
  coef(netempo)
  sjPlot:: tab_model(netempo)
  sjPlot:: tab_model(netempo2)
lm(data=(diffRATIO_data$Ne_N_alpha),  htz_slope~1+mean_Ne_N_alpha)->nenalpha
sjPlot:: tab_model(nenalpha)
lm(data=(diffRATIO_data$Ne_N_alpha),  htz_slope~1+frac)->nenalpha
sjPlot:: tab_model(nenalpha)


# compare effect htz&Ne by Fp 
  diffRATIO_data$Nb$F_p <- as.factor(paste0("Fmax_", diffRATIO_data$Nb$Fmax, "_Fstart_", diffRATIO_data$Nb$Fstart))
  lmer(data=(diffRATIO_data$Nb),  htz_slope~1+mean_Nb+(1+mean_Nb|Species))->nb
  lmer(data=(diffRATIO_data$Nb),  htz_slope~1+mean_Nb+(1+mean_Nb|F_p))->nb
  lmer(data=(diffRATIO_data$Nb),  htz_slope~1+mean_Nb+(1|Species))->nb
  lmer(data=(diffRATIO_data$Nb),  htz_slope~1+mean_Nb+(1|F_p))->nb
  sjPlot:: tab_model(nb)
  coef(nb)
diffRATIO_data$Ne_N_1$F_p <- as.factor(paste0("Fmax_", diffRATIO_data$Ne_N_1$Fmax, "_Fstart_", diffRATIO_data$Ne_N_1$Fstart))
lmer(data=(diffRATIO_data$Ne_N_1),  htz_slope~1+mean_Ne_N_1+(1|F_p))-> nen1
lmer(data=(diffRATIO_data$Ne_N_1),  htz_slope~1+mean_Ne_N_1+(1|Species))-> nen1
lmer(data=(diffRATIO_data$Ne_N_1),  htz_slope~1+mean_Ne_N_1+(1+mean_Ne_N_1|F_p))-> nen1
lmer(data=(diffRATIO_data$Ne_N_1),  htz_slope~1+mean_Ne_N_1+(1+mean_Ne_N_1|Species))-> nen1
  diffRATIO_data$Ne_temporal$Fp <- paste0(diffRATIO_data$Ne_temporal$Fmax, "_", diffRATIO_data$Ne_temporal$Fstart)
  # diffRATIO_data$Ne_temporal$frac <- diffRATIO_data$Ne_temporal$mean_Ne_temporal/diffRATIO_data$Ne_temporal$mean_NT
  netempoa <- lmer(data=(diffRATIO_data$Ne_temporal),  htz_slope~1+mean_Ne_temporal+(1|F_p))
  netempob <- lmer(data=(diffRATIO_data$Ne_temporal),  htz_slope~1+mean_Ne_temporal+(1|Species))
  netempoc <- lmer(data=(diffRATIO_data$Ne_temporal),  htz_slope~1+mean_Ne_temporal+(1+mean_Ne_temporal|F_p))
  netempod <- lmer(data=(diffRATIO_data$Ne_temporal),  htz_slope~1+mean_Ne_temporal+(1+mean_Ne_temporal|Species))
  coef(netempo)
  sjPlot:: tab_model(netempoa)
  sjPlot:: tab_model(netempob)
  sjPlot:: tab_model(netempoc)
  sjPlot:: tab_model(netempod)
  AIC(netempoa, netempob, netempoc, netempod)
lmer(data=(diffRATIO_data$Ne_N_alpha),  htz_slope~1+mean_Ne_N_alpha+(1|F_p))->nenalpha
lmer(data=(diffRATIO_data$Ne_N_alpha),  htz_slope~1+mean_Ne_N_alpha+(1|Species))->nenalpha
lmer(data=(diffRATIO_data$Ne_N_alpha),  htz_slope~1+mean_Ne_N_alpha+(1+mean_Ne_N_alpha|F_p))->nenalpha
lmer(data=(diffRATIO_data$Ne_N_alpha),  htz_slope~1+mean_Ne_N_alpha+(1+mean_Ne_N_alpha|Species))->nenalpha
sjPlot:: tab_model(nenalpha)


# predation mortality by COHORT explaines heterozygosity slope better than by year
# fishing mortality by YEAR explaines heterozygosity slope better than by cohort

data_htz_year_cohort <- full_join(diffRATIO_data$Ne_N_alpha, diffRATIO_data$Ne_temporal, by=c("replicate", "Fmax", "Fstart", "Species", "htz_slope"))
data_htz_year_cohort$F_p <- paste0("Fmax_", data_htz_year_cohort$Fmax, "_Fstart_", data_htz_year_cohort$Fstart)
data_htz_year_cohort <- data_htz_year_cohort[-which(is.na(data_htz_year_cohort$year_MpredbyY)), ]

  Pred_cohort <- lmer(data=data_htz_year_cohort,  htz_slope~1+cohort_mean_of_Mpred_per_age+(1+cohort_mean_of_Mpred_per_age|F_p))
  Pred_cohort2 <- lmer(data=data_htz_year_cohort, htz_slope~1+cohort_mean_of_Mpred_per_age+(1|F_p))
  sjPlot:: tab_model(Pred_cohort)
  sjPlot:: tab_model(Pred_cohort2)
  coef(Pred_cohort)
  coef(Pred_cohort2)

  Pred_year <- lmer(data=data_htz_year_cohort,  htz_slope~1+year_MpredbyY+(1+year_MpredbyY|F_p))
  Pred_year2 <- lmer(data=data_htz_year_cohort, htz_slope~1+year_MpredbyY+(1|F_p))
  sjPlot:: tab_model(Pred_year)
  sjPlot:: tab_model(Pred_year2)
  
  Pred_cohort3 <- lmer(data=data_htz_year_cohort,  htz_slope~1+cohort_mean_of_Mpred_per_age+(1+cohort_mean_of_Mpred_per_age|Species)) #THIS ONE
  Pred_cohort4 <- lmer(data=data_htz_year_cohort, htz_slope~1+cohort_mean_of_Mpred_per_age+(1|Species))
  sjPlot:: tab_model(Pred_cohort3)
  sjPlot:: tab_model(Pred_cohort4)
  coef(Pred_cohort3)
  coef(Pred_cohort4)
  
  Pred_year3 <- lmer(data=data_htz_year_cohort,  htz_slope~1+year_MpredbyY+(1+year_MpredbyY|Species))
  Pred_year4 <- lmer(data=data_htz_year_cohort, htz_slope~1+year_MpredbyY+(1|Species))
  sjPlot:: tab_model(Pred_year3)
  sjPlot:: tab_model(Pred_year4)

AIC(Pred_cohort, Pred_cohort2, Pred_year, Pred_year2, Pred_cohort3, Pred_cohort4, Pred_year3, Pred_year4)


  Pred_cohort <- lmer(data=data_htz_year_cohort,  htz_slope~1+cohort_mean_of_F_per_age+(1+cohort_mean_of_F_per_age|F_p)) #THIS ONE
  Pred_cohort2 <- lmer(data=data_htz_year_cohort, htz_slope~1+cohort_mean_of_F_per_age+(1|F_p))
  sjPlot:: tab_model(Pred_cohort)
  sjPlot:: tab_model(Pred_cohort2)
  coef(Pred_cohort)
  coef(Pred_cohort2)

  Pred_year <- lmer(data=data_htz_year_cohort,  htz_slope~1+year_FbyY+(1+year_FbyY|F_p)) 
  Pred_year2 <- lmer(data=data_htz_year_cohort, htz_slope~1+year_FbyY+(1|F_p))# this one
  sjPlot:: tab_model(Pred_year)
  sjPlot:: tab_model(Pred_year2)
  
  Pred_cohort3 <- lmer(data=data_htz_year_cohort,  htz_slope~1+cohort_mean_of_F_per_age+(1+cohort_mean_of_F_per_age|Species)) #THIS ONE
  Pred_cohort4 <- lmer(data=data_htz_year_cohort, htz_slope~1+cohort_mean_of_F_per_age+(1|Species))
  sjPlot:: tab_model(Pred_cohort3)
  sjPlot:: tab_model(Pred_cohort4)
  coef(Pred_cohort3)
  coef(Pred_cohort4)
  
  Pred_year3 <- lmer(data=data_htz_year_cohort,  htz_slope~1+year_FbyY+(1+year_FbyY|Species))
  Pred_year4 <- lmer(data=data_htz_year_cohort, htz_slope~1+year_FbyY+(1|Species))
  sjPlot:: tab_model(Pred_year3)
  sjPlot:: tab_model(Pred_year4)

  AIC(Pred_cohort, Pred_cohort2, Pred_year, Pred_year2, Pred_cohort3, Pred_cohort4, Pred_year3, Pred_year4)


# what about total mortality x F_p to explain the intercept of the Ne ~ Mortality
# explication :     x$natural = x$pred + x$starv + x$other + x$out
#                   x$for = metabolisme oxydatif, uniquement avec Bioen
diffRATIO_data$Ne_N_alpha$F_p <- paste0("Fmax_", diffRATIO_data$Ne_N_alpha$Fmax, "_Fstart_", diffRATIO_data$Ne_N_alpha$Fstart)

Mort_cohort <- lmer(data=diffRATIO_data$Ne_N_alpha, mean_Ne_N_alpha~1+cohort_sum_of_mortalityByAge_per_age+(1|F_p)) 
Mort_cohort0 <- lm(data=diffRATIO_data$Ne_N_alpha, #THIS ONE
                    cohort_sum_of_mortalityByAge_per_age~1+cohort_mean_of_Mstar_per_age+cohort_mean_of_Mfor_per_age+cohort_mean_of_F_per_age+cohort_mean_of_Mpred_per_age) 
Mort_cohort1 <- lmer(data=diffRATIO_data$Ne_N_alpha,
                    cohort_sum_of_mortalityByAge_per_age~1+cohort_mean_of_Mstar_per_age+cohort_mean_of_Mfor_per_age+cohort_mean_of_F_per_age+cohort_mean_of_Mpred_per_age+(1|F_p))

AIC(Mort_cohort0, Mort_cohort1, Mort_cohort)

sjPlot:: tab_model(Mort_cohort)
sjPlot:: tab_model(Mort_cohort0)
sjPlot:: tab_model(Mort_cohort1)
coef(Pred_cohort)
coef(Pred_cohort2)

diffRATIO_data$Nb$F_p <- paste0("Fmax_", diffRATIO_data$Nb$Fmax, "_Fstart_", diffRATIO_data$Nb$Fstart)

Mort_year <- lmer(data=diffRATIO_data$Nb, mean_Nb~1+year_mortalityByY+(1|F_p))
Mort_year <- lmer(data=diffRATIO_data$Nb, htz_slope~1+year_mortalityByY+(1|F_p))
Mort_year0 <- lm(data=diffRATIO_data$Nb,
                   year_mortalityByY~1+year_MstarbyY+year_MforbyY+year_FbyY+year_MnatbyY+year_MpredbyY) #THIS ONE
Mort_year1 <- lmer(data=diffRATIO_data$Nb,
                 year_mortalityByY~1+year_MstarbyY+year_MforbyY+year_FbyY+year_MnatbyY+year_MpredbyY+(1|F_p)) 
Mort_year2 <- lmer(data=diffRATIO_data$Nb,
                   year_mortalityByY~1+year_MstarbyY+year_MforbyY+year_FbyY+year_MnatbyY+year_MpredbyY+(1|Species)) 
AIC(Mort_year0, Mort_year1, Mort_year,Mort_year2)

sjPlot:: tab_model(Mort_year)
sjPlot:: tab_model(Mort_year0)
sjPlot:: tab_model(Mort_year1)
coef(Mort_cohort)
coef(Mort_cohort2)




# what about variance in reproductive success x F_p to explain htz or Ne
# + compare the htz ~Ne and htz~raw formula Ne
diffRATIO_data$Ne_N_alpha$F_p <- as.factor(paste0("Fmax_", diffRATIO_data$Ne_N_alpha$Fmax, "_Fstart_", diffRATIO_data$Ne_N_alpha$Fstart))
diffRATIO_data$Ne_N_alpha$test_Ne <- (4*diffRATIO_data$Ne_N_alpha$mean_Numerator_N_1*diffRATIO_data$Ne_N_alpha$mean_generation_length)/(diffRATIO_data$Ne_N_alpha$mean_Vk_point+2)

test_Nea <- lm(data=diffRATIO_data$Ne_N_alpha, htz_slope~1+test_Ne) 
test_Neb <- lmer(data=diffRATIO_data$Ne_N_alpha, htz_slope~1+test_Ne+(1|Species)) 
test_Ne2a <- lm(data=diffRATIO_data$Ne_N_alpha, htz_slope~1+mean_Ne_N_alpha) 
test_Ne2b <- lmer(data=diffRATIO_data$Ne_N_alpha, htz_slope~1+mean_Ne_N_alpha+(1|Species)) 
sjPlot:: tab_model(test_Nea)
sjPlot:: tab_model(test_Neb)
sjPlot:: tab_model(test_Ne2a)
sjPlot:: tab_model(test_Ne2b)
AIC(test_Nea, test_Neb, test_Ne2a, test_Ne2b)

Vkptcohort <- lmer(data=diffRATIO_data$Ne_N_alpha, mean_Ne_N_alpha~1+mean_Vk_point+(1|F_p)) 
Vkptcohort <- lmer(data=diffRATIO_data$Ne_N_alpha, htz_slope~1+mean_Vk_point+(1+mean_Vk_point|Species)) 
VkptcohortMort <- lmer(data=diffRATIO_data$Ne_N_alpha, htz_slope~1+mean_Vk_point+
                         cohort_sum_of_mortalityByAge_per_age+(1+mean_Vk_point|Species)) 

  # la mortalité par pêche et par prédation explique la différence de htz entre les pressions de pêche
VkptcohortMort0a <- lmer(data=diffRATIO_data$Ne_N_alpha, htz_slope~1+mean_k_point+
                          mean_Vk_point+cohort_mean_of_Mfor_per_age+
                          cohort_mean_of_Mstar_per_age+cohort_mean_of_F_per_age+
                          cohort_mean_of_Mpred_per_age+(1|F_p)) 
  # la mortalité non par pêche explique la différence de htz entre les espèces
VkptcohortMort0b <- lmer(data=diffRATIO_data$Ne_N_alpha, htz_slope~1+mean_k_point+
                          mean_Vk_point+cohort_mean_of_Mfor_per_age+
                          cohort_mean_of_Mstar_per_age+cohort_mean_of_F_per_age+
                          cohort_mean_of_Mpred_per_age+(1|Species)) 
  # variabilité induite par la pression de pêche + species = le succès reproducteur est mis sur le meme pied d'égalité que la mortalité naturelle
VkptcohortMort0c <- lmer(data=diffRATIO_data$Ne_N_alpha, htz_slope~1+mean_k_point+
                           mean_Vk_point+cohort_mean_of_Mfor_per_age+
                           cohort_mean_of_Mstar_per_age+cohort_mean_of_F_per_age+
                           cohort_mean_of_Mpred_per_age+(1|Species)+(1|F_p)) # best
sjPlot:: tab_model(VkptcohortMort0a)
sjPlot:: tab_model(VkptcohortMort0b)
sjPlot:: tab_model(VkptcohortMort0c)
AIC(VkptcohortMort0a, VkptcohortMort0b,VkptcohortMort0c)
cor(diffRATIO_data$Ne_N_alpha$mean_k_point, diffRATIO_data$Ne_N_alpha$cohort_mean_of_Mpred_per_age)
diffRATIO_data$Ne_N_alpha$Fstart <- as.numeric(diffRATIO_data$Ne_N_alpha$Fstart)
diffRATIO_data$Ne_N_alpha$Fmax <- as.numeric(diffRATIO_data$Ne_N_alpha$Fmax)
pca1 <- PCA(diffRATIO_data$Ne_N_alpha[, which(colnames(diffRATIO_data$Ne_N_alpha) %in% 
                c("mean_Vk_point", "cohort_sum_of_mortalityByAge_per_age", "mean_k_point", 
                  "cohort_mean_of_F_per_age", "cohort_mean_of_Mpred_per_age", "mean_generation_length", "Fmax"))],
          scale.unit = TRUE,
          graph = FALSE)
fviz_pca_biplot(pca1, 
                col.ind = (diffRATIO_data$Ne_N_alpha$htz_slope+0.01+abs(min(diffRATIO_data$Ne_N_alpha$htz_slope))),
                alpha.ind = 0.5
                )

VkptcohortMort1 <- lmer(data=diffRATIO_data$Ne_N_alpha, htz_slope~1+mean_Vk_point+cohort_mean_of_Mfor_per_age+cohort_mean_of_Mstar_per_age+cohort_mean_of_F_per_age+cohort_mean_of_Mpred_per_age+(1+mean_Vk_point|Species)) 
Vkptcohort0 <- lm(data=diffRATIO_data$Ne_N_alpha, 
                   htz_slope~1+mean_Vk_point+mean_k_point+cohort_mean_of_Mfor_per_age+cohort_mean_of_Mstar_per_age+cohort_mean_of_F_per_age+cohort_mean_of_Mnat_per_age+cohort_mean_of_Mpred_per_age)
Vkptcohort2 <- lmer(data=diffRATIO_data$Ne_N_alpha, 
                  htz_slope~1+mean_Vk_point+mean_k_point+cohort_mean_of_Mfor_per_age+cohort_mean_of_Mstar_per_age+cohort_mean_of_F_per_age+cohort_mean_of_Mnat_per_age+cohort_mean_of_Mpred_per_age+
                    (1|F_p)) # not the same mortalities count when F_p or Species 
Vkptcohort3 <- lmer(data=diffRATIO_data$Ne_N_alpha, 
                    htz_slope~1+mean_Vk_point+mean_k_point+cohort_mean_of_Mfor_per_age+cohort_mean_of_Mstar_per_age+cohort_mean_of_F_per_age+cohort_mean_of_Mnat_per_age+cohort_mean_of_Mpred_per_age+
                      mean_generation_length+(1|Species)) 
Vkptcohort4 <- lmer(data=diffRATIO_data$Ne_N_alpha, 
                    htz_slope~1+mean_Vk_point+mean_k_point+cohort_mean_of_Mfor_per_age+cohort_mean_of_Mstar_per_age+cohort_mean_of_Mnat_per_age+cohort_mean_of_F_per_age+cohort_mean_of_Mnat_per_age+cohort_mean_of_Mpred_per_age+
                      mean_generation_length+(1|Species)+(0+mean_Vk_point+cohort_sum_of_mortalityByAge_per_age|F_p)) #*****#THIS ONE => Not surprising : intercept depends on the species, and slope on the fishing pressure
AIC(Vkptcohort3, Vkptcohort2, Vkptcohort4, Vkptcohort)

sjPlot:: tab_model(Vkptcohort)
sjPlot:: tab_model(VkptcohortMort)
sjPlot:: tab_model(VkptcohortMort0)
sjPlot:: tab_model(VkptcohortMort1)
sjPlot:: tab_model(Vkptcohort0)
sjPlot:: tab_model(Vkptcohort2)
sjPlot:: tab_model(Vkptcohort3)
sjPlot:: tab_model(Vkptcohort4)
coef(Pred_cohort)
coef(Pred_cohort2)



diffRATIO_data$Nb$F_p <- as.factor(paste0("Fmax_", diffRATIO_data$Nb$Fmax, "_Fstart_", diffRATIO_data$Nb$Fstart))
diffRATIO_data$Nb$test_Nb <- (diffRATIO_data$Nb$mean_k_point*diffRATIO_data$Nb$mean_NT-2)/(diffRATIO_data$Nb$mean_k_point-
                                                    1+diffRATIO_data$Nb$mean_V_k_point/diffRATIO_data$Nb$mean_k_point)

test_Nba <- lm(data=diffRATIO_data$Nb, htz_slope~1+test_Nb) 
test_Nba <- lm(data=diffRATIO_data$Nb, htz_slope~1+mean_Nb) 
test_Nbb <- lmer(data=diffRATIO_data$Nb, htz_slope~1+test_Nb+(1|F_p)) 
test_Nb2b <- lmer(data=diffRATIO_data$Nb, htz_slope~1+mean_Nb_N_alpha+(1+F_p|Species)) 
sjPlot:: tab_model(test_Nba)
sjPlot:: tab_model(test_Nbb)
sjPlot:: tab_model(Vkptyear1)

diffRATIO_data$Nb$F_p <- paste0("Fmax_", diffRATIO_data$Nb$Fmax, "_Fstart_", diffRATIO_data$Nb$Fstart)
diffRATIO_data$Nb <- left_join(diffRATIO_data$Nb, data_full_merged_age[[1]]$Ne_temporal[, c(1:4, grep(pattern="egg", colnames(data_full_merged_age[[1]]$Ne_temporal)))])
diffRATIO_data$Nb$year_mean_of_eggNumber_per_age <- scale(diffRATIO_data$Nb$year_mean_of_eggNumber_per_age)
Vkptyear <- lmer(data=diffRATIO_data$Nb, mean_Nb~1+year_mortalityByY+year_mean_of_eggNumber_per_age+(1|F_p))
Vkptyear <- lmer(data=diffRATIO_data$Nb, htz_slope~1+year_mean_of_eggNumber_per_age+(1|F_p))
Vkptyear0 <- lm(data=diffRATIO_data$Nb,
                 year_mortalityByY~1+year_MstarbyY+year_MforbyY+year_FbyY+year_MnatbyY+year_MpredbyY) #THIS ONE
Vkptyear1 <- lmer(data=diffRATIO_data$Nb,
                   year_mortalityByY~1+year_MstarbyY+year_MforbyY+year_FbyY+year_MnatbyY+year_MpredbyY+(1|F_p)) 

AIC(Vkptyear0, Vkptyear1, Vkptyear)

sjPlot:: tab_model(Vkptyear)
sjPlot:: tab_model(Vkptyear0)
sjPlot:: tab_model(Vkptyear1)
coef(Vkptcohort)
coef(Vkptcohort2)

AIC(Pred_cohort, Pred_cohort2, Pred_year, Pred_year2)



diffRATIO_data$Ne_temporal$F_p <- as.factor(paste0("Fmax_", diffRATIO_data$Ne_temporal$Fmax, "_Fstart_", diffRATIO_data$Ne_temporal$Fstart))




a <- rbind(list_data_traits$genet_drift_data$Ne_N_alpha, list_data_traits$trait_data$data_Ne, 
           list_data_traits$trait_data$htz_list[, -which(colnames(list_data_traits$trait_data$htz_list) == "SP_ref")])
a <- as.data.frame(pivot_wider(a, names_from = "variable", values_from = "value"))
a[, -c(1:4)] <- scale(a[, -c(1:4)])
lm_t <- lm(data=a, htz_slope~1+mean_Ne_N_alpha)
summary(lm_t)
# 
# 
# 
# # relationship between Ne and predation mortality is not linear and depends on the functional traits (?)
#   predmNe_N_alpha2 <- lmer(data=data_mean_sc$Ne_N_alpha,  cohort_mean_of_Mpred_per_age~1+mean_Ne_N_alpha+(1+mean_Ne_N_alpha|Species))
#   sjPlot:: tab_model(predmNe_N_alpha2)
#   coef(predmNe_N_alpha2)
# predmNe_N_12 <- lmer(data=data_mean_sc$Ne_N_1,  cohort_mean_of_Mpred_per_age~1+mean_Ne_N_1+(1+mean_Ne_N_1|Species))
# sjPlot:: tab_model(predmNe_N_12)
# coef(predmNe_N_12)
#   predmNb2 <- lmer(data=data_mean_sc$Nb,  year_MpredbyY~1+mean_Nb+(1+mean_Nb|Species))
#   sjPlot:: tab_model(predmNb2)
#   coef(predmNb2)
# predmNe_temporal2 <- lmer(data=data_mean_sc$Ne_temporal, year_MpredbyY ~1+mean_Ne_temporal+year_biomass+(1+year_biomass|Species))
# sjPlot:: tab_model(predmNe_temporal2)
# coef(predmNe_temporal2)
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# htzm2 <- lmer(data=left_join(data_mean_sc$Nb,data_mean_sc$Ne_temporal),  htz_slope~1+Nb+mean_of_Mpred_per_age+generation_length+
#                 (1+Nb+mean_of_Mpred_per_age+generation_length|Species))
# coef(htzm2)
# sjPlot:: tab_model(htzm2)
# 
# 
# modelmer <- buildmer(as.formula(paste0("htz_slope", " ~",
#                         paste(colnames(data_mean_sc$Ne_N_1)[-which(colnames(data_mean_sc$Ne_N_1) %in%
#                         c("Species", "replicate", "Fmax", "Fstart", "htz_slope", "Fp"))], collapse="+"),
#                         "+ (1|Fp) + (1|replicate) + (1|Species)")), data=data_mean_sc$Ne_N_1)
# 
# 
# 
# 
# 
# mod_rs <- lmer(
#   Ne_N_1 ~ mean_of_F_per_age * mean_of_Mpred_per_age +
#     ( mean_of_F_per_age * mean_of_Mpred_per_age | Species) +
#     ( mean_of_F_per_age * mean_of_Mpred_per_age | Fmax) +
#     ( mean_of_F_per_age * mean_of_Mpred_per_age | Fstart),
#   data = data_mean_sc$Ne_N_1
# )
# coef(mod_rs)
