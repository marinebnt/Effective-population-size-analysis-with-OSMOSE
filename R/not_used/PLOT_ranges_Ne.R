# make sure that the Ne are coherent with expected values :
# Test multiple school size multipliers and then explore the effect on Ne & else
# 
# marine beneat
# 02/06/25
# 




source(file.path("osmose-ns-v4.4.1_clim/R/functions/source.R"))
repo        = "Simu_LONG"
multip_name = "m"
multipliers = na.omit(as.numeric(stringr::str_remove(list.files("osmose-ns-v4.4.1_clim/output/Simu_MULTI"), "output_m")))


# function to prepare environment : launch simulation if repo missing
function_INIT1(MULTIP=1, REPO=repo)


# prepare environment for this particular plot
rep=1
m=1
incr_multiplier <- 0
k <- 0
dt_Ne_tempo <- c(20)




# initialisation
simulation = paste0(repo, "/output_m", m)
Ne_version = paste0("School_multiplier_", m)
if (!(m %in% c(1,2,3,4,5,6,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9))){simulation = paste0(m)}
cat("\n\t---> folder ", simulation, "\n")
incr_multiplier <- incr_multiplier + 1
ns1 = NULL
outputDir_Simu  = file.path(configDir, "output", simulation)
outputDir_Ne    = file.path(configDir, "output", "Ne","by_multiplier", Ne_version)
dir.create(outputDir_Ne)
dim_vec  <- c("abundanceByAge", "eggNumber", "mortalityByAge", "schoolNumber")
ns1 = read_osmose(path = outputDir_Simu)

# function to prepare environment : 
INIT2 <- function_INIT2(outputDir_Simu, outputDir_Ne, dim_vec, species, nb_rep, OUTPUT_data=ns1 , conf, seeding_y, REPO=repo)
data_Ne_raw  <- INIT2[[1]]
data_Ne      <- INIT2[[2]]
tempoNe      <- INIT2[[3]]




#* th : for the ne estimated acc to Agene method
#* tp : for the Ne estimated acc to Variance Ne method (Marty et al 2015)
#* h  : harmonic mean
#* Nc : cohort abundance of the mature individuals (intracohort reproductors)
#* Nr : abundance of the mature individuals of the population at each time step (intercohort reproductors)
Ne_th_abs <- list()
Ne_thbyNc  <- list()
Ne_th_h   <- list()
Ne_th_hbyNc_h <- list()
Ne_tp_abs <- list()
Ne_tpbyNr  <- list()
Ne_tpbyNr_h <- list()
H_abs     <- list()
# Nec_abs <- list()
# Nec_NT  <- list()
m1=c()
m2=c()
m3=c()
m4=c()
m5=c()
m6=c()
m7=c()
m8=c()
m9=c()
m10=c()
m11=c()


for (i in species){#"AtlanticMackerel"){#
  cat("m=", m, "\n")
  cat("sp=", i, "\n")
  k = k+1 # increment
  
  # to have cohort birth time steps 
  cohort_time_corr <- data_Ne[["cohort_time_step_corresponding"]]
  co_to_include <- colnames(data_Ne[["longevity"]][[i]][[paste0("rep_", rep)]])[which(colSums(data_Ne[["longevity"]][[i]][[paste0("rep_", rep)]] )>0)]
  
  # load useful data
  raw_Ne_i       <- as.data.frame(data_Ne$Ne[[i]][[rep]])
  raw_Ne_i$age   <- as.numeric(rownames(raw_Ne_i))
  raw_ab_c_i     <- as.data.frame(data_Ne$abundanceByAge[[i]])
  raw_ab_c_i$age <- as.numeric(rownames(raw_ab_c_i))
  tempo_Ne       <- tempoNe[[i]]
  tempo_Ne_time_seq <- as.numeric(colnames(tempo_Ne[[1]]))
  tempo_Ne_dt    <- 20
  new_age_one    <- as.numeric(data_Ne$new_first_age[i])

  
# ZERO : get Nr and Nc : population abundance of the reproducing adults, or cohort abundance of reproducing adults
#        also, estimate the harmonic means of these vectors         
  
  # load data
  # tab_Nc_Nr <- get(load(paste0(outputDir_Ne, "/abundance_mature.RData")))
  age_mat <- floor(min(ns1$ageMature[,i,rep], na.rm = T))
  tab_Nc_Nr <- ns1$abundanceByAge[[i]][,,rep]
  tab_Nc_Nr <- tab_Nc_Nr[,-which(as.numeric(colnames(tab_Nc_Nr)) < age_mat)]
  
  #Nc : by cohort
  Nc   <-c()
  Nc_h <- c()
  for (co in names(table(cohort_time_corr[which(cohort_time_corr>=0, arr.ind = T)]))){
    coho_right_dt <- cohort_time_corr[which(rownames(cohort_time_corr) %in% rownames(tab_Nc_Nr)),]
    ID_coho_abund_mature <- which(coho_right_dt==co, arr.ind = T)
    if(length(ID_coho_abund_mature)>0) {
      if(any(as.numeric(ID_coho_abund_mature[,"row"]) <= nrow(tab_Nc_Nr) & 
             as.numeric(ID_coho_abund_mature[,"col"]) <= ncol(tab_Nc_Nr))){
        id_keep <- which(as.numeric(ID_coho_abund_mature[,"row"]) <= nrow(tab_Nc_Nr) & 
                         as.numeric(ID_coho_abund_mature[,"col"]) <= ncol(tab_Nc_Nr))
        Nc[co] <- sum(tab_Nc_Nr[ID_coho_abund_mature[as.numeric(id_keep),]])
      }else { Nc[co] <- 0}
    }else {Nc[co] <- 0}
  }
  
  Nc = as.numeric(colSums(raw_ab_c_i[,-which(colnames(raw_ab_c_i) == "age")]))
  MAT_TOT <-   data.frame("Nc"=c(Nc)) #, rep(0, years-length(Nc)))
  rownames(MAT_TOT) <- 0:c(years-1)
  MAT_TOT$age <- 0:c(years-1)
  
  for (l in rownames(MAT_TOT)[1]:c(years-dt_Ne_tempo)){
    MAT_TOT[l+1, "Nc_h"] <- psych::harmonic.mean(MAT_TOT$Nc[which(rownames(MAT_TOT) %in% c(l:c(l+20)))])
  }

  # Nr : by year
  MAT_TOT[rownames(tab_Nc_Nr),"Nr"]   <- c(rowSums(tab_Nc_Nr))
  Nr_h <- c()
  for (l in rownames(MAT_TOT)[1]:c(years-dt_Ne_tempo)){
    MAT_TOT[l+1, "Nr_h"] <- psych::harmonic.mean(MAT_TOT$Nr[which(rownames(MAT_TOT) %in% c(l:c(l+20)))])
  }
  
# FIRST : estimate Ne, Ne harmonic mean and weighted Ne by abundance
  init_age     <- new_age_one
  raw_Ne       <- pivot_longer(raw_Ne_i, cols = -age, names_to = "Cohort", values_to = "Ne")
  raw_ab_c     <- pivot_longer(raw_ab_c_i, cols = -age, names_to = "Cohort", values_to = "Abundance")
  raw_ab_Ne    <- full_join(raw_Ne, raw_ab_c)
  # match the cohort age with the time step to which corresponding
  for (co in 0:c(ncohort-1)){ 
    coho_tab <- cohort_time_corr[-which(as.numeric(rownames(cohort_time_corr))<init_age),]
    ID_co <- rownames(coho_tab)[which(coho_tab==co, arr.ind = T)[,1]]
    if(length(ID_co)==0){next}
    raw_ab_Ne[grep(paste0("\\bcohort_", co, "\\b"), raw_ab_Ne$Cohort)[1:length(ID_co)], "Time"] <- ID_co
  }
  raw_ab_Ne    <- raw_ab_Ne[which(raw_ab_Ne$Cohort %in% co_to_include),]
  raw_ab_Ne_a1 <- raw_ab_Ne[which(raw_ab_Ne$age == min(raw_ab_Ne$age)),]
  MAT_TOT[raw_ab_Ne_a1$Time, "Ne_th_abs"] = raw_ab_Ne_a1$Ne
  for (l in min(rownames(MAT_TOT), na.rm = T):c(years-dt_Ne_tempo)){
    MAT_TOT$Ne_th_h[l+1] <- psych::harmonic.mean(MAT_TOT$Ne_th_abs[which(rownames(MAT_TOT) 
                                                                     %in% c(l:c(l+20)))])
  }# compare output with : aggregate(raw_ab_Ne, Ne~Time, mean) to make sure result is realistic
  
# SECOND : devide that FIRST output by the population size
  MAT_TOT$Ne_thbyNc     <- MAT_TOT$Ne_th_abs/MAT_TOT$Nc
  MAT_TOT$Ne_th_hbyNc_h <- MAT_TOT$Ne_th_abs/MAT_TOT$Nc_h
  
# THIRD : get Ne temporal
  MAT_TOT[as.numeric(names(tempo_Ne))+1, "Ne_tp_abs"] <- as.numeric(tempo_Ne)
  
# FOURTH : devide Ne tempo by Nc
  MAT_TOT$Ne_tpbyNr <- MAT_TOT$Ne_tp_abs/MAT_TOT$Nr
  MAT_TOT$Ne_tpbyNr_h <- MAT_TOT$Ne_tp_abs/MAT_TOT$Nr_h
  
  # compare output with : aggregate(raw_ab_Ne, Ne~Time, mean) to make sure result is realistic

# # FIFTH : get Ne by cohort
#   Ne_c <- raw_Ne_i[1,-ncol(raw_Ne_i)]
#   
# # SIXTH : devide Ne by the cohort total size (NT)
#   Ne_cbyNT <- Ne_c/colSums(raw_ab_c_i[,-ncol(raw_ab_c_i)])
#   NT       <- colSums(raw_ab_c_i[,-ncol(raw_ab_c_i)])

# SEVENTH : collect data of HTZ
  # if (!file.exists(paste0(outputDir_Ne, "/htz", i, ".RData"))){
        cat("saving htz data\n")
        data <- ncdf4::nc_open(paste0(outputDir_Simu, "/Genetic/ns_observedHtz-", i, "Simu0.nc"))
        data_nc <- ncdf4::ncvar_get(data, "observed_heterozygosity")
        trait1_data <- data_nc[, 1, ]
        ncdf4::nc_close(data)
        mean_htz <- colMeans(trait1_data)
        time <- 0:c(dim(data_nc)[3]-1)
        data_htz <- data.frame(age=time, HTZ=mean_htz)
        save(data_htz, file=paste0(outputDir_Ne, "/htz", i, ".RData"))
      # }  else {
        data_htz <- get(load(paste0(outputDir_Ne, "/htz", i, ".RData")))
        # cat("=> loading already existing htz info\n")}
  MAT_TOT <- dplyr::full_join(MAT_TOT, data_htz)
  
  
# COMBINE dataset  
  MAT_TOT$Time <- as.numeric(rownames(MAT_TOT))
  df <- as.data.frame(MAT_TOT[complete.cases(MAT_TOT), ])
  if(any(is.infinite(as.matrix(df)))) {df <- df[-which(df == Inf, arr.ind = T)[,1], ]}

# plot all and store plots
  mi   = max(c(df$Nc), na.rm = T)/max(df$Ne_th_abs, na.rm = T)
  mii  = max(c(df$Nc), na.rm = T)/max(df$Ne_thbyNc, na.rm = T)
  miii = max(c(df$Nc_h), na.rm = T)/max(df$Ne_th_h, na.rm = T)
  miv  = max(c(df$Nc_h), na.rm = T)/max(df$Ne_th_hbyNc_h, na.rm = T)
  mvii = max(c(df$Nr), na.rm = T)/max(df$Ne_tp_abs, na.rm = T)
  mviii= max(c(df$Nr), na.rm = T)/max(df$Ne_tpbyNr, na.rm = T)
  mix  = max(c(df$Nr_h), na.rm = T)/max(df$Ne_tpbyNr_h, na.rm = T)
  mx  = max(c(df$Ne_th_abs), na.rm = T)/max(df$H_abs, na.rm = T)
  mxi = max(c(df$Ne_tp_abs), na.rm = T)/max(df$H_abs, na.rm = T)
  # mv   = max(c(df$NT), na.rm = T)/max(df$Nec_abs, na.rm = T)
  # mvi  = max(c(df$NT), na.rm = T)/max(df$Nec_NT, na.rm = T)
  
  m1 = c(m1, mi)
  m2 = c(m2, mii)
  m3 = c(m3, miii)
  m4 = c(m4, miv)
  m7 = c(m7, mvii)
  m8 = c(m8, mviii)
  m9 = c(m9, mix)
  m10= c(m10, mx)
  m11= c(m11, mxi)
  # m5 = c(m5, mv)
  # m6 = c(m6, mvi)
  
  Ne_th_abs[[k]] <- ggplot(df)+
    geom_point(aes(x=Time, y=Ne_th_abs))+
    geom_line(aes(x=Time, y=Nc/mi))+
    labs(title = paste0("Ne_th and Nc on [50:100]"))+
    scale_y_continuous(
      # Features of the first axis
      name = "Ne_th (points)",
      # Add a second axis and specify its features
      sec.axis = sec_axis(~.*mi, name=paste0("Nc (line)"))
    )+ 
    theme(plot.title = element_text(size=12))

  Ne_thbyNc[[k]] <- ggplot(df)+
    geom_point(aes(x=Time, y=Ne_thbyNc))+
    geom_line(aes(x=Time, y=Nc/mii))+
    # geom_line(aes(x=Time, y=Ne/multip2), col="red")+
    labs(title = paste0("Ne_th/Nc and Nc on [50:100]"))+
    scale_y_continuous(
      # Features of the first axis
      name = "Ne_th/Nc (points)",
      # Add a second axis and specify its features
      sec.axis = sec_axis(~.*mii, name=paste0("Nc (line)"))
    )+ 
    theme(plot.title = element_text(size=12))
  Ne_th_h[[k]]  <- ggplot(df)+
    geom_point(aes(x=Time, y=Ne_th_h))+
    geom_line(aes(x=Time, y=Nc_h/miii))+
    # geom_line(aes(x=Time, y=Ne/multip2), col="red")+
    labs(title = paste0("h means of Ne_th and Nc on [50:100]"))+
    scale_y_continuous(
      # Features of the first axis
      name = "h mean of Ne_th (points)",
      # Add a second axis and specify its features
      sec.axis = sec_axis(~.*miii, name=paste0("h mean of Nc (line)"))
    )+ 
    theme(plot.title = element_text(size=12))
  Ne_th_hbyNc_h[[k]] <- ggplot(df)+
    geom_point(aes(x=Time, y=Ne_th_hbyNc_h))+
    geom_line(aes(x=Time, y=Nc/miv))+
    # geom_line(aes(x=Time, y=Ne/multip2), col="red")+
    labs(title = paste0("h means of Ne_th/Nc and Nc on [50:100]"))+
    scale_y_continuous(
      # Features of the first axis
      name = "h mean of Ne_th/Nc (points)",
      # Add a second axis and specify its features
      sec.axis = sec_axis(~.*miv, name=paste0("h mean of Nc (line)"))
    )+ 
    theme(plot.title = element_text(size=12))
  Ne_tp_abs[[k]] <- ggplot(df)+
    geom_point(aes(x=Time, y=Ne_tp_abs))+
    geom_line(aes(x=Time, y=Nr/mvii))+
    # geom_line(aes(x=Time, y=Ne/multip2), col="red")+
    labs(title = paste0("Ne_tp and Nr on [50:100]"))+
    scale_y_continuous(
      # Features of the first axis
      name = "Ne_tp (points)",
      # Add a second axis and specify its features
      sec.axis = sec_axis(~.*mvii, name=paste0("Nr (line)"))
    )+ 
    theme(plot.title = element_text(size=12))
  Ne_tpbyNr[[k]] <- ggplot(df)+
    geom_point(aes(x=Time, y=Ne_tpbyNr))+
    geom_line(aes(x=Time, y=Nr/mviii))+
    # geom_line(aes(x=Time, y=Ne/multip2), col="red")+
    labs(title = paste0("Ne_tp/Nr and Nr on [50:100]"))+
    scale_y_continuous(
      # Features of the first axis
      name = "Ne_tp/Nr (points)",
      # Add a second axis and specify its features
      sec.axis = sec_axis(~.*mviii, name=paste0("Nr (line)"))
    )+ 
    theme(plot.title = element_text(size=12))
  Ne_tpbyNr_h[[k]] <- ggplot(df)+
    geom_point(aes(x=Time, y=Ne_tpbyNr_h))+
    geom_line(aes(x=Time, y=Nr/mix))+
    # geom_line(aes(x=Time, y=Ne/multip2), col="red")+
    labs(title = paste0("Ne_tp/mean h Nr and mean h Nr on [50:100]"))+
    scale_y_continuous(
      # Features of the first axis
      name = "Ne_tp/mean h Nr (points)",
      # Add a second axis and specify its features
      sec.axis = sec_axis(~.*mix, name=paste0("mean h Nr (line)"))
    )+ 
    theme(plot.title = element_text(size=12))
  H_abs[[k]] <-  ggplot(df)+
    geom_point(aes(x=Time, y=HTZ))+
    # geom_line(aes(x=Time, y=Ne/multip2), col="red")+
    labs(title = paste0("HTZ"))+
    theme(plot.title = element_text(size=12))
    # ggplot(df)+
    # geom_point(aes(x=Time, y=HTZ))+
    # geom_line(aes(x=Time, y=Ne_th_abs/mx))+
    # geom_line(aes(x=Time, y=Ne_tp_abs/mxi))+
    # # geom_line(aes(x=Time, y=Ne/multip2), col="red")+
    # labs(title = paste0("SCALE ISSUE Htz, Ne_th, Ne_tp on [50:100]"))+
    # scale_y_continuous(
    #   # Features of the first axis
    #   name = "Htz (points)",
    #   # Add a second axis and specify its features
    #   sec.axis = sec_axis(~.*mx, name=paste0("Ne_tp, Ne_th (line)"))
    # )+ 
    # theme(plot.title = element_text(size=12))
}

# STORE THE PLOTS
dir.create(paste0(configDir, "/output/Ne"))
pdf(file = paste0(configDir, "/output/Ne/rangeNeLONG_",m,".pdf"), width = 9.3, height = 6.0)
# pdf(file = paste0(paste0(pathNCDF, "_FINAL_plot/", "tot", ".pdf")), width=25, height=15)
par(mfrow=c(1,1))
for (k in 1:length(Ne_th_abs)){
  i = species[k]
  mi   = m1[k]
  mii  = m2[k]
  miii = m3[k]
  miv  = m4[k]
  mvii = m7[k]
  mviii= m8[k]
  mix  = m9[k]
  mx   = m10[k]
  mxi  = m11[k]
  # mv   = m5[k]
  # mvi  = m6[k]
  print(cowplot::plot_grid(Ne_th_abs[[k]],
                           Ne_thbyNc[[k]],
                           Ne_th_h[[k]],
                           Ne_th_hbyNc_h[[k]],
                           Ne_tp_abs[[k]],
                           Ne_tpbyNr[[k]],
                           Ne_tpbyNr_h[[k]], ncol=3, nrow=3,
                           labels = paste0(i, "\n")))
}
dev.off()

plot(df$Nc~df$Nc_h)
abline(a =0, b=1)
plot(df$Nr~df$Nr_h)
abline(a =0, b=1)
plot(df$Nr~df$Nc)
abline(a =0, b=1)
plot(df$Ne_th_abs~df$Ne_th_h)
abline(a =0, b=1)
plot(df$Ne_tp_abs~df$Ne_th_abs)
abline(a =0, b=1)
plot(df$Nr_h~df$Nc_h)
abline(a =0, b=1)
plot(df$Nr_h~df$Nc_h)
abline(a =0, b=1)












#PREVIOUS VERSION

# NebyNc <- (data_Ne$Ne[[i]][[rep]]/[1,50:100]
# NetempobyNc <- (tempoNe[[i]]/rowSums(ns1$abundance)[which(rownames(ns1$abundance)%in%names(tempoNe[[i]]))])
# df <- data.frame("Ne_by_Nc" = NebyNc, "Time"=50:100, 
#                  "Nc"=data_Ne$abundanceByAge[[i]][,,1][1,50:100],
#                  "Ne" =data_Ne$Ne[[i]][[rep]][1,50:100],
#                  "Netempo_by_Nc" = NetempobyNc)
# multip = max(c(df$Nc), na.rm = T)/max(df$Ne_by_Nc, na.rm = T)
# multip2 = max(c(df$Nc), na.rm = T)/max(df$Netempo_by_Nc, na.rm = T)
# m = c(m, multip)
# m2 = c(m2, multip2)
# p[[k]] <-ggplot(df)+
#   geom_point(aes(x=Time, y=Ne_by_Nc))+
#   geom_line(aes(x=Time, y=Nc/multip))+
#   # geom_line(aes(x=Time, y=Ne/multip2), col="red")+
#   labs(title = paste0("Ne/Nc on time [50:100] for species of age 1 ", i))+
#   scale_y_continuous(
#     # Features of the first axis
#     name = "Ne/Nc (points)",
#     # Add a second axis and specify its features
#     sec.axis = sec_axis(~.*multip, name=paste0("Nc (line)"))
#   )
# 
# q[[k]]<-ggplot(df, aes(y=Netempo_by_Nc))+
#   geom_point(aes(x=Time, y=Netempo_by_Nc))+
#   geom_line(aes(x=Time, y=Nc/multip2))+
#   labs(title = paste0("Netempo/Nc on time [50:100]", i))+
#   scale_y_continuous(
#     # Features of the first axis
#     name = "Ne/Nc (points)",
#     # Add a second axis and specify its features
#     sec.axis = sec_axis(~.*multip2, name=paste0("Nc (line)"))
#   )



# Nec_abs[[k]] <- ggplot(df)+
#   geom_point(aes(x=Time, y=Nec_abs))+
#   geom_line(aes(x=Time, y=NT/mv))+
#   # geom_line(aes(x=Time, y=Ne/multip2), col="red")+
#   labs(title = paste0(i, "Ne and NT on time [50:100]"))+
#   scale_y_continuous(
#     # Features of the first axis
#     name = "Ne cohort (points)",
#     # Add a second axis and specify its features
#     sec.axis = sec_axis(~.*mv, name=paste0("NT (line)"))
#   )
# Nec_NT[[k]] <- ggplot(df)+
#   geom_point(aes(x=Time, y=Nec_NT))+
#   geom_line(aes(x=Time, y=NT/mvi))+
#   # geom_line(aes(x=Time, y=Ne/multip2), col="red")+
#   labs(title = paste0(i, "Ne of cohort/NT and NT on time [50:100]"))+
#   scale_y_continuous(
#     # Features of the first axis
#     name = "Ne/NT (points)",
#     # Add a second axis and specify its features
#     sec.axis = sec_axis(~.*mvi, name=paste0("NT (line)"))
#   )