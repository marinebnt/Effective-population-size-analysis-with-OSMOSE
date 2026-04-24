# make sure that the Ne are coherent with expected values :
# Test multiple school size multipliers and then explore the effect on Ne & else
# 
# marine beneat
# 02/06/25
# 



main = "New_gpu_conf.R"
source(file.path("osmose-ns-v4.4.1_clim/R/functions/source.R"))
repo        = "RESTART"
multip_name = "y"
#multipliers = na.omit(as.numeric(stringr::str_remove(list.files(paste0("osmose-ns-v4.4.1_clim/output/")), paste0("output_", y))))
m = get_par(conf, "simulation.time.nyear")
multipliers=m

# function to prepare environment : launch simulation if repo missing
function_INIT1(MULTIP=m, REPO=repo)


# prepare environment for this particular plot
rep=1
incr_multiplier <- 0
dt_Ne_tempo <- c(20)




# initialisation
simulation = paste0(repo, "/output_", multip_name, m)
Ne_version = paste0("School_", repo, m)
if (!(m %in% c(1,2,3,4,5,6,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9))){simulation = paste0(m)}
cat("\n\t---> folder ", simulation, "\n")
incr_multiplier <- incr_multiplier + 1
ns1 = NULL
outputDir_Simu  = file.path(configDir, "output", simulation)
outputDir_Ne    = file.path(configDir, "output", "Ne","by_multiplier", Ne_version)
dir.create(outputDir_Ne, recursive = T)
dim_vec  <- c("abundanceByAge", "eggNumber", "mortalityByAge", "schoolNumber")
ns1 = read_osmose(path = outputDir_Simu)

k <- 0
H_abs     <- list()
for (i in species){#"AtlanticMackerel"){#
  cat("m=", m, "\n")
  cat("sp=", i, "\n")
  k = k+1 # increment
  
  # # to have cohort birth time steps 
  # cohort_time_corr <- data_Ne[["cohort_time_step_corresponding"]]
  # co_to_include <- colnames(data_Ne[["longevity"]][[i]][[paste0("rep_", rep)]])[which(colSums(data_Ne[["longevity"]][[i]][[paste0("rep_", rep)]] )>0)]
  # 
  # load useful data
  # raw_Ne_i       <- as.data.frame(data_Ne$Ne[[i]][[rep]])
  # raw_Ne_i$age   <- as.numeric(rownames(raw_Ne_i))
  # raw_ab_c_i     <- as.data.frame(data_Ne$abundanceByAge[[i]])
  # raw_ab_c_i$age <- as.numeric(rownames(raw_ab_c_i))
  # tempo_Ne       <- tempoNe[[i]]
  # tempo_Ne_time_seq <- as.numeric(colnames(tempo_Ne[[1]]))
  # tempo_Ne_dt    <- 20
  # new_age_one    <- as.numeric(data_Ne$new_first_age[i])
  
  
  # ZERO : get Nr and Nc : population abundance of the reproducing adults, or cohort abundance of reproducing adults
  #        also, estimate the harmonic means of these vectors         
  
  # load data
  # tab_Nc_Nr <- get(load(paste0(outputDir_Ne, "/abundance_mature.RData")))
  # age_mat <- floor(min(ns1$ageMature[,i,rep], na.rm = T))
  # tab_Nc_Nr <- ns1$abundanceByAge[[i]][,,rep]
  # tab_Nc_Nr <- tab_Nc_Nr[,-which(as.numeric(colnames(tab_Nc_Nr)) < age_mat)]
  
  # #Nc : by cohort
  # Nc   <-c()
  # Nc_h <- c()
  # for (co in names(table(cohort_time_corr[which(cohort_time_corr>=0, arr.ind = T)]))){
  #   coho_right_dt <- cohort_time_corr[which(rownames(cohort_time_corr) %in% rownames(tab_Nc_Nr)),]
  #   ID_coho_abund_mature <- which(coho_right_dt==co, arr.ind = T)
  #   if(length(ID_coho_abund_mature)>0) {
  #     if(any(as.numeric(ID_coho_abund_mature[,"row"]) <= nrow(tab_Nc_Nr) & 
  #            as.numeric(ID_coho_abund_mature[,"col"]) <= ncol(tab_Nc_Nr))){
  #       id_keep <- which(as.numeric(ID_coho_abund_mature[,"row"]) <= nrow(tab_Nc_Nr) & 
  #                          as.numeric(ID_coho_abund_mature[,"col"]) <= ncol(tab_Nc_Nr))
  #       Nc[co] <- sum(tab_Nc_Nr[ID_coho_abund_mature[as.numeric(id_keep),]])
  #     }else { Nc[co] <- 0}
  #   }else {Nc[co] <- 0}
  # }
  # 
  # Nc = as.numeric(colSums(raw_ab_c_i[,-which(colnames(raw_ab_c_i) == "age")]))
  # MAT_TOT <-   data.frame("Nc"=c(Nc)) #, rep(0, years-length(Nc)))

  
  # collect data of HTZ
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
  MAT_TOT          <- data_htz
  rownames(MAT_TOT) <- 0:c(years-1)
  
  
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
  
  H_abs[[k]] <- ggplot(df)+
    geom_point(aes(x=Time, y=HTZ), size=1)+
    labs(title = paste0(i, " HTZ"))+
    theme(plot.title = element_text(size=11))
}

# STORE THE PLOTS
library(gridExtra)
dir.create(paste0(configDir, "/output/Ne"))
pdf(file = paste0(configDir, "/output/Ne/rangeHTZLONG_",m,".pdf"), width = 9.3, height = 6.0)
grid.arrange(grobs = H_abs, ncol = 4, nrow = 4)
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







