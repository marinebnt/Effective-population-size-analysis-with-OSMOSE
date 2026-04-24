main = "New_gpu_conf.R"
source(file.path("osmose-ns-v4.4.1_clim/R/functions/source.R"))
repo        = "Simu_LONG"
multip_name = "y"
multipliers = na.omit(as.numeric(stringr::str_remove(list.files("osmose-ns-v4.4.1_clim/output/Simu_MULTI"), "output_m")))


# function to prepare environment : launch simulation if repo missing

multipliers = c(1)
function_INIT1(MULTIP=multipliers, REPO=repo, FORCE=T)
# 
# 
# # prepare environment for this particular plot
# rep=1
# m=1
# 
# 
# cat("load simulation output")
# simulation = paste0(repo, "/output_", multip_name, m)
# Ne_version = paste0("School_", multip_name, m)
# if (!(m %in% c(1,2,3,4,5,6,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9))){simulation = paste0(m)}
# cat("\n\t---> folder ", simulation, "\n")
# ns1 = NULL
# outputDir_Simu  = file.path(configDir, "output", simulation)
# outputDir_Ne =file.path(configDir, "output", "Ne","by_multiplier", Ne_version)
# ns1 = read_osmose(path = outputDir_Simu)
# if(!dir.exists(outputDir_Ne)) {dir.create(outputDir_Ne)}
# 
# # function to prepare environment : 
# INIT2 <- function_INIT2(outputDir_Simu, outputDir_Ne, dim_vec, species, nb_rep, OUTPUT_data=ns1, conf, seeding_y, REPO=repo, TEMPORAL = F)
# data_Ne_raw  <- INIT2[[1]]
# data_Ne      <- INIT2[[2]]
# tempoNe      <- INIT2[[3]]
# 
# 
# 
# 
# 
# 
# species
# k=0
# par(mfrow=c(4,4))
# mapply(dx=data_Ne$b_rate_sc, y=species, FUN = function(dx, y) {
#   mapply(df=dx, y=y, FUN = function(df, y) {
#     print(plot(df[1:max(which(df[,"cohort_80"] != 0)),"cohort_80"]~as.numeric(rownames(df)[1:max(which(df[,"cohort_80"] != 0))]),
#           xlab="age", ylab="Birth rate scaled"),
#           title(main=y))
#     })
#   })
# 
# 
# 
# 
# 
# species
# k=0
# par(mfrow=c(4,4))
# sapply(species, FUN = function(y) {
#     print(plot(data_Ne$b_rate_sc[[y]][[1]][1:max(which(data_Ne$b_rate_sc[[y]][[1]][,"cohort_80"] != 0)),"cohort_80"]~
#                  data_Ne$meanSizeByAge[[y]][1:max(which(data_Ne$b_rate_sc[[y]][[1]][,"cohort_80"] != 0)), "cohort_80"],
#                xlab="mean length at age", ylab="Birth rate scaled"),
#           title(main=y))
#   })
# 
# 
# 
# species
# k=0
# par(mfrow=c(4,4))
# mapply(dx=data_Ne$b_rate, y=species, FUN = function(dx, y) {
#   mapply(df=dx, y=y, FUN = function(df, y) {
#     print(plot(df[1:max(which(df[,"cohort_80"] != 0)),"cohort_80"]~as.numeric(rownames(df)[1:max(which(df[,"cohort_80"] != 0))]),
#                xlab="age", ylab="Birth rate not scaled"),
#           title(main=y))
#   })
# })
# 
# 
# 
# 
# 
# species
# k=0
# par(mfrow=c(4,4))
# sapply(species, FUN = function(y) {
#   print(plot(data_Ne$b_rate[[y]][[1]][1:max(which(data_Ne$b_rate[[y]][[1]][,"cohort_80"] != 0)),"cohort_80"]~
#                data_Ne$meanSizeByAge[[y]][1:max(which(data_Ne$b_rate[[y]][[1]][,"cohort_80"] != 0)), "cohort_80"],
#              xlab="mean length at age", ylab="Birth rate not scaled"),
#         title(main=y))
# })
