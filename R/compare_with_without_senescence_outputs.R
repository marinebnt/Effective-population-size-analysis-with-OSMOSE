library(dplyr)


path_sc <- c("osmose-ns-v4.4.1_clim/output/PLOTS/FINAL_Ne_Htz_noSenescence/Complete_INIT2AtlanticHerring.RData")
path <- c("osmose-ns-v4.4.1_clim/output/PLOTS/FINAL_Ne_Htz2/Complete_INIT2AtlanticHerring.RData")


data_sc <- get(load(path_sc))
data <- get(load(path))

head(data)
head(data_sc)

colnames(data_sc[[1]])[3] <- "value_sc"
colnames(data_sc[[2]])[3] <- "value_sc"






full <- full_join(data[[1]], data_sc[[1]])

full$minus <- abs(full$value - full$value_sc)
View(full)

full_a <- full_join(data[[2]], data_sc[[2]])
fulla2 <- table(full_a[which(full_a$variable=="abundanceByAge"),][which(full_a[which(full_a$variable=="abundanceByAge"), "value"] == 0), "age"])
fulla2 <- names(fulla2)[which(fulla2 > quantile(fulla2, 0.5))]
fulla2 <- full_a[-which(as.numeric(full_a$age) %in% fulla2),]

full_a$minus <- abs(full_a$value - full_a$value_sc)
View(full_a)
fulla2$minus <- abs(fulla2$value - fulla2$value_sc)
View(fulla2)


table(full_a[which(full_a$minus==0),"variable"])*100/table(full_a$variable)
table(fulla2[which(fulla2$minus==0),"variable"])*100/table(fulla2$variable)


View(fulla2[which(fulla2$variable == "Numerator_N_1"),])








path_sc <- c("osmose-ns-v4.4.1_clim/output/PLOTS/FINAL_Ne_Htz_noSenescence/INIT/Fmax1_Fstart0.5_spAtlanticHerring/INIT2.RData")
path <- c("osmose-ns-v4.4.1_clim/output/PLOTS/FINAL_Ne_Htz2/INIT/Fmax0.5_Fstart1_spAtlanticHerring/INIT2.RData")


data_sc <- get(load(path_sc))
data <- get(load(path))
# 
# head(data)
# head(data_sc)
# 
# colnames(data_sc[[1]])[3] <- "value_sc"
# colnames(data_sc[[2]])[3] <- "value_sc"


data$data_Ne$AtlanticHerring$Numerator_N_alpha 
data_sc$data_Ne$AtlanticHerring$Numerator_N_alpha
range(data$data_Ne$AtlanticCod$abundanceByAge$rep_1) - range(data_sc$data_Ne$AtlanticCod$abundanceByAge$rep_20)



full <- full_join(data[[1]], data_sc[[1]])

full$minus <- abs(full$value - full$value_sc)
View(full)

full_a <- full_join(data[[2]], data_sc[[2]])
fulla2 <- table(full_a[which(full_a$variable=="abundanceByAge"),][which(full_a[which(full_a$variable=="abundanceByAge"), "value"] == 0), "age"])
fulla2 <- names(fulla2)[which(fulla2 > quantile(fulla2, 0.5))]
fulla2 <- full_a[-which(as.numeric(full_a$age) %in% fulla2),]

full_a$minus <- abs(full_a$value - full_a$value_sc)
View(full_a)
fulla2$minus <- abs(fulla2$value - fulla2$value_sc)
View(fulla2)


table(full_a[which(full_a$minus==0),"variable"])*100/table(full_a$variable)
table(fulla2[which(fulla2$minus==0),"variable"])*100/table(fulla2$variable)


View(fulla2[which(fulla2$variable == "Numerator_N_1"),])







path_sc2 <- c("osmose-ns-v4.4.1_clim/output/PLOTS/FINAL_Ne_Htz_noSenescence/INIT/INIT2_sensec2.RData")
path_sc <- c("osmose-ns-v4.4.1_clim/output/PLOTS/FINAL_Ne_Htz_noSenescence/INIT/Fmax0.5_Fstart1_spAtlanticHerring/INIT2.RData")


data_sc2 <- get(load(path_sc2))
data_sc <- get(load(path_sc))
# 
# head(data)
# head(data_sc)
# 
# colnames(data_sc[[1]])[3] <- "value_sc"
# colnames(data_sc[[2]])[3] <- "value_sc"


data_sc2$data_Ne$AtlanticHerring$abundanceByAge
data_sc$data_Ne$AtlanticHerring$abundanceByAge
View((data_sc2$data_Ne$AtlanticHerring$abundanceByAge$rep_1) - (data_sc$data_Ne$AtlanticHerring$abundanceByAge$rep_1))

full <- full_join(data_sc2[[1]], data_sc[[1]])

full$minus <- abs(full$value - full$value_sc)
View(full)

full_a <- full_join(data[[2]], data_sc[[2]])
fulla2 <- table(full_a[which(full_a$variable=="abundanceByAge"),][which(full_a[which(full_a$variable=="abundanceByAge"), "value"] == 0), "age"])
fulla2 <- names(fulla2)[which(fulla2 > quantile(fulla2, 0.5))]
fulla2 <- full_a[-which(as.numeric(full_a$age) %in% fulla2),]

full_a$minus <- abs(full_a$value - full_a$value_sc)
View(full_a)
fulla2$minus <- abs(fulla2$value - fulla2$value_sc)
View(fulla2)


table(full_a[which(full_a$minus==0),"variable"])*100/table(full_a$variable)
table(fulla2[which(fulla2$minus==0),"variable"])*100/table(fulla2$variable)


View(fulla2[which(fulla2$variable == "Numerator_N_1"),])










