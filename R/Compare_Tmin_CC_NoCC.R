library(ncdf4)
library(s2dv)
library(osmose)
chap2 <- "C:/Users/mbeneat/Documents/osmose/Chap2/osmose-ns/osmose-ns-v4.4.1_clim/physical_ns.nc"
chap2ltl <- "C:/Users/mbeneat/Documents/osmose/Chap2/osmose-ns/osmose-ns-v4.4.1_clim/LTL_ns_corrected_one.nc"
chap3 <- paste0("C:/Users/mbeneat/Documents/osmose/Chap2/osmose-ns/osmose-ns-v4.4.1_clim/Forcages_RCP/Physics_RCP85/20", 90:99, "_physic.nc")
chap3ltl <- paste0("C:/Users/mbeneat/Documents/osmose/Chap2/osmose-ns/osmose-ns-v4.4.1_clim/Forcages_RCP/LTL_RCP85/LTL/LTL_20", 90:99, ".nc")

conf <- osmose::read_osmose("C:/Users/mbeneat/Documents/osmose/Chap2/osmose-ns/osmose-ns-v4.4.1_clim/param_Eco_Bioen_neutral.R")
z_lay <- as.matrix(get_par(conf, "lay"))
rownames(z_lay) <- stringr::str_remove(rownames(z_lay), pattern = "species.zlayer.sp")
spec_name <- as.matrix(get_par(conf, "name"))
rownames(spec_name) <- stringr::str_remove(rownames(spec_name), pattern = "species.name.sp")


# TEMPERATURE

## NO CC ------------
d2 <- ncdf4::nc_open(chap2)
temp_lay_2 <- ncvar_get(d2, "T")
ncdf4::nc_close(d2)
dim(temp_lay_2)
id_time <- which(dim(temp_lay_2) == 24) #dimension of the time steps

temp_lay_2_notime <- MeanDims(data = temp_lay_2, dims = 4)
matrix_mean_T_sp <- data.frame(matrix(nrow=length(species), ncol=2))
k=0
for (sp in species){
  k=k+1
  id_sp_name <- rownames(spec_name)[which(spec_name == sp)]
  z_lay_sp   <- as.numeric(z_lay[which(rownames(z_lay) == id_sp_name)])
  if(z_lay_sp == 0) z_lay_sp <- c(d2$dim$nz$vals)
  
  temp_lay_2_sp <- temp_lay_2_notime[,,z_lay_sp]
  if(length(dim(temp_lay_2_sp))>2) temp_lay_2_sp <- MeanDims(data = temp_lay_2, dims = 3)
  matrix_mean_T_sp[k,] <- c(sp, mean(temp_lay_2_sp, na.rm=T))
}
colnames(matrix_mean_T_sp) <- c("Species", "No CC")
matrix_mean_T_sp[, "No CC"] <- as.numeric(matrix_mean_T_sp[, "No CC"])


## WITH CC ------------
array_chap3 <- array(dim = c(25,52,3,10))
k=0
for (path_3 in chap3){
  k=k+1
  d3 <- ncdf4::nc_open(path_3)
  temp_lay_3 <- ncvar_get(d3, "T")
  ncdf4::nc_close(d3)
  dim(temp_lay_3)
  id_time <- which(dim(temp_lay_3) == 24) #dimension of the time steps
  array_chap3[,,,k] <- MeanDims(data = temp_lay_3, dims = id_time)
}

temp_lay_3_notime <- MeanDims(data = array_chap3, dims = 4)

k=0
for (sp in species){
  k=k+1
  id_sp_name <- rownames(spec_name)[which(spec_name == sp)]
  z_lay_sp   <- as.numeric(z_lay[which(rownames(z_lay) == id_sp_name)])
  if(z_lay_sp == 0) z_lay_sp <- c(d3$dim$nz$vals)
  
  temp_lay_3_sp <- temp_lay_3_notime[,,z_lay_sp]
  if(length(dim(temp_lay_3_sp))>2) temp_lay_3_sp <- MeanDims(data = temp_lay_3, dims = 3)
  matrix_mean_T_sp[k,"With CC"] <-mean(temp_lay_3_sp, na.rm=T)
}
colnames(matrix_mean_T_sp)[3] <- c("With CC")
matrix_mean_T_sp[, "With CC"] <- as.numeric(matrix_mean_T_sp[, "With CC"])


# SAVE as RDATA ------------
Scenario = "Chap3"
PATHOUT  = paste0("output/PLOTS/", toupper(Scenario), "/KEY_TRAITS_MATURE_removed0inNedata/")
if (!dir.exists(PATHOUT)) dir.create(PATHOUT, recursive = TRUE)
if (!dir.exists(paste0(PATHOUT, "/Clean_for_paper"))) dir.create(paste0(PATHOUT, "/Clean_for_paper"), recursive = TRUE)
pathoutputRDATA = paste0(PATHOUT, "RDATA")

save(matrix_mean_T_sp, file=paste0(pathoutputRDATA, "/Tmean_by_species.RData"))



# OXYGEN

## NO CC ------------
d2 <- ncdf4::nc_open(chap2)
ox_lay_2 <- ncvar_get(d2, "O2")
ncdf4::nc_close(d2)
dim(ox_lay_2)
id_time <- which(dim(ox_lay_2) == 24) #dimension of the time steps

ox_lay_2_notime <- MeanDims(data = ox_lay_2, dims = 4)
matrix_mean_O2_sp <- data.frame(matrix(nrow=length(species), ncol=2))
k=0
for (sp in species){
  k=k+1
  id_sp_name <- rownames(spec_name)[which(spec_name == sp)]
  z_lay_sp   <- as.numeric(z_lay[which(rownames(z_lay) == id_sp_name)])
  if(z_lay_sp == 0) z_lay_sp <- c(d2$dim$nz$vals)
  
  ox_lay_2_sp <- ox_lay_2_notime[,,z_lay_sp]
  if(length(dim(ox_lay_2_sp))>2) ox_lay_2_sp <- MeanDims(data = ox_lay_2, dims = 3)
  matrix_mean_O2_sp[k,] <- c(sp, mean(ox_lay_2_sp, na.rm=T))
}
colnames(matrix_mean_O2_sp) <- c("Species", "No CC")
matrix_mean_O2_sp[, "No CC"] <- as.numeric(matrix_mean_O2_sp[, "No CC"])


## WITH CC ------------
array_chap3 <- array(dim = c(25,52,3,10))
k=0
for (path_3 in chap3){
  k=k+1
  d3 <- ncdf4::nc_open(path_3)
  ox_lay_3 <- ncvar_get(d3, "O2")
  ncdf4::nc_close(d3)
  dim(ox_lay_3)
  id_time <- which(dim(ox_lay_3) == 24) #dimension of the time steps
  array_chap3[,,,k] <- MeanDims(data = ox_lay_3, dims = id_time)
}

ox_lay_3_notime <- MeanDims(data = array_chap3, dims = 4)

k=0
for (sp in species){
  k=k+1
  id_sp_name <- rownames(spec_name)[which(spec_name == sp)]
  z_lay_sp   <- as.numeric(z_lay[which(rownames(z_lay) == id_sp_name)])
  if(z_lay_sp == 0) z_lay_sp <- c(d3$dim$nz$vals)
  
  ox_lay_3_sp <- ox_lay_3_notime[,,z_lay_sp]
  if(length(dim(ox_lay_3_sp))>2) ox_lay_3_sp <- MeanDims(data = ox_lay_3, dims = 3)
  matrix_mean_O2_sp[k,"With CC"] <-mean(ox_lay_3_sp, na.rm=T)
}
colnames(matrix_mean_O2_sp)[3] <- c("With CC")
matrix_mean_O2_sp[, "With CC"] <- as.numeric(matrix_mean_O2_sp[, "With CC"])


# SAVE as RDATA ------------
Scenario = "Chap3"
PATHOUT  = paste0("output/PLOTS/", toupper(Scenario), "/KEY_TRAITS_MATURE_removed0inNedata/")
if (!dir.exists(PATHOUT)) dir.create(PATHOUT, recursive = TRUE)
if (!dir.exists(paste0(PATHOUT, "/Clean_for_paper"))) dir.create(paste0(PATHOUT, "/Clean_for_paper"), recursive = TRUE)
pathoutputRDATA = paste0(PATHOUT, "RDATA")

save(matrix_mean_O2_sp, file=paste0(pathoutputRDATA, "/Oxygen_by_species.RData"))



