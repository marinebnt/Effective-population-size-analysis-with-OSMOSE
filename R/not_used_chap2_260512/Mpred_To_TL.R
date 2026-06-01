library(tibble)
library(dplyr)
library(stringr)
library(tidyr)

main = "New_gpu_conf.R"
source("osmose-ns-v4.4.1_clim/R/functions/source.R")
path_RData <- "C:/Users/mbeneat/Documents/osmose/Chap2/osmose-ns/osmose-ns-v4.4.1_clim/output/PLOTS/Trait_variables_by_year_by_age2/RData"
list_f     <- list.files(path_RData, recursive = T, full.names = T)


if (!file.exists("osmose-ns-v4.4.1_clim/output/RDATA/mortalities.RData")){
  for (sp_f in species){
    
    # load and keep trait by age values
    file_sp  <- grep(list_f, pattern=sp_f)
    rdata_sp <- get(load(list_f[file_sp]))
    rdata_sp$REF_sp <- sp_f
    rdata_sp <- rdata_sp[-which(rdata_sp$age == 0),]
    rdata_sp <- rdata_sp[which(rdata_sp$TRAIT %in% c("MpredbyY", "FbyY", "MstarbyY", "MnatbyY", "abundanceByAge", "abundanceMature")),]
    
    # select maturation age data from different rdata
    file_year <- get(load(stringr::str_replace(list_f[file_sp], pattern = "Trait_variables_by_year_by_age2", replacement = "Trait_variables_by_year2")))
    age_mat   <- file_year[which(file_year$TRAIT == "age_Mature"),c("Species", "value", "time")]
    age_mat   <- age_mat %>%
          group_by(Species) %>%
          summarize("age_mat" = mean(value, na.rm=T))
    
    # bind traits and mature age by species
    rdata_sp <- pivot_wider(rdata_sp, names_from = "TRAIT", values_from = "value")
    rdata_sp <- full_join(rdata_sp, age_mat)

    ## MnatbyY is the Madditionnal in reality. 
    if (sp_f==species[1]){
      mortalities <- rdata_sp
    } else {
      mortalities <- rbind(mortalities, rdata_sp)
    }
  }
  save(mortalities, file="osmose-ns-v4.4.1_clim/output/RDATA/mortalities.RData")
} else {
  mortalities <- get(load(file="osmose-ns-v4.4.1_clim/output/RDATA/mortalities.RData"))
}


for (Fp in unique(rdata_sp$Fish_press)){
  
  mortalities_F <- mortalities[which(mortalities$Fish_press == Fp), ]
  
  # select only non mature individuals
  id_remove_juv <- sapply(species, FUN = function(x) {
    age_m <- mortalities_F$age_mature[which(mortalities_F$Species == x)[1]]
    ID_remove <- which((mortalities_F$Species == x) & (mortalities_F$age>age_m))
  })
  data_juv <- mortalities_F[id_remove_juv,]
  
  # sum of the abundance by species : to ponderate the mortalities by abundance
  abundanceByAge <- data_juv %>%
    group_by(Species, rep) %>%
    summarize(sum_abund = sum(abundanceByAge))
  data_juv <- full_join(data_juv, abundanceByAge)
  
  mpred < - 
  
  
}
  






last_Mpred <- data.frame(matrix(NA, 
                                ncol=length(unique(data_complete_y$Fish_press)), 
                                nrow = length(species)))
colnames(last_Mpred) <- unique(data_complete_y$Fish_press)
rownames(last_Mpred) <- species
for (Fi in unique(data_complete_y$Fish_press)){
  Fi_chr <- as.character(Fi)
  data_F <- data_complete_y[which(data_complete_y$Fish_press == Fi),]
  data_F <- data_F[-which(data_F$age == 0),]
  data_F$age   <- as.numeric(data_F$age)
  data_F$value <- as.numeric(data_F$value)
  data_F = data_F[which(data_F$TRAIT == "MpredbyY"),]
  data_F <- data_F%>%
    group_by(Species, age, TRAIT)%>%
    summarise(mean=mean(value, na.rm=T),
              sd=sd(value, na.rm=T))
  last_Mpred[, Fi_chr] <- unlist(sapply(species, FUN = function(x) {
    data_F_sp <- data_F[which(data_F$Species == x),]
    data_F_sp <- data_F_sp[-which(data_F_sp$mean == 0),]
    return(data_F_sp[which.max(data_F_sp$age),"mean"])
  }))
}

