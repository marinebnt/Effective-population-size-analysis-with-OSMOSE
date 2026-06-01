library(tidyr)
library(dplyr)
library(ggplot2)
Time1=150
Time2=160

path <- c("osmose-ns-v4.4.1_clim/output/PLOTS/Trait_variables_by_cohort_by_age2/RDATA/Fmax1.5_Fstart1.5_spAtlanticCod/data_traits_summary.RData")
data <- get(load(path))

data$Fish_press <- as.factor(data$Fish_press)
data$replicate <- as.factor(data$replicate)
data$age <- as.numeric(data$age)
data$value <- as.numeric(data$value)

ggplot(data[which(data$TRAIT == "mean_of_Mfor_per_age" & data$replicate == 1 & data$Species == "AtlanticCod"),], 
       aes(color=replicate, group=(replicate)))+
  geom_line(aes(x=age, y=value))+
  facet_wrap(~Fish_press)


ggplot(data[which(data$TRAIT == "mean_of_Mnat_per_age" & data$replicate == 1 & data$Species == "AtlanticCod"),], 
       aes(color=replicate, group=(replicate)))+
  geom_line(aes(x=age, y=value))+
  facet_wrap(~Fish_press)
