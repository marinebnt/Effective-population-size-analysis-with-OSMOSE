#" plot for presentation"
library(ggplot2)
library(dplyr)
library(tidyr)

# Parameters
time <- 1:40
N0 <- 100       # initial abundance
M <- 0.01        # natural mortality rate per time step
mu <- 10         # peak fishing pressure
sigma <- 2       # spread of fishing pressure
q <- 0.2         # max fishing mortality
threshold <- 10 # time when fishing starts
F_rate <- 0.08   # fishing mortality applied after threshold

# Knife-edge fishing mortality
F1 <- ifelse(time >= threshold, F_rate, 0)
F1b<- ifelse(time >= threshold, F_rate*0.5, 0)
F1c<- ifelse(time >= threshold, F_rate*1.5, 0)

F_rate <- 0.04 

# Gaussian fishing mortality
F2 <- dnorm(time, mean = mu, sd = sigma)
F2 <- q * F2 / max(F2)   # scale so max fishing mortality = q
F2b<- 0.5*q * F2 / max(F2)   # scale so max fishing mortality = q
F2c<- 1.5*q * F2 / max(F2)   # scale so max fishing mortality = q

# Abundance dynamics with natural + fishing mortality
abund_curve <- function (Fisheries){
  N <- numeric(length(time))
  N[1] <- N0
  for(t in 2:length(time)){
    N[t] <- N[t-1] * exp(-(M + Fisheries[t]))  # exponential survival model
  }
  return(N)
}

N1 <- abund_curve(F1)
N1b <- abund_curve(F1b)
N1c <- abund_curve(F1c)
N2 <- abund_curve(F2)
N2b <- abund_curve(F2b)
N2c <- abund_curve(F2c)


# Baseline abundance without fishing
N_nofish <- N0 * exp(-M * (time - 1))

# Long dataframe for abundance
df_abundance <- tibble(
  time = time,
  `No fishing` = N_nofish,
  `With fishing` = N1
) %>%
  pivot_longer(-time, names_to = "Scenario", values_to = "Abundance")

# Data frames
df_abundance <- tibble(
  time = time,
  abundance_fished = N1,
  abundance_fished_b = N1b,
  abundance_fished_c = N1c,
  abundance_nofish = N_nofish
)

df_fishing <- tibble(
  time = time,
  fishing_pressure = F1*500, # rescale for area visibility
  fishing_pressure_b = F1b*500, # rescale for area visibility
  fishing_pressure_c = F1c*500 # rescale for area visibility
)

# Plot
B <- ggplot() +
  
  # Abundance as filled areas
  geom_area(data = df_abundance,
            aes(x = time, y = abundance_nofish, fill = "Abundance no fishing"),
            alpha = 0.1, color = NA) +
   geom_area(data = df_abundance,
             aes(x = time, y = abundance_fished_c, fill = "Abundance with fishing"),
             alpha = 0.2, color = NA) +
   geom_area(data = df_abundance,
             aes(x = time, y = abundance_fished, fill = "Abundance with fishing"),
             alpha = 0.2, color = NA) +
   geom_area(data = df_abundance,
             aes(x = time, y = abundance_fished_b, fill = "Abundance with fishing"),
             alpha = 0.2, color = NA) +
  
  # Fishing pressure as gray area
   geom_area(data = df_fishing,
             aes(x = time, y = fishing_pressure_c),
             fill = "gray80", alpha = 0.6) +
   geom_area(data = df_fishing,
             aes(x = time, y = fishing_pressure, fill="Fishing pressure"),
             alpha = 0.6) +
   geom_area(data = df_fishing,
             aes(x = time, y = fishing_pressure_b),
             fill = "gray60", alpha = 0.6) +
   geom_line(data = df_fishing, aes(x = time, y = fishing_pressure),
             color = "firebrick3", size = 1.2, alpha=0.6, linetype="dotted") +
   geom_line(data = df_fishing, aes(x = time, y = fishing_pressure_c),
             color = "firebrick4", size = 1.2, alpha=0.6, linetype="longdash") +
   geom_line(data = df_fishing, aes(x = time, y = fishing_pressure_b),
             color = "firebrick1", size = 1.2, alpha=0.6) +
  
  # Lines on top for clarity
  geom_line(data = df_abundance, aes(x = time, y = abundance_nofish),
            color = "darkgreen", size = 1.2, alpha=0.5) +
   geom_line(data = df_abundance, aes(x = time, y = abundance_fished),
             color = "firebrick3", size = 1.2, linetype="dotted") +
   geom_line(data = df_abundance, aes(x = time, y = abundance_fished_b),
             color = "firebrick1", size = 1.2) +
   geom_line(data = df_abundance, aes(x = time, y = abundance_fished_c),
             color = "firebrick4", size = 1.2, linetype="longdash") +
  
  scale_fill_manual(values = c("Abundance no fishing" = "darkgreen",
                               "Abundance with fishing" = "firebrick",
                               "Fishing pressure" = "grey70")) +
  theme_classic(base_size = 14) +
  labs(title = "Abundance decline under knife-edge fishing mortality",
       x = "Time", y = "Abundance / Pressure", fill = "") +
  scale_y_continuous(
      name = "Abundance",
      sec.axis = sec_axis(~./500,
                          name = "Fishing mortality rate"))






# Data frames
df_abundance <- tibble(
  time = time,
  abundance_fished = N2,
  abundance_fished_b = N2b,
  abundance_fished_c = N2c,
  abundance_nofish = N_nofish
)

df_fishing <- tibble(
  time = time,
  fishing_pressure = F2*300, # rescale for area visibility
  fishing_pressure_b = F2b*300, # rescale for area visibility
  fishing_pressure_c = F2c*300 # rescale for area visibility
)

# Plot
A <- ggplot() +
  # Abundance as filled areas
  geom_area(data = df_abundance,
            aes(x = time, y = abundance_nofish, fill = "Abundance no fishing"),
            alpha = 0.1, color = NA) +
   geom_area(data = df_abundance,
             aes(x = time, y = abundance_fished_c, fill = "Abundance with fishing"),
             alpha = 0.2, color = NA) +
   geom_area(data = df_abundance,
             aes(x = time, y = abundance_fished, fill = "Abundance with fishing"),
             alpha = 0.2, color = NA) +
   geom_area(data = df_abundance,
             aes(x = time, y = abundance_fished_b, fill = "Abundance with fishing"),
             alpha = 0.2, color = NA) +
  
  # Fishing pressure as gray area
   geom_area(data = df_fishing,
             aes(x = time, y = fishing_pressure_c),
             fill = "gray80", alpha = 0.6) +
   geom_area(data = df_fishing,
             aes(x = time, y = fishing_pressure),
             fill = "gray70", alpha = 0.6) +
   geom_area(data = df_fishing,
             aes(x = time, y = fishing_pressure_b),
             fill = "gray60", alpha = 0.6) +
   geom_line(data = df_fishing, aes(x = time, y = fishing_pressure, linetype="Fmax=0.04"),
             color = "firebrick3", size = 1.2, alpha=0.6) +
   geom_line(data = df_fishing, aes(x = time, y = fishing_pressure_c, linetype="Fmax=0.06"),
             color = "firebrick4", size = 1.2, alpha=0.6) +
   geom_line(data = df_fishing, aes(x = time, y = fishing_pressure_b, linetype="Fmax=0.02"),
             color = "firebrick1", size = 1.2, alpha=0.6) +
  

  
  # Lines on top for clarity
  geom_line(data = df_abundance, aes(x = time, y = abundance_nofish),
            color = "darkgreen", size = 1.2, alpha=0.5) +
   geom_line(data = df_abundance, aes(x = time, y = abundance_fished),
             color = "firebrick3", size = 1.2, linetype="dotted") +
   geom_line(data = df_abundance, aes(x = time, y = abundance_fished_b),
             color = "firebrick1", size = 1.2) +
   geom_line(data = df_abundance, aes(x = time, y = abundance_fished_c),
             color = "firebrick4", size = 1.2, linetype="longdash") +
  
  scale_fill_manual(values = c("Abundance no fishing" = "darkgreen",
                               "Abundance with fishing" = "firebrick")) +
  scale_linewidth_manual(values = c("Fmax=0.04"="dotted",
                            "Fmax=0.02"="longdash",
                            "Fmax=0.06"="solid"))+
  
  theme_classic(base_size = 14) +
  labs(title = "Abundance decline under Gaussian fishing mortality",
       x = "Time", y = "Abundance / Pressure", fill = "", linetype="") +
  scale_y_continuous(
    name = "Abundance",
    sec.axis = sec_axis(~./300,
                        name = "Fishing mortality rate"))



ggpubr::ggarrange(B, A, ncol=1, nrow=2,
                   labels = c("A", "B"), common.legend = TRUE, legend = "right")


# 
# 
# 
# # Dataframe
# df <- tibble(
#   time = time,
#   `Abundance (with fishing) F_rate=0.2` = N1,
#   `Abundance (with fishing) F_rate=0.1` = N1b,
#   `Abundance (with fishing) F_rate=0.3` = N1c,
#   `Abundance (no fishing)` = N_nofish,
#   `Fishing mortality F_rate=0.2` = F1 * (max(N0)/1),  # rescale fishing curve for plotting
#   `Fishing mortality F_rate=0.1` = F1b * (max(N0)/1),  # rescale fishing curve for plotting
#   `Fishing mortality F_rate=0.3` = F1c * (max(N0)/1)  # rescale fishing curve for plotting
# )
# 
# df_long <- df %>%
#   pivot_longer(-time, names_to = "Variable", values_to = "Value")
# 
# # Plot
# ggplot(df_long, aes(x = time, y = Value, color = Variable, linetype = Variable)) +
#   geom_area(size=1.2) +
#   scale_color_manual(values = c("Abundance (with fishing) F_rate=0.1" = "yellow2",
#                      "Abundance (with fishing) F_rate=0.2" = "orange",
#                      "Abundance (with fishing) F_rate=0.3" = "red3",
#                      "Abundance (no fishing)" = "forestgreen",
#                      "Fishing mortality F_rate=0.1" = "lightblue",  # rescale fishing curve for plotting
#                      "Fishing mortality F_rate=0.2" = "royalblue",  # rescale fishing curve for plotting
#                      "Fishing mortality F_rate=0.3" = "darkblue"  # rescale fishing curve for plotting
#   )) +
#   scale_y_continuous(
#     name = "Abundance",
#     sec.axis = sec_axis(~./(max(N0)/1),
#                         name = "Fishing mortality")
#   ) +
#   theme_classic(base_size = 14) +
#   labs(title = "Cohort decline with natural and fishing mortality",
#        x = "Time", color = "", linetype = "")
# 
# 
# 
# 
# 
# 
# 
# # Data frames
# df_abundance <- tibble(
#   time = time,
#   abundance_fished = N1,
#   abundance_nofish = N_nofish
# )
# 
# df_fishing <- tibble(
#   time = time,
#   fishing_pressure = F1 * max(N0) # rescale for area visibility
# )
# 
# ggplot(df) +
#   # Fishing pressure as gray area
#   geom_area(aes(x = time, y = ),
#             fill = "gray60", alpha = 0.5) +
#   
#   # Abundance as filled areas
#   geom_area(data = df_abundance,
#             aes(x = time, y = abundance_nofish, fill = "No fishing"),
#             alpha = 0.3, color = NA) +
#   geom_area(data = df_abundance,
#             aes(x = time, y = abundance_fished, fill = "With fishing"),
#             alpha = 0.6, color = NA) +
#   
#   # Lines on top for clarity
#   geom_line(data = df_abundance, aes(x = time, y = abundance_nofish),
#             color = "darkgreen", size = 1.2) +
#   geom_line(data = df_abundance, aes(x = time, y = abundance_fished),
#             color = "firebrick", size = 1.2) +
#   
#   scale_fill_manual(values = c("No fishing" = "darkgreen",
#                                "With fishing" = "firebrick")) +
#   theme_classic(base_size = 14) +
#   labs(title = "Cohort decline under knife-edge fishing mortality",
#        x = "Time", y = "Abundance / Pressure", fill = "") 
