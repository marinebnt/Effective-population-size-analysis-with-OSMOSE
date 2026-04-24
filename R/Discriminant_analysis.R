library(MASS)
library(tidyverse)
library(caret)
library(mda)
library(klaR)
library(ggrepel)
# library(ggord)
theme_set(theme_classic())


main = "New_gpu_conf.R" # config file
source("osmose-ns-v4.4.1_clim/R/functions/source.R") # load config and key parameters and functions
source("osmose-ns-v4.4.1_clim/R/functions/functions_key_traits.R") # plots against Fmax
pathout <- paste0("osmose-ns-v4.4.1_clim/output/PLOTS/Trait_variables_by_cohort_by_age2/otherplots")
if (!dir.exists(pathout)) dir.create(pathout, recursive = T)


# prepare dataset
# df <-as.data.frame(sapply(c("relativefecundity", "species.k.", "length2weight.condition", 
#                          "length2weight.allometric", 
#                "linf", "simulation.nschool", "bioen.maint.energy.c_m", 
#                ".bioen.maint.e.maint", 
#                "maturity.m1", "maturity.m0", "maturity.r", "predation.c", 
#                ), 
df <-as.data.frame(sapply(c("lifespan", "mortality.additional.larva.rate", "predation.ingestion", 
                            "mortality.additional.rate",
                            "egg.weight", "length2weight.allometric", 
                            "length2weight.condition", 
                            "increase.ratio", 
                            # "ratio.min", "ratio.max", 
                            "oxygen.c2", "oxygen.c1",
                            "mobilized.tp", "mobilized.e.d", "mobilized.e.mobi", 
                            "maint.e.maint", "maint.energy.c_m", 
                            "forage.k2", "forage.k1", 
                            "maturity.r", "m0", "m1", "predation.c" 
               ), 
             FUN = function(x) (unlist(osmose::get_par(conf, x)))))
rownames(df) <- get_par(conf, "species.name")[1:nrow(df)]
df$"Ne/N_1" <- ifelse(rownames(df) %in% 
                               c("Hake", "Whiting", "Saithe"), "increasing",
                      ifelse(rownames(df) %in% 
                               c("AtlanticMackerel", "CommonDab", "GreyGurnard", 
                                 "EuropeanPlaice"), "decreasing",
                      ifelse(rownames(df) %in% 
                              c("NorwayPout", "Sandeel", "CommonSole", 
                                "AtlanticHerring"), "constantVkdecreasing", "constant"))) # "HorseMackerel", "CommonSole", "Haddock"
df$"Ne/N_2" <- ifelse(rownames(df) %in% 
                        c("Hake", "Whiting", "Saithe"), "increasing",
                      ifelse(rownames(df) %in% 
                               c("AtlanticMackerel", "CommonDab", "GreyGurnard", 
                                 "EuropeanPlaice"), "decreasing",
                             ifelse(rownames(df) %in% 
                                      c("NorwayPout", "Sandeel", "EuropeanSprat", 
                                        "AtlanticHerring", "HorseMackerel", "CommonSole", "Haddock", "AtlanticCod"), "constant", "out"))) # "HorseMackerel", "CommonSole", "Haddock"


df <- df[-which(rownames(df) == "Shrimp"),] #remove the shrimp

# Split the data into training (80%) and test set (20%)
training.samples <- df$"Ne/N_1" %>%
  createDataPartition(p = 0.8, list = FALSE)
train.data <- df[training.samples, ]
test.data <- df[-training.samples, ]
# Normalize the data
preproc.param <- train.data %>% 
  preProcess(method = c("center", "scale"))
# Transform the data using the estimated parameters
train.transformed <- preproc.param %>% predict(train.data)
train.transformed <- df %>% 
  preProcess(method = c("center", "scale")) %>% 
  predict(df)
test.transformed <- preproc.param %>% predict(test.data)



# METHOD : https://www.sthda.com/english/articles/36-classification-methods-essentials/146-discriminant-analysis-essentials-in-r/
  
### Linear discriminant analysis - LDA ---------------

# The LDA algorithm starts by finding directions that maximize the separation between classes, then use these 
# directions to predict the class of individuals. These directions, called linear discriminants, are a linear 
# combinations of predictor variables.
# LDA assumes that predictors are normally distributed (Gaussian distribution) and that the different classes 
# have class-specific means and equal variance/covariance.

# Before performing LDA, consider:
# Inspecting the univariate distributions of each variable and make sure that they are normally distribute. 
# If not, you can transform them using log and root for exponential distributions and Box-Cox for skewed distributions.
# removing outliers from your data and standardize the variables to make their scale comparable.
# The linear discriminant analysis can be easily computed using the function lda() [MASS package].

# Before performing LDA, consider:
# Inspecting the univariate distributions of each variable and make sure that they are normally distribute. If not, you can transform them using log and root for exponential distributions and Box-Cox for skewed distributions.
# removing outliers from your data and standardize the variables to make their scale comparable.
# Diagramme de densité
ggplotdata <- train.data[which(sapply(train.data, is.numeric))]
ggplotdata$id <- 1:nrow(ggplotdata)
ggplotdata <- pivot_longer(ggplotdata, cols=-c("id"), 
                           values_to = "density", names_to = "variable")
ggplot(ggplotdata)+
  geom_density(aes(density))+
  facet_wrap(~variable, scales="free")
# QQ plot
ggplot(ggplotdata, aes(sample=density))+
  stat_qq()+
  stat_qq_line()+
  facet_wrap(~variable, scales="free")

# Fit the model
model <- lda(`Ne/N_1`~., data = train.transformed[, -which(colnames(train.transformed) == "Ne/N_2")])
plot(model)
# Make predictions
predictions <- model %>% predict(test.transformed[, -which(colnames(train.transformed) == "Ne/N_2")])
head(predictions$class, 6)
head(predictions$posterior, 6) 
head(predictions$x, 3) 
names(predictions)
# plot model 
lda.data <- cbind(train.transformed, predict(model)$x)
ggplot(lda.data, aes(LD1, LD2, label=rownames(lda.data))) +
  geom_point(aes(color = `Ne/N_1`), size=2)+
  geom_text(hjust=0, vjust=0)
# Model accuracy
mean(predictions$class==test.transformed$`Ne/N_1`) # 0 = bad, 1=Good

axis_expl <- model$scaling # a matrix which transforms observations to discriminant functions, normalized so that within groups covariance matrix is spherical.
pi1 <- ggplot(lda.data, aes(`LD1`, `LD2`)) +
  geom_jitter(aes(color = `Ne/N_1`), size=2)+
  geom_text_repel(hjust=0, vjust=0, max.overlaps = 15, aes(label=rownames(lda.data)))+
  geom_segment(data=axis_expl, aes(x = 0, y = 0, xend = (LD1/0.5),
                                   yend = (LD2/0.5)), 
               arrow = arrow(length = unit(1/2, "picas")),
               color="grey51") +
  geom_text_repel(data=axis_expl, aes(x = (LD1/0.5),
                                      y = (LD2/0.5), 
                                      label=rownames(axis_expl)), 
                  color="grey51", max.overlaps = 100)+
  labs(title="LDA analysis explaining the species' Ne/N variation with their traits")
ggsave(pi1, file="C:/Users/mbeneat/Documents/osmose/Chap2/osmose-ns/osmose-ns-v4.4.1_clim/output/PLOTS/LDA1.pdf")



# same but reordered the Ne/N ---

# Fit the model
model <- lda(`Ne/N_2`~., data = train.transformed[, -which(colnames(train.transformed) == "Ne/N_1")])
plot(model)
# Make predictions
predictions <- model %>% predict(test.transformed[, -which(colnames(train.transformed) == "Ne/N_2")])
head(predictions$class, 6)
head(predictions$posterior, 6) 
head(predictions$x, 3) 
names(predictions)
# plot model 
lda.data <- cbind(train.transformed, predict(model)$x)
ggplot(lda.data, aes(LD1, LD2, label=rownames(lda.data))) +
  geom_point(aes(color = `Ne/N_2`), size=2)+
  geom_text(hjust=0, vjust=0)
# Model accuracy
mean(predictions$class==test.transformed$`Ne/N_2`) # 0 = bad, 1=Good

axis_expl <- model$scaling # a matrix which transforms observations to discriminant functions, normalized so that within groups covariance matrix is spherical.
pi2 <- ggplot(lda.data, aes(`LD1`, `LD2`)) +
  geom_jitter(aes(color = `Ne/N_2`), size=2)+
  geom_text_repel(hjust=0, vjust=0, max.overlaps = 15, aes(label=rownames(lda.data)))+
  geom_segment(data=axis_expl, aes(x = 0, y = 0, xend = (LD1/0.5),
                                   yend = (LD2/0.5)), 
               arrow = arrow(length = unit(1/2, "picas")),
               color="grey51") +
  geom_text_repel(data=axis_expl, aes(x = (LD1/0.5),
                                      y = (LD2/0.5), 
                                      label=rownames(axis_expl)), 
                  color="grey51", max.overlaps = 100)+
  labs(title="LDA analysis explaining the species' Ne/N variation with their traits")
ggsave(pi2, file="C:/Users/mbeneat/Documents/osmose/Chap2/osmose-ns/osmose-ns-v4.4.1_clim/output/PLOTS/LDA2.pdf")

pi1
pi2

### Quadratic discriminant analysis - QDA ---------------


# NOT USED : for large datasets
# QDA is little bit more flexible than LDA, in the sense that it does not assumes the equality of variance/covariance. 
# In other words, for QDA the covariance matrix can be different for each class.
# LDA tends to be a better than QDA when you have a small training set. ********
# In contrast, QDA is recommended if the training set is very large, so that the variance of the classifier is not a major 
# issue, or if the assumption of a common covariance matrix for the K classes is clearly untenable (James et al. 2014).
# QDA can be computed using the R function qda() [MASS package]

### Mixture discriminant analysis - MDA ---------------

# The LDA classifier assumes that each class comes from a single normal (or Gaussian) distribution. 
# This is too restrictive.
# For MDA, there are classes, and each class is assumed to be a Gaussian mixture of subclasses, 
# where each data point has a probability of belonging to each class. Equality of covariance matrix, 
# among classes, is still assumed.
# MDA might outperform LDA and QDA is some situations, as illustrated below. In this example data, 
# we have 3 main groups of individuals, each having 3 no adjacent subgroups. 
# The solid black lines on the plot represent the decision boundaries of LDA, QDA and MDA. 
# It can be seen that the MDA classifier have identified correctly the subclasses compared to LDA and QDA, 
# which were not good at all in modeling this data.

# Fit the model
model <- mda(`Ne/N`~., data = train.transformed)
plot(model)
# Make predictions
predicted.classes <- model %>% predict(test.transformed)
# Model accuracy
mean(predicted.classes == test.transformed$`Ne/N`)
# plot model 
mda.data <- cbind(train.transformed, (model$fit$fitted.values))
ggplot(mda.data, aes(`1`, `2`, label=rownames(mda.data))) +
  geom_jitter(aes(color = `Ne/N`), size=2)+
  geom_text_repel(hjust=0, vjust=0, max.overlaps = 15)
ggplot(mda.data, aes(`1`, `2`, label=rownames(mda.data))) +
  geom_point(aes(color = `Ne/N`), size=2)+
  geom_text_repel(hjust=0, vjust=0, max.overlaps = 15)

### Flexible discriminant analysis - FDA ---------------


# FDA is a flexible extension of LDA that uses non-linear combinations of predictors such as splines. 
# FDA is useful to model multivariate non-normality or non-linear relationships among variables within each group, 
# allowing for a more accurate classification.


# Fit the model
model <- fda(`Ne/N`~., data = train.transformed)
plot(model)
model$percent.explained
model$values
model$fit$coefficients
# Make predictions
predicted.classes <- model %>% predict(test.transformed)
# Model accuracy
mean(predicted.classes == test.transformed$`Ne/N`)
# plot model 
fda.data <- cbind(train.transformed, (model$fit$fitted.values))
axis_expl <- model$fit$coefficients
colnames(axis_expl) <- c("Axis1", "Axis2")
p1 <- ggplot(fda.data, aes(`1`, `2`, label=rownames(fda.data))) +
  geom_point(aes(color = `Ne/N`), size=2)+
  geom_text(hjust=0, vjust=0)
ggsave(p1, file="C:/Users/mbeneat/Documents/osmose/Chap2/osmose-ns/osmose-ns-v4.4.1_clim/output/PLOTS/FDA.raw.pdf")
p2 <- ggplot(fda.data, aes(`1`, `2`)) +
  geom_jitter(aes(color = `Ne/N`), size=2)+
  geom_text_repel(hjust=0, vjust=0, max.overlaps = 15, aes(label=rownames(fda.data)))+
  ylim(c(-3, +3))+
  xlim(c(-3, +3))+
  geom_segment(data=axis_expl, aes(x = 0, y = 0, xend = (Axis1/1.5),
                                            yend = (Axis2/1.5)), 
               arrow = arrow(length = unit(1/2, "picas")),
               color="grey51") +
  geom_text_repel(data=axis_expl, aes(x = (Axis1/1.5),
                                y = (Axis2/1.5), 
                                label=rownames(axis_expl)), 
                  color="grey51", max.overlaps = 100)+
  labs(title="FDA analysis explaining the species' Ne/N variation with their traits")
ggsave(p2, file="C:/Users/mbeneat/Documents/osmose/Chap2/osmose-ns/osmose-ns-v4.4.1_clim/output/PLOTS/FDA.pdf")

# devtools::install_github("bedeffinian/utilities.package")

### Regularized discriminant analysis - RDA---------------

# RDA builds a classification rule by regularizing the group covariance matrices (Friedman 1989) 
# allowing a more robust model against multicollinearity in the data. This might be very useful for a 
# large multivariate data set containing highly correlated predictors.
# Regularized discriminant analysis is a kind of a trade-off between LDA and QDA. Recall that, 
# in LDA we assume equality of covariance matrix for all of the classes. QDA assumes different 
# covariance matrices for all the classes. Regularized discriminant analysis is an intermediate between LDA and QDA.
# RDA shrinks the separate covariances of QDA toward a common covariance as in LDA. 
# This improves the estimate of the covariance matrices in situations where the number of predictors 
# is larger than the number of samples in the training data, potentially leading to an improvement of the model accuracy.

# Fit the model
model <- rda(`Ne/N`~., data = train.transformed)
plot(model)
# Make predictions
predictions <- model %>% predict(test.transformed)
# Model accuracy
mean(predictions$class == test.transformed$`Ne/N`)
