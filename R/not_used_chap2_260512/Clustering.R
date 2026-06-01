main = "New_gpu_conf.R" # config file
source("osmose-ns-v4.4.1_clim/R/functions/source.R") # load config and key parameters and functions
source("osmose-ns-v4.4.1_clim/R/functions/functions_key_traits.R") # plots against Fmax
pathout <- paste0("osmose-ns-v4.4.1_clim/output/PLOTS/Trait_variables_by_cohort_by_age2/otherplots")
if (!dir.exists(pathout)) dir.create(pathout, recursive = T)

# load dataset
data_complete_sp = NULL
data_complete_c = NULL
base_path = "osmose-ns-v4.4.1_clim/output/PLOTS/Trait_variables_by_cohort_by_age2"
folder_name <- basename(base_path)
suffix <- gsub("Trait_variables_by_", "", folder_name)
suffix <- gsub("_", "_", suffix)  # keep underscores
save_dir <- file.path(dirname(dirname(base_path)), "RDATA")
save_file <- file.path(save_dir, paste0("summary_trait_", suffix, ".RData"))
if (!dir.exists(save_dir)){dir.create(save_dir, recursive = T)}
FORCE=F
if (!file.exists(save_file) || FORCE==T){
  data_complete_c <- data_traits_summary(base_path)
} else {
  data_complete_c = get(load(file=save_file))
}

# dimensions dataset required : [time, Sp, value, rep, Fish_press, TRAIT]
head(data_complete_c)
data_complete_c <- data_complete_c[, which(colnames(data_complete_c) %in% c("age", "time", "Species", "value", "TRAIT", "replicate", "Fish_press"))]
data_complete_c[grep(data_complete_c$TRAIT, pattern="fec"), "TRAIT"] <- "mean_of_eggNumber_per_age"
unique(data_complete_c$TRAIT)

# adding natural mortality dimension
data_complete_c[which((data_complete_c$TRAIT) == "mean_of_Mnat_per_age"), "TRAIT"] <- "mean_of_Madd_per_age"
data_complete_c <- pivot_wider(data = data_complete_c, names_from = "TRAIT", values_from = "value")
data_complete_c <- data_complete_c[,-grep(colnames(data_complete_c), pattern="mean_of_Dis_per_age")]
data_complete_c <- data_complete_c[,-grep(colnames(data_complete_c), pattern="mean_of_Z_per_age")]
data_complete_c$mean_of_Mnat_per_age <- data_complete_c$mean_of_Madd_per_age + data_complete_c$mean_of_Mpred_per_age + data_complete_c$mean_of_Mstar_per_age
data_complete_c <- pivot_longer(data = data_complete_c, names_to = "TRAIT", values_to = "value", cols=-colnames(data_complete_c)[1:4])

# changing colnames to match expectations
colnames(data_complete_c)[which(colnames(data_complete_c) == "TRAIT")] <- "variable"
data_complete_c$Fmax   <- as.numeric(unlist(stringr::str_remove(unlist(stringr::str_split(data_complete_c$Fish_press, " Fstart"))[seq(1, 2*length(data_complete_c$Fish_press), by=2)], "Fmax=")))
data_complete_c$Fstart <- as.numeric(unlist(stringr::str_split(data_complete_c$Fish_press, " Fstart="))[seq(2, 2*length(data_complete_c$Fish_press), by=2)])
data_complete_c <- (data_complete_c[,-which(colnames(data_complete_c) == "Fish_press")])
data_complete_c <- data_complete_c[-which(data_complete_c$age == 0), ]
PATHOUT  = "osmose-ns-v4.4.1_clim/output/PLOTS/KEY_TRAITS/"


# adding variables or transforming variables
data_complete_supp = NULL
for (sp in species){
  cat(sp)
  
  df_trait_sp <- data_complete_c %>% filter(Species == sp)
  df_trait_sp <- df_trait_sp %>% filter(Fstart == 1)
  
  df_trait_t    <- estimateSum(trait_data_full_i = df_trait_sp, var=paste0("mean_of_F_per_age", "\nsum", "mean_of_Mnat_per_age"))
  df_trait_t$variable    <- "sum_mortality_rates"
  df_trait_sp   <- rbind(df_trait_sp, df_trait_t)
  
  df_trait_t    <- estimateSurv(trait_data_full_i = df_trait_sp)
  df_trait_sp   <- rbind(df_trait_sp, df_trait_t)
  
  df_trait_t    <- estimateDiv(trait_data_full_i = df_trait_sp, var="mean_of_eggNumber_per_age\nDivmean_of_abundanceByAge_per_age")
  df_trait_sp   <- rbind(df_trait_sp, df_trait_t)
  
  df_trait_t    <- estimateDiv(trait_data_full_i = df_trait_sp, var="mean_of_Mpred_per_age\nDivsum_mortality_rates")
  df_trait_sp   <- rbind(df_trait_sp, df_trait_t)
  
  var_i    <- unique(df_trait_sp$variable)
  # compare all mortality rates to their value *abundance to have the number of individuals to understand where is the increase coming from
  for (v in var_i){
    if (any(grepl(pattern="mean_of_M", v) | grepl(pattern="mean_of_F", v) | grepl(pattern="sum_mortality_rates", v))) {
      df_trait_t    <- estimateMnumber(trait_data_full_i = df_trait_sp, var=paste0(v))
      df_trait_sp   <- rbind(df_trait_sp, df_trait_t)
    }
  }
  var_i    <- unique(df_trait_sp$variable)
  for (v in var_i){
    if (any(grepl(pattern="deaths", v))) {
      df_trait_t    <- estimatecolSum(trait_data_full_i = df_trait_sp, var=paste0("colSum", v))
      df_trait_sp   <- rbind(df_trait_sp, df_trait_t)
    }
  }
  
  if (is.null(data_complete_supp)) {
    data_complete_supp <- df_trait_sp
  } else(data_complete_supp   <- rbind(data_complete_supp, df_trait_sp))
}


data_w <- pivot_wider(data_complete_supp, names_from = "variable", values_from = "value")
cat("removing the NANs, don't know if good move")
#normalize
# data_w_norm <- data_w[, c("mean_of_F_per_age", "mean_of_Mpred_per_age", "mean_of_Mfor_per_age", "mean_of_Mstar_per_age")]
data_w_norm <- data_w
data_w_norm <- data_w_norm[, -which(colnames(data_w_norm) == c("colSumMpred\nDivsum_mortality_rates_deaths"))]
data_w_norm <- na.omit(data_w_norm) # remove rows with NA values
# data_w_norm[is.na(data_w_norm)] <- 0
data_w_norm[, -c(1:5)] =  scale(data_w_norm[, -c(1:5)])
data_w_norm[is.na(data_w_norm)] <- 0
# data_w[, c("mean_of_F_per_age", "mean_of_Mpred_per_age", "mean_of_Mfor_per_age", "mean_of_Mstar_per_age")] =
#   scale(data_w[, c("mean_of_F_per_age", "mean_of_Mpred_per_age", "mean_of_Mfor_per_age", "mean_of_Mstar_per_age")])
data_w_norm$`Ne/N` <- ifelse(data_w_norm$Species %in% 
                               c("Hake", "Whiting", "AtlanticHerring", "AtlanticCod"), "increasing",
                             ifelse(data_w_norm$Species %in% 
                                      c("EuropeanSprat"), "increasingCST",
                                    ifelse(data_w_norm$Species %in% 
                                             c("AtlanticMackerel", "Saithe", "CommonDab", "GreyGurnard", "Sandeel", 
                                               "CommonSole"), "decreasing",
                                           ifelse(data_w_norm$Species %in% 
                                                    c("EuropeanPlaice", "NorwayPout"), "decreasingCST", "CST"))))

# keeping only Fmax == 1 and Fstart ==1, because we are trying to characterize the species, not to see the consequences of fishing 
# data_w_norm <- data_w_norm[data_w_norm$Fmax==1, ]
# data_w_norm <- data_w_norm[data_w_norm$Fstart==1, ]
# data_w_norm <- data_w_norm[, -which(colnames(data_w_norm) %in% c("Fstart", "Fmax"))]
  
##### PCA --------------- 
#https://www.datacamp.com/fr/tutorial/pca-analysis-r?dc_referrer=https%3A%2F%2Fwww.google.com%2F
library('corrr')
library(ggcorrplot)
library("FactoMineR")
library("factoextra")
data.pca <- princomp(data_w_norm[,-c(3,17,40)])
summary(data.pca)
data.pca$loadings[, 1:2]
fviz_eig(data.pca, addlabels = TRUE)
# Graph of the variables
fviz_pca_var(data.pca, col.var = "black")
fviz_cos2(data.pca, choice = "var", axes = 1:2)
fviz_pca_var(data.pca, col.var = "cos2",
             gradient.cols = c("black", "orange", "green"),
             repel = TRUE)
# fviz_pca_biplot(data.pca, #col.var = "cos2", # trop long
#              # gradient.cols.var = c("black", "orange", "green"),
#              repel = TRUE, geom = "point", 
#              col.ind = as.factor(data_w_norm$Species))
# fviz_pca_biplot(data.pca, axes = c(1,6),  # trop long
#            # col.var = "cos2", gradient.cols.var = c("black", "orange", "green"),
#            repel = TRUE, label.val="none", habillage=data_w_norm$`Ne/N`#, 
#            # addEllipses=TRUE, ellipse.level=0.95
#            ) 
  
fviz_pca_biplot(data.pca, axes = c(4,5), 
           # col.var = "cos2", gradient.cols.var = c("black", "orange", "green"),
           repel = TRUE, label="var")+
           # addEllipses=TRUE, ellipse.level=0.95) 
  geom_point(aes(colour=data_w_norm$`Ne/N`)) + 
  scale_color_discrete("Set1")+
  stat_ellipse(aes(colour=data_w_norm$`Ne/N`))
fviz_cos2(data.pca, choice = "var", axes = c(3,4))
fviz_cos2(data.pca, choice = "var", axes = c(4,5))
fviz_cos2(data.pca, choice = "var", axes = c(5,6))
fviz_cos2(data.pca, choice = "var", axes = c(6,7))



data_plot_selected <- data_w_norm[,c(grep(colnames(data_w_norm), pattern="Mpred"), 
                                     grep(colnames(data_w_norm), pattern="Madd"),
                                     grep(colnames(data_w_norm), pattern="Mnat"),
                                     grep(colnames(data_w_norm), pattern="F"),
                                     grep(colnames(data_w_norm), pattern="Mfor"),
                                     grep(colnames(data_w_norm), pattern="biomass"),
                                     grep(colnames(data_w_norm), pattern="abundanceByAge"),
                                     grep(colnames(data_w_norm), pattern="mortality"),
                                     grep(colnames(data_w_norm), pattern="death"),
                                     grep(colnames(data_w_norm), pattern="surv"),
                                     grep(colnames(data_w_norm), pattern="egg"),
                                     which(colnames(data_w_norm) == "age"))]
data_plot_selected <- data_plot_selected[,-grep(colnames(data_plot_selected), pattern="Div")]
data_plot_selected <- data_plot_selected[,-which(duplicated(colnames(data_plot_selected)))]
data.pca <- princomp(data_plot_selected[,-9])
summary(data.pca)
data.pca$loadings[, 1:2]
fviz_eig(data.pca, addlabels = TRUE)
# Graph of the variables
fviz_pca_var(data.pca, col.var = "black")
fviz_cos2(data.pca, choice = "var", axes = 1:2)
fviz_pca_var(data.pca, col.var = "cos2",
             gradient.cols = c("black", "orange", "green"),
             repel = TRUE)
fviz_pca_biplot(data.pca, axes = c(2,1), 
                # col.var = "cos2", gradient.cols.var = c("black", "orange", "green"),
                repel = TRUE, label="var")+
  # addEllipses=TRUE, ellipse.level=0.95) 
  geom_point(aes(colour=data_w_norm$`Ne/N`)) + 
  scale_color_discrete("Set1") +
  stat_ellipse(aes(colour=data_w_norm$`Ne/N`))



##### K-MEANS ------------

# fitting and evaluating the model

set.seed(123)
data_plot_selected2 <- data_w_norm[,which(colnames(data_w_norm) %in%
                          c("Fmax", "mean_of_eggNumber_per_age", "mean_of_Mpred_per_age", 
                            "mean_of_Madd_per_age", "mean_of_Mnat_per_age", "mean_of_F_per_age", 
                            "mean_of_schoolNumber_per_age", 
                            "colSumall_deaths", "sum_mortality_rates",
                            "surv_curve",
                          "mean_of_biomassByAge_per_age"))]
km.out <- kmeans(data_plot_selected2, centers = 3, nstart = 20)
km.out

# Decide how many clusters to look at
n_clusters <- 5

# Initialize total within sum of squares error: wss
wss <- numeric(n_clusters)

set.seed(123)

# Look over 1 to n possible clusters
for (i in 1:n_clusters) {
  # Fit the model: km.out
  km.out <- kmeans(data_plot_selected2, centers = i, nstart = 20)
  # Save the within cluster sum of squares
  wss[i] <- km.out$tot.withinss
}

# Produce a scree plot
wss_df <- tibble(clusters = 1:n_clusters, wss = wss)

scree_plot +
  geom_hline(
    yintercept = wss, 
    linetype = 'dashed', 
    col = ifelse(cumsum(wss/sum(wss))>0.8, "red", "black")
  )

library(ggplot2)

df_plot <- data.frame("Species"=data_w_norm$Species, "cluster"=km.out$cluster)

ggplot(df_plot) +
  geom_point(aes(x = Species,
                   y = cluster))+
  ggtitle("Clusters by species") +
  theme(axis.text.x = element_text(angle = 90))

km.out$centers

# Clsuter 2 very well defines the group of fish with decreasing Ne/N 
#  + Haddock & Horse mackerel == all sp that do not have increasing Ne/N, except from common sole 
# (& includes european sprat that should be excluded)
# > km.out$centers
# mean_of_Mpred_per_age surv_curve mean_of_eggNumber_per_age mean_of_meanSizeByAge_per_age 
# 1            -0.1690100 -0.1027748               -0.03453899                    -0.2920973
# 2             4.2447954  2.2717659               -0.04909751                    -1.2938177
# 3            -0.3390315 -0.1096254              138.38791562                     0.1010426
# 4            -0.3399962 -0.1094883                0.17312938                     2.1251898



##### K-MEDOIDS ------------

# fitting and evaluating the model
library(ggplot2)
library(cluster)
set.seed(123)

data_plot_selected2 <- data_w_norm[, -c(1,2,3,38)] #[,which(colnames(data_w_norm) %in%
                                            # c("mean_of_eggNumber_per_age",
                                            #   "mean_of_Mpred_per_age",
                                            #   "surv_curve",
                                            #   "mean_of_meanSizeByAge_per_age"))]
data_plot_selected2 <- data_w_norm[,which(colnames(data_w_norm) %in%
                                            c("Fmax", "mean_of_eggNumber_per_age", "mean_of_Mpred_per_age", 
                                              "mean_of_Madd_per_age", "mean_of_Mnat_per_age", "mean_of_F_per_age", 
                                              "mean_of_schoolNumber_per_age", 
                                              "colSumall_deaths", "sum_mortality_rates",
                                              "surv_curve",
                                              "mean_of_biomassByAge_per_age"))]
n_clusters <- 6
kmed.out <- pam(data_plot_selected2, n_clusters)


df_plot <- data.frame("Species"=data_w_norm$Species, "cluster"=kmed.out$cluster)

ggplot(df_plot) +
  geom_point(aes(x = Species,
                 y = cluster))+
  ggtitle("Clusters by species") +
  theme(axis.text.x = element_text(angle = 90))

kmed.out$medoids







##### SUPERVISED CLASSIFICATION: Predicting Ne/N from trait variables -----

# The goal is to predict data_w_norm$`Ne/N` (the target) using all variables in data_w (non-normalized)
library(caret)
library(randomForest)
library(e1071)

# Prepare the classification dataset using non-normalized data
data_classif0 <- data_w

# Add the target variable from normalized data FIRST (before removing columns)
data_classif0$Ne_N  <- (ifelse(data_classif0$Species %in% 
                               c("Hake", "Whiting", "AtlanticHerring", "AtlanticCod"), "increasing",
                             ifelse(data_classif0$Species %in% 
                                      c("EuropeanSprat"), "increasingCST",
                                    ifelse(data_classif0$Species %in% 
                                             c("AtlanticMackerel", "Saithe", "CommonDab", "GreyGurnard", "Sandeel", 
                                               "CommonSole"), "decreasing",
                                           ifelse(data_classif0$Species %in% 
                                                    c("EuropeanPlaice", "NorwayPout"), "decreasingCST", "CST")))))

# keep only life history traits 
data_classif0 <- data_classif0[,which(colnames(data_classif0) %in% 
                                        c("mean_of_abundanceByAge_by_age", "mean_of_biomassByAge_by_age", 
                                        "mean_of_meanSizeByAge_per_age", "mean_of_Mpred_per_age", 
                                         "mean_of_eggNumber_per_age", "age", "replicate", "Fstart", "Fmax", "Species", "Ne_N"))]

# Ensure finite numeric predictors for all models
num_cols <- sapply(data_classif0, is.numeric)
cat("Columns ", colnames(data_classif0)[which(num_cols==F)], " are not numeric.\n")

# Handle missing values (keep only complete rows)
data_classif <- data_classif0[complete.cases(data_classif0), , drop = FALSE]

# Summary of target distribution
cat("Target variable (Ne/N) distribution:\n")
print(table(data_classif0$Ne_N))

set.seed(123)

# Create training and test sets (70-30 split)
trainIndex <- createDataPartition(data_classif0$Ne_N, p = 0.7, list = FALSE)
data_train0 <- data_classif0[trainIndex, ]
data_test0 <- data_classif0[-trainIndex, ]

# Remove problematic and metadata columns that shouldn't be features
data_classif <- data_classif0[, -which(colnames(data_classif0) %in% 
                                        c("age", "Species", "replicate", "Fstart", "Fmax"))]
# data_classif <- data_classif[, -which(colnames(data_classif) == "colSumMpred\nDivsum_mortality_rates_deaths"), drop = FALSE] # because all is NA
# data_classif <- data_classif[, -grep(colnames(data_classif), pattern = "Div"), drop = FALSE] # because name has \n that fails in the analysis

# Ensure finite numeric predictors for all models
num_cols <- sapply(data_classif, is.numeric)
cat("Columns ", colnames(data_classif)[which(num_cols==F)], " are not numeric.\n")

# Handle missing values (keep only complete rows)
data_classif <- data_classif[complete.cases(data_classif), , drop = FALSE]

# Summary of target distribution
cat("Target variable (Ne/N) distribution:\n")
print(table(data_classif$Ne_N))

set.seed(123)

# Create training and test sets (70-30 split)
trainIndex <- createDataPartition(data_classif$Ne_N, p = 0.7, list = FALSE)
data_train <- data_classif[trainIndex, ]
data_test <- data_classif[-trainIndex, ]

cat("\nTraining set size:", nrow(data_train), "\n")
cat("Test set size:", nrow(data_test), "\n")


# ---- Single Decision Tree Classification ----
cat("\n===== SINGLE DECISION TREE =====\n")
library(rpart)
library(rpart.plot)

# Train a single decision tree
tree_model <- rpart(Ne_N ~ ., 
                    data = data_train, 
                    method = "class",
                    control = rpart.control(maxdepth = 5, cp = 0.01))

# Predictions on test set
tree_pred <- predict(tree_model, data_test, type = "class")
tree_cm <- confusionMatrix(tree_pred, as.factor(data_test$Ne_N))
print(tree_cm)

# Plot the decision tree
rpart.plot(tree_model, main = "Single Decision Tree for Ne/N Classification")


# ---- Random Forest Classification ----
cat("\n===== RANDOM FOREST =====\n")
data_train$Ne_N <- as.factor(data_train$Ne_N)
rf_model <- randomForest(Ne_N ~ ., 
                         data = data_train, 
                         ntree = 500,
                         mtry = sqrt(ncol(data_train)-1),
                         importance = TRUE,
                         verbose = TRUE)

# Predictions on test set
rf_pred <- predict(rf_model, data_test)
rf_cm <- confusionMatrix(rf_pred, as.factor(data_test$Ne_N))
print(rf_cm)

# Feature importance
rf_importance <- importance(rf_model)
rf_importance_sorted <- rf_importance[order(rf_importance[, "MeanDecreaseGini"], decreasing = TRUE), ]
cat("\nTop 15 Most Important Features (Random Forest):\n")
print(head(rf_importance_sorted, 15))

# Visualize feature importance
par(mar = c(5, 8, 2, 2))
varImpPlot(rf_model, n.var = ncol(data_train)-1, main = "Feature Importance - Random Forest")


# ---- Support Vector Machine (SVM) Classification ----
cat("\n===== SUPPORT VECTOR MACHINE (SVM) =====\n")

# Prepare data for SVM (standardize features)
data_train_svm <- data_train
data_test_svm <- data_test

# Extract features only
X_train <- data_train[, -which(colnames(data_train) == "Ne_N")]
X_test <- data_test[, -which(colnames(data_test) == "Ne_N")]

# Standardize
preproc <- preProcess(X_train, method = c("center", "scale"))
X_train_scaled <- predict(preproc, X_train)
X_test_scaled <- predict(preproc, X_test)

# Train SVM with radial kernel
svm_model <- svm(x = X_train_scaled, 
                 y = data_train$Ne_N,
                 kernel = "radial",
                 cost = 1,
                 gamma = 0.1,
                 probability = TRUE)

# Predictions
svm_pred <- predict(svm_model, X_test_scaled)
svm_cm <- confusionMatrix(svm_pred, as.factor(data_test$Ne_N))
print(svm_cm)


# ---- Logistic Regression (baseline) - using multinomial for multiple classes ----
cat("\n===== MULTINOMIAL LOGISTIC REGRESSION (Baseline) =====\n")

# Using caret's train function for multinomial logistic regression
lr_model <- train(Ne_N ~ ., 
                  data = data_train,
                  method = "multinom",
                  trControl = trainControl(method = "cv", number = 5),
                  trace = FALSE)

# Predictions
lr_pred <- predict(lr_model, data_test)
lr_cm <- confusionMatrix(lr_pred, as.factor(data_test$Ne_N))
print(lr_cm)






# ---- GAM Multinomial Classification ----
cat("\n===== GAM MULTINOMIAL =====\n")
library(mgcv)

# mgcv::multinom requires integer classes 0..K and often works better with syntactically valid names
data_train_gam <- data_train0
data_test_gam <- data_test0

colnames(data_train_gam) <- make.names(colnames(data_train_gam), unique = TRUE)
colnames(data_test_gam) <- make.names(colnames(data_test_gam), unique = TRUE)

class_levels <- levels(as.factor(data_train_gam$Ne_N))
data_train_gam$Ne_N <- factor(data_train_gam$Ne_N, levels = class_levels)
data_test_gam$Ne_N <- factor(data_test_gam$Ne_N, levels = class_levels)

data_train_gam$Ne_N_num <- as.numeric(data_train_gam$Ne_N) - 1
data_train_gam <- data_train_gam[, -which(colnames(data_train_gam) == "Ne_N")]
K_multi <- length(class_levels) - 1

feature_cols_gam <- setdiff(colnames(data_train_gam), c("Ne_N", "Ne_N_num"))
smooth_terms <- paste0("s(", feature_cols_gam, ", bs='cs', k=", min(10, nrow(data_train_gam)-1), ")")
rhs_terms <- paste(smooth_terms, collapse = " + ")

gam_formula_list <- c(
  list(as.formula(paste("Ne_N_num ~", rhs_terms))),
  replicate(K_multi - 1, as.formula(paste("~", rhs_terms)), simplify = FALSE)
)

# Final GAM safety filter against NA/NaN/Inf in predictors
gam_num_cols <- setdiff(names(which(sapply(data_train_gam, is.numeric))), "Ne_N_num")
if (length(gam_num_cols) > 0) {
  data_train_gam[, gam_num_cols] <- lapply(data_train_gam[, gam_num_cols, drop = FALSE], function(x) {
    x[!is.finite(x)] <- NA
    x
  })
  data_test_gam[, gam_num_cols] <- lapply(data_test_gam[, gam_num_cols, drop = FALSE], function(x) {
    x[!is.finite(x)] <- NA
    x
  })
}
data_train_gam <- na.omit(data_train_gam)
data_test_gam <- na.omit(data_test_gam)
data_train_gam2 <- data_train_gam
# data_train_gam2 <- data_train_gam2[,c(1, ncol(data_train_gam2))]#+0.0001
gam_multi_model <- gam(
  Ne_N_num ~ s(mean_of_meanSizeByAge_per_age, by=interaction(Species, Fmax, age)), #+s(mean_of_Mpred_per_age) + s(mean_of_eggNumber_per_age),
  data = data_train_gam2,
  # family = multinom(K = 4),
  method = "REML"
)

gam_prob <- predict(gam_multi_model, newdata = data_test_gam, type = "response")
if (is.null(dim(gam_prob))) {
  gam_prob <- cbind(gam_prob)
}
if (ncol(gam_prob) == K_multi) {
  gam_prob <- cbind(1 - rowSums(gam_prob), gam_prob)
}

gam_pred_idx <- max.col(gam_prob, ties.method = "first")
gam_pred <- factor(class_levels[gam_pred_idx], levels = class_levels)
gam_cm <- confusionMatrix(gam_pred, as.factor(data_test_gam$Ne_N))
print(gam_cm)

library(mgcv)
library(ggplot2)
data(mtcars)
head(mtcars)
# Fit a GAM to model the relationship between mpg and hp
gam_model <- gam(mpg ~ s(hp), data = mtcars)

# Summary of the GAM
summary(gam_model)


# ---- Model Comparison ----
cat("\n===== MODEL COMPARISON =====\n")

model_comparison <- data.frame(
  Model = c("Random Forest", "SVM", "Logistic Regression", "GAM Multinomial"),
  Accuracy = c(rf_cm$overall["Accuracy"], 
               svm_cm$overall["Accuracy"],
               lr_cm$overall["Accuracy"],
               gam_cm$overall["Accuracy"]),
  Kappa = c(rf_cm$overall["Kappa"],
            svm_cm$overall["Kappa"],
            lr_cm$overall["Kappa"],
            gam_cm$overall["Kappa"]),
  row.names = NULL
)
print(model_comparison)

# Plot accuracy comparison
barplot(model_comparison$Accuracy,
        names.arg = model_comparison$Model,
        ylim = c(0, 1),
        ylab = "Accuracy",
        main = "Model Accuracy Comparison for Ne/N Prediction",
  col = c("forestgreen", "steelblue", "orange", "firebrick3"))
legend("topright", "test set", bty = "n")


# ---- Prediction on full dataset (if needed) ----
# If you want to apply the best model to the entire dataset:

best_model <- rf_model  # or choose another model
full_predictions <- predict(best_model, 
                            data_w[, -which(colnames(data_w) %in% 
                                           c("age", "time", "Species", "replicate", "Fstart", "Fmax"))])

cat("\nFull dataset predictions distribution:\n")
print(table(full_predictions))
