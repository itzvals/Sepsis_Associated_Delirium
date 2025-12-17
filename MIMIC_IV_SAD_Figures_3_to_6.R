#########################################################################################################
# ZHANG PAPER REPRODUCTION
# FIGURE 3-6 
#########################################################################################################

# Clean the working directory
rm(list = ls())
# Load packages
library(glmnet)
library(caret)
library(ggplot2)
library(cowplot)
library(data.table)
library(RColorBrewer)
library(e1071)
library(xgboost)
library(randomForest)
library(rpart)
library(kknn)
library(naivebayes)
library(pROC)
library(reshape2)
library(cowplot)
library(SHAPforxgboost)
library(patchwork)
library(dplyr)
library(grid)

# Load final SAD cohort
load("sepsis_cohort_complete.RData")  

# Rename DF
sepsis_cohort <- complete_cohort 

# Verify dataset
cat(sprintf("✓ Loaded: %d patients, %d variables\n", 
            nrow(sepsis_cohort), ncol(sepsis_cohort)))

# Convert to data.table if needed
setDT(sepsis_cohort)

# Remove time related variables (related to figure 2)
sepsis_cohort <- subset(sepsis_cohort, select = -c(icu_admit_time, 
                                                   icu_discharge_time, death_time,
                                                   cam_positive, sad_group))

# View distribution
cat(sprintf("SAD distribution: %d cases (%.1f%%)\n", 
            sum(sepsis_cohort$sad_binary, na.rm = TRUE),
            mean(sepsis_cohort$sad_binary, na.rm = TRUE) * 100))

# Define what not to use as predictors
exclude_vars <- c(
  # Any known factor that contributed to SAD
  "sad_binary", "icd_delirium", "cam_status"
)

# Get all predictors: everything except excluded vars
predictor_cols <- setdiff(names(sepsis_cohort), exclude_vars)

# Check number 
cat(sprintf("Initial predictor pool: %d variables\n", length(predictor_cols)))

# Keep only numeric/binary predictors (LASSO needs numeric)
predictor_data <- sepsis_cohort[, ..predictor_cols]

numeric_cols <- names(predictor_data)[sapply(predictor_data, function(x) {
  is.numeric(x) || is.integer(x)  
})]

predictor_data <- predictor_data[, ..numeric_cols]

# Remove zero-variance columns since they can't help predicitions
predictor_data <- predictor_data[, which(sapply(predictor_data, function(x) {
  !all(is.na(x)) && (sd(x, na.rm = TRUE) > 0)
})), with = FALSE]

# Check final number
cat(sprintf("Final predictors: %d variables\n", ncol(predictor_data)))

# Add the outcome back for the modeling
data_for_modeling <- predictor_data
data_for_modeling$sad_binary <- sepsis_cohort$sad_binary

# Remove any rows with missing outcome
data_for_modeling <- data_for_modeling[!is.na(sad_binary)]

# View final dataset
cat(sprintf("\nFinal dataset: %d patients, %d predictors + 1 outcome\n", 
            nrow(data_for_modeling), ncol(data_for_modeling) - 1))

# Verify outcome is present
cat(sprintf("Outcome (sad) distribution: %d SAD (%.1f%%)\n",
            sum(data_for_modeling$sad_binary),
            mean(data_for_modeling$sad_binary) * 100))


#########################################################################################################
# TRAIN-TEST SPLIT
#########################################################################################################

set.seed(123)
index <- createDataPartition(data_for_modeling$sad_binary, p = 0.7, list = FALSE)
train <- data_for_modeling[index, ]
test <- data_for_modeling[-index, ]

cat(sprintf("\nTrain set: %d patients (%d SAD, %.1f%%)\n", 
            nrow(train), sum(train$sad_binary), mean(train$sad_binary)*100))
cat(sprintf("Test set:  %d patients (%d SAD, %.1f%%)\n\n", 
            nrow(test), sum(test$sad_binary), mean(test$sad_binary)*100))

#########################################################################################################
# LASSO REGRESSION - CODE BASED ON Zhang et al. GITHUB
# https://github.com/bbycat927/SAD/blob/main/MLcode.R
#########################################################################################################

# Prepare matrices
X_train <- as.matrix(train[, -"sad_binary"])
Y_train <- as.matrix(train$sad_binary)

X_test <- as.matrix(test[, -"sad_binary"])
Y_test <- as.matrix(test$sad_binary)

# Run cross-validated LASSO
# NOTE: glmnet will scale by default: standardize = TRUE
set.seed(123)
lasso.cv <- cv.glmnet(
  x = X_train,
  y = Y_train,
  family = "binomial",
  alpha = 1,  # LASSO (alpha=1), Ridge would be alpha=0
  nfolds = 10,
  type.measure = "deviance"
)

# Get lambda values
lambda_min <- lasso.cv$lambda.min
lambda_1se <- lasso.cv$lambda.1se

# View lambda values
cat(sprintf("Lambda min: %.6f\n", lambda_min))
cat(sprintf("Lambda 1se: %.6f\n", lambda_1se))

# Fit full LASSO model
lasso_model <- glmnet(
  x = X_train,
  y = Y_train,
  family = "binomial",
  alpha = 1
)

# Get coefficients
lasso_coef_1se <- coef(lasso.cv, s = "lambda.1se")
# Identify non-zero features
selected_features <- rownames(lasso_coef_1se)[lasso_coef_1se[,1] != 0]
# Filter to selected features 
selected_features <- selected_features[selected_features != "(Intercept)"]
# View features using lambda.1se:
print(selected_features)

#########################################################################################################
# CORSS-VALIDATION (CV) PLOT
#########################################################################################################

# Extract CV results
cv_results <- data.frame(
  lambda = lasso.cv$lambda,
  cvm = lasso.cv$cvm,  # Mean cross-validated error
  cvsd = lasso.cv$cvsd,  # Standard error
  nzero = lasso.cv$nzero  # Number of non-zero coefficients
)

cv_results$cvup <- cv_results$cvm + cv_results$cvsd
cv_results$cvlo <- cv_results$cvm - cv_results$cvsd

# Plot CV 
fig3a <- ggplot(cv_results, aes(x = log(lambda), y = cvm)) +
  geom_point(color = "red", size = 1) +
  geom_errorbar(aes(ymin = cvlo, ymax = cvup), 
                width = 0.01, color = "darkgray", alpha = 0.5) +
  geom_vline(xintercept = log(lambda_min), 
             linetype = "dashed", color = "black", size = 0.8) +
  geom_vline(xintercept = log(lambda_1se), 
             linetype = "dashed", color = "black", size = 0.8) +
  labs(
    x = "Log(λ)",
    y = "Binomial Deviance",
    title = ""
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    plot.margin = ggplot2::margin(t = 30, r = 5, b = 5, l = 5)  
  ) +
  annotate("text", x = log(lambda_min), y = max(cv_results$cvup) * 0.95,
           label = "lambda.min", hjust = -0.1, size = 3.5) +
  annotate("text", x = log(lambda_1se), y = max(cv_results$cvup) * 0.95,
           label = "lambda.1se", hjust = -0.1, size = 3.5) +
  # Add number of variables as secondary x-axis
  scale_x_continuous(
    sec.axis = sec_axis(
      ~ .,
      name = " ",
      breaks = log(cv_results$lambda)[seq(1, nrow(cv_results), length.out = 10)],
      labels = cv_results$nzero[seq(1, nrow(cv_results), length.out = 10)]
    )
  )

# View
print(fig3a)
# Save
ggsave(
  "Figure3_CV.png", fig3a, width = 12, height = 5, dpi = 300)

#########################################################################################################
# LASSO Coefficient Paths 
#########################################################################################################

# Extract coefficient paths
coef_matrix <- as.matrix(coef(lasso_model))[-1, ]  
lambda_seq <- lasso_model$lambda

# Count non-zero coefficients at each lambda
nzero_at_lambda <- colSums(coef_matrix != 0)

# Convert to long format for plotting
coef_long <- data.frame(
  lambda = rep(lambda_seq, each = nrow(coef_matrix)),
  coefficient = as.vector(coef_matrix),
  variable = rep(rownames(coef_matrix), times = length(lambda_seq))
)

# Create data frame for number of variables
nzero_df <- data.frame(
  lambda = lambda_seq,
  nzero = nzero_at_lambda
)

n_vars <- length(unique(coef_long$variable))

#########################################################################################################
# FIGURE 3
#########################################################################################################

# Chose a palette with many colors 
color_palette <- colorRampPalette(brewer.pal(12, "Paired"))(n_vars)

fig3b <- ggplot(coef_long, aes(x = log(lambda), y = coefficient, 
                               group = variable, color = variable)) +
  geom_line(alpha = 0.7, size = 0.7) +
  scale_color_manual(values = color_palette) +
  geom_vline(xintercept = log(lambda_min), 
             linetype = "dashed", color = "red", size = 0.8) +
  geom_vline(xintercept = log(lambda_1se), 
             linetype = "dashed", color = "blue", size = 0.8) +
  labs(
    x = "Log Lambda",
    y = "Coefficients",
    title = ""
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    legend.position = "none",
    plot.margin = ggplot2::margin(t = 30, r = 5, b = 5, l = 5)  
  ) +
  # Add number of variables as secondary x-axis
  scale_x_continuous(
    sec.axis = sec_axis(
      ~ .,
      name = " ",
      breaks = log(nzero_df$lambda)[seq(1, nrow(nzero_df), length.out = 10)],
      labels = nzero_df$nzero[seq(1, nrow(nzero_df), length.out = 10)]
    )
  ) +
  annotate("text", x = log(lambda_min), y = max(coef_long$coefficient) * 0.95,
           label = "lambda.min", hjust = -0.1, size = 3.5) +
  annotate("text", x = log(lambda_1se), y = max(coef_long$coefficient) * 0.95,
           label = "lambda.1se", hjust = -0.1, size = 3.5)

# View
print(fig3b)
# Save
ggsave(
  "Figure3_LASSO_Coef.png", fig3b, width = 12, height = 5, dpi = 300)

# Combine into one figure
fig3_combined <- plot_grid(
  fig3a,
  fig3b,
  labels = c("A", "B"),
  label_size = 16,
  label_fontface = "bold",
  ncol = 2
)

# View
print(fig3_combined)
# Save
ggsave(
  "Figure3_CV_and_LASSO_Coef.png", fig3_combined, width = 12, height = 5, dpi = 300)

#########################################################################################################
# EXTRACT LASSO-SELECTED FEATURES
#########################################################################################################

lasso_coef_1se <- coef(lasso.cv, s = "lambda.1se")
# Convert to regular vector
coef_vec <- as.numeric(lasso_coef_1se)
# Extract variable names
var_names <- rownames(lasso_coef_1se)

# Select features with non-zero coefficients
selected_features <- var_names[coef_vec != 0]

# Remove intercept
selected_features <- selected_features[selected_features != "(Intercept)"]

# View the selected predictors (total 35)
cat(sprintf("\n✓ LASSO selected %d features\n", length(selected_features)))
for(i in seq_along(selected_features)) {
  cat(sprintf("  %2d. %s\n", i, selected_features[i]))
}

# Create coefficient summary
selected_coefs <- coef_vec[coef_vec != 0 & var_names != "(Intercept)"]
lasso_summary <- data.frame(
  Rank = 1:length(selected_features),
  Feature = selected_features,
  Coefficient = selected_coefs,
  AbsCoef = abs(selected_coefs)
)
lasso_summary <- lasso_summary[order(-lasso_summary$AbsCoef), ]
lasso_summary$Rank <- 1:nrow(lasso_summary)

#
print(lasso_summary, row.names = FALSE)
# Save
write.csv(lasso_summary, "LASSO_Selected_Features.csv", row.names = FALSE)
cat("\n✓ Saved: LASSO_Selected_Features.csv\n")

#########################################################################################################
# PREPARE TRAIN/TEST DATA WITH SELECTED FEATURES
#########################################################################################################

# Convert back to data.frame
train <- as.data.frame(train)
test <- as.data.frame(test)

# Filter to selected features + outcome
train <- train[, c(selected_features, "sad_binary"), drop = FALSE]
test <- test[, c(selected_features, "sad_binary"), drop = FALSE]

cat(sprintf("\n✓ Train: %d rows × %d features + outcome\n", 
            nrow(train), ncol(train) - 1))
cat(sprintf("✓ Test:  %d rows × %d features + outcome\n", 
            nrow(test), ncol(test) - 1))

# NOTE: The code for the 7 models was also modeled after Zhang et al code (same repository as above)
#########################################################################################################
# MODEL 1: LOGISTIC REGRESSION
#########################################################################################################

lm_model <- glm(sad_binary ~ ., data = train, family = binomial(link = "logit"))
summary(lm_model)

lm_pred_prob <- predict(lm_model, test, type = "response")
threshold <- 0.5
predictions_binary <- ifelse(lm_pred_prob > threshold, 1, 0)

confusionMatrix(as.factor(predictions_binary), as.factor(test$sad_binary),positive = "1")

#########################################################################################################
# MODEL 2: SUPPORT VECTOR MACHINE
#########################################################################################################

# Convert to factor for SVM
train$sad_binary <- as.factor(train$sad_binary)
test$sad_binary <- as.factor(test$sad_binary)

svm_model <- svm(sad_binary ~ ., data = train, probability = TRUE)
svm_pred <- predict(svm_model, test, probability = TRUE)
svm_pred_prob <- attr(svm_pred, "probabilities")[, "1"]

confusionMatrix(data = factor(svm_pred, levels = c("0", "1")), reference = test$sad_binary, positive = "1")

#########################################################################################################
# MODEL 3: XGBOOST
#########################################################################################################

# Convert back to numeric for XGBoost
train$sad_binary <- as.numeric(as.character(train$sad_binary)) 
test$sad_binary <- as.numeric(as.character(test$sad_binary)) 

train_matrix <- xgb.DMatrix(data.matrix(train[, !names(train) %in% "sad_binary"]), label = train$sad_binary)
test_matrix <- xgb.DMatrix(data.matrix(test[, !names(test) %in% "sad_binary"]), label = test$sad_binary)

params <- list(
  objective = "binary:logistic",
  eval_metric = "logloss",
  max_depth = 3, 
  eta = 0.1, 
  gamma = 0.5, 
  colsample_bytree = 1, 
  min_child_weight = 1,
  subsample = 0.5
)

watchlist <- list(train = train_matrix, val = test_matrix)

xgb_model <- xgb.train(
  params = params,
  data = train_matrix,
  nrounds = 125,
  watchlist = watchlist,
  early_stopping_rounds = 10, 
  print_every_n = 10,
  maximize = FALSE
)

xgb_pred_prob <- predict(xgb_model, test_matrix)
xgb_pred <- ifelse(xgb_pred_prob > 0.5, 1, 0)
xgb_pred_factor <- factor(xgb_pred, levels = c(0, 1))
test_sad_factor <- factor(test$sad_binary, levels = c(0, 1))

confusionMatrix(data = xgb_pred_factor, reference = test_sad_factor, positive = "1")
saveRDS(xgb_model, "zhang_xgboost_model.rds")

#########################################################################################################
# MODEL 4: RANDOM FOREST
#########################################################################################################

# Convert to factor for RF
train$sad_binary <- as.factor(train$sad_binary)
test$sad_binary <- as.factor(test$sad_binary)

# Clean column names for RF
names(train) <- make.names(names(train))
names(test) <- make.names(names(test))

rf_model <- randomForest(sad_binary ~ ., data = train, ntree = 500, mtry = 6)
rf_pred <- predict(rf_model, newdata = test)
rf_pred_prob <- predict(rf_model, newdata = test, type = "prob")[, "1"]
rf_pred_binary <- ifelse(rf_pred_prob > 0.5, 1, 0)

confusionMatrix(data = rf_pred, reference = test$sad_binary, positive = "1")

#########################################################################################################
# MODEL 5: DECISION TREE
#########################################################################################################

dt_model <- rpart(sad_binary ~ ., data = train, method = "class")
dt_pred_prob <- predict(dt_model, newdata = test, type = "prob")[, "1"]
dt_pred <- ifelse(dt_pred_prob > 0.5, 1, 0)

confusionMatrix(factor(dt_pred, levels = c("0", "1")), test$sad_binary, positive = "1")

#########################################################################################################
# MODEL 6: NAIVE BAYES
#########################################################################################################

nb_model <- naiveBayes(sad_binary ~ ., data = train)
nb_pred_prob <- predict(nb_model, newdata = test, type = "raw")[, "1"]
nb_pred <- ifelse(nb_pred_prob > 0.5, 1, 0)

confusionMatrix(factor(nb_pred, levels = c("0", "1")), test$sad_binary, positive = "1")

#########################################################################################################
# MODEL 7: K-NEAREST NEIGHBORS
#########################################################################################################

knn_model <- kknn(sad_binary ~ ., train, test, k = 10, distance = 2, kernel = "rectangular")
knn_pred_prob <- predict(knn_model, newdata = test, type = "prob") 
knn_pred_prob <- knn_pred_prob[, "1"]
knn_pred_prob <- as.numeric(knn_pred_prob)
threshold <- 0.5 
knn_pred <- ifelse(knn_pred_prob > threshold, 1, 0)

confusionMatrix(factor(knn_pred, levels = c("0", "1")), test$sad_binary, positive = "1")

#########################################################################################################
# FIGURE 4A: ROC CURVES TRAINING COHORT
#########################################################################################################

# Build ROC list 
roc_list <- list(
  LR = roc(test$sad_binary, lm_pred_prob),
  SVM = roc(test$sad_binary, svm_pred_prob),
  XGBoost = roc(test$sad_binary, xgb_pred_prob),
  RF = roc(as.numeric(test$sad_binary), rf_pred_prob),
  KNN = roc(as.numeric(test$sad_binary), knn_pred_prob),
  DT = roc(as.numeric(test$sad_binary), dt_pred_prob),
  NB = roc(as.numeric(test$sad_binary), nb_pred_prob)
)

auc_values <- sapply(roc_list, function(x) round(auc(x), 3))
print(auc_values)

plot_colors <- brewer.pal(length(roc_list), "Set1")

# Build ROC plot data
roc_df <- do.call(rbind, lapply(seq_along(roc_list), function(i) {
  data.frame(
    Model = names(roc_list)[i],
    TPR = roc_list[[i]]$sensitivities,
    FPR = 1 - roc_list[[i]]$specificities
  )
}))


#########################################################################################################
# FIGURE 4: ROC CURVES
#########################################################################################################

fig4A <- ggplot(roc_df, aes(x = FPR, y = TPR, color = Model)) +
  geom_line(size = 1) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  scale_color_manual(values = plot_colors) +
  theme_bw() +
  labs(
    title = "Figure 4A: ROC Curves (Internal Validation)",
    x = "False Positive Rate",
    y = "True Positive Rate"
  ) +
  annotate("text", x = 0.6, y = seq(0.4, 0.1, length.out = length(auc_values)),
           label = paste(names(auc_values), "AUC=", auc_values),
           color = plot_colors, size = 3.5, hjust = 0)

# Convert test$sad_binary to numeric once for consistency
test_numeric <- as.numeric(as.character(test$sad_binary))

# Create ROC objects for all models
roc_list <- list(
  LR = roc(test_numeric, lm_pred_prob, quiet = TRUE),
  SVM = roc(test_numeric, svm_pred_prob, quiet = TRUE),
  XGBoost = roc(test_numeric, xgb_pred_prob, quiet = TRUE),
  RF = roc(test_numeric, rf_pred_prob, quiet = TRUE), 
  DT = roc(test_numeric, dt_pred_prob, quiet = TRUE),
  NB = roc(test_numeric, nb_pred_prob, quiet = TRUE),
  KNN = roc(test_numeric, knn_pred_prob, quiet = TRUE)
)

# Extract AUC values
auc_values <- sapply(roc_list, function(x) auc(x))

# View
print(round(auc_values, 4))

# Prepare data for ggplot
roc_data <- data.frame()
for(model_name in names(roc_list)) {
  roc_obj <- roc_list[[model_name]]
  
  temp_df <- data.frame(
    Model = rep(model_name, length(roc_obj$specificities)),
    Specificity = roc_obj$specificities,
    Sensitivity = roc_obj$sensitivities,
    AUC = rep(auc(roc_obj), length(roc_obj$specificities))
  )
  roc_data <- rbind(roc_data, temp_df)
}

# Define colors
model_colors <- c(
  "LR" = "#4DAF4A",
  "SVM" = "#377EB8",
  "XGBoost" = "#E41A1C",
  "RF" = "#984EA3",
  "DT" = "#FF7F00",
  "NB" = "#FFFF33",
  "KNN" = "#A65628"
)

# Create ROC plot with custom legend labels
auc_labels <- sapply(names(roc_list), function(model) {
  sprintf("%s (AUC=%.3f)", model, auc(roc_list[[model]]))
})

# Create the base plot w/out label - done to add the letter 
fig4A_base <- ggplot(roc_data, aes(x = 1 - Specificity, y = Sensitivity, 
                                   color = Model, group = Model)) +
  geom_line(linewidth = 1.2, alpha = 0.8) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", 
              color = "black", linewidth = 0.8) +
  scale_color_manual(
    values = model_colors,
    labels = auc_labels,  
    name = ""  
  ) +
  labs(
    x = "False Positive Rate",
    y = "True Positive Rate"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    legend.position = c(0.85, 0.25),
    legend.background = element_rect(fill = "white", color = "white", linewidth = 0.5),
    legend.title = element_blank(),  
    legend.text = element_text(size = 9),
    legend.key.size = unit(0.8, "cm"),  
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    axis.line = element_blank()
  ) +
  coord_equal()

# Add "A" letter
fig4A <- ggdraw() +
  draw_plot(fig4A_base) +
  draw_label("A", x = 0.02, y = 0.98, fontface = "bold", size = 16, hjust = 0, vjust = 1)

# View
print(fig4A)
#Save
ggsave("Figure4_ROC_Curves.png", fig4A, width = 8, height = 8, dpi = 300)


#########################################################################################################
# FIGURE 4B: EXTERNAL VALIDATION (MIMIC → eICU) - FIXED WITH FEATURE MAPPING
#########################################################################################################

# Load eICU data
if (file.exists("sepsis_cohort_eicu_complete.RData")) {
  
  load("sepsis_cohort_eicu_complete.RData")
  eicu_data <- complete_cohort
  setDT(eicu_data)
  
  cat(sprintf(" eICU loaded: %d patients, %d variables\n", 
              nrow(eicu_data), ncol(eicu_data)))
  
  # FEATURE MAPPING: MIMIC → eICU
  # Define feature name mappings (MIMIC name → eICU name)
  feature_map <- list(
    "anchor_age" = "age_numeric",
    "hemoglobin" = NULL,        # Not available in eICU
    "platelets" = NULL,         # Not available in eICU
    "magnesium" = NULL,         # Not available in eICU
    "phosphate" = NULL,         # Not available in eICU
    "gcs_score" = NULL,         # Not available in eICU
    "mbp" = NULL,               # Not available in eICU
    "spo2" = NULL,              # Not available in eICU
    "admit_elective" = NULL,    # Not available in eICU
    "admit_observation" = NULL, # Not available in eICU
    "admit_surgical" = NULL     # Not available in eICU
  )
  
  # Create mapping function
  map_feature_name <- function(mimic_name) {
    if (mimic_name %in% names(feature_map)) {
      eicu_name <- feature_map[[mimic_name]]
      if (is.null(eicu_name)) {
        return(NULL)  
      }
      return(eicu_name)
    }
    return(mimic_name)  
  }
  
  # Map all selected features
  eicu_feature_names <- sapply(selected_features, map_feature_name)
  
  # Remove NULL 
  available_mimic_features <- selected_features[!sapply(eicu_feature_names, is.null)]
  available_eicu_features <- unlist(eicu_feature_names[!sapply(eicu_feature_names, is.null)])
  
  cat(sprintf("\nFeature mapping results:\n"))
  cat(sprintf("  MIMIC features: %d\n", length(selected_features)))
  cat(sprintf("  Available in eICU: %d\n", length(available_eicu_features)))
  cat(sprintf("  Missing: %d\n", length(selected_features) - length(available_eicu_features)))
  
  if (length(available_eicu_features) < length(selected_features)) {
    missing_features <- setdiff(selected_features, available_mimic_features)
    cat("\n Missing features in eICU:\n")
    for (feat in missing_features) {
      cat(sprintf("  - %s\n", feat))
    }
  }
  
  # PREPARE eICU TEST SET

  # Remove variables not needed
  vars_to_remove <- c(
    "unitdischargeoffset", "cam_positive", "sad_group",
    "icd_delirium", "cam_status", "nursing_delirium", "age"
  )
  vars_to_remove <- intersect(vars_to_remove, names(eicu_data))
  
  eicu_clean <- eicu_data[, !vars_to_remove, with = FALSE]
  
  # Select only available features + outcome
  eicu_test <- as.data.frame(eicu_clean[, c(available_eicu_features, "sad_binary"), with = FALSE])
  eicu_test <- eicu_test[!is.na(eicu_test$sad_binary), ]
  
  # Rename eICU columns to match MIMIC model expectations
  names(eicu_test)[match(available_eicu_features, names(eicu_test))] <- available_mimic_features
  
  # Clean column names for consistency
  names(eicu_test) <- make.names(names(eicu_test))
  
  cat(sprintf("\neICU test set: %d patients (%d SAD, %.1f%%)\n",
              nrow(eicu_test),
              sum(eicu_test$sad_binary),
              mean(eicu_test$sad_binary) * 100))
  
  cat(sprintf("Features: %d\n\n", ncol(eicu_test) - 1))
  
  # HANDLE MISSING FEATURES (Impute with median/mode from MIMIC training)

  # Get features that MIMIC expects but eICU doesn't have
  missing_in_eicu <- setdiff(selected_features, available_mimic_features)
  
  if (length(missing_in_eicu) > 0) {
    cat(sprintf("\n Imputing %d missing features with MIMIC training set values:\n", 
                length(missing_in_eicu)))
    
    # Get MIMIC training data values for imputation
    train_for_impute <- as.data.frame(train)
    
    for (feat in missing_in_eicu) {
      if (feat %in% names(train_for_impute)) {
        if (is.numeric(train_for_impute[[feat]])) {
          # Numeric: use median
          impute_val <- median(train_for_impute[[feat]], na.rm = TRUE)
          eicu_test[[feat]] <- impute_val
          cat(sprintf("  - %s: %.2f (median)\n", feat, impute_val))
        } else {
          # Categorical: use mode
          impute_val <- names(sort(table(train_for_impute[[feat]]), decreasing = TRUE))[1]
          eicu_test[[feat]] <- as.integer(impute_val)
          cat(sprintf("  - %s: %s (mode)\n", feat, impute_val))
        }
      }
    }
    cat("\n")
  }
  
  # VERIFY FEATURE ALIGNMENT
  
  # Get expected features from each model
  train_features <- setdiff(names(train), "sad_binary")
  test_features <- setdiff(names(eicu_test), "sad_binary")
  
  cat(sprintf("MIMIC train features: %d\n", length(train_features)))
  cat(sprintf("eICU test features:   %d\n", length(test_features)))
  
  # Ensure eICU has all features MIMIC expects (in same order)
  for (feat in train_features) {
    if (!feat %in% test_features) {
      cat(sprintf("  Adding missing feature: %s\n", feat))
      # Add with median from training
      if (is.numeric(train[[feat]])) {
        eicu_test[[feat]] <- median(train[[feat]], na.rm = TRUE)
      } else {
        eicu_test[[feat]] <- 0
      }
    }
  }
  
  # Reorder columns to match training data
  eicu_test <- eicu_test[, c(train_features, "sad_binary")]
  
  cat(" Feature alignment complete\n\n")
  
  # PREDICT ON eICU WITH MIMIC MODELS
  
  # LR Predictions
  lr_pred_prob_eicu <- predict(lm_model, newdata = eicu_test, type = "response")
  cat(" LR predictions\n")
  
  # SVM Predictions
  eicu_test$sad_binary <- as.factor(eicu_test$sad_binary)
  svm_pred_eicu <- predict(svm_model, newdata = eicu_test, probability = TRUE)
  svm_pred_prob_eicu <- attr(svm_pred_eicu, "probabilities")[, "1"]
  cat(" SVM predictions\n")
  
  # XGBoost Predictions
  eicu_test$sad_binary <- as.numeric(as.character(eicu_test$sad_binary))
  test_matrix_eicu <- xgb.DMatrix(
    data.matrix(eicu_test[, !names(eicu_test) %in% "sad_binary"])
  )
  xgb_pred_prob_eicu <- predict(xgb_model, test_matrix_eicu)
  cat(" XGBoost predictions\n")
  
  # RF Predictions 
  eicu_test$sad_binary <- as.factor(eicu_test$sad_binary)
  rf_pred_prob_raw <- predict(rf_model, newdata = eicu_test, type = "prob")
  # Check column names and get the probability for class "1"
  if ("1" %in% colnames(rf_pred_prob_raw)) {
    rf_pred_prob_eicu <- rf_pred_prob_raw[, "1"]
  } else if (ncol(rf_pred_prob_raw) >= 2) {
    rf_pred_prob_eicu <- rf_pred_prob_raw[, 2]  
  } else {
    stop("Cannot extract RF probabilities - unexpected format")
  }
  cat(" RF predictions\n")
  
  # DT Predictions 
  dt_pred_prob_raw <- predict(dt_model, newdata = eicu_test, type = "prob")
  # Check if it's a matrix or vector
  if (is.matrix(dt_pred_prob_raw) || is.data.frame(dt_pred_prob_raw)) {
    if ("1" %in% colnames(dt_pred_prob_raw)) {
      dt_pred_prob_eicu <- dt_pred_prob_raw[, "1"]
    } else if (ncol(dt_pred_prob_raw) >= 2) {
      dt_pred_prob_eicu <- dt_pred_prob_raw[, 2]
    } else {
      dt_pred_prob_eicu <- as.numeric(dt_pred_prob_raw)
    }
  } else {
    dt_pred_prob_eicu <- as.numeric(dt_pred_prob_raw)
  }
  cat(" DT predictions\n")
  
  # NB Predictions
  setDT(eicu_test)
  
  nb_pred_prob_raw <- tryCatch({
    # Try naivebayes package style (type = "prob")
    predict(nb_model, newdata = eicu_test, type = "prob")
  }, error = function(e) {
    # Fallback to e1071 style or class prediction
    predict(nb_model, newdata = eicu_test)
  })
  
  # Extract probabilities
  if (is.matrix(nb_pred_prob_raw) || is.data.frame(nb_pred_prob_raw)) {
    # It's a probability matrix
    if ("1" %in% colnames(nb_pred_prob_raw)) {
      nb_pred_prob_eicu <- nb_pred_prob_raw[, "1"]
    } else if (ncol(nb_pred_prob_raw) >= 2) {
      nb_pred_prob_eicu <- nb_pred_prob_raw[, 2]
    } else {
      nb_pred_prob_eicu <- as.numeric(nb_pred_prob_raw)
    }
  } else {
    # It's class predictions - convert to probabilities (0 or 1)
    nb_pred_prob_eicu <- as.numeric(as.character(nb_pred_prob_raw))
  }
  
  cat(" NB predictions\n")
  
  # KNN Predictions 
  eicu_test_knn <- as.data.frame(eicu_test)
  train_knn <- as.data.frame(train)
  
  # Ensure both have same factor levels
  eicu_test_knn$sad_binary <- factor(eicu_test_knn$sad_binary, levels = c("0", "1"))
  train_knn$sad_binary <- factor(train_knn$sad_binary, levels = c("0", "1"))
  
  knn_model_eicu <- kknn(
    sad_binary ~ .,
    train = train_knn,
    test = eicu_test_knn,
    k = 10,
    distance = 2,
    kernel = "rectangular"
  )
  
  # Get probabilities for class "1"
  if ("1" %in% colnames(knn_model_eicu$prob)) {
    knn_pred_prob_eicu <- as.numeric(knn_model_eicu$prob[, "1"])
  } else if (ncol(knn_model_eicu$prob) >= 2) {
    knn_pred_prob_eicu <- as.numeric(knn_model_eicu$prob[, 2])
  } else {
    knn_pred_prob_eicu <- as.numeric(knn_model_eicu$prob)
  }
  cat(" KNN predictions\n\n")

  # CALCULATE ROC CURVES FOR EXTERNAL VALIDATION

  # Convert outcome to numeric
  test_outcome_eicu <- as.numeric(as.character(eicu_test$sad_binary))
  
  # Create ROC objects
  roc_list_eicu <- list(
    LR = roc(test_outcome_eicu, lr_pred_prob_eicu, quiet = TRUE),
    SVM = roc(test_outcome_eicu, svm_pred_prob_eicu, quiet = TRUE),
    XGBoost = roc(test_outcome_eicu, xgb_pred_prob_eicu, quiet = TRUE),
    RF = roc(test_outcome_eicu, rf_pred_prob_eicu, quiet = TRUE),
    DT = roc(test_outcome_eicu, dt_pred_prob_eicu, quiet = TRUE),
    NB = roc(test_outcome_eicu, nb_pred_prob_eicu, quiet = TRUE),
    KNN = roc(test_outcome_eicu, knn_pred_prob_eicu, quiet = TRUE)
  )
  
  # Extract AUC values
  auc_values_eicu <- sapply(roc_list_eicu, function(x) auc(x))
  
  cat("External Validation AUCs (MIMIC → eICU) \n")
  for (model in names(auc_values_eicu)) {
    cat(sprintf("  %-10s: %.4f\n", model, auc_values_eicu[model]))
  }
  cat("\n")
  
  # FIGURE 4B
  # Prepare for ggplot
  roc_data_eicu <- data.frame()
  for (model_name in names(roc_list_eicu)) {
    roc_obj <- roc_list_eicu[[model_name]]
    
    temp_df <- data.frame(
      Model = rep(model_name, length(roc_obj$specificities)),
      Specificity = roc_obj$specificities,
      Sensitivity = roc_obj$sensitivities,
      AUC = rep(auc(roc_obj), length(roc_obj$specificities))
    )
    roc_data_eicu <- rbind(roc_data_eicu, temp_df)
  }
  
  # Create custom legend labels with AUC values
  auc_labels_eicu <- sapply(names(roc_list_eicu), function(model) {
    sprintf("%s (AUC=%.3f)", model, auc(roc_list_eicu[[model]]))
  })
  
  # Create ROC plot
  fig4B_base <- ggplot(roc_data_eicu, aes(x = 1 - Specificity, y = Sensitivity, 
                                          color = Model, group = Model)) +
    geom_line(linewidth = 1.2, alpha = 0.8) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", 
                color = "black", linewidth = 0.8) +
    scale_color_manual(
      values = model_colors,
      labels = auc_labels_eicu,
      name = ""
    ) +
    labs(
      title = " ",
      x = "False Positive Rate",
      y = "True Positive Rate"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
      legend.position = c(0.85, 0.25),
      legend.background = element_rect(fill = "white", color = "white", linewidth = 0.5),
      legend.title = element_blank(),
      legend.text = element_text(size = 9),
      legend.key.size = unit(0.8, "cm"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
      axis.line = element_blank()
    ) +
    coord_equal()
  
  # Add "B" label 
  fig4B <- ggdraw() +
    draw_plot(fig4B_base) +
    draw_label("B", x = 0.02, y = 0.98, fontface = "bold", size = 16, hjust = 0, vjust = 1)
  
  print(fig4B)
  
  # Save figure
  ggsave("Figure4B_ROC_External_Validation.png", fig4B, 
         width = 8, height = 8, dpi = 300)
  
  cat(" Saved: Figure4B_ROC_External_Validation.png\n\n")
  
  # SAVE COMPARISON TABLE

  comparison_df <- data.frame(
    Model = names(auc_values),
    AUC_Internal_MIMIC = round(auc_values, 4),
    AUC_External_eICU = round(auc_values_eicu, 4),
    Difference = round(auc_values - auc_values_eicu, 4)
  )
  comparison_df <- comparison_df[order(-comparison_df$AUC_External_eICU), ]
  
  cat(" AUC COMPARISON: Internal vs External \n")
  print(comparison_df, row.names = FALSE)
  
  write.csv(comparison_df, "AUC_Comparison_Internal_vs_External.csv", row.names = FALSE)
  cat("\n✓ Saved: AUC_Comparison_Internal_vs_External.csv\n")
  
  # Save feature mapping for reference
  mapping_df <- data.frame(
    MIMIC_Feature = selected_features,
    eICU_Feature = sapply(selected_features, function(f) {
      mapped <- map_feature_name(f)
      if (is.null(mapped)) return("NOT AVAILABLE")
      return(mapped)
    }),
    Status = sapply(selected_features, function(f) {
      if (f %in% available_mimic_features) "Available" else "Imputed"
    })
  )
  
  write.csv(mapping_df, "Feature_Mapping_MIMIC_to_eICU.csv", row.names = FALSE)
  cat(" Saved: Feature_Mapping_MIMIC_to_eICU.csv\n")
  
} else {
  cat(" eICU data not found. Skipping external validation.\n")
  cat(" Please ensure 'sepsis_cohort_eicu_complete.RData' exists.\n\n")
}

#########################################################################################################
# FIGURE 4C: CALIBRATION CURVES FOR TOP 3 MODELS 
#########################################################################################################

# Prepare calibration data for top 3 models (XGBoost, RF, SVM)
calibration_data <- data.frame(
  Actual = as.numeric(as.character(test$sad_binary)),
  XGBoost = xgb_pred_prob,
  RF = rf_pred_prob,
  SVM = svm_pred_prob
)

# Convert to long format
cal_data <- reshape2::melt(
  calibration_data,
  id.vars = "Actual",
  variable.name = "Model",
  value.name = "Predicted"
)

# Create decile bins within each model
cal_data <- cal_data %>%
  group_by(Model) %>%
  mutate(bin = ntile(Predicted, 10)) %>%
  ungroup()

# Summarize per bin
calibration_plot <- cal_data %>%
  group_by(Model, bin) %>%
  summarise(
    Predicted_Probability = mean(Predicted, na.rm = TRUE),
    Actual_Probability = mean(Actual, na.rm = TRUE),
    n = n(),
    .groups = "drop"
  )


# Create the base calibration plot WITHOUT label
fig4C_base <- ggplot(
  calibration_plot,
  aes(x = Predicted_Probability, y = Actual_Probability, color = Model)
) +
  geom_line(linewidth = 1.2) +
  geom_abline(
    intercept = 0,
    slope = 1,
    linewidth = 0.8,
    color = "gray50",
    linetype = "dashed"
  ) +
  # Manually add "Ideal" to legend by creating a dummy geom
  geom_line(
    data = data.frame(
      Predicted_Probability = c(0, 1),
      Actual_Probability = c(0, 1),
      Model = "Ideal"
    ),
    linetype = "dashed",
    linewidth = 0.8
  ) +
  scale_color_manual(
    values = c(
      "XGBoost" = "#E41A1C",
      "RF" = "#984EA3",
      "SVM" = "#377EB8",
      "Ideal" = "#020c0d"
    )
  ) +
  coord_equal() +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  labs(
    x = "Predicted Probability",
    y = "Actual Probability"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    legend.title = element_blank(),
    legend.position = c(0.85, 0.25), 
    legend.background = element_rect(fill = "white", color = "white", linewidth = 0.5),
    legend.text = element_text(size = 10),
    axis.line = element_blank()
  )

# Add "C" label 
fig4C <- ggdraw() +
  draw_plot(fig4C_base) +
  draw_label("C", x = 0.02, y = 0.98, fontface = "bold", size = 16, hjust = 0, vjust = 1)

print(fig4C)

# Save figure
ggsave("Figure4C_Calibration_Curves.png", fig4C, 
       width = 8, height = 8, dpi = 300)



#########################################################################################################
# FIGURE 4D: DECISION CURVE ANALYSIS FOR TOP 3 MODELS
#########################################################################################################

# Function to compute net benefit
compute_dca <- function(pred_prob, true_outcome, thresholds) {
  keep <- is.finite(pred_prob) & is.finite(true_outcome)
  p <- pred_prob[keep]
  y <- true_outcome[keep]
  
  y <- as.integer(y == 1)
  n <- length(y)
  
  sapply(thresholds, function(th) {
    if(th == 0 || th == 1) return(NA)  # Avoid division by zero
    pred_positive <- p >= th
    tp <- sum(pred_positive & y == 1)
    fp <- sum(pred_positive & y == 0)
    (tp / n) - (fp / n) * (th / (1 - th))
  })
}

# Prepare DCA data
dca_data <- data.frame(
  outcome = as.numeric(as.character(test$sad_binary)),  
  XGBoost = xgb_pred_prob,
  RF = rf_pred_prob,
  SVM = svm_pred_prob
)

# Calculate prevalence
prevalence <- mean(dca_data$outcome)
cat(sprintf("Prevalence: %.3f\n", prevalence))

# Define threshold range (1% to 99%)
thresholds <- seq(0.01, 0.99, by = 0.01)

# Compute net benefit for all strategies
dca_results <- data.frame(
  threshold = thresholds,
  XGBoost = compute_dca(dca_data$XGBoost, dca_data$outcome, thresholds),
  RF = compute_dca(dca_data$RF, dca_data$outcome, thresholds),
  SVM = compute_dca(dca_data$SVM, dca_data$outcome, thresholds),
  TreatAll = prevalence - (1 - prevalence) * thresholds / (1 - thresholds),
  TreatNone = 0
)

# Convert to long format
dca_long <- reshape2::melt(dca_results, id.vars = "threshold",
                           variable.name = "Model", value.name = "NetBenefit")

# Remove NA values
dca_long <- dca_long[!is.na(dca_long$NetBenefit), ]

# Define colors
dca_colors <- c(
  "XGBoost" = "#E41A1C",
  "RF" = "#984EA3", 
  "SVM" = "#377EB8",
  "TreatAll" = "#4DAF4A",
  "TreatNone" = "#d3be14"
)

# Create the base DCA plot WITHOUT label
fig4D_base <- ggplot(dca_long, aes(x = threshold, y = NetBenefit, 
                                   color = Model, linetype = Model)) +
  geom_line(linewidth = 1.2) +
  geom_hline(yintercept = 0, linetype = "dotted", color = "white") +
  scale_color_manual(
    values = dca_colors,
    name = ""
  ) +
  scale_linetype_manual(
    values = c(
      "XGBoost" = "solid",
      "RF" = "solid",
      "SVM" = "solid", 
      "TreatAll" = "solid",
      "TreatNone" = "solid"
    ),
    name = ""
  ) +
  scale_x_continuous(
    labels = scales::percent_format(),
    breaks = seq(0, 1, 0.25),
    limits = c(0, 1),
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    limits = c(-0.05, 0.35),
    breaks = seq(0, 0.35, 0.05),
    expand = c(0, 0)
  ) +
  labs(
    x = "Threshold Probability",
    y = "Net Benefit"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    panel.grid.major = element_line(color = "white", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    legend.position = c(0.85, 0.75),
    legend.background = element_rect(fill = "white", color = "white", linewidth = 0.5),
    legend.title = element_blank(),
    legend.text = element_text(size = 10),
    axis.line = element_blank()
  )

# Add "D" label 
fig4D <- ggdraw() +
  draw_plot(fig4D_base) +
  draw_label("D", x = 0.02, y = 0.98, fontface = "bold", size = 16, hjust = 0, vjust = 1)

print(fig4D)

# Save
ggsave("Figure4D_DCA.png", fig4D, 
       width = 8, height = 6, dpi = 300)


#########################################################################################################
# FIGURE 5: SHAP ANALYSIS 
#########################################################################################################

#########################################################################################################
#  LOAD  XGBOOST MODEL/ RECREATE TRAINING DATA
#########################################################################################################

# Load saved XGBoost model
xgb_model <- readRDS("zhang_xgboost_model.rds")

# Get model features
model_features <- xgb_model$feature_names
cat(sprintf("Model expects %d features\n", length(model_features)))

# RECREATE train/test split with LASSO-selected features
# Load the original complete data
load("sepsis_cohort_complete.RData")
sepsis_cohort <- complete_cohort
setDT(sepsis_cohort)

# Remove time-related variables
sepsis_cohort <- subset(sepsis_cohort, select = -c(icu_admit_time, 
                                                   icu_discharge_time, 
                                                   death_time,
                                                   cam_positive, 
                                                   sad_group))

# Get predictor columns (everything except outcome-related)
exclude_vars <- c("sad_binary", "icd_delirium", "cam_status")
predictor_cols <- setdiff(names(sepsis_cohort), exclude_vars)

# Keep only numeric predictors
predictor_data <- sepsis_cohort[, ..predictor_cols]
numeric_cols <- names(predictor_data)[sapply(predictor_data, function(x) {
  is.numeric(x) || is.integer(x)
})]
predictor_data <- predictor_data[, ..numeric_cols]

# Remove zero-variance columns
predictor_data <- predictor_data[, which(sapply(predictor_data, function(x) {
  !all(is.na(x)) && (sd(x, na.rm = TRUE) > 0)
})), with = FALSE]

# Add outcome back
data_for_modeling <- predictor_data
data_for_modeling$sad_binary <- sepsis_cohort$sad_binary
data_for_modeling <- data_for_modeling[!is.na(sad_binary)]

# Recreate train/test split with SAME seed
set.seed(123)
index <- createDataPartition(data_for_modeling$sad_binary, p = 0.7, list = FALSE)
train_original <- data_for_modeling[index, ]
test_original <- data_for_modeling[-index, ]

# Filter to model features only
train_shap <- train_original[, c(model_features, "sad_binary"), with = FALSE]

cat(sprintf("✓ Recreated training data: %d rows × %d features\n", 
            nrow(train_shap), length(model_features)))

#########################################################################################################
# DATA FOR SHAP
#########################################################################################################

# Prepare feature matrix (numeric only, no outcome)
train_x <- as.matrix(train_shap[, ..model_features])
storage.mode(train_x) <- "numeric"

# View dimensions
cat(sprintf(" Feature matrix: %d rows × %d features\n", 
            nrow(train_x), ncol(train_x)))

#########################################################################################################
# GET SHAP VALUES
#########################################################################################################

# Quick access 
dtrain <- xgb.DMatrix(data = train_x)

# Compute SHAP contributions
shap_contrib <- predict(xgb_model, dtrain, predcontrib = TRUE)
shap_contrib <- as.data.table(shap_contrib)
shap_contrib[, ID := .I]

# Get feature columns (exclude BIAS)
feature_cols <- setdiff(colnames(shap_contrib), c("BIAS", "ID"))

cat(sprintf(" SHAP computed for %d features\n", length(feature_cols)))

#########################################################################################################
# FIGURE 5A: SHAP FEATURE IMPORTANCE (TOP 15 FEATURES)
#########################################################################################################

# Calculate mean absolute SHAP values
shap_importance <- data.table(
  feature = feature_cols,
  importance = colMeans(abs(as.matrix(shap_contrib[, ..feature_cols])))
)

shap_importance <- shap_importance[order(-importance)]

# Select top 15
top_n <- 15
top_15_features <- shap_importance$feature[1:top_n]
shap_importance_top15 <- shap_importance[1:top_n, ]

cat(sprintf("✓ Top %d features selected\n", top_n))
print(shap_importance_top15)

# Create bar plot
fig5a <- ggplot(shap_importance_top15,
                aes(x = reorder(feature, importance), y = importance)) +
  geom_bar(stat = "identity", fill = "#f0b70f") +
  geom_text(aes(label = sprintf("%.4f", importance)), 
            hjust = -0.1, 
            size = 3,
            color = "black") +
  coord_flip() +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(
    x = "Feature", 
    y = "Mean (|SHAP Value|)"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.y = element_text(size = 11),
    panel.grid.major.y = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
  )

print(fig5a)

ggsave("Figure5A_SHAP_Importance_Top15.png", fig5a, 
       width = 8, height = 6, dpi = 300)

#########################################################################################################
# FIGURE 5B: SHAP SUMMARY PLOT 
#########################################################################################################

# Use the SAME top 15 features from Panel A
feature_cols_top15 <- as.character(top_15_features)

# Prepare SHAP values (long format)
shap_long <- melt(
  shap_contrib,
  id.vars = "ID",
  measure.vars = feature_cols_top15,
  variable.name = "variable",
  value.name = "value"
)

# Get raw feature values
train_x_top15 <- train_x[, feature_cols_top15, drop = FALSE]

train_long <- as.data.table(train_x_top15)
train_long[, ID := .I]
train_long <- melt(
  train_long,
  id.vars = "ID",
  variable.name = "variable",
  value.name = "rfvalue"
)

setDT(train_long)
train_long$variable <- as.character(train_long$variable)
shap_long$variable <- as.character(shap_long$variable)

# Merge SHAP and feature values
shap_plot_data <- merge(shap_long, train_long, by = c("ID", "variable"), all.x = TRUE)
setDT(shap_plot_data)

# Remove NAs
shap_plot_data <- shap_plot_data[complete.cases(shap_plot_data), ]

# Add SHAP importance as mean_value column (this is what gets displayed)
shap_importance_lookup <- setNames(shap_importance_top15$importance, 
                                   shap_importance_top15$feature)

shap_plot_data[, mean_value := shap_importance_lookup[as.character(variable)]]

# Standardize feature values for color gradient
shap_plot_data[, stdfvalue := {
  feat_mean <- mean(rfvalue, na.rm = TRUE)
  feat_sd <- sd(rfvalue, na.rm = TRUE)
  if(feat_sd > 0) {
    (rfvalue - feat_mean) / feat_sd
  } else {
    0
  }
}, by = variable]

# Set factor levels in SAME ORDER as Panel A
shap_plot_data$variable <- factor(
  shap_plot_data$variable,
  levels = top_15_features
)

cat(sprintf("✓ %d data points for plotting\n", nrow(shap_plot_data)))
# Use the built-in shap.plot.summary function
fig5b <- shap.plot.summary(
  data_long = shap_plot_data,
  scientific = FALSE
)

# Customize colors and legend
fig5b_final <- fig5b +
  scale_color_gradient(
    low = "#7c0aba",
    high = "#f9c306",
    limits = c(-2, 2),
    oob = scales::squish,
    breaks = c(-2, 0, 2),
    labels = c("Low", "", "High")
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 11, face = "bold"),
    legend.text = element_text(size = 10),
    legend.key.height = unit(3, "cm"),
    legend.key.width = unit(0.4, "cm"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    panel.grid.major.x = element_line(color = "gray90", linewidth = 0.3),
    panel.grid.minor = element_blank()
  ) +
  guides(
    color = guide_colorbar(
      title = "Feature value",
      title.position = "top",
      title.hjust = 0.5,
      barheight = unit(3, "cm"),
      barwidth = unit(0.4, "cm"),
      ticks = FALSE 
    )
  )

# Override the x-axis scale to match the paper
fig5b_final <- fig5b_final +
  ylab("SHAP value")

# View
print(fig5b_final)

ggsave("Figure5B_SHAP_Summary_Top15.png", fig5b_final, 
       width = 10, height = 7, dpi = 300)


#########################################################################################################
# Combine 5A/B 
#########################################################################################################

# Combine the two plots
combined_fig5 <- fig5a + fig5b_final + 
  plot_layout(ncol = 2, widths = c(1, 1.2)) +
  plot_annotation(
    title = " ",
    tag_levels = "A"
  )

# Display in RStudio
print(combined_fig5)

# Save to file
ggsave(
  filename = "Figure5_Combined_Top15_3.png",
  plot = combined_fig5,
  width = 18, height = 7, dpi = 300  # ✅ Adjusted for 15 features
)

#########################################################################################################
# CREATE FXN FOR SHAP-BASED PARTIAL DEPENDENCE PLOTS
#########################################################################################################
make_shap_pdp <- function(feature_name, feature_data, shap_data) {
  
  plot_data <- data.frame(
    feature_value = feature_data[[feature_name]],
    shap_value = shap_data[[feature_name]]
  )
  
  plot_data <- plot_data[complete.cases(plot_data), ]
  
  # Detect binary feature
  unique_vals <- sort(unique(plot_data$feature_value))
  is_binary <- length(unique_vals) == 2 && all(unique_vals %in% c(0, 1))
  
  p <- ggplot(plot_data, aes(x = feature_value, y = shap_value)) +
    geom_hline(
      yintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.5
    )
  
  if (is_binary) {
    # Binary feature: jitter + mean line
    p <- p +
      geom_jitter(
        width = 0.05, height = 0, alpha = 0.3, size = 0.6, color = "#5624ac"
      ) +
      stat_summary(fun = mean, geom = "line", aes(group = 1),
                   color = "red", linewidth = 1.2
      )
  } else {
    # Continuous feature: scatter + LOESS
    p <- p +
      geom_point(alpha = 0.3, size = 0.5, color = "#5624ac") +
      geom_smooth(
        method = "loess",
        span = 0.75,
        se = FALSE,
        color = "red",
        linewidth = 1.2
      )
  }
  
  p +
    labs(
      title = '',
      x = feature_name,
      y = "SHAP value"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 11),
      axis.title.x = element_text(size = 9),
      axis.title.y = element_text(size = 9),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)  # Add border
    )
}

#########################################################################################################
# GET SHAP PDP FOR TOP 15 FEATURES 
#########################################################################################################

# Extract SHAP values and feature values for top 15
shap_values_df <- as.data.frame(shap_contrib[, ..feature_cols_top15])
feature_values_df <- as.data.frame(train_x_top15)

cat("\nGenerating SHAP-based PDPs\n")

# Use feature_cols_top15 instead of all_top15_features
pdp_list <- lapply(
  feature_cols_top15,  # ✅ FIXED: was all_top15_features
  make_shap_pdp, 
  feature_data = feature_values_df,
  shap_data = shap_values_df
)


#########################################################################################################
# MAKE A 5 X 3 LAYOUT LIKE THE PAPER 
#########################################################################################################

fig6 <- plot_grid(plotlist = pdp_list, ncol = 5, nrow = 3, align = "hv")

fig6_final <- ggdraw() +
  draw_label(
    " ",
    fontface = "bold",
    x = 0.5, y = 0.98, size = 14, hjust = 0.5
  ) +
  draw_label(
    "",
    x = 0.5, y = 0.01, size = 8, hjust = 0.5, color = "gray30"
  ) +
  draw_plot(fig6, y = 0.04, height = 0.92)

print(fig6_final)

# Save
ggsave(filename = "Figure6_SHAP_PDPs_Top15.png", plot = fig6_final, width = 16, height = 9, dpi = 300)

