#########################################################################################################
# ZHANG PAPER REPRODUCTION
# NOVEL EXTENSION
#########################################################################################################

# Clean work space
rm(list = ls())

# Load packages
library(data.table)
library(caret)
library(glmnet)
library(e1071)
library(xgboost)
library(randomForest)
library(rpart)
library(kknn)
library(naivebayes)
library(pROC)
library(ggplot2)
library(reshape2)
library(cowplot)
library(patchwork)
library(dplyr)
library(tidyr)

# Load finalized SAD cohort (post MICE) 
load("sepsis_cohort_complete.RData")  
# Rename
sepsis_cohort <- complete_cohort
# Make table
setDT(sepsis_cohort)

# Remove time related variables (related to figure 2)
sepsis_cohort <- subset(sepsis_cohort, select = -c(icu_admit_time, 
                                                   icu_discharge_time, death_time,
                                                   cam_positive, sad_group))

# Load LASSO-selected features
lasso_features <- read.csv("LASSO_Selected_Features.csv", stringsAsFactors = FALSE)
selected_features <- lasso_features$Feature

cat(sprintf("Loaded %d LASSO-selected features\n", length(selected_features)))

model_data <- sepsis_cohort[, c(selected_features, "sad_binary"), with = FALSE]
model_data <- as.data.frame(model_data)

cat(sprintf("Model data: %d rows × %d columns\n", 
            nrow(model_data), ncol(model_data)))

# Outcome check
print(table(model_data$sad_binary))
#########################################################################################################
# TRAIN/TEST SPLIT
#########################################################################################################
# Seed was the same as the original analysis for the reproduction portion of the analysis

set.seed(123)  

train_index <- createDataPartition(model_data$sad_binary, p = 0.7, list = FALSE)
train <- model_data[train_index, ]
test <- model_data[-train_index, ]

cat(sprintf("Train: %d patients (%.1f%% SAD)\n",
  nrow(train), mean(train$sad_binary) * 100
))
cat(sprintf("Test:  %d patients (%.1f%% SAD)\n",
  nrow(test), mean(test$sad_binary) * 100
))

# Factor versions (for SVM, NB, RF)
train_factor <- train
test_factor  <- test
train_factor$sad_binary <- factor(train_factor$sad_binary, levels = c(0, 1))
test_factor$sad_binary <- factor(test_factor$sad_binary, levels = c(0, 1))

# Verification check
stopifnot(
  all(selected_features %in% names(sepsis_cohort)),
  "sad_binary" %in% names(sepsis_cohort)
)

cat("Feature and outcome alignment verified\n")

#########################################################################################################
# TRAIN ALL BASE MODELS (WITH CLASS WEIGHTS)
#########################################################################################################

# Ensure outcome numeric 0/1
train$sad_binary <- as.integer(train$sad_binary)
test$sad_binary  <- as.integer(test$sad_binary)
stopifnot(all(train$sad_binary %in% c(0,1)))
stopifnot(all(test$sad_binary %in% c(0,1)))

# Calculate class weights
prev1 <- mean(train$sad_binary)    
w0 <- 1 / (1 - prev1)               
w1 <- 1 / prev1   

# Normalize so average weight ~1 
w_sum <- w0 + w1
w0 <- 2 * w0 / w_sum
w1 <- 2 * w1 / w_sum

cat(sprintf("Class 1 prevalence: %.1f%% | w1=%.3f, w0=%.3f\n", prev1*100, w1, w0))

train_weights <- ifelse(train$sad_binary == 1, w1, w0)

#########################################################################################################
# MODEL 1: LOGISTIC REGRESSION (WITH WEIGHTS)
#########################################################################################################

lr_model <- glm(sad_binary ~ ., data = train, family = binomial, weights = train_weights)
lr_pred_prob <- predict(lr_model, newdata = test, type = "response")
saveRDS(lr_model, "lr_model_weighted.rds")

lr_roc <- roc(test$sad_binary, lr_pred_prob, levels = c(0,1), direction = "<", quiet = TRUE)
lr_auc <- auc(lr_roc)
cat(sprintf("LR (weighted): AUC = %.4f\n", lr_auc))

#########################################################################################################
# MODEL 2: SVM (WITH CLASS WEIGHTS)
#########################################################################################################

# Ensure factor outcome
train_factor$sad_binary <- factor(train_factor$sad_binary, levels = c(0, 1))
test_factor$sad_binary  <- factor(test_factor$sad_binary,  levels = c(0, 1))

svm_class_weights <- c("0" = w0, "1" = w1)

svm_model <- svm(
  sad_binary ~ ., 
  data = train_factor,
  probability = TRUE,
  class.weights = svm_class_weights
)

svm_pred <- predict(svm_model, newdata = test_factor, probability = TRUE)
svm_pred_prob <- attr(svm_pred, "probabilities")[, "1"]

saveRDS(svm_model, "svm_model_weighted.rds")

svm_auc <- auc(roc(test$sad_binary, svm_pred_prob, quiet = TRUE))
cat(sprintf("SVM (weighted): AUC = %.4f\n", svm_auc))


#########################################################################################################
# MODEL 3: XGBOOST (WITH SCALE_POS_WEIGHT)
#########################################################################################################

scale_pos_weight <- sum(train$sad_binary == 0) / sum(train$sad_binary == 1)
cat(sprintf("   scale_pos_weight: %.3f\n", scale_pos_weight))

train_x <- data.matrix(train[, selected_features])
test_x  <- data.matrix(test[, selected_features])

train_matrix <- xgb.DMatrix(data = train_x, label = train$sad_binary)
test_matrix  <- xgb.DMatrix(data = test_x,  label = test$sad_binary)

params <- list(
  objective = "binary:logistic",
  eval_metric = "logloss",
  max_depth = 3,
  eta = 0.1,
  gamma = 0.5,
  colsample_bytree = 1,
  min_child_weight = 1,
  subsample = 0.5,
  scale_pos_weight = scale_pos_weight
)

xgb_model <- xgb.train(
  params = params,
  data = train_matrix,
  nrounds = 125,
  watchlist = list(train = train_matrix, val = test_matrix),
  early_stopping_rounds = 10,
  print_every_n = 10,
  maximize = FALSE
)

saveRDS(xgb_model, "xgb_model_weighted.rds")

xgb_pred_prob <- predict(xgb_model, test_matrix)

xgb_auc <- auc(roc(test$sad_binary, xgb_pred_prob, quiet = TRUE))
cat(sprintf("XGBoost (weighted): AUC = %.4f\n", xgb_auc))


#########################################################################################################
# MODEL 4: RANDOM FOREST (WITH CLASS WEIGHTS)
#########################################################################################################

rf_class_weights <- c(w0, w1)
names(rf_class_weights) <- c("0", "1")

n_features <- ncol(train) - 1
mtry_value <- max(1, floor(sqrt(n_features)))

rf_model <- randomForest(
  sad_binary ~ .,
  data = train_factor,
  ntree = 500,
  mtry = mtry_value,
  classwt = rf_class_weights
)

rf_pred_prob <- predict(rf_model, newdata = test_factor, type = "prob")[, "1"]

saveRDS(rf_model, "rf_model_weighted.rds")

rf_auc <- auc(roc(test$sad_binary, rf_pred_prob, quiet = TRUE))
cat(sprintf("RF (weighted): AUC = %.4f\n", rf_auc))

#########################################################################################################
# MODEL 5: DECISION TREE (WITH LOSS MATRIX)
#########################################################################################################

loss_matrix <- matrix(c(
  0,  w1,   # FP cost
  w0, 0     # FN cost
), nrow = 2, byrow = TRUE)

dt_model <- rpart(
  sad_binary ~ .,
  data = train_factor,
  method = "class",
  parms = list(loss = loss_matrix),
  control = rpart.control(cp = 0.01)
)

dt_pred_prob <- predict(dt_model, newdata = test_factor, type = "prob")[, "1"]
saveRDS(dt_model, "dt_model_weighted.rds")

dt_splits <- if (!is.null(dt_model$splits)) nrow(dt_model$splits) else 0

if (dt_splits == 0) {
  cat("️ DT failed to split — excludine\n")
  dt_failed <- TRUE
  dt_auc <- 0.5
} else {
  dt_failed <- FALSE
  dt_auc <- auc(roc(test$sad_binary, dt_pred_prob, quiet = TRUE))
  cat(sprintf(" DT (weighted): AUC = %.4f\n", dt_auc))
}

# NOTE: DT failed, so excluding it from the rest of the models
#########################################################################################################
# MODEL 6: NAIVE BAYES (WITH CLASS PRIORS)
#########################################################################################################

# Define priors from training distribution (or inverse prevalence)
nb_priors <- c(
  "0" = w0 / (w0 + w1),
  "1" = w1 / (w0 + w1)
)

nb_model <- naive_bayes(
  sad_binary ~ .,
  data  = train_factor,
  prior = nb_priors
)

nb_pred_prob <- predict(nb_model, newdata = test_factor, type = "prob")[, "1"]

saveRDS(nb_model, "nb_model_weighted_priors.rds")

nb_auc <- auc(roc(test$sad_binary, nb_pred_prob, quiet = TRUE))

cat(sprintf(" NB (with priors): AUC = %.4f\n", nb_auc))

#########################################################################################################
# MODEL 7: KNN (DISTANCE-WEIGHTED, SCALED FEATURES)
#########################################################################################################

# Scale features (REQUIRED for KNN)
preproc <- preProcess(train_factor[, selected_features], method = c("center", "scale"))

train_knn <- train_factor
test_knn  <- test_factor

train_knn[, selected_features] <- predict(preproc, train_factor[, selected_features])
test_knn[, selected_features]  <- predict(preproc, test_factor[, selected_features])

# Distance-weighted KNN (closer neighbors matter more)
knn_model <- kknn(
  sad_binary ~ .,
  train  = train_knn,
  test   = test_knn,
  k      = 10,
  distance = 2,
  kernel = "triangular"   # distance-weighted voting
)

knn_pred_prob <- knn_model$prob[, "1"]

saveRDS(knn_model, "knn_model_distance_weighted.rds")

knn_auc <- auc(roc(test$sad_binary, knn_pred_prob, quiet = TRUE))

cat(sprintf(" KNN (distance-weighted): AUC = %.4f\n", knn_auc))

#########################################################################################################
# COMBINED ROC CURVES (WEIGHTED MODELS)
#########################################################################################################

# Recall only 6/7 models worked, DT is excluded
roc_list <- list(
  "Logistic Regression" = roc(test$sad_binary, lr_pred_prob, quiet = TRUE),
  "SVM"                 = roc(test$sad_binary, svm_pred_prob, quiet = TRUE),
  "XGBoost"             = roc(test$sad_binary, xgb_pred_prob, quiet = TRUE),
  "Random Forest"       = roc(test$sad_binary, rf_pred_prob,  quiet = TRUE),
  "Naive Bayes"         = roc(test$sad_binary, nb_pred_prob,  quiet = TRUE),
  "KNN"                 = roc(test$sad_binary, knn_pred_prob, quiet = TRUE)
)

auc_values <- sapply(roc_list, auc)

cat("\n Individual Model AUCs:\n")
for (i in seq_along(auc_values)) {
  cat(sprintf("   %s: %.4f\n", names(auc_values)[i], auc_values[i]))
}

# Convert to df with AUC replicated for each row
roc_df <- do.call(rbind, lapply(names(roc_list), function(m) {
  r <- roc_list[[m]]
  auc_val <- auc(r)  # Calculate once
  n_points <- length(r$specificities)  # Number of points in ROC curve
  
  data.frame(
    FPR = 1 - r$specificities,
    TPR = r$sensitivities,
    Model = m,                           # Just model name
    AUC = rep(auc_val, n_points)         # Replicate AUC for each point
  )
}))

# Create labels with AUC
roc_df$Model_Label <- sprintf("%s (AUC = %.3f)", roc_df$Model, roc_df$AUC)

# Define colors (only include the 6 models you're using)
model_colors <- c(
  "Logistic Regression" = "#4DAF4A",
  "SVM"                 = "#377EB8",
  "XGBoost"             = "#E41A1C",
  "Random Forest"       = "#984EA3",
  "Naive Bayes"         = "#FFFF33",
  "KNN"                 = "#A65628"
)

# Plot
roc_plot <- ggplot(roc_df, aes(x = FPR, y = TPR, color = Model)) +
  geom_line(linewidth = 1.2) +
  geom_abline(
    slope = 1,
    intercept = 0,
    linetype = "dashed",
    color = "black"
  ) +
  scale_color_manual(
    values = model_colors,
    name = "Model"
  ) +
  coord_equal() +
  labs(
    title = "ROC Curves for Class-Weighted Models (SAD)",
    x = "False Positive Rate",
    y = "True Positive Rate"
  ) +
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5),   
    legend.position = c(0.75, 0.25),
    legend.background = element_rect(fill = "white", color = "white", linewidth = 0.5),
    legend.title = element_blank(),  
    legend.text = element_text(size = 9),
    legend.key.width = unit(1.2, "lines"),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
  )

print(roc_plot)

ggsave("Figure_Novel_weighted_models_ROC.png", roc_plot, width = 8, height = 6, dpi = 300)

#########################################################################################################
# RISK STRATIFICATION 
#########################################################################################################

# Using best model: XGBoost

risk_df <- data.frame(
  true_outcome = test$sad_binary,
  pred_prob    = xgb_pred_prob
)

# Define clinically interpretable bins
risk_df$risk_group <- cut(
  risk_df$pred_prob,
  breaks = c(0, 0.20, 0.40, 0.60, 1.00),
  labels = c("Low risk", "Intermediate risk", "High risk", "Very high risk"),
  include.lowest = TRUE
)

# Summarize risk by group
risk_summary <- risk_df %>%
  group_by(risk_group) %>%
  summarise(
    N = n(),
    Observed_SAD_rate = mean(true_outcome),
    .groups = "drop"
  )

# View
print(risk_summary)

# View stratification in a barplot
Fig1 <- ggplot(risk_summary, aes(x = risk_group, y = Observed_SAD_rate, fill = risk_group)) +
  geom_col(width = 0.7, color = "black") +
  geom_text(aes(label = sprintf("%.1f%%", Observed_SAD_rate * 100)),
            vjust = -0.5, size = 5) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                     expand = expansion(mult = c(0, 0.1))) +
  scale_fill_brewer(palette = "YlOrRd") +
  labs(
    title = "Observed Delirium Risk Across Model-Predicted Risk Groups",
    x = "Predicted Risk Group (XGBoost)",
    y = "Observed SAD Incidence"
  ) +
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "none",
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
  )

print(Fig1)
# Save
ggsave("Figure_Novel1_barplot.png", Fig1, width = 8, height = 6, dpi = 300)

#########################################################################################################
# CALIBRATION ANALYSIS 
#########################################################################################################

calib_df <- data.frame(
  pred_prob = xgb_pred_prob,
  true_outcome = test$sad_binary
)

# Bin 
calib_df$bin <- cut(
  calib_df$pred_prob,
  breaks = quantile(calib_df$pred_prob, probs = seq(0, 1, 0.1)),
  include.lowest = TRUE,
  labels = FALSE
)

calib_summary <- calib_df %>%
  group_by(bin) %>%
  summarise(
    Mean_predicted_risk = mean(pred_prob),
    Observed_risk = mean(true_outcome),
    N = n(),
    .groups = "drop"
  )

# View
print(calib_summary)

# Calibration Plot
Fig2 <- ggplot(calib_summary, aes(x = Mean_predicted_risk, y = Observed_risk)) +
  geom_point(size = 3, color = "#377EB8") +
  geom_line(color = "#377EB8", linewidth = 1) +
  geom_abline(
    intercept = 0,
    slope = 1,
    linetype = "dashed",
    color = "gray50"
  ) +
  coord_equal() +
  scale_x_continuous(limits = c(0, 1), labels = scales::percent_format()) +
  scale_y_continuous(limits = c(0, 1), labels = scales::percent_format()) +
  labs(
    title = "Calibration Plot for Weighted XGBoost Model",
    x = "Mean Predicted Probability",
    y = "Observed Probability"
  ) +
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
  )

print(Fig2)
# Save
ggsave("Figure_Novel2_callibration.png", Fig2, width = 8, height = 6, dpi = 300)

# Brier score
brier_score <- mean((xgb_pred_prob - test$sad_binary)^2)
cat(sprintf("Brier score: %.4f\n", brier_score))

#########################################################################################################
# CONFORMAL PREDICTION (SPLIT-CONFORMAL)
#########################################################################################################

set.seed(222)

# True labels and predicted probabilities from XGBoost
y_test <- test$sad_binary
p_test <- xgb_pred_prob   # from weighted XGBoost

n_test <- length(y_test)

# Split test set into calibration and evaluation (50/50)
cal_idx  <- sample(seq_len(n_test), size = floor(0.5 * n_test))
eval_idx <- setdiff(seq_len(n_test), cal_idx)

y_cal  <- y_test[cal_idx]
p_cal  <- p_test[cal_idx]

y_eval <- y_test[eval_idx]
p_eval <- p_test[eval_idx]

# Nonconformity scores on calibration set
s_cal <- ifelse(y_cal == 1, 1 - p_cal, p_cal)

#coverage level and thershold
alpha <- 0.10   # 90% marginal coverage

# Conformal threshold
qhat <- quantile(s_cal,probs = ceiling((length(s_cal) + 1) * (1 - alpha)) / length(s_cal),
  type = 1)

cat(sprintf("Conformal threshold q̂ (α = %.2f): %.4f\n", alpha, qhat))

# Prediction sets on evaluation data
# Include class 1 if its nonconformity <= qhat
include_1 <- (1 - p_eval) <= qhat

# Include class 0 if its nonconformity <= qhat
include_0 <- p_eval <= qhat

# Prediction set size
set_size <- include_0 + include_1

# Evalutate conformal properties
# coverage
covered <- ifelse(y_eval == 1, include_1, include_0)
coverage <- mean(covered)

cat(sprintf("Observed coverage: %.3f\n", coverage))
# confidence
prop_singleton <- mean(set_size == 1)
prop_ambiguous <- mean(set_size == 2)

cat(sprintf("Confident (singleton) predictions: %.1f%%\n", 100 * prop_singleton))
cat(sprintf("Uncertain (abstain) cases: %.1f%%\n", 100 * prop_ambiguous))

# Only perofrom on confident cases
singleton_idx <- which(set_size == 1)

pred_singleton <- ifelse(include_1[singleton_idx], 1, 0)
y_singleton    <- y_eval[singleton_idx]

auc_singleton <- auc(roc(y_singleton, p_eval[singleton_idx], quiet = TRUE))

cat(sprintf(
  "AUC on confident predictions only: %.4f (n = %d)\n",
  auc_singleton,
  length(singleton_idx)
))

# Create risk stratification abstention table
conformal_summary <- data.frame(
  RiskGroup = c("Confident Non-SAD", "Confident SAD", "Uncertain"),
  Count = c(
    sum(set_size == 1 & include_0),
    sum(set_size == 1 & include_1),
    sum(set_size == 2)
  )
)

conformal_summary$Percent <- 100 * conformal_summary$Count / sum(conformal_summary$Count)
print(conformal_summary)

# Plot conformal summary 
conf_sum <- ggplot(conformal_summary,
                   aes(x = RiskGroup, y = Percent, fill = RiskGroup)) +
  geom_bar(stat = "identity", width = 0.7) +
  geom_text(aes(label = sprintf("%.1f%%", Percent)),
            vjust = -0.5, size = 5) +
  scale_fill_manual(values = c(
    "Confident Non-SAD" = "#4DAF4A",
    "Confident SAD"     = "#E41A1C",
    "Uncertain"         = "#999999"
  )) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  labs(
    title = "Conformal Risk Stratification for SAD Prediction",
    subtitle = "90% marginal coverage (α = 0.10)",
    x = NULL,
    y = "Percentage of Patients"
  ) +
  theme_bw(base_size = 14) +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
  )

print(conf_sum)

# Save
ggsave("Figure_Novel3_conformal_suummary.png", conf_sum, width = 8, height = 6, dpi = 300)

# Coverage vs alpha plot
alphas <- seq(0.05, 0.25, by = 0.05)

coverage_df <- bind_rows(lapply(alphas, function(a) {
  q <- quantile(
    s_cal,
    probs = ceiling((length(s_cal) + 1) * (1 - a)) / length(s_cal),
    type = 1
  )
  inc1 <- (1 - p_eval) <= q
  inc0 <- p_eval <= q
  covered <- ifelse(y_eval == 1, inc1, inc0)
  
  data.frame(
    alpha = a,
    coverage = mean(covered)
  )
}))

cov <- ggplot(coverage_df, aes(x = alpha, y = coverage)) +
  geom_line(linewidth = 1.2, color = "#377EB8") +
  geom_point(size = 3, color = "#377EB8") +
  geom_abline(intercept = 1, slope = -1,
              linetype = "dashed", color = "gray50") +
  scale_x_continuous(labels = scales::percent_format()) +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(
    title = "Empirical Coverage vs Nominal Level",
    x = "Significance Level (α)",
    y = "Observed Coverage"
  ) +
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
  )

print(cov)
# Save
ggsave("Figure_Novel4_coverage_alpha.png", cov, width = 8, height = 6, dpi = 300)

#########################################################################################################
# NOVEL METHOD 4: DISTRIBUTION SHIFT ANALYSIS
#########################################################################################################
# Helper fxns

#Find a column w/ multiple possible names
pick_col <- function(df, candidates) {
  hit <- intersect(candidates, names(df))
  if (length(hit) == 0) return(NULL)
  hit[1]
}

# Calibration model (logistic recalibration): y ~ logit(p)
calibration_stats <- function(y, p) {
  eps <- 1e-6
  p2 <- pmin(pmax(p, eps), 1 - eps)
  lp <- qlogis(p2)
  fit <- glm(y ~ lp, family = binomial())
  co <- coef(fit)
  list(
    intercept = unname(co[1]),
    slope = unname(co[2])
  )
}

# Brier score
brier_score_fn <- function(y, p) mean((p - y)^2)

# Expected Calibration Error (ECE)
ece_score <- function(y, p, n_bins = 10) {
  eps <- 1e-6
  p2 <- pmin(pmax(p, eps), 1 - eps)
  
  bins <- cut(
    p2,
    breaks = seq(0, 1, length.out = n_bins + 1),
    include.lowest = TRUE
  )
  
  df <- data.frame(y = y, p = p2, bin = bins)
  
  agg <- aggregate(
    cbind(y, p) ~ bin,
    df,
    function(x) c(mean = mean(x), n = length(x))
  )
  
  # Defensive unpacking
  ybar <- sapply(agg$y, function(z) ifelse(is.null(z["mean"]), NA, z["mean"]))
  pbar <- sapply(agg$p, function(z) ifelse(is.null(z["mean"]), NA, z["mean"]))
  n    <- sapply(agg$y, function(z) ifelse(is.null(z["n"]),    0,  z["n"]))
  
  # Drop empty bins
  keep <- n > 0
  if (sum(keep, na.rm = TRUE) == 0) return(NA_real_)
  
  sum((n[keep] / sum(n[keep])) * abs(ybar[keep] - pbar[keep]))
}

# AUROC
auroc <- function(y, p) as.numeric(auc(roc(y, p, quiet = TRUE)))

# Nonconformity for binary 
nonconformity <- function(y, p) ifelse(y == 1, 1 - p, p)

# Quantile used in split conformal
conformal_qhat <- function(s_cal, alpha) {
  n <- length(s_cal)
  idx <- ceiling((n + 1) * (1 - alpha)) / n
  as.numeric(quantile(s_cal, probs = idx, type = 1))
}

# Evaluate conformal prediction sets for a given qhat
eval_conformal <- function(y, p, qhat) {
  include_1 <- (1 - p) <= qhat
  include_0 <- p <= qhat
  set_size  <- include_0 + include_1
  covered   <- ifelse(y == 1, include_1, include_0)
  
  list(
    coverage = mean(covered),
    prop_singleton = mean(set_size == 1),
    prop_ambiguous = mean(set_size == 2),
    set_size = set_size
  )
}

# Prepare evaluation dataframe
y_test_int <- as.integer(test$sad_binary)
p_test_num <- as.numeric(xgb_pred_prob)

eval_df <- test
eval_df$y <- y_test_int
eval_df$p <- p_test_num

# Split test into calibration vs evaluation
set.seed(222)

n <- nrow(eval_df)
cal_idx_shift  <- sample(seq_len(n), size = floor(0.5 * n))
evl_idx_shift  <- setdiff(seq_len(n), cal_idx_shift)

df_cal <- eval_df[cal_idx_shift, ]
df_evl <- eval_df[evl_idx_shift, ]

# 90% target marginal coverage
alpha_shift <- 0.10   

s_cal_shift <- nonconformity(df_cal$y, df_cal$p)
qhat_global <- conformal_qhat(s_cal_shift, alpha_shift)

cat(sprintf("Global conformal qhat (alpha=%.2f): %.4f\n", alpha_shift, qhat_global))

# Define realistic within MIMIC distribution shifts

# Try to locate common columns (adjust candidate lists if needed)
col_SOFA <- pick_col(eval_df, c("sofa_total", "SOFA", "sofa", "sofa_score"))
col_MV   <- pick_col(eval_df, c("ventilated", "MV", "mech_vent", "ventilation"))
col_SED  <- pick_col(eval_df, c("any_sedative_used", "sedation", "sedative"))

cat("\n📊 Identified columns for distribution shifts:\n")
if (!is.null(col_SOFA)) cat(sprintf("   SOFA: %s\n", col_SOFA))
if (!is.null(col_MV)) cat(sprintf("   Mechanical ventilation: %s\n", col_MV))
if (!is.null(col_SED)) cat(sprintf("   Sedation: %s\n", col_SED))

# Define shifts
shifts <- list()

if (!is.null(col_SOFA)) {
  sofa_q75 <- quantile(df_evl[[col_SOFA]], 0.75, na.rm = TRUE)
  shifts[["High severity (SOFA Q4)"]] <- function(d) {
    if (!col_SOFA %in% names(d)) rep(FALSE, nrow(d))
    else d[[col_SOFA]] >= sofa_q75
  }
}

if (!is.null(col_MV)) {
  shifts[["Mechanically ventilated (MV=1)"]] <- function(d) {
    if (!col_MV %in% names(d)) rep(FALSE, nrow(d))
    else d[[col_MV]] == 1
  }
}

if (!is.null(col_SED)) {
  shifts[["Sedated (Sedation=1)"]] <- function(d) {
    if (!col_SED %in% names(d)) rep(FALSE, nrow(d))
    else d[[col_SED]] == 1
  }
}

# Always include baseline
shifts[["Baseline (all eval)"]] <- function(d) rep(TRUE, nrow(d))

cat("\n Defined distribution shifts:\n")
for (sname in names(shifts)) {
  cat(sprintf("   - %s\n", sname))
}

# Run robust evaluation
cat("\nEvaluating model performance across distribution shifts...\n")

results <- lapply(names(shifts), function(sname) {
  
  mask_evl <- shifts[[sname]](df_evl)
  d_e <- df_evl[mask_evl, ]
  n_e <- nrow(d_e)
  
  # Skip tiny groups
  if (n_e < 100) {
    return(data.frame(
      Shift = sname, n = n_e,
      AUROC = NA, Brier = NA, Cal_Intercept = NA, Cal_Slope = NA, ECE = NA,
      Cov_Global = NA, Sing_Global = NA, Amb_Global = NA,
      Cov_Group = NA,  Sing_Group = NA,  Amb_Group = NA
    ))
  }
  
  # Base metrics
  A <- auroc(d_e$y, d_e$p)
  B <- brier_score_fn(d_e$y, d_e$p)
  C <- calibration_stats(d_e$y, d_e$p)
  E <- ece_score(d_e$y, d_e$p, n_bins = 10)
  
  # Conformal using GLOBAL qhat
  conf_global <- eval_conformal(d_e$y, d_e$p, qhat_global)
  
  # Group-conditional conformal: compute qhat within same subgroup using calibration data
  mask_cal <- tryCatch(
    shifts[[sname]](df_cal),
    error = function(e) rep(FALSE, nrow(df_cal))
  )
  
  # Ensure logical vector of correct length
  if (!is.logical(mask_cal) || length(mask_cal) != nrow(df_cal)) {
    mask_cal <- rep(FALSE, nrow(df_cal))
  }
  
  d_c <- df_cal[mask_cal, , drop = FALSE]
  if (nrow(d_c) >= 100) {
    s_c <- nonconformity(d_c$y, d_c$p)
    qhat_group <- conformal_qhat(s_c, alpha_shift)
    conf_group <- eval_conformal(d_e$y, d_e$p, qhat_group)
  } else {
    conf_group <- list(coverage = NA, prop_singleton = NA, prop_ambiguous = NA)
  }
  
  data.frame(
    Shift = sname,
    n = n_e,
    AUROC = A,
    Brier = B,
    Cal_Intercept = C$intercept,
    Cal_Slope = C$slope,
    ECE = E,
    Cov_Global = conf_global$coverage,
    Sing_Global = conf_global$prop_singleton,
    Amb_Global  = conf_global$prop_ambiguous,
    Cov_Group = conf_group$coverage,
    Sing_Group = conf_group$prop_singleton,
    Amb_Group  = conf_group$prop_ambiguous
  )
})

shift_report <- do.call(rbind, results)

cat("\n Distribution Shift Analysis Results:\n")
print(shift_report)

# Save 
write.csv(shift_report, "distribution_shift_analysis.csv", row.names = FALSE)
cat("Saved: distribution_shift_analysis.csv\n")

# Create plot_df before using it
plot_df <- shift_report %>%
  select(Shift, AUROC, ECE, Cov_Global, Sing_Global) %>%
  pivot_longer(
    cols = c(AUROC, ECE, Cov_Global, Sing_Global),
    names_to = "Metric",
    values_to = "Value"
  ) %>%
  filter(!is.na(Value))

plot_df$Metric <- factor(
  plot_df$Metric,
  levels = c("AUROC", "ECE", "Cov_Global", "Sing_Global"),
  labels = c(
    "AUROC\n(Discrimination)", 
    "ECE\n(Calibration Error)", 
    "Coverage\n(90% target)",
    "Confident Predictions\n(%)"
  )
)

# Shift labels
plot_df$Shift <- dplyr::recode(
  plot_df$Shift,
  "Baseline (all eval)"                = "All Patients",
  "Sedated (Sedation=1)"               = "Sedated",
  "Mechanically ventilated (MV=1)"     = "Mechanically Ventilated",
  "High severity (SOFA Q4)"             = "High Severity"
)

# Legend order 
plot_df$Shift <- factor(
  plot_df$Shift,
  levels = c(
    "All Patients",
    "High Severity",
    "Mechanically Ventilated",
    "Sedated"
  )
)

# Add colors
shift_colors <- c(
  "Sedated"                 = "#F28E2B",
  "Mechanically Ventilated" = "#59A14F",
  "High Severity"           = "#E15759",
  "All Patients"            = "#4E79A7"
)

# Plot
shift_plot <- ggplot(plot_df, aes(x = Shift, y = Value, fill = Shift)) +
  geom_col(width = 0.6) +
  facet_wrap(
    ~ Metric,
    scales = "free_y",
    nrow = 1
  ) +
  coord_flip() +
  scale_fill_manual(
    values = shift_colors,
    breaks = c(
      "Sedated",
      "Mechanically Ventilated",
      "High Severity",
      "All Patients"
    ),
    name = "Clinical Subgroup"
  ) +
  labs(
    title = "Model Reliability Under Clinical Distribution Shifts",
    x = NULL,
    y = NULL
  ) +
  theme_bw(base_size = 14) +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 11),
    strip.text = element_text(size = 11),
    plot.title = element_text(hjust = 0.5, size = 15),
    plot.subtitle = element_text(hjust = 0.5, size = 12, color = "gray30"),
    panel.spacing.x = unit(2, "lines"),
    axis.text.y = element_text(size = 11),
    plot.margin = ggplot2::margin(10, 30, 10, 10, unit = "pt")
  )

print(shift_plot)
# Save
ggsave(
  "Figure_Novel5_distribution_shifts.png", shift_plot, width = 16,
  height = 6, dpi = 300)

