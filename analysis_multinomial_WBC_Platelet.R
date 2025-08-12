# author: Abhishek Sawalkar (adapted)
# lasso_multinomial_analysis.R
# Description: keep preprocessing same as original script, then run multinomial LASSO (glmnet).
# ------------------------------------------------------------------------------

# Libraries
library(dplyr)
library(readr)
library(glmnet)
library(pROC)
library(ggplot2)

# ---------------------------------------------------------------------------
# Working directory (same as your original script)
# ---------------------------------------------------------------------------
setwd("~/Library/CloudStorage/Box-Box/KrishnanA-Stats-Share/LucyExtendedTermProjects/Project5-WBC_RNAseq-Analysis/WBC-Platelet-Coop-Learning-AS/Cooperative-Learning-WBC-Platelet-Data")
output_dir <- "results/lasso/multinomial_platelet_wbc"
dir.create(output_dir, showWarnings = FALSE)
# ---------------------------------------------------------------------------
# Load & preprocess (keeps your original preprocessing)
# ---------------------------------------------------------------------------
platelet_data <- read.csv("data/NormCountsPlateletAdjusted.csv", stringsAsFactors = FALSE)
wbc_data      <- read.csv("data/NormCountsWBCAdjusted.csv", stringsAsFactors = FALSE)
outcome_data  <- read.csv("data/DeSeq2_AllData_n120_PTPRCadj.csv", stringsAsFactors = FALSE)

clean_sample_names <- function(data) {
  # keep behavior same as your original: use first column for rownames, remove ".x" suffix,
  # then drop the first column.
  rownames(data) <- gsub("\\.x$", "", data[[1]])
  data <- data[, -1, drop = FALSE]
  return(data)
}

platelet_data <- clean_sample_names(platelet_data)
wbc_data      <- clean_sample_names(wbc_data)
outcome_data  <- clean_sample_names(outcome_data)

cat("Sample names (platelet):", head(rownames(platelet_data)), "\n")
cat("Sample names (wbc):", head(rownames(wbc_data)), "\n")
cat("Sample names (outcome):", head(rownames(outcome_data)), "\n")

# find common patients and subset
platelet_patients <- rownames(platelet_data)
wbc_patients      <- rownames(wbc_data)
outcome_patients  <- rownames(outcome_data)

common_patients <- intersect(intersect(platelet_patients, wbc_patients), outcome_patients)
cat("Number of common patients:", length(common_patients), "\n")
if (length(common_patients) == 0) stop("No common patients found - check sample names.")

X <- platelet_data[common_patients, , drop = FALSE]
Z <- wbc_data[common_patients, , drop = FALSE]
Y <- outcome_data[common_patients, , drop = FALSE]

# Sort by rownames to ensure consistent order
X <- X[order(rownames(X)), , drop = FALSE]
Z <- Z[order(rownames(Z)), , drop = FALSE]
Y <- Y[order(rownames(Y)), , drop = FALSE]

# Extract disease column (first column, as in your script)
MPNDisease <- Y[, 1, drop = FALSE]

# Prefix columns to preserve view identity
colnames(X) <- paste0("X_", colnames(X))
colnames(Z) <- paste0("Z_", colnames(Z))

# Normalize entire X and Z (same as your script)
X_normalized <- scale(X)
Z_normalized <- scale(Z)
X <- X_normalized
Z <- Z_normalized

# ---------------------------------------------------------------------------
# 50/50 random split with stratified sampling
# ---------------------------------------------------------------------------
set.seed(123)  # reproducible

train_ratio <- 0.5
min_train_per_class <- 8  # adjust if needed

classes <- unique(MPNDisease[, 1])

train_indices <- unlist(
  lapply(split(seq_len(nrow(X)), MPNDisease[, 1]), function(idx) {
    n_cl <- length(idx)
    # ensure at least min_train_per_class and leave at least 1 for test
    n_train_cl <- max(min_train_per_class, floor(n_cl * train_ratio))
    n_train_cl <- min(n_train_cl, n_cl - 1)
    sample(idx, n_train_cl)
  })
)

test_indices <- setdiff(seq_len(nrow(X)), train_indices)

# Create datasets
X_train <- X[train_indices, , drop = FALSE]
Z_train <- Z[train_indices, , drop = FALSE]
X_test  <- X[test_indices, , drop = FALSE]
Z_test  <- Z[test_indices,  , drop = FALSE]

# Keep factor levels consistent for multinomial classification
Y_train <- factor(MPNDisease[train_indices, 1], levels = classes)
Y_test  <- factor(MPNDisease[test_indices,  1], levels = classes)

# Print distributions
cat("Class levels in Y_train:", levels(Y_train), "\n")
cat("Y_train distribution:\n"); print(table(Y_train))
cat("Class levels in Y_test:", levels(Y_test), "\n")
cat("Y_test distribution:\n");  print(table(Y_test))







library(glmnet)

# Combine views into design matrices
train_mat <- as.matrix(cbind(X_train, Z_train))
test_mat  <- as.matrix(cbind(X_test,  Z_test))

# Sanity check
stopifnot(nrow(train_mat) == length(Y_train))
stopifnot(nrow(test_mat)  == length(Y_test))

# Train multinomial LASSO logistic regression with cross-validation
set.seed(123)
cvfit <- cv.glmnet(
  x = train_mat,
  y = Y_train,
  family = "multinomial",
  alpha = 1,
  type.measure = "class",
  nfolds = 5
)

best_lambda <- cvfit$lambda.min
best_lambda
# Save model object
saveRDS(cvfit, file = file.path(output_dir, "cvfit_glmnet_multinomial.rds"))

# Predict classes and probabilities on test set
pred_class_test_raw <- predict(cvfit, newx = test_mat, s = best_lambda, type = "class")
predicted_class_test <- as.character(as.vector(pred_class_test_raw))

prob_array <- predict(cvfit, newx = test_mat, s = best_lambda, type = "response")

# Correctly convert to matrix of probs (samples x classes)
prob_matrix <- prob_array[, 1, , drop = TRUE]

# Add row and column names
rownames(prob_matrix) <- rownames(test_mat)
colnames(prob_matrix) <- dimnames(prob_array)[[3]]

print(dim(prob_matrix))  # Should be 35 x 4 for your 4 classes
print(head(prob_matrix))

# Now select classes of interest
classes_of_interest <- c("ET", "PV", "MF")

prob_cols <- lapply(classes_of_interest, function(cl) {
  if (cl %in% colnames(prob_matrix)) {
    prob_matrix[, cl]
  } else {
    rep(NA_real_, nrow(prob_matrix))
  }
})

prob_df_interest <- as.data.frame(prob_cols)
colnames(prob_df_interest) <- classes_of_interest

# Create final dataframe
test_predictions_df <- data.frame(
  Sample_ID = rownames(test_mat),
  prob_df_interest,
  Predicted = predicted_class_test,
  True = as.character(Y_test),
  stringsAsFactors = FALSE
)

print(head(test_predictions_df))

prob_array <- predict(cvfit, newx = test_mat, s = best_lambda, type = "response")
if (length(dim(prob_array)) == 3) {
  prob_matrix <- prob_array[, 1, , drop=TRUE]
  colnames(prob_matrix) <- dimnames(prob_array)[[3]]
} else {
  prob_matrix <- as.matrix(prob_array)
  colnames(prob_matrix) <- levels(Y_train)
}

# Ensure all classes ET, PV, MF exist (CTRL not required for predictions CSV)
classes_of_interest <- c("ET", "PV", "MF")
classes_of_interest <- c("ET", "PV", "MF")
# Reorder columns of prob_matrix to classes_of_interest (fill NA if missing)
prob_cols <- lapply(classes_of_interest, function(cl) {
  if (cl %in% colnames(prob_matrix)) {
    prob_matrix[, cl]
  } else {
    rep(NA_real_, nrow(prob_matrix))
  }
})
prob_df_interest <- as.data.frame(prob_cols)
colnames(prob_df_interest) <- classes_of_interest


# Create final predictions dataframe
test_predictions_df <- data.frame(
  Sample_ID = rownames(test_mat),
  prob_df_interest,
  Predicted = predicted_class_test,
  True = as.character(Y_test),
  stringsAsFactors = FALSE
)

# Save predictions to CSV for ROC plotting
write.csv(test_predictions_df, file = file.path(output_dir, "test_predictions_multinomial.csv"), row.names = FALSE)

# Print a summary of predictions
print(head(test_predictions_df))