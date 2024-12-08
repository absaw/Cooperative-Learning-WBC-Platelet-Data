# author : abhishek Sawalkar
# importing libraries
library(dplyr)
library(readr)
library(multiview)
library(ggplot2)
setwd("~/Library/CloudStorage/Box-Box/KrishnanA-Stats-Share/LucyExtendedTermProjects/Project5-WBC_RNAseq-Analysis/WBC-Platelet-Coop-Learning-AS/Cooperative-Learning-WBC-Platelet-Data")
#===============================================================================
# PART 1 :: Pre-processing data
#===============================================================================
# Load the datasets
# Load the datasets
# Load data
platelet_data <- read.csv("data/NormCountsPlatelet.csv", stringsAsFactors = FALSE)
wbc_data <- read.csv("data/NormCountsWBC.csv", stringsAsFactors = FALSE)
outcome_data <- read.csv("data/DeSeq2_AllData_n120_PTPRCadj.csv", stringsAsFactors = FALSE)

# Function to clean sample names
clean_sample_names <- function(data) {
  # Remove rows with ".y" suffix in the sample name
  data <- data[!grepl("\\.y$", data[[1]]), ]
  # Set row names and remove ".x" suffix
  rownames(data) <- gsub("\\.x$", "", data[[1]])
  data <- data[, -1]  # Remove the first column
  return(data)
}

# Apply cleaning function to each dataset
platelet_data <- clean_sample_names(platelet_data)
wbc_data <- clean_sample_names(wbc_data)
outcome_data <- clean_sample_names(outcome_data)

# Verify changes
cat("Updated sample names in platelet_data:\n", head(rownames(platelet_data)), "\n")
cat("Updated sample names in wbc_data:\n", head(rownames(wbc_data)), "\n")
cat("Updated sample names in outcome_data:\n", head(rownames(outcome_data)), "\n")
# Extract patient IDs
platelet_patients <- rownames(platelet_data)
wbc_patients <- rownames(wbc_data)
outcome_patients <- rownames(outcome_data)

# Find common patients
common_patients <- intersect(intersect(platelet_patients, wbc_patients), outcome_patients)
cat("Number of common patients:", length(common_patients), "\n")




# platelet_data <- read.csv("data/NormCountsPlatelet.csv")
# wbc_data <- read_csv("data/NormCountsWBC.csv")
# outcome_data <- read_csv("data/DeSeq2_AllData_n120_PTPRCadj.csv")
# 
# # Remove ".x" suffix from sample names in X and Z
# rownames(platelet_data) <- gsub("\\.x$", "", rownames(platelet_data))  # For X (Platelet data)
# rownames(wbc_data) <- gsub("\\.x$", "", rownames(wbc_data))  # For Z (WBC data)
# rownames(outcome_data) <- gsub("\\.x$", "", rownames(outcome_data)) 
# # Verify the changes
# cat("Updated sample names in X:\n", head(rownames()), "\n")
# cat("Updated sample names in Z:\n", head(rownames(Z)), "\n")
# # Extract patient IDs
# platelet_patients <- platelet_data[[1]]
# wbc_patients <- wbc_data[[1]]
# outcome_patients <- outcome_data[[1]]
# # Find common patients
# common_patients <- intersect(intersect(platelet_patients, wbc_patients), outcome_patients)
# length(common_patients)


# Subset the data for common patients
X <- platelet_data[common_patients, , drop = FALSE]  # Subset rows using rownames
Z <- wbc_data[common_patients, , drop = FALSE]       # Subset rows using rownames
Y <- outcome_data[common_patients, , drop = FALSE]   # Subset rows using rownames

# Sort the datasets by rownames to ensure consistent order
X <- X[order(rownames(X)), ]
Z <- Z[order(rownames(Z)), ]
Y <- Y[order(rownames(Y)), ]

# Extract the MPNDisease column
MPNDisease <- Y[, 1, drop = FALSE]  # Keep it as a data frame

# Rename columns to distinguish views
colnames(X) <- paste0("X_", colnames(X))  # Add "X_" prefix to platelet data
colnames(Z) <- paste0("Z_", colnames(Z))  # Add "Z_" prefix to WBC data
#MPNDisease<-Y %>% select(1,2)
# Combine X, Z, and MPNDisease into a single data frame for splitting
# Normalize X and Z datasets
X_normalized <- scale(X)  # Normalize platelet data
Z_normalized <- scale(Z)  # Normalize WBC data
X <- X_normalized
Z <- Z_normalized
#combined_data <- data.frame(X, Z, MPNDisease = MPNDisease[[1]])

#===============================================================================
# Random sampling
#===============================================================================
# Set seed for reproducibility
set.seed(123)

# Define split ratios
train_ratio <- 0.5
val_ratio <- 0
test_ratio <- 0.5

# Number of samples
n <- nrow(X)
train_size = floor(train_ratio * n)
train_size
val_size= floor(val_ratio*n)
val_size
test_size = n-train_size-val_size
test_size
# Generate random indices for training, validation, and test sets
# Generate a shuffled sequence of row indices
shuffled_indices <- sample(seq_len(n))

# Assign the first `train_size` indices to the training set
train_indices <- shuffled_indices[1:train_size]

# Assign the next `val_size` indices to the validation set
#val_indices <- shuffled_indices[(train_size + 1):(train_size + val_size)]

# Assign the remaining indices to the test set
test_indices <- shuffled_indices[(train_size + val_size + 1):n]
# For stratified sampling
library(caret)  
set.seed(123)
# Stratified sampling to preserve class proportions
n <- nrow(X)
MPNDisease <- Y[, 1, drop = TRUE]
train_indices <- createDataPartition(MPNDisease, p = train_ratio, list = FALSE)
test_indices <- setdiff(seq_len(n), train_indices)
# Training set
X_train <- X[train_indices, ]
Z_train <- Z[train_indices, ]
Y_train <- MPNDisease[train_indices]

# Validation set
#X_val <- X[val_indices, ]
#Z_val <- Z[val_indices, ]
#Y_val <- MPNDisease[val_indices, ]

# Test set
X_test <- X[test_indices, ]
Z_test <- Z[test_indices, ]
Y_test <- MPNDisease[test_indices]

# Prepare data for multiview package
# Convert MPNDisease to a factor
Y_train <- as.factor(Y_train)
#Y_val <- as.factor(Y_val[[1]])
Y_test <- as.factor(Y_test)

# Create predictor lists for the multiview package
x_train <- list(X_train, Z_train)
#x_val <- list(X_val, Z_val)
x_test <- list(X_test, Z_test)

# Check class distributions
table(Y_train)
#table(Y_val)
table(Y_test)

#===============================================================================
# PART 2 :: Model Fit
#===============================================================================

# multiview package doesnt support mutlinomial classification
# so for our problem, with 4 different subtypes possible, one v/s all method
# needs to be used. Hence binomial() family then can be used for binary
# classifiaction

# Prepare Binary Labels for PV
# Create binary labels: 1 for PV, 0 for others
Y_binary_train <- as.numeric(Y_train == "ET")
Y_binary_train
#Y_binary_val <- as.numeric(Y_val == "PV")
Y_binary_test <- as.numeric(Y_test == "ET")
Y_binary_test
# Check the distribution of binary labels
cat("Training set distribution:\n")
table(Y_binary_train)
cat("Validation set distribution:\n")
#table(Y_binary_val)
cat("Test set distribution:\n")
table(Y_binary_test)

# Perform cross-validation to find the optimal lambda
cvfit <- cv.multiview(x_train, Y_binary_train, family = binomial(), type.measure = "class")
best_lambda <- cvfit$lambda.min

# Lambda with the most regularization that is within 1 standard error of the minimum error
lambda_1se <- cvfit$lambda.1se

# Plot cross-validation results
plot(cvfit)
# Highlight the selected lambda values
abline(v = log(best_lambda), col = "blue", lty = 2)  # Optimal lambda
abline(v = log(lambda_1se), col = "red", lty = 2)    # 1-SE lambda

# Refit the model with the optimal lambda
fit_optimized <- multiview(x_train, Y_binary_train, family = binomial(), lambda = best_lambda)
best_rho <- fit_optimized$rho
best_rho
# Optimal lambda minimizing cross-validation error
# Make predictions on the test set
test_probabilities <- predict(fit_optimized, newx = x_test, s = best_lambda, type = "response")

# Convert probabilities to binary predictions
test_predictions <- as.numeric(test_probabilities > 0.5)

# Calculate test accuracy
test_accuracy <- mean(test_predictions == Y_binary_test)
cat("Test Accuracy with Optimized Lambda:", test_accuracy, "\n")
# Confusion Matrix for Test Set
cat("Confusion Matrix for Test Set:\n")
print(table(Predicted = test_predictions, Actual = Y_binary_test))
# Analyze Feature Importance
# Extract coefficients for X (platelet data) and Z (WBC data)
coef <- coef(fit_optimized, s = best_lambda)
coef_matrix <- as.matrix(coef)
# Remove the intercept (row 1)
coef_matrix <- coef_matrix[-1, , drop = FALSE]
# Get feature names and their coefficients
feature_names <- rownames(coef_matrix)
coefficients <- as.numeric(coef_matrix)
# Combine into a data frame for sorting
coef_df <- data.frame(Feature = feature_names, Coefficient = coefficients)
# Rank by absolute value of the coefficients
top_features <- coef_df %>%
  filter(Coefficient != 0) %>%
  arrange(desc(abs(Coefficient))) %>%
  head(20)
# Display the top 20 features
print(top_features)

# Create a bar plot of the top 10 features
plt <-ggplot(top_features, aes(x = reorder(Feature, Coefficient), y = Coefficient)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +
  labs(title = "PV : Top 20 Features by Coefficient", x = "Feature", y = "Coefficient") +
  theme_minimal()
plt <- plt + theme(
  panel.background = element_rect(fill = "white"),  # Light background
  plot.background = element_rect(fill = "white"),   # Light outer background
  axis.text = element_text(color = "black"),        # Black axis text
  axis.title = element_text(color = "black"),       # Black axis titles
  plot.title = element_text(color = "black", face = "bold")  # Black title
)
ggsave('images/PVTopFeatures.png')

# Define rho (agreement penalty)
rho <- 0.3

# Prepare input data
view_list <- list(X = X_train, Z = Z_train)  # Platelet (X) and WBC (Z) data

# Evaluate view contribution
contribution <- view.contribution(
  x_list = view_list,      # List of views
  y = Y_binary_train,      # Target variable
  rho = rho,               # Agreement penalty
  family = binomial(),     # Model family (binomial for classification)
  eval_data = list(X_train, Z_train)  # Optional: Data for evaluation
)

# Print contribution results
print(contribution)
# Save the Model
# Save the fitted model for later use
#saveRDS(fit_pv, file = "model/pv_cooperative_learning_model.rds")

# To load the model in future:
# fit_pv <- readRDS("pv_cooperative_learning_model.rds")

# rho_grid <- c(0.01, 0.1, 0.5, 1)
# 
# results <- lapply(rho_grid, function(rho) {
#   cv.multiview(x_train, Y_binary_train, family = binomial(), rho = rho, type.measure = "class")
# })
# 
# # Extract the best rho and corresponding performance
# errors <- sapply(results, function(cvfit) min(cvfit$cvm))
# best_rho <- rho_grid[which.min(errors)]
# cat("Best rho:", best_rho, "\n")
#Debugging 
cat("Number of rows in X_test:", nrow(X_test), "\n")
cat("Number of rows in Z_test:", nrow(Z_test), "\n")
cat("Length of Y_binary_test:", length(Y_binary_test), "\n")
cat("Length of test_predictions:", length(test_predictions), "\n")



# Function for all classes
# Define a function to process a given class
process_class <- function(class_name, x_train, Y_train, x_test, Y_test, output_dir) {
  # Binary labeling for the target class
  Y_binary_train <- as.numeric(Y_train == class_name)
  Y_binary_test <- as.numeric(Y_test == class_name)
  
  # Check the distribution of binary labels
  cat("Training set distribution for", class_name, ":\n")
  print(table(Y_binary_train))
  cat("Test set distribution for", class_name, ":\n")
  print(table(Y_binary_test))
  
  # Perform cross-validation to find the optimal lambda
  cvfit <- cv.multiview(x_train, Y_binary_train, family = binomial(), type.measure = "class")
  best_lambda <- cvfit$lambda.min
  lambda_1se <- cvfit$lambda.1se
  
  # Plot cross-validation results
  cvplot<-plot(cvfit)+
      abline(v = log(best_lambda), col = "blue", lty = 2)+
      abline(v = log(lambda_1se), col = "red", lty = 2)    # 1-SE lambda
  
  # Save the plot
  ggsave(filename = paste0(output_dir, "/", class_name, "_CVPlot.png"))
  
  # Refit the model with the optimal lambda
  fit_optimized <- multiview(x_train, Y_binary_train, family = binomial(), lambda = best_lambda)
  
  # Extract the best rho
  best_rho <- fit_optimized$rho
  cat("Best rho for", class_name, ":", best_rho, "\n")
  
  # Make predictions on the test set
  test_probabilities <- predict(fit_optimized, newx = x_test, s = best_lambda, type = "response")
  test_predictions <- as.numeric(test_probabilities > 0.5)
  
  # Calculate test accuracy
  test_accuracy <- mean(test_predictions == Y_binary_test)
  cat("Test Accuracy for", class_name, ":", test_accuracy, "\n")
  # Collect model details
  model_details <- data.frame(
    Class = class_name,
    BestLambda = best_lambda,
    Lambda1SE = lambda_1se,
    BestRho = best_rho,
    TestAccuracy = test_accuracy
  )
  
  # Save or append results to the CSV file
  if (!file.exists(results_file)) {
    write.csv(model_details, results_file, row.names = FALSE)  # Create a new file
  } else {
    write.table(model_details, results_file, row.names = FALSE, col.names = FALSE, sep = ",", append = TRUE)  # Append to the file
  }
  # Confusion Matrix
  cat("Confusion Matrix for Test Set for", class_name, ":\n")
  print(table(Predicted = test_predictions, Actual = Y_binary_test))
  
  # Analyze Feature Importance
  coef <- coef(fit_optimized, s = best_lambda)
  coef_matrix <- as.matrix(coef)[-1, , drop = FALSE]  # Remove intercept
  feature_names <- rownames(coef_matrix)
  coefficients <- as.numeric(coef_matrix)
  coef_df <- data.frame(Feature = feature_names, Coefficient = coefficients)
  
  # Get top 20 features
  top_features <- coef_df %>%
    filter(Coefficient != 0) %>%
    arrange(desc(abs(Coefficient))) %>%
    head(20)
  cat("Top Features for", class_name, ":\n")
  print(top_features)
  
  # Create a bar plot of the top 20 features
  plt <- ggplot(top_features, aes(x = reorder(Feature, Coefficient), y = Coefficient)) +
    geom_bar(stat = "identity", fill = "steelblue") +
    coord_flip() +
    labs(title = paste0(class_name, ": Top Features by Coefficient"), x = "Feature", y = "Coefficient") +
    theme_minimal()+ 
    theme(
      panel.background = element_rect(fill = "white"),  # Light background
      plot.background = element_rect(fill = "white"),   # Light outer background
      axis.text = element_text(color = "black"),        # Black axis text
      axis.title = element_text(color = "black"),       # Black axis titles
      plot.title = element_text(color = "black", face = "bold")  # Black title
    )
  
  # Display and save the plot
  plt
  print(plt)
  ggsave(filename = paste0(output_dir, "/", class_name, "_TopFeatures.png"))
}

# 2 ======================
process_class <- function(class_name, x_train, Y_train, x_test, Y_test, output_dir) {
  # Binary labeling for the target class
  Y_binary_train <- as.numeric(Y_train == class_name)
  Y_binary_test <- as.numeric(Y_test == class_name)
  
  # Check the distribution of binary labels
  cat("Training set distribution for", class_name, ":\n")
  print(table(Y_binary_train))
  cat("Test set distribution for", class_name, ":\n")
  print(table(Y_binary_test))
  
  # Perform cross-validation to find the optimal lambda
  cvfit <- cv.multiview(x_train, Y_binary_train, family = binomial(), type.measure = "class")
  best_lambda <- cvfit$lambda.min
  lambda_1se <- cvfit$lambda.1se
  
  # Plot cross-validation results
png(filename = paste0(output_dir, "/", class_name, "_CVPlot.png"), width = 800, height = 600)
plot(cvfit)
abline(v = log(best_lambda), col = "blue", lty = 2)
abline(v = log(lambda_1se), col = "red", lty = 2)
dev.off()  # Close the graphics device
  
  # Refit the model with the optimal lambda
  fit_optimized <- multiview(x_train, Y_binary_train, family = binomial(), lambda = lambda_1se)
  
  # Extract the best rho
  best_rho <- fit_optimized$rho
  cat("Best rho for", class_name, ":", best_rho, "\n")
  
  # Make predictions on the test set
  test_probabilities <- predict(fit_optimized, newx = x_test, s = lambda_1se, type = "response")
  test_predictions <- as.numeric(test_probabilities > 0.5)
  
  # Calculate test accuracy
  test_accuracy <- mean(test_predictions == Y_binary_test)
  cat("Test Accuracy for", class_name, ":", test_accuracy, "\n")
  
  # Collect model details
  model_details <- data.frame(
    Class = class_name,
    BestLambda = best_lambda,
    Lambda1SE = lambda_1se,
    BestRho = best_rho,
    TestAccuracy = test_accuracy
  )
  
  # Save or append results to the CSV file
  if (!file.exists(results_file)) {
    write.csv(model_details, results_file, row.names = FALSE)  # Create a new file
  } else {
    write.table(model_details, results_file, row.names = FALSE, col.names = FALSE, sep = ",", append = TRUE)  # Append to the file
  }
  
  # Confusion Matrix
  cat("Confusion Matrix for Test Set for", class_name, ":\n")
  print(table(Predicted = test_predictions, Actual = Y_binary_test))
  
  # Plot and save the confusion matrix
  confusion_data <- as.data.frame(table(Predicted = test_predictions, Actual = Y_binary_test))
  confusion_plot <- ggplot(confusion_data, aes(x = Actual, y = Predicted)) +
    geom_tile(aes(fill = Freq), color = "white") +
    geom_text(aes(label = Freq), size = 5) +
    scale_fill_gradient(low = "white", high = "blue") +
    labs(title = paste("Confusion Matrix for", class_name), x = "Actual", y = "Predicted") +
    theme_minimal()
  
  ggsave(
    filename = paste0(output_dir, "/", class_name, "_ConfusionMatrix.png"),
    plot = confusion_plot, dpi = 300, width = 8, height = 6
  )
  
  # Analyze Feature Importance
  coef <- coef(fit_optimized, s = best_lambda)
  coef_matrix <- as.matrix(coef)[-1, , drop = FALSE]  # Remove intercept
  feature_names <- rownames(coef_matrix)
  coefficients <- as.numeric(coef_matrix)
  coef_df <- data.frame(Feature = feature_names, Coefficient = coefficients)
  
  # Get top 20 features
  top_features <- coef_df %>%
    filter(Coefficient != 0) %>%
    arrange(desc(abs(Coefficient))) %>%
    head(20)
  cat("Top Features for", class_name, ":\n")
  print(top_features)
  
  # Save top features to a CSV file
  write.csv(
    top_features,
    file = paste0(output_dir, "/", class_name, "_TopFeatures.csv"),
    row.names = FALSE
  )
  
  # Create a bar plot of the top 20 features
  plt <- ggplot(top_features, aes(x = reorder(Feature, Coefficient), y = Coefficient)) +
    geom_bar(stat = "identity", fill = "steelblue") +
    coord_flip() +
    labs(title = paste0(class_name, ": Top Features by Coefficient"), x = "Feature", y = "Coefficient") +
    theme_minimal() + 
    theme(
      panel.background = element_rect(fill = "white"),  # Light background
      plot.background = element_rect(fill = "white"),   # Light outer background
      axis.text = element_text(color = "black"),        # Black axis text
      axis.title = element_text(color = "black"),       # Black axis titles
      plot.title = element_text(color = "black", face = "bold")  # Black title
    )
  
  # Save the plot
  ggsave(
    filename = paste0(output_dir, "/", class_name, "_TopFeatures.png"),
    plot = plt, dpi = 300, width = 8, height = 6
  )
}

# Classes to process
classes <- c("PV","ET", "MF", "CTRL")

# Output directory to save results
output_dir <- "results_stratified_sampled"
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}
# CSV file to save results
results_file <- "results_stratified_sampled/model_details.csv"
if (file.exists(results_file)) {
  file.remove(results_file)  # Remove old file if exists
}
process_class(classes[4], x_train, Y_train, x_test, Y_test, output_dir)
# Iterate over classes and process each one
for (class_name in classes) {
  process_class(class_name, x_train, Y_train, x_test, Y_test, output_dir)
}



# Classes to process
classes <- c("PV", "ET", "MF", "CTRL")

# Initialize a matrix to store probabilities for all classes
test_probabilities_matrix <- matrix(0, nrow = length(Y_test), ncol = length(classes))
colnames(test_probabilities_matrix) <- classes

# Iterate over classes to train binary classifiers and get probabilities
for (i in seq_along(classes)) {
  class_name <- classes[i]
  
  # Binary labeling for the target class
  Y_binary_train <- as.numeric(Y_train == class_name)
  Y_binary_test <- as.numeric(Y_test == class_name)
  
  # Perform cross-validation to find the optimal lambda
  cvfit <- cv.multiview(x_train, Y_binary_train, family = binomial(), type.measure = "class")
  best_lambda <- cvfit$lambda.min
  
  # Refit the model with the optimal lambda
  fit_optimized <- multiview(x_train, Y_binary_train, family = binomial(), lambda = best_lambda)
  
  # Get probabilities for the test set
  test_probabilities <- predict(fit_optimized, newx = x_test, s = best_lambda, type = "response")
  
  # Store probabilities for this class
  test_probabilities_matrix[, i] <- test_probabilities
}

# Convert probabilities to final predictions
# Assign each sample to the class with the highest probability
final_predictions <- colnames(test_probabilities_matrix)[max.col(test_probabilities_matrix)]

# Calculate accuracy of multinomial classification
multinomial_accuracy <- mean(final_predictions == Y_test)
cat("Multinomial Classification Accuracy:", multinomial_accuracy, "\n")

# Confusion matrix for multinomial classification
multinomial_confusion_matrix <- table(Predicted = final_predictions, Actual = Y_test)
cat("Confusion Matrix for Multinomial Classification:\n")
print(multinomial_confusion_matrix)
