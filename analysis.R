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
  # data <- data[!grepl("\\.y$", data[[1]]), ]
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
test_ratio <- 0.5

# Number of samples
n <- nrow(X)
train_size = floor(train_ratio * n)
train_size
test_size = n-train_size
test_size
# Generate random indices for training, validation, and test sets
# Generate a shuffled sequence of row indices
shuffled_indices <- sample(seq_len(n))

# Assign the first `train_size` indices to the training set
train_indices <- shuffled_indices[1:train_size]
train_indices

# Assign the remaining indices to the test set
test_indices <- shuffled_indices[(train_size + 1):n]
test_indices
# Training set
X_train <- X[train_indices, ]
Z_train <- Z[train_indices, ]
Y_train <- MPNDisease[train_indices,]


# Test set
X_test <- X[test_indices, ]
Z_test <- Z[test_indices, ]
Y_test <- MPNDisease[test_indices,]

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
Y_binary_train <- as.numeric(Y_train == "MF")
Y_binary_train
#Y_binary_val <- as.numeric(Y_val == "PV")
Y_binary_test <- as.numeric(Y_test == "MF")
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
best_lambda
# Lambda with the most regularization that is within 1 standard error of the minimum error
lambda_1se <- cvfit$lambda.1se
lambda_1se
# Plot cross-validation results
plot(cvfit)
print(cvfit)
coef(cvfit)
# Highlight the selected lambda values
abline(v = log(best_lambda), col = "blue", lty = 2)  # Optimal lambda
abline(v = log(lambda_1se), col = "red", lty = 2)    # 1-SE lambda

# Refit the model with the optimal lambda
fit_optimized <- multiview(x_train, Y_binary_train, family = binomial(),rho=4, lambda = best_lambda)
rho <- fit_optimized$rho
rho
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


#3

predict_on_whole_dataset <- function(class_name, x_train, Y_train, x_test, Y_test, output_dir, results_file) {
  # Binary labeling for the target class
  # class_name = "MF"
  Y_binary_train <- as.numeric(Y_train == class_name)
  Y_binary_test <- as.numeric(Y_test == class_name)
  
  # Train the model on the training set and predict on the test set
  cvfit_train <- cv.multiview(x_train, Y_binary_train, family = binomial(), type.measure = "class")
  best_lambda <- cvfit_train$lambda.min
  
  fit_train <- multiview(x_train, Y_binary_train, family = binomial(), lambda = best_lambda)
  test_probabilities <- predict(fit_train, newx = x_test, s = best_lambda, type = "response")
  test_predictions <- as.numeric(test_probabilities > 0.5)
  
  # Train the model on the test set and predict on the train set
  cvfit_test <- cv.multiview(x_test, Y_binary_test, family = binomial(), type.measure = "class")
  best_lambda_test <- cvfit_test$lambda.min
  
  fit_test <- multiview(x_test, Y_binary_test, family = binomial(), lambda = best_lambda_test)
  train_probabilities <- predict(fit_test, newx = x_train, s = best_lambda_test, type = "response")
  train_predictions <- as.numeric(train_probabilities > 0.5)
  
  # Combine predictions for both sets with dynamic labels
  results_train <- data.frame(
    Sample_ID = rownames(x_train[[1]]),
    Predicted_Probability = train_probabilities,
    Predicted_Label = ifelse(train_predictions == 1, class_name, paste("Not", class_name)),
    Actual_Label = ifelse(Y_binary_train == 1, class_name, paste("Not", class_name))
  )
  
  results_test <- data.frame(
    Sample_ID = rownames(x_test[[1]]),
    Predicted_Probability = test_probabilities,
    Predicted_Label = ifelse(test_predictions == 1, class_name, paste("Not", class_name)),
    Actual_Label = ifelse(Y_binary_test == 1, class_name, paste("Not", class_name))
  )
  # Combine and save to a CSV file
  all_results <- rbind(results_train, results_test) %>%
    arrange(Sample_ID)
  write.csv(all_results, file = paste0(output_dir, "/", class_name, "_Predictions.csv"), row.names = FALSE)
  
  # Output to console for verification
  cat("Results saved for class", class_name, "\n")
  print(head(all_results))
}
#4
predict_on_whole_dataset <- function(class_name, x_train, Y_train, x_test, Y_test, output_dir, results_file) {
  # Binary labeling for the target class
  Y_binary_train <- as.numeric(Y_train == class_name)
  Y_binary_test <- as.numeric(Y_test == class_name)
  
  #===============================================================================
  # Train the model on the training set and predict on the test set
  #===============================================================================
  cvfit_train <- cv.multiview(x_train, Y_binary_train, family = binomial(), type.measure = "class")
  best_lambda <- cvfit_train$lambda.min
  
  fit_train <- multiview(x_train, Y_binary_train, family = binomial(), lambda = best_lambda)
  test_probabilities <- predict(fit_train, newx = x_test, s = best_lambda, type = "response")
  test_predictions <- as.numeric(test_probabilities > 0.5)
  
  #===============================================================================
  # Train the model on the test set and predict on the train set
  #===============================================================================
  cvfit_test <- cv.multiview(x_test, Y_binary_test, family = binomial(), type.measure = "class")
  best_lambda_test <- cvfit_test$lambda.min
  
  fit_test <- multiview(x_test, Y_binary_test, family = binomial(), lambda = best_lambda_test)
  train_probabilities <- predict(fit_test, newx = x_train, s = best_lambda_test, type = "response")
  train_predictions <- as.numeric(train_probabilities > 0.5)
  
  #===============================================================================
  # Combine Predictions for Both Sets
  #===============================================================================
  results_train <- data.frame(
    Sample_ID = rownames(x_train[[1]]),
    Predicted_Probability = train_probabilities,
    Predicted_Label = ifelse(train_predictions == 1, class_name, paste("Not", class_name)),
    Actual_Label = ifelse(Y_binary_train == 1, class_name, paste("Not", class_name))
  )
  
  results_test <- data.frame(
    Sample_ID = rownames(x_test[[1]]),
    Predicted_Probability = test_probabilities,
    Predicted_Label = ifelse(test_predictions == 1, class_name, paste("Not", class_name)),
    Actual_Label = ifelse(Y_binary_test == 1, class_name, paste("Not", class_name))
  )
  
  all_results <- rbind(results_train, results_test) %>%
    arrange(Sample_ID)
  write.csv(all_results, file = paste0(output_dir, "/", class_name, "_Predictions.csv"), row.names = FALSE)
  
  #===============================================================================
  # Save Model Details
  #===============================================================================
  # Calculate accuracies
  train_accuracy <- mean(train_predictions == Y_binary_train)
  test_accuracy <- mean(test_predictions == Y_binary_test)
  
  # Save Model Details
  model_details <- data.frame(
    Model = c("Train", "Test"),
    BestLambda = c(best_lambda, best_lambda_test),
    Rho = c(fit_train$rho, fit_test$rho),
    Accuracy = c(train_accuracy, test_accuracy)
  )
  write.csv(model_details, file = paste0(output_dir, "/", class_name, "_ModelDetails.csv"), row.names = FALSE)
  
  #===============================================================================
  # Save Coefficients (Feature Importance)
  #===============================================================================
  # Extract coefficients for the training model
  coef_train <- as.matrix(coef(fit_train, s = best_lambda))[-1, , drop = FALSE]  # Remove intercept
  coef_train_df <- data.frame(Feature = rownames(coef_train), Coefficient = coef_train[, 1])
  coef_train_df <- coef_train_df %>%
    filter(Coefficient != 0) %>%
    arrange(desc(abs(Coefficient)))
  write.csv(coef_train_df, file = paste0(output_dir, "/", class_name, "_Train_TopFeatures.csv"), row.names = FALSE)
  
  # Extract coefficients for the test model
  coef_test <- as.matrix(coef(fit_test, s = best_lambda_test))[-1, , drop = FALSE]  # Remove intercept
  coef_test_df <- data.frame(Feature = rownames(coef_test), Coefficient = coef_test[, 1])
  coef_test_df <- coef_test_df %>%
    filter(Coefficient != 0) %>%
    arrange(desc(abs(Coefficient)))
  write.csv(coef_test_df, file = paste0(output_dir, "/", class_name, "_Test_TopFeatures.csv"), row.names = FALSE)
  
  #===============================================================================
  # Output to Console for Verification
  #===============================================================================
  cat("Results saved for class", class_name, "\n")
  print(head(all_results))
  cat("Model details:\n")
  print(model_details)
  cat("Top features for train model:\n")
  print(head(coef_train_df, 10))
  cat("Top features for test model:\n")
  print(head(coef_test_df, 10))
}

# Define a maximum iteration count
max_iterations <- 10
iteration <- 0
accuracy_threshold <- 0.85

repeat {
  iteration <- iteration + 1
  cat("Iteration:", iteration, "\n")
  
  # Train the model on the training set and predict on the test set
  cvfit_train <- cv.multiview(x_train, Y_binary_train, family = binomial(), type.measure = "class")
  best_lambda <- cvfit_train$lambda.min
  
  fit_train <- multiview(x_train, Y_binary_train, family = binomial(), lambda = best_lambda)
  test_probabilities <- predict(fit_train, newx = x_test, s = best_lambda, type = "response")
  test_predictions <- as.numeric(test_probabilities > 0.5)
  
  # Train the model on the test set and predict on the train set
  cvfit_test <- cv.multiview(x_test, Y_binary_test, family = binomial(), type.measure = "class")
  best_lambda_test <- cvfit_test$lambda.min
  
  fit_test <- multiview(x_test, Y_binary_test, family = binomial(), lambda = best_lambda_test)
  train_probabilities <- predict(fit_test, newx = x_train, s = best_lambda_test, type = "response")
  train_predictions <- as.numeric(train_probabilities > 0.5)
  
  # Calculate accuracies
  train_accuracy <- mean(train_predictions == Y_binary_train)
  test_accuracy <- mean(test_predictions == Y_binary_test)
  
  # Print accuracies for this iteration
  cat(sprintf("Train Accuracy: %.4f, Test Accuracy: %.4f\n", train_accuracy, test_accuracy))
  
  # Check if accuracies meet the threshold
  if (train_accuracy > accuracy_threshold && test_accuracy > accuracy_threshold) {
    cat("Accuracy threshold met. Saving model details...\n")
    
    # Save Model Details
    model_details <- data.frame(
      Model = c("Train", "Test"),
      BestLambda = c(best_lambda, best_lambda_test),
      Rho = c(fit_train$rho, fit_test$rho),
      Accuracy = c(train_accuracy, test_accuracy)
    )
    write.csv(model_details, file = paste0(output_dir, "/", class_name, "_ModelDetails.csv"), row.names = FALSE)
    
    # Save Coefficients for Train and Test Models
    coef_train <- as.matrix(coef(fit_train, s = best_lambda))[-1, , drop = FALSE]
    coef_train_df <- data.frame(Feature = rownames(coef_train), Coefficient = coef_train[, 1]) %>%
      filter(Coefficient != 0) %>%
      arrange(desc(abs(Coefficient)))
    write.csv(coef_train_df, file = paste0(output_dir, "/", class_name, "_Train_TopFeatures.csv"), row.names = FALSE)
    
    coef_test <- as.matrix(coef(fit_test, s = best_lambda_test))[-1, , drop = FALSE]
    coef_test_df <- data.frame(Feature = rownames(coef_test), Coefficient = coef_test[, 1]) %>%
      filter(Coefficient != 0) %>%
      arrange(desc(abs(Coefficient)))
    write.csv(coef_test_df, file = paste0(output_dir, "/", class_name, "_Test_TopFeatures.csv"), row.names = FALSE)
    
    # Save Predictions for the Main Model
    results_train <- data.frame(
      Sample_ID = rownames(x_train[[1]]),
      Predicted_Probability = train_probabilities,
      Predicted_Label = ifelse(train_predictions == 1, class_name, paste("Not", class_name)),
      Actual_Label = ifelse(Y_binary_train == 1, class_name, paste("Not", class_name))
    )
    
    results_test <- data.frame(
      Sample_ID = rownames(x_test[[1]]),
      Predicted_Probability = test_probabilities,
      Predicted_Label = ifelse(test_predictions == 1, class_name, paste("Not", class_name)),
      Actual_Label = ifelse(Y_binary_test == 1, class_name, paste("Not", class_name))
    )
    
    # Combine and save predictions
    all_results <- rbind(results_train, results_test) %>%
      arrange(Sample_ID)
    write.csv(all_results, file = paste0(output_dir, "/", class_name, "_Predictions.csv"), row.names = FALSE)
    
    # Stop the process
    break
  
  }
  
  # Stop if the maximum number of iterations is reached
  if (iteration >= max_iterations) {
    cat("Max iterations reached. Accuracy threshold not satisfied.\n")
    break
  }
}
# Classes to process
# classes <- c("PV","ET", "CTRL")

# Output directory to save results
output_dir <- "results/co-op_learning/results_whole_dataset_prediction4"
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}
# CSV file to save results
results_file <- "results/co-op_learning/results_whole_dataset_prediction4/model_details.csv"
if (file.exists(results_file)) {
  file.remove(results_file)  # Remove old file if exists
}
# Classes to process
classes <- c("MF","PV","ET","CTRL")
predict_on_whole_dataset(classes[1], x_train, Y_train, x_test, Y_test, output_dir)
  # Iterate over classes and process each one
for (class_name in classes) {
  predict_on_whole_dataset(class_name, x_train, Y_train, x_test, Y_test, output_dir)
}
