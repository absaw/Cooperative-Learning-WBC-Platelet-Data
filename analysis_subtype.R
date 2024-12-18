# Load necessary libraries
library(dplyr)
library(readr)
library(multiview)
library(ggplot2)

#===============================================================================
# PART 1 :: Load and preprocess data
#===============================================================================
# Load the datasets
platelet_data <- read.csv("data/NormCountsPlatelet.csv", stringsAsFactors = FALSE)
wbc_data <- read.csv("data/NormCountsWBC.csv", stringsAsFactors = FALSE)
outcome_data <- read.csv("data/DeSeq2_AllData_n120_PTPRCadj.csv", stringsAsFactors = FALSE)

# Function to clean sample names
clean_sample_names <- function(data) {
  # Remove rows with ".y" suffix in the sample name
  data <- data[!grepl("\\.y$", data[[1]]), ]
  data <- data[!grepl("\\.z$", data[[1]]), ]
  # Set row names and remove ".x" suffix
  rownames(data) <- gsub("\\.x$", "", data[[1]])
  data <- data[, -1]  # Remove the first column
  return(data)
}

clean_sample_names <- function(data) {

  # Remove rows with ".y" suffix in the sample name
  # data <- data[!grepl("\\.y$", data[[1]]), ]
  # data <- data[!grepl("\\.z$", data[[1]]), ]
  # Set row names and remove ".x" suffix
  rownames(data) <- gsub("\\.x$", "", data[[1]])
  data <- data[, -1]  # Remove the first column
  return(data)
}

# Apply cleaning function to each dataset
platelet_data <- clean_sample_names(platelet_data)
wbc_data <- clean_sample_names(wbc_data)
outcome_data <- clean_sample_names(outcome_data)

# Find common patients across datasets
common_patients <- intersect(intersect(rownames(platelet_data), rownames(wbc_data)), rownames(outcome_data))
cat("Number of common patients:", length(common_patients), "\n")

# Subset data for common patients
X <- platelet_data[common_patients, , drop = FALSE]  # Platelet data
Z <- wbc_data[common_patients, , drop = FALSE]       # WBC data
Y <- outcome_data[common_patients, , drop = FALSE]   # Outcome data (MPNDisease)

# Sort datasets to ensure consistent order
X <- X[order(rownames(X)), ]
Z <- Z[order(rownames(Z)), ]
Y <- Y[order(rownames(Y)), ]
# Rename columns to distinguish views
colnames(X) <- paste0("X_", colnames(X))  # Add "X_" prefix to platelet data
colnames(Z) <- paste0("Z_", colnames(Z))  # Add "Z_" prefix to WBC data
# Normalize X and Z datasets
X <- scale(X)
Z <- scale(Z)

#===============================================================================
# PART 2 :: Subtype vs Control Analysis
#===============================================================================
# Function to preprocess data for subtype vs control classification
preprocess_subtype_vs_control <- function(subtype, control = "CTRL", outcome_data, X, Z) {
  # subtype = "MF"
  # control = "CTRL"
  
  # Subset rows for the specific subtype and control
  selected_samples <- rownames(outcome_data)[outcome_data[[1]] %in% c(subtype, control)]
  # Ensure selected_samples exists in both X and Z i.e. common patients
  selected_samples <- intersect(selected_samples, common_patients)
  cat("Selected number of samples for subtype ", subtype," is ", length(selected_samples))
  # selected_samples_2 <- intersect(selected_samples, rownames(Z))
  # Filter the datasets for the selected samples
  X_selected <- X[selected_samples, , drop = FALSE]
  Z_selected <- Z[selected_samples, , drop = FALSE]
  Y_selected <- outcome_data[selected_samples, , drop = FALSE]
  
  # Convert the outcome to binary: 1 for subtype, 0 for control
  Y_selected <- as.numeric(Y_selected[[1]] == subtype)
  
  return(list(X = X_selected, Z = Z_selected, Y = Y_selected))
}

# Function to perform multiview classification for subtype vs control
multiview_subtype_vs_control <- function(subtype, X, Z, Y, output_dir, seed = 123) {
  set.seed(seed)
  
  # Split the data into train-test (50-50 random split)
  n <- nrow(X)
  indices <- sample(seq_len(n))
  train_indices <- indices[1:(n / 2)]
  test_indices <- indices[(n / 2 + 1):n]
  
  X_train <- X[train_indices, ]
  Z_train <- Z[train_indices, ]
  Y_train <- Y[train_indices]
  
  X_test <- X[test_indices, ]
  Z_test <- Z[test_indices, ]
  Y_test <- Y[test_indices]
  
  # Prepare data for multiview package
  x_train <- list(X_train, Z_train)
  x_test <- list(X_test, Z_test)
  
  # Perform cross-validation to find the optimal lambda
  cvfit <- cv.multiview(x_train, Y_train, family = binomial(), type.measure = "class")
  best_lambda <- cvfit$lambda.min
  lambda_1se <- cvfit$lambda.1se
  selected_lambda <- best_lambda
  # Plot cross-validation results
  png(filename = paste0(output_dir, "/", subtype, "_CVPlot.png"), width = 800, height = 600)
  plot(cvfit)
  abline(v = log(best_lambda), col = "blue", lty = 2)
  abline(v = log(lambda_1se), col = "red", lty = 2)
  dev.off()
  
  # Refit the model with the optimal lambda
  fit_optimized <- multiview(x_train, Y_train, family = binomial(), lambda = selected_lambda)
  
  # Make predictions on the test set
  test_probabilities <- predict(fit_optimized, newx = x_test, s = selected_lambda, type = "response")
  test_predictions <- as.numeric(test_probabilities > 0.5)
  
  # Calculate test accuracy
  test_accuracy <- mean(test_predictions == Y_test)
  cat("Test Accuracy for", subtype, "vs CTRL:", test_accuracy, "\n")
  # Collect model details
  model_details <- data.frame(
    Class = subtype,
    BestLambda = best_lambda,
    Lambda1SE = lambda_1se,
    TestAccuracy = test_accuracy
  )
  
  # Save or append results to the CSV file
  if (!file.exists(results_file)) {
    write.csv(model_details, results_file, row.names = FALSE)  # Create a new file
  } else {
    write.table(model_details, results_file, row.names = FALSE, col.names = FALSE, sep = ",", append = TRUE)  # Append to the file
  }
  # Save confusion matrix
  confusion_data <- as.data.frame(table(Predicted = test_predictions, Actual = Y_test))
  confusion_plot <- ggplot(confusion_data, aes(x = Actual, y = Predicted)) +
    geom_tile(aes(fill = Freq), color = "white") +
    geom_text(aes(label = Freq), size = 5) +
    scale_fill_gradient(low = "white", high = "blue") +
    labs(title = paste("Confusion Matrix for", subtype, "vs CTRL"), x = "Actual", y = "Predicted") +
    theme_minimal()
  
  ggsave(
    filename = paste0(output_dir, "/", subtype, "_ConfusionMatrix.png"),
    plot = confusion_plot, dpi = 300, width = 8, height = 6
  )
  
  # Analyze Feature Importance
  coef <- coef(fit_optimized, s = selected_lambda)
  coef_matrix <- as.matrix(coef)[-1, , drop = FALSE]  # Remove intercept
  feature_names <- rownames(coef_matrix)
  coefficients <- as.numeric(coef_matrix)
  coef_df <- data.frame(Feature = feature_names, Coefficient = coefficients)
  
  # Get top 20 features
  top_features <- coef_df %>%
    filter(Coefficient != 0) %>%
    arrange(desc(abs(Coefficient))) %>%
    head(20)
  
  # Save top features
  write.csv(
    top_features,
    file = paste0(output_dir, "/", subtype, "_TopFeatures.csv"),
    row.names = FALSE
  )
  
  # Create a bar plot of the top 20 features
  plt <- ggplot(top_features, aes(x = reorder(Feature, Coefficient), y = Coefficient)) +
    geom_bar(stat = "identity", fill = "steelblue") +
    coord_flip() +
    labs(title = paste0(subtype, " vs CTRL: Top Features by Coefficient"), x = "Feature", y = "Coefficient") +
    theme_minimal()
  
  ggsave(
    filename = paste0(output_dir, "/", subtype, "_TopFeatures.png"),
    plot = plt, dpi = 300, width = 8, height = 6
  )
}

# List of subtypes
subtypes <- c("PV", "ET", "MF")

# Output directory
output_dir <- "results/co-op_learning/results_subtype_vs_control_12_17"
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}
results_file <- "results/co-op_learning/results_subtype_vs_control_12_17/model_details.csv"
if (file.exists(results_file)) {
  file.remove(results_file)  # Remove old file if exists
}
data <- preprocess_subtype_vs_control(subtypes[1], control = "CTRL", outcome_data, X, Z)
multiview_subtype_vs_control(subtypes[1], data$X, data$Z, data$Y, output_dir)
# Iterate over subtypes and perform analysis
for (subtype in subtypes) {
  # Preprocess data for the current subtype vs control
  data <- preprocess_subtype_vs_control(subtype, control = "CTRL", outcome_data, X, Z)
  
  # Perform multiview analysis
  multiview_subtype_vs_control(subtype, data$X, data$Z, data$Y, output_dir)
}
