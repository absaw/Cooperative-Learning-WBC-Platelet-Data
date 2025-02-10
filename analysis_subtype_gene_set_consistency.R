# Load necessary libraries
library(dplyr)
library(readr)
library(multiview)
library(ggplot2)

#===============================================================================
# PART 1 :: Load and preprocess data
#===============================================================================
# Load the datasets
platelet_data <- read.csv("data/NormCountsPlateletAdjusted.csv", stringsAsFactors = FALSE)
wbc_data <- read.csv("data/NormCountsWBCAdjusted.csv", stringsAsFactors = FALSE)
outcome_data <- read.csv("data/DeSeq2_AllData_n120_PTPRCadj.csv", stringsAsFactors = FALSE)

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
# PART 2: Subtype vs Control Analysis Functions
#===============================================================================
# Function to preprocess data for subtype vs control classification
preprocess_subtype_vs_control <- function(subtype, control = "CTRL", outcome_data, X, Z) {
  selected_samples <- rownames(outcome_data)[outcome_data[[1]] %in% c(subtype, control)]
  selected_samples <- intersect(selected_samples, common_patients)
  
  X_selected <- X[selected_samples, , drop = FALSE]
  Z_selected <- Z[selected_samples, , drop = FALSE]
  Y_selected <- outcome_data[selected_samples, , drop = FALSE]
  Y_selected <- as.numeric(Y_selected[[1]] == subtype)  # Convert outcome to binary
  
  return(list(X = X_selected, Z = Z_selected, Y = Y_selected))
}

# Function to perform multiview classification for subtype vs control
multiview_subtype_vs_control <- function(subtype, X, Z, Y, output_dir, seed = 123, iterations = 100) {
  set.seed(seed)
  
  # Initialize a list to store top features and model details across iterations
  all_features <- list()
  model_details_list <- list()
  
  for (iteration in seq_len(iterations)) {
    cat("Iteration:", iteration, "for Subtype:", subtype, "\n")
    
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
    
    # Refit the model with the optimal lambda
    fit_optimized <- multiview(x_train, Y_train, family = binomial(), lambda = best_lambda)
    
    # Make predictions on the test set
    test_probabilities <- predict(fit_optimized, newx = x_test, s = best_lambda, type = "response")
    test_predictions <- as.numeric(test_probabilities > 0.5)
    
    # Calculate test accuracy
    test_accuracy <- mean(test_predictions == Y_test)
    cat("Test Accuracy for", subtype, "vs CTRL in Iteration", iteration, ":", test_accuracy, "\n")
    
    # Analyze Feature Importance
    coef <- coef(fit_optimized, s = best_lambda)
    coef_matrix <- as.matrix(coef)[-1, , drop = FALSE]  # Remove intercept
    coef_df <- data.frame(
      Iteration = iteration,  # Add iteration number
      Feature = rownames(coef_matrix),
      Coefficient = as.numeric(coef_matrix)
    ) %>%
      filter(Coefficient != 0) %>%
      arrange(desc(abs(Coefficient)))
    
    # Save top features for this iteration
    iteration_file <- paste0(output_dir, "/", subtype, "_TopFeatures_Iteration_", iteration, ".csv")
    write.csv(coef_df, file = iteration_file, row.names = FALSE)
    
    # Store the features for comparison
    all_features[[iteration]] <- coef_df$Feature
    
    # Save model details for this iteration
    model_details_list[[iteration]] <- data.frame(
      Iteration = iteration,
      Subtype = subtype,
      TestAccuracy = test_accuracy,
      BestLambda = best_lambda
    )
  }
  
  # Find common features across all iterations
  common_features <- Reduce(intersect, all_features)
  
  # Save common features
  common_features_file <- paste0(output_dir, "/", subtype, "_CommonFeatures.csv")
  write.csv(data.frame(Feature = common_features), file = common_features_file, row.names = FALSE)
  cat("Common features saved for Subtype:", subtype, "\n")
  
  # Save model details for all iterations
  model_details_df <- do.call(rbind, model_details_list)
  model_details_file <- paste0(output_dir, "/", subtype, "_ModelDetails.csv")
  write.csv(model_details_df, file = model_details_file, row.names = FALSE)
  
  # Count occurrences of each feature across all iterations
  feature_counts <- unlist(all_features) %>%
    table() %>%
    as.data.frame()
  
  # Rename columns for clarity
  colnames(feature_counts) <- c("Feature", "Occurrences")
  
  # Sort by most frequently occurring features
  feature_counts <- feature_counts %>%
    arrange(desc(Occurrences))
  
  # Save feature counts to a CSV file
  feature_counts_file <- paste0(output_dir, "/", subtype, "_FeatureCounts.csv")
  write.csv(feature_counts, file = feature_counts_file, row.names = FALSE)
  
  cat("Feature counts saved for Subtype:", subtype, "\n")
}
#===============================================================================
# PART 3: Run Analysis
#===============================================================================
# List of subtypes
subtypes <- c("PV", "ET", "MF")

# Output directory
output_dir <- "results/co-op_learning/results_adjusted_subtype_vs_control_consistent"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}
data <- preprocess_subtype_vs_control(subtypes[1], control = "CTRL", outcome_data, X, Z)

# Perform multiview analysis with multiple iterations
multiview_subtype_vs_control(subtypes[1], data$X, data$Z, data$Y, output_dir)
# Perform analysis for each subtype
for (subtype in subtypes) {
  # Preprocess data for the current subtype vs control
  data <- preprocess_subtype_vs_control(subtype, control = "CTRL", outcome_data, X, Z)
  
  # Perform multiview analysis with multiple iterations
  multiview_subtype_vs_control(subtype, data$X, data$Z, data$Y, output_dir)
}
