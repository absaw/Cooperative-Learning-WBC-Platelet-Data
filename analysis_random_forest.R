install.packages("randomForest")
library(randomForest)
library(dplyr)
library(readr)
# library(multiview)
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

test_size = n-train_size
test_size
# Generate random indices for training, validation, and test sets
# Generate a shuffled sequence of row indices
shuffled_indices <- sample(seq_len(n))

# Assign the first `train_size` indices to the training set
train_indices <- shuffled_indices[1:train_size]

# Assign the remaining indices to the test set
test_indices <- shuffled_indices[(train_size + 1):n]
# For stratified sampling
# library(caret)  
# set.seed(123)
# # Stratified sampling to preserve class proportions
# n <- nrow(X)
# MPNDisease <- Y[, 1, drop = TRUE]
# train_indices <- createDataPartition(MPNDisease, p = train_ratio, list = FALSE)
# test_indices <- setdiff(seq_len(n), train_indices)
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
Y_test <- as.factor(Y_test)

# Check class distributions
table(Y_train)
table(Y_test)



# Combine normalized X and Z features into a single data frame
train_data <- data.frame(X_train, Z_train, MPNDisease = Y_train)
test_data <- data.frame(X_test, Z_test, MPNDisease = Y_test)

# Train Random Forest model for multinomial classification
set.seed(123)
rf_model <- randomForest(
  MPNDisease ~ .,            # Formula: MPNDisease as the target, all other columns as predictors
  data = train_data,         # Training data
  ntree = 100,               # Number of trees
  mtry = floor(sqrt(ncol(train_data) - 1)),  # Number of variables randomly sampled at each split
  importance = TRUE          # Calculate variable importance
)

# Print model summary
print(rf_model)

# Evaluate the model on the test set
test_predictions <- predict(rf_model, newdata = test_data)

# Confusion matrix
conf_matrix <- table(Predicted = test_predictions, Actual = test_data$MPNDisease)
cat("Confusion Matrix:\n")
print(conf_matrix)

# Calculate overall accuracy
accuracy <- sum(diag(conf_matrix)) / sum(conf_matrix)
cat("Test Accuracy:", accuracy, "\n")

# Plot feature importance
importance_values <- importance(rf_model)
importance_df <- data.frame(
  Feature = rownames(importance_values),
  MeanDecreaseAccuracy = importance_values[, "MeanDecreaseAccuracy"],
  MeanDecreaseGini = importance_values[, "MeanDecreaseGini"]
)

# Plot MeanDecreaseAccuracy
library(ggplot2)
importance_plot <- ggplot(importance_df, aes(x = reorder(Feature, MeanDecreaseAccuracy), y = MeanDecreaseAccuracy)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +
  labs(title = "Feature Importance (Mean Decrease Accuracy)", x = "Features", y = "Mean Decrease Accuracy") +
  theme_minimal()

# Save and display plot
ggsave("Feature_Importance.png", plot = importance_plot, dpi = 300, width = 8, height = 6)
print(importance_plot)