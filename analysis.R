# author : abhishek Sawalkar
# importing libraries
library(dplyr)
library(readr)
library(multiview)
setwd("~/Library/CloudStorage/Box-Box/KrishnanA-Stats-Share/LucyExtendedTermProjects/Project5-WBC_RNAseq-Analysis/WBC-Platelet-Coop-Learning-AS/Cooperative-Learning-WBC-Platelet-Data")
#===============================================================================
# PART 1 :: Pre-processing data
#===============================================================================
# Load the datasets

platelet_data <- read.csv("data/NormCountsPlatelet.csv")
wbc_data <- read_csv("data/NormCountsWBC.csv")
outcome_data <- read_csv("data/DeSeq2_AllData_n120_PTPRCadj.csv")

# Extract patient IDs
platelet_patients <- platelet_data[[1]]
wbc_patients <- wbc_data[[1]]
outcome_patients <- outcome_data[[1]]

# Find common patients
common_patients <- intersect(intersect(platelet_patients, wbc_patients), outcome_patients)
length(common_patients)

# Subset the data for common patients
# Filter rows for common patients and remove the first column (patient IDs)
#X <- platelet_data %>% filter(.[[1]] %in% common_patients)
X <- platelet_data %>% filter(.[[1]] %in% common_patients)
Z <- wbc_data %>% filter(.[[1]] %in% common_patients) 
Y <- outcome_data %>% filter(.[[1]] %in% common_patients)

# Ensure the order of rows is consistent across datasets
X <- X %>% arrange(.[[1]])
Y <- Y %>% arrange(.[[1]])
Z <- Z %>% arrange(.[[1]])

# to remove first column of sampleids
X <- X %>% select(-1)
Z <- Z %>% select(-1)
MPNDisease<-Y %>% select(2)

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
combined_data <- data.frame(X, Z, MPNDisease = MPNDisease[[1]])

# Set seed for reproducibility
set.seed(123)

# Define split ratios
train_ratio <- 0.8
val_ratio <- 0
test_ratio <- 0.2

# Number of samples
n <- nrow(combined_data)
train_size = floor(train_ratio * n)
val_size= floor(val_ratio*n)
test_size = n-train_size-val_size
# Generate random indices for training, validation, and test sets
# Generate a shuffled sequence of row indices
shuffled_indices <- sample(seq_len(n))

# Assign the first `train_size` indices to the training set
train_indices <- shuffled_indices[1:train_size]

# Assign the next `val_size` indices to the validation set
#val_indices <- shuffled_indices[(train_size + 1):(train_size + val_size)]

# Assign the remaining indices to the test set
test_indices <- shuffled_indices[(train_size + val_size + 1):n]


# Training set
X_train <- X[train_indices, ]
Z_train <- Z[train_indices, ]
Y_train <- MPNDisease[train_indices, ]

# Validation set
#X_val <- X[val_indices, ]
#Z_val <- Z[val_indices, ]
#Y_val <- MPNDisease[val_indices, ]

# Test set
X_test <- X[test_indices, ]
Z_test <- Z[test_indices, ]
Y_test <- MPNDisease[test_indices, ]

# Prepare data for multiview package
# Convert MPNDisease to a factor
Y_train <- as.factor(Y_train[[1]])
#Y_val <- as.factor(Y_val[[1]])
Y_test <- as.factor(Y_test[[1]])

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
# multiview package doesn't support mutlinomial classification
# so for our problem, with 4 different subtypes possible, one v/s all method 
# needs to be used. Hence binomial() family then can be used for binary 
#classifiaction

# Prepare Binary Labels for PV
# Create binary labels: 1 for PV, 0 for others
Y_binary_train <- as.numeric(Y_train == "PV")
#Y_binary_val <- as.numeric(Y_val == "PV")
Y_binary_test <- as.numeric(Y_test == "PV")

# Check the distribution of binary labels
cat("Training set distribution:\n")
table(Y_binary_train)
cat("Validation set distribution:\n")
#table(Y_binary_val)
cat("Test set distribution:\n")
table(Y_binary_test)


# Train the Binary Classifier for PV
# Fit the cooperative learning model
fit_pv <- multiview(x_train, Y_binary_train, family = binomial(), rho = 0.5)
?multiview
plot(fit_pv)

# Evaluate the Model on the Test Set
# Predict probabilities for the test set
test_probabilities <- predict(fit_pv, newx = list(X_test,Z_test),s = c(0.01), type = "response")
test_probabilities
# Convert probabilities to binary predictions (threshold = 0.5)
test_predictions <- as.numeric(test_probabilities > 0.5)

# Calculate test accuracy
test_accuracy <- mean(test_predictions == Y_binary_test)
cat("Test Accuracy for PV:", test_accuracy, "\n")

# Confusion Matrix for Test Set
cat("Confusion Matrix for Test Set:\n")
print(table(Predicted = test_predictions, Actual = Y_binary_test))

# Analyze Feature Importance
# Extract coefficients for X (platelet data) and Z (WBC data)
coef <- coef(fit_pv,s=0.1)
coef
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
library(ggplot2)

# Create a bar plot of the top 10 features
ggplot(top_features, aes(x = reorder(Feature, Coefficient), y = Coefficient)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +
  labs(title = "PV Top 20 Features by Coefficient: Binomial Distribution", x = "Feature", y = "Coefficient") +
  theme_minimal()


# Perform cross-validation to find the optimal lambda
cvfit <- cv.multiview(x_train, Y_binary_train, family = binomial(), type.measure = "class")

# Optimal lambda minimizing cross-validation error
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

# Make predictions on the test set
test_probabilities <- predict(fit_optimized, newx = x_test, s = best_lambda, type = "response")

# Convert probabilities to binary predictions
test_predictions <- as.numeric(test_probabilities > 0.5)

# Calculate test accuracy
test_accuracy <- mean(test_predictions == Y_binary_test)
cat("Test Accuracy with Optimized Lambda:", test_accuracy, "\n")
# Save the Model
# Save the fitted model for later use
#saveRDS(fit_pv, file = "model/pv_cooperative_learning_model.rds")

# To load the model in future:
# fit_pv <- readRDS("pv_cooperative_learning_model.rds")


#Debugging 
cat("Number of rows in X_test:", nrow(X_test), "\n")
cat("Number of rows in Z_test:", nrow(Z_test), "\n")
cat("Length of Y_binary_test:", length(Y_binary_test), "\n")
cat("Length of test_predictions:", length(test_predictions), "\n")
