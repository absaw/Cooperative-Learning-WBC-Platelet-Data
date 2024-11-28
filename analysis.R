# author : abhishek Sawalkar
# importing libraries
library(dplyr)
library(readr)

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
train_ratio <- 0.7
val_ratio <- 0.15
test_ratio <- 0.15

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
val_indices <- shuffled_indices[(train_size + 1):(train_size + val_size)]

# Assign the remaining indices to the test set
test_indices <- shuffled_indices[(train_size + val_size + 1):n]


# Training set
X_train <- X[train_indices, ]
Z_train <- Z[train_indices, ]
Y_train <- MPNDisease[train_indices, ]

# Validation set
X_val <- X[val_indices, ]
Z_val <- Z[val_indices, ]
Y_val <- MPNDisease[val_indices, ]

# Test set
X_test <- X[test_indices, ]
Z_test <- Z[test_indices, ]
Y_test <- MPNDisease[test_indices, ]

# Prepare data for multiview package
# Convert MPNDisease to a factor
Y_train <- as.factor(Y_train[[1]])
Y_val <- as.factor(Y_val[[1]])
Y_test <- as.factor(Y_test[[1]])

# Create predictor lists for the multiview package
x_train <- list(X_train, Z_train)
x_val <- list(X_val, Z_val)
x_test <- list(X_test, Z_test)

# Check class distributions
table(Y_train)
table(Y_val)
table(Y_test)
