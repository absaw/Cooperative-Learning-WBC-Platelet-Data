# Load required libraries
library(tidyverse)
library(ggplot2)
library(PRROC)

#===============================================================================
# PART 1: Load the Results Files
#===============================================================================
# Load datasets for predictions
wbc_multi_res <- read.csv("../Results/WBC_MultiLasso_Group_Prediction_OverlapWBC-Platelet_n~72.csv")
plate_multi_res <- read.csv("../Results/Platelet_MultiLasso_Group.csv")
base_res <- read.csv("../Results/BaseMultinomial.csv")
coop_learning_res <- read.csv("results/co-op_learning/results_whole_dataset_prediction_final/MF_Predictions.csv")

#===============================================================================
# PART 2: Prepare ROC Data for Each Model
#===============================================================================
# Function to prepare ROC data
generate_roc_data <- function(fg, bg) {
  roc <- PRROC::roc.curve(scores.class0 = fg, scores.class1 = bg, curve = TRUE)
  data <- data.frame(FP = roc$curve[, 1], TP = roc$curve[, 2])
  list(roc_data = data, auc = roc$auc)
}

# WBC Model
roc_wbc <- generate_roc_data(
  fg = wbc_multi_res$Prob.MF[wbc_multi_res$True.Label == "MF"],
  bg = wbc_multi_res$Prob.MF[wbc_multi_res$True.Label != "MF"]
)

# Platelet Model
roc_plate <- generate_roc_data(
  fg = plate_multi_res$Prob.MF[plate_multi_res$True.Label == "MF"],
  bg = plate_multi_res$Prob.MF[plate_multi_res$True.Label != "MF"]
)

# Baseline Model
roc_base <- generate_roc_data(
  fg = base_res$Prob.MF[base_res$True.Label == "MF"],
  bg = base_res$Prob.MF[base_res$True.Label != "MF"]
)

# Co-op Learning Model
roc_coop <- generate_roc_data(
  fg = coop_learning_res$Predicted_Probability[coop_learning_res$Actual_Label == "MF"],
  bg = coop_learning_res$Predicted_Probability[coop_learning_res$Actual_Label != "MF"]
)

#===============================================================================
# PART 3: Combine and Plot ROC Curves
#===============================================================================
# Combine all ROC data with model labels
roc_data_combined <- bind_rows(
  mutate(roc_wbc$roc_data, Model = "WBC Lasso Logistic"),
  mutate(roc_plate$roc_data, Model = "Platelet Lasso Logistic"),
  mutate(roc_base$roc_data, Model = "Baseline Logistic"),
  mutate(roc_coop$roc_data, Model = "Co-op Learning")
)

# Define colors for each model
colors <- c(
  "WBC Lasso Logistic" = "red",
  "Platelet Lasso Logistic" = "blue",
  "Baseline Logistic" = "black",
  "Co-op Learning" = "green"
)

# Create the ROC plot
p <- ggplot() +
  geom_line(data = roc_data_combined, aes(x = FP, y = TP, color = Model), size = 1) +
  geom_abline(slope = 1, intercept = 0, color = "gray60", linetype = "twodash", size = 0.8) +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.1)) +
  scale_y_continuous(limits = c(0, 1.005), breaks = seq(0, 1, 0.1)) +
  scale_color_manual(name = "", values = colors) +
  labs(
    title = "Combined ROC Curves",
    x = "False Positive Rate",
    y = "True Positive Rate"
  ) +
  theme_bw() +
  theme(
    text = element_text(size = 14, family = "sans"),
    axis.title = element_text(size = 20, face = "bold"),
    axis.text.y = element_text(size = 16),
    axis.text.x = element_text(size = 16),
    legend.position = "bottom"
  ) +
  annotate("text", x = 0.85, y = 0.05, label = paste0(
    "AUC:\n",
    "WBC = ", round(roc_wbc$auc, 3), ", ",
    "Platelet = ", round(roc_plate$auc, 3), ", ",
    "Baseline = ", round(roc_base$auc, 3), ", ",
    "Co-op = ", round(roc_coop$auc, 3)
  ), size = 4, color = "black", hjust = 1)

# Display the plot
print(p)

# Save the plot as a PNG file
output_dir <- "results/combined_roc"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
ggsave(
  filename = file.path(output_dir, "Combined_ROC_Curve.png"),
  plot = p, width = 8, height = 6, dpi = 300
)
