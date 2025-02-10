# Load required libraries
library(tidyverse)
library(ggplot2)
library(PRROC)

#===============================================================================
# PART 1: Load the Results File
#===============================================================================
# Read the results file for predictions
setwd("~/Library/CloudStorage/Box-Box/KrishnanA-Stats-Share/LucyExtendedTermProjects/Project5-WBC_RNAseq-Analysis/WBC-Platelet-Coop-Learning-AS/Cooperative-Learning-WBC-Platelet-Data")

results <- read.csv("results/co-op_learning/results_adjusted_mf/MF_Predictions.csv")

# Verify the structure of the results file
head(results)

# Columns in the file:
# 1. Sample_ID: Unique patient IDs
# 2. Predicted_Probability: Probability of being MF
# 3. Predicted_Label: Predicted binary label (MF/Not MF)
# 4. Actual_Label: Actual binary label (MF/Not MF)

#===============================================================================
# PART 2: Prepare Data for ROC Curve
#===============================================================================
# Extract prediction probabilities and actual labels
fg <- results$s1[results$Actual_Label == "MF"]  # Positive class (MF)
bg <- results$s1[results$Actual_Label != "MF"]  # Negative class (Not MF)

#===============================================================================
# PART 3: Generate the ROC Curve
#===============================================================================
# Calculate the ROC curve
roc_results <- PRROC::roc.curve(scores.class0 = fg, scores.class1 = bg, curve = TRUE)

# Print the AUROC score
cat("AUROC Score for MF vs Not MF:", roc_results$auc, "\n")

#===============================================================================
# PART 4: Save ROC Data for Plotting
#===============================================================================
# # Convert ROC curve to a data frame for custom plotting
# roc_df <- data.frame(FPR = roc_results$curve[,1],  # False Positive Rate
#                      TPR = roc_results$curve[,2]) # True Positive Rate
#===============================================================================
# PART 5: Plot the ROC Curve for MF Only
#===============================================================================
# Save the ROC curve data to a data frame
roc_mf_tab <- data.frame(FP = roc_results$curve[, 1],  # False Positive Rate
                         TP = roc_results$curve[, 2]) # True Positive Rate

# Define colors for the MF ROC curve
colors <- c("Co-op Learning" = "red")

# Create the ROC plot for MF
p <- ggplot() + 
  # Plot the MF ROC curve
  geom_point(data = roc_mf_tab, aes(x = FP, y = TP, col = "Co-op Learning"), size = 0.2) + 
  geom_line(data = roc_mf_tab, aes(x = FP, y = TP, col = "Co-op Learning"), size = 1, lty = 1) + 
  # Add diagonal reference line
  geom_abline(slope = 1, intercept = 0, col = "gray60", linetype = "twodash", size = 0.8) +
  # Format axes
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.1)) + 
  scale_y_continuous(limits = c(0, 1.005), breaks = seq(0, 1, 0.1)) + 
  # Add color legend
  scale_color_manual(name = "", values = colors) + 
  # Labels and theme
  labs(title = "ROC Curve for MF Classification",
       x = "False Positive Rate", 
       y = "True Positive Rate") + 
  theme_bw() + 
  theme(text = element_text(size = 14, family = "sans"), 
        axis.title = element_text(size = 20, face = "bold"), 
        axis.text.y = element_text(size = 16), 
        axis.text.x = element_text(size = 16), 
        legend.position = "bottom")
# Add AUC annotation to the plot
p <- p + annotate("text", x = 0.85, y = 0.05, 
                  label = paste("AUC =", round(roc_results$auc, 3)), 
                  size = 6, color = "black", fontface = "bold")
# Display the plot
print(p)

output_dir <- "results/co-op_learning/results_adjusted_mf"
# Save the ROC plot as a PNG file
ggsave(paste0(output_dir,"/MF_ROC_Curve.png"), plot = p, width = 8, height = 6, dpi = 300)
