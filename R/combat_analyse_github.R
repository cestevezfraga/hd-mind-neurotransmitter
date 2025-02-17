#Script to review the results of COMBAT

library(effectsize)
library(tidyverse)
library(ggm)
library(ggplot2)
library(tidyr)
library(dplyr)
library(readr)
library(lsr)
library(ggplot2)
library(tidyr)
library(dplyr)
library (data.table)
library(readr)
library(reshape2)

rm(list = ls())

#Get the combat / non combat (raw) node degree data
combat<-read.csv("/Users/charlie/Desktop/my_projects/neurotransmitter/github/data/data_combat.csv",  header = FALSE)

raw<-read.csv("/Users/charlie/Desktop/my_projects/neurotransmitter/github/data/mind_node.csv", header = FALSE)

t_combat<- transpose(combat)
t_raw<-transpose(raw)

ct_data <- read.csv("/Users/charlie/Desktop/my_projects/neurotransmitter/github/data/track_fsv7.csv")
colnames <-colnames(ct_data[8:75])
colnames(t_combat)<-colnames
colnames(t_raw)<-colnames

demos<-read.csv("/Users/charlie/Desktop/my_projects/neurotransmitter/github/data/demos_long.csv")

merged_combat<-cbind(demos,t_combat)
merged_raw<-cbind(demos,t_raw)


t_test_results <- list()
# Initialize an empty vector to store adjusted p-values
adjusted_p_values <- numeric()

for (i in 9:76) {
  # Extract cortical thickness from the ith ROI in merged_raw
  mind_raw <- merged_raw[, i]
  
  # Extract cortical thickness from the ith ROI in merged_combat
  mind_combat <- merged_combat[, i]
  
  # Perform t-test
  t_test_result <- t.test(mind_raw, mind_combat)
  
  # Save t-test result in the list
  t_test_results[[i]] <- t_test_result
  
  # Store p-value for FDR correction
  adjusted_p_values <- c(adjusted_p_values, t_test_result$p.value)
}




uncorrected_p_values<-adjusted_p_values
# Perform FDR correction
adjusted_p_values <- p.adjust(adjusted_p_values, method = "fdr")

# Open a file for writing
file <- file("/Users/charlie/Desktop/my_projects/neurotransmitter/github/results/all_combat_vs_raw.txt", "w")

# Write ROI names, uncorrected p-values, and adjusted p-values to the file
for (i in 9:76) {
  # Get ROI name
  roi_name <- colnames(merged_raw)[i]
  
  # Write ROI name, uncorrected p-value, and adjusted p-value to the file
  writeLines(paste(roi_name, ": Uncorrected p-value =", t_test_results[[i]]$p.value, ", Adjusted p-value =", adjusted_p_values[i - 8]), file)
}

# Close the file
close(file)


results_df <- data.frame(
  ROI = colnames(merged_raw)[9:76], # ROI names
  Uncorrected_p_value = uncorrected_p_values, # Uncorrected p-values
  Adjusted_p_value = adjusted_p_values # Adjusted p-values
)

# Save the results as a table
write.table(results_df, file = "/Users/charlie/Desktop/my_projects/neurotransmitter/github/results/all_combat_vs_raw.csv", sep = ",", row.names = FALSE, quote = FALSE)





######### FIGURES #########

combat<-read.csv("/Users/charlie/Desktop/my_projects/neurotransmitter/github/data_combat.csv",  header = FALSE)

raw<-read.csv("/Users/charlie/Desktop/my_projects/neurotransmitter/github/results/mind_node.csv", header = FALSE)



t_combat<- transpose(combat)
t_raw<-transpose(raw)
num_boxplots <- ncol(t_combat)
x_positions <- seq(1, by = 2.1, length.out = num_boxplots)

ct_data <- read.csv("/Users/charlie/Desktop/my_projects/neurotransmitter/github/data/track_fsv7.csv")
colnames <-colnames(ct_data[8:75])
colnames(t_combat)<-colnames
colnames(t_raw)<-colnames

col_names <- colnames(t_combat)
col_names <- gsub("_thickness", "", col_names)



# Set up plotting area
pdf_file <- "/Users/charlie/Desktop/my_projects/neurotransmitter/github/results/boxplot_combat_vs_raw_all.pdf"
pdf(pdf_file, width = 20, height = 8)

par(mar = c(6, 3, 3, 3) + 0.1)

color_combat <- rgb(0, 0, 1, alpha = 0.4)  # Transparent blue for t_combat
color_raw <- rgb(0, 1, 0, alpha = 0.4)   

# Plot t_combat
boxplot(t_combat, 
        at = x_positions - 1,
        horizontal = FALSE,  # Horizontal orientation
        col = color_combat,   # Color of the boxes for t_combat
        border = "darkblue",    # Border color
        notch = FALSE,       # Disable notches
        outline = FALSE,     # Show outliers
        las = 2,             # Rotate axis labels
        names = col_names,  # Use column names as labels
        cex.axis = 0.5,      # Adjust the size of axis labels
        main = "All HD" # Main title
)

# Overlay t_raw boxplots
boxplot(t_raw, 
        at = x_positions,
        horizontal = FALSE,  # Horizontal orientation
        col = color_raw,   # Color of the boxes for t_raw
        border = "darkgreen",    # Border color
        notch = FALSE,       # Disable notches
        outline = FALSE,     # Show outliers
        las = 2,             # Rotate axis labels
        names = col_names,  # Use column names as labels
        cex.axis = 0.5,      # Adjust the size of axis labels
        add = TRUE          # Add the boxplot to the existing plot
)

# Add legend
legend("topright", legend = c("Combat", "Raw"), fill = c(color_combat, color_raw))

# Close the PDF device
dev.off()





