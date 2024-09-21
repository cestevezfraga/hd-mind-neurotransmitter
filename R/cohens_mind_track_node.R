
#Script to obtain Cohen's D from a MIND similarity matrix
#Two methods:
##Node degree: Averaging ROI data followed by Cohen's D
##Edge weight: Estimate differences between patients and controls first 


library(tidyverse)
library(ggm)
library(ggplot2)
library(tidyr)
library(dplyr)
library(readr)
library(ppcor)

rm(list = ls())

#####Track#####

##Obtaining degree
data_directory <- "/Users/charlie/Desktop/my_projects/neurotransmitter/imaging_metrics/thickness/track/carlos"


# List of subject names

subject_ids <- c(
'XXX-XXX-XXX', 'YYY-YYY-YYY'
)


# Initialize an empty data frame to store the merged data
merged_data <- data.frame()

# Loop through each subject and process their data
for (subject_id in subject_ids) {
  # Construct the file path
  file_path <- file.path(data_directory, subject_id, "mri", "mind.csv")
  
  # Check if the file exists
  if (file.exists(file_path)) {
    # Read the CSV file
    subject_data <- read.csv(file_path)
    subject_data <- subject_data[-1]  # Remove the first column
    
    # Calculate the column-wise average
    subject_avg <- data.frame(matrix(0, ncol = ncol(subject_data), nrow = 1))
    colnames(subject_avg) <- colnames(subject_data)
    subject_avg[1, ] <- colMeans(subject_data)
    #subject_avg <- colMeans(subject_data, na.rm = TRUE)
    
    
    # Add the subject ID as a column
    subject_avg <- data.frame(ID = subject_id, subject_avg)
    
    # Append the subject's averaged data to the merged_data data frame
    merged_data <- rbind(merged_data, subject_avg)
  } else {
    cat("File not found for subject:", subject_id, "\n")
  }
}


# Reset column names if needed (assuming the columns in all files are the same)
colnames(merged_data) <- c("ID", colnames(merged_data)[-1])

ct_data <- read.csv("/Users/charlie/Desktop/my_projects/neurotransmitter/imaging_metrics/thickness/track/reprocessed_results/track_fsv7.csv")


merged_data <- merge(merged_data, ct_data[, c("ID", "site", "sex", "cag", "group", "age1", "burden1", "eTIV")], by = "ID", all.x = TRUE)

desired_order <- c("ID", "site", "sex", "cag", "group", "age1", "burden1", "eTIV", names(merged_data)[!names(merged_data) %in% c("ID", "site", "sex", "cag", "group", "age1", "burden1", "eTIV")])

# Reorder the columns in merged_data
merged_data <- merged_data[, desired_order]



patients <- merged_data %>%filter(group == 2)
controls <- merged_data %>%filter(group == 0)

colnames <-colnames(merged_data[9:76])
k=9


test=0.12
df<-data.frame(test)
write_csv(df,'/Users/charlie/Desktop/my_projects/neurotransmitter/R/mind/track/results/correlate/correlation.csv')

test=0.12
df<-data.frame(test)
write_csv(df,'/Users/charlie/Desktop/my_projects/neurotransmitter/R/mind/track/results/correlate/p_correlation.csv')



for (i in colnames) {
  setwd('/Users/charlie/Desktop/my_projects/neurotransmitter/R/mind/track/results/correlate/')
  
  # Access the column in the dataframe using the column name
  print(i)
  roi<-patients[, i]
  
  m<-pcor.test(roi, patients$burden1, patients[,c("site","sex","eTIV","age1")])
  

  coefficient<-m$estimate
  pval<-m$p

  z = data.frame(coefficient = coefficient, pval=pval)
  rownames(z)<-c(print(i))
  write.csv(z, paste0(i))
  
  df_d<-read.csv('/Users/charlie/Desktop/my_projects/neurotransmitter/R/mind/track/results/correlate/correlation.csv')
  df_d2<-data.frame(coefficient)
  colnames(df_d2)[1] = c(print(i))
  df_d3<-cbind(df_d, df_d2)
  write.csv (df_d3, '/Users/charlie/Desktop/my_projects/neurotransmitter/R/mind/track/results/correlate/correlation.csv')
  
  
  df_e<-read.csv('/Users/charlie/Desktop/my_projects/neurotransmitter/R/mind/track/results/correlate/p_correlation.csv')
  df_e2<-data.frame(pval)
  colnames(df_d2)[1] = c(print(i))
  df_e3<-cbind(df_e, df_e2)
  write.csv (df_e3, '/Users/charlie/Desktop/my_projects/neurotransmitter/R/mind/track/results/correlate/p_correlation.csv')
  
  
  p <- ggplot(patients, aes_string(x = roi, y = patients$burden1)) +
    geom_point() +
    geom_smooth(method = "lm", se = FALSE, color = 'red') +  # Add linear regression line
    labs(title = paste("Correlation MIND - DBS", i, "\nCoefficient:", round(m$estimate, 2), "\nP-value:", round(m$p, 4)))
  
  ggsave(paste0("/Users/charlie/Desktop/my_projects/neurotransmitter/R/mind/track/results/correlate/","dbs", i, "_scatterplot.pdf"), p, width = 8, height = 6, dpi = 300)
  
  
  
  if (k<75) {k = k+1} else print('analysis done') 
  
  
}




#Remove columns that start with "X" or have non cortical items
h<- read.csv("/Users/charlie/Desktop/my_projects/neurotransmitter/R/mind/track/results/correlate/correlation.csv")
b<- h %>% dplyr::select(-(contains("X") | contains('Mean')|contains('Caudate')|contains('Putamen')|contains('cortex')|contains('icv')|contains('aparc')|contains('test')))
b_trans<-t(b)

yy<- read.csv("/Users/charlie/Desktop/my_projects/neurotransmitter/R/mind/track/results/correlate/p_correlation.csv")
yyy<- yy %>% dplyr::select(-(contains("X") | contains('Mean')|contains('Caudate')|contains('Putamen')|contains('cortex')|contains('icv')|contains('aparc')|contains('test')))
yyy_trans<-data.frame(t(yyy))
yyy_trans$P_FDR <- p.adjust(yyy_trans$t.yyy., method = "BH")
yyy_trans <- yyy_trans %>% rename(Pval = t.yyy.)

correlate<-cbind(b_trans,yyy_trans)
correlate<-correlate %>% rename (coefficient = b_trans)
write.csv(correlate, '/Users/charlie/Desktop/my_projects/neurotransmitter/R/mind/track/results/correlate/correlation_transp.csv')


hhh<-data.frame(b_trans)
gg<-ggplot(hhh, aes(x= hhh$b_trans)) +
  geom_density()
ggsave(gg, filename = 'coefficients.pdf', device = 'pdf', width = 11.69, height = 8.27, dpi = 300, path ='/Users/charlie/Desktop/my_projects/neurotransmitter/R/mind/track/results/correlate/')


hhh<-data.frame(yyy_trans$Pval)
gg<-ggplot(hhh, aes(x= hhh$yyy_trans)) +
  geom_density()
ggsave(gg, filename = 'Pvals.pdf', device = 'pdf', width = 11.69, height = 8.27, dpi = 300, path ='/Users/charlie/Desktop/my_projects/neurotransmitter/R/mind/track/results/correlate/')



hhh<-data.frame(yyy_trans$P_FDR)
gg<-ggplot(hhh, aes(x= hhh$yyy_trans)) +
  geom_density()
ggsave(gg, filename = 'FDR_Pvals.pdf', device = 'pdf', width = 11.69, height = 8.27, dpi = 300, path ='/Users/charlie/Desktop/my_projects/neurotransmitter/R/mind/track/results/correlate/')











######PLASMA NFL######

rm(list = ls())


data_directory <- "/Users/charlie/Desktop/my_projects/neurotransmitter/imaging_metrics/thickness/track/carlos"

subject_ids <- c(
  'XXX-XXX-XXX', 'YYY-YYY-YYY'
)

# Initialize an empty data frame to store the merged data
merged_data <- data.frame()


# Initialize an empty data frame to store the merged data
merged_data <- data.frame()

# Loop through each subject and process their data
for (subject_id in subject_ids) {
  # Construct the file path
  file_path <- file.path(data_directory, subject_id, "mri", "mind.csv")
  
  # Check if the file exists
  if (file.exists(file_path)) {
    # Read the CSV file
    subject_data <- read.csv(file_path)
    subject_data <- subject_data[-1]  # Remove the first column
    
    # Calculate the column-wise average
    subject_avg <- data.frame(matrix(0, ncol = ncol(subject_data), nrow = 1))
    colnames(subject_avg) <- colnames(subject_data)
    subject_avg[1, ] <- colMeans(subject_data)
    #subject_avg <- colMeans(subject_data, na.rm = TRUE)
    
    
    # Add the subject ID as a column
    subject_avg <- data.frame(ID = subject_id, subject_avg)
    
    # Append the subject's averaged data to the merged_data data frame
    merged_data <- rbind(merged_data, subject_avg)
  } else {
    cat("File not found for subject:", subject_id, "\n")
  }
}


# Reset column names if needed (assuming the columns in all files are the same)
colnames(merged_data) <- c("ID", colnames(merged_data)[-1])

ct_data <- read.csv("/Users/charlie/Desktop/my_projects/neurotransmitter/imaging_metrics/thickness/track/reprocessed_results/track_fsv7.csv")




merged_data <- merge(merged_data, ct_data[, c("ID", "site", "sex", "cag", "group", "age1", "burden1", "eTIV", "plasma_nfl","lognfl" )], by = "ID", all.x = TRUE)

desired_order <- c("ID", "site", "sex", "cag", "group", "age1", "burden1", "eTIV", "plasma_nfl","lognfl", names(merged_data)[!names(merged_data) %in% c("ID", "site", "sex", "cag", "group", "age1", "burden1", "eTIV","plasma_nfl","lognfl")])

# Reorder the columns in merged_data
merged_data <- merged_data[, desired_order]



patients <- merged_data %>%filter(group == 2)
patients <- na.omit(patients)
controls <- merged_data %>%filter(group == 0)

colnames <-colnames(merged_data[11:78])
k=9




test=0.12
df<-data.frame(test)
write_csv(df,'/Users/charlie/Desktop/my_projects/neurotransmitter/R/mind/track/results/correlate/correlation_nfl_plasma.csv')

test=0.12
df<-data.frame(test)
write_csv(df,'/Users/charlie/Desktop/my_projects/neurotransmitter/R/mind/track/results/correlate/p_correlation_nfl_plasma.csv')



for (i in colnames) {
  setwd('/Users/charlie/Desktop/my_projects/neurotransmitter/R/mind/track/results/correlate/')
  
  # Access the column in the dataframe using the column name
  print(i)
  roi<-patients[, i]
  
  m<-pcor.test(roi, patients$plasma_nfl, patients[,c("site","sex","eTIV","age1")])
  
  
  coefficient<-m$estimate
  pval<-m$p
  
  z = data.frame(coefficient = coefficient, pval=pval)
  rownames(z)<-c(print(i))
  write.csv(z, paste0(i))
  
  df_d<-read.csv('/Users/charlie/Desktop/my_projects/neurotransmitter/R/mind/track/results/correlate/correlation_nfl_plasma.csv')
  df_d2<-data.frame(coefficient)
  colnames(df_d2)[1] = c(print(i))
  df_d3<-cbind(df_d, df_d2)
  write.csv (df_d3, '/Users/charlie/Desktop/my_projects/neurotransmitter/R/mind/track/results/correlate/correlation_nfl_plasma.csv')
  
  
  df_e<-read.csv('/Users/charlie/Desktop/my_projects/neurotransmitter/R/mind/track/results/correlate/p_correlation_nfl_plasma.csv')
  df_e2<-data.frame(pval)
  colnames(df_d2)[1] = c(print(i))
  df_e3<-cbind(df_e, df_e2)
  write.csv (df_e3, '/Users/charlie/Desktop/my_projects/neurotransmitter/R/mind/track/results/correlate/p_correlation_nfl_plasma.csv')
  
  
  p <- ggplot(patients, aes_string(x = roi, y = patients$plasma_nfl)) +
    geom_point() +
    geom_smooth(method = "lm", se = FALSE, color = 'red') +  # Add linear regression line
    labs(title = paste("Correlation MIND - plasma NFL", i, "\nCoefficient:", round(m$estimate, 2), "\nP-value:", round(m$p, 4)))
  
  ggsave(paste0("/Users/charlie/Desktop/my_projects/neurotransmitter/R/mind/track/results/correlate/","nfl_plasma", i, "_scatterplot.pdf"), p, width = 8, height = 6, dpi = 300)
  
  
  
  if (k<75) {k = k+1} else print('analysis done') 
  
  
}


patients_long <- gather(patients, ROI, mind, starts_with("lh_"), starts_with("rh_"))

p <- ggplot(patients_long, aes(x = mind, y = plasma_nfl)) +
  geom_point(alpha = 0.3) +  # Add points
  geom_smooth(method = "lm", se = FALSE, color = "red") +  # Add a linear regression line
  labs(title = "Correlation between nodal degree and plasma NFL",
       x = "Nodal degree",
       y = "Plasma NFL ") +
  theme_minimal()

# Print the plot

ggsave("/Users/charlie/Desktop/my_projects/neurotransmitter/R/mind/track/results/correlate/ALL_ROI_nfl_plasma_scatterplot.pdf", p, width = 8, height = 6, dpi = 300)




#Remove columns that start with "X" or have non cortical items
h<- read.csv("/Users/charlie/Desktop/my_projects/neurotransmitter/R/mind/track/results/correlate/correlation_nfl_plasma.csv")
b<- h %>% dplyr::select(-(contains("X") | contains('Mean')|contains('Caudate')|contains('Putamen')|contains('cortex')|contains('icv')|contains('aparc')|contains('test')))
b_trans<-t(b)

yy<- read.csv("/Users/charlie/Desktop/my_projects/neurotransmitter/R/mind/track/results/correlate/p_correlation_nfl_plasma.csv")
yyy<- yy %>% dplyr::select(-(contains("X") | contains('Mean')|contains('Caudate')|contains('Putamen')|contains('cortex')|contains('icv')|contains('aparc')|contains('test')))
yyy_trans<-data.frame(t(yyy))
yyy_trans$P_FDR <- p.adjust(yyy_trans$t.yyy., method = "BH")
yyy_trans <- yyy_trans %>% rename(Pval = t.yyy.)

correlate<-cbind(b_trans,yyy_trans)
correlate<-correlate %>% rename (coefficient = b_trans)
write.csv(correlate, '/Users/charlie/Desktop/my_projects/neurotransmitter/R/mind/track/results/correlate/correlation_transp_nfl_plasma.csv')


hhh<-data.frame(b_trans)
gg<-ggplot(hhh, aes(x= hhh$b_trans)) +
  geom_density()
ggsave(gg, filename = 'coefficients_nfl_plasma.pdf', device = 'pdf', width = 11.69, height = 8.27, dpi = 300, path ='/Users/charlie/Desktop/my_projects/neurotransmitter/R/mind/track/results/correlate/')


hhh<-data.frame(yyy_trans$Pval)
gg<-ggplot(hhh, aes(x= hhh$yyy_trans)) +
  geom_density()
ggsave(gg, filename = 'Pvals_nfl_plasma.pdf', device = 'pdf', width = 11.69, height = 8.27, dpi = 300, path ='/Users/charlie/Desktop/my_projects/neurotransmitter/R/mind/track/results/correlate/')



hhh<-data.frame(yyy_trans$P_FDR)
gg<-ggplot(hhh, aes(x= hhh$yyy_trans)) +
  geom_density()
ggsave(gg, filename = 'FDR_Pvals_nfl_plasma.pdf', device = 'pdf', width = 11.69, height = 8.27, dpi = 300, path ='/Users/charlie/Desktop/my_projects/neurotransmitter/R/mind/track/results/correlate/')








