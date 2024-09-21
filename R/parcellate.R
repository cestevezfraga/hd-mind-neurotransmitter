
#Script to compare different parcellations
#Using MIND similarity matrices
#For mHD from TrackHD

##Node degree: Averaging ROI data followed by Cohen's D
##Edge weight: Estimate differences between patients and controls first 

library(effectsize)
library(tidyverse)
library(ggm)
library(ggplot2)
library(tidyr)
library(dplyr)
library(readr)
library(lsr)

rm(list = ls())

#Obtaining degree
data_directory <- "/Users/charlie/Desktop/my_projects/neurotransmitter/imaging_metrics/thickness/track/carlos"
ct_data <- read.csv("/Users/charlie/Desktop/my_projects/neurotransmitter/imaging_metrics/thickness/track/reprocessed_results/track_fsv7.csv")


# List of subject names

subject_ids <- c(
  'XXX-XXX-XXX', 'YYY-YYY-YYY'
)



###### DK-68 ######

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

merged_data_DK68<-merged_data
write_csv(merged_data_DK68,'/Users/charlie/Desktop/my_projects/neurotransmitter/R/parcellate/results/merged_data_DK68.csv')

avg_rois <- rowMeans(merged_data[, -1])
avg_sample_DK68<-mean(avg_rois)

colnames(merged_data)[colnames(merged_data) == "ID"] <- "ID"
merged_data <- merge(merged_data, ct_data[, c("ID", "site", "sex", "cag", "group", "age1", "burden1", "eTIV")], by = "ID", all.x = TRUE)
desired_order <- c("ID", "site", "sex", "cag", "group", "age1", "burden1", "eTIV", names(merged_data)[!names(merged_data) %in% c("ID", "site", "sex", "cag", "group", "age1", "burden1", "eTIV")])
merged_data <- merged_data[, desired_order]

##Estimating Cohen's D##
colnames <-colnames(merged_data[9:76])
k=9

test=0.12
df<-data.frame(test)
write_csv(df,'/Users/charlie/Desktop/my_projects/neurotransmitter/R/parcellate/results/cohens_d_DK68.csv')

for (i in colnames) {
  setwd('/Users/charlie/Desktop/my_projects/neurotransmitter/R/parcellate/results/')
  
  # Access the column in the dataframe using the column name
  print(i)
  roi<-merged_data[, i]
  
  # Linear mixed effects model that describes cortical thickness
  # The T value here is the adjusted T 
  model <- lm( roi ~ merged_data$group + merged_data$age1 + merged_data$sex + merged_data$site + merged_data$eTIV)
  df.residual(model)
  m<-summary(model)
  print(m)
  
  
  estimate<-m$coefficients[ , 1]
  se<-m$coefficients[ , 2]
  adjusted_tval<-m$coefficients[ , 3]
  p<-m$coefficients[ , 4]
  
  #Get the adjusted T value for the ROI
  adjusted_tval_cxthick<-adjusted_tval[2]
  print(adjusted_tval_cxthick)
  
  
  # To estimate  the partial D value (adjusted for covariantes)
  # We will need the adjusted T from  lm()
  adjusted_d<-t_to_d(adjusted_tval_cxthick, df_error=df.residual(model))[[1]]
  adjusted_d_ci_low<-t_to_d(adjusted_tval_cxthick, df_error=df.residual(model))[[3]]
  adjusted_d_ci_high<-t_to_d(adjusted_tval_cxthick, df_error=df.residual(model))[[4]]
  
  #Now to estimate unadjusted T and D
  # First get the data in patients and controls
  x<- merged_data  %>% dplyr::select(k)
  patients<-subset(x, merged_data$group == "2")
  controls<-subset(x, merged_data$group == "0")
  
  patients_mean<-mean(patients[[1]])
  controls_mean<-mean(controls[[1]])
  
  s1<-sd(patients[[1]])
  s2<-sd(controls[[1]])
  n1<-length(patients[[1]])
  n2<-length(controls[[1]])
  pooled <- sqrt(((n1-1)*s1^2 + (n2-1)*s2^2) / (n1+n1-2))
  
  #Plain T value
  unadj_Tval <- (patients_mean - controls_mean) / pooled
  
  #Plain Cohen's D
  unadj_d<-cohensD(x = patients[[1]], y = controls[[1]])
  
  z = data.frame(D_adjusted = adjusted_d, adjusted_d_ci_low = adjusted_d_ci_low, 
                 adjusted_d_ci_high, D_unadjusted = unadj_d, 
                 T_adjusted = adjusted_tval_cxthick,T_unadjusted = unadj_Tval )
  rownames(z)<-c(print(i))
  write.csv(z, paste0(i))
  
  df_d<-read.csv('/Users/charlie/Desktop/my_projects/neurotransmitter/R/parcellate/results/cohens_d_DK68.csv')
  df_d2<-data.frame(adjusted_d)
  colnames(df_d2)[1] = c(print(i))
  df_d3<-cbind(df_d, df_d2)
  write.csv (df_d3, '/Users/charlie/Desktop/my_projects/neurotransmitter/R/parcellate/results/cohens_d_DK68.csv')
  
  if (k<75) {k = k+1} else print('analysis done') 
  
}
#Remove columns that start with "X" or have non cortical items
h<- read.csv("/Users/charlie/Desktop/my_projects/neurotransmitter/R/parcellate/results/cohens_d_DK68.csv")
b<- h %>% dplyr::select(-(contains("X") | contains('Mean')|contains('Caudate')|contains('Putamen')|contains('cortex')|contains('icv')|contains('aparc')|contains('test')))
b_trans<-t(b)
write.csv(b_trans, '/Users/charlie/Desktop/my_projects/neurotransmitter/R/parcellate/results/cohens_d_transp_DK68.csv')

















###### DK-318 #####

# Initialize an empty data frame to store the merged data
merged_data <- data.frame()

# Loop through each subject and process their data
for (subject_id in subject_ids) {
  # Construct the file path
  file_path <- file.path(data_directory, subject_id, "mri", "mind_500.sym.aparc.csv")
  
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


merged_data_DK318<-merged_data
write_csv(merged_data_DK318,'/Users/charlie/Desktop/my_projects/neurotransmitter/R/parcellate/results/merged_data_DK318.csv')

avg_rois <- rowMeans(merged_data[, -1])
avg_sample_DK318<-mean(avg_rois)


colnames(merged_data)[colnames(merged_data) == "ID"] <- "ID"
merged_data <- merge(merged_data, ct_data[, c("ID", "site", "sex", "cag", "group", "age1", "burden1", "eTIV")], by = "ID", all.x = TRUE)
desired_order <- c("ID", "site", "sex", "cag", "group", "age1", "burden1", "eTIV", names(merged_data)[!names(merged_data) %in% c("ID", "site", "sex", "cag", "group", "age1", "burden1", "eTIV")])
merged_data <- merged_data[, desired_order]

##Estimating Cohen's D##
colnames <-colnames(merged_data[9:326])
k=9

test=0.12
df<-data.frame(test)
write_csv(df,'/Users/charlie/Desktop/my_projects/neurotransmitter/R/parcellate/results/cohens_d_DK318.csv')

for (i in colnames) {
  setwd('/Users/charlie/Desktop/my_projects/neurotransmitter/R/parcellate/results/')
  
  # Access the column in the dataframe using the column name
  print(i)
  roi<-merged_data[, i]
  
  # Linear mixed effects model that describes cortical thickness
  # The T value here is the adjusted T 
  model <- lm( roi ~ merged_data$group + merged_data$age1 + merged_data$sex + merged_data$site + merged_data$eTIV)
  df.residual(model)
  m<-summary(model)
  print(m)
  
  
  estimate<-m$coefficients[ , 1]
  se<-m$coefficients[ , 2]
  adjusted_tval<-m$coefficients[ , 3]
  p<-m$coefficients[ , 4]
  
  #Get the adjusted T value for the ROI
  adjusted_tval_cxthick<-adjusted_tval[2]
  print(adjusted_tval_cxthick)
  
  
  # To estimate  the partial D value (adjusted for covariantes)
  # We will need the adjusted T from  lm()
  adjusted_d<-t_to_d(adjusted_tval_cxthick, df_error=df.residual(model))[[1]]
  adjusted_d_ci_low<-t_to_d(adjusted_tval_cxthick, df_error=df.residual(model))[[3]]
  adjusted_d_ci_high<-t_to_d(adjusted_tval_cxthick, df_error=df.residual(model))[[4]]
  
  #Now to estimate unadjusted T and D
  # First get the data in patients and controls
  x<- merged_data  %>% dplyr::select(k)
  patients<-subset(x, merged_data$group == "2")
  controls<-subset(x, merged_data$group == "0")
  
  patients_mean<-mean(patients[[1]])
  controls_mean<-mean(controls[[1]])
  
  s1<-sd(patients[[1]])
  s2<-sd(controls[[1]])
  n1<-length(patients[[1]])
  n2<-length(controls[[1]])
  pooled <- sqrt(((n1-1)*s1^2 + (n2-1)*s2^2) / (n1+n1-2))
  
  #Plain T value
  unadj_Tval <- (patients_mean - controls_mean) / pooled
  
  #Plain Cohen's D
  unadj_d<-cohensD(x = patients[[1]], y = controls[[1]])
  
  z = data.frame(D_adjusted = adjusted_d, adjusted_d_ci_low = adjusted_d_ci_low, 
                 adjusted_d_ci_high, D_unadjusted = unadj_d, 
                 T_adjusted = adjusted_tval_cxthick,T_unadjusted = unadj_Tval )
  rownames(z)<-c(print(i))
  write.csv(z, paste0(i))
  
  df_d<-read.csv('/Users/charlie/Desktop/my_projects/neurotransmitter/R/parcellate/results/cohens_d_DK318.csv')
  df_d2<-data.frame(adjusted_d)
  colnames(df_d2)[1] = c(print(i))
  df_d3<-cbind(df_d, df_d2)
  write.csv (df_d3, '/Users/charlie/Desktop/my_projects/neurotransmitter/R/parcellate/results/cohens_d_DK318.csv')
  
  if (k<327) {k = k+1} else print('analysis done') 
  
}
#Remove columns that start with "X" or have non cortical items
h<- read.csv("/Users/charlie/Desktop/my_projects/neurotransmitter/R/parcellate/results/cohens_d_DK318.csv")
b<- h %>% dplyr::select(-(contains("X") | contains('Mean')|contains('Caudate')|contains('Putamen')|contains('cortex')|contains('icv')|contains('aparc')|contains('test')))
b_trans<-t(b)
write.csv(b_trans, '/Users/charlie/Desktop/my_projects/neurotransmitter/R/parcellate/results/cohens_d_transp_DK318.csv')

















###### SCHAEFFER 100 #####

# Initialize an empty data frame to store the merged data
merged_data <- data.frame()


# Loop through each subject and process their data
for (subject_id in subject_ids) {
  # Construct the file path
  file_path <- file.path(data_directory, subject_id, "mri", "mind_Schaefer2018_100Parcels_7Networks_order.csv")
  
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



merged_data_Schaeffer_100<-merged_data
write_csv(merged_data_Schaeffer_100,'/Users/charlie/Desktop/my_projects/neurotransmitter/R/parcellate/results/merged_data_Schaeffer_100.csv')

avg_rois <- rowMeans(merged_data[, -1])
avg_sample_Schaeffer_100<-mean(avg_rois)


colnames(merged_data)[colnames(merged_data) == "ID"] <- "ID"
merged_data <- merge(merged_data, ct_data[, c("ID", "site", "sex", "cag", "group", "age1", "burden1", "eTIV")], by = "ID", all.x = TRUE)
desired_order <- c("ID", "site", "sex", "cag", "group", "age1", "burden1", "eTIV", names(merged_data)[!names(merged_data) %in% c("ID", "site", "sex", "cag", "group", "age1", "burden1", "eTIV")])
merged_data <- merged_data[, desired_order]

##Estimating Cohen's D##
colnames <-colnames(merged_data[9:108])
k=9

test=0.12
df<-data.frame(test)
write_csv(df,'/Users/charlie/Desktop/my_projects/neurotransmitter/R/parcellate/results/cohens_d_Schaeffer_100.csv')

for (i in colnames) {
  setwd('/Users/charlie/Desktop/my_projects/neurotransmitter/R/parcellate/results/')
  
  # Access the column in the dataframe using the column name
  print(i)
  roi<-merged_data[, i]
  
  # Linear mixed effects model that describes cortical thickness
  # The T value here is the adjusted T 
  model <- lm( roi ~ merged_data$group + merged_data$age1 + merged_data$sex + merged_data$site + merged_data$eTIV)
  df.residual(model)
  m<-summary(model)
  print(m)
  
  
  estimate<-m$coefficients[ , 1]
  se<-m$coefficients[ , 2]
  adjusted_tval<-m$coefficients[ , 3]
  p<-m$coefficients[ , 4]
  
  #Get the adjusted T value for the ROI
  adjusted_tval_cxthick<-adjusted_tval[2]
  print(adjusted_tval_cxthick)
  
  
  # To estimate  the partial D value (adjusted for covariantes)
  # We will need the adjusted T from  lm()
  adjusted_d<-t_to_d(adjusted_tval_cxthick, df_error=df.residual(model))[[1]]
  adjusted_d_ci_low<-t_to_d(adjusted_tval_cxthick, df_error=df.residual(model))[[3]]
  adjusted_d_ci_high<-t_to_d(adjusted_tval_cxthick, df_error=df.residual(model))[[4]]
  
  #Now to estimate unadjusted T and D
  # First get the data in patients and controls
  x<- merged_data  %>% dplyr::select(k)
  patients<-subset(x, merged_data$group == "2")
  controls<-subset(x, merged_data$group == "0")
  
  patients_mean<-mean(patients[[1]])
  controls_mean<-mean(controls[[1]])
  
  s1<-sd(patients[[1]])
  s2<-sd(controls[[1]])
  n1<-length(patients[[1]])
  n2<-length(controls[[1]])
  pooled <- sqrt(((n1-1)*s1^2 + (n2-1)*s2^2) / (n1+n1-2))
  
  #Plain T value
  unadj_Tval <- (patients_mean - controls_mean) / pooled
  
  #Plain Cohen's D
  unadj_d<-cohensD(x = patients[[1]], y = controls[[1]])
  
  z = data.frame(D_adjusted = adjusted_d, adjusted_d_ci_low = adjusted_d_ci_low, 
                 adjusted_d_ci_high, D_unadjusted = unadj_d, 
                 T_adjusted = adjusted_tval_cxthick,T_unadjusted = unadj_Tval )
  rownames(z)<-c(print(i))
  write.csv(z, paste0(i))
  
  df_d<-read.csv('/Users/charlie/Desktop/my_projects/neurotransmitter/R/parcellate/results/cohens_d_Schaeffer_100.csv')
  df_d2<-data.frame(adjusted_d)
  colnames(df_d2)[1] = c(print(i))
  df_d3<-cbind(df_d, df_d2)
  write.csv (df_d3, '/Users/charlie/Desktop/my_projects/neurotransmitter/R/parcellate/results/cohens_d_Schaeffer_100.csv')
  
  if (k<107) {k = k+1} else print('analysis done') 
  
}
#Remove columns that start with "X" or have non cortical items
h<- read.csv("/Users/charlie/Desktop/my_projects/neurotransmitter/R/parcellate/results/cohens_d_Schaeffer_100.csv")
b<- h %>% dplyr::select(-(contains("X") | contains('Mean')|contains('Caudate')|contains('Putamen')|contains('cortex')|contains('icv')|contains('aparc')|contains('test')))
b_trans<-t(b)
write.csv(b_trans, '/Users/charlie/Desktop/my_projects/neurotransmitter/R/parcellate/results/cohens_d_transp_Schaeffer_100.csv')












###### SCHAEFFER 200 #####

# Initialize an empty data frame to store the merged data
merged_data <- data.frame()


# Loop through each subject and process their data
for (subject_id in subject_ids) {
  # Construct the file path
  file_path <- file.path(data_directory, subject_id, "mri", "mind_Schaefer2018_200Parcels_7Networks_order.csv")
  
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

merged_data_Schaeffer_200<-merged_data
write_csv(merged_data_Schaeffer_200,'/Users/charlie/Desktop/my_projects/neurotransmitter/R/parcellate/results/merged_data_Schaeffer_200.csv')


avg_rois <- rowMeans(merged_data[, -1])
avg_sample_Schaeffer_200<-mean(avg_rois)



colnames(merged_data)[colnames(merged_data) == "ID"] <- "ID"
merged_data <- merge(merged_data, ct_data[, c("ID", "site", "sex", "cag", "group", "age1", "burden1", "eTIV")], by = "ID", all.x = TRUE)
desired_order <- c("ID", "site", "sex", "cag", "group", "age1", "burden1", "eTIV", names(merged_data)[!names(merged_data) %in% c("ID", "site", "sex", "cag", "group", "age1", "burden1", "eTIV")])
merged_data <- merged_data[, desired_order]

##Estimating Cohen's D##
colnames <-colnames(merged_data[9:208])
k=9

test=0.12
df<-data.frame(test)
write_csv(df,'/Users/charlie/Desktop/my_projects/neurotransmitter/R/parcellate/results/cohens_d_Schaeffer_200.csv')

for (i in colnames) {
  setwd('/Users/charlie/Desktop/my_projects/neurotransmitter/R/parcellate/results/')
  
  # Access the column in the dataframe using the column name
  print(i)
  roi<-merged_data[, i]
  
  # Linear mixed effects model that describes cortical thickness
  # The T value here is the adjusted T 
  model <- lm( roi ~ merged_data$group + merged_data$age1 + merged_data$sex + merged_data$site + merged_data$eTIV)
  df.residual(model)
  m<-summary(model)
  print(m)
  
  
  estimate<-m$coefficients[ , 1]
  se<-m$coefficients[ , 2]
  adjusted_tval<-m$coefficients[ , 3]
  p<-m$coefficients[ , 4]
  
  #Get the adjusted T value for the ROI
  adjusted_tval_cxthick<-adjusted_tval[2]
  print(adjusted_tval_cxthick)
  
  
  # To estimate  the partial D value (adjusted for covariantes)
  # We will need the adjusted T from  lm()
  adjusted_d<-t_to_d(adjusted_tval_cxthick, df_error=df.residual(model))[[1]]
  adjusted_d_ci_low<-t_to_d(adjusted_tval_cxthick, df_error=df.residual(model))[[3]]
  adjusted_d_ci_high<-t_to_d(adjusted_tval_cxthick, df_error=df.residual(model))[[4]]
  
  #Now to estimate unadjusted T and D
  # First get the data in patients and controls
  x<- merged_data  %>% dplyr::select(k)
  patients<-subset(x, merged_data$group == "2")
  controls<-subset(x, merged_data$group == "0")
  
  patients_mean<-mean(patients[[1]])
  controls_mean<-mean(controls[[1]])
  
  s1<-sd(patients[[1]])
  s2<-sd(controls[[1]])
  n1<-length(patients[[1]])
  n2<-length(controls[[1]])
  pooled <- sqrt(((n1-1)*s1^2 + (n2-1)*s2^2) / (n1+n1-2))
  
  #Plain T value
  unadj_Tval <- (patients_mean - controls_mean) / pooled
  
  #Plain Cohen's D
  unadj_d<-cohensD(x = patients[[1]], y = controls[[1]])
  
  z = data.frame(D_adjusted = adjusted_d, adjusted_d_ci_low = adjusted_d_ci_low, 
                 adjusted_d_ci_high, D_unadjusted = unadj_d, 
                 T_adjusted = adjusted_tval_cxthick,T_unadjusted = unadj_Tval )
  rownames(z)<-c(print(i))
  write.csv(z, paste0(i))
  
  df_d<-read.csv('/Users/charlie/Desktop/my_projects/neurotransmitter/R/parcellate/results/cohens_d_Schaeffer_200.csv')
  df_d2<-data.frame(adjusted_d)
  colnames(df_d2)[1] = c(print(i))
  df_d3<-cbind(df_d, df_d2)
  write.csv (df_d3, '/Users/charlie/Desktop/my_projects/neurotransmitter/R/parcellate/results/cohens_d_Schaeffer_200.csv')
  
  if (k<207) {k = k+1} else print('analysis done') 
  
}
#Remove columns that start with "X" or have non cortical items
h<- read.csv("/Users/charlie/Desktop/my_projects/neurotransmitter/R/parcellate/results/cohens_d_Schaeffer_200.csv")
b<- h %>% dplyr::select(-(contains("X") | contains('Mean')|contains('Caudate')|contains('Putamen')|contains('cortex')|contains('icv')|contains('aparc')|contains('test')))
b_trans<-t(b)
write.csv(b_trans, '/Users/charlie/Desktop/my_projects/neurotransmitter/R/parcellate/results/cohens_d_transp_Schaeffer_200.csv')






###### SCHAEFFER 400 #####

# Initialize an empty data frame to store the merged data
merged_data <- data.frame()


# Loop through each subject and process their data
for (subject_id in subject_ids) {
  # Construct the file path
  file_path <- file.path(data_directory, subject_id, "mri", "mind_Schaefer2018_400Parcels_7Networks_order.csv")
  
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


merged_data_Schaeffer_400<-merged_data
write_csv(merged_data_Schaeffer_400,'/Users/charlie/Desktop/my_projects/neurotransmitter/R/parcellate/results/merged_data_Schaeffer_400.csv')

avg_rois <- rowMeans(merged_data[, -1])
avg_sample_Schaeffer_400<-mean(avg_rois)


colnames(merged_data)[colnames(merged_data) == "ID"] <- "ID"
merged_data <- merge(merged_data, ct_data[, c("ID", "site", "sex", "cag", "group", "age1", "burden1", "eTIV")], by = "ID", all.x = TRUE)
desired_order <- c("ID", "site", "sex", "cag", "group", "age1", "burden1", "eTIV", names(merged_data)[!names(merged_data) %in% c("ID", "site", "sex", "cag", "group", "age1", "burden1", "eTIV")])
merged_data <- merged_data[, desired_order]

##Estimating Cohen's D##
colnames <-colnames(merged_data[9:408])
k=9

test=0.12
df<-data.frame(test)
write_csv(df,'/Users/charlie/Desktop/my_projects/neurotransmitter/R/parcellate/results/cohens_d_Schaeffer_400.csv')

for (i in colnames) {
  setwd('/Users/charlie/Desktop/my_projects/neurotransmitter/R/parcellate/results/')
  
  # Access the column in the dataframe using the column name
  print(i)
  roi<-merged_data[, i]
  
  # Linear mixed effects model that describes cortical thickness
  # The T value here is the adjusted T 
  model <- lm( roi ~ merged_data$group + merged_data$age1 + merged_data$sex + merged_data$site + merged_data$eTIV)
  df.residual(model)
  m<-summary(model)
  print(m)
  
  
  estimate<-m$coefficients[ , 1]
  se<-m$coefficients[ , 2]
  adjusted_tval<-m$coefficients[ , 3]
  p<-m$coefficients[ , 4]
  
  #Get the adjusted T value for the ROI
  adjusted_tval_cxthick<-adjusted_tval[2]
  print(adjusted_tval_cxthick)
  
  
  # To estimate  the partial D value (adjusted for covariantes)
  # We will need the adjusted T from  lm()
  adjusted_d<-t_to_d(adjusted_tval_cxthick, df_error=df.residual(model))[[1]]
  adjusted_d_ci_low<-t_to_d(adjusted_tval_cxthick, df_error=df.residual(model))[[3]]
  adjusted_d_ci_high<-t_to_d(adjusted_tval_cxthick, df_error=df.residual(model))[[4]]
  
  #Now to estimate unadjusted T and D
  # First get the data in patients and controls
  x<- merged_data  %>% dplyr::select(k)
  patients<-subset(x, merged_data$group == "2")
  controls<-subset(x, merged_data$group == "0")
  
  patients_mean<-mean(patients[[1]])
  controls_mean<-mean(controls[[1]])
  
  s1<-sd(patients[[1]])
  s2<-sd(controls[[1]])
  n1<-length(patients[[1]])
  n2<-length(controls[[1]])
  pooled <- sqrt(((n1-1)*s1^2 + (n2-1)*s2^2) / (n1+n1-2))
  
  #Plain T value
  unadj_Tval <- (patients_mean - controls_mean) / pooled
  
  #Plain Cohen's D
  unadj_d<-cohensD(x = patients[[1]], y = controls[[1]])
  
  z = data.frame(D_adjusted = adjusted_d, adjusted_d_ci_low = adjusted_d_ci_low, 
                 adjusted_d_ci_high, D_unadjusted = unadj_d, 
                 T_adjusted = adjusted_tval_cxthick,T_unadjusted = unadj_Tval )
  rownames(z)<-c(print(i))
  write.csv(z, paste0(i))
  
  df_d<-read.csv('/Users/charlie/Desktop/my_projects/neurotransmitter/R/parcellate/results/cohens_d_Schaeffer_400.csv')
  df_d2<-data.frame(adjusted_d)
  colnames(df_d2)[1] = c(print(i))
  df_d3<-cbind(df_d, df_d2)
  write.csv (df_d3, '/Users/charlie/Desktop/my_projects/neurotransmitter/R/parcellate/results/cohens_d_Schaeffer_400.csv')
  
  if (k<407) {k = k+1} else print('analysis done') 
  
}
#Remove columns that start with "X" or have non cortical items
h<- read.csv("/Users/charlie/Desktop/my_projects/neurotransmitter/R/parcellate/results/cohens_d_Schaeffer_400.csv")
b<- h %>% dplyr::select(-(contains("X") | contains('Mean')|contains('Caudate')|contains('Putamen')|contains('cortex')|contains('icv')|contains('aparc')|contains('test')))
b_trans<-t(b)
write.csv(b_trans, '/Users/charlie/Desktop/my_projects/neurotransmitter/R/parcellate/results/cohens_d_transp_Schaeffer_400.csv')






###### SCHAEFFER 100 - 17 Networks- #####

# Initialize an empty data frame to store the merged data
merged_data <- data.frame()


# Loop through each subject and process their data
for (subject_id in subject_ids) {
  # Construct the file path
  file_path <- file.path(data_directory, subject_id, "mri", "mind_Schaefer2018_100Parcels_17Networks_order.csv")
  
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


merged_data_Schaeffer_100_17<-merged_data
write_csv(merged_data_Schaeffer_100_17,'/Users/charlie/Desktop/my_projects/neurotransmitter/R/parcellate/results/merged_data_Schaeffer_100_17.csv')


avg_rois <- rowMeans(merged_data[, -1])
avg_sample_Schaeffer_100_17<-mean(avg_rois)



colnames(merged_data)[colnames(merged_data) == "ID"] <- "ID"
merged_data <- merge(merged_data, ct_data[, c("ID", "site", "sex", "cag", "group", "age1", "burden1", "eTIV")], by = "ID", all.x = TRUE)
desired_order <- c("ID", "site", "sex", "cag", "group", "age1", "burden1", "eTIV", names(merged_data)[!names(merged_data) %in% c("ID", "site", "sex", "cag", "group", "age1", "burden1", "eTIV")])
merged_data <- merged_data[, desired_order]

##Estimating Cohen's D##
colnames <-colnames(merged_data[9:108])
k=9

test=0.12
df<-data.frame(test)
write_csv(df,'/Users/charlie/Desktop/my_projects/neurotransmitter/R/parcellate/results/cohens_d_Schaeffer_100_17.csv')

for (i in colnames) {
  setwd('/Users/charlie/Desktop/my_projects/neurotransmitter/R/parcellate/results/')
  
  # Access the column in the dataframe using the column name
  print(i)
  roi<-merged_data[, i]
  
  # Linear mixed effects model that describes cortical thickness
  # The T value here is the adjusted T 
  model <- lm( roi ~ merged_data$group + merged_data$age1 + merged_data$sex + merged_data$site + merged_data$eTIV)
  df.residual(model)
  m<-summary(model)
  print(m)
  
  
  estimate<-m$coefficients[ , 1]
  se<-m$coefficients[ , 2]
  adjusted_tval<-m$coefficients[ , 3]
  p<-m$coefficients[ , 4]
  
  #Get the adjusted T value for the ROI
  adjusted_tval_cxthick<-adjusted_tval[2]
  print(adjusted_tval_cxthick)
  
  
  # To estimate  the partial D value (adjusted for covariantes)
  # We will need the adjusted T from  lm()
  adjusted_d<-t_to_d(adjusted_tval_cxthick, df_error=df.residual(model))[[1]]
  adjusted_d_ci_low<-t_to_d(adjusted_tval_cxthick, df_error=df.residual(model))[[3]]
  adjusted_d_ci_high<-t_to_d(adjusted_tval_cxthick, df_error=df.residual(model))[[4]]
  
  #Now to estimate unadjusted T and D
  # First get the data in patients and controls
  x<- merged_data  %>% dplyr::select(k)
  patients<-subset(x, merged_data$group == "2")
  controls<-subset(x, merged_data$group == "0")
  
  patients_mean<-mean(patients[[1]])
  controls_mean<-mean(controls[[1]])
  
  s1<-sd(patients[[1]])
  s2<-sd(controls[[1]])
  n1<-length(patients[[1]])
  n2<-length(controls[[1]])
  pooled <- sqrt(((n1-1)*s1^2 + (n2-1)*s2^2) / (n1+n1-2))
  
  #Plain T value
  unadj_Tval <- (patients_mean - controls_mean) / pooled
  
  #Plain Cohen's D
  unadj_d<-cohensD(x = patients[[1]], y = controls[[1]])
  
  z = data.frame(D_adjusted = adjusted_d, adjusted_d_ci_low = adjusted_d_ci_low, 
                 adjusted_d_ci_high, D_unadjusted = unadj_d, 
                 T_adjusted = adjusted_tval_cxthick,T_unadjusted = unadj_Tval )
  rownames(z)<-c(print(i))
  write.csv(z, paste0(i))
  
  df_d<-read.csv('/Users/charlie/Desktop/my_projects/neurotransmitter/R/parcellate/results/cohens_d_Schaeffer_100_17.csv')
  df_d2<-data.frame(adjusted_d)
  colnames(df_d2)[1] = c(print(i))
  df_d3<-cbind(df_d, df_d2)
  write.csv (df_d3, '/Users/charlie/Desktop/my_projects/neurotransmitter/R/parcellate/results/cohens_d_Schaeffer_100_17.csv')
  
  if (k<107) {k = k+1} else print('analysis done') 
  
}
#Remove columns that start with "X" or have non cortical items
h<- read.csv("/Users/charlie/Desktop/my_projects/neurotransmitter/R/parcellate/results/cohens_d_Schaeffer_100_17.csv")
b<- h %>% dplyr::select(-(contains("X.9") | contains("X.8") |contains("X.7") |contains("X.6") |contains("X.5") |contains("X.4") |contains("X.3") |contains("X.2") |contains("X.1") |contains('Mean')|contains('Caudate')|contains('Putamen')|contains('cortex')|contains('icv')|contains('aparc')|contains('test')))%>%
  select(-X)
b_trans<-t(b)
write.csv(b_trans, '/Users/charlie/Desktop/my_projects/neurotransmitter/R/parcellate/results/cohens_d_transp_Schaeffer_100_17.csv')








###### SCHAEFFER 500 - 17 Networks- #####


# Initialize an empty data frame to store the merged data
merged_data <- data.frame()


# Loop through each subject and process their data
for (subject_id in subject_ids) {
  # Construct the file path
  file_path <- file.path(data_directory, subject_id, "mri", "mind_Schaefer2018_500Parcels_17Networks_order.csv")
  
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


merged_data_Schaeffer_500_17<-merged_data
write_csv(merged_data_Schaeffer_500_17,'/Users/charlie/Desktop/my_projects/neurotransmitter/R/parcellate/results/merged_data_Schaeffer_500_17.csv')

avg_rois <- rowMeans(merged_data[, -1])
avg_sample_Schaeffer_500_17<-mean(avg_rois)



colnames(merged_data)[colnames(merged_data) == "ID"] <- "ID"
merged_data <- merge(merged_data, ct_data[, c("ID", "site", "sex", "cag", "group", "age1", "burden1", "eTIV")], by = "ID", all.x = TRUE)
desired_order <- c("ID", "site", "sex", "cag", "group", "age1", "burden1", "eTIV", names(merged_data)[!names(merged_data) %in% c("ID", "site", "sex", "cag", "group", "age1", "burden1", "eTIV")])
merged_data <- merged_data[, desired_order]

##Estimating Cohen's D##
colnames <-colnames(merged_data[9:508])
k=9

test=0.12
df<-data.frame(test)
write_csv(df,'/Users/charlie/Desktop/my_projects/neurotransmitter/R/parcellate/results/cohens_d_Schaeffer_500_17.csv')

for (i in colnames) {
  setwd('/Users/charlie/Desktop/my_projects/neurotransmitter/R/parcellate/results/')
  
  # Access the column in the dataframe using the column name
  print(i)
  roi<-merged_data[, i]
  
  # Linear mixed effects model that describes cortical thickness
  # The T value here is the adjusted T 
  model <- lm( roi ~ merged_data$group + merged_data$age1 + merged_data$sex + merged_data$site + merged_data$eTIV)
  df.residual(model)
  m<-summary(model)
  print(m)
  
  
  estimate<-m$coefficients[ , 1]
  se<-m$coefficients[ , 2]
  adjusted_tval<-m$coefficients[ , 3]
  p<-m$coefficients[ , 4]
  
  #Get the adjusted T value for the ROI
  adjusted_tval_cxthick<-adjusted_tval[2]
  print(adjusted_tval_cxthick)
  
  
  # To estimate  the partial D value (adjusted for covariantes)
  # We will need the adjusted T from  lm()
  adjusted_d<-t_to_d(adjusted_tval_cxthick, df_error=df.residual(model))[[1]]
  adjusted_d_ci_low<-t_to_d(adjusted_tval_cxthick, df_error=df.residual(model))[[3]]
  adjusted_d_ci_high<-t_to_d(adjusted_tval_cxthick, df_error=df.residual(model))[[4]]
  
  #Now to estimate unadjusted T and D
  # First get the data in patients and controls
  x<- merged_data  %>% dplyr::select(k)
  patients<-subset(x, merged_data$group == "2")
  controls<-subset(x, merged_data$group == "0")
  
  patients_mean<-mean(patients[[1]])
  controls_mean<-mean(controls[[1]])
  
  s1<-sd(patients[[1]])
  s2<-sd(controls[[1]])
  n1<-length(patients[[1]])
  n2<-length(controls[[1]])
  pooled <- sqrt(((n1-1)*s1^2 + (n2-1)*s2^2) / (n1+n1-2))
  
  #Plain T value
  unadj_Tval <- (patients_mean - controls_mean) / pooled
  
  #Plain Cohen's D
  unadj_d<-cohensD(x = patients[[1]], y = controls[[1]])
  
  z = data.frame(D_adjusted = adjusted_d, adjusted_d_ci_low = adjusted_d_ci_low, 
                 adjusted_d_ci_high, D_unadjusted = unadj_d, 
                 T_adjusted = adjusted_tval_cxthick,T_unadjusted = unadj_Tval )
  rownames(z)<-c(print(i))
  write.csv(z, paste0(i))
  
  df_d<-read.csv('/Users/charlie/Desktop/my_projects/neurotransmitter/R/parcellate/results/cohens_d_Schaeffer_500_17.csv')
  df_d2<-data.frame(adjusted_d)
  colnames(df_d2)[1] = c(print(i))
  df_d3<-cbind(df_d, df_d2)
  write.csv (df_d3, '/Users/charlie/Desktop/my_projects/neurotransmitter/R/parcellate/results/cohens_d_Schaeffer_500_17.csv')
  
  if (k<507) {k = k+1} else print('analysis done') 
  
}
#Remove columns that start with "X" or have non cortical items
h<- read.csv("/Users/charlie/Desktop/my_projects/neurotransmitter/R/parcellate/results/cohens_d_Schaeffer_500_17.csv")
b<- h %>% dplyr::select(-(contains("X.9") | contains("X.8") |contains("X.7") |contains("X.6") |contains("X.5") |contains("X.4") |contains("X.3") |contains("X.2") |contains("X.1") |contains('Mean')|contains('Caudate')|contains('Putamen')|contains('cortex')|contains('icv')|contains('aparc')|contains('test')))%>%
  select(-X)
b_trans<-t(b)
write.csv(b_trans, '/Users/charlie/Desktop/my_projects/neurotransmitter/R/parcellate/results/cohens_d_transp_Schaeffer_500_17.csv')











###### SCHAEFFER 1000 - 17 Networks- #####

# Initialize an empty data frame to store the merged data
merged_data <- data.frame()


# Loop through each subject and process their data
for (subject_id in subject_ids) {
  # Construct the file path
  file_path <- file.path(data_directory, subject_id, "mri", "mind_Schaefer2018_1000Parcels_17Networks_order.csv")
  
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


merged_data_Schaeffer_1000_17<-merged_data
write_csv(merged_data_Schaeffer_1000_17,'/Users/charlie/Desktop/my_projects/neurotransmitter/R/parcellate/results/merged_data_Schaeffer_1000_17.csv')


avg_rois <- rowMeans(merged_data[, -1])
avg_sample_Schaeffer_1000_17<-mean(avg_rois)




colnames(merged_data)[colnames(merged_data) == "ID"] <- "ID"
merged_data <- merge(merged_data, ct_data[, c("ID", "site", "sex", "cag", "group", "age1", "burden1", "eTIV")], by = "ID", all.x = TRUE)
desired_order <- c("ID", "site", "sex", "cag", "group", "age1", "burden1", "eTIV", names(merged_data)[!names(merged_data) %in% c("ID", "site", "sex", "cag", "group", "age1", "burden1", "eTIV")])
merged_data <- merged_data[, desired_order]

##Estimating Cohen's D##
colnames <-colnames(merged_data[9:1008])
k=9

test=0.12
df<-data.frame(test)
write_csv(df,'/Users/charlie/Desktop/my_projects/neurotransmitter/R/parcellate/results/cohens_d_Schaeffer_1000_17.csv')

for (i in colnames) {
  setwd('/Users/charlie/Desktop/my_projects/neurotransmitter/R/parcellate/results/')
  
  # Access the column in the dataframe using the column name
  print(i)
  roi<-merged_data[, i]
  
  # Linear mixed effects model that describes cortical thickness
  # The T value here is the adjusted T 
  model <- lm( roi ~ merged_data$group + merged_data$age1 + merged_data$sex + merged_data$site + merged_data$eTIV)
  df.residual(model)
  m<-summary(model)
  print(m)
  
  
  estimate<-m$coefficients[ , 1]
  se<-m$coefficients[ , 2]
  adjusted_tval<-m$coefficients[ , 3]
  p<-m$coefficients[ , 4]
  
  #Get the adjusted T value for the ROI
  adjusted_tval_cxthick<-adjusted_tval[2]
  print(adjusted_tval_cxthick)
  
  
  # To estimate  the partial D value (adjusted for covariantes)
  # We will need the adjusted T from  lm()
  adjusted_d<-t_to_d(adjusted_tval_cxthick, df_error=df.residual(model))[[1]]
  adjusted_d_ci_low<-t_to_d(adjusted_tval_cxthick, df_error=df.residual(model))[[3]]
  adjusted_d_ci_high<-t_to_d(adjusted_tval_cxthick, df_error=df.residual(model))[[4]]
  
  #Now to estimate unadjusted T and D
  # First get the data in patients and controls
  x<- merged_data  %>% dplyr::select(k)
  patients<-subset(x, merged_data$group == "2")
  controls<-subset(x, merged_data$group == "0")
  
  patients_mean<-mean(patients[[1]])
  controls_mean<-mean(controls[[1]])
  
  s1<-sd(patients[[1]])
  s2<-sd(controls[[1]])
  n1<-length(patients[[1]])
  n2<-length(controls[[1]])
  pooled <- sqrt(((n1-1)*s1^2 + (n2-1)*s2^2) / (n1+n1-2))
  
  #Plain T value
  unadj_Tval <- (patients_mean - controls_mean) / pooled
  
  #Plain Cohen's D
  unadj_d<-cohensD(x = patients[[1]], y = controls[[1]])
  
  z = data.frame(D_adjusted = adjusted_d, adjusted_d_ci_low = adjusted_d_ci_low, 
                 adjusted_d_ci_high, D_unadjusted = unadj_d, 
                 T_adjusted = adjusted_tval_cxthick,T_unadjusted = unadj_Tval )
  rownames(z)<-c(print(i))
  write.csv(z, paste0(i))
  
  df_d<-read.csv('/Users/charlie/Desktop/my_projects/neurotransmitter/R/parcellate/results/cohens_d_Schaeffer_1000_17.csv')
  df_d2<-data.frame(adjusted_d)
  colnames(df_d2)[1] = c(print(i))
  df_d3<-cbind(df_d, df_d2)
  write.csv (df_d3, '/Users/charlie/Desktop/my_projects/neurotransmitter/R/parcellate/results/cohens_d_Schaeffer_1000_17.csv')
  
  if (k<107) {k = k+1} else print('analysis done') 
  
}
#Remove columns that start with "X" or have non cortical items
h<- read.csv("/Users/charlie/Desktop/my_projects/neurotransmitter/R/parcellate/results/cohens_d_Schaeffer_1000_17.csv")
b<- h %>% dplyr::select(-(contains("X.9") | contains("X.8") |contains("X.7") |contains("X.6") |contains("X.5") |contains("X.4") |contains("X.3") |contains("X.2") |contains("X.1") |contains('Mean')|contains('Caudate')|contains('Putamen')|contains('cortex')|contains('icv')|contains('aparc')|contains('test')))%>%
  select(-X)
b_trans<-t(b)
write.csv(b_trans, '/Users/charlie/Desktop/my_projects/neurotransmitter/R/parcellate/results/cohens_d_transp_Schaeffer_1000_17.csv')









###### HCP #####

# Initialize an empty data frame to store the merged data
merged_data <- data.frame()


# Loop through each subject and process their data
for (subject_id in subject_ids) {
  # Construct the file path
  file_path <- file.path(data_directory, subject_id, "mri", "mind_HCP.csv")
  
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


merged_data_HCP<-merged_data
write_csv(merged_data_HCP,'/Users/charlie/Desktop/my_projects/neurotransmitter/R/parcellate/results/merged_data_HCP.csv')


avg_rois <- rowMeans(merged_data[, -1])
avg_sample_HCP<-mean(avg_rois)



colnames(merged_data)[colnames(merged_data) == "ID"] <- "ID"
merged_data <- merge(merged_data, ct_data[, c("ID", "site", "sex", "cag", "group", "age1", "burden1", "eTIV")], by = "ID", all.x = TRUE)
desired_order <- c("ID", "site", "sex", "cag", "group", "age1", "burden1", "eTIV", names(merged_data)[!names(merged_data) %in% c("ID", "site", "sex", "cag", "group", "age1", "burden1", "eTIV")])
merged_data <- merged_data[, desired_order]

##Estimating Cohen's D##
colnames <-colnames(merged_data[9:368])
k=9

test=0.12
df<-data.frame(test)
write_csv(df,'/Users/charlie/Desktop/my_projects/neurotransmitter/R/parcellate/results/cohens_d_HCP.csv')

for (i in colnames) {
  setwd('/Users/charlie/Desktop/my_projects/neurotransmitter/R/parcellate/results/')
  
  # Access the column in the dataframe using the column name
  print(i)
  roi<-merged_data[, i]
  
  # Linear mixed effects model that describes cortical thickness
  # The T value here is the adjusted T 
  model <- lm( roi ~ merged_data$group + merged_data$age1 + merged_data$sex + merged_data$site + merged_data$eTIV)
  df.residual(model)
  m<-summary(model)
  print(m)
  
  
  estimate<-m$coefficients[ , 1]
  se<-m$coefficients[ , 2]
  adjusted_tval<-m$coefficients[ , 3]
  p<-m$coefficients[ , 4]
  
  #Get the adjusted T value for the ROI
  adjusted_tval_cxthick<-adjusted_tval[2]
  print(adjusted_tval_cxthick)
  
  
  # To estimate  the partial D value (adjusted for covariantes)
  # We will need the adjusted T from  lm()
  adjusted_d<-t_to_d(adjusted_tval_cxthick, df_error=df.residual(model))[[1]]
  adjusted_d_ci_low<-t_to_d(adjusted_tval_cxthick, df_error=df.residual(model))[[3]]
  adjusted_d_ci_high<-t_to_d(adjusted_tval_cxthick, df_error=df.residual(model))[[4]]
  
  #Now to estimate unadjusted T and D
  # First get the data in patients and controls
  x<- merged_data  %>% dplyr::select(k)
  patients<-subset(x, merged_data$group == "2")
  controls<-subset(x, merged_data$group == "0")
  
  patients_mean<-mean(patients[[1]])
  controls_mean<-mean(controls[[1]])
  
  s1<-sd(patients[[1]])
  s2<-sd(controls[[1]])
  n1<-length(patients[[1]])
  n2<-length(controls[[1]])
  pooled <- sqrt(((n1-1)*s1^2 + (n2-1)*s2^2) / (n1+n1-2))
  
  #Plain T value
  unadj_Tval <- (patients_mean - controls_mean) / pooled
  
  #Plain Cohen's D
  unadj_d<-cohensD(x = patients[[1]], y = controls[[1]])
  
  z = data.frame(D_adjusted = adjusted_d, adjusted_d_ci_low = adjusted_d_ci_low, 
                 adjusted_d_ci_high, D_unadjusted = unadj_d, 
                 T_adjusted = adjusted_tval_cxthick,T_unadjusted = unadj_Tval )
  rownames(z)<-c(print(i))
  write.csv(z, paste0(i))
  
  df_d<-read.csv('/Users/charlie/Desktop/my_projects/neurotransmitter/R/parcellate/results/cohens_d_HCP.csv')
  df_d2<-data.frame(adjusted_d)
  colnames(df_d2)[1] = c(print(i))
  df_d3<-cbind(df_d, df_d2)
  write.csv (df_d3, '/Users/charlie/Desktop/my_projects/neurotransmitter/R/parcellate/results/cohens_d_HCP.csv')
  
  if (k<75) {k = k+1} else print('analysis done') 
  
}
#Remove columns that start with "X" or have non cortical items
h<- read.csv("/Users/charlie/Desktop/my_projects/neurotransmitter/R/parcellate/results/cohens_d_HCP.csv")
b<- h %>% dplyr::select(-(contains("X") | contains('Mean')|contains('Caudate')|contains('Putamen')|contains('cortex')|contains('icv')|contains('aparc')|contains('test')))
b_trans<-t(b)
write.csv(b_trans, '/Users/charlie/Desktop/my_projects/neurotransmitter/R/parcellate/results/cohens_d_transp_HCP.csv')









# Create a vector containing the names of the variables you want to save
variables <- c("avg_sample_DK318", "avg_sample_DK68", "avg_sample_HCP", 
               "avg_sample_Schaeffer_100", "avg_sample_Schaeffer_100_17", "avg_sample_Schaeffer_1000_17", 
               "avg_sample_Schaeffer_200", "avg_sample_Schaeffer_400", 
               "avg_sample_Schaeffer_500_17")

# Initialize an empty string to store the content
content_string <- ""

# Loop through each variable and concatenate its name and content
for (variable in variables) {
  # Get the content of the variable
  content <- get(variable)
  # Concatenate the name and content
  content_string <- paste0(content_string, paste0(variable, " = ", content, "\n"))
}

# Set the directory where you want to save the file
directory <- "/Users/charlie/Desktop/my_projects/neurotransmitter/R/parcellate/results"

# Define the file path
file_path <- file.path(directory, "all_variables.txt")

# Write the content string to the text file
writeLines(content_string, file_path)
