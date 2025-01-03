
#Script to obtain Cohen's D from individual MIND similarity matrices


library(effectsize)
library(tidyverse)
library(ggm)
library(ggplot2)
library(tidyr)
library(dplyr)
library(readr)
library(lsr)

rm(list = ls())

##########NODE DEGREE##############

#####Track#####

##Obtaining degree
data_directory <- "/Users/charlie/Desktop/my_projects/neurotransmitter/github/data"


# List of subject names

subject_ids <- c(
  '000-000-001', '000-000-002','000-000-122', '000-000-033','000-000-101','000-000-178', '000-000-055'
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
    subject_avg <- data.frame(hdid = subject_id, subject_avg)
    
    # Append the subject's averaged data to the merged_data data frame
    merged_data <- rbind(merged_data, subject_avg)
  } else {
    cat("File not found for subject:", subject_id, "\n")
  }
}


# Reset column names if needed (assuming the columns in all files are the same)
colnames(merged_data) <- c("ID", colnames(merged_data)[-1])

ct_data <- read.csv("/Users/charlie/Desktop/my_projects/neurotransmitter/github/data/track_fsv7.csv")


merged_data <- merge(merged_data, ct_data[, c("ID", "site", "sex", "cag", "group", "age1", "burden1", "eTIV")], by = "ID", all.x = TRUE)

desired_order <- c("ID", "site", "sex", "cag", "group", "age1", "burden1", "eTIV", names(merged_data)[!names(merged_data) %in% c("ID", "site", "sex", "cag", "group", "age1", "burden1", "eTIV")])

# Reorder the columns in merged_data
merged_data <- merged_data[, desired_order]

write_csv(merged_data,'/Users/charlie/Desktop/my_projects/neurotransmitter/github/results/nodal_degree.csv')

##Estimating Cohen's D##

colnames <-colnames(merged_data[9:76])
k=9

test=0.12
df<-data.frame(test)
write_csv(df,'/Users/charlie/Desktop/my_projects/neurotransmitter/github/results/cohens_d.csv')


for (i in colnames) {
  setwd('/Users/charlie/Desktop/my_projects/neurotransmitter/github/results')
  
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
  
  df_d<-read.csv('/Users/charlie/Desktop/my_projects/neurotransmitter/github/results/cohens_d.csv')
  df_d2<-data.frame(adjusted_d)
  colnames(df_d2)[1] = c(print(i))
  df_d3<-cbind(df_d, df_d2)
  write.csv (df_d3, '/Users/charlie/Desktop/my_projects/neurotransmitter/github/results/cohens_d.csv')
  
  if (k<75) {k = k+1} else print('analysis done') 
  
}

#Remove columns that start with "X" or have non cortical items
h<- read.csv("/Users/charlie/Desktop/my_projects/neurotransmitter/github/results/cohens_d.csv")
b<- h %>% dplyr::select(-(contains("X") | contains('Mean')|contains('Caudate')|contains('Putamen')|contains('cortex')|contains('icv')|contains('aparc')|contains('test')))
b_trans<-t(b)
write.csv(b_trans, '/Users/charlie/Desktop/my_projects/neurotransmitter/github/results/cohens_d_right_transp.csv')



hhh<-data.frame(b_trans)
gg<-ggplot(hhh, aes(x= hhh$b_trans)) +
  geom_density(color="darkblue", fill="lightblue")+ theme_classic()+theme(
    text = element_text(size = 40),  # Set text size for all elements
    axis.title = element_text(size = 30),  # Set title size
    axis.text = element_text(size = 25)  # Set axis labels size
  ) + xlab("Cohen's d") + ylab ("Density") + xlim(-1, 0.5) +  ylim(0, 2.5)     
ggsave(gg, filename = 'cohensD_track_scale.pdf', device = 'pdf', width = 11.69, height = 8.27, dpi = 300, path ='/Users/charlie/Desktop/my_projects/neurotransmitter/github/results')



# Z-score all the cohen's D values
z_scores<-(b_trans-mean(b_trans))/sd(b_trans)  
write.csv(z_scores, '/Users/charlie/Desktop/my_projects/neurotransmitter/github/results/Z_cohens_d.csv')  

hh<-data.frame(z_scores)
g<-ggplot(hh, aes(x= hh$z_scores)) +
  geom_density()
ggsave(g, filename = 'z_scores_track_scale.pdf', device = 'pdf', width = 11.69, height = 8.27, dpi = 300, path ='/Users/charlie/Desktop/my_projects/neurotransmitter/github/results')




