

library(ggplot2)
library(tidyr)
library(dplyr)
library (data.table)
library(readr)

rm(list = ls())


#####YAS#####

##Obtaining degree
data_directory <- "/Users/charlie/Desktop/my_projects/neurotransmitter/imaging_metrics/thickness/yas/data/fsv7"


# List of subject names

subject_ids <- c(
  'XXX-XXX-XXX', 'YYY-YYY-YYY'
)

# Initialize an empty data frame to store the merged data
merged_data_yas <- data.frame()

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
    merged_data_yas <- rbind(merged_data_yas, subject_avg)
  } else {
    cat("File not found for subject:", subject_id, "\n")
  }
} 



transposed_data_yas <- transpose(merged_data_yas)
transposed_data_yas_noname <- transposed_data_yas[-1, ]  # Remove the first row
colnames(transposed_data_yas) <- transposed_data_yas[1, ]  # Use the first row as column names
transposed_data_yas_name <- transposed_data_yas[-1, ]  # Remove the first row

# Reset column names if needed (assuming the columns in all files are the same)
colnames(merged_data_yas) <- c("HDID", colnames(merged_data_yas)[-1])
ct_data_yas <- read.csv("/Users/charlie/Desktop/my_projects/neurotransmitter/imaging_metrics/thickness/yas/results/by_harry/fsv7/yas.csv")
merged_data_yas <- merge(merged_data_yas, ct_data_yas[, c("HDID", "site", "sex", "cag", "group", "age1", "burden1", "eTIV")], by = "HDID", all.x = TRUE)
desired_order <- c("HDID", "site", "sex", "cag", "group", "age1", "burden1", "eTIV", names(merged_data_yas)[!names(merged_data_yas) %in% c("HDID", "site", "sex", "cag", "group", "age1", "burden1", "eTIV")])
# Reorder the columns in merged_data
merged_data_yas <- merged_data_yas[, desired_order]








#####Trackon#####

##Obtaining degree
data_directory <- "/Users/charlie/Desktop/my_projects/neurotransmitter/imaging_metrics/thickness/trackon/carlos"


# List of subject names

subject_ids <- c(
  'XXX-XXX-XXX', 'YYY-YYY-YYY'
)


# Initialize an empty data frame to store the merged data
merged_data_trackon <- data.frame()

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
    merged_data_trackon <- rbind(merged_data_trackon, subject_avg)
  } else {
    cat("File not found for subject:", subject_id, "\n")
  }
} 



transposed_data_trackon <- transpose(merged_data_trackon)
transposed_data_trackon_noname <- transposed_data_trackon[-1, ]  # Remove the first row
colnames(transposed_data_trackon) <- transposed_data_trackon[1, ]  # Use the first row as column names
transposed_data_trackon_name <- transposed_data_trackon[-1, ]  # Remove the first row

# Reset column names if needed (assuming the columns in all files are the same)
colnames(merged_data_trackon) <- c("HDID", colnames(merged_data_trackon)[-1])
ct_data_trackon <- read.csv("/Users/charlie/Desktop/my_projects/neurotransmitter/imaging_metrics/thickness/trackon/reprocessed_results/trackon_fsv7.csv")
merged_data_trackon <- merge(merged_data_trackon, ct_data_trackon[, c("HDID", "site", "sex", "cag", "group", "age1", "DBS_1", "eTIV")], by = "HDID", all.x = TRUE)
desired_order <- c("HDID", "site", "sex", "cag", "group", "age1", "DBS_1", "eTIV", names(merged_data_trackon)[!names(merged_data_trackon) %in% c("HDID", "site", "sex", "cag", "group", "age1", "burden1", "eTIV")])
# Reorder the columns in merged_data
merged_data_trackon <- merged_data_trackon[, desired_order]
merged_data_trackon <- merged_data_trackon %>% 
  rename(burden1 = DBS_1)








#####Track#####

##Obtaining degree
data_directory <- "/Users/charlie/Desktop/my_projects/neurotransmitter/imaging_metrics/thickness/track/carlos"


# List of subject names

subject_ids <- c(
   'XXX-XXX-XXX', 'YYY-YYY-YYY'
)


# Initialize an empty data frame to store the merged data
merged_data_track <- data.frame()

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
    merged_data_track <- rbind(merged_data_track, subject_avg)
  } else {
    cat("File not found for subject:", subject_id, "\n")
  }
} 



transposed_data_track <- transpose(merged_data_track)
transposed_data_track_noname <- transposed_data_track[-1, ]  # Remove the first row
colnames(transposed_data_track) <- transposed_data_track[1, ]  # Use the first row as column names
transposed_data_track_name <- transposed_data_track[-1, ]  # Remove the first row

# Reset column names if needed (assuming the columns in all files are the same)
colnames(merged_data_track) <- c("HDID", colnames(merged_data_track)[-1])
ct_data_track <- read.csv("/Users/charlie/Desktop/my_projects/neurotransmitter/imaging_metrics/thickness/track/reprocessed_results/track_fsv7.csv")
merged_data_track <- merge(merged_data_track, ct_data_track[, c("HDID", "site", "sex", "cag", "group", "age1", "burden1", "eTIV")], by = "HDID", all.x = TRUE)
desired_order <- c("HDID", "site", "sex", "cag", "group", "age1", "burden1", "eTIV", names(merged_data_track)[!names(merged_data_track) %in% c("HDID", "site", "sex", "cag", "group", "age1", "burden1", "eTIV")])
# Reorder the columns in merged_data
merged_data_track <- merged_data_track[, desired_order]



#####Merged demos#####

select_merged_data_yas <- merged_data_yas[, 1:8]
select_merged_data_trackon <- merged_data_trackon[, 1:8]
select_merged_data_track <- merged_data_track[, 1:8]


select_merged_data_trackon$site <- as.character(select_merged_data_trackon$site)
select_merged_data_trackon <- select_merged_data_trackon %>%
  mutate(site = case_when(
    site == "0" ~ "2",
    site == "1" ~ "3",
    site == "2" ~ "4",
    site == "3" ~ "5",
    TRUE ~ site  # Keep other values unchanged
  ))



select_merged_data_track$site <- as.character(select_merged_data_track$site)
select_merged_data_track <- select_merged_data_track %>%
  mutate(site = case_when(
    site == "0" ~ "6",
    site == "1" ~ "7",
    site == "2" ~ "8",
    site == "3" ~ "9",
    TRUE ~ site  # Keep other values unchanged
  ))

appended_demos <- rbind(select_merged_data_yas, select_merged_data_trackon, select_merged_data_track)
write.csv(appended_demos,file = '/Users/charlie/Desktop/my_projects/neurotransmitter/R/combat/results/demos.csv', row.names = FALSE)

transposed_appended_demos <- transpose(appended_demos)
colnames(transposed_appended_demos) <- transposed_appended_demos[1, ]  # Use the first row as column names
transposed_appended_demos <- transposed_appended_demos[-1, ] 
write.csv(appended_demos,file = '/Users/charlie/Desktop/my_projects/neurotransmitter/R/combat/results/demos_long.csv', row.names = FALSE)

demos_yas <- appended_demos[appended_demos$site %in% c(1), ]
write.csv(demos_yas,file = '/Users/charlie/Desktop/my_projects/neurotransmitter/R/combat/results/demos_long_yas.csv', row.names = FALSE)

demos_trackon <- appended_demos[appended_demos$site %in% c(2, 3, 4, 5), ]
write.csv(demos_trackon,file = '/Users/charlie/Desktop/my_projects/neurotransmitter/R/combat/results/demos_long_trackon.csv', row.names = FALSE)

demos_track <- appended_demos[appended_demos$site %in% c(6, 7, 8, 9), ]
write.csv(demos_track,file = '/Users/charlie/Desktop/my_projects/neurotransmitter/R/combat/results/demos_long_track.csv', row.names = FALSE)






write.csv(transposed_appended_demos,file = '/Users/charlie/Desktop/my_projects/neurotransmitter/R/combat/results/demos.csv', row.names = FALSE)


#####Merged MIND#####
appended_data_name <- cbind( transposed_data_yas_name,transposed_data_trackon_name,transposed_data_track_name)
write.csv(appended_data_name,file = '/Users/charlie/Desktop/my_projects/neurotransmitter/R/combat/results/mind_node_hdid.csv', row.names = FALSE)

write.csv(transposed_data_yas_name,file = '/Users/charlie/Desktop/my_projects/neurotransmitter/R/combat/results/mind_node_hdid_yas.csv', row.names = FALSE)
write.csv(transposed_data_trackon_name,file = '/Users/charlie/Desktop/my_projects/neurotransmitter/R/combat/results/mind_node_hdid_trackon.csv', row.names = FALSE)
write.csv(transposed_data_track_name,file = '/Users/charlie/Desktop/my_projects/neurotransmitter/R/combat/results/mind_node_hdid_track.csv', row.names = FALSE)





appended_data_noname <- cbind( transposed_data_yas_noname,transposed_data_trackon_noname,transposed_data_track_noname)

colnames(appended_data_noname) <- NULL
write.csv(appended_data_noname,file = '/Users/charlie/Desktop/my_projects/neurotransmitter/R/combat/results/mind_node.csv', row.names = FALSE)

write.csv(transposed_data_yas_noname,file = '/Users/charlie/Desktop/my_projects/neurotransmitter/R/combat/results/mind_node_yas.csv', row.names = FALSE)
write.csv(transposed_data_trackon_noname,file = '/Users/charlie/Desktop/my_projects/neurotransmitter/R/combat/results/mind_node_trackon.csv', row.names = FALSE)
write.csv(transposed_data_track_noname,file = '/Users/charlie/Desktop/my_projects/neurotransmitter/R/combat/results/mind_node_track.csv', row.names = FALSE)

