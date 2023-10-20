install.packages("ghql")
library("ghql")
### set folder path to location of your kallisto outputs
folder_path <- "H:/Project/downsteam/kallisto_counts/"
folders <- list.dirs(path = folder_path, recursive=TRUE, full.names = FALSE)
# remove the root folder name
folders <- folders[- 1]
# check
print(folders)

count_matrix <- NULL
for (folder in folders) {
  # specify the full file path to the abundance data
  file_path <- paste0("H:/Project/downsteam/kallisto_counts/", folder, "/abundance.tsv")
  
  # read in the abundance data
  data <- read.table(file_path, header = TRUE, sep = "\t")
  
  # add the abundance data to the count matrix
  count_matrix <- cbind(count_matrix, data[, "est_counts"])
}
print(count_matrix)
rownames(count_matrix) <- data[, "target_id"]
colnames(count_matrix) <- gsub(".*/", " ", folders)
print(count_matrix)
write.csv(count_matrix, file = "H:/Project/downsteam/kallisto_counts/merged_kallisto_counts.csv", row.names = TRUE)

### Acknowledgements:
### These R codes are based on the  website: "https://zbengt.github.io/2023-03-09-Mergin-Kallisto_Abundance/"
### and make some adjustments to make it more suitable for my data. (line 5-7)
