#
# Create heat map from a differential expression count table.
#
# Load the library.
#install.packages("gplots")
library(gplots)

# The name of the file that contains the counts.
count_file = "H:\\Project\\downsteam\\different_expression\\result_VsaFsa_remove0.csv"

# The name of the output file.
output_file = "H:\\Project\\downsteam\\different_expression\\result_VsaFsa_remove0_heatmap.pdf"

# Inform the user.
print("# Tool: Create Heatmap ")
print(paste("# Input: ", count_file))
print(paste("# Output: ", output_file))

# FDR cutoff.
MIN_FDR = 0.05

# Plot width
WIDTH = 12

# Plot height.
HEIGHT = 13

# Set the margins
MARGINS = c(9, 12)

# Relative heights of the rows in the plot.
LHEI = c(1, 5)

# Read normalized counts from the standard input.
data = read.csv(count_file, header=T, as.is=TRUE)

# Subset data for values under a treshold.
data = subset(data, data$FDR <= MIN_FDR)

# The heatmap row names will be taken from the first column.
row_names = data[, 1]

# The code assumes that the normalized data matrix is listed to the right of the falsePos column.
idx = which(colnames(data) == "falsePos") + 1

# The normalized counts are on the right size.
counts = data[, idx : ncol(data)]

# Load the data from the second column on.
values = as.matrix(counts)

# Adds a little noise to each element to avoid the
# clustering function failure on zero variance rows.
values = jitter(values, factor = 1, amount = 0.00001)

# Normalize each row to a z-score
zscores = NULL
for (i in 1 : nrow(values)) {
  row = values[i,]
  zrow = (row - mean(row)) / sd(row)
  zscores = rbind(zscores, zrow)
}

# Set the row names on the zscores.
row.names(zscores) = row_names

# Turn the data into a matrix for heatmap2.
zscores = as.matrix(zscores)

# Set the color palette.
col = greenred



heatmap.2(zscores, 
               col=greenred, 
               density.info="none", Colv=NULL,
    dendrogram="row", trace="none", margins=c(9,12), lhei=c(1,5))
# Create a PDF device to save the plot
pdf(output_file, width = 12, height = 13)

#dev.off()



### Citation&&Acknowledgement:
### This R code is based on the Website of "https://www.biostarhandbook.com/books/rnaseq/differential-expression.html"
### and I make some adjustments to suit for my data.
