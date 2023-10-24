
# Install the DESeq2 package
BiocManager::install("DESeq2",force = TRUE)
# Load the library.
library("DESeq2")
# The pathway to the file that contains the counts
counts_file <- "H:\\Project\\downsteam\\different_expression\\merged_kallisto_counts.csv"
# The pathway to the file that contains design infomations
design_file <- "H:\\Project\\downsteam\\different_expression\\design.csv"
# The final result file
output_file <- "H:\\Project\\downsteam\\different_expression\\result_VsaFsa.csv"
# Read the sample file
original_colData <- read.csv(design_file, stringsAsFactors = F)
# Extract some groups that want to focus on
filter_condition <- original_colData$condition %in% c("VSA", "FSA")
colData <- original_colData[filter_condition,]
# Turn conditions into factors.
colData$condition = factor(colData$condition)

# The first level should correspond to the first entry in the file!
# Required later when building a model.
colData$condition = relevel(colData$condition, toString(colData$condition[1]))

# Isolate the sample names.
sample_names <- colData$sample

# Read the data from the standard input.
original_df = read.csv(counts_file, header=TRUE, row.names=1 )
# Extract some columns that want to focus on
df = subset(original_df, select = c("KS_PF_1","KS_PF_2","KS_PF_3",
                                    "SH_PF_1","SH_PF_2","SH_PF_3"))

# Created rounded integers for the count data
countData = round(df[, sample_names])

# Other columns in the dataframe that are not sample information. 
otherCols = df[!(names(df) %in% sample_names)]
#
# Running DESeq2
#

# Create DESEq2 dataset.
dds = DESeqDataSetFromMatrix(countData=countData, colData=colData, design = ~condition)

# Run deseq
dse = DESeq(dds)

# Format the results.
res = results(dse)

#
# The rest of the code is about formatting the output dataframe.
#

# Turn the DESeq2 results into a data frame.
data = cbind(otherCols, data.frame(res))

# Create the foldChange column.
data$foldChange = 2 ^ data$log2FoldChange

# Rename columns to better reflect reality.
names(data)[names(data)=="pvalue"] <-"PValue"
names(data)[names(data)=="padj"] <- "FDR"

# Create a real adjusted pvalue
data$PAdj = p.adjust(data$PValue, method="hochberg")

# Sort the data by PValue to compute false discovery counts.
data = data[with(data, order(PValue, -foldChange)), ]

# Compute the false discovery counts on the sorted table.
data$falsePos = 1:nrow(data) * data$FDR

# Create the additional columns that we wish to present.
data$baseMeanA = 1
data$baseMeanB = 1


# Get the normalized counts.
normed = counts(dse, normalized=TRUE)

# Round normalized counts to a single digit.
normed = round(normed, 1)

# Merge the two datasets by row names.
total <- merge(data, normed, by=0)

# Sort again for output.
total = total[with(total, order(PValue, -foldChange)), ]

# Sample names for condition A:KS_PF
col_names_A = data.frame(split(colData, colData$condition)[1])[,1]

# Sample names for condition B:SH_PF
col_names_B = data.frame(split(colData, colData$condition)[2])[,1]


# Create the individual baseMean columns.
total$baseMeanA = rowMeans(total[, col_names_A])
total$baseMeanB = rowMeans(total[, col_names_B])


# Bringing some sanity to numbers. Round columns to fewer digits.
total$foldChange = round(total$foldChange, 3)
total$log2FoldChange = round(total$log2FoldChange, 1)
total$baseMean  = round(total$baseMean, 1)
total$baseMeanA = round(total$baseMeanA, 1)
total$baseMeanB =  round(total$baseMeanB, 1)
total$lfcSE = round(total$lfcSE, 2)
total$stat = round(total$stat, 2)
total$FDR = round(total$FDR, 4)
total$falsePos = round(total$falsePos, 0)

# Reformat these columns as string.
total$PAdj = formatC(total$PAdj, format = "e", digits = 1)
total$PValue = formatC(total$PValue, format = "e", digits = 1)

# Rename the first column.
colnames(total)[1] <- "name"

# Reorganize columns names to make more sense.
new_cols = c("name", names(otherCols), "baseMean","baseMeanA","baseMeanB","foldChange",
             "log2FoldChange","lfcSE","stat","PValue","PAdj", "FDR","falsePos",col_names_A, col_names_B)

# Slice the dataframe with new columns.
total = total[, new_cols]
names(total)[names(total)=="baseMeanA"] <- "VSA"
names(total)[names(total)=="baseMeanB"] <- "FSA"
total_VsaFsa <- total

# Remove the specific row whose "baseMean" is 0
total_VsaFsa_Remove0 = subset(total_VsaFsa, total_VsaFsa$baseMean != 0)


# Write the results to the standard output.
write.csv(total_VsaFsa_Remove0,
          file="H:\\Project\\downsteam\\different_expression\\result_VsaFsa_remove0.csv",
          row.names=FALSE, quote=FALSE)



###Acknowlegement:
###This R code is based on Website "https://www.biostarhandbook.com/books/rnaseq/differential-expression.html"
###and I make some adjustments to suit for my data.
