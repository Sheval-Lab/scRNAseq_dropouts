library(SAVER)
library(Seurat)

# Command line arguments
args <- commandArgs(TRUE)
path_to_data <- args[1] # Raw un-normalized counts  
path_to_result <- args[2] # Imputed counts

# Read the data
dataset <- read.table(path_to_data, header = TRUE, row.names = 1)

# Run the imputation algorithm
result <- saver(dataset, ncores = 20, estimates.only = TRUE)

# Re-normalize data so that total sum of counts is scaled up to 10000
total_sums <- colSums(result)
result_norm <- log1p(sweep(result, 2, total_sums, FUN = "/") * 10000)

# Log-transform adjusted count matrix
imputed_dataset <- log1p(result_norm) 

# Save adjusted count matrix
write.table(imputed_dataset, gzfile(paste0(path_to_result, ".saver.txt.gz")), sep = "\t", quote = FALSE)
