library(Rmagic)

# Command line arguments
args <- commandArgs(TRUE)
path_to_data <- args[1] # Log-normalized counts
t <- as.integer(args[2])
path_to_result <- args[3] # Imputed counts

# Read and prepare the data
dataset <- read.table(path_to_data, header = TRUE, row.names = 1)
dataset <- t(as.matrix(dataset))

# Run the imputation algorithm
result <- magic(dataset, genes = "all_genes", t = t)

# Get adjusted count matrix
imputed_dataset <- result$result
imputed_dataset <- t(imputed_dataset)

# Save adjusted count matrix
write.table(imputed_dataset, gzfile(paste0(path_to_result, ".magic_", t, ".txt.gz")), sep = "\t", quote = FALSE)
