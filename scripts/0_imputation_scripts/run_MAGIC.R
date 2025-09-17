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

# Re-normalize data so that total sum of counts is scaled up to 10000
imputed_dataset_unlog <- exp(imputed_dataset) - 1
total_sums <- colSums(imputed_dataset_unlog)
imputed_dataset_lognorm <- log1p(sweep(imputed_dataset_unlog, 2, total_sums, FUN = "/") * 10000)

# Save adjusted count matrix
write.table(imputed_dataset_lognorm, gzfile(paste0(path_to_result, ".magic_", t, ".txt.gz")), sep = "\t", quote = FALSE)
