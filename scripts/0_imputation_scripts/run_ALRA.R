library(ALRA)

# Command line arguments
args <- commandArgs(TRUE)
path_to_data <- args[1] # Log-normalized counts 
path_to_result <- args[2] # Imputed counts

# Read and prepare the data
dataset <- read.table(path_to_data, header = TRUE, row.names = 1)
dataset <- t(as.matrix(dataset))

# Run the imputation algorithm
result <- alra(dataset)

# Get adjusted count matrix
imputed_dataset <- result[[3]]
rownames(imputed_dataset) <- rownames(dataset)
imputed_dataset <- t(imputed_dataset)

# Save adjusted count matrix
write.table(imputed_dataset, gzfile(paste0(path_to_result, ".alra.txt.gz")), sep = "\t", quote = FALSE)
