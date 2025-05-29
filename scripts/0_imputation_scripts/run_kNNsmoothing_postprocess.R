library(Seurat)

# Command line arguments
args <- commandArgs(TRUE)
path_to_data <- args[1] # Imputed un-normalized counts  
k <- args[2] # The number of neighbors used for smoothing
path_to_result <- args[3] # Imputed log-normalized counts

# Read in the data
imputed_unnorm_dataset <- read.table(path_to_data, header = TRUE, row.names=1)

# Log-normalize adjusted count matrix
imputed_dataset <- CreateSeuratObject(counts = imputed_unnorm_dataset)
imputed_dataset <- NormalizeData(imputed_dataset, normalization.method = "LogNormalize", scale.factor = 10000)

# Save adjusted count matrix
write.table(as.matrix(GetAssayData(object = imputed_dataset, assay = "RNA", slot = "data")), gzfile(paste0(path_to_result, ".knn_", k, ".txt.gz")), sep = "\t", quote = FALSE)
