library(Seurat)

# Command line arguments
args <- commandArgs(TRUE)
res_name_prefix <- args[1] # Imputed results
dca_model <- args[2] # Shortened name of DCA model

# Read in the data
imputed_unnorm_dataset <- read.table(paste0(res_name_prefix, ".dca_", dca_model, "_res_dir/mean.tsv"), header = TRUE, row.names=1)

# Log-normalize adjusted count matrix
imputed_dataset <- CreateSeuratObject(counts = imputed_unnorm_dataset)
imputed_dataset <- NormalizeData(imputed_dataset, normalization.method = "LogNormalize", scale.factor = 10000)

# Save adjusted count matrix
write.table(as.matrix(GetAssayData(object = imputed_dataset, assay = "RNA", slot = "data")), gzfile(paste0(res_name_prefix, ".dca_", dca_model, ".txt.gz")), sep = "\t", quote = FALSE)

