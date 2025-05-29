library(scImpute)
library(Seurat)

# Command line arguments
args <- commandArgs(TRUE)
path_to_data <- args[1] # Raw un-normalized counts  
path_to_intermediate_files <- args[2] # Full path of the output directory, which is used to store all intermdediate and final outputs
t <- args[3] # Threshold set on dropout probability
path_to_result <- args[4] # Imputed counts

# Run the imputation algorithm
scimpute(
  count_path = path_to_data,               # full path to raw count matrix
  infile     = "txt",                      # format of input file
  outfile    = "txt",                      # format of output file
  out_dir    = path_to_intermediate_files, # full path to output directory
  labeled    = FALSE,                      # cell type labels not available
  drop_thre  = t,                          # threshold set on dropout probability
  Kcluster   = 1,                          # 1 cell subpopulation
  ncores     = 20)                         # number of cores used in parallel computation

# Log-normalize adjusted count matrix
imputed_unnorm_dataset <- read.table(file.path(path_to_intermediate_files, "scimpute_count.txt"), header = TRUE, row.names=1)
imputed_dataset <- CreateSeuratObject(counts = imputed_unnorm_dataset)
imputed_dataset <- NormalizeData(imputed_dataset, normalization.method = "LogNormalize", scale.factor = 10000)

# Save adjusted count matrix
write.table(as.matrix(GetAssayData(object = imputed_dataset, assay = "RNA", slot = "data")), gzfile(paste0(path_to_result, ".scimpute_", t, ".txt.gz")), sep = "\t", quote = FALSE)

