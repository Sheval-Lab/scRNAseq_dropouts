library(furrr)
library(Seurat)

data_dir <- "data/2_k562_dataset/cellranger_outs"
results_dir <- "results/2_k562_analysis/gene_expression_filtered"


preprocessing <- function(path_to_data, path_to_result){
  # Get sample ID
  sample_name = basename(path_to_data)
  
  # Read in gene x cell matrix 
  dataset = Read10X(data.dir = file.path(path_to_data, "outs", "filtered_feature_bc_matrix"))
  
  # Perform non-expressed genes and poor-quality cells fltering
  dataset = CreateSeuratObject(counts = dataset, min.cells = 3, min.features = 200)
  
  # Log-normalized data 
  dataset = NormalizeData(dataset, normalization.method = "LogNormalize", scale.factor = 10000)
  
  # Save 
  saveRDS(dataset, file.path(path_to_result, paste0(sample_name, ".sobj.rds"))) # Seurat object
  write.table(as.matrix(GetAssayData(object = dataset, assay = "RNA", slot = "counts")), gzfile(file.path(path_to_result, paste0(sample_name, ".raw_counts.txt.gz"))), sep = "\t", quote = FALSE) # Un-normalized counts: used by kNN-smoothing, SAVER, scImpute, scVI
  write.table(as.matrix(GetAssayData(object = dataset, assay = "RNA", slot = "data")), gzfile(file.path(path_to_result, paste0(sample_name, ".norm_counts.txt.gz"))), sep = "\t", quote = FALSE) # Log-normalized counts: used by ALRA, MAGIC 
}


# Run preprocessing in parallel
plan(multisession, workers = 4)
future_walk(list.dirs(data_dir, recursive = FALSE), preprocessing, results_dir)

