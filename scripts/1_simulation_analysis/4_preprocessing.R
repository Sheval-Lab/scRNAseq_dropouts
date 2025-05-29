library(tidyverse)
library(furrr)
library(Seurat)

data_dir <- "data/1_simulated_datasets"
results_dir <- "results/1_simulation_analysis/gene_expression_filtered"


preprocessing <- function(path_to_data, path_to_result){
  # Get sample ID
  sample_name = str_remove(basename(path_to_data), ".txt.gz")
  
  # Read in gene x cell matrix 
  dataset = read.table(path_to_data, header = TRUE)
  rownames(dataset) = paste0("Gene", 1:nrow(dataset))
  colnames(dataset) = paste0("Cell", 1:ncol(dataset))
  
  # Perform non-expressed genes and poor-quality cells fltering
  dataset = CreateSeuratObject(counts = dataset, min.cells = 3, min.features = 200)
  
  # Log-normalized data 
  dataset = NormalizeData(dataset, normalization.method = "LogNormalize", scale.factor = 10000)
  
  # Save 
  saveRDS(dataset, file.path(path_to_result, paste0(sample_name, ".sobj.rds"))) # Seurat object
  write.table(as.matrix(GetAssayData(object = dataset, assay = "RNA", slot = "counts")), gzfile(file.path(path_to_result, paste0(sample_name, ".raw.txt.gz"))), sep = "\t", quote = FALSE) # Un-normalized counts: used by kNN-smoothing, SAVER, scImpute, scVI
  write.table(as.matrix(GetAssayData(object = dataset, assay = "RNA", slot = "data")), gzfile(file.path(path_to_result, paste0(sample_name, ".norm.txt.gz"))), sep = "\t", quote = FALSE) # Log-normalized counts: used by ALRA, MAGIC 
}


# Run preprocessing in parallel
plan(multisession, workers = 4)
future_walk(list.files(data_dir, pattern = "counts.txt.gz", full.names = TRUE, recursive = TRUE), preprocessing, results_dir)
