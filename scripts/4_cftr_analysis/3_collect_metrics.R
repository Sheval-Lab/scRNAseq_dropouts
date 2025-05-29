library(tidyverse)
source("scripts/functions.R")

data_dir <- "results/4_cftr_analysis"
results_dir <- "results/4_cftr_analysis/metrics"


# Collect data on the fraction of cells with non-zero expression of a gene and mean normalized expression of a gene ----
## Get list of files with normalized counts matrices
raw_counts_files <- list.files(file.path(data_dir, "gene_expression_filtered"), pattern = "norm.txt.gz$", full.names = TRUE)
imputed_counts_files <- list.files(file.path(data_dir, "gene_expression_imputed"), pattern = ".txt.gz$", full.names = TRUE)

## Calculate and aggregate metrics tables
metrics <- collect_metrics(c(raw_counts_files, imputed_counts_files), replace_na = FALSE)
dataset_size <- metrics$dataset_size
dropout_rate <- metrics$dropout_rate
nonzero_fraction <- metrics$nonzero_fraction
avg_norm_counts <- metrics$avg_norm_counts


## Save metrics into files -----------------------------------------------------
write_tsv(dataset_size, file.path(results_dir, "dataset_size.txt"))
write_tsv(dropout_rate, file.path(results_dir, "dropout_rate.txt"))
write_tsv(nonzero_fraction, file.path(results_dir, "nonzero_fraction.txt"))
write_tsv(avg_norm_counts, file.path(results_dir, "avg_norm_counts.txt"))


# Calculate nonzero fraction of CFTR for specific cell populations -------------
calculate_metrics_in_unsplit_ds <- function(path_to_matrix, genes_upper, groupping_vars, metadata_df) {
  # Get imputation run ID
  run_id = str_remove(basename(path_to_matrix), ".txt.gz")
  
  # Read in counts matrix
  data = read.table(path_to_matrix, header = TRUE, row.names = 1)
  
  # Subset gene expression matrix by gene of interest, and add cells metadata
  data_flt = data %>% 
    rownames_to_column("GeneID") %>% 
    filter(str_detect(str_to_upper(GeneID), genes_upper)) %>% 
    pivot_longer(where(is.numeric), names_to = "BarcodeID", values_to = "count") %>% 
    left_join(metadata_df, by = "BarcodeID") 
  
  # For each combination of airway section and cell type calculate nonzero fraction 
  data_flt %>% 
    group_by(across(all_of(c(groupping_vars, "GeneID")))) %>% 
    summarise(
      avg_count = mean(count),
      nonzero_fraction = mean(count > 0)) %>% 
    ungroup() %>% 
    mutate(run = run_id)
}


### Select files - imputation results for mixed cell type populations
okuda_imputed_counts_files_subset <- imputed_counts_files %>% str_subset("all_")

### Cells metadata
okuda_cells_meta <- read_tsv(file.path(data_dir, "gene_expression_filtered", "okuda2021_cells_metadata.txt"))

#### Replace '-' with '.' in BarcodeID because read.table will change colnames that way...
okuda_cells_meta <- okuda_cells_meta %>% mutate(BarcodeID = str_replace(BarcodeID, "-", "."))

okuda_metrics_in_unsplit_ds <- map(okuda_imputed_counts_files_subset, calculate_metrics_in_unsplit_ds, genes_upper = "CFTR", groupping_vars = c("dataset", "Airway", "CellGroup"), okuda_cells_meta) %>% 
  list_rbind()

write_tsv(okuda_metrics_in_unsplit_ds, file.path(results_dir, "okuda2021_metrics_in_unsplit_ds.txt"))

