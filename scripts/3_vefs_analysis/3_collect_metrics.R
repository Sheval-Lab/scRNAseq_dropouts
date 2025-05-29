library(tidyverse)
source("scripts/functions.R")

data_dir <- "results/3_vefs_analysis"
results_dir <- "results/3_vefs_analysis/metrics"


# Collect data on the fraction of cells with non-zero expression of a gene and mean normalized expression of a gene ----
## Get list of files with normalized counts matrices
ref_counts_files <- file.path("data", "3_montoro2018", "GSE103354_Trachea_fullLength_TPM.txt.gz") # path to reference Smart-seq2 dataset
raw_counts_files <- list.files(file.path(data_dir, "gene_expression_filtered"), pattern = "norm.txt.gz$", full.names = TRUE)
imputed_counts_files <- list.files(file.path(data_dir, "gene_expression_imputed"), pattern = ".txt.gz$", full.names = TRUE)

## Calculate and aggregate metrics tables
metrics <- collect_metrics(c(ref_counts_files, raw_counts_files, imputed_counts_files), replace_na = FALSE)
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
## Select files - imputation results for mixed cell type populations
counts_files_subset <- c(ref_counts_files, str_subset(c(raw_counts_files, imputed_counts_files), "Total"))

calculate_metrics_in_unsplit_ds <- function(path_to_matrix) {
  # Get imputation run ID
  run_id = str_remove(basename(path_to_matrix), ".txt.gz")

  # Read in counts matrix
  data = read.table(path_to_matrix, header = TRUE, row.names = 1)

  # Subset gene expression matrix by gene of interest, and add cells metadata
  data_flt = data[c("Ace2", "Tmprss2", "Furin", "Ctsl"), ] %>%
    rownames_to_column("GeneID") %>%
    pivot_longer(where(is.numeric), names_to = "BarcodeID", values_to = "count") %>%
    mutate(CellGroup = str_extract(BarcodeID, "[A-Za-z]+$")) %>% 
    select(-BarcodeID)

  # For each cell type calculate nonzero fraction
  data_flt %>%
    group_by(CellGroup, GeneID) %>%
    summarise(
      avg_count = mean(count),
      nonzero_fraction = mean(count > 0)) %>%
    ungroup() %>%
    mutate(run = run_id)
}


metrics_in_unsplit_ds <- map(counts_files_subset, calculate_metrics_in_unsplit_ds) %>%
  list_rbind()

write_tsv(metrics_in_unsplit_ds, file.path(results_dir, "montoro2018_metrics_in_unsplit_ds.txt"))

