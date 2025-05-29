library(tidyverse)
source("scripts/functions.R")

data_dir <- "results/2_k562_analysis/gene_expression_imputed"
results_dir <- "results/2_k562_analysis/imputed_metrics"


# Collect data on the fraction of cells with non-zero expression of a gene and mean normalized expression of a gene ----
## Get list of files with normalized counts matrices
imputed_counts_files <- list.files(data_dir, pattern = ".txt.gz$", full.names = TRUE)

## Calculate and aggregate metrics tables
imputed_metrics <- collect_metrics(imputed_counts_files)
imputed_nonzero_fraction <- imputed_metrics$nonzero_fraction
imputed_avg_norm_counts <- imputed_metrics$avg_norm_counts


## Save metrics into files -----------------------------------------------------
write_tsv(imputed_nonzero_fraction, gzfile(file.path(results_dir, "imputed_nonzero_fraction.txt.gz")))
write_tsv(imputed_avg_norm_counts, gzfile(file.path(results_dir, "imputed_avg_norm_counts.txt.gz")))


