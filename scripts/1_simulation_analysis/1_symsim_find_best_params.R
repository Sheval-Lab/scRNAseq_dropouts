library(tidyverse)
library(SymSim)

data_dir <- "data/3_montoro2018"
results_dir <- "data/1_simulated_datasets/symsim"


# Load reference UMI scRNA-seq data -  from Montoro et al., 2018 ---------------
montoro_umi <- read.table(file.path(data_dir, "GSE103354_Trachea_droplet_UMIcounts.txt.gz"), header = TRUE, row.names = 1)


# Specify cell types to choose and sizes of cell populations -------------------
celltypes <- c("basal", "club", "ciliated")
population_sizes <- c(3000, 1000, 300)


# Search for optimal simulation parameters -------------------------------------
find_params <- function(cell_type, cell_number, ref_data) {
  set.seed(42)
  
  # Prepare reference data
  sampled_data = ref_data %>% 
    # Subset reference scRNA-seq data by chosen cell type
    dplyr::select(contains(str_to_title(cell_type))) %>% 
    # Select desired number of cells
    dplyr::select(sample(1:ncol(.), cell_number)) %>% 
    # Filter out all-zero genes
    filter(rowSums(.) !=0 ) %>% 
    as.matrix()
  
  # Calculate fraction of cells with non-zero expression of a gene
  nonzero_fraction = data.frame(
    GeneID = paste0("Gene", 1:nrow(sampled_data)),
    nonzero_fraction = rowMeans(sampled_data > 0)) 
  write_tsv(nonzero_fraction, file.path(data_dir, paste0("reference_nonzero_fraction.", cell_type, ".txt")))
  
  # Estimate parameters
  best_params = BestMatchParams("UMI", sampled_data, file.path(results_dir, paste0("symsim_", cell_type, "_best_params.qqplot")), n_optimal=5) # depth_range = c(20e3,500e3)

  # Add genes and cells numbers to best_params table
  best_params$gene_number = nrow(sampled_data)
  best_params$cell_number = ncol(sampled_data)

  # Save estimated parameters
  write_tsv(best_params, file.path(results_dir, paste0("symsim_", cell_type, "_best_params.txt")))
}


## Run generation of datasets... -----------------------------------------------
walk2(celltypes, population_sizes, find_params, montoro_umi)

