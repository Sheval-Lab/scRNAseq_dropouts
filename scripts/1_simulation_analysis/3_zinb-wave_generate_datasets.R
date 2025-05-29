library(tidyverse)
library(splatter)
# library(MASS)

data_dir <- "data/3_montoro2018"
results_dir <- "data/1_simulated_datasets/zinb-wave"


# Load reference UMI scRNA-seq data -  from Montoro et al., 2018 ---------------
montoro_umi <- read.table(file.path(data_dir, "GSE103354_Trachea_droplet_UMIcounts.txt.gz"), header = TRUE, row.names = 1)


# Specify cell types to choose and sizes of cell populations -------------------
celltypes <- c("basal", "club", "ciliated")
population_sizes <- c(3000, 1000, 300)


# Generate synthetic datasets --------------------------------------------------
run_zinbwave <- function(cell_type, cell_number, ref_data) {
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
  
  # Estimate parameters
  params = zinbEstimate(sampled_data)
  
  # Generate synthetic data
  synthetic_data = zinbSimulate(params = params, sparsify = TRUE, verbose = TRUE)
  
  # Save simulated true and observed expression matrices
  MASS::write.matrix(assay(synthetic_data, "TrueCounts"), gzfile(file.path(results_dir, paste0("zinb-wave_", cell_type, "_true_counts.txt.gz"))), sep="\t")
  MASS::write.matrix(assay(synthetic_data, "counts"), gzfile(file.path(results_dir, paste0("zinb-wave_", cell_type, "_obs_counts.txt.gz"))), sep="\t")
}


## Run generation of datasets... -----------------------------------------------
walk2(celltypes, population_sizes, run_zinbwave, montoro_umi)
