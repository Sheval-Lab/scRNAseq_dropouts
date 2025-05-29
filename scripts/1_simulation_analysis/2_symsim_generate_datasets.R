library(tidyverse)
library(SymSim)

data_dir <- "data/1_simulated_datasets/symsim"


# Load estimated best parameters for simulations -------------------------------
params_files <- list.files(data_dir, pattern = "best_params.txt$", full.names = TRUE)
params_df <- purrr::map(params_files, read_tsv, id = "celltype") %>% 
  list_rbind() %>% 
  mutate(celltype = str_split_i(basename(celltype), pattern = "_", i = 2)) %>% 
  group_by(celltype) %>% 
  dplyr::mutate(n_param_set = 1:n(), .after = celltype) %>% 
  ungroup()


# Select one of the parameters set for each cell type --------------------------
## Use Q-Q plots to make the choice
selected_params_df <- params_df %>% 
  mutate(celltype_n_set = paste(celltype, n_param_set, sep = "_")) %>% 
  # basal: 2nd best parameters set
  # club: 2nd best parameters set
  # ciliated: 2nd best parameters set
  filter(celltype_n_set %in% c("basal_2", "club_3", "ciliated_2"))


# Generate synthetic datasets --------------------------------------------------
run_symsim <- function(cell_type, params_table) {
  # Select parameters obtained for specified cell type population
  params = params_table %>% filter(celltype == cell_type)
  

  # Generate true expression matrix
  true_counts = SimulateTrueCounts(
    ncells_total = params$cell_number, 
    ngenes = params$gene_number, 
    evf_type = "one.population",
    Sigma = params$Sigma, 
    randseed = 42,
    gene_effect_prob = params$gene_effect_prob, 
    gene_effects_sd = params$gene_effects_sd,
    scale_s = params$scale_s, 
    mean_hge = params$mean_hge,
    prop_hge = params$prop_hge)
  
  # Generate observed expression matrix
  data(gene_len_pool)
  gene_lengths = sample(gene_len_pool, params$gene_number, replace = FALSE)
  obs_counts = True2ObservedCounts(
    true_counts = true_counts[[1]],
    meta_cell = true_counts[[3]], 
    protocol = "UMI",
    alpha_mean = params$alpha_mean, 
    alpha_sd = params$alpha_sd,
    gene_len = gene_lengths,
    depth_mean = params$depth_mean, 
    depth_sd = params$depth_sd,
    nPCR1 = params$nPCR1)
  
  # Save simulated true and observed expression matrices
  write.table(true_counts[[1]], gzfile(file.path(data_dir, paste0("symsim_", cell_type, "_true_counts.txt.gz"))), sep = "\t", quote = FALSE)
  write.table(obs_counts[[1]], gzfile(file.path(data_dir, paste0("symsim_", cell_type, "_obs_counts.txt.gz"))), sep = "\t", quote = FALSE)
}


## Run generation of datasets... -----------------------------------------------
walk(selected_params_df$celltype, run_symsim, selected_params_df)
