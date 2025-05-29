calculate_metrics <- function(path_to_matrix) {
  # Get sample ID
  sample_name = str_split_i(basename(path_to_matrix), pattern = "\\.", i = 1)
  
  # Get method name
  method = basename(path_to_matrix) %>% 
    # Remove file format
    str_remove(".txt.gz") %>%
    # Remove sample ID
    str_remove_all(paste0(sample_name, "."))
  
  # Read in counts matrix
  data = read.table(path_to_matrix, header = TRUE, row.names = 1)
  
  # Save dataset size
  dataset_size = data.frame(
    sample = sample_name,
    method = method,
    n_genes = nrow(data),
    n_cells = ncol(data))
  
  # Calculate dropout rate 
  dropout_rate = data.frame(
    sample = sample_name,
    method = method,
    dropout_rate = mean(data == 0))
  
  # Calculate fraction of cells with non-zero expression of a gene
  nonzero_fraction = as.data.frame(rowMeans(data > 0)) 
  
  # Calculate mean normalized expression
  avg_norm_counts = as.data.frame(rowMeans(data))
  
  # Add sample ID as column name
  colnames(nonzero_fraction) = paste(sample_name, method, sep = ".")
  colnames(avg_norm_counts) = paste(sample_name, method, sep = ".")
  
  # Add a column with gene names
  nonzero_fraction = rownames_to_column(nonzero_fraction, "GeneID")
  avg_norm_counts = rownames_to_column(avg_norm_counts, "GeneID")
  
  return(list(dataset_size, dropout_rate, nonzero_fraction, avg_norm_counts))
}

collect_metrics <- function(counts_files, replace_na = TRUE) {
  # Calculate metrics and store them in lists
  metrics_list = map(counts_files, calculate_metrics) %>% list_transpose()
  
  # Separate lists to get dataframes with dataset size, dropout rates, nonzero fraction, and average expression metrics 
  dataset_size = metrics_list[[1]]
  dropout_rate = metrics_list[[2]]
  nonzero_fraction = metrics_list[[3]] %>% reduce(full_join, by = "GeneID")
  avg_norm_counts = metrics_list[[4]] %>% reduce(full_join, by = "GeneID")
  
  # Replace NA's with zeros 
  if (replace_na == TRUE) {
    nonzero_fraction = nonzero_fraction %>% mutate(across(where(is.numeric), \(x) replace_na(x, 0))) 
    avg_norm_counts = avg_norm_counts %>% mutate(across(where(is.numeric), \(x) replace_na(x, 0))) 
  }
  
  return(list("dataset_size" = dataset_size, "dropout_rate" = dropout_rate, "nonzero_fraction" = nonzero_fraction, "avg_norm_counts" = avg_norm_counts))
}


collect_ds_metrics <- function(sobj, dataset_name) {
  sobj@meta.data %>% 
    select(nCount_RNA, nFeature_RNA) %>% 
    pivot_longer(everything(), names_to = "QC") %>% 
    group_by(QC) %>% 
    summarise(
      min_value = min(value),
      mean_value = mean(value),
      median_value = median(value),
      max_value = max(value)) %>% 
    mutate(dataset = dataset_name)
}


# Calculate nonzero fraction and Q-group for a set of genes
extract_gene_expr_q <- function(sobj, groupping_vars, genes) {
  # Extract all genes expression values from Seurat object
  genes_expression = FetchData(sobj, vars = c(groupping_vars, rownames(sobj)), slot = "data")
  
  # Annotate gene expression quartiles
  gene_q = genes_expression %>% 
    pivot_longer(where(is.numeric), names_to = "GeneID", values_to = "count") %>% 
    group_by(across(all_of(c(groupping_vars, "GeneID")))) %>% 
    mutate(n_cells = sum(count > 0)) %>% 
    filter(n_cells >= 3) %>% 
    summarise(
      n_cells = n(),
      avg_count = mean(count),
      nonzero_fraction = mean(count > 0)) %>% 
    mutate(gene_quartile = as.factor(paste0("Q", ntile(avg_count, 4)))) %>% 
    ungroup() %>% 
    filter(str_to_upper(GeneID) %in% genes)
  
  gene_q
}
