library(tidyverse)
library(Seurat)
source("scripts/functions.R")

data_dir <- "data/4_okuda2021"
results_dir <- "results/4_cftr_analysis/gene_expression_filtered"

# Okuda2021 datasets -----------------------------------------------------------
## Load Seurat RDS object with all cell x gene expression ----------------------
### Drop-seq dataset
ds_sobj <- readRDS(file.path(data_dir, "GSE160673_Ken_7scRNAseq_mt20_merged16labeled_seurat.rds"))

#### Rename Ionocyte/Neuroendocrine cells into simple name
ds_sobj$CellGroup <- if_else(ds_sobj$CellGroup == "Ionocyte/NE", "Ionocyte", ds_sobj$CellGroup) %>% as.vector()

### 10x Chromium dataset
ch_sobj <- readRDS(file.path(data_dir, "GSE160664_Cedars_all_seurat.rds"))

#### Aggregate fine cell types into larger groups (as in ds_sobj)
ch_sobj$CellGroup <- case_when(
  str_detect(str_to_lower(ch_sobj$predicted.id), "basal") ~ "Basal",
  str_detect(ch_sobj$predicted.id, "Ciliated") ~ "Ciliated",
  str_detect(ch_sobj$predicted.id, "Cycling") ~ "Cycling",
  str_detect(ch_sobj$predicted.id, "Ionocyte") ~ "Ionocyte",
  str_detect(ch_sobj$predicted.id, "Secretory") ~ "Secretory",
  str_detect(ch_sobj$predicted.id, "Undefined") ~ "Undefined",
  TRUE ~ "Others") 


## Datasets basic QC metrics ---------------------------------------------------
okuda_datasets_qc <- bind_rows(
  collect_ds_metrics(ds_sobj, "Drop-seq"),
  collect_ds_metrics(ch_sobj, "10x Chromium"))

write_tsv(okuda_datasets_qc, file.path(results_dir, "okuda2021_datasets_qc.txt"))


## Create metadata table -------------------------------------------------------
ds_meta <- ds_sobj@meta.data %>% 
  rownames_to_column("BarcodeID") %>% 
  select(BarcodeID, Airway, CellGroup)

ch_meta <- ch_sobj@meta.data %>% 
  rownames_to_column("BarcodeID") %>% 
  select(BarcodeID, Airway, CellGroup)

bind_rows(
  mutate(ds_meta, dataset = "Drop-seq"),
  mutate(ch_meta, dataset = "10x Chromium")) %>% 
  write_tsv(file.path(results_dir, "okuda2021_cells_metadata.txt"))


## Inspect CFTR expression levels ----------------------------------------------
ds_cftr <- extract_gene_expr_q(ds_sobj, groupping_vars = c("Airway", "CellGroup"), genes = str_to_upper(c("CFTR")))
ch_cftr <- extract_gene_expr_q(ch_sobj, groupping_vars = c("Airway", "CellGroup"), genes = str_to_upper(c("CFTR")))

bind_rows(
  mutate(ds_cftr, dataset = "Drop-seq"),
  mutate(ch_cftr, dataset = "10x Chromium")) %>% 
  write_tsv(file.path(results_dir, "okuda2021_cftr_expression_groups.txt"))


## Extract cell populations and save the data ----------------------------------
save_data_from_sobj <- function(sobj, airway = "both", celltype = "all", dataset, path_to_result) {
  # Subset airway section or cell type (or both)
  if (airway != "both") { sobj <- subset(sobj, subset = Airway == airway) }
  if (celltype != "all") { sobj <- subset(sobj, subset = CellGroup == celltype) }
  
  # Determine dataset name
  dataset_name = paste(celltype, airway, dataset, sep = "_")
  
  # Use fixed threshold to filter out genes with undetected expression
  counts = GetAssayData(sobj, assay="RNA", slot="counts")   
  detected_genes = rowSums(counts>0 ) >= 3
  
  write.table(as.matrix(GetAssayData(object = sobj, assay = "RNA", slot = "counts"))[names(detected_genes), ], gzfile(file.path(path_to_result, paste0(dataset_name, ".raw.txt.gz"))), sep = "\t", quote = FALSE) # Un-normalized counts: used by kNN-smoothing, SAVER, scImpute, scVI
  write.table(as.matrix(GetAssayData(object = sobj, assay = "RNA", slot = "data"))[names(detected_genes), ], gzfile(file.path(path_to_result, paste0(dataset_name, ".norm.txt.gz"))), sep = "\t", quote = FALSE) # Log-normalized counts: used by ALRA, MAGIC 
}


### Data frame with all options to walk over
okuda_datasets <- data.frame(
  airway = rep(rep(c("both", "LAE", "SAE"), each = 5), times = 2),
  celltype = rep(c("all", "Basal", "Ciliated", "Ionocyte", "Secretory"), times = 6),
  dataset = rep(c("dropseq", "10x"), each = 15))

### Drop-seq dataset
okuda_datasets %>% 
  filter(dataset == "dropseq") %>% 
  pwalk(save_data_from_sobj, sobj = ds_sobj, path_to_result = results_dir)

### 10x Chromium dataset
okuda_datasets %>% 
  filter(dataset == "10x") %>% 
  pwalk(save_data_from_sobj, sobj = ch_sobj, path_to_result = results_dir)


