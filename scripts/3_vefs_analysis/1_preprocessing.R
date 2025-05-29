library(tidyverse)
library(Seurat)
source("scripts/functions.R")

data_dir <- "data/3_montoro2018"
results_dir <- "results/3_vefs_analysis/gene_expression_filtered"


# Load scRNA-seq data (Montoro et al., 2018) -----------------------------------
## 10x Chromium dataset
montoro_umi <- read.table(file.path(data_dir, "GSE103354_Trachea_droplet_UMIcounts.txt.gz"), header = TRUE, row.names = 1)

## Smart-seq2 dataset
montoro_fl <- read.table(file.path(data_dir, "GSE103354_Trachea_fullLength_TPM.txt.gz"), header = TRUE, row.names = 1)


# Dataset basic QC metrics -----------------------------------------------------
umi_sobj <- CreateSeuratObject(counts = montoro_umi, min.cells = 3, min.features = 200)
fl_sobj <- CreateSeuratObject(counts = montoro_fl, min.cells = 3, min.features = 200)

umi_sobj$CellGroup <- Cells(umi_sobj) %>% str_split_i("_", 3)
umi_sobj$dataset <- "10x Chromium"
umi_sobj$percent.mt <- PercentageFeatureSet(umi_sobj, pattern = "^mt-")

fl_sobj$CellGroup <- Cells(fl_sobj) %>% str_split_i("_", 4)
fl_sobj$dataset <- "Smart-seq2"
fl_sobj$percent.mt <- PercentageFeatureSet(fl_sobj, pattern = "^mt-")

datasets_qc <- bind_rows(
  collect_ds_metrics(umi_sobj, "10x Chromium"),
  collect_ds_metrics(fl_sobj, "Smart-seq2"))

write_tsv(datasets_qc, file.path(results_dir, "montoro2018_datasets_qc.txt"))


# Create metadata table --------------------------------------------------------
umi_meta <- umi_sobj@meta.data %>% 
  rownames_to_column("BarcodeID") %>% 
  select(BarcodeID, CellGroup, dataset)

fl_meta <- fl_sobj@meta.data %>% 
  rownames_to_column("BarcodeID") %>% 
  select(BarcodeID, CellGroup, dataset)

bind_rows(fl_meta, umi_meta) %>% 
  write_tsv(file.path(results_dir, "montoro2018_cells_metadata.txt"))


# Cell type population sizes ---------------------------------------------------
umi_cells_count <- umi_meta %>% 
  count(CellGroup, dataset)

fl_cells_count <- fl_meta %>% 
  count(CellGroup, dataset)

bind_rows(umi_cells_count, fl_cells_count) %>% 
  write_tsv(file.path(results_dir, "montoro2018_cells_count.txt"))


# Inspect viral entry factors (VEFs) expression levels -------------------------
## Log-normalized UMI-count data 
umi_sobj <- NormalizeData(umi_sobj, normalization.method = "LogNormalize", scale.factor = 10000)

umi_vefs <- extract_gene_expr_q(umi_sobj, groupping_vars = c("CellGroup"), genes = str_to_upper(c("Ace2", "Tmprss2", "Furin", "Ctsl"))) %>% 
  mutate(dataset = "10x Chromium")
fl_vefs <- extract_gene_expr_q(fl_sobj, groupping_vars = c("CellGroup"), genes = str_to_upper(c("Ace2", "Tmprss2", "Furin", "Ctsl"))) %>% 
  mutate(dataset = "Smart-seq2")

bind_rows(umi_vefs, fl_vefs) %>% 
  write_tsv(file.path(results_dir, "montoro2018_vefs_expression_groups.txt"))


# Inspect CFTR expression levels -----------------------------------------------
umi_cftr <- extract_gene_expr_q(umi_sobj, groupping_vars = c("CellGroup"), genes = str_to_upper(c("Cftr"))) %>% 
  mutate(dataset = "10x Chromium")
fl_cftr <- extract_gene_expr_q(fl_sobj, groupping_vars = c("CellGroup"), genes = str_to_upper(c("Cftr"))) %>% 
  mutate(dataset = "Smart-seq2")

bind_rows(umi_cftr, fl_cftr) %>% 
  write_tsv(file.path(results_dir, "montoro2018_cftr_expression_groups.txt"))


# Extract cell populations, normalize and save the data ------------------------
preprocessing <- function(counts_mtx, dataset_name, path_to_result){
  # Perform non-expressed genes and poor-quality cells filtering
  dataset = CreateSeuratObject(counts = counts_mtx, min.cells = 3, min.features = 200)
  
  # Log-normalized data 
  dataset = NormalizeData(dataset, normalization.method = "LogNormalize", scale.factor = 10000)
  
  # Save 
  write.table(as.matrix(GetAssayData(object = dataset, assay = "RNA", slot = "counts")), gzfile(file.path(path_to_result, paste0(dataset_name, ".raw.txt.gz"))), sep = "\t", quote = FALSE) # Un-normalized counts: used by kNN-smoothing, SAVER, scImpute, scVI
  write.table(as.matrix(GetAssayData(object = dataset, assay = "RNA", slot = "data")), gzfile(file.path(path_to_result, paste0(dataset_name, ".norm.txt.gz"))), sep = "\t", quote = FALSE) # Log-normalized counts: used by ALRA, MAGIC 
}


## Total dataset ---------------------------------------------------------------
preprocessing(montoro_umi, "all", results_dir)


## Subset cell populations by cell type ----------------------------------------
### Basal cells ----------------------------------------------------------------
montoro_umi %>% select(ends_with("Basal")) %>% preprocessing("Basal", results_dir)

### Club cells ----------------------------------------------------------------
montoro_umi %>% select(ends_with("Club")) %>% preprocessing("Club", results_dir)

### Ciliated cells ----------------------------------------------------------------
montoro_umi %>% select(ends_with("Ciliated")) %>% preprocessing("Ciliated", results_dir)

