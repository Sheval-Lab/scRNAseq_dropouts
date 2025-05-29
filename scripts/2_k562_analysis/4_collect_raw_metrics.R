library(tidyverse)
library(Seurat)
library(RSpectra) # RunUMAP doesn't work without it...
library(ragg)
source("scripts/functions.R")

cellragre_outs <- "data/2_k562_dataset/cellranger_outs"
data_dir <- "results/2_k562_analysis/gene_expression_filtered"
results_dir <- "results/2_k562_analysis/raw_metrics"


# Sequencing depth and number of detected genes dependency ---------------------
## Collect data ----------------------------------------------------------------
collect_cellranger_stats <- function(path_to_data) {
  # Get sample ID
  sample_name = basename(path_to_data)
  
  # Read in metrics summary file
  sample_metrics = read_csv(file.path(path_to_data, "outs", "metrics_summary.csv"))
  
  # Select required columns
  sample_metrics = select(sample_metrics, "Estimated Number of Cells", "Mean Reads per Cell", "Median Genes per Cell", "Number of Reads", "Total Genes Detected", "Median UMI Counts per Cell")
  
  # Add sample ID to metrics table
  sample_metrics = mutate(sample_metrics, Sample = sample_name, .before = 1)
  
  return(sample_metrics)
}


cellranger_stats <- map(list.dirs(cellragre_outs, recursive = FALSE), collect_cellranger_stats) %>% list_rbind()
write_tsv(cellranger_stats, file.path(results_dir, "cellranger_stats.txt"))


## Plot `Median Genes per Cell` ~ `Median UMI Counts per Cell` -----------------
cellranger_stats <- read_tsv(file.path(results_dir, "cellranger_stats.txt"))

ggplot(cellranger_stats, aes(x = `Median UMI Counts per Cell`, y = `Median Genes per Cell`)) +
  geom_point() +
  scale_x_continuous(
    breaks = seq(0, 200, 40) * 1000,
    labels = scales::number_format(big.mark = ","),
    limits = c(0, 180000)) +
  scale_y_continuous(
    breaks = seq(4, 10, 2) * 1000,
    labels = scales::number_format(big.mark = ",")) +
  labs(x = "Median UMI Counts per Cell", y = "Median\ndetected genes per cell") +
  theme_classic() +
  theme(
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(angle = 30, hjust = 1))

ggsave(file.path(results_dir, "cellranger_stats.umi.png"), device = agg_png, width = 5.5, height = 4.5, units = "cm", dpi = 300, scaling = 0.7) 
ggsave(file.path(results_dir, "cellranger_stats.umi.pdf"), width = 5.5, height = 4.5, units = "cm", scale = 1/0.7) 


ggplot(cellranger_stats, aes(x = `Mean Reads per Cell`, y = `Median UMI Counts per Cell`)) +
  geom_point() +
  scale_x_continuous(
    breaks = seq(0, 600, 100) * 1000,
    labels = scales::number_format(big.mark = ",")) +
  scale_y_continuous(
    breaks = seq(0, 200, 40) * 1000,
    labels = scales::number_format(big.mark = ","),
    limits = c(0, 180000)) +
  theme_classic() +
  theme(
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(angle = 30, hjust = 1))

ggsave(file.path(results_dir, "cellranger_stats.umi_reads.png"), device = agg_png, width = 5.5, height = 4.5, units = "cm", dpi = 300, scaling = 0.7) 
ggsave(file.path(results_dir, "cellranger_stats.umi_reads.pdf"), width = 5.5, height = 4.5, units = "cm", scale = 1/0.7) 


# Collect data on the fraction of cells with non-zero expression of a gene and mean normalized expression of a gene ----
## Get list of files with normalized counts matrices
norm_counts_files <- list.files(data_dir, pattern = "norm_counts", full.names = TRUE)

## Calculate and aggregate metrics tables
raw_metrics <- collect_metrics(norm_counts_files)
raw_dropout_rate <- raw_metrics$dropout_rate
raw_nonzero_fraction <- raw_metrics$nonzero_fraction
raw_avg_norm_counts <- raw_metrics$avg_norm_counts


## Save metrics into files -----------------------------------------------------
write_tsv(raw_dropout_rate, file.path(results_dir, "raw_dropout_rate.txt"))
write_tsv(raw_nonzero_fraction, file.path(results_dir, "raw_nonzero_fraction.txt"))
write_tsv(raw_avg_norm_counts, file.path(results_dir, "raw_avg_norm_counts.txt"))


# Show cells homogeneity of original dataset -----------------------------------
## Ggplot theme for UMAP plot --------------------------------------------------
theme_umap <- function() {
  list(
    labs(x = "UMAP 1", y = "UMAP 2"),
      theme(
        axis.text = element_blank(),
        axis.ticks = element_blank()))
}

## Preprocess dataset with Seurat ----------------------------------------------
orig_sobj <- readRDS(file.path(data_dir, "sample_100.sobj.rds"))
orig_sobj[["percent.mt"]] <- PercentageFeatureSet(orig_sobj, pattern = "^MT-")
orig_sobj <- FindVariableFeatures(orig_sobj) %>% 
  ScaleData() %>% 
  RunPCA()


## Plot PCA ------------------------------------------------------------------
var_explained <- round(orig_sobj@reductions$pca@stdev ^ 2 / orig_sobj@reductions$pca@misc$total.variance * 100, 2)[1:2]

DimPlot(orig_sobj, reduction = "pca", cols = alpha("dodgerblue3", 0.5), pt.size = 0.1) + NoLegend() +
  labs(
    x = paste0("PC1: ", var_explained[1], "% variance"),
    y = paste0("PC2: ", var_explained[2], "% variance"))

ggsave(file.path(results_dir, "pca.png"), device = agg_png, width = 5.5, height = 4.5, units = "cm", dpi = 300, scaling = 0.7)
ggsave(file.path(results_dir, "pca.pdf"), width = 5.5, height = 4.5, units = "cm", scale = 1/0.7)

ElbowPlot(orig_sobj)


## Plot UMAPs ------------------------------------------------------------------
orig_sobj <- RunUMAP(orig_sobj, dims = 1:10)

### UMAP with no color identity ------------------------------------------------
DimPlot(orig_sobj, reduction = "umap") + NoLegend() + theme_umap()

ggsave(file.path(results_dir, "umap.png"), width = 6, height = 4.5, units = "cm", dpi = 300, scale = 2)
ggsave(file.path(results_dir, "umap.pdf"), width = 6, height = 4.5, units = "cm", scale = 2)

### UMAP showing mitochondrial genes percent -----------------------------------
FeaturePlot(orig_sobj, features = "percent.mt", reduction = "umap") + theme_umap()

ggsave(file.path(results_dir, "umap_mt.png"), width = 6, height = 4.5, units = "cm", dpi = 300, scale = 2)
ggsave(file.path(results_dir, "umap_mt.pdf"), width = 6, height = 4.5, units = "cm", scale = 2)

### UMAP showing cell cycle phase (predicted from expression of cell cycle gene markers) ----
orig_sobj <- CellCycleScoring(orig_sobj, s.features = cc.genes$s.genes, g2m.features = cc.genes$g2m.genes)
DimPlot(orig_sobj, group.by = "Phase", reduction = "umap") + theme_umap()

ggsave(file.path(results_dir, "umap_cc.png"), width = 6, height = 4.5, units = "cm", dpi = 300, scale = 2)
ggsave(file.path(results_dir, "umap_cc.pdf"), width = 6, height = 4.5, units = "cm", scale = 2)

