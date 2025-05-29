library(tidyverse)
library(ggtext)
library(patchwork)
library(ggh4x)
library(ragg)

data_dir <- "results/4_cftr_analysis"
results_dir <- "results/4_cftr_analysis/metrics"

theme_set(
  theme_classic() +
    theme(
      axis.text = element_text(color = "black")))

cols_dataset <- c("scRNA ISH" = "red3", "Drop-seq" = "#1B9E77", "10x Chromium" = "#a8dadc")

celltypes2plot <- c("Basal", "Ciliated", "Ionocyte", "Secretory")


# Cell type population sizes ---------------------------------------------------
okuda_cells_meta <- read_tsv(file.path(data_dir, "gene_expression_filtered", "okuda2021_cells_metadata.txt"))

okuda_cells_count <- okuda_cells_meta %>% 
  count(Airway, CellGroup, dataset) %>% 
  filter(CellGroup %in% c("Basal", "Ciliated", "Ionocyte", "Secretory"))

write_tsv(cells_count, file.path(results_dir, "okuda2021_cells_count.txt"))


# scRNA ISH measurments of the proportion of CFTR+ cells -----------------------
## From Fig. 4C (Okuda et al., 2021)
okuda_scrnaish_df <- data.frame(
  Airway = rep(c("LAE", "SAE"), each = 4),
  CellGroup = rep(c("Basal", "Ciliated", "Ionocyte", "Secretory"), 2),
  # GeneID = rep("CFTR", 8),
  # avg_count = rep(NA, 8),
  nonzero_fraction = c(0.68, 0.265, 0.955, 0.852, 0.489, 0.376, 0.743, 0.894),
  # gene_quartile = rep(NA, 8),
  dataset = rep("scRNA ISH", 8),
  run = rep("raw", 8))


# Load scRNA-seq data ----------------------------------------------------------
## Nonzero fraction of CFTR before imputation
okuda_raw_metrics <- read_tsv(file.path(data_dir, "gene_expression_filtered", "okuda2021_cftr_expression_groups.txt"))

## Nonzero fraction of CFTR after imputation
okuda_imputed_metrics <- read_tsv(file.path(data_dir, "metrics", "okuda2021_metrics_in_unsplit_ds.txt"))

## Combine metrics
okuda_metrics_combined <- bind_rows(
  okuda_scrnaish_df,
  transmute(okuda_raw_metrics, Airway, CellGroup, nonzero_fraction, dataset, run = "raw"),
  select(okuda_imputed_metrics, Airway, CellGroup, nonzero_fraction, dataset, run)) %>% 
  mutate(type = if_else(run == "raw", "Raw", "Imputed"))

write_tsv(okuda_metrics_combined, file.path(results_dir, "okuda2021_combined_metrics.txt"))

okuda_metrics2plot <- okuda_metrics_combined %>% 
  filter(str_detect(run, "raw|all_both")) %>% 
  filter(CellGroup %in% celltypes2plot) %>% 
  mutate(
    CellGroup = fct_rev(as.factor(CellGroup)),
    Airway_full = if_else(Airway == "LAE", "Large airway epithelial", "Small airway epithelial"),
    dataset = factor(dataset, levels = c("scRNA ISH", "Drop-seq", "10x Chromium")),
    type = factor(type, levels = c("Raw", "Imputed"))) 

# Barplot ----------------------------------------------------------------------
okuda_metrics2plot %>% 
  ggplot(aes(x = fct_rev(dataset), y = nonzero_fraction, fill = dataset, color = type)) +
  geom_col(position = position_dodge2(preserve = "single")) +
  facet_grid2(Airway_full ~ fct_rev(CellGroup), axes = "x", remove_labels = "x", scales = "free_x") +
  scale_color_manual(values = c("gray80", "black"), name = "") +
  scale_fill_manual(values = cols_dataset, breaks = names(cols_dataset), name = "Dataset") +
  guides(color = guide_legend(override.aes = list(fill=NA)), order = 1) +
  labs(x = "", y = "<i>CFTR</i> non-zero fraction") +
  theme(
    axis.title.y = element_markdown(),
    axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(results_dir, paste0("barplot_cftr_full.png")), device = agg_png, width = 14, height = 8, units = "cm", dpi = 300, scaling = 0.7)
ggsave(file.path(results_dir, paste0("barplot_cftr_full.pdf")), width = 14, height = 7, units = "cm", scale = 1/0.7)

