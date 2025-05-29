library(tidyverse)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(ggpubr)
library(patchwork)
library(ragg)

theme_set(
  theme_classic() +
    theme(
      panel.spacing.x = unit(1, "lines"),
      axis.title = element_text(size = 10), 
      axis.text = element_text(color = "black", size = 7),
      axis.ticks = element_line(color = "black"),
      strip.text = element_text(size = 8))
)


data_dir <- "results/2_k562_analysis"
results_dir <- "results/2_k562_analysis/imputed_metrics"


# Create the mapping of downsampled datasets names and sequencing depth values ----
cellranger_stats <- read_tsv(file.path(data_dir, "raw_metrics", "cellranger_stats.txt"))

## Create a named vector
samples_naming <- cellranger_stats %>% pull(`Median UMI Counts per Cell`) 
samples_naming <- floor(samples_naming / 1000)
names(samples_naming) <- cellranger_stats$Sample


# Split genes into 4 groups based on mean expression level (estimated with original dataset) ----
raw_nonzero_fraction <- read_tsv(file.path(data_dir, "raw_metrics", "raw_nonzero_fraction.txt"))
raw_avg_norm_counts <- read_tsv(file.path(data_dir, "raw_metrics", "raw_avg_norm_counts.txt"))

## Split dataset into 4 quartiles 
genes_in_4_groups_df <- raw_avg_norm_counts %>% 
  select(GeneID, sample_100.norm_counts) %>% 
  mutate(gene_quartile = as.factor(paste0("Q", ntile(sample_100.norm_counts, 4))))

genes_in_4_groups <- genes_in_4_groups_df$gene_quartile
names(genes_in_4_groups) <- genes_in_4_groups_df$GeneID


## Plot expression level in 4 groups of genes ----------------------------------
ggplot(genes_in_4_groups_df, aes(x = gene_quartile, y = sample_100.norm_counts, fill = gene_quartile)) +
  geom_violin(alpha = 0.6) +
  geom_boxplot(width = 0.15) +
  scale_y_log10(labels = scales::number) +
  scale_fill_brewer(palette = "Dark2") +
  labs(x = "", y = "Mean expression level\n(mean logcounts)") +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text = element_text(color = "black"))

ggsave(file.path(results_dir, "gene_expression_quartiles.png"), width = 6, height = 4.5, units = "cm", dpi = 300, scale = 2)  
ggsave(file.path(results_dir, "gene_expression_quartiles.pdf"), width = 6, height = 4.5, units = "cm", scale = 2)  


## Plot heatmap of mean gene expression level in downsampled datasets ---------- 
raw_avg_norm_counts_mtx <- raw_avg_norm_counts %>% 
  select(where(is.numeric)) %>% 
  as.matrix()
rownames(raw_avg_norm_counts_mtx) <- raw_avg_norm_counts$GeneID

hm_expression <- Heatmap(
  raw_avg_norm_counts_mtx, 
  show_column_dend = FALSE,
  show_row_names = FALSE, 
  cluster_rows = FALSE, 
  cluster_columns = FALSE, 
  row_split = genes_in_4_groups[rownames(raw_avg_norm_counts_mtx)], 
  column_labels = samples_naming,
  column_names_rot = 0, 
  row_names_rot = 90, 
  row_title_rot = 0, 
  name = "Expression\nlevel", 
  row_title_gp = gpar(fontsize = c(14)),
  column_names_gp = gpar(fontsize = c(10)), 
  column_names_centered = TRUE, 
  # column_title = "Sequencing depth, mean reads per cell (in thousands)", 
  column_title = "Sequencing depth, median UMI counts per cell (in thousands)", 
  column_title_side = "bottom", 
  column_title_gp = gpar(fontsize = 14),
  col = colorRamp2(quantile(raw_avg_norm_counts$sample_100.norm_counts, probs = seq(0, 1, 0.1)), rev(brewer.pal(11, "Spectral"))),
  heatmap_legend_param = list(
    legend_width = unit(10, "cm"),
    break_dist = rep(1, 5),
    at = round(quantile(raw_avg_norm_counts$sample_100.norm_counts, probs = seq(0, 1, 0.2)), 3)))

agg_png(file.path(results_dir, "hm_raw_avg_norm_counts.umi.png"), width = 7.5, height = 4.5, units = "cm", res = 300, scaling = 0.5)
hm_expression
dev.off()

pdf(file.path(results_dir, "hm_raw_avg_norm_counts.umi.pdf"), width = 7.5, height = 4.5)
hm_expression
dev.off()


# Plot boxplots of nonzero fraction of raw and imputed data --------------------
plt_box_fraction <- function(df, method2plot = "", nrows = 1) {
  # Calculate median in original dataset to plot it as a red line
  ref_df = df %>% 
    filter(sample == "sample_100") %>% 
    group_by(gene_quartile) %>% 
    summarise(fraction = median(fraction))
  
  plt = df %>%
    filter(method %in% c("norm_counts", method2plot)) %>% 
    ggplot(aes(x = sample, y = fraction, fill = gene_quartile)) +
    geom_boxplot(aes(alpha = sample), linewidth = 0.25, outlier.size = 0.25) +
    # Add median of original dataset as a red line
    geom_hline(
      aes(yintercept = fraction),
      ref_df,
      color = "red3",
      linetype = "dashed") +
    facet_wrap(vars(gene_quartile), nrow = nrows) +
    scale_x_discrete(labels = samples_naming) +
    scale_y_continuous(breaks = seq(0, 1, by = 0.25)) +
    scale_fill_brewer(palette = "Dark2") +
    labs(
      # x = "Sequencing depth, mean reads per cell (in thousands)", 
      x = "Sequencing depth, median UMI counts per cell (in thousands)", 
      y = "Non-zero fraction") +
    theme(
      # aspect.ratio = 0.4,
      legend.position = "none",
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
  
  # Save plots
  if (method2plot == "") {
    file_name_suffix = "norm_counts"
  } else {
    file_name_suffix = method2plot
  }
  
  plt 
}


## Boxplots of raw nonzero fraction --------------------------------------------
raw_nonzero_fraction_lng <- raw_nonzero_fraction %>% 
  # Add gene group ID
  left_join(select(genes_in_4_groups_df, GeneID, gene_quartile), by = "GeneID") %>% 
  pivot_longer(starts_with("sample_"), values_to = "fraction") %>% 
  separate(name, into = c("sample", "method"), sep = "\\.", extra = "merge")


### Make boxplots
raw_nonzero_fraction_lng %>% plt_box_fraction(nrows = 2)

ggsave(file.path(results_dir, "boxplot.umi.png"), device = agg_png, width = 7.5, height = 4.5, units = "cm", dpi = 300, scaling = 0.7) 
ggsave(file.path(results_dir, "boxplot.umi.pdf"), width = 7.5, height = 4.5, units = "cm", scale = 1/0.7) 


## Boxplots of imputed nonzero fraction ----------------------------------------
imputed_nonzero_fraction <- read_tsv(file.path(results_dir, "imputed_nonzero_fraction.txt.gz"))

imputed_nonzero_fraction_lng <- imputed_nonzero_fraction %>% 
  # Add gene group ID
  left_join(select(genes_in_4_groups_df, GeneID, gene_quartile), by = "GeneID") %>% 
  # Add nonzero fraction values observed in original dataset
  left_join(select(raw_nonzero_fraction, GeneID, sample_100.norm_counts), by = "GeneID") %>% 
  pivot_longer(starts_with("sample_"), values_to = "fraction") %>% 
  separate(name, into = c("sample", "method"), sep = "\\.", extra = "merge")

methods2plot <- c("magic_3", "saver", "scvi", "scimpute_0.5", "knn_10", "alra")

plt_box_fraction <- function(df, method2plot = "", nrows = 1) {
  # Calculate median in original dataset to plot it as a red line
  ref_df = df %>% 
    filter(sample == "sample_100") %>% 
    group_by(gene_quartile) %>% 
    summarise(fraction = median(fraction))
  
  df = df %>%
    filter(method %in% c("norm_counts", method2plot)) %>% 
    mutate(method = if_else(method == "norm_counts", paste(method2plot, collapse = ","), method)) %>% 
    separate_rows(method, sep = ",")
  
  plt = df %>%
    filter(method %in% c(method2plot)) %>% 
    ggplot(aes(x = sample, y = fraction, fill = gene_quartile)) +
    geom_boxplot(aes(alpha = sample), linewidth = 0.25, outlier.size = 0.25) +
    # Add median of original dataset as a red line
    geom_hline(
      aes(yintercept = fraction),
      ref_df,
      color = "red3",
      linetype = "dashed") +
    # facet_wrap(vars(gene_quartile), nrow = nrows) +
    facet_grid(
      method ~ gene_quartile, 
      labeller = as_labeller(c(
        "Q1" = "Q1", 
        "Q2" = "Q2", 
        "Q3" = "Q3", 
        "Q4" = "Q4", 
        "magic_3" = "MAGIC", 
        "saver" =  "SAVER", 
        "scvi" = "scVI", 
        "scimpute_0.5" = "scImpute", 
        "knn_10" = "kNN-\nsmoothing", 
        "alra" = "ALRA"))) +
    scale_x_discrete(labels = samples_naming) +
    scale_y_continuous(breaks = seq(0, 1, by = 0.25)) +
    scale_fill_brewer(palette = "Dark2") +
    labs(
      # x = "Sequencing depth, mean reads per cell (in thousands)", 
      x = "Sequencing depth, median UMI counts per cell (in thousands)", 
      y = "Non-zero fraction") +
    theme(
      # aspect.ratio = 0.4,
      legend.position = "none",
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
  
  # Save plots
  if (method2plot == "") {
    file_name_suffix = "norm_counts"
  } else {
    file_name_suffix = method2plot
  }
  
  plt 
  
}

boxpls_nzf <- map(methods2plot, plt_box_fraction, df = imputed_nonzero_fraction_lng)

boxpls_nzf[[1]] <- boxpls_nzf[[1]] + 
  theme(axis.text.x = element_blank(), title = element_text(size = rel(0.7)))

boxpls_nzf[[2]] <- boxpls_nzf[[2]] + 
  theme(axis.text.x = element_blank(), title = element_text(size = rel(0.7)))

boxpls_nzf[[3]] <- boxpls_nzf[[3]] + 
  theme(axis.text.x = element_blank(), title = element_text(size = rel(0.7)))

boxpls_nzf[[4]] <- boxpls_nzf[[4]] + 
  theme(axis.text.x = element_blank(), title = element_text(size = rel(0.7)))

boxpls_nzf[[5]] <- boxpls_nzf[[5]] + 
  theme(axis.text.x = element_blank(), title = element_text(size = rel(0.7)))

boxpls_nzf[[6]] <- boxpls_nzf[[6]] + 
  theme(title = element_text(size = rel(0.7)))

wrap_plots(boxpls_nzf, nrow = 6) + plot_layout(axis_titles = "collect")

ggsave(file.path(results_dir, "boxplot_all2.umi.png"), device = agg_png, width = 14, height = 13, units = "cm", dpi = 300, scaling = 0.7) 
ggsave(file.path(results_dir, "boxplot_all2.umi.pdf"), width = 14, height = 13, units = "cm", scale = 1/0.7)


# Plot medians of nonzero fractions in imputed data when imputation parameters were varied ----
params_labels <- data.frame(
  method = c("magic", "scimpute", "knn"),
  param_label = paste("Parameter", c("t", "drop_thre", "k")))

plt_dot_params <- function(df, method2plot = "") {
  # Calculate median in original dataset to plot it as a red line
  ref_df = df %>% 
    filter(sample == "sample_100") %>% 
    group_by(gene_quartile) %>% 
    summarise(fraction = median(fraction))
  
  # Select parameter name
  param_label = params_labels %>% 
    filter(method == method2plot) %>% 
    pull(param_label)
  
  df = df %>%
    filter(str_detect(method, paste("norm_counts", method2plot, sep = "|"))) %>% 
    mutate(
      method_param = str_extract(method, "[\\d\\.]+$"),
      method_param = replace_na(method_param, "reference\nmedian"),
      method_param = fct_inseq(method_param))
  
  n_params = df %>% 
    filter(method_param != "reference\nmedian") %>% 
    pull(method_param) %>% 
    n_distinct()
  
  plt = df %>% 
    ggplot(aes(x = sample, y = fraction)) +
    stat_summary(aes(color = method_param), geom = "point", fun = "median", alpha = 0.7) +
    # Add median of original dataset as a red line
    geom_hline(
      aes(yintercept = fraction),
      ref_df,
      color = "red3",
      linetype = "dashed") +
    facet_wrap(vars(gene_quartile), nrow = 1) +
    scale_x_discrete(labels = samples_naming) +
    scale_y_continuous(breaks = seq(0, 1, by = 0.25)) +
    scale_color_manual(values = c(brewer.pal(n_params, "Set2"), "red3")) + 
    guides(color = guide_legend(override.aes = list(size = 2))) +
    labs(
      # x = "Sequencing depth, mean reads per cell (in thousands)", 
      x = "Sequencing depth, median UMI counts per cell (in thousands)",  
      y = "Non-zero fraction",
      color = param_label) + 
    theme(
      legend.text = element_text(size = rel(0.6)),
      legend.key.size = unit(0.5, "lines"),
      axis.title = element_text(size = rel(1.2)),
      axis.text.x = element_text(angle = 90, hjust = 1, size = rel(0.8)))
  
  # Save plots
  file_name_suffix = method2plot

  plt 
}

dotplts <- map(params_labels$method, plt_dot_params, df = imputed_nonzero_fraction_lng)

dotplts[[1]] <- dotplts[[1]] +
  labs(title = "MAGIC") +
  theme(axis.text.x = element_blank(), title = element_text(size = rel(0.7)))

dotplts[[2]] <- dotplts[[2]] +
  labs(title = "scImpute") +
  theme(axis.text.x = element_blank(), title = element_text(size = rel(0.7)))

dotplts[[3]] <- dotplts[[3]] +
  labs(title = "kNN-smoothing") +
  theme(title = element_text(size = rel(0.7)))

wrap_plots(dotplts, nrow = 3) + plot_layout(axis_titles = "collect", guides = "keep")

ggsave(file.path(results_dir, "dotplot_all.umi.png"), device = agg_png, width = 14, height = 9, units = "cm", dpi = 300, scaling = 0.7) 
ggsave(file.path(results_dir, "dotplot_all.umi.pdf"), width = 14, height = 9, units = "cm", scale = 1/0.7)

