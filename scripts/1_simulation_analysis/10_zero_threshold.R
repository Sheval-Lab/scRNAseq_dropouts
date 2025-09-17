library(tidyverse)
library(patchwork)
library(ragg)

theme_set(
  theme_classic() +
    theme(
      aspect.ratio = 1,
      panel.spacing.x = unit(1, "lines"),
      axis.title = element_text(size = 10), 
      axis.text = element_text(color = "black", size = 7),
      axis.ticks = element_line(color = "black"),
      strip.text = element_text(size = 8))
)


data_dir <- "results/1_simulation_analysis"
results_dir <- "results/1_simulation_analysis/metrics"


methods_labels <- c("magic" = "MAGIC", "saver" = "SAVER", "scvi" = "scVI", "dca" = "DCA", "scbig" = "scBiG")
cell_type_labels <- c("basal" = "Basal", "ciliated" = "Ciliated", "club" = "Club")


# Plot density of close to zero imputed values ---------------------------------
datasets <- data.frame(
  sym_method = rep(c("symsim", "zinb-wave"), each = 3),
  cell_type = rep(c("basal", "club", "ciliated"), times = 2))

combine_and_plot <- function(sym_method, cell_type) {
  print(paste(sym_method, cell_type))
  
  # # Load observed and imputed (5 methods) data
  obs = read.table(file.path(data_dir, "gene_expression_filtered", paste0(sym_method, "_", cell_type, "_obs_counts.norm.txt.gz")), header = TRUE, row.names = 1)
  magic = read.table(file.path(data_dir, "gene_expression_imputed", paste0(sym_method, "_", cell_type, ".magic_3.txt.gz")), header = TRUE, row.names = 1)
  saver = read.table(file.path(data_dir, "gene_expression_imputed", paste0(sym_method, "_", cell_type, ".saver.txt.gz")), header = TRUE, row.names = 1)
  scvi = read.table(file.path(data_dir, "gene_expression_imputed", paste0(sym_method, "_", cell_type, ".scvi.txt.gz")), header = TRUE, row.names = 1)
  dca = read.table(file.path(data_dir, "gene_expression_imputed", paste0(sym_method, "_", cell_type, ".dca_nb.txt.gz")), header = TRUE, row.names = 1)
  scbig = read.table(file.path(data_dir, "gene_expression_imputed", paste0(sym_method, "_", cell_type, ".scbig.txt.gz")), header = TRUE, row.names = 1)

  # Find minimal non-zero value in observed data
  min_obs = min(obs[obs != 0])
  print(min_obs)

  # Combine imputed data
  imp = reduce(
    list(
      magic %>%
        rownames_to_column("GeneID") %>%
        pivot_longer(-GeneID, names_to = "CellID", values_to = "magic") %>%
        transmute(GeneID_CellID = paste(GeneID, CellID, sep = "_"), magic),
      saver %>%
        rownames_to_column("GeneID") %>%
        pivot_longer(-GeneID, names_to = "CellID", values_to = "saver") %>%
        transmute(GeneID_CellID = paste(GeneID, CellID, sep = "_"), saver),
      scvi %>%
        rownames_to_column("GeneID") %>%
        pivot_longer(-GeneID, names_to = "CellID", values_to = "scvi") %>%
        transmute(GeneID_CellID = paste(GeneID, CellID, sep = "_"), scvi),
      dca %>%
        rownames_to_column("GeneID") %>%
        pivot_longer(-GeneID, names_to = "CellID", values_to = "dca") %>%
        transmute(GeneID_CellID = paste(GeneID, CellID, sep = "_"), dca),
      scbig %>%
        rownames_to_column("GeneID") %>%
        pivot_longer(-GeneID, names_to = "CellID", values_to = "scbig") %>%
        transmute(GeneID_CellID = paste(GeneID, CellID, sep = "_"), scbig)),
    inner_join,
    by = "GeneID_CellID") %>%
    pivot_longer(-GeneID_CellID, names_to = "method") %>%
    mutate(method = factor(method, levels = names(methods_labels)))

  # Save combined imputed data
  write_tsv(imp, gzfile(file.path(results_dir, paste0(sym_method, "_", cell_type, ".imp_counts.txt.gz"))))
  # imp = read_tsv(file.path(results_dir, paste0(sym_method, "_", cell_type, ".imp_counts.txt.gz")))
  
  # Density plot
  dplt = imp %>% 
    mutate(
      value = if_else(value >= 2, 2, value),
      celltype = cell_type) %>% 
    ggplot(aes(x = value)) +
    geom_density(n = 256) +
    facet_grid(celltype ~ method, labeller = as_labeller(c(methods_labels, cell_type_labels))) +
    scale_x_continuous(
      breaks = seq(0, 2, 0.5),
      labels = c(seq(0, 1.5, 0.5), "≥2")) +
    labs(x = "Expression level", y = "Density") +
    theme(aspect.ratio = 1)
  
  ggsave(file.path(results_dir, paste0(sym_method, "_", cell_type, ".imp_counts.dplt.png")), dplt, device = agg_png, width = 14, height = 4, units = "cm", dpi = 300, scaling = 0.7) 
  ggsave(file.path(results_dir, paste0(sym_method, "_", cell_type, ".imp_counts.dplt.pdf")), dplt, width = 14, height = 4, units = "cm", scale = 1/0.7)
  
  dplt_line = dplt +
    geom_vline(xintercept = min_obs, color = "red3", linewidth = 0.5, linetype = "dashed") 
    
  ggsave(file.path(results_dir, paste0(sym_method, "_", cell_type, ".imp_counts.dplt_line.png")), dplt_line, device = agg_png, width = 14, height = 4, units = "cm", dpi = 300, scaling = 0.7) 
  ggsave(file.path(results_dir, paste0(sym_method, "_", cell_type, ".imp_counts.dplt_line.pdf")), dplt_line, width = 14, height = 4, units = "cm", scale = 1/0.7)
}

pwalk(datasets, combine_and_plot)


# Test varying non-zero cutoffs ------------------------------------------------
## SymSim - basal cells as reference for simulation ----------------------------
ss_basal_imp <-  read_tsv(file.path(results_dir, "symsim_basal.imp_counts.txt.gz"))

ss_basal_imp_thresholds <- ss_basal_imp %>% 
  mutate(GeneID = str_extract(GeneID_CellID, "Gene\\d+")) %>% 
  mutate(
    value_t0.01 = if_else(value < 0.01, 0, 1),
    value_t0.05 = if_else(value < 0.05, 0, 1),
    value_t0.1 = if_else(value < 0.1, 0, 1),
    value_t0.2 = if_else(value < 0.2, 0, 1),
    value_t0.5 = if_else(value < 0.5, 0, 1)) %>% 
  select(-c(GeneID_CellID, value)) %>% 
  pivot_longer(starts_with("value_"), names_to = "threshold", values_to = "bin") %>% 
  group_by(GeneID, method, threshold) %>% 
  summarise(nonzero_fraction = mean(bin)) %>% 
  ungroup()
  
write_tsv(ss_basal_imp_thresholds, gzfile(file.path(results_dir, "ss_basal_imp_thresholds.txt.gz")))


## Get nonzero fraction values before and after imputation ---------------------
nonzero_fraction <- read_tsv(file.path(results_dir, "nonzero_fraction.txt"))


## Plot imputed vs true nonzero fraction ---------------------------------------
ss_basal_imp_thresholds %>% 
  left_join(select(nonzero_fraction, GeneID, symsim_basal_true_counts.norm), by = "GeneID") %>% 
  mutate(method = factor(method, levels = names(methods_labels))) %>% 
  ggplot(aes(x = symsim_basal_true_counts.norm, y = nonzero_fraction)) +
  geom_point(size = 0.5, alpha = 0.1) +
  # Diagonal line
  geom_abline(intercept = 0, slope = 1, color = "red3", linewidth = 0.5, linetype = "dashed") +
  facet_grid(
    threshold ~ method,
    labeller = as_labeller(c(
      methods_labels, 
      "value_t0.01" = "0: <0.01",
      "value_t0.05" = "0: <0.05",
      "value_t0.1" = "0: <0.1",
      "value_t0.2" = "0: <0.2",
      "value_t0.5" = "0: <0.5"))) + 
  labs(
    x = "True non-zero fraction", 
    y = "Non-zero fraction after imputation") 


ggsave(file.path(results_dir, "ss_basal_imp_thresholds.png"), device = agg_png, width = 14, height = 12, units = "cm", dpi = 300, scaling = 0.7) 
ggsave(file.path(results_dir, "ss_basal_imp_thresholds.pdf"), width = 14, height = 12, units = "cm", scale = 1/0.7)

