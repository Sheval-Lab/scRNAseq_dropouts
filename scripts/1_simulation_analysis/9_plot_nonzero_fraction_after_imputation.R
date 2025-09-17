library(tidyverse)
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

results_dir <- "results/1_simulation_analysis/metrics"

methods_labels <- c("MAGIC", "SAVER", "scVI", "DCA", "scBiG", "scImpute", "kNN-smoothing", "ALRA")


# Load nonzero fraction data ---------------------------------------------------
## Get nonzero fraction values before and after imputation ---------------------
nonzero_fraction <- read_tsv(file.path(results_dir, "nonzero_fraction.txt"))

## Reshape dataframe for plotting ----------------------------------------------
nonzero_fraction_lng <- nonzero_fraction %>% 
  pivot_longer(-GeneID, names_to = "dataset_method", values_to = "nonzero_fraction") %>% 
  mutate(
    dataset = str_extract(dataset_method, "[a-z-]+_[a-z]+"),
    method = case_when(
      str_detect(dataset_method, "obs") ~ "observed",
      str_detect(dataset_method, "true") ~ "true",
      TRUE ~ str_remove(str_extract(dataset_method, "\\.[a-z0-9_.]+$"), "^\\."))) %>% 
  select(-dataset_method) %>% 
  pivot_wider(names_from = method, values_from = nonzero_fraction) %>% 
  pivot_longer(-c(GeneID, dataset, observed, true), names_to = "method", values_to = "imputed") %>% 
  drop_na() %>% 
  mutate(method = factor(method, levels = c("magic_3", "saver", "scvi", "dca_nb", "scbig", "scimpute_0.5", "knn_10", "alra")))

## Create vector of datasets names
datasets2plot <- unique(nonzero_fraction_lng$dataset)


# Plot imputed vs true nonzero fraction ----------------------------------------
plt_nonzero_fraction_corr <- function(dataset2plot, df_nonzero_fraction) {
  # Filter nonzero fraction values for selected dataset
  df_nonzero_fraction = df_nonzero_fraction %>% 
    filter(dataset == dataset2plot)
  
  # Prepare methods labels to change facet labels
  names(methods_labels) = sort(unique(df_nonzero_fraction$method))
  
  # Plot imputed vs true values
  plt_true = df_nonzero_fraction %>% 
    ggplot(aes(x = true, y = imputed)) +
    geom_point(size = 0.5, alpha = 0.1) +
    # Diagonal line
    geom_abline(intercept = 0, slope = 1, color = "red3", linewidth = 0.5, linetype = "dashed") +
    facet_wrap(vars(method), labeller = as_labeller(methods_labels)) +
    labs(
      x = "True non-zero fraction", 
      y = "Non-zero fraction after imputation") 
  
  plt_true
  
  ggsave(file.path(results_dir, paste0("imp_vs_true_nzf_", dataset2plot, ".png")), device = agg_png, width = 10, height = 6.5, units = "cm", dpi = 300, scaling = 0.8) 
  ggsave(file.path(results_dir, paste0("imp_vs_true_nzf_", dataset2plot, ".pdf")), width = 10, height = 6.5, units = "cm", scale = 1/0.8) 
  
  ## Add line - median per bin
  plt_true +
    stat_summary_bin(
      fun = "median", bins = 100,
      color = "#f77f00", linewidth = 0.5, 
      geom = "line")
  
  ggsave(file.path(results_dir, paste0("imp_vs_true_nzf_wline_", dataset2plot, ".png")), device = agg_png, width = 10, height = 6.5, units = "cm", dpi = 300, scaling = 0.8) 
  ggsave(file.path(results_dir, paste0("imp_vs_true_nzf_wline_", dataset2plot, ".pdf")), width = 10, height = 6.5, units = "cm", scale = 1/0.8) 
  
  
  # Plot imputed vs observed values 
  plt_obs = df_nonzero_fraction %>% 
    ggplot(aes(x = observed, y = imputed)) +
    geom_point(size = 0.5, alpha = 0.1) +
    # Diagonal line
    geom_abline(intercept = 0, slope = 1, color = "red3", linewidth = 0.5, linetype = "dashed") +
    stat_summary_bin(
      fun = "median", bins = 100,
      color = "#f77f00", linewidth = 0.5, 
      geom = "line")+
    facet_wrap(vars(method), labeller = as_labeller(methods_labels)) +
    labs(
      x = "Observed non-zero fraction", 
      y = "Non-zero fraction after imputation") 
  
  plt_obs
  
  ggsave(file.path(results_dir, paste0("imp_vs_obs_nzf_", dataset2plot, ".png")), device = agg_png, width = 10, height = 6.5, units = "cm", dpi = 300, scaling = 0.8) 
  ggsave(file.path(results_dir, paste0("imp_vs_obs_nzf_", dataset2plot, ".pdf")), width = 10, height = 6.5, units = "cm", scale = 1/0.8) 
  
  ## Add line - median per bin
  plt_obs +
    stat_summary_bin(
      fun = "median", bins = 100,
      color = "#f77f00", linewidth = 0.5, 
      geom = "line")
  
  ggsave(file.path(results_dir, paste0("imp_vs_obs_nzf_wline_", dataset2plot, ".png")), device = agg_png, width = 10, height = 6.5, units = "cm", dpi = 300, scaling = 0.8) 
  ggsave(file.path(results_dir, paste0("imp_vs_obs_nzf_wline_", dataset2plot, ".pdf")), width = 10, height = 6.5, units = "cm", scale = 1/0.8) 
  
}

plt_nonzero_fraction_corr <- function(dataset2plot, df_nonzero_fraction) {
  # Filter nonzero fraction values for selected dataset
  df_nonzero_fraction = df_nonzero_fraction %>% 
    filter(dataset == dataset2plot)
  
  # Prepare methods labels to change facet labels
  names(methods_labels) = sort(unique(df_nonzero_fraction$method))
  
  point_color = case_when(
    str_detect(dataset2plot, "symsim") ~ "dodgerblue3",
    str_detect(dataset2plot, "zinb-wave") ~ "palevioletred")
  
  # Plot imputed vs true values
  plt_true = df_nonzero_fraction %>% 
    ggplot(aes(x = true, y = imputed)) +
    geom_point(shape = 16, size = 0.5, color = point_color, alpha = 0.05) +
    # Diagonal line
    geom_abline(intercept = 0, slope = 1, color = "black", linewidth = 0.5, linetype = "dashed") +
    facet_wrap(vars(method), labeller = as_labeller(methods_labels), ncol = 4) +
    labs(
      x = "True non-zero fraction", 
      y = "Non-zero fraction after imputation") 
  
  plt_true
  
  ggsave(file.path(results_dir, paste0("col_imp_vs_true_nzf_", dataset2plot, ".png")), device = agg_png, width = 11.5, height = 6.5, units = "cm", dpi = 300, scaling = 0.8) 
  ggsave(file.path(results_dir, paste0("col_imp_vs_true_nzf_", dataset2plot, ".pdf")), width = 11.5, height = 6.5, units = "cm", scale = 1/0.8) 
}

## Make plots
walk(datasets2plot, plt_nonzero_fraction_corr, nonzero_fraction_lng)
