library(tidyverse)

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

methods_labels <- c("MAGIC", "SAVER", "scVI", "scImpute", "kNN-smoothing", "ALRA")


# Load data to plot deltaNRMSE(expression level) vs observed nonzero fraction ----
## Get observed nonzero fraction values in synthetic datasets ------------------
obs_nonzero_fraction <- read_tsv(file.path(results_dir, "nonzero_fraction.txt")) %>% 
  select(GeneID, contains("obs")) %>% 
  pivot_longer(-GeneID, names_to = "dataset", values_to = "obs_nonzero_fraction") %>% 
  mutate(dataset = str_remove(dataset, "_obs_counts.norm"))


# Load NRMSE values ------------------------------------------------------------
nrmse_df <- read_tsv(file.path(results_dir, "nrmse.txt"))

## Create vector of datasets names
datasets2plot <- select(nrmse_df, -GeneID) %>% 
  colnames() %>% 
  str_split_i(pattern = "\\.", i = 1) %>% 
  unique()


plt_deltanrmse_vs_obs_nonzero_fraction <- function(dataset2plot, df_nrmse, df_nonzero_fraction) {
  # Calculate deltaNRMSE = NRMSE(observed) - NRMSE(imputed)
  df_deltanrmse = df_nrmse %>% 
    select(GeneID, contains(dataset2plot)) %>% 
    rename_with(~str_remove(.x, paste0(dataset2plot, "."))) %>% 
    pivot_longer(-c(GeneID, observed), names_to = "method", values_to = "imputed") %>% 
    transmute(GeneID, method, deltanrmse = observed - imputed) %>% 
    mutate(method = factor(method, levels = c("magic_3", "saver", "scvi", "scimpute_0.5", "knn_10", "alra")))
  
  # Filter nonzero fraction values observed in selected dataset
  df_nonzero_fraction = df_nonzero_fraction %>% 
    filter(dataset == dataset2plot)
  
  # Prepare methods labels to change facet labels
  names(methods_labels) = sort(unique(df_deltanrmse$method))
  
  plt = left_join(df_deltanrmse, df_nonzero_fraction, by = "GeneID") %>% 
    drop_na() %>% 
    ggplot(aes(x = obs_nonzero_fraction, y = deltanrmse)) +
    geom_point(size = 1, alpha = 0.1) +
    geom_hline(yintercept = 0, color = "red3", linewidth = 0.5, linetype = "dashed") +
    stat_summary_bin(
      fun = "median", bins = 100,
      color = "#f77f00", linewidth = 0.5, 
      geom = "line") +
    facet_wrap(vars(method), labeller = as_labeller(methods_labels)) + 
    labs(
      x = "Observed nonzero fraction", 
      y = expression(ΔNRMSE == NRMSE[obs] - NRMSE[imp])) 
  
  plt
  
  ggsave(file.path(results_dir, paste0("deltanrmse_", dataset2plot, ".png")), width = 6, height = 4, units = "cm", dpi = 300, scale = 2) 
  ggsave(file.path(results_dir, paste0("deltanrmse_", dataset2plot, ".pdf")), width = 6, height = 4, units = "cm", scale = 2) 
  
  # Limit Y-axis
  plt +
    coord_cartesian(ylim = c(-20, 20))
  
  ggsave(file.path(results_dir, paste0("framed_deltanrmse_", dataset2plot, ".png")), width = 6, height = 4, units = "cm", dpi = 300, scale = 2) 
  ggsave(file.path(results_dir, paste0("framed_deltanrmse_", dataset2plot, ".pdf")), width = 6, height = 4, units = "cm", scale = 2) 
}


## Make plots
walk(datasets2plot, plt_deltanrmse_vs_obs_nonzero_fraction, nrmse_df, obs_nonzero_fraction)

