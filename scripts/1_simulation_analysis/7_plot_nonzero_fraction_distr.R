library(tidyverse)
library(ragg)

theme_set(
  theme_classic() +
    theme(
      legend.position = "none",
      axis.title = element_text(size = 10), 
      axis.text = element_text(color = "black", size = 7),
      axis.ticks = element_line(color = "black"))
)


data_dir <- "data/3_montoro2018"
results_dir <- "results/1_simulation_analysis/metrics"


# Load tables with nonzero fractions values ------------------------------------
## Get nonzero fraction values in the reference dataset (Montoro et al., 2018) ----
reference_nonzero_fraction_files <- list.files(data_dir, pattern = "reference_nonzero_fraction", full.names = TRUE)

reference_nonzero_fraction <- map(reference_nonzero_fraction_files, read_tsv, id = "sample") %>% 
  list_rbind() %>% 
  mutate(
    sample = basename(sample),
    type = str_remove(str_split_i(sample, pattern = "\\.", i = 1), "_nonzero_fraction"), 
    celltype = str_split_i(sample, pattern = "\\.", i = 2)) %>% 
  select(type, celltype, nonzero_fraction)


## Get nonzero fraction values in synthetic datasets ---------------------------
sym_nonzero_fraction <- read_tsv(file.path(results_dir, "nonzero_fraction.txt")) %>% 
  select(matches("obs|true")) %>% 
  pivot_longer(everything(), names_to = "sample", values_to = "nonzero_fraction") %>% 
  mutate(
    type = paste(str_split_i(sample, pattern = "_", i = 1), str_split_i(sample, pattern = "_", i = 3), sep = "_"),
    celltype = str_split_i(sample, pattern = "_", i = 2)) %>% 
  select(type, celltype, nonzero_fraction)


## Combine tables with nonzero fraction -----------------------------------------
nonzero_fraction <- bind_rows(reference_nonzero_fraction, sym_nonzero_fraction) %>% 
  drop_na() 

write_tsv(nonzero_fraction, file.path(results_dir, "plt_nonzero_fraction_distr.txt"))


# Plot nonzero fraction distributions in synthetic datasets compared to reference ----
plt_vln_fraction <- function(dataset2plot, df) {
  df %>% 
    # Select data on one of 3 cell types
    filter(celltype == dataset2plot) %>% 
    ggplot(aes(x = type, y = nonzero_fraction, fill = type)) +
    geom_violin() +
    scale_x_discrete(
      labels = c("Reference\ndateset", "SymSim\nObserved", "SymSim\nTrue", "ZINB-WaVE\nObserved", "ZINB-WaVE\nTrue")) +
    scale_fill_manual(values = c("red3", "slategray1", "dodgerblue3", "mistyrose1", "palevioletred")) +
    labs(x = "", y = "Non-zero fraction") 
  
  ggsave(file.path(results_dir, paste0("vln_", dataset2plot, ".png")), device = agg_png, width = 7.5, height = 4, units = "cm", dpi = 300, scaling = 0.8) 
  ggsave(file.path(results_dir, paste0("vln_", dataset2plot, ".pdf")), width = 7.5, height = 4, units = "cm", scale = 1/0.8) 
}


## Make violin plots
walk(unique(nonzero_fraction$celltype), plt_vln_fraction, nonzero_fraction)
