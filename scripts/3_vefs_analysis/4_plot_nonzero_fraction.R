library(tidyverse)
library(ggtext)
library(patchwork)
library(ggh4x)
library(ragg)

data_dir <- "results/3_vefs_analysis"
results_dir <- "results/3_vefs_analysis/metrics"

theme_set(
  theme_classic() +
    theme(
      axis.text = element_text(color = "black")))


cols_dataset <- c("Smart-seq2" = "red3", "10x Chromium" = "#a8dadc")


# Load scRNA-seq data ----------------------------------------------------------
## Get cell populations sizes in 10x Chromium dataset
cells_numbers <- read_tsv(file.path(data_dir, "gene_expression_filtered", "montoro2018_cells_count.txt")) 
  

## Nonzero fraction of VEFs genes in reference dataset (Smart-seq2) and 
## in experimental (10x Chromium) dataset before and after imputation
metrics_combined <- read_tsv(file.path(results_dir, "montoro2018_metrics_in_unsplit_ds.txt"))

metrics_combined <- metrics_combined %>% 
  transmute(
    GeneID,
    CellGroup, 
    nonzero_fraction, 
    dataset = if_else(str_detect(run, "TPM"), "Smart-seq2", "10x Chromium"),
    run, 
    type = if_else(str_detect(run, "norm|TPM"), "Raw", "Imputed")) %>% 
  

write_tsv(metrics_combined, file.path(results_dir, "combined_metrics.txt"))


metrics2plot <- metrics_combined %>% 
  filter(CellGroup != "Goblet") %>% # goblet cells are absent in Smart-seq2 dataset
  mutate(
    dataset = factor(dataset, levels = c("Smart-seq2", "10x Chromium")),
    type = factor(type, levels = c("Raw", "Imputed")))


# Imputation on datasets split by cell type population -------------------------
ref_nonzero_fraction <- read_tsv(file.path(data_dir, "gene_expression_filtered", "montoro2018_vefs_expression_groups.txt"))
ref_nonzero_fraction <- ref_nonzero_fraction %>% 
  select(-c(avg_count, gene_quartile)) %>% 
  mutate(type = "Raw") %>% 
  filter(dataset == "Smart-seq2")

nonzero_fraction <- read_tsv(file.path(results_dir, "nonzero_fraction.txt"))


split_ds_metrcis2plot <- nonzero_fraction %>%
  filter(GeneID %in% c("Ace2", "Tmprss2", "Furin", "Ctsl")) %>% 
  pivot_longer(where(is.numeric), names_to = "run", values_to = "nonzero_fraction") %>% 
  filter(str_detect(run, "Total|TPM", negate = TRUE)) %>% 
  transmute(
    CellGroup = str_extract(run, "^[A-Za-z]+"),
    GeneID, 
    nonzero_fraction,
    dataset = "10x Chromium",
    type = if_else(str_detect(run, "norm"), "Raw", "Imputed")) %>% 
  bind_rows(ref_nonzero_fraction) %>% 
  mutate(type = factor(type, levels = c("Raw", "Imputed"))) 


# Combine two imputed set of values into one barplot ---------------------------
both_metrics2plot <- bind_rows(
  metrics2plot %>% 
    select(-run) %>% 
    mutate(type = if_else(type == "Imputed", "Imputed total", type)),
  split_ds_metrcis2plot %>% 
    filter(type == "Imputed") %>% 
    mutate(type = "Imputed split") %>% 
    select(-n_cells)) %>% 
  left_join(cells_numbers, by = c("CellGroup", "dataset")) %>% 
  mutate(type = factor(type, levels = c("Raw", "Imputed total", "Imputed split")))

both_metrics2plot_CellGroup_labels <- both_metrics2plot %>% 
  select(CellGroup, dataset, n) %>%
  distinct() %>%
  arrange(dataset) %>% 
  group_by(CellGroup) %>% 
  mutate(
    CellGroup_label = paste(n, collapse = " and "),
    CellGroup_label = paste0(CellGroup, "<br><span style='font-size:8pt'>(", CellGroup_label, " cells)</span>")) %>% 
  ungroup() %>% 
  distinct(CellGroup, CellGroup_label)

both_metrics2plot_dataset_labels <- both_metrics2plot %>% 
  select(CellGroup, dataset, n) %>%
  distinct() %>%
  mutate(dataset_label = paste0(dataset, "<br><span style='font-size:7pt'>(", n, " cells)</span>")) %>% 
  arrange(CellGroup, dataset) %>% 
  select(-n)


## With cell number info in X-labels -------------------------------------------
both_metrics2plot %>% 
  left_join(both_metrics2plot_dataset_labels, by = c("CellGroup", "dataset")) %>%
  ggplot(aes(x = dataset_label, y = nonzero_fraction, color = type, fill = dataset)) +
  geom_col(position = position_dodge2(preserve = "single")) +
  facet_grid2(GeneID ~ CellGroup, axes = "x", remove_labels = "x", scales = "free_x") +
  scale_color_manual(
    values = c("gray80", "gray40", "black"), 
    labels = c("Raw", "Imputed\n(total dataset)", "Imputed\n(cell type populations)"),
    name = "") +
  scale_fill_manual(values = cols_dataset, breaks = names(cols_dataset), name = "Dataset") +
  guides(color = guide_legend(override.aes = list(fill=NA), order = 1)) +
  labs(x = "", y = "Non-zero fraction") +
  theme(
    strip.text.x = element_text(size = rel(0.88)),
    strip.text.y = element_text(face = "italic"),
    axis.text.x = element_markdown(angle = 45, hjust = 1))

ggsave(file.path(results_dir, paste0("barplot_vefs_full_ncellsx.png")), device = agg_png, width = 14, height = 9, units = "cm", dpi = 300, scaling = 0.7)
ggsave(file.path(results_dir, paste0("barplot_vefs_full_ncellsx.pdf")), width = 14, height = 9, units = "cm", scale = 1/0.7)


