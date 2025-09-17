# Based on https://divingintogeneticsandgenomics.com/post/negative-bionomial-distribution-in-single-cell-rnaseq/

library(tidyverse)
library(MASS)
library(pscl)
library(ragg)
library(patchwork)
library(scales)

theme_set(
  theme_classic() +
    theme(
      panel.spacing.x = unit(1, "lines"),
      axis.title = element_text(size = 10), 
      axis.text = element_text(color = "black", size = 7),
      axis.ticks = element_line(color = "black"),
      strip.text = element_text(size = 8)))


data_dir <- "results/1_simulation_analysis/gene_expression_filtered"
results_dir <- "results/1_simulation_analysis/metrics"

# Load synthetic 'observed' data ----------------------------------------------- 
## ZINB-WaVE
zw_obs <- read.table(file.path(data_dir, "zinb-wave_basal_obs_counts.raw.txt.gz"), header = TRUE, row.names = 1)


# Try fit to NB and inspect zero inflation -------------------------------------
zw_obs_summary <- data.frame(
  gene = rownames(zw_obs),
  mean_expr = rowMeans(zw_obs),
  gene_var = rowVars(as.matrix(zw_obs)),
  dropout_rate = rowMeans(zw_obs == 0))

ggplot(zw_obs_summary, aes(x = mean_expr, y = dropout_rate)) +
  geom_point(alpha = 0.6) +
  scale_x_log10(labels = label_log(digits = 1)) +
  xlab("Mean expression") +
  ylab("Dropout rate") +
  theme_minimal() 

zw_obs_results <- data.frame(
  gene = character(),
  nb_loglik = numeric(),
  zinb_loglik = numeric(),
  lrt_stat = numeric(),
  p_value = numeric(),
  stringsAsFactors = FALSE)

for (i in 1:nrow(zw_obs)) {
  gene_counts <- as.matrix(zw_obs)[i, ]
  gene_name <- rownames(zw_obs)[i]
  print(gene_name)
  
  # Fit NB model
  nb_fit <- try(glm.nb(gene_counts ~ 1), silent = TRUE)
  if (inherits(nb_fit, "try-error")) next
  
  # Fit ZINB model
  zinb_fit <- try(zeroinfl(gene_counts ~ 1 | 1, dist = "negbin"), silent = TRUE)
  if (inherits(zinb_fit, "try-error")) next
  
  # Extract log-likelihoods
  nb_ll <- logLik(nb_fit)
  zinb_ll <- logLik(zinb_fit)
  
  # LRT: ZINB has 1 extra parameter (zero-inflation probability)
  lrt_stat <- 2 * (zinb_ll - nb_ll)
  p_val <- pchisq(as.numeric(lrt_stat), df = 1, lower.tail = FALSE)
  
  zw_obs_results <- rbind(zw_obs_results, data.frame(
    gene = gene_name,
    nb_loglik = as.numeric(nb_ll),
    zinb_loglik = as.numeric(zinb_ll),
    lrt_stat = as.numeric(lrt_stat),
    p_value = p_val
  ))
}

zw_obs_results$adj_pval <- p.adjust(zw_obs_results$p_value, method = "BH")
zw_obs_results$significant_zi <- zw_obs_results$adj_pval < 0.05

zw_zi_results <- left_join(zw_obs_summary, zw_obs_results, by = "gene")
write_tsv(zw_zi_results, file.path(results_dir, "zw_obs_zero_inflation.txt"))


## Get NB and ZINB models parameters -------------------------------------------
zw_zi_results <- read_tsv(file.path(results_dir, "zw_obs_zero_inflation.txt"))


### NB -------------------------------------------------------------------------
model <- lm(gene_var ~  1* mean_expr + I(mean_expr^2) + 0, data = zw_zi_results)
summary(model)

size_nb <- 1 / coef(model)

zw_zi_results <- zw_zi_results %>% 
  mutate(zeros_nb = (size_nb / (mean_expr + size_nb)) ^ size_nb)


### ZINB -----------------------------------------------------------------------
zw_zi_results <- zw_zi_results %>% 
  mutate(
    psi_hat = pmax(0, pmin(1, (dropout_rate - zeros_nb) / (1 - zeros_nb))),
    psi_global = median(psi_hat, na.rm = TRUE),
    psi_global = pmin(psi_global, 0.95),
    zeros_zinb = psi_global + (1 - psi_global) * zeros_nb)

## Plot ------------------------------------------------------------------------
zw_zi_results %>% 
  filter(!is.na(significant_zi)) %>% 
  ggplot(aes(x = mean_expr, y = dropout_rate)) +
  geom_point(aes(color = significant_zi), alpha = 0.6, size = 1) +
  geom_line(aes(y = zeros_nb), color = "red3", linetype = "solid", size = 1) +
  scale_x_log10(labels = label_log(digits = 1)) +
  scale_color_manual(
    values = c("FALSE" = "black", "TRUE" = "dodgerblue3"),
    breaks = c("TRUE", "FALSE"),
    labels = c("Zero inflation", "No zero inflation")) +
  labs(x = "Mean expression", y = "Dropout rate", color = "") +
  theme(
    aspect.ratio = 2/3,
    legend.position = "top")

ggsave(file.path(results_dir, "zw_obs_zero_inflation.png"), device = agg_png, width = 7, height = 5.5, units = "cm", dpi = 300, scaling = 0.7) 



