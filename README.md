# Low expressed RNA analysis using single cell RNA-sequencing data 

The repository contains the code to reproduce all analyses and figures of the paper "Evaluating imputation methods for accurate estimation of cell population fractions in single-cell RNA sequencing" by Valyaeva A.A., Tochilkina M.S., and Sheval E.V. 



## Scripts 

### Auxiliary scripts for running imputation algorithms and collecting metrics

- `scripts/functions.R` contains functions to calculate metrics: sample size (number of cells and genes), dropout rate per sample, nonzero fraction per gene (fraction of cells expressing the gene), mean normalized expression per gene, and collects this metrics across the samples to store them in a list of data frames.

- `scripts/0_imputation_scripts` directory contains scripts to run scRNA-seq imputation methods ALRA (`run_ALRA.R`), DCA (`run_DCA_postprocess.R`), kNN-smoothing (`knn_smooth.py` and `run_kNNsmoothing_postprocess.R`), MAGIC (`run_MAGIC.R`), SAVER (`run_SAVER.R`), scBiG (`run_scBiG.py`), scImpute (`run_scImpute.R`), and scVI (`run_scVI.R`). DCA is executed from command line.


### Testing imputation algorithms on simulated scRNA-seq data

The scripts used in this part of the anlysis are stored in `scripts/1_simulation_analysis` directory.

- `1_symsim_find_best_params.R` and `2_symsim_generate_datasets.R` runs [SimSym](https://doi.org/10.1038/s41467-019-10500-w) to estimate parameters for data simulation from real scRNA-seq dataset of mouse tracheal epithelial cells [(Montoro et al. 2018)](https://doi.org/10.1038/s41586-018-0393-7) and generates 3 synthetics datasets based on the populations of basal, club, and ciliated cells.

- `3_zinb-wave_generate_datasets.R` runs [ZINB-WaVE](https://doi.org/10.1038/s41467-017-02554-5) to estimate parameters for data simulation from real scRNA-seq dataset of mouse tracheal epithelial cells [(Montoro et al. 2018)](https://doi.org/10.1038/s41586-018-0393-7) and generates 3 synthetic datasets based on the populations of basal, club, and ciliated cells.

- `4_preprocessing.R` performs filtering of non-expressed genes (min.cells = 3) and poor-quality cells (min.features = 200) and exports datasets as Seurat objects and text files with un-normalized and log-normalized gene expression counts. 

- `5_run_imputations.sh` applies all 6 imputation methods to synthetic datasets. 

- `6_collect_metrics.R` performs aggregation of per-sample metrics of the fraction of cells with non-zero expression of a gene and mean normalized expression of a gene before and after scRNA-seq data imputation. It also computes NRMSE values for observed and imputed count matrices.

- `7_plot_nonzero_fraction_distr.R`, `8_plot_delta_nrmse.R`, and `9_plot_nonzero_fraction_after_imputation.R` are used to make plots: distributions of nonzero fraction values in the reference and simulated datasets, scatter plots of ΔNRMSE *versus* observed nonzero fraction, and scatter plots of nonzero fraction after imputation *versus* true nonzero fraction.

- `10_zero_threshold.R` explores the impact of different non-zero expression thresholds on 5 imputation methods (MAGIC, SAVER, scVI, DCA, and scBiG) performance in the recovery of non-zero fraction. 

- `11_check_nb_or_zinb.R` fits ZINB-WaVE-generated data to NB and ZINB distribution to find out what model describes the synthetic data better.


### Testing imputation algorithms on real scRNA-seq data of K562 cells expressing dCas9-KRAB (GEO ID: [GSE204750](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE204750))

The scripts used in this part of the anlysis are stored in `scripts/2_k562_analysis` directory.

- `1_downsampling.sh` performs downsampling of K562 cells dataset. As input the script requires a number [0:1], which corresponds to reads proportion to sample (-p parameter of seqkit sample command). Multiple runs of downsampling with varying -p parameter settings can be run in parallel: `parallel bash sampling.sh ::: 0.02 0.04 0.06`.

- `2_cellranger.sh` runs Cell Ranger pipeline (10x Genomics) that includes aligning reads to human reference transcriptome (version GRCh38-2020-A) and gene expression quantification. As input the script requires proportion value, used to downsample the full dataset. All samples can be processed in a cycle: `for i in 0.02 0.04 0.06 0.08 0.10 0.12 0.14 0.16 0.18 0.20 0.30 0.40 0.50 0.60 0.80; do bash 2_cellranger.sh $i; done`. 

- `3_preprocessing.R` performs filtering of non-expressed genes (min.cells = 3) and poor-quality cells (min.features = 200) and exports datasets as Seurat objects and text files with un-normalized and log-normalized gene expression counts. 

- `4_collect_raw_metrics.R` performs aggregation of Cell Ranger outputs metrics (e.g. estimated number of cells, mean reads per cell, and median genes per cell) and per-sample metrics of the fraction of cells with non-zero expression of a gene and mean normalized expression of a gene. Creates a plot showing median genes per cell ~ mean reads per cell dependency.

- `5_run_imputations.sh` applies all 6 imputation methods to downsampled datasets. 

- `6_collect_imputed_metrics.R` performs aggregation of per-sample metrics of the fraction of cells with non-zero expression of a gene and mean normalized expression of a gene after scRNA-seq data imputation. 

- `7_K562_analysis.R` visualizes the results of imputation of downsampled datasets and performs the comparison between imputed and original deep-sequenced data. 


### Example of adjusting gene expression of SARS-CoV-2 entry factors (*Ace2*, *Tmprss2*, *Furin*, and *Ctsl*) in populations of mouse tracheal epithelial cells [(Montoro et al. 2018)](https://doi.org/10.1038/s41586-018-0393-7) (GEO ID: [GSE103354](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE103354)) and *CFTR* gene expression in populations of human airway epithelial cells [(Okuda et al. 2021)](https://doi.org/10.1164/rccm.202008-3198oc) (GEO IDs: [GSE160673](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE160673) and [GSE160664](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE160664))

The scripts used in this part of the anlysis are stored in `scripts/3_vefs_analysis`  and `scripts/4_cftr_analysis` directories. Both directories contain the same set of four scripts:

- `1_preprocessing.R` collects datasets' basic QC metrics, splits datasets by cell populations, calculates nonzero fraction values for the genes of interest, performs filtering of non-expressed genes (min.cells = 3) and poor-quality cells (min.features = 200) and exports datasets as text files with un-normalized and log-normalized gene expression counts. 

- `2_run_imputations.sh` applies ALRA imputation methods to datasets. 

- `3_collect_metrics.R` performs aggregation of per-sample metrics of the fraction of cells with non-zero expression of a gene and mean normalized expression of a gene before and after scRNA-seq data imputation. It also calculates nonzero fraction and average expression for the genes of interest for every cell population in full datasets (unsplitted by cell populations).

- `4_plot_nonzero_fraction.R` is used to make barplots illustrating the comparison of the reference, raw, and imputed values of nonzero fraction for the genes of interest.
 

## Data

Synthetic scRNA-seq data generated in this study are available on Zenodo: https://doi.org/10.5281/zenodo.15546550. All other datasets used in the study are publicly available in [NCBI GEO database](https://www.ncbi.nlm.nih.gov/geo/). 