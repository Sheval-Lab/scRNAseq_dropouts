#!/bin/bash 


# Scripts to run imputation methods
alra_script="../0_imputation_scripts/run_ALRA.R"
magic_script="../0_imputation_scripts/run_MAGIC.R"
knn_script="../0_imputation_scripts/knn_smooth.py"
knn_script_postprocess="../0_imputation_scripts/run_kNNsmoothing_postprocess.R"
saver_script="../0_imputation_scripts/run_SAVER.R"
scimpute_script="../0_imputation_scripts/run_scImpute.R"
scvi_script="../0_imputation_scripts/run_scVI.R"

# Input and output directories
data_dir="../../results/4_cftr_analysis/gene_expression_filtered"
results_dir="../../results/4_cftr_analysis/gene_expression_imputed"

# Run imputation algorithms
for dataset in all_both_dropseq Basal_both_dropseq Ciliated_both_dropseq Ionocyte_both_dropseq Secretory_both_dropseq all_LAE_dropseq Basal_LAE_dropseq Ciliated_LAE_dropseq Ionocyte_LAE_dropseq Secretory_LAE_dropseq all_SAE_dropseq Basal_SAE_dropseq Ciliated_SAE_dropseq Ionocyte_SAE_dropseq Secretory_SAE_dropseq all_both_10x Basal_both_10x Ciliated_both_10x Ionocyte_both_10x Secretory_both_10x all_LAE_10x Basal_LAE_10x Ciliated_LAE_10x Ionocyte_LAE_10x Secretory_LAE_10x all_SAE_10x Basal_SAE_10x Ciliated_SAE_10x Ionocyte_SAE_10x Secretory_SAE_10x
do
  raw_counts_file=${data_dir}/${dataset}.raw.txt.gz
  lognorm_counts_file=${data_dir}/${dataset}.norm.txt.gz
  imputed_counts_file_prefix=${results_dir}/${dataset}
  
  echo "Start ${dataset} - ALRA"
  Rscript $alra_script $lognorm_counts_file $imputed_counts_file_prefix

  echo "End ${dataset}"
done

