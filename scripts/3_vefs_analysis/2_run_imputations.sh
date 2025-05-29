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
data_dir="../../results/3_vefs_analysis/gene_expression_filtered"
results_dir="../../results/3_vefs_analysis/gene_expression_imputed"

# Run imputation algorithms
for dataset in Total Basal Ciliated Club 
do
  raw_counts_file=${data_dir}/${dataset}.raw.txt.gz
  lognorm_counts_file=${data_dir}/${dataset}.norm.txt.gz
  imputed_counts_file_prefix=${results_dir}/${dataset}
  
  echo "Start ${dataset} - ALRA"
  Rscript $alra_script $lognorm_counts_file $imputed_counts_file_prefix

  echo "End ${dataset}"
done

