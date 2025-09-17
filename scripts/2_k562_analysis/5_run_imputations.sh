#!/bin/bash 


# Scripts to run imputation methods
alra_script="../0_imputation_scripts/run_ALRA.R"
magic_script="../0_imputation_scripts/run_MAGIC.R"
knn_script="../0_imputation_scripts/knn_smooth.py"
knn_script_postprocess="../0_imputation_scripts/run_kNNsmoothing_postprocess.R"
saver_script="../0_imputation_scripts/run_SAVER.R"
scimpute_script="../0_imputation_scripts/run_scImpute.R"
scvi_script="../0_imputation_scripts/run_scVI.R"
scbig_script="../0_imputation_scripts/run_scBiG.py"
dca_script_postprocess="../0_imputation_scripts/run_DCA_postprocess.R"

# Input and output directories
data_dir="../../results/2_k562_analysis/gene_expression_filtered"
results_dir="../../results/2_k562_analysis/gene_expression_imputed"

# Run imputation algorithms
for prop in 002 004 006 008 010 012 014 016 018 020 030 040 050 060 080
do
  raw_counts_file=${data_dir}/sample_${prop}.raw_counts.txt.gz
  lognorm_counts_file=${data_dir}/sample_${prop}.norm_counts.txt.gz
  imputed_counts_file_prefix=${results_dir}/sample_${prop}
  
  echo "Start sample_${prop} - ALRA"
  Rscript $alra_script $lognorm_counts_file $imputed_counts_file_prefix
  
  echo "Start sample_${prop} - MAGIC"
  for t in 1 2 3 5 7 10
  do
    Rscript $magic_script $lognorm_counts_file $t $imputed_counts_file_prefix
  done
  
  #echo "Start sample_${prop} - kNN-smoothing"
  for k in 2 3 4 6 8 10 16 32 
  do
    python3 $knn_script -k $k -f ${raw_counts_file} -o ${imputed_counts_file_prefix}.knn_${k}.raw.txt
    Rscript $knn_script_postprocess ${imputed_counts_file_prefix}.knn_${k}.raw.txt $k $imputed_counts_file_prefix
    rm ${imputed_counts_file_prefix}.knn_${k}.raw.txt
  done 
  
  echo "Start sample_${prop} - SAVER"
  Rscript $saver_script $raw_counts_file $imputed_counts_file_prefix
  
  #echo "Start sample_${prop} - scImpute"
  for t in 0.3 0.5 0.9 
  do
    Rscript $scimpute_script $raw_counts_file $results_dir/sample_${prop}_scimpute_${t}/ $t $imputed_counts_file_prefix
    rm -r $results_dir/sample_${prop}_scimpute_${t}/
  done
  
  echo "Start sample_${prop} - scVI"
  Rscript $scvi_script ${raw_counts_file} ${imputed_counts_file_prefix}.scvi.txt.gz
  
  echo "Start sample_${prop} - scBiG"
  python $scbig_script ${raw_counts_file} ${imputed_counts_file_prefix}.scbig.txt.gz
  
  echo "Start sample_${prop} - DCA"
  dca --threads 6 ${raw_counts_file} ${imputed_counts_file_prefix}.dca_nb_res_dir # default nb-conddisp model
  Rscript $dca_script_postprocess $imputed_counts_file_prefix nb
  # rm -r ${imputed_counts_file_prefix}.dca_nb_res_dir
  
  echo "End sample_${prop}"
done

