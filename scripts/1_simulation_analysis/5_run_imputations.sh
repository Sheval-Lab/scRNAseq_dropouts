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
data_dir="../../results/1_simulation_analysis/gene_expression_filtered"
results_dir="../../results/1_simulation_analysis/gene_expression_imputed"

# Run imputation algorithms
for sim_celltype in symsim_basal symsim_club symsim_ciliated zinb-wave_basal zinb-wave_club zinb-wave_ciliated
do
  raw_counts_file=${data_dir}/${sim_celltype}_obs_counts.raw.txt.gz
  lognorm_counts_file=${data_dir}/${sim_celltype}_obs_counts.norm.txt.gz
  imputed_counts_file_prefix=${results_dir}/${sim_celltype}
  
  echo "Start ${sim_celltype} - ALRA"
  Rscript $alra_script $lognorm_counts_file $imputed_counts_file_prefix

  echo "Start ${sim_celltype} - MAGIC"
  Rscript $magic_script $lognorm_counts_file 3 $imputed_counts_file_prefix

  echo "Start ${sim_celltype} - kNN-smoothing"
  python3 $knn_script -k 10 -f ${raw_counts_file} -o ${imputed_counts_file_prefix}.knn_10.raw.txt
  Rscript $knn_script_postprocess ${imputed_counts_file_prefix}.knn_10.raw.txt 10 $imputed_counts_file_prefix
  rm ${imputed_counts_file_prefix}.knn_10.raw.txt
  
  echo "Start ${sim_celltype} - SAVER"
  Rscript $saver_script $raw_counts_file $imputed_counts_file_prefix
  
  echo "Start ${sim_celltype} - scImpute"
  Rscript $scimpute_script $raw_counts_file $results_dir/${sim_celltype}_scimpute_0.5/ 0.5 $imputed_counts_file_prefix
  rm -r $results_dir/${sim_celltype}_scimpute_0.5/

  echo "Start ${sim_celltype} - scVI"
  Rscript $scvi_script ${raw_counts_file} ${imputed_counts_file_prefix}.scvi.txt.gz
  
  echo "Start ${sim_celltype} - scBiG"
  python $scbig_script ${raw_counts_file} ${imputed_counts_file_prefix}.scbig.txt.gz

  echo "Start ${sim_celltype} - DCA"
  dca --threads 6 ${raw_counts_file} ${imputed_counts_file_prefix}.dca_nb_res_dir
  Rscript $dca_script_postprocess $imputed_counts_file_prefix nb
  rm -r ${imputed_counts_file_prefix}.dca_nb_res_dir
  
  echo "End ${sim_celltype}"
done

