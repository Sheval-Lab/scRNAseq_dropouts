#!/bin/bash 
# Reads proportion to sample: 0.02 0.04 0.06 0.08 0.10 0.12 0.14 0.16 0.18 0.20 0.30 0.40 0.50 0.60 0.80
prop=$1

data_dir="../../data/2_k562_dataset/fastqs"
cd $data_dir
mkdir sample_${prop}

seqkit sample -p ${prop} -s 10 -2 sample_1.00/K562_S1_L001_R1_001.fastq.gz -o sample_${prop}/K562_S1_L001_R1_001.fastq.gz
seqkit sample -p ${prop} -s 10 -2 sample_1.00/K562_S1_L001_R2_001.fastq.gz -o sample_${prop}/K562_S1_L001_R2_001.fastq.gz