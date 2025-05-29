#!/bin/bash 
prop=$1

# Navigate to working directory 
out_dir="../../data/2_k562_dataset/cellranger_outs"
cd $out_dir

# Relative to out_dir
ref="../refdata-gex-GRCh38-2020-A"

cellranger count --id=sample_$(echo ${prop} | sed -e 's/\.//g') --transcriptome=$ref --fastqs=../fastqs/sample_${prop} --sample=K562 --localcores=8 --localmem=64
