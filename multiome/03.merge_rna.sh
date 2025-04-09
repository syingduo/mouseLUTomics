#!/bin/bash

## enviroment and software
#export TMPDIR=/diskmnt/Projects/Users/rliu/TMP/

eval "$(conda shell.bash hook)"
conda activate /diskmnt/Projects/Users/s.yingduo/miniconda3/envs/scrna

## run
merging_name=all_18_samples_new
sample_info_path=/diskmnt/Projects/Users/s.yingduo/05gudmap/02multiome/LUT/${merging_name}/sample_rds_path.txt
output_path=/diskmnt/Projects/Users/s.yingduo/05gudmap/02multiome/LUT/${merging_name}/

cd ${output_path}
mkdir logs

Rscript /diskmnt/Projects/Users/s.yingduo/05gudmap/02multiome/script/merging_conserve.memory_sample_id_not_in_regression.R \
    --merging_name=${merging_name} \
    --sample_path_table=${sample_info_path} \
    --output=${output_path} \
    --marker.toplot=NULL 1>logs/merge_rna.log 2>&1
