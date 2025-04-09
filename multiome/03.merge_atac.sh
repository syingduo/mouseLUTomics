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

# step1
Rscript /diskmnt/Projects/Users/s.yingduo/05gudmap/02multiome/script/merge_all_peaks_from_individual_objects_v1.1.0.R \
    --sample_info_path=${sample_info_path} \
    --merging_name=${merging_name} \
    --output_path=${output_path} 1>logs/merge_peaks.log 2>&1

# step2
Rscript /diskmnt/Projects/Users/s.yingduo/05gudmap/02multiome/script/atac_pipeline_v1.1.1.R \
    --sample_info_path=${sample_info_path} \
    --merging_name=${merging_name} \
    --output_path=${output_path} \
    --pc_first=2 1>logs/atac_processing.log 2>&1
