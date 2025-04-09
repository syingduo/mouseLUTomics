#!/bin/bash

## enviroment and software
# export TMPDIR=/diskmnt/Projects/Users/rliu/Projects/TMP
eval "$(conda shell.bash hook)"
conda activate /diskmnt/Projects/Users/rliu/Software/miniconda/envs/Seurat_RNA_ATAC_ST_20220615
macs2_path=/diskmnt/Projects/Users/rliu/Software/miniconda/envs/Seurat_RNA_ATAC_ST/bin/macs2
pipeline=/diskmnt/Projects/Users/rliu/Projects/GUDMAP/Seurat_combo/src/multiome_pipeline_v1.2_mouse.R
atac_data=/diskmnt/Projects/GUDMAP_analysis/combo_snRNA_snATAC

## run
WORKING_DIR=/diskmnt/Projects/Users/s.yingduo/05gudmap/02multiome/snRNA_snATAC
cd $WORKING_DIR
rm -rf logs && mkdir logs

Rscript $pipeline --sample=GM005MwP21 --atac_data=$atac_data --macs2_path=$macs2_path --prf_min=1000 >>logs/GM005MwP21.log 2>&1 ## 20230424
Rscript $pipeline --sample=GM006MwP21 --atac_data=$atac_data --macs2_path=$macs2_path --prf_min=1000 >>logs/GM006MwP21.log 2>&1 ## 20230424
