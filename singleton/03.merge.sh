#!/bin/bash
eval "$(conda shell.bash hook)"
conda activate /diskmnt/Projects/Users/rliu/Software/miniconda/envs/Seurat_RNA_ATAC_ST_20220615

Rscript /diskmnt/Projects/Users/s.yingduo/05gudmap/01singleton/merging.R \
    --merging_name=GM_P21 \
    --sample_path_table=/diskmnt/Projects/Users/s.yingduo/05gudmap/01singleton/sample_info_summary.txt \
    --output=/diskmnt/Projects/Users/s.yingduo/05gudmap/01singleton/two_samples_P21/ \
    --marker.toplot=/diskmnt/Projects/Users/rliu/Projects/PKD/celltype_markers/KidneyCellMarkers.txt
