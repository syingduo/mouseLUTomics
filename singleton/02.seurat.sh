#!/bin/bash
eval "$(conda shell.bash hook)"
conda activate /diskmnt/Projects/Users/rliu/Software/miniconda/envs/Seurat_RNA_ATAC_ST_20220615

Rscript /diskmnt/Projects/Users/rliu/Projects/GUDMAP/Seurat/Seurat_v4_analysis_auto_v0.4.R \
    --input=/diskmnt/Projects/Users/s.yingduo/05gudmap/01singleton/GPAD-GM006MwP21-MuVN1ZthV1_1Bn1_1/outs/raw_feature_bc_matrix \
    --pre_filter=300 \
    --nfeature_min=1000 --nfeature_max=7500 --ncount_min=1000 --ncount_max=40000 --mito_max=0.1 \
    --output=/diskmnt/Projects/Users/s.yingduo/05gudmap/01singleton/GPAD-GM006MwP21-MuVN1ZthV1_1Bn1_1/nfeature_1000_7500_ncount_1000_40000_mito_0.1/ \
    --sample_id=GPAD-GM006MwP21-MuVN1ZthV1_1Bn1_1 \
    --project=GUDMAP
