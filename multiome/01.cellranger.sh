#!/bin/bash

mkdir logs
/diskmnt/Projects/Users/rliu/Software/cellranger-arc-2.0.0/cellranger-arc count \
    --id=20220130-E16_5F \
    --reference=/diskmnt/Datasets/Reference/Cellranger-ARC/refdata-cellranger-arc-mm10-2020-A-2.0.0 \
    --libraries=/diskmnt/Projects/GUDMAP_analysis/combo_snRNA_snATAC/FASTQ/20220130-E16_5F/library.csv \
    --localcores=48 \
    --localmem=128 >logs/20220130-E16_5F.log
