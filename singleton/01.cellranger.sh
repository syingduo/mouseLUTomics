#!/bin/bash

## run
/diskmnt/Projects/Users/rliu/Software/cellranger-7.1.0/cellranger count \
    --id=GPAD-GM006MwP21-MuVN1ZthV1_1Bn1_1 \
    --fastqs=/diskmnt/Projects/GUDMAP_analysis/FASTQ/GPAD-GM006MwP21-MuVN1ZthV1_1Bn1_1 \
    --transcriptome=/diskmnt/Datasets/Reference/Cellranger-2020-A/refdata-gex-mm10-2020-A \
    --chemistry=threeprime \
    --jobmode=local \
    --localcores=40 \
    --localmem=128 \
    --include-introns=true >logs/GPAD-GM006MwP21-MuVN1ZthV1_1Bn1_1.worklog
