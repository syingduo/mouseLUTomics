## 20230214; changed the sample annotation column from rds_object to Path (so that the sample input information from snRNA and ATAC merging pipeline would be consistent)

#### Ruiyang Liu, 20230124
# Adopted from: /diskmnt/Projects/Users/rliu/Projects/ccRCC/cell_of_origin_study/snATAC/snATAC_sample_merging_20220801/scripts/Merge_ATAC_samples_auto.v.3.2.R

library(optparse)

option_list <- list(
    make_option(c("-s", "--sample_info_path"), ## /diskmnt/Projects/Users/rliu/Projects/normal_mouse_kidney_development/snATAC/Merging/NMK3_and_old_male/sample_info.txt
        type = "character",
        default = NULL,
        help = "samples_for_merging,use \",\" to separate multiple input", # e.g. HTHK1301462,HTHK1301463
        metavar = "character"
    ),
    make_option(c("-M", "--merging_name"), ## case_id <- "2418"
        type = "character",
        default = "",
        help = "merging  name",
        metavar = "character"
    ),
    make_option(c("-o", "--output_path"),
        type = "character",
        default = "./",
        help = "path to output directory",
        metavar = "characer"
    ),
    make_option(c("--pc_num"),
        type = "integer",
        default = 30,
        help = "number of principal components to use",
        metavar = "integer"
    ),
    make_option(c("--pc_first"),
        type = "integer",
        default = 1,
        help = "first principal components to use (should be 1 or 2)",
        metavar = "integer"
    )
)



library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(ggplot2)
library(RColorBrewer)

require(magrittr)
require(readr)
require(Matrix)
require(tidyr)
require(dplyr)
set.seed(1234)

library(plyr)
library(tibble)
library(data.table)
library(reshape)

library(BSgenome)
library(BSgenome.Mmusculus.UCSC.mm10)
library(EnsDb.Mmusculus.v79)

library(future)

### get input parameters
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)


## specify the additional options:
# opt$sample_info_path <- "/diskmnt/Projects/Users/rliu/Projects/normal_mouse_kidney_development/snRNA_ATAC_combo/merging_ATAC/test_NMK_W3MX3/sample_rds_path.txt"
# opt$merging_name <- "test_NMK_W3MX3"
# opt$output_path <- "/diskmnt/Projects/Users/rliu/Projects/normal_mouse_kidney_development/snRNA_ATAC_combo/merging_ATAC/test_NMK_W3MX3/"
# opt$pc_first <- 2

merging <- opt$merging_name

out_path <- opt$output_path
if (!file.exists(out_path)) {
    dir.create(out_path)
}
if (!file.exists(paste0(out_path, "/peaks"))) {
    dir.create(paste0(out_path, "/peaks"))
}

sample_info <- opt$sample_info_path %>% read.table(head = TRUE, sep = "\t")


### some parallelization-solution from the tutorial:
plan("multicore", workers = 4)
options(future.globals.maxSize = 200 * 1024^3) # for 100 Gb RAM

samples <- sample_info$sample_id %>% as.character()

if (!file.exists(paste0(out_path, "peaks/", merging, "_snATAC_Merged.rds.gz"))) {
    atac <- vector(mode = "list", length = length(samples))
    for (i in 1:nrow(sample_info)) {
        cat("read in ", sample_info[i, "sample_id"] %>% as.character(), "...\n")
        atac[[i]] <- readRDS(sample_info[i, "Path"] %>% as.character())
        DefaultAssay(atac[[i]]) <- "X500peaksMACS2"
        # DefaultAssay(atac[[i]]) <- "ATAC_MACS2"

        for (assay in Assays(atac[[i]])) {
            if (assay != "X500peaksMACS2") {
                # if (assay != "ATAC_MACS2") {
                atac[[i]][[assay]] <- NULL
            }
        }
    }

    ##### To obtain the best results - use ALL peaks!
    combined.peaks <- UnifyPeaks(object.list = atac, mode = "reduce")
    peakwidths <- width(combined.peaks)
    combined.peaks <- combined.peaks[peakwidths < 10000 & peakwidths > 20]
    combined.peaks
    peaks.use <- combined.peaks


    # We don't filter cells like in the tutorial, because we use already filtered matrices. And all cells are pass those filters in the tutorial.

    matrix.counts <- vector(mode = "list", length = length(samples))

    for (i in 1:length(samples)) {
        matrix.counts[[i]] <- FeatureMatrix(
            fragments = Fragments(atac[[i]]@assays$X500peaksMACS2),
            # fragments = Fragments(atac[[i]]@assays$ATAC_MACS2),
            features = peaks.use,
            sep = c("-", "-"),
            cells = colnames(atac[[i]])
        )
    }

    for (i in 1:length(samples)) {
        atac[[i]][["peaksinters"]] <- CreateChromatinAssay(counts = matrix.counts[[i]], fragments = Fragments(atac[[i]]@assays$X500peaksMACS2))
        # atac[[i]][["peaksinters"]] <- CreateChromatinAssay(counts = matrix.counts[[i]], fragments = Fragments(atac[[i]]@assays$ATAC_MACS2))
        atac[[i]]$dataset <- samples[i]
        DefaultAssay(atac[[i]]) <- "peaksinters"
        atac[[i]]$X500peaksMACS2 <- NULL
        # atac[[i]]$ATAC_MACS2 <- NULL
    }

    #### Merging:
    combined <- merge(x = atac[[1]], y = atac[2:length(samples)], add.cell.ids = samples)
    saveRDS(combined, paste0(out_path, "peaks/", merging, "_snATAC_Merged.rds.gz"), compress = TRUE)
    rm(combined)
    rm(atac)
    gc()
}


##########################################################
#### Now create Feature Matrix with the new set of peaks:###
#### new peaks are generated based on the script here: /diskmnt/Projects/Users/rliu/Projects/normal_mouse_kidney_development/snRNA_ATAC_combo/merging_ATAC/src/merge_all_peaks_from_individual_objects_v1.0.R
##########################################################
recentered_final <- read.table(paste0(out_path, "peaks/recentered_final.filtered.tsv"), sep = "\t", header = T)

combined <- readRDS(paste0(out_path, "peaks/", merging, "_snATAC_Merged.rds.gz"))
atac <- combined

recentered_p <- StringToGRanges(recentered_final$new_peak, sep = c("-", "-"))
matrix.counts <- FeatureMatrix(
    fragments = Fragments(atac@assays$peaksinters),
    features = recentered_p,
    sep = c("-", "-"),
    cells = colnames(atac)
)

### next is here:
atac[["peaksMACS2"]] <- CreateChromatinAssay(
    counts = matrix.counts,
    fragments = Fragments(atac@assays$peaksinters)
)
DefaultAssay(atac) <- "peaksMACS2"

# atac[['X500peaksMACS2']]<-NULL
atac[["peaksinters"]] <- NULL

### Overlapping ranges supplied -- check this later:


#### filter mt genes
atac <- subset(atac, percent.mito <= 1)

### Remove some assays
saveRDS(atac, paste0(out_path, "peaks/", merging, "_snATAC_Merged_BasedOnSelectedPeaks.rds.gz"), compress = TRUE)


###################
# ends here for now#
###################

atac <- FindTopFeatures(atac, min.cutoff = "q0")
atac <- RunTFIDF(atac)
atac <- RunSVD(atac)
atac <- RunUMAP(object = atac, reduction = "lsi", dims = opt$pc_first:50)
atac <- FindNeighbors(object = atac, reduction = "lsi", dims = opt$pc_first:50)

### this helps:
options(future.globals.maxSize = 891289600, future.seed = TRUE)
atac <- FindClusters(object = atac, verbose = FALSE, algorithm = 3)
saveRDS(atac, paste0(out_path, "peaks/", merging, "_snATAC_Merged_BasedOnSelectedPeaks_Normalized.rds.gz"), compress = TRUE)


### Add ATAC Gene activity
options(future.globals.maxSize = 200 * 1024^3) # for 100 Gb RAM
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
seqlevelsStyle(annotations) <- "UCSC"
genome(annotations) <- "mm10"

Annotation(atac) <- annotations

gene.activities <- GeneActivity(atac)
atac[["ATACGeneActivity"]] <- CreateAssayObject(counts = gene.activities)

atac <- NormalizeData(object = atac, assay = "ATACGeneActivity", normalization.method = "LogNormalize", scale.factor = median(atac$nCount_RNA))

saveRDS(atac, paste0(out_path, "peaks/", merging, "_snATAC_Merged_BasedOnSelectedPeaks_Normalized_With_GeneActivity.rds.gz"), compress = TRUE)
