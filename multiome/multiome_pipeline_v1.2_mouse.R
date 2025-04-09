### 2021-10-19, Nadezhda V. Terekhanova: read here for the changes made/diskmnt/Projects/HTAN_analysis/snATAC/BR_HTAN/Signac.v.1.3/Scripts/Multiome/README

# References based on
# https://satijalab.org/signac/articles/pbmc_multiomic.html

library(optparse)
set.seed(1234)
library(future)
plan("multiprocess", workers = 20)
options(future.globals.maxSize = 100 * 1024^3)

option_list <- list(
    make_option(c("-s", "--sample"),
        type = "character",
        default = NULL,
        help = "sample_name",
        metavar = "character"
    ),
    make_option(c("-a", "--atac_data"),
        type = "character",
        default = NULL,
        help = "path to data folder (e.g. cellranger output's raw matrices folder)",
        metavar = "character"
    ),
    make_option(c("-m", "--macs2_path"),
        type = "character",
        default = NULL,
        help = "path to installed MACS2",
        metavar = "character"
    ),

    # RNA QC metrix (from /diskmnt/Projects/HTAN_analysis/snRNA/Pipeline_development/Seurat_v3_analysis_auto_v0.6.R)
    make_option(c("--pre_filter"),
        type = "integer",
        default = 300,
        help = "min number of reads per cell to prefilter",
        metavar = "integer"
    ),
    make_option(c("--nfeature_min"),
        type = "integer",
        default = 200,
        help = "nFeature_RNA min value for filtering",
        metavar = "integer"
    ),
    make_option(c("--nfeature_max"),
        type = "integer",
        default = 10000,
        help = "nFeature_RNA max value for filtering",
        metavar = "integer"
    ),
    make_option(c("--ncount_min"),
        type = "integer",
        default = 1000,
        help = "nCount_RNA min value for filtering",
        metavar = "integer"
    ),
    make_option(c("--ncount_max"),
        type = "integer",
        default = 80000,
        help = "nCount_RNA max value for filtering",
        metavar = "integer"
    ),
    make_option(c("--mito_max"),
        type = "double",
        default = 0.1,
        help = "maximum allowed mitochondrial fraction",
        metavar = "double"
    ),
    make_option(c("--rna_pc_num"),
        type = "integer",
        default = 30,
        help = "number of principal components to use",
        metavar = "integer"
    ),
    # ATAC QC metrics
    make_option(c("--prf_min"),
        type = "integer",
        default = 3000,
        help = "peak_region_fragments_minimum value for filtering",
        metavar = "integer"
    ),
    make_option(c("--prf_max"),
        type = "integer",
        default = 20000,
        help = "peak_region_fragments_maximum value for filtering",
        metavar = "integer"
    ),
    make_option(c("--pct_min"),
        type = "integer",
        default = 15,
        help = "pct_reads_in_peaks_minimum value for filtering",
        metavar = "integer"
    ),
    make_option(c("--bl_ratio"),
        type = "double",
        default = 0.05,
        help = "blacklist_ratio_minimum value for filtering",
        metavar = "double"
    ),
    # Changed to default=4, based on the latest Signac-vignette
    make_option(c("--ns_max"),
        type = "integer",
        default = 4,
        help = "nucleosome_signal_maximum value for filtering",
        metavar = "integer"
    ),
    make_option(c("--tss"),
        type = "integer",
        default = 2,
        help = "tss_enrichment_minimum value for filtering",
        metavar = "integer"
    ),
    make_option(c("--atac_pc_num"),
        type = "integer",
        default = 30,
        help = "number of principal components to use",
        metavar = "integer"
    ),
    make_option(c("--atac_pc_first"),
        type = "integer",
        default = 1,
        help = "first principal components to use (should be 1 or 2)",
        metavar = "integer"
    )
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$sample) | is.null(opt$atac_data)) {
    print_help(opt_parser)
    stop("At least two arguments must be supplied (sample_name,atac_data).n", call. = FALSE)
}

## input data
sample <- opt$sample
atac_data_folder <- opt$atac_data

print("Input parameters")
print(paste("--peak_region_fragments_min:", opt$prf_min, sep = ""))
print(paste("--peak_region_fragments_max:", opt$prf_max, sep = ""))
print(paste("--pct_reads_in_peaks_minimum:", opt$pct_min, sep = ""))
print(paste("--blacklist_ratio_minimum:", opt$bl_ratio, sep = ""))
print(paste("--nucleosome_signal_maximum:", opt$ns_max, sep = ""))
print(paste("--tss_enrichment_minimum:", opt$tss, sep = ""))
# print(paste("--pc_first:",opt$pc_first,sep=""))
# print(paste("--pc_num:",opt$pc_num,sep=""))


print(paste("ATAC data:", atac_data_folder, sep = ""))

## output data
print(sample)

##### LOAD REQUIRED PACKAGES##########
library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(GenomicRanges)
library(EnsDb.Mmusculus.v79)
library(ggplot2)
library(patchwork)
library(BSgenome.Mmusculus.UCSC.mm10)

# atac_data_folder='/diskmnt/Projects/cptac_scratch/CPTAC3_analysis/CPTAC3_EC_snRNA_Analysis/Cellranger_arc/'
# sample='CPT1541DU-T1N1'


###########################
######## LOAD IN DATA#######
###########################

dir.create("out")
outputpath <- paste("out/", sample, "/", sep = "")
dir.create(outputpath)

counts <- Read10X_h5(paste(atac_data_folder, "/", sample, "/outs/raw_feature_bc_matrix.h5", sep = ""),
    use.names = TRUE, unique.features = TRUE
)
metadata <- read.csv(
    file = paste(atac_data_folder, "/", sample, "/outs/per_barcode_metrics.csv", sep = ""),
    header = TRUE, row.names = 1
)
fragment.path <- paste(atac_data_folder, "/", sample, "/outs/atac_fragments.tsv.gz", sep = "")

# get gene annotations for mm10
annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
genome(annotation) <- "NA"
seqlevelsStyle(annotation) <- "UCSC"
genome(annotation) <- "mm10"

# Use solution from the post (https://github.com/timoast/signac/issues/570), comment Apr. 22
# Create Seurat object containing chromatin assay, and add RNA assay to it
chrom_assay <- CreateChromatinAssay(
    counts = counts$Peaks,
    sep = c(":", "-"),
    fragments = fragment.path,
    annotation = annotation,
    min.features = 200,
    min.cells = -1
)

pbmc <- CreateSeuratObject(
    counts = chrom_assay,
    assay = "ATAC",
    meta.data = metadata
)


# Create RNA assay and add it to the object
pbmc[["RNA"]] <- CreateAssayObject(
    counts = counts$"Gene Expression"[, colnames(pbmc)]
)


#### 2021-03-20: Change Peaks to MACS2
###########################################################
############ MACS2 peak calling#############################
###########################################################
peaks <- CallPeaks(
    object = pbmc,
    macs2.path = opt$macs2_path
)
# remove peaks on nonstandard chromosomes and in genomic blacklist regions
peaks <- keepStandardChromosomes(peaks, pruning.mode = "coarse")
peaks <- subsetByOverlaps(x = peaks, ranges = blacklist_mm10, invert = TRUE)

p <- as.data.frame(peaks)
write.table(p, paste(outputpath, "MACS2_peaks.", sample, ".tsv", sep = ""),
    sep = "\t", quote = FALSE,
    row.names = FALSE
)

p$peak_center <- p$start + p$relative_summit_position
p$recentered_start <- p$peak_center - 250
p$recentered_end <- p$peak_center + 250

#### Now check that new start and end don't go beyond the chromosome boundaries
chr_size <- read.table("/diskmnt/Projects/Users/rliu/Projects/normal_mouse_kidney_development/Resources/genome_reference/mm10.chrom.sizes.txt", sep = "\t", header = FALSE)
colnames(chr_size) <- c("seqnames", "chr_length")
p1 <- merge(p, chr_size, all.x = TRUE)

p1 <- p1[p1$recentered_end <= p1$chr_length & p1$recentered_start >= 0, ]
p1$length <- p1$recentered_end - p1$recentered_start + 1
p1$new_peak <- paste(p1$seqnames, p1$recentered_start, p1$recentered_end, sep = "-")

#### Add step of removing peaks with "N" and peaks on chrY:
peaks_gr <- StringToGRanges(p1$new_peak, sep = c("-", "-"))
seq <- getSeq(BSgenome.Mmusculus.UCSC.mm10, peaks_gr) # extract fasta sequence
names(seq) <- p1$new_peak
peaks.match.pattern <- vmatchPattern("N", seq) # match peak sequence with N in them
peaks.withN <- names(peaks.match.pattern)[elementNROWS(peaks.match.pattern) > 0] # remove peaks with "N"
p1 <- p1[!(p1$new_peak %in% peaks.withN), ]
p1 <- p1[p1$seqnames != "chrY", ]
print('removed peaks with "N" and the ones on chrY')

##### Change row.names -- this is needed to be changed, since some peaks were removed, and we use row.names ids
rownames(p1) <- c(1:nrow(p1))

recentered_p <- StringToGRanges(p1$new_peak, sep = c("-", "-"))

olap <- as.data.frame(findOverlaps(recentered_p, recentered_p))
olap1 <- olap[olap$queryHits != olap$subjectHits, ]

recentered_non_olap <- p1[-olap1$queryHits, ]
recentered_olap <- p1[olap1$queryHits, ]


pairs <- cbind(p1[olap1$queryHits, c(1, 13, 14, 10)], olap1$queryHits, p1[olap1$subjectHits, c(1, 13, 14, 10)], olap1$subjectHits)
colnames(pairs) <- c("chr_1", "st_1", "en_1", "score_1", "row_1", "chr_2", "st_2", "en_2", "score_2", "row_2")

### Remove this step to keep weak peak in clusters of peaks with N>=3 (like suggested by ArchR-team)
# 	 pairs=pairs[pairs$score_1>=pairs$score_2,]
pairs <- pairs[order(-pairs$score_1), ]
all_st <- NULL
for (i in 1:nrow(pairs)) {
    if (nrow(pairs) > 0) {
        all_st <- rbind(all_st, p1[rownames(p1) == pairs[1, 5], ])
        pairs <- pairs[pairs$row_1 != pairs[1, 10], ]
        pairs <- pairs[-1, ]
    }
}

all_st <- as.data.frame(all_st)
all_st <- all_st[!duplicated(all_st), ]

recentered_final <- rbind(recentered_non_olap, all_st)
write.table(recentered_final, paste(outputpath, "recentered_final.filtered", sample, ".tsv", sep = ""),
    sep = "\t", quote = FALSE, row.names = FALSE
)
recentered_p <- StringToGRanges(recentered_final$new_peak, sep = c("-", "-"))
matrix.counts <- FeatureMatrix(
    fragments = Fragments(pbmc@assays$ATAC),
    features = recentered_p,
    sep = c("-", "-"),
    cells = colnames(pbmc)
)

pbmc[["X500peaksMACS2"]] <- CreateChromatinAssay(
    counts = matrix.counts,
    fragments = Fragments(pbmc@assays$ATAC)
)
DefaultAssay(pbmc) <- "X500peaksMACS2"

peak.data <- GetAssayData(object = pbmc, assay = "X500peaksMACS2", slot = "counts")

### use atac_fragments instead of passed_filters
total_fragments_cell <- pbmc$atac_fragments
peak.counts <- colSums(x = peak.data)
frip <- peak.counts / total_fragments_cell
pbmc <- AddMetaData(object = pbmc, metadata = frip, col.name = "frip_500MACS2")
pbmc <- AddMetaData(object = pbmc, metadata = peak.counts, col.name = "peak_RF_500MACS2")

annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
genome(annotation) <- "NA"
seqlevelsStyle(annotation) <- "UCSC"
genome(annotation) <- "mm10"
Annotation(pbmc) <- annotation


###########################################################
############ Quality Control OF SC-ATAC DATA################
###########################################################
# https://satijalab.org/signac/articles/pbmc_vignette.html

# compute nucleosome signal score per cell
pbmc <- NucleosomeSignal(object = pbmc)

# compute TSS enrichment score per cell
pbmc <- TSSEnrichment(object = pbmc, fast = FALSE)

pbmc$pct_reads_in_peaks <- pbmc$peak_RF_500MACS2 / pbmc$atac_fragments * 100

# add blacklist ratio and fraction of reads in peaks
blacklist_counts <- CountsInRegion(pbmc, assay = "X500peaksMACS2", region = blacklist_mm10)
blacklist_counts <- blacklist_counts[names(peak.counts)]

pbmc <- AddMetaData(object = pbmc, metadata = blacklist_counts, col.name = "blacklist_region_fragments")
pbmc$blacklist_ratio <- pbmc$blacklist_region_fragments / pbmc$peak_RF_500MACS2

# inspecting TSS-enrichment scores
pbmc$high.tss <- ifelse(pbmc$TSS.enrichment > 2, "High", "Low")
tss_plot <- TSSPlot(pbmc, group.by = "high.tss") + NoLegend()

# inspecting fragment length periodicity
pbmc$nucleosome_group <- ifelse(pbmc$nucleosome_signal > opt$ns_max, "NS > opt$ns_max", "NS < opt$ns_max")
if (length(table(pbmc$nucleosome_group)) > 1) {
    fragment_period_plot <- FragmentHistogram(object = pbmc, group.by = "nucleosome_group")
} else {
    fragment_period_plot <- NULL
}

QC_plot <- VlnPlot(
    object = pbmc,
    features = c(
        "pct_reads_in_peaks", "peak_RF_500MACS2",
        "TSS.enrichment", "blacklist_ratio", "nucleosome_signal"
    ),
    pt.size = 0.1,
    ncol = 5
)

pdf(paste(outputpath, "/", sample, "_0_QC.pdf", sep = ""), height = 6, width = 12)
print(tss_plot)
print(fragment_period_plot)
print(QC_plot)
dev.off()

# Add percent mito:
pbmc$percent.mito <- PercentageFeatureSet(pbmc, pattern = "^mt-", assay = "RNA")

# remove cells that are outliers for these QC metrics
pbmc <- subset(
    x = pbmc,
    # ATAC-filters:
    subset = peak_RF_500MACS2 > opt$prf_min &
        peak_RF_500MACS2 < opt$prf_max &
        pct_reads_in_peaks > opt$pct_min &
        blacklist_ratio < opt$bl_ratio &
        nucleosome_signal < opt$ns_max &
        TSS.enrichment > opt$tss &
        # RNA-filters:
        nFeature_RNA > opt$nfeature_min &
        nFeature_RNA < opt$nfeature_max &
        nCount_RNA > opt$ncount_min &
        nCount_RNA < opt$ncount_max &
        percent.mito < opt$mito_max * 100
)

##################################################
######## Gene expression data processing###########
##################################################
# Use params as in the Dan's script:
DefaultAssay(pbmc) <- "RNA"
pbmc <- SCTransform(pbmc, vars.to.regress = c("nCount_RNA", "percent.mito"), return.only.var.genes = F)
pbmc <- RunPCA(pbmc, npcs = opt$rna_pc_num, verbose = FALSE)
pbmc <- RunUMAP(
    pbmc,
    reduction = "pca",
    dims = 1:opt$rna_pc_num,
    reduction.name = "rna.umap",
    reduction.key = "rnaUMAP_"
)


##################################################
######## DNA accessibility data processing#########
##################################################

DefaultAssay(pbmc) <- "X500peaksMACS2"
### Will take all features,
pbmc <- FindTopFeatures(pbmc, min.cutoff = "q0")
pbmc <- RunTFIDF(pbmc)
pbmc <- RunSVD(pbmc)

pbmc <- RunUMAP(
    pbmc,
    reduction = "lsi",
    dims = opt$atac_pc_first:opt$atac_pc_num,
    reduction.name = "atac.umap",
    reduction.key = "atacUMAP_"
)

# Check if first LSI-component correlated with the sequencibg depth. If it is, then re-run using LSI components starting from 2 (for example, 2:30 instead of 1:30)
depth_corr_plot <- DepthCor(pbmc)
pdf(paste(outputpath, "/", sample, "_DepthCorrelation_1_QC.pdf", sep = ""), height = 6, width = 12)
print(depth_corr_plot)
dev.off()

# ATAC gene activity
DefaultAssay(pbmc) <- "X500peaksMACS2"
message("[ATAC] infer gene activity to ATACGeneActivity assay...")
gene.activities <- GeneActivity(pbmc)
# add the gene activity matrix to the Seurat object as a new assay and normalize it
pbmc[["ATACGeneActivity"]] <- CreateAssayObject(counts = gene.activities)
pbmc <- NormalizeData(
    object = pbmc,
    assay = "ATACGeneActivity",
    normalization.method = "LogNormalize",
    scale.factor = median(pbmc$nCount_ATACGeneActivity)
)


# Joint clustering
message("[Joint] cluster cells based on joint RNA and ATAC...")
DefaultAssay(pbmc) <- "X500peaksMACS2"
# build a joint neighbor graph using both assays
pbmc <- FindMultiModalNeighbors(
    object = pbmc,
    reduction.list = list("pca", "lsi"),
    dims.list = list(1:opt$rna_pc_num, opt$atac_pc_first:opt$atac_pc_num),
    verbose = TRUE
)
# build a joint UMAP visualization
pbmc <- RunUMAP(
    object = pbmc,
    nn.name = "weighted.nn", # weighted nearest neighbor
    reduction.name = "wnn.umap",
    reduction.key = "wnnUMAP_",
    verbose = TRUE
)
pbmc <- FindClusters(
    pbmc,
    graph.name = "wsnn",
    algorithm = 3, # SLM
    #    	resolution = opt$joint_cluster_res,
    verbose = TRUE
)

# plot joint clustering result
message("[Joint] ... plot joint clustering result")
p1 <- DimPlot(pbmc, label = TRUE, repel = TRUE, reduction = "wnn.umap") +
    NoLegend() +
    ggtitle("Joint WNN")
p2 <- DimPlot(pbmc,
    reduction = "rna.umap", label = TRUE,
    repel = TRUE, label.size = 2.5
) + NoLegend() + ggtitle("RNA")
p3 <- DimPlot(pbmc,
    reduction = "atac.umap", label = TRUE,
    repel = TRUE, label.size = 2.5
) + NoLegend() + ggtitle("ATAC")
p <- p1 + p2 + p3
pdf(paste(outputpath, "/", sample, "_2_Dimplots.pdf", sep = ""), height = 4, width = 12)
print(p)
dev.off()

DefaultAssay(pbmc) <- "RNA"

# Save object
saveRDS(pbmc, file = paste(outputpath, "/", sample, "_processed_multiome.rds", sep = ""))
