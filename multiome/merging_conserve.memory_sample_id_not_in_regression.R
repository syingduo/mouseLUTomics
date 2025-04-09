library(optparse)
library(dplyr)
library(Seurat)
library(tibble)
library(cowplot)
library(ggplot2)

library(data.table)

option_list <- c(
    make_option(c("--merging_name"),
        type = "character",
        default = "integration",
        help = "specify the name of the integration",
        metavar = "character"
    ),
    make_option(c("--sample_path_table"),
        type = "character",
        default = NULL,
        help = "path to a table which stores rds path for all the samples needed for integration, additional columns are allowed",
        metavar = "character"
    ),
    make_option(c("--npcs"),
        type = "integer",
        default = 30,
        help = "top # of PCs for downstream analysis",
        metavar = "integer"
    ),
    make_option(c("--n_neighbors"),
        type = "integer",
        default = 30,
        help = "# of neighbors used for UMAP and clustering",
        metavar = "integer"
    ),
    make_option(c("--cluster_resolution"),
        type = "double",
        default = 0.3,
        help = "resolution for clustering",
        metavar = "double"
    ),
    make_option(c("-o", "--output"),
        type = "character",
        default = "./",
        help = "output path",
        metavar = "character"
    ),
    make_option(c("--force.reprocess"),
        type = "logical",
        default = FALSE,
        help = "whether to re-process the data when rds object already exists",
        metavar = "character"
    ),
    make_option(c("--marker.toplot"),
        type = "character",
        default = NULL,
        help = "path to the file where markers are organized",
        metavar = "character"
    )
)
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)
merging_name <- opt$merging_name
sample_path_table <- opt$sample_path_table
out_path <- opt$output

if (!file.exists(out_path)) {
    dir.create(out_path)
}

# sink(paste0(out_path, "integration_processing.txt"))

sample_info_summary <- read.table(opt$sample_path_table, head = TRUE, sep = "\t")
n_sample <- nrow(sample_info_summary)
if ((!file.exists(paste0(out_path, "merged_", merging_name, ".rds"))) | opt$force.reprocess == TRUE) {
    obj_list <- list()
    sample_ids <- c()
    for (i in 1:n_sample) {
        sample_id <- sample_info_summary[i, "sample_id"] %>% as.character()

        cat("read in sample ", sample_id, "\n")
        tmp <- readRDS(sample_info_summary[i, "Path"] %>% as.character())
        tmp$sample_id <- sample_id
        tmp$original_barcode <- rownames(tmp@meta.data)
        DefaultAssay(tmp) <- "RNA"
        for (assay in setdiff(Assays(tmp), c("RNA", "SCT"))) {
            tmp[[assay]] <- NULL
        }
        obj_list[[sample_id]] <- tmp
        sample_ids <- c(sample_ids, sample_id)
    }

    cat("merging...\n")
    merged <- merge(obj_list[[sample_ids[1]]],
        y = unlist(obj_list)[-1],
        add.cell.ids = sample_ids,
        project = merging_name
    )

    ### filter mt genes
    merged <- subset(merged, percent.mito <= 1)

    # cat("SCTransform...\n")
    # merged <- SCTransform(merged, vars.to.regress = c("nCount_RNA", "percent.mito"), return.only.var.genes = F, conserve.memory = TRUE)

    merged <- NormalizeData(merged, normalization.method = "LogNormalize", scale.factor = 10000)
    merged <- FindVariableFeatures(merged, selection.method = "vst")
    merged <- ScaleData(obj.ref)

    cat("PCA...\n")
    merged <- RunPCA(merged, npcs = opt$npcs)

    cat("Run UMAP...\n")
    merged <- RunUMAP(merged, reduction = "pca", dims = 1:opt$n_neighbors)

    cat("Find Neighbors\n")
    merged <- FindNeighbors(merged, reduction = "pca", dims = 1:opt$n_neighbors)

    cat("Find clusters\n")
    merged <- FindClusters(merged, resolution = opt$cluster_resolution, algorithm = 2) #
    saveRDS(merged, paste0(out_path, "merged_", merging_name, ".rds"))
} else {
    merged <- readRDS(paste0(out_path, "merged_", merging_name, ".rds"))
}



sink()


if ((!is.null(opt$marker.toplot)) & file.exists(opt$marker.toplot)) {
    top.markers <- fread(opt$marker.toplot, head = TRUE, sep = "\t")$gene_symbol %>% as.character()
    pdf(paste0(out_path, "Marker_Expression_in_", merging_name, ".pdf"), useDingbats = FALSE, width = 20, height = 50)
    FeaturePlot(merged, features = top.markers, cols = c("grey", "red"), pt.size = 0.1) %>% print()
    dev.off()
}

# Visualization
pdf(paste0(out_path, "DimReductionInCase_", merging_name, ".pdf"), width = 24, useDingbats = FALSE)
p1 <- DimPlot(object = merged, reduction = "umap", group.by = "sample_id", pt.size = 0.1, raster = FALSE) + theme(aspect.ratio = 1) + ggtitle(merging_name)
p2 <- DimPlot(object = merged, reduction = "umap", label = TRUE, raster = FALSE) + theme(aspect.ratio = 1)
print(plot_grid(p1, p2))
dev.off()


pdf(paste0(out_path, "DimReductionInCase_", merging_name, "_split_by_cases.pdf"), width = 7.2 * n_sample, useDingbats = FALSE)
DimPlot(object = merged, reduction = "umap", split.by = "sample_id", label = TRUE, pt.size = 0.1, raster = FALSE)
# DimPlot(object = merged, reduction = "umap", split.by = "sample_id",
#          group.by = "cell_type",pt.size=0.05,label=TRUE)
dev.off()
