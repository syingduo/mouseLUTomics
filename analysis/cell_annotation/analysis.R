#!/bin/Rscript

## set up diretory and load libraries
setwd("/diskmnt/Projects/Users/s.yingduo/05gudmap/02multiome/LUT/analysis/all_18_samples_new")
library(Seurat)
library(SingleCellExperiment)
library(ComplexHeatmap)
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(viridis)
library(reshape2)

## read in object
joint <- readRDS("/diskmnt/Projects/Users/s.yingduo/05gudmap/02multiome/LUT/all_18_samples_new/results_res.0.3/snRNA_ATAC_jointly_analyzed_for_all_18_samples_new.rds")
meta <- joint@meta.data
dim(meta)

## add clincal information and doublets dectection to meta file
tmp <- word(rownames(meta), sep = fixed("_"), 1)

samples <- names(table(tmp))
age <- c("P21", "P21", "E16", "E16", "P21", "P21", "P1", "P1", "E16", "E16", "P21", "P21", "P1", "P1", "P1", "P1", "E16", "E16")
sex <- c("F", "M", "M", "F", "F", "M", "M", "F", "F", "M", "M", "F", "F", "M", "M", "F", "M", "F")
sample_id <- c("GM005MwP21_F", "GM006MwP21_M", "GM007MwE16_M", "GM008MwE16_F", "GM011MwP21_F", "GM012MwP21_M", "GM013MwP1_M", "GM014MwP1_F", "GM015MwE16_F", "GM016MwE16_M", "GM017MwP21_M", "GM018MwP21_F", "GM019MwP1_F", "GM020MwP1_M", "GM001MwP1_M", "GM002MwP1_F", "GM003MwE16_M", "GM004MwE16_F")
dt <- data.frame(samples, age, sex, sample_id)

tmp1 <- tmp
for (j in dt$samples) {
    tmp1 <- gsub(j, dt$age[which(dt$samples == j)], tmp1)
}
age <- tmp1
tmp1 <- tmp
for (j in dt$samples) {
    tmp1 <- gsub(j, dt$sex[which(dt$samples == j)], tmp1)
}
sex <- tmp1
tmp1 <- tmp
for (j in dt$samples) {
    tmp1 <- gsub(j, dt$sample_id[which(dt$samples == j)], tmp1)
}
sample_id <- tmp1
meta <- cbind(meta, sample_id, age, sex)
joint@meta.data <- meta

doublets <- read.table("/diskmnt/Projects/Users/s.yingduo/05gudmap/02multiome/LUT/doublets_detection/doublets_all_18_samples.tsv", header = T, row.names = 1)
joint <- AddMetaData(joint, doublets)
meta <- joint@meta.data
dim(meta)
saveRDS(joint, "/diskmnt/Projects/Users/s.yingduo/05gudmap/02multiome/LUT/all_18_samples_new/results_res.0.3/snRNA_ATAC_jointly_analyzed_for_all_18_samples_new_add.rds")


### snRNA-seq analysis
### i. Find markers
joint <- readRDS("/diskmnt/Projects/Users/s.yingduo/05gudmap/02multiome/LUT/all_18_samples_new/results_res.0.3/snRNA_ATAC_jointly_analyzed_for_all_18_samples_new_add.rds")
dim(joint)

DefaultAssay(joint) <- "RNA"
joint.rna.markers <- FindAllMarkers(joint,
    assay = "RNA", only.pos = TRUE, min.pct = 0.25,
    logfc.threshold = 0.25
)

rna.markers.table <- joint.rna.markers %>%
    group_by(cluster) %>%
    top_n(n = 60, wt = avg_log2FC)

save(joint.rna.markers, rna.markers.table, file = "rna.markers.RData")
load("rna.markers.RData")

cluster <- 0:37
topGene <- character()
cell_type <- rep("NULL", length(cluster))
for (i in 1:38) {
    topGene[i] <- paste(rna.markers.table$gene[(i * 60 - 59):(i * 60)], collapse = ",")
}
markers.check <- data.frame(cluster, topGene, cell_type)
tmp <- paste("cluster", cluster, sep = "_")
marker_genes <- paste(tmp, topGene, sep = ":")
write.table(marker_genes, "marker_genes1.tsv", sep = "\t", quote = F, col.names = F, row.names = F)


### ii. manual annotation
markers.manual <- read.table("cell_annotation/cell_annotation.tsv", header = T, sep = "\t")
markers.manual <- markers.manual[order(markers.manual$Cluster), ]
new.cluster.id <- markers.manual$Cell_type
names(new.cluster.id) <- levels(joint)
joint <- RenameIdents(joint, new.cluster.id)

### plot annotated umamp
pdf("rna.umap.cell.annotation.pdf", width = 10, height = 10)
DefaultAssay(joint) <- "RNA"
p1 <- DimPlot(joint, label = TRUE, repel = T, raster = FALSE, label.size = 4, reduction = "rna.umap") + NoLegend()
p1
dev.off()

### feature plot
marker_genes <- as.data.frame(t(read.csv("cell_annotation/cellxgene.txt", header = F)))
v1 <- marker_genes$V1
v2 <- marker_genes$V4
v2[v2 != ""]
intersect(v1, v2)
genes <- c("Gsta4", "Wfdc2", "Ly6d")

pdf("cell_annotation/tmp.pdf", height = 5, width = 15)
DefaultAssay(joint) <- "RNA"
FeaturePlot(joint, genes, pt.size = 0.2, ncol = 3, raster = FALSE) + theme(
    # legend.position = "none",
    axis.title.x = element_blank(), axis.title.y = element_blank()
)
dev.off()


### iii. Dot plot
markers <- unique(unlist(strsplit(markers.manual$Makers, ",")))
pdf("dot_plot.pdf", height = 10, width = 30)
DefaultAssay(joint) <- "RNA"
p <- DotPlot(object = joint, features = markers, cols = c("lightblue", "red"), dot.scale = 2.5) +
    theme(
        axis.text.x = element_text(size = 12, angle = 70, hjust = 1, vjust = 1),
        axis.text.y = element_text(size = 14)
    )
p
dev.off()


### iv. make the stack plot
dt <- meta %>%
    count(dataset, seurat_clusters) %>%
    group_by(seurat_clusters)
percent <- c(dt[1:25, ]$n / 5224 * 100, dt[26:55, ]$n / 7054 * 100, dt[56:76, ]$n / 5754 * 100, dt[77:103, ]$n / 7211 * 100, dt[104:133, ]$n / 6580 * 100, dt[134:157, ]$n / 7647 * 100)
dt <- data.frame(dt, percent)

theme_bar <- function(..., bg = "white") {
    require(grid)
    theme_classic(...) +
        theme(
            rect = element_rect(fill = bg),
            plot.margin = unit(rep(0.5, 4), "lines"),
            panel.background = element_rect(fill = "transparent", color = "black"),
            panel.border = element_rect(fill = "transparent", color = "transparent"),
            panel.grid = element_blank(), # 去网格线
            axis.title.x = element_blank(), # 去x轴标签,
            axis.title.y = element_blank(),
            # axis.title.y = element_text(face = "bold", size = 8), # y轴标签加粗及字体大小
            axis.text.y = element_text(face = "bold", size = 8), # 坐标轴刻度标签加粗
            axis.text.x = element_text(size = 8, face = "bold", angle = 60, hjust = 1, vjust = 1),
            # axis.ticks = element_line(color='black'),#坐标轴刻度线
            # axis.ticks.margin = unit(0.8,"lines"),
            legend.title = element_blank(), # 去除图例标题
            # legend.justification=c(1,0),#图例在画布的位置(绘图区域外)
            # legend.position = c(0.76, 0.86), # 图例在绘图区域的位置
            # legend.position='top',#图例放在顶部
            # legend.direction = "horizontal", # 设置图例水平放置
            # legend.spacing.x = unit(2, 'cm'),
            legend.text = element_text(face = "bold", size = 10, margin = margin(r = 20)),
            legend.background = element_rect(linetype = "solid", colour = "black"),
            # legend.margin=margin(0,0,-7,0)#图例与绘图区域边缘的距离
            # legend.box.margin =margin(-10,0,0,0)
            plot.title = element_text(size = 18, hjust = 0.5, face = "bold")
        )
}

pdf("tmp2.pdf", height = 8, width = 8)
ggplot(dt, aes(x = dataset, y = percent, fill = seurat_clusters)) +
    geom_bar(stat = "identity", position = "stack") +
    scale_y_continuous(expand = c(0, 0)) +
    theme_bar()
dev.off()


## integration
library(Seurat)
library(cowplot)
library(harmony)

pdf(paste0(out_path, "/RunHarmony.pdf"), width = 12)
options(repr.plot.height = 2.5, repr.plot.width = 6)
combined <- RunHarmony(object = combined, group.by.vars = "batch", reduction = "pca", assay.use = "SCT")

combined <- RunHarmony(object = combined, group.by.vars = "batch", reduction = "pca", assay.use = "RNA")
dev.off()

harmony_embeddings <- Embeddings(combined, "harmony")
harmony_embeddings[1:5, 1:5]

pdf(paste0(out_path, "/corrected_PC.pdf"), width = 12)
p1 <- DimPlot(object = combined, reduction = "harmony", pt.size = .1, group.by = "batch")
p2 <- VlnPlot(object = combined, features = "harmony_1", group.by = "batch", pt.size = .1)
plot_grid(p1, p2, align = "h", nrow = 1, rel_widths = c(1 / 2, 1 / 2))
dev.off()
