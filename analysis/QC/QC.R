#!/bin/Rscript

## set up directory
setwd("/diskmnt/Projects/Users/s.yingduo/05gudmap/02multiome/LUT/qc")
library(ggplot2)
library(cowplot)
library(stringr)

## make table
path <- "/diskmnt/Projects/GUDMAP_analysis/combo_snRNA_snATAC"
sample <- read.table("/diskmnt/Projects/Users/s.yingduo/05gudmap/02multiome/LUT/all_18_samples/sample_rds_path.txt", sep = "\t", header = T)
sample <- sample$sample_id

meta <- list()
for (i in sample) {
    file <- paste0(path, "/", i, "/outs/summary.csv")
    tmp <- read.csv(file, header = T)
    tmp <- tmp[, c("Sample.ID", "Estimated.number.of.cells", "ATAC.Mean.raw.read.pairs.per.cell", "ATAC.Median.high.quality.fragments.per.cell", "GEX.Mean.raw.reads.per.cell", "GEX.Median.genes.per.cell")]
    meta <- rbind(meta, tmp)
}
dim(meta)

## make plots
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
            axis.text.y = element_text(size = 12 face = "bold"), # 坐标轴刻度标签加粗
            axis.ticks.x = element_blank(), # 坐标轴刻度线
            # axis.ticks.margin = unit(0.8,"lines"),
            axis.text.x = element_text(size = 10, face = "bold", angle = 65, hjust = 1, vjust = 1),
            legend.title = element_blank(), # 去除图例标题
            # legend.justification=c(1,0),#图例在画布的位置(绘图区域外)
            legend.position = c(0.39, 0.76), # 图例在绘图区域的位置
            # legend.position='top',#图例放在顶部
            legend.direction = "horizontal", # 设置图例水平放置
            # legend.spacing.x = unit(2, 'cm'),
            legend.text = element_text(face = "bold", size = 11, margin = margin(r = 20)),
            legend.background = element_rect(linetype = "solid", colour = "black"),
            # legend.margin=margin(0,0,-7,0)#图例与绘图区域边缘的距离
            # legend.box.margin =margin(-10,0,0,0)
            plot.title = element_text(size = 14, face = "bold", hjust = 0.5, vjust = -6.5)
        )
}

p1 <- ggplot(meta, mapping = aes(x = reorder(Sample.ID, -Estimated.number.of.cells), y = Estimated.number.of.cells)) +
    geom_bar(
        stat = "identity",
        fill = "#00AFBB",
        width = 0.6
    ) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, 20 * 1000), breaks = seq(0, 20 * 1000, 2000)) +
    labs(title = "Number of cells") +
    annotate("text",
        x = 15, y = 18000, label = paste0("Mean: ", round(mean(meta$Estimated.number.of.cells), 2)),
        color = "purple", size = 5, fontface = "bold"
    ) +
    annotate("text",
        x = 15.1, y = 17400, label =
            paste0("SD: ", round(sd(meta$Estimated.number.of.cells), 2)),
        color = "purple", size = 5, fontface = "bold"
    ) +
    theme_bar() +
    theme(axis.text.x = element_blank())

p2 <- ggplot(meta, mapping = aes(x = reorder(Sample.ID, -Estimated.number.of.cells), y = ATAC.Mean.raw.read.pairs.per.cell)) +
    geom_bar(
        stat = "identity",
        fill = "#FFB300",
        width = 0.6
    ) +
    scale_y_continuous(
        expand = c(0, 0),
        limits = c(0, 1000 * 60),
        breaks = seq(0, 1000 * 60, 6000)
    ) +
    labs(title = "ATAC Mean reads per cell") +
    annotate("text",
        x = 15, y = 54000, label = paste0("Mean: ", round(mean(meta$ATAC.Mean.raw.read.pairs.per.cell), 2)),
        color = "purple", size = 5, fontface = "bold"
    ) +
    annotate("text",
        x = 15.2, y = 52500, label =
            paste0("SD: ", round(sd(meta$ATAC.Mean.raw.read.pairs.per.cell), 2)),
        color = "purple", size = 5, fontface = "bold"
    ) +
    theme_bar() +
    theme(axis.text.x = element_blank())

p3 <- ggplot(meta, mapping = aes(x = reorder(Sample.ID, -Estimated.number.of.cells), y = ATAC.Median.high.quality.fragments.per.cell)) +
    geom_bar(
        stat = "identity",
        fill = "#FA2017",
        width = 0.6
    ) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, 1000 * 25), breaks = seq(0, 1000 * 25, 2500)) +
    labs(title = "ATAC Median genes per cell") +
    annotate("text",
        x = 15, y = 23000, label = paste0("Mean: ", round(mean(meta$ATAC.Median.high.quality.fragments.per.cell), 2)),
        color = "purple", size = 5, fontface = "bold"
    ) +
    annotate("text",
        x = 15.2, y = 22200, label =
            paste0("SD: ", round(sd(meta$ATAC.Median.high.quality.fragments.per.cell), 2)),
        color = "purple", size = 5, fontface = "bold"
    ) +
    theme_bar() +
    theme(axis.text.x = element_blank())

p4 <- ggplot(meta, mapping = aes(x = reorder(Sample.ID, -Estimated.number.of.cells), y = GEX.Mean.raw.reads.per.cell)) +
    geom_bar(
        stat = "identity",
        fill = "hotpink",
        width = 0.6
    ) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, 1000 * 70), breaks = seq(0, 1000 * 70, 7000)) +
    labs(title = "RNA Mean reads per cell") +
    annotate("text",
        x = 15, y = 63000, label = paste0("Mean: ", round(mean(meta$GEX.Mean.raw.reads.per.cell), 2)),
        color = "purple", size = 5, fontface = "bold"
    ) +
    annotate("text",
        x = 15.2, y = 61000, label =
            paste0("SD: ", round(sd(meta$GEX.Mean.raw.reads.per.cell), 2)),
        color = "purple", size = 5, fontface = "bold"
    ) +
    theme_bar() +
    theme(axis.text.x = element_blank())


p5 <- ggplot(meta, mapping = aes(x = reorder(Sample.ID, -Estimated.number.of.cells), y = GEX.Median.genes.per.cell)) +
    geom_bar(
        stat = "identity",
        fill = "#00C5FF",
        width = 0.6
    ) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, 100 * 25), breaks = seq(0, 100 * 25, 250)) +
    labs(title = "RNA Median genes per cell") +
    annotate("text",
        x = 13, y = 2200, label = paste0("Mean: ", round(mean(meta$GEX.Median.genes.per.cell), 2)),
        color = "purple", size = 5, fontface = "bold"
    ) +
    annotate("text",
        x = 13, y = 2100, label =
            paste0("SD: ", round(sd(meta$GEX.Median.genes.per.cell), 2)),
        color = "purple", size = 5, fontface = "bold"
    ) +
    theme_bar()

pdf("QC of LUT.pdf", width = 12, height = 35)
plot_grid(p1, p2, p3, p4, p5, ncol = 1, nrow = 5, align = "v")
dev.off()

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


## analysis
joint <- readRDS("/diskmnt/Projects/Users/s.yingduo/05gudmap/02multiome/LUT/all_18_samples_new/results_res.0.3/snRNA_ATAC_jointly_analyzed_for_all_18_samples_new_add.rds")
dim(joint)

### meta information
meta <- joint@meta.data
dim(meta)

### UMAP plot
table(meta$seurat_clusters)
pdf("rna.umap.pdf")
DefaultAssay(joint) <- "RNA"
p1 <- DimPlot(joint, label = TRUE, raster = FALSE) #+ NoLegend()
p1
dev.off()

pdf("rna.umap_sample.pdf")
DefaultAssay(joint) <- "RNA"
p1 <- DimPlot(joint, label = F, raster = FALSE, group.by = "sample_id")
p1
dev.off()

pdf("rna.umap_sex.pdf")
DefaultAssay(joint) <- "RNA"
p1 <- DimPlot(joint, label = F, raster = FALSE, group.by = "sex")
p1
dev.off()

pdf("rna.umap_doublets.pdf")
DefaultAssay(joint) <- "RNA"
FeaturePlot(joint, features = "doublet_score_rna", raster = FALSE) & theme(plot.title = element_text(size = 10))
dev.off()

pdf("rna.umap_mito.pdf")
DefaultAssay(joint) <- "RNA"
FeaturePlot(joint, features = "percent.mito", raster = FALSE) & theme(plot.title = element_text(size = 10))
dev.off()

