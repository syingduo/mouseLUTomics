#!/bin/Rscript

library(Seurat)
library(Signac)
library(dplyr)
library(tibble)
library(ggplot2)
library(patchwork)

merging_name <- "all_18_samples_new"
out_path <- paste0("/diskmnt/Projects/Users/s.yingduo/05gudmap/02multiome/LUT/", merging_name, "/results_res.0.15/")
dir.create(out_path, recursive = T)

obj_ATAC <- readRDS(paste0("/diskmnt/Projects/Users/s.yingduo/05gudmap/02multiome/LUT/", merging_name, "/peaks/", merging_name, "_snATAC_Merged_BasedOnSelectedPeaks_Normalized.rds.gz"))
obj_RNA <- readRDS(paste0("/diskmnt/Projects/Users/s.yingduo/05gudmap/02multiome/LUT/", merging_name, "/merged_", merging_name, ".rds"))


obj_joint <- obj_ATAC
obj_joint[["RNA"]] <- CreateAssayObject(counts = GetAssayData(obj_RNA, assay = "RNA"))
obj_joint[["SCT"]] <- CreateAssayObject(data = GetAssayData(obj_RNA, assay = "SCT"))
obj_joint[["umap"]] <- NULL

DefaultAssay(obj_joint) <- "SCT"
obj_joint[["pca"]] <- obj_RNA[["pca"]]
obj_joint[["rna.umap"]] <- obj_RNA[["umap"]]

DefaultAssay(obj_joint) <- "peaksMACS2"
obj_joint[["atac.umap"]] <- obj_ATAC[["umap"]]

obj_joint@meta.data$SCT_snn_res.0.5 <- obj_RNA@meta.data[obj_joint@meta.data %>% rownames(), "SCT_snn_res.0.5"]


#### Now we ready to perform joint Dim-reduction:
# Joint clustering
message("[Joint] cluster cells based on joint RNA and ATAC...")
DefaultAssay(obj_joint) <- "peaksMACS2"

# build a joint neighbor graph using both assays
obj_joint <- FindMultiModalNeighbors(
    object = obj_joint, reduction.list = list("pca", "lsi"),
    dims.list = list(1:30, 2:30), verbose = TRUE
)

# build a joint UMAP visualization
obj_joint <- RunUMAP(
    object = obj_joint, nn.name = "weighted.nn", # weighted nearest neighbor
    reduction.name = "wnn.umap", reduction.key = "wnnUMAP_", verbose = TRUE
)

#### Run this later, once you are happy with the results:
obj_joint <- FindClusters(obj_joint,
    graph.name = "wsnn", algorithm = 3, # SLM
    resolution = 0.15,
    verbose = TRUE
)


### Now make plots:
p1 <- DimPlot(obj_joint, reduction = "wnn.umap", label = TRUE, repel = TRUE, label.size = 2.5, raster = FALSE) + NoLegend() + ggtitle("Joint WNN") + theme(aspect.ratio = 1)
p2 <- DimPlot(obj_joint, reduction = "rna.umap", label = TRUE, repel = TRUE, label.size = 2.5, raster = FALSE) + NoLegend() + ggtitle("RNA") + theme(aspect.ratio = 1)
p3 <- DimPlot(obj_joint, reduction = "atac.umap", label = TRUE, repel = TRUE, label.size = 2.5, raster = FALSE) + NoLegend() + ggtitle("ATAC") + theme(aspect.ratio = 1)
p_wsnn <- (p1 + p2 + p3) + plot_annotation(title = "wsnn based cluster")

p1 <- DimPlot(obj_joint, group.by = "peaksMACS2_snn_res.0.8", label = TRUE, repel = TRUE, reduction = "wnn.umap", raster = FALSE) + NoLegend() + ggtitle("Joint WNN") + theme(aspect.ratio = 1)
p2 <- DimPlot(obj_joint, group.by = "peaksMACS2_snn_res.0.8", reduction = "rna.umap", label = TRUE, repel = TRUE, label.size = 2.5, raster = FALSE) + NoLegend() + ggtitle("RNA") + theme(aspect.ratio = 1)
p3 <- DimPlot(obj_joint, group.by = "peaksMACS2_snn_res.0.8", reduction = "atac.umap", label = TRUE, repel = TRUE, label.size = 2.5, raster = FALSE) + NoLegend() + ggtitle("ATAC") + theme(aspect.ratio = 1)
p_peaksMACS2 <- (p1 + p2 + p3) + plot_annotation(title = "peaksMACS2 based cluster")

# SCT_snn_res.0.5
p1 <- DimPlot(obj_joint, group.by = "SCT_snn_res.0.5", label = TRUE, repel = TRUE, reduction = "wnn.umap", raster = FALSE) + NoLegend() + ggtitle("Joint WNN") + theme(aspect.ratio = 1)
p2 <- DimPlot(obj_joint, group.by = "SCT_snn_res.0.5", reduction = "rna.umap", label = TRUE, repel = TRUE, label.size = 2.5, raster = FALSE) + NoLegend() + ggtitle("RNA") + theme(aspect.ratio = 1)
p3 <- DimPlot(obj_joint, group.by = "SCT_snn_res.0.5", reduction = "atac.umap", label = TRUE, repel = TRUE, label.size = 2.5, raster = FALSE) + NoLegend() + ggtitle("ATAC") + theme(aspect.ratio = 1)
p_SCT <- (p1 + p2 + p3) + plot_annotation(title = "SCT based cluster")

pdf(paste0(out_path, "UMAP_for_", merging_name, ".pdf", sep = ""), height = 6, width = 20)
print(p_wsnn)
print(p_peaksMACS2)
# print(p_SCT)
dev.off()

saveRDS(obj_joint, paste0(out_path, "snRNA_ATAC_jointly_analyzed_for_", merging_name, ".rds"))















### Check for doublets
cluster_doublet_proportion <- obj_joint@meta.data %>%
    group_by(wsnn_res.0.8, predicted_doublet) %>%
    count() %>%
    tidyr::pivot_wider(names_from = "predicted_doublet", values_from = "n") %>%
    mutate(`TRUE` = ifelse(is.na(`TRUE`), 0, `TRUE`)) %>%
    mutate(doublet_proportion = `TRUE` / (`TRUE` + `FALSE`)) %>%
    arrange(desc(doublet_proportion), na.rm = TRUE)
tmp <- kmeans(cluster_doublet_proportion$doublet_proportion, 2)
doublet_enriched_clusters <- cluster_doublet_proportion[which(tmp$cluster == 1), "wsnn_res.0.8", drop = TRUE] %>%
    as.character() %>%
    as.numeric()

BC_doublet_1 <- obj_joint@meta.data %>%
    rownames_to_column() %>%
    filter(wsnn_res.0.8 %in% doublet_enriched_clusters) %>%
    .$rowname
BC_doublet_2 <- obj_joint@meta.data %>%
    rownames_to_column() %>%
    filter(is.na(predicted_doublet) | `predicted_doublet` == TRUE) %>%
    .$rowname

BC_doublet_all <- c(BC_doublet_1, BC_doublet_2) %>% unique()
saveRDS(BC_doublet_all, paste0(out_path, "doublet_to_remove.rds"))


p1 <- DimPlot(obj_joint, reduction = "wnn.umap", cells.highlight = BC_doublet_all) + ggtitle("Joint WNN") + theme(aspect.ratio = 1)
p2 <- DimPlot(obj_joint, reduction = "rna.umap", cells.highlight = BC_doublet_all) + ggtitle("RNA") + theme(aspect.ratio = 1)
p3 <- DimPlot(obj_joint, reduction = "atac.umap", cells.highlight = BC_doublet_all) + ggtitle("ATAC") + theme(aspect.ratio = 1)
p <- (p1 + p2 + p3) + plot_annotation(title = "highlighted doublets")
pdf(paste0(out_path, "UMAP_for_", merging_name, "_highlighting_doublets.pdf", sep = ""), height = 6, width = 20)
print(p)
dev.off()
