# Clustering and annotation on the B2MKO and Mix subset

# Load required packages
library(Seurat)
library(ggplot2)
library(readxl)
library(tidyverse)

# -------------------------------------------------------------------
# Load unfiltered data and set plotting order for conditions
# -------------------------------------------------------------------

dat <- readRDS("./objects/b2mko_mix_subset_from_all_samples_QCed_object.rds")
condition_levels <- c("B2MKO", "Mix")
dat@meta.data <- dat@meta.data %>%
  mutate(condition = factor(condition, levels = condition_levels))

# -------------------------------------------------------------------
# Initial normalization, feature selection, scaling, and PCA
# -------------------------------------------------------------------

# Per-cell log normalization
dat <- NormalizeData(dat)

# Feature selection
dat <- FindVariableFeatures(dat, selection.method = "vst", nfeatures = 2000)
top10_hvg <- head(VariableFeatures(dat), 10)
plot1 <- VariableFeaturePlot(dat)
plot2 <- LabelPoints(plot = plot1, points = top10_hvg, repel = TRUE)
ggsave("./plots/B2MKO_mix_hvg_scatter.png", plot = plot2)

# Data scaling
all.genes <- rownames(dat)
dat <- ScaleData(dat, features = all.genes)

# PCA
dat <- RunPCA(dat, features = VariableFeatures(object = dat), npcs = 100)
ElbowPlot(dat, ndims = 100)
ggsave(filename = "./plots/B2MKO_mix_PCA_elbowplot.png")

# -------------------------------------------------------------------
# Neighborhood identification, clustering, and UMAP
# -------------------------------------------------------------------

set.seed(00101)
resolutions <- c(0.1, 0.2, 0.5, 0.8, 1.0, 1.5, 2.0)

# Neighborhood identification
dat <- FindNeighbors(dat, reduction = "pca", dims = 1:30)

# Clustering across multiple resolutions
for (res in resolutions) {
  cluster_name <- paste0("seurat_cluster_res_", res)
  dat <- FindClusters(dat, cluster.name = cluster_name, resolution = res)
}

# UMAP
dat <- RunUMAP(dat, dims = 1:30)

for (res in resolutions) {
  cluster_name <- paste0("seurat_cluster_res_", res)
  Idents(dat) <- cluster_name
  DimPlot(dat)
  plot_name <- paste0("./plots/B2MKO_mix_UMAP_", cluster_name, ".png")
  ggsave(filename = plot_name)
}

# Find cluster markers for resolution 0.5
Idents(dat) <- "seurat_cluster_res_0.5"
seurat05_markers <- FindAllMarkers(dat, only.pos = TRUE)
write.csv(seurat05_markers, file = "./numerical_results/B2MKO_Mix_seurat05_cluster_markers.csv")

saveRDS(dat, file = "./objects/B2MKO_Mix_pca_umap_clustered.rds")

# -------------------------------------------------------------------
# Remove low-quality cluster and re-run clustering workflow
# -------------------------------------------------------------------

# NOTE: Seurat resolution 0.5, cluster 12 showed high mt-gene content and no distinct markers. Remove the cells in the cluster

dat <- readRDS("./objects/B2MKO_Mix_pca_umap_clustered.rds")

Idents(dat) <- "seurat_cluster_res_0.5"
dat_qc <- subset(dat, subset = seurat_cluster_res_0.5 != 12)

# Drop previous cluster columns from metadata
cluster_columns_to_drop <- c(
  "seurat_cluster_res_0.1",
  "seurat_cluster_res_0.2",
  "seurat_cluster_res_0.5",
  "seurat_cluster_res_0.8",
  "seurat_cluster_res_1",
  "seurat_cluster_res_1.5",
  "seurat_cluster_res_2"
)

dat_qc@meta.data <- dat_qc@meta.data %>% select(-all_of(cluster_columns_to_drop))

# Feature selection for filtered object
dat_qc <- FindVariableFeatures(dat_qc, selection.method = "vst", nfeatures = 2000)
top10_hvg <- head(VariableFeatures(dat_qc), 10)
plot1 <- VariableFeaturePlot(dat_qc)
plot2 <- LabelPoints(plot = plot1, points = top10_hvg, repel = TRUE)
ggsave("./plots/B2MKO_mix_QC2_hvg_scatter.png", plot = plot2)

# Data scaling
all.genes <- rownames(dat_qc)
dat_qc <- ScaleData(dat_qc, features = all.genes)

# PCA
dat_qc <- RunPCA(dat_qc, features = VariableFeatures(object = dat_qc), npcs = 100)
ElbowPlot(dat_qc, ndims = 100)
ggsave(filename = "./plots/B2MKO_mix_QC2_PCA_elbowplot.png")

set.seed(00101)

# Neighborhood identification
dat_qc <- FindNeighbors(dat_qc, reduction = "pca", dims = 1:30)

# Clustering across multiple resolutions
for (res in resolutions) {
  cluster_name <- paste0("seurat_cluster_res_", res)
  dat_qc <- FindClusters(dat_qc, cluster.name = cluster_name, resolution = res)
}

# UMAP
dat_qc <- RunUMAP(dat_qc, dims = 1:30)

for (res in resolutions) {
  cluster_name <- paste0("seurat_cluster_res_", res)
  Idents(dat_qc) <- cluster_name
  DimPlot(dat_qc)
  plot_name <- paste0("./plots/B2MKO_mix_QC2_UMAP_", cluster_name, ".png")
  ggsave(filename = plot_name)
}

saveRDS(dat_qc, file = "./objects/B2MKO_Mix_QC2_REpca_umap_clustered.rds")

# -------------------------------------------------------------------
# Cell-type annotation using cluster resolution 0.5
# -------------------------------------------------------------------

dat_qc <- readRDS("./objects/B2MKO_Mix_QC2_REpca_umap_clustered.rds")

# Marker set used for level-1 annotation
# DC: Flt3, Zbtb46
# pDC: Siglech, Bcl11a
# NK: Ncr1, Klrb1c
# T cell: Cd3g, Trac
# CD8: Cd8a
# CD8 activated: Pdcd1, Havcr2, Lag3, Cd69
# CD4: Cd4
# Cd4 Treg: Foxp3, Tnfrsf4 
# gd T: Trgc1, Il23r 
# Macrophage: Csf1r
# Neutrophil: S100a8, S100a9
# B cell: Mzb1, Pou2af1
# Mast cell: Cpa3, Tpsab1
# Proliferation: Mki67, Top2a, Cdk1

marker_genes <- c(
  "Flt3", "Zbtb46", "Siglech", "Bcl11a", "Ncr1", "Klrb1c", "Cd3g", "Trac",
  "Cd8a", "Pdcd1", "Havcr2", "Lag3", "Cd69", "Cd4", "Foxp3", "Tnfrsf4",
  "Trgc1", "Il23r", "Csf1r", "S100a8", "S100a9", "Mzb1", "Pou2af1",
  "Cpa3", "Tpsab1", "Mki67", "Top2a", "Cdk1"
)

FeaturePlot(dat_qc, features = marker_genes, pt.size = 0.3)
ggsave("./plots/B2MKO_Mix_QC2_marker_gene_plot2.png", width = 11, height = 17)

Idents(dat_qc) <- "seurat_cluster_res_0.5"
DimPlot(dat_qc, label = TRUE)
ggsave(filename = "./plots/B2MKO_mix_QC2_UMAP_seurat_cluster_res_0.5_labeled.png")


# Find markers of seurat_cluster_res_0.5
Idents(dat_qc) <- "seurat_cluster_res_0.5"
dat_qc_05_markers <- FindAllMarkers(dat_qc, only.pos = TRUE)
write.csv(dat_qc_05_markers, file = "./numerical_results/B2MKO_Mix_QC2_seurat05_cluster_markers.csv")

# -------------------------------------------------------------------
# Finalized Cell-type annotation based on seurat_cluster_res_0.5
# -------------------------------------------------------------------

# 0: NK
# 1: pDC
# 2: Treg
# 3: T_CD4
# 4: T_proliferating
# 5: T_CD4
# 6: T_CD8
# 7: NKT
# 8: gdT
# 9: B
# 10: B
# 11: T_CD8_activated/exhausted
# 12: cDC2
# 13: cDC1
# 14: Macrophages
# 15: pDC
# 16: B
# 17: Mast_cell
# 18: Plasma_cell

# Assign annotation to the Seurat object
dat_qc@meta.data$sub_celltype_l1 <- "Unknown"

dat_qc@meta.data$sub_celltype_l1[which(dat_qc@meta.data$seurat_cluster_res_0.5 %in% c(0))] <- "NK"
dat_qc@meta.data$sub_celltype_l1[which(dat_qc@meta.data$seurat_cluster_res_0.5 %in% c(1, 15))] <- "pDC"
dat_qc@meta.data$sub_celltype_l1[which(dat_qc@meta.data$seurat_cluster_res_0.5 %in% c(2))] <- "Treg"
dat_qc@meta.data$sub_celltype_l1[which(dat_qc@meta.data$seurat_cluster_res_0.5 %in% c(3, 5))] <- "CD4_T"
dat_qc@meta.data$sub_celltype_l1[which(dat_qc@meta.data$seurat_cluster_res_0.5 %in% c(4))] <- "Proliferating_T"
dat_qc@meta.data$sub_celltype_l1[which(dat_qc@meta.data$seurat_cluster_res_0.5 %in% c(6))] <- "CD8_T"
dat_qc@meta.data$sub_celltype_l1[which(dat_qc@meta.data$seurat_cluster_res_0.5 %in% c(11))] <- "CD8_T_activated"
dat_qc@meta.data$sub_celltype_l1[which(dat_qc@meta.data$seurat_cluster_res_0.5 %in% c(7))] <- "NKT"
dat_qc@meta.data$sub_celltype_l1[which(dat_qc@meta.data$seurat_cluster_res_0.5 %in% c(8))] <- "gdT"
dat_qc@meta.data$sub_celltype_l1[which(dat_qc@meta.data$seurat_cluster_res_0.5 %in% c(9, 10, 16, 18))] <- "B"
dat_qc@meta.data$sub_celltype_l1[which(dat_qc@meta.data$seurat_cluster_res_0.5 %in% c(12))] <- "cDC2"
dat_qc@meta.data$sub_celltype_l1[which(dat_qc@meta.data$seurat_cluster_res_0.5 %in% c(13))] <- "cDC1"
dat_qc@meta.data$sub_celltype_l1[which(dat_qc@meta.data$seurat_cluster_res_0.5 %in% c(14))] <- "Macrophage"
dat_qc@meta.data$sub_celltype_l1[which(dat_qc@meta.data$seurat_cluster_res_0.5 %in% c(17))] <- "Mast_cell"

# Final annotated UMAP
Idents(dat_qc) <- "sub_celltype_l1"
DimPlot(dat_qc, label = TRUE)
ggsave(filename = "./plots/B2MKO_mix_QC2_UMAP_celltype_label.png", width = 8)

saveRDS(dat_qc, file = "./objects/B2MKO_Mix_QC2_REpca_umap_clustered_annotated.rds")

