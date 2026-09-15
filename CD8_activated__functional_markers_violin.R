##### Assess activated CD8 T cell functional markers #####

# Load required packages
library(Seurat)
library(ggplot2)
library(SeuratExtend)
library(patchwork)

# Source helper script
source("./scripts/subset_reprocessing.R") # script used to subset Seurat object

# Load the processed dataset
dat_qc <- readRDS("./objects/B2MKO_Mix_QC2_REpca_umap_clustered_annotated.rds")
Idents(dat_qc) <- "sub_celltype_l1" # use annotated cell types as grouping identity

# Restrict to activated CD8 T cells
cells <- colnames(dat_qc)[dat_qc$sub_celltype_l1 %in% c("CD8_T_activated")]

# Marker genes to visualize
marker_genes <- c("Ifng", "Entpd1", "Cd69", "Tbx21", "Tnf", "Klrg1", "Gzmb", "Prf")
condition_cols <- c("#A8C2D0", "#0B1F2A80")

# Build the violin plot
p <- VlnPlot2(
  dat_qc,
  features = marker_genes,
  split.by = "condition",
  cells = cells,
  stat.method = "wilcox.test",
  hide.ns = TRUE,
  cols = condition_cols,
  box = TRUE,
  label = "p.signif"
) +
  facet_wrap(~feature, ncol = 8, nrow = 1, strip.position = "bottom") +
  theme(strip.text.x = element_text(angle = 0, hjust = 0.5))

print(p)

ggsave(
  "./plots/CD8_T_activated_functional_markers.pdf",
  width = 8,
  height = 3
)
