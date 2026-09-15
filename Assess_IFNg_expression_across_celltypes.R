##### Assess IFNg expression across annotated cell types #####

# Load required packages
library(Seurat)
library(ggplot2)
library(SeuratExtend)

# Load processed dataset
 dat_qc <- readRDS("./objects/B2MKO_Mix_QC2_REpca_umap_clustered_annotated.rds")

# Use annotated cell types as grouping identity
Idents(dat_qc) <- "sub_celltype_l1"

# Dot plot across all annotated cell types
p_dot <- DotPlot(dat_qc, features = "Ifng")
print(p_dot)

ggsave(
  filename = "./plots/B2MKO_mix_QC2_IFNg_across_celltypes.pdf"
)

# -------------------------------------------------------------
# Violin plot for IFNg expression in activated CD8 T cells
# -------------------------------------------------------------

cells <- colnames(dat_qc)[dat_qc$sub_celltype_l1 %in% c("CD8_T_activated")]

p_vln <- VlnPlot2(
  dat_qc,
  features = "Ifng",
  split.by = "condition",
  cells = cells,
  stat.method = "wilcox.test",
  hide.ns = TRUE,
  cols = c("#A8C2D0", "#0B1F2A80"),
  box = TRUE
)

print(p_vln)

ggsave(
  "./plots/B2MKO_mix_QC2_IFNg_CD8_activated_B2MKO_vs_Mix_vlnplot.pdf",
  width = 3
)


