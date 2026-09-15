# NK chemokine expression across cell types

# Load required packages
library(Seurat)
library(ggplot2)
library(SeuratExtend)
library(patchwork)

# Load annotated object
 dat_qc <- readRDS("./objects/B2MKO_Mix_QC2_REpca_umap_clustered_annotated.rds")
Idents(dat_qc) <- "sub_celltype_l1"


# Assess NK chemokine CXCL9 and CXCL10 expression across annotated cell types

VlnPlot2(dat_qc, features = c("Cxcl9", "Cxcl10"), split.by = "condition", ncol = 1, stat.method = "wilcox.test", 
         hide.ns = TRUE, cols = c("#A8C2D0", "#0B1F2A80"), box = TRUE)

ggsave(
    filename = "./plots/B2MKO_mix_QC2_NK_trafficking_chemokines_by_condition.pdf", 
    height = 5, 
    width = 8)





