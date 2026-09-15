# Pseudobulk differential gene analysis between immune cell types in B2MKO vs. Mix tumors

library(Seurat)
library(dplyr)
library(tidyverse)


# -------------------------------------------------------------------
# Load data and set cell type identities
# -------------------------------------------------------------------
dat_qc <- readRDS("./objects/B2MKO_Mix_QC2_REpca_umap_clustered_annotated.rds")

Idents(dat_qc) <- "sub_celltype_l1"


# -------------------------------------------------------------------
# Pseudobulking based on condition, MouseID, and sub_celltype_l1
# -------------------------------------------------------------------

pseudo_dat_qc <- AggregateExpression(dat_qc, assays = "RNA", return.seurat = T, group.by = c("condition", "MouseID", "sub_celltype_l1"))

Cells(pseudo_dat_qc)

pseudo_dat_qc$celltype.condition <- paste(pseudo_dat_qc$sub_celltype_l1, pseudo_dat_qc$condition, sep = "_")

# -------------------------------------------------------------------
# Perform DEG on NK cells
# -------------------------------------------------------------------

Idents(pseudo_dat_qc) <- "celltype.condition"

bulk.nk.de <- FindMarkers(object = pseudo_dat_qc, 
                            ident.1 = "NK_Mix", 
                            ident.2 = "NK_B2MKO",
                            test.use = "DESeq2")

write.csv(bulk.nk.de, file = "./numerical_results/Mix_vs_B2MKO_nk_aggregated_counts_deg.csv")


# -------------------------------------------------------------------
# Perform DEG on activated CD8 T cells
# -------------------------------------------------------------------

Idents(pseudo_dat_qc) <- "celltype.condition"

bulk.cd8_t_activated.de <- FindMarkers(object = pseudo_dat_qc, 
                                       ident.1 = "CD8-T-activated_Mix", 
                                       ident.2 = "CD8-T-activated_B2MKO",
                                       test.use = "DESeq2")

write.csv(bulk.cd8_t_activated.de, file = "./numerical_results/Mix_vs_B2MKO_cd8_t_activated_aggregated_counts_deg.csv")


# -------------------------------------------------------------------
# Perform DEG on non activated CD8 T cells
# -------------------------------------------------------------------

Idents(pseudo_dat_qc) <- "celltype.condition"

bulk.cd8_t_activated.de <- FindMarkers(object = pseudo_dat_qc, 
                                       ident.1 = "CD8-T_Mix", 
                                       ident.2 = "CD8-T_B2MKO",
                                       test.use = "DESeq2")

write.csv(bulk.cd8_t_activated.de, file = "./numerical_results/Mix_vs_B2MKO_cd8_aggregated_counts_deg.csv")
  
  


  
  
  
  
  