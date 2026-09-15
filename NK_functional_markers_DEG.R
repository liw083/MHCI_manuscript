# Assess NK cell functional markers and perform differential gene expression on all NK based on tumor types

library(Seurat)
library(ggplot2)
library(SeuratExtend)
library(patchwork)
Idents(dat_qc) <- "sub_celltype_l1"

cells <- colnames(dat_qc)[dat_qc$sub_celltype_l1 %in% c("NK")]

# -------------------------------------------------------------------
# Assess NK cell activation markers
# -------------------------------------------------------------------

VlnPlot2(dat_qc, features = c("Klrg1", "Cd69", "Cd226", "Il2ra", "Icos", "Ly6a"), 
         split.by = "condition", cells = cells, stat.method = "wilcox.test", 
         hide.ns = TRUE, cols = c("#A8C2D0", "#0B1F2A80"), box = TRUE, label = "p.signif") + 
  facet_wrap(~feature, ncol = 6, nrow = 1, strip.position = "bottom") + 
  theme(strip.text.x = element_text(angle = 0, hjust = 0.5))

ggsave("./plots/NK_subset__activation_markers.pdf", width = 8, height = 3)


# -------------------------------------------------------------------
# Assess NK cell perforin/granzymes #####
# -------------------------------------------------------------------

VlnPlot2(dat_qc, features = c("Prf1", "Gzmk","Gzma", "Gzmb"), 
         split.by = "condition", cells = cells, stat.method = "wilcox.test", p.adjust.method = "holm",
         hide.ns = TRUE, cols = c("#A8C2D0", "#0B1F2A80"), box = TRUE, label = "p.signif", show.mean = FALSE) + 
         facet_wrap(~feature, ncol = 6, nrow = 1, strip.position = "bottom") + 
  theme(strip.text.x = element_text(angle = 0, hjust = 0.5)
)

ggsave("./plots/NK_subset__perforin_granzymes.pdf", width = 6, height = 3)


# -------------------------------------------------------------------
# Assess NK cell proliferation #####
# -------------------------------------------------------------------

VlnPlot2(dat_qc, features = c("Mki67", "Top2a", "Cdk1"), 
         split.by = "condition", cells = cells, stat.method = "wilcox.test", p.adjust.method = "holm",
         hide.ns = TRUE, cols = c("#A8C2D0", "#0B1F2A80"), box = TRUE, label = "p.signif", show.mean = FALSE) + 
  facet_wrap(~feature, ncol = 6, nrow = 1, strip.position = "bottom") + 
  theme(strip.text.x = element_text(angle = 0, hjust = 0.5)
  )

ggsave("./plots/NK_subset__proliferation_markers.pdf", width = 4, height = 3)


# -------------------------------------------------------------------
# Assess NK cell immune checkpoints #####
# -------------------------------------------------------------------
VlnPlot2(dat_qc, features = c("Tigit", "Lag3", "Ctla4", "Pdcd1"), 
         split.by = "condition", cells = cells, stat.method = "wilcox.test", p.adjust.method = "holm",
         hide.ns = TRUE, cols = c("#A8C2D0", "#0B1F2A80"), box = TRUE, label = "p.signif", show.mean = FALSE) + 
  facet_wrap(~feature, ncol = 6, nrow = 1, strip.position = "bottom") + 
  theme(strip.text.x = element_text(angle = 0, hjust = 0.5)
  )

ggsave("./plots/NK_subset__immune_checkpoints.pdf", width = 6, height = 3)
