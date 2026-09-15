# Initial QC and subsetting for an object containing B2MKO and Mix samples only

# Load required packages
library(Seurat)
library(ggplot2)

sobj <- readRDS("./objects/seurat_object_filtered_with_metadata.rds")


# Filter cells based on feature counts and mitochondrial content

sobj[["percent.mt"]] <- PercentageFeatureSet(sobj, pattern = "^mt")
VlnPlot(sobj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
ggsave(filename = "./plots/sobj_all_labeled_cells_nfeaturerna_ncountrna_percentmt.png", width = 8, height = 5)

plot1 <- FeatureScatter(sobj, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(sobj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

ggsave(filename = "./plots/sobj_all_labeled_cells_ncountrna_vs_percentmt.png", plot = plot1)
ggsave(filename = "./plots/sobj_all_labeled_cells_ncountrna_vs_nfeaturerna.png", plot = plot2)


# Filter out cells based on RNA feature count and percent mitochondria content. Save the filtered cells to a new object called sobj_qc
sobj_qc <- subset(sobj, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 10)
saveRDS(sobj_qc, file = "./objects/All_samples_QCed_object.rds")


# Subset for B2MKO and Mix #
sub_b2mko_mix <- subset(sobj_qc, subset = condition %in% c("B2MKO", "Mix"))
saveRDS(sub_b2mko_mix, file = "./objects/b2mko_mix_subset_from_all_samples_QCed_object.rds")


