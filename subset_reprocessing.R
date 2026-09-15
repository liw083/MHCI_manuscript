# Define function for processing cell subsets

# Re-select features, rescale, and PCA
subset_hvg_rescale_pca <- function(sobj, subset_colname, subset_to_include, output_folder_name, output_subset_name){
  
  # subset based on cells of interest
  sub <- subset(sobj, cells = rownames(sobj@meta.data[sobj@meta.data[[subset_colname]] %in% subset_to_include, ]))
  print(unique(sub@meta.data[[subset_colname]]))
  
  # Identify HVG
  sub <- FindVariableFeatures(sub, selection.method = "vst", nfeatures = 2000)
  top10 <- head(VariableFeatures(sub), 10)
  plot1 <- VariableFeaturePlot(sub)
  plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
  plot1 + plot2
  output_path <- paste0(output_folder_name, output_subset_name)
  ggsave(filename = paste0(output_path,"_hvg_scatter.png"), plot = plot2)
  
  # Rescale
  all.genes <- rownames(sub)
  sub <- ScaleData(sub, features = all.genes)
  
  # Rerun pca
  sub <- RunPCA(sub, features = VariableFeatures(sub))
  VizDimLoadings(sub, dims = 1:2, reduction = "pca")
  ggsave(filename = paste0(output_path,"_pca_VizDimLoading.png"))
  
  DimPlot(sub, reduction = "pca", group.by = "condition")
  ggsave(filename = paste0(output_path,"_pca_scatter_condition.png"))
  
  ElbowPlot(sub, ndims = 50)
  ggsave(filename = paste0(output_path, "_filtered_pca_elbow.png"))
  
  return(sub)
}


# Reperform neighorhood identification, cluster identification, and UMAP
subset_recluster_umap <- function(sobj, n_pc, umap_res, output_folder_name, output_subset_name){
  
  set.seed(10101)
  
  # recluster cells
  sobj <- FindNeighbors(sobj, reduction = "pca", dims = 1:n_pc)
  
  for (res in c(0.1, 0.2, 0.5, 0.8, 1.0, 1.5, 2.0)){
    
    cluster_name <- paste0("subcluster_", res)
    
    sobj <- FindClusters(sobj, cluster.name = cluster_name, resolution = res)
  }
  
  # rerun umap
  sobj <- RunUMAP(sobj, dims = 1:n_pc)
  
  Idents(sobj) <- paste0("subcluster_", umap_res)
  output_path <- paste0(output_folder_name, output_subset_name)
  
  DimPlot(sobj, reduction = "umap")
  ggsave(filename = paste0(output_path, "_umap_res", umap_res, ".png"), width = 5, height = 5)
  
  DimPlot(sobj, reduction = "umap", label = TRUE)
  ggsave(filename = paste0(output_path, "_umap_res", umap_res, "_labeled.png"), width = 5, height = 5)
  
  DimPlot(sobj, reduction = "umap", group.by = "MouseID")
  ggsave(filename = paste0(output_path, "_umap_MouseID_merged.png"))
  
  DimPlot(sobj, reduction = "umap", split.by = "MouseID")
  ggsave(filename = paste0(output_path, "_umap_MouseID_split.png"), width = 30, height = 4)
  
  DimPlot(sobj, reduction = "umap", group.by = "condition")
  ggsave(filename = paste0(output_path, "_umap_condition_merged.png"))
  
  DimPlot(sobj, reduction = "umap", split.by = "condition")
  ggsave(filename = paste0(output_path, "_umap_condition_split.png"), width = 20, height = 4)
  
  return(sobj)
}
