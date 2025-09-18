library(Seurat)        
library(harmony)       
library(MAESTRO)      
# Normalize UMI counts 
  combined_seurat <- NormalizeData(combined_seurat, normalization.method = "LogNormalize", scale.factor = 1e4)
  
  # Identify top 2,000 highly variable genes for PCA
  combined_seurat <- FindVariableFeatures(combined_seurat, selection.method = "vst", nfeatures = 2000)
  
  # Scale data
  all.genes <- rownames(combined_seurat)
  combined_seurat <- ScaleData(combined_seurat, features = all.genes)
  
  # PCA with top 50 principal components
  combined_seurat <- RunPCA(combined_seurat, features = VariableFeatures(combined_seurat), npcs = 50, verbose = FALSE)
  combined_seurat <- RunHarmony(combined_seurat, group.by.vars = "sample")
  
  # Build the SNN graph and cluster using the top 25 PCs; resolution = 0.4
  combined_seurat <- FindNeighbors(combined_seurat, reduction = "pca", dims = 1:25)
  combined_seurat <- FindClusters(combined_seurat, resolution = 0.4)
  combined_seurat <- RunUMAP(combined_seurat, reduction = "pca", dims = 1:25)
  
  # Identify marker genes per cluster using MAESTRO (logFC > 1, padj < 0.05)
  markers <- MAESTRO::FindAllMarkersMAESTRO(
    combined_seurat,
    logfc.threshold = 1.0
  )
