# Required package versions
# Matrix >= 1.6.5
# spatstat.explore >= 3.3.2
# Seurat / SeuratObject >= 5.0.0

# ---------- Dimension Reduction & Clustering ----------

# Load pre-filtered ArchRProject object
proj.filter <- readRDS("Save-ArchR-Project.rds")

# Step 1: Iterative LSI (Latent Semantic Indexing)
proj.filter <- addIterativeLSI(
  ArchRProj = proj.filter,
  useMatrix = "TileMatrix",          # Use the accessibility matrix
  name = "IterativeLSI",             # Name of the reduced dimension
  iterations = 4,                    # Number of LSI iterations
  clusterParams = list(
    resolution = 4,
    sampleCells = 10000,
    n.start = 10
  ),
  varFeatures = 50000,              # Number of variable features used
  dimsToUse = 1:20,
  force = TRUE,
  seed = 10
)

# Step 2: Batch Correction using Harmony
proj.filter <- addHarmony(
  ArchRProj = proj.filter,
  reducedDims = "IterativeLSI",
  name = "Harmony",                  # Name of harmony output
  groupBy = "Sample",                # Batch correction by sample
  force = TRUE,
  theta = 2                          # Regularization parameter (higher = stronger correction)
)

# Step 3: Clustering
proj.filter <- addClusters(
  input = proj.filter,
  reducedDims = "Harmony",          # Use harmony corrected dimensions
  method = "Seurat",                # Clustering method (e.g., Seurat, Louvain, etc.)
  name = "Clusters",                # Name of metadata column
  resolution = 0.5,                 # Resolution for clustering granularity
  force = TRUE,
  seed = 11
)

# Step 4: Cluster composition by sample (confusion matrix)
cM <- confusionMatrix(
  paste0(proj.filter$Clusters),
  paste0(proj.filter$Sample)
)
cM <- cM / Matrix::rowSums(cM)
p <- pheatmap::pheatmap(
  mat = as.matrix(cM),
  color = paletteContinuous("whiteBlue"),
  border_color = "black"
)

# Step 5: UMAP embedding
proj.filter <- addUMAP(
  ArchRProj = proj.filter,
  reducedDims = "Harmony",          # Use harmony for UMAP embedding
  name = "UMAP",
  nNeighbors = 3,
  minDist = 1e-11,                  # Extremely small value for tight clustering
  metric = "cosine",
  force = TRUE,
  seed = 12
)

# Save ArchR project
saveArchRProject(proj.filter, outputDirectory = "./")

# Step 6: Visualize UMAP by Clusters and Sample
p7 <- plotEmbedding(
  ArchRProj = proj.filter,
  colorBy = "cellColData",
  name = "Clusters",
  embedding = "UMAP",
  size = 0.2
)

p8 <- plotEmbedding(
  ArchRProj = proj.filter,
  colorBy = "cellColData",
  name = "Sample",
  embedding = "UMAP",
  size = 0.2
)

# Combine plots side by side
p7 | p8
