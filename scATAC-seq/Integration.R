# ============================ #
# scATAC-seq and scRNA-seq Integration (ArchR + Seurat v4)
# ============================ #

# IMPORTANT: Seurat and SeuratObject must be downgraded to v4 for compatibility

# Integrate scATAC-seq with scRNA-seq reference
proj.filter <- addGeneIntegrationMatrix(
  ArchRProj = proj.filter,                        # ArchR Project (scATAC)
  useMatrix = "GeneScoreMatrix",                  # Use gene scores
  matrixName = "GeneIntegrationMatrix",           # Output matrix name
  reducedDims = "Harmony",                        # Harmony-reduced dimensions
  seRNA = OSCC,                                   # Reference scRNA-seq Seurat object
  addToArrow = TRUE,
  groupRNA = "celltype",                          # scRNA annotation column
  nameCell = "predictedCell",
  nameGroup = "predictedGroup",
  nameScore = "predictedScore",
  force = TRUE
)

# Plot predicted groups on UMAP
plotEmbedding(proj.filter, colorBy = "cellColData", name = "predictedGroup", size = 0.05)

# Match clusters to predicted groups
cM <- confusionMatrix(proj.filter$Clusters, proj.filter$predictedGroup)
cM <- cM / Matrix::rowSums(cM)
labelNew <- colnames(cM)[apply(cM, 1, which.max)]
mapLabs <- cbind(rownames(cM), labelNew)

# Visualize confusion matrix heatmap
pheatmap::pheatmap(
  mat = as.matrix(cM),
  color = paletteContinuous("whiteBlue"),
  border_color = "black"
)

# Assign predicted group labels
proj.filter$Clusters_sub <- mapLabs[match(proj.filter$Clusters, mapLabs[,1]), 2]

# Check distribution by sample
table(proj.filter$Clusters_sub, proj.filter$Sample)

# Plot labeled clusters
plotEmbedding(proj.filter, colorBy = "cellColData", name = "Clusters_sub", embedding = "UMAP")

# ==== Feature plots for marker genes ====
getAvailableMatrices(proj.filter)
sc_marker <- c("MAST4", "PECAM1", "KRT14", "DCN", "CPA3", "RGS5", "LYZ", "NKG7",
               "MZB1", "CD3D", "CD4", "CD8A")

featurePlots <- plotEmbedding(
  ArchRProj = proj.filter,
  colorBy = "GeneIntegrationMatrix",
  name = sc_marker,
  embedding = "UMAP",
  imputeWeights = getImputeWeights(proj.filter)
)

# Customize ggplot themes for all feature plots
featurePlotsFormatted <- lapply(featurePlots, function(x) {
  x + guides(color = FALSE, fill = FALSE) +
    theme_ArchR(baseSize = 6.5) +
    theme(
      plot.margin = unit(c(0, 0, 0, 0), "cm"),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank()
    ) +
    scale_fill_viridis_c()
})

# Combine into a grid
do.call(cowplot::plot_grid, c(list(ncol = 4), featurePlotsFormatted))

# Save project
save(proj.filter, file = "projfilter20241112.RData")

# ==== Custom visualization using ggplot2 ====
# UMAPs colored by predicted group, sample, and cluster
atac.archr <- plotEmbedding(proj.filter, colorBy = "cellColData", name = "predictedGroup_ArchR")
atac.archr.emb <- as.data.frame(atac.archr$data)
atac.archr.emb$cell.type.archr <- gsub(".*:", "", sub("-", ":", atac.archr.emb$color))
atac.archr.emb$cell.type.archr <- factor(atac.archr.emb$cell.type.archr, levels = levels(as.factor(OSCC$celltype)))

# Sample
atac.sample <- plotEmbedding(proj.filter, colorBy = "cellColData", name = "Sample")
atac.sample.emb <- as.data.frame(atac.sample$data)
atac.sample.emb$sample <- gsub(".*:", "", sub("-", ":", atac.sample.emb$color))

# Clusters
atac.cluster <- plotEmbedding(proj.filter, colorBy = "cellColData", name = "Clusters")
atac.cluster.emb <- as.data.frame(atac.cluster$data)
atac.cluster.emb$cluster <- gsub(".*:", "", sub("-", ":", atac.cluster.emb$color))

# Merge into one dataframe
atac.emb.all <- cbind(
  atac.archr.emb[, c("x", "y", "cell.type.archr")],
  atac.sample.emb$sample,
  atac.cluster.emb$cluster
)
colnames(atac.emb.all) <- c("UMAP1", "UMAP2", "Predicted.Group.ArchR", "Sample", "ATAC_clusters")

# Plot: Predicted Group
p1 <- ggplot(atac.emb.all, aes(x = UMAP1, y = UMAP2, color = Predicted.Group.ArchR)) +
  geom_point(size = .1) +
  theme_classic() +
  ggtitle("scATAC-seq: Predicted Group (ArchR)") +
  theme(plot.title = element_text(face = "bold")) +
  theme(legend.key.size = unit(0.2, "cm"))

# Plot: Sample
p2 <- ggplot(atac.emb.all, aes(x = UMAP1, y = UMAP2, color = Sample)) +
  geom_point(size = .1) +
  theme_classic() +
  ggtitle("scATAC-seq: Sample") +
  theme(plot.title = element_text(face = "bold"))

# Plot: ATAC Clusters
p3 <- ggplot(atac.emb.all, aes(x = UMAP1, y = UMAP2, color = ATAC_clusters)) +
  geom_point(size = .1) +
  theme_classic() +
  ggtitle("scATAC-seq: ATAC Clusters") +
  theme(plot.title = element_text(face = "bold"))

# Combine all three
p1 + p2 + p3

# ============================ #
# GroupList Mapping for Refined Integration
# ============================ #

groupList <- SimpleList(
  T = SimpleList(
    ATAC = proj.filter$cellNames[proj.filter$Clusters %in% c("C5", "C7", "C8", "C9", "C10")],
    RNA = rownames(OSCC@meta.data[OSCC@meta.data$celltype %in% c("T"), ])
  ),
  B = SimpleList(
    ATAC = proj.filter$cellNames[proj.filter$Clusters %in% c("C3", "C4")],
    RNA = rownames(OSCC@meta.data[OSCC@meta.data$celltype %in% c("B cells"), ])
  ),
  Epithelial = SimpleList(
    ATAC = proj.filter$cellNames[proj.filter$Clusters %in% c("C15")],
    RNA = rownames(OSCC@meta.data[OSCC@meta.data$celltype %in% c("Epithelial cells"), ])
  ),
  Fibroblasts = SimpleList(
    ATAC = proj.filter$cellNames[proj.filter$Clusters %in% c("C19", "C20")],
    RNA = rownames(OSCC@meta.data[OSCC@meta.data$celltype %in% c("Fibroblasts"), ])
  ),
  Mast = SimpleList(
    ATAC = proj.filter$cellNames[proj.filter$Clusters %in% c("C11", "C12")],
    RNA = rownames(OSCC@meta.data[OSCC@meta.data$celltype %in% c("Mast cells"), ])
  ),
  Mural = SimpleList(
    ATAC = proj.filter$cellNames[proj.filter$Clusters %in% c("C18")],
    RNA = rownames(OSCC@meta.data[OSCC@meta.data$celltype %in% c("Mural cells"), ])
  ),
  Myeloid = SimpleList(
    ATAC = proj.filter$cellNames[proj.filter$Clusters %in% c("C13", "C14")],
    RNA = rownames(OSCC@meta.data[OSCC@meta.data$celltype %in% c("Myeloid cells"), ])
  ),
  NK = SimpleList(
    ATAC = proj.filter$cellNames[proj.filter$Clusters %in% c("C6")],
    RNA = rownames(OSCC@meta.data[OSCC@meta.data$celltype %in% c("NK"), ])
  ),
  Plasma = SimpleList(
    ATAC = proj.filter$cellNames[proj.filter$Clusters %in% c("C1", "C2")],
    RNA = rownames(OSCC@meta.data[OSCC@meta.data$celltype %in% c("Plasma cells"), ])
  ),
  Endothelial = SimpleList(
    ATAC = proj.filter$cellNames[proj.filter$Clusters %in% c("C16", "C17")],
    RNA = rownames(OSCC@meta.data[OSCC@meta.data$celltype %in% c("Endothelial cells"), ])
  )
)

# Apply integration with group list
sc <- addGeneIntegrationMatrix(
  ArchRProj = proj.filter,
  useMatrix = "GeneScoreMatrix",
  matrixName = "GeneIntegrationMatrix",
  reducedDims = "Harmony",
  seRNA = OSCC,
  addToArrow = TRUE,
  force = TRUE,
  groupList = groupList,
  groupRNA = "celltype",
  nameCell = "predictedCell_ArchR",
  nameGroup = "predictedGroup_ArchR",
  nameScore = "predictedScore_ArchR",
  transferParams = list(dims = 1:50)
)

# Visualization after refined integration
plotEmbedding(sc, colorBy = "cellColData", name = "predictedGroup_ArchR")

# Label clusters with predicted groups
cM <- confusionMatrix(sc$Clusters, sc$predictedGroup_ArchR)
cM <- cM / Matrix::rowSums(cM)
labelNew <- colnames(cM)[apply(cM, 1, which.max)]
mapLabs <- cbind(rownames(cM), labelNew)

sc$Clusters_cell <- mapLabs[match(sc$Clusters, mapLabs[,1]), 2]

# Final UMAP plots
p1 <- plotEmbedding(sc, colorBy = "cellColData", name = "Clusters_cell", embedding = "UMAP")
p2 <- plotEmbedding(proj.filter, colorBy = "cellColData", name = "predictedGroup", embedding = "UMAP")

# Combine side-by-side
p1 | p2
