# ArchR scATAC-seq Full Preprocessing & QC Pipeline

# ===== Step 0: Environment Setup =====
# Version Requirements:
# Matrix >= 1.6.5, spatstat.explore == 3.3.2, Seurat >= 5.0, ggplot2 >= 3.4

library(ArchR)
library(pheatmap)
library(Rsamtools)
library(scran)
library(scater)
library(dplyr)
library(Seurat)
library(patchwork)
library(SingleCellExperiment)
library(ComplexHeatmap)
library(ggplot2)
library(stringr)
library(EnsDb.Hsapiens.v86)
library(viridis)
library(scuttle)

# ===== Step 1: Set genome & working directory =====
addArchRGenome("hg38")
setwd("/your/project/path")  # <-- Change this!

# ===== Step 2: Input fragments and sample names =====
input.file.list <- c(
  "path/to/P05_OLK/outs/fragments.tsv.gz",
  "path/to/P05_Tumor/outs/fragments.tsv.gz",
  "path/to/P10_OLK/outs/fragments.tsv.gz",
  "path/to/P10_Tumor/outs/fragments.tsv.gz",
  "path/to/P11_Tumor/outs/fragments.tsv.gz",
  "path/to/P12_OLP/outs/fragments.tsv.gz",
  "path/to/P12_Tumor/outs/fragments.tsv.gz"
)

sampleNames <- c("P05_OLK", "P05_Tumor", "P10_OLK", "P10_Tumor", "P11_Tumor", "P12_OLP", "P12_Tumor")

# ===== Step 3: Create Arrow files =====
ArrowFiles <- createArrowFiles(
  inputFiles = input.file.list,
  sampleNames = sampleNames,
  minTSS = 4,
  minFrags = 1000,
  addTileMat = TRUE,
  addGeneScoreMat = TRUE,
  excludeChr = c("chrM", "chrY", "chrX")
)

# ===== Step 4: Doublet Filtering =====
doubScores <- addDoubletScores(
  input = ArrowFiles,
  k = 20,
  knnMethod = "UMAP",
  useMatrix = "TileMatrix",
  nTrials = 5,
  LSIMethod = 1,
  scaleDims = FALSE,
  corCutOff = 0.75,
  UMAPParams = list(n_neighbors = 30, min_dist = 0.3, metric = "cosine", verbose = TRUE)
)

# ===== Step 5: Create ArchR Project =====
projncov <- ArchRProject(
  ArrowFiles = ArrowFiles,
  outputDirectory = "ATAC_out",
  copyArrows = TRUE
)

# ===== Step 6: QC Plots (TSSEnrichment & log10(nFrags)) =====
p1 <- plotGroups(projncov, groupBy = "Sample", colorBy = "cellColData", name = "TSSEnrichment", plotAs = "ridges")
p2 <- plotGroups(projncov, groupBy = "Sample", colorBy = "cellColData", name = "TSSEnrichment", plotAs = "violin", alpha = 0.4, addBoxPlot = TRUE)
p3 <- plotGroups(projncov, groupBy = "Sample", colorBy = "cellColData", name = "log10(nFrags)", plotAs = "ridges")
p4 <- plotGroups(projncov, groupBy = "Sample", colorBy = "cellColData", name = "log10(nFrags)", plotAs = "violin", alpha = 0.4, addBoxPlot = TRUE)

combined_plot <- (p1 + p2) / (p3 + p4)

# ===== Step 7: Sample-specific Density QC Plots =====
QC_plot <- list()
for (i in seq_along(sampleNames)) {
  sc <- projncov[projncov$Sample == sampleNames[i], ]
  df <- getCellColData(sc, select = c("nFrags", "TSSEnrichment"))
  df <- data.frame("nFrags" = log10(df$nFrags), "TSSEnrichment" = df$TSSEnrichment)

  p <- ggPoint(
    x = df[,1], y = df[,2],
    colorDensity = TRUE,
    continuousSet = "sambaNight",
    xlabel = "Log10 Unique Fragments",
    ylabel = "TSS Enrichment",
    xlim = c(log10(500), quantile(df[,1], 0.99)),
    ylim = c(0, quantile(df[,2], 0.99))
  ) +
    geom_hline(yintercept = 4, lty = "dashed") +
    geom_vline(xintercept = 3, lty = "dashed") +
    ggtitle(sampleNames[i])

  QC_plot[[i]] <- p
}
wrap_plots(QC_plot, ncol = 4)
plotPDF(QC_plot, name = "QC-Sample-density.pdf", ArchRProj = projncov, addDOC = FALSE, width = 12, height = 6)

# ===== Step 8: Fragment Size Histogram & TSS Profile =====
p5 <- plotFragmentSizes(projncov) + ggtitle("Fragment Size Histogram")
p6 <- plotTSSEnrichment(projncov) + ggtitle("TSS Enrichment")
plotPDF(p5, p6, name = "QC-Sample-FragSizes-TSSProfile.pdf", ArchRProj = projncov, addDOC = FALSE)

# ===== Step 9: Add Custom Metadata (optional) =====
atac_samples <- data.frame(
  Sample = c("P05_OLK", "P10_OLK", "P12_OLP"),
  Age = c(6, 10, 12),
  Sex = c("F", "F", "M")
)
projncov$Age <- as.numeric(atac_samples$Age[match(projncov$Sample, atac_samples$Sample)])
projncov$Sex <- atac_samples$Sex[match(projncov$Sample, atac_samples$Sample)])

# ===== Step 10: Save/Reload Project =====
saveArchRProject(projncov, outputDirectory = "./", load = FALSE)
projncov <- loadArchRProject("projncov")

# ===== Step 11: Remove Doublets =====
proj.filter <- filterDoublets(projncov)

# ===== Step 12: Filter TSS Outliers =====
tss_outliers <- list()
tss_outliers_names <- list()
for (i in seq_along(sampleNames)) {
  sample_i <- sampleNames[i]
  tss_enrich <- proj.filter$TSSEnrichment[proj.filter$Sample == sample_i]
  tss_outliers[[sample_i]] <- isOutlier(tss_enrich, nmads = 1, type = "lower")
  tss_outliers_names[[sample_i]] <- proj.filter$cellNames[proj.filter$Sample == sample_i][tss_outliers[[sample_i]]]
}
tss_outliers_names <- unlist(tss_outliers_names)
proj.filter <- proj.filter[!proj.filter$cellNames %in% tss_outliers_names, ]
