# ================================
# scATAC-seq: Peak Calling & Marker Peak Analysis
# ================================

library(ArchR)
library(ggplot2)
library(GenomicRanges)
library(rtracklayer)
library(grid)

# ===== Step 1: Create Group Coverages =====
sc <- addGroupCoverages(
  ArchRProj = projT,
  groupBy = "celltype",
  minCells = 50,
  force = TRUE
)
table(sc$celltype)

# ===== Step 2: Call Reproducible Peaks (MACS2 required) =====
pathToMacs2 <- "/home/Shared/miniconda3/envs/macs2/bin/macs2"  # Path to MACS2 binary
sc <- addReproduciblePeakSet(
  ArchRProj = sc,
  groupBy = "celltype",
  peaksPerCell = 500,
  pathToMacs2 = pathToMacs2,
  force = TRUE
)

# ===== Step 3: Add Peak Matrix =====
sc <- addPeakMatrix(ArchRProj = sc, force = TRUE)
getAvailableMatrices(sc)
# Expected output:
# [1] "GeneIntegrationMatrix" "GeneScoreMatrix" "PeakMatrix" "TileMatrix"

# ===== Step 4: Identify Marker Peaks =====
use_groups <- unique(sc$celltype)
markerPeaks <- getMarkerFeatures(
  ArchRProj = sc,
  useMatrix = "PeakMatrix",
  groupBy = "celltype",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)

# Retrieve markers as a list
markerList <- getMarkers(markerPeaks, cutOff = "FDR <= 0.05 & Log2FC >= 0.5")
A <- as.data.frame(markerList)
marker_Peaks <- markerList@listData

# Format marker peak data
for (x in seq_along(marker_Peaks)) {
  marker_Peaks[[x]] <- as.data.frame(marker_Peaks[[x]])
}
for (x in seq_along(marker_Peaks)) {
  if (nrow(marker_Peaks[[x]]) == 0) {
    marker_Peaks[[x]][1, ] <- names(marker_Peaks[x])
  } else {
    marker_Peaks[[x]]$Clusters_cell <- names(marker_Peaks[x])
  }
}

# ===== Step 5: Marker Peak Heatmap =====
heatmapPeaks <- plotMarkerHeatmap(
  seMarker = markerPeaks,
  cutOff = "FDR <= 0.05 & Log2FC >= 3",
  labelMarkers = NULL,
  binaryClusterRows = TRUE,
  transpose = FALSE
)

plotPDF(
  heatmapPeaks,
  name = "Peak-Marker-Heatmap",
  width = 10, height = 15,
  ArchRProj = sc,
  addDOC = FALSE
)

# ===== Step 6: Extract Peaks for Specific Cell Types =====
Peaks_CD8ISG <- A[A$group_name %in% c("CD8T-ISG15") & A$FDR < 0.05, ]
Peaks_CD4ISG <- A[A$group_name %in% c("CD4T-ISG15") & A$FDR < 0.05, ]

# ===== Step 7: MA Plot for CD8T-ISG15 =====
MAplot <- markerPlot(
  seMarker = markerPeaks,
  name = "CD8T-ISG15",
  cutOff = "FDR <= 0.05 & Log2FC >= 3"
)
print(MAplot)

plotPDF(
  MAplot,
  name = "CD8T_Peak-Marker-MAplot",
  width = 6, height = 6,
  ArchRProj = sc,
  addDOC = FALSE
)

ggsave(
  filename = "MAplot.pdf",
  plot = print(MAplot),
  width = 6, height = 6
)

# ===== Step 8: Browser Track for Specific Genes (e.g., CBLB) =====
p <- plotBrowserTrack(
  ArchRProj = sc,
  groupBy = "celltype",
  geneSymbol = "CBLB",      # Gene of interest
  upstream = 200000,
  downstream = 50000,
  tileSize = 250,
  minCells = 250
)

plotPDF(
  plotList = p,
  name = "Plot-Tracks-Marker-Genes.pdf",
  ArchRProj = sc,
  width = 5, height = 5,
  addDOC = FALSE
)

# ===== Step 9: Plot by Genomic Region (Custom GRanges) =====
region_gr <- GRanges(
  seqnames = "chr14",
  ranges = IRanges(start = 71345821, end = 71505821)  # Example region
)

p <- plotBrowserTrack(
  ArchRProj = sc,
  groupBy = "celltype",
  region = region_gr,
  tileSize = 250,
  minCells = 250
)

grid.draw(p)

# ===== Step 10: Optional - Load External Gene Annotation (e.g., GTF) =====
gtf_file <- "/home/guile/cellranger/genes.gtf"
gtf_data <- import(gtf_file)
feature_annotations <- getFeatureAnnotations(sc)

# View a few rows
head(as.data.frame(gtf_data))
head(feature_annotations)

# Final filtered peak tables:
# - Peaks_CD8ISG
# - Peaks_CD4ISG
# ===============================================
# Step 11: Annotate Marker Peaks Using GTF Genes
# ===============================================

library(rtracklayer)
library(GenomicRanges)
library(dplyr)

# ----- 1. Load GTF annotation -----
gtf_file <- "/home/guile/cellranger/genes.gtf"  # Replace with correct path
gtf_data <- import(gtf_file)  # Read GTF file

# Filter for gene-level annotations only
gtf_genes <- gtf_data[gtf_data$type == "gene"]

# Create GRanges object for gene annotation
gtf_gr <- GRanges(
  seqnames = seqnames(gtf_genes),
  ranges = IRanges(start = start(gtf_genes), end = end(gtf_genes)),
  strand = strand(gtf_genes),
  symbol = mcols(gtf_genes)$gene_name  # Extract gene symbol
)

# ----- 2. Convert Peaks to GRanges -----
peaks_gr <- GRanges(
  seqnames = Peaks_CD4ISG$seqnames,
  ranges = IRanges(start = Peaks_CD4ISG$start, end = Peaks_CD4ISG$end),
  Log2FC = Peaks_CD4ISG$Log2FC,
  FDR = Peaks_CD4ISG$FDR
)

# ----- 3. Find Overlaps Between Peaks and Genes -----
overlaps <- findOverlaps(peaks_gr, gtf_gr)

# Extract overlapping peak-gene pairs
peak_symbols <- data.frame(
  seqnames = seqnames(peaks_gr[queryHits(overlaps)]),
  start = start(peaks_gr[queryHits(overlaps)]),
  end = end(peaks_gr[queryHits(overlaps)]),
  Log2FC = mcols(peaks_gr[queryHits(overlaps)])$Log2FC,
  FDR = mcols(peaks_gr[queryHits(overlaps)])$FDR,
  gene_symbol = mcols(gtf_gr[subjectHits(overlaps)])$symbol
)

# Optional: preview
head(peak_symbols)

# ----- 4. Merge Annotation Back to Peaks Table -----
Peaks_CD4ISG_annotated <- dplyr::left_join(
  Peaks_CD4ISG,
  peak_symbols,
  by = c("seqnames", "start", "end")
)

# Remove unmatched (NA) entries
Peaks_CD4ISG_annotated <- na.omit(Peaks_CD4ISG_annotated)

# Final output: Annotated Peaks with Gene Symbols
head(Peaks_CD4ISG_annotated)
