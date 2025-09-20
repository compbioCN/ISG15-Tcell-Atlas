# -------------------------------
# 1. Identify marker peaks
# -------------------------------

# Filter marker peaks: FDR <= 0.05 and Log2FC >= 0.5
markerList <- getMarkers(markerPeaks, cutOff = "FDR <= 0.05 & Log2FC >= 0.5")

# Convert full list to one flat data frame for downstream extraction
all_markers_df <- as.data.frame(markerList)

# Extract marker peak tables by cluster
marker_Peaks <- markerList@listData
for (i in seq_along(marker_Peaks)) {
  marker_Peaks[[i]] <- as.data.frame(marker_Peaks[[i]])
}

# Annotate each marker table with corresponding cluster name
for (i in seq_along(marker_Peaks)) {
  if (nrow(marker_Peaks[[i]]) == 0) {
    marker_Peaks[[i]][1, ] <- names(marker_Peaks[i])
  } else {
    marker_Peaks[[i]]$Clusters_cell <- names(marker_Peaks[i])
  }
}

# -------------------------------
# 2. Plot heatmap of marker peaks
# -------------------------------
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
  width = 10,
  height = 15,
  ArchRProj = sc,
  addDOC = FALSE
)

# -------------------------------
# 3. Subset peaks for specific clusters
# -------------------------------
Peaks_CD8ISG <- all_markers_df[all_markers_df$group_name == "CD8T-ISG15", ]
Peaks_CD4ISG <- all_markers_df[all_markers_df$group_name == "CD4T-ISG15", ]

# -------------------------------
# 4. MA plot for one group
# -------------------------------
MAplot <- markerPlot(
  seMarker = markerPeaks,
  name = "CD8T-ISG15",
  cutOff = "FDR <= 0.05 & Log2FC >= 3"
)

# Display and save
print(MAplot)

plotPDF(
  MAplot,
  name = "CD8T_Peak-Marker-MAplot",
  width = 6,
  height = 6,
  ArchRProj = sc,
  addDOC = FALSE
)

ggsave(filename = "MAplot.pdf", plot = MAplot, width = 6, height = 6)

# -------------------------------
# 5. Genome browser-style track plot (for single gene)
# -------------------------------

# Plot accessibility around a selected gene (e.g., "CBLB")
browser_track <- plotBrowserTrack(
  ArchRProj = sc,
  groupBy = "celltype",
  geneSymbol = "CBLB",
  upstream = 200000,
  downstream = 50000,
  tileSize = 250,
  minCells = 250
)

# Save to PDF
plotPDF(
  plotList = browser_track,
  name = "Plot-Tracks-Marker-Genes.pdf",
  ArchRProj = sc,
  addDOC = FALSE,
  width = 5,
  height = 5
)

# -------------------------------
# 6. Custom genomic region track (optional)
# -------------------------------
library(GenomicRanges)
region_gr <- GRanges(
  seqnames = "chr14",
  ranges = IRanges(start = 71345821, end = 71505821)  # Define region of interest
)

# Draw custom genomic region accessibility track
region_track <- plotBrowserTrack(
  ArchRProj = sc,
  groupBy = "celltype",
  region = region_gr,
  tileSize = 250,
  minCells = 250
)

# Display
grid.draw(region_track)
