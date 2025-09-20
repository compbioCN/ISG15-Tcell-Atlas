# =======================================
# CellChat Analysis: Secreted Signaling
# =======================================

# Load required libraries
library(CellChat)
library(ggplot2)

# ----------------------------
# Step 1: Prepare CellChat Data
# ----------------------------

# Extract data from existing CellChat object
data.input <- cellchat[["RNA"]]@data              # Gene expression matrix
meta <- cellchat@meta.data                        # Cell annotations

# Create new CellChat object
cellchatC <- createCellChat(object = data.input, meta = meta, group.by = "celltype")

# Load human signaling database
CellChatDB <- CellChatDB.human
showDatabaseCategory(CellChatDB)                  # Display available categories

# Focus on secreted signaling only
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling")
cellchatC@DB <- CellChatDB.use

# ----------------------------
# Step 2: Preprocessing
# ----------------------------

cellchatC <- subsetData(cellchatC, features = NULL)            # Use all genes
cellchatC <- identifyOverExpressedGenes(cellchatC)             # Find over-expressed genes
cellchatC <- identifyOverExpressedInteractions(cellchatC)      # Find over-expressed interactions
cellchatC <- projectData(cellchatC, PPI.mouse)                 # Project onto PPI (mouse default)

# ----------------------------
# Step 3: Compute Communication Probabilities
# ----------------------------

cellchatC <- computeCommunProb(cellchatC, raw.use = TRUE)      # Raw probability matrix
cellchatC <- filterCommunication(cellchatC, min.cells = 10)    # Remove small clusters
cellchatC <- computeCommunProbPathway(cellchatC)               # Per-pathway probabilities
cellchatC <- aggregateNet(cellchatC)                           # Aggregate network

# ----------------------------
# Step 4: Visualization
# ----------------------------

# Group size (used for plotting node size)
groupSize <- as.numeric(table(cellchatC@idents))

# Circle plot: number of interactions between cell types
netVisual_circle(
  cellchatC@net$count, 
  vertex.weight = groupSize, 
  weight.scale = TRUE,
  label.edge = FALSE, 
  title.name = "Number of Interactions"
)

# ----------------------------
# Step 5: Centrality Analysis
# ----------------------------

cellchatC <- netAnalysis_computeCentrality(cellchatC, slot.name = "netP")

# ----------------------------
# Step 6: Bubble Plot for Specific Interactions
# ----------------------------

bubble_plot <- netVisual_bubble(
  cellchatC, 
  sources.use = c("CD4T-ISG15"),
  remove.isolate = TRUE,
  angle.x = 45
)
