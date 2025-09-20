# Required Libraries
library(Seurat)
library(scRepertoire)
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggalluvial)

#-------------------------#
# Step 1: Load scRNA Data #
#-------------------------#
# Replace this path with your actual Seurat object
seurat_obj <- readRDS("GSE200996_scRNA.rds")

# OPTIONAL: Cell type annotations if not already present
# seurat_obj$celltype <- Idents(seurat_obj) or other method

#------------------------------------------#
# Step 2: Load and Process TCR Annotation Files
#------------------------------------------#
# Load all filtered_contig_annotations.csv.gz in a folder

# Path to directory containing TCR files
vdj_dir <- "/home/guile/CvI/IMTHE/GSE200996 scTCR"
tcr_files <- list.files(vdj_dir, pattern = "filtered_contig_annotations.*.csv.gz", full.names = TRUE)

# Auto-generate sample names from file names
sample_ids <- gsub("_filtered_contig_annotations.*", "", basename(tcr_files))

# Read TCR files
TCR_list <- lapply(tcr_files, read.csv)
names(TCR_list) <- sample_ids

# Combine TCR data
combined_tcr <- combineTCR(
  TCR_list,
  samples = sample_ids,
  ID = ifelse(grepl("pre", sample_ids, ignore.case = TRUE), "pre", "post"),
  cells = "T-AB"
)

#-----------------------------------------#
# Step 3: Format TCR barcodes for matching
#-----------------------------------------#
# Strip barcode prefix then re-add consistent format
for (i in seq_along(combined_tcr)) {
  combined_tcr[[i]] <- stripBarcode(combined_tcr[[i]], column = 1, connector = "_", num_connects = 3)
  combined_tcr[[i]]$barcode <- paste0(sample_ids[i], "_", combined_tcr[[i]]$barcode)
}

#--------------------------------------------------#
# Step 4: Integrate TCR data into Seurat object
#--------------------------------------------------#
seurat_tcr <- combineExpression(
  combined_tcr,
  seurat_obj,
  cloneCall = "aa",
  proportion = FALSE,
  cloneTypes = c(Single=1, Small=3, Medium=10, Large=30, Hyperexpanded=100)
)

#--------------------------------------------------#
# Step 5: Extract CD8T-ISG15 clones from Pre-treatment
#--------------------------------------------------#
cd8_pre_clones <- seurat_tcr %>%
  subset(orig.ident %in% sample_ids[grepl("pre", sample_ids, ignore.case = TRUE)] & celltype == "CD8T-ISG15" & !is.na(CTaa)) %>%
  pull(CTaa) %>%
  unique()

#--------------------------------------------------#
# Step 6: Trace clone fate in Post-treatment samples
#--------------------------------------------------#
clone_flow_df <- seurat_tcr %>%
  subset(orig.ident %in% sample_ids & CTaa %in% cd8_pre_clones) %>%
  mutate(group = ifelse(grepl("pre", orig.ident, ignore.case = TRUE), "Pre", "Post")) %>%
  select(CTaa, group, celltype)

# Pivot for alluvial input
alluvial_input <- clone_flow_df %>%
  group_by(CTaa, group, celltype) %>%
  summarise(n = n(), .groups = "drop") %>%
  pivot_wider(names_from = group, values_from = celltype) %>%
  filter(!is.na(Pre) & !is.na(Post))

# Summarize flows
plot_data <- alluvial_input %>%
  group_by(CTaa, Pre, Post) %>%
  summarise(n = n(), .groups = "drop")

#--------------------------------------------------#
# Step 7: Draw alluvial plot
#--------------------------------------------------#
ggplot(plot_data, aes(axis1 = Pre, axis2 = Post, y = n)) +
  geom_alluvium(aes(fill = Post), width = 1/12) +
  geom_stratum(width = 1/12, fill = "white", color = "black") +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 3.5, fontface = "bold") +
  theme_minimal() +
  theme(
    axis.title = element_blank(),
    panel.grid = element_blank(),
    legend.position = "none"
  ) +
  scale_fill_manual(values = c(
    "CD8Tem" = "#C2A5CF",
    "CD8Trm" = "#92C5DE",
    "Tprolif" = "#F4A582"
  )) +
  ggtitle("Clonal Flow: CD8T-ISG15 (Pre) to CD8 States (Post)")
