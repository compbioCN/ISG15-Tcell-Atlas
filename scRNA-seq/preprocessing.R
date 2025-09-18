library(Seurat)        # v4.3.0 recommended
library(harmony)       # loaded if you later want batch correction
library(MAESTRO)       # v1.5.1 for FindAllMarkersMAESTRO
library(reticulate)    # to call Scrublet (Python)

# Root directory containing dataset folders
datasets_dir <- "/home/guile/CvI/data raw"

# Get all dataset folder paths 
dataset_dirs <- list.dirs(datasets_dir, recursive = FALSE)

# Helper: run Scrublet on a Seurat object (per sequencing library) and return filtered object
run_scrublet <- function(seu, expected_doublet_rate = 0.06) {
  # Scrublet expects a cells-by-genes raw count matrix (numpy array). Seurat stores genes x cells.
  if (!py_module_available("scrublet")) {
    warning("Python module 'scrublet' not found. Skipping doublet removal for this sample.")
    return(seu)
  }
  scrublet <- import("scrublet")
  np <- import("numpy")
  
  # Extract counts and transpose to cells x genes
  mat <- as.matrix(GetAssayData(seu, assay = "RNA", slot = "counts"))
  if (ncol(mat) == 0 || nrow(mat) == 0) return(seu)
  m_cells_by_genes <- t(mat)
  
  # Initialize Scrublet and call doublets
  scr <- scrublet$Scrublet(np$array(m_cells_by_genes, dtype = "float32"))
  res <- scr$call_doublets(expected_doublet_rate = expected_doublet_rate)
  
  # Predicted doublets is a boolean array aligned to cells
  pred <- as.logical(res[[2]])
  if (length(pred) != nrow(m_cells_by_genes)) return(seu)
  
  # Keep only cells predicted as singlets
  keep_cells <- colnames(seu)[!pred]
  subset(seu, cells = keep_cells)
}

# Container for processed Seurat objects per dataset
processed_datasets <- list()

# Iterate through each dataset folder
for (dataset_dir in dataset_dirs) {
  # Collect all sample subfolders (each is a sequencing library)
  sample_dirs <- list.dirs(dataset_dir, recursive = FALSE)
  
  # Read, annotate, QC, and doublet-filter each sample/library
  seurat_list <- lapply(sample_dirs, function(sample_dir) {
    sample_name <- basename(sample_dir)  # use subfolder name as patient/library id
    
    # Read10X matrices (Cell Ranger v3.1.0 output assumed; prefitered matrices acceptable)
    counts <- Read10X(data.dir = sample_dir)
    
    # Create Seurat object with minimal prefilter; apply study QC below
    seu <- CreateSeuratObject(
      counts = counts,
      min.cells = 0,   # do not prefilter by genes expressed across cells
      min.features = 0 # apply nFeature_RNA >= 500 later per spec
    )
    seu$patient <- sample_name
    
    # Compute mitochondrial percentage (human GRCh38: genes start with "MT-")
    seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")
    
    # Cell-level filtering: keep cells with >=500 detected genes and <=15% mitochondrial UMIs
    seu <- subset(seu, subset = nFeature_RNA >= 500 & percent.mt <= 15)
    
    # Doublet removal with Scrublet per library (expected_doublet_rate = 0.06)
    seu <- run_scrublet(seu, expected_doublet_rate = 0.06)
    
    return(seu)
  })
  
  # Merge all samples within the current dataset
  if (length(seurat_list) == 1) {
    combined_seurat <- seurat_list[[1]]
  } else {
    combined_seurat <- Reduce(function(x, y) merge(x, y), seurat_list)
  }
