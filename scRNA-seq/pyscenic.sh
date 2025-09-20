# ------------------------------------------
# Step 1-3: pySCENIC workflow (bash commands)
# ------------------------------------------

# 1. Inference gene regulatory network (GRN)
nohup pyscenic grn \
  --num_workers 16 \
  --sparse \
  --method grnboost2 \
  --output grn.csv \
  sce.loom allTFs_hg38.txt &

# 2. Prune GRN using motif enrichment
nohup pyscenic ctx \
  --num_workers 16 \
  --output regulons.csv \
  --expression_mtx_fname sce.loom \
  --mode "custom_multiprocessing" \
  --annotations_fname motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl \
  grn.csv \
  hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather &

# 3. Score regulon activity per cell
nohup pyscenic aucell \
  --num_workers 16 \
  --output sample_SCENIC.loom \
  sce.loom \
  regulons.csv &
