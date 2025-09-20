# Spatiotemporal Dynamics and Immunotherapy Relevance of Interferon-Responsive ISG15⁺ T Cells in Cancer

This repository hosts the analysis code, data processing workflows, and supporting scripts for the manuscript:

> **Spatiotemporal Dynamics and Immunotherapy Relevance of Interferon-Responsive ISG15⁺ T Cells in Cancer**

---

## 📖 Overview
Immunotherapy has transformed the treatment landscape across solid tumors, yet the heterogeneity of T cell states remains a barrier to efficacy.  
In this study, we constructed a **pan-cancer single-cell atlas of >2.3 million T cells from 41 cancer types**, identifying a previously underappreciated **interferon-responsive ISG15⁺ T cell subset**.  

Key findings:
- ISG15⁺ T cells are conserved across cancers, progressively enriched in HNSCC, and localize to tumor cores.
- ISG15⁺ CD8⁺ T cells display **transcriptional plasticity**, co-expressing cytotoxic and exhaustion programs.
- Multi-omics and spatial analyses revealed **distinct chromatin accessibility** and **epithelial interaction patterns (LGALS9–CD44 axis)**.
- ISG15⁺ T cells decline upon ICI treatment in both **human cohorts** and **murine models**, with pre-treatment abundance correlating with response.
- Clonal tracing indicates their potential to transition toward effector states, positioning them as a **reversible transitional state** linked to immunotherapy outcomes.

---

## 📂 Repository Contents
- `scRNA-seq/` — Scripts for preprocessing, QC, clustering, and integration of scRNA-seq datasets  
- `spatial/` — Code for Visium HD and Xenium spatial transcriptomics analyses  
- `scATAC-seq/` — Workflows for scATAC-seq preprocessing and regulatory analysis using ArchR  
- `bulk/` — Bulk RNA-seq GSVA scoring and survival analysis scripts 
---

## ⚙️ Requirements
- R (≥ 4.1) with Seurat, Harmony, survminer, GSVA, CellChat, MAESTRO  
- Python (≥ 3.9) with Scanpy, BBKNN, Squidpy, scFEA, PySCENIC  
- 10x Genomics Cell Ranger / Space Ranger (for raw data processing)  
- ArchR (≥ 1.0.2) for scATAC-seq analysis  

Installation of required packages can be managed via `conda` or `renv` (R).  
A detailed `environment.yml` and `requirements.txt` are provided.

---

## 🚀 Usage
1. Clone the repository:
   ```bash
   git clone https://github.com/compbioCN/ISG15-Tcell-Atlas.git
   cd ISG15-Tcell-Atlas
