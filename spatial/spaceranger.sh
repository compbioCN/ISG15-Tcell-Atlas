#!/bin/bash

# ================================
# Batch Spaceranger Count Script
# ================================
# This script runs `spaceranger count` for multiple samples,
# assuming the following folder structure:
#
# project_root/
# ├── raw_fastq/               # contains one subfolder per sample
# ├── TIF/                     # contains one .tif per sample
# ├── TIF_CYTA/                # contains one cyta .tif per sample (optional)
# ├── JSON/                    # alignment JSON file(s)
# ├── spaceranger_resources/   # probe-set and related files
# ├── ref/                     # transcriptome reference
# └── run_spaceranger_batch.sh # this script
# ================================

# -------- Configuration --------
SPACERANGER="spaceranger"  # Make sure spaceranger is in your PATH
TRANSCRIPTOME="ref/refdata-gex-GRCh38-2020-A"
PROBESET="spaceranger_resources/Visium_Human_Transcriptome_Probe_Set_v2.0_GRCh38-2020-A.csv"
CYTA_FOLDER="TIF_CYTA"
IMAGE_FOLDER="TIF"
JSON_FILE="JSON/alignment.json"     # Set to one fixed file or make dynamic per sample
SLIDE="SLIDE_ID"                    # Replace with your slide ID or customize per sample
AREA="A1"
FASTQ_BASE="raw_fastq"
OUT_BASE="spaceranger_output"

# -------- Batch Execution --------
mkdir -p "$OUT_BASE"

for SAMPLE_DIR in "$FASTQ_BASE"/*; do
  SAMPLE_NAME=$(basename "$SAMPLE_DIR")
  FASTQ_PATH="$SAMPLE_DIR"
  IMAGE_PATH="${IMAGE_FOLDER}/${SAMPLE_NAME}.tif"
  CYTA_PATH="${CYTA_FOLDER}/CYTA_${SAMPLE_NAME}.tif"
  OUTPUT_DIR="${OUT_BASE}/${SAMPLE_NAME}_outs"

  echo ">>> Processing sample: $SAMPLE_NAME"

  "$SPACERANGER" count \
    --id="${SAMPLE_NAME}_outs" \
    --transcriptome="$TRANSCRIPTOME" \
    --fastqs="$FASTQ_PATH" \
    --probe-set="$PROBESET" \
    --cytaimage="$CYTA_PATH" \
    --image="$IMAGE_PATH" \
    --loupe-alignment="$JSON_FILE" \
    --slide="$SLIDE" \
    --area="$AREA" \
    --create-bam=false \
    --jobmode=local

  echo "✓ Completed: $SAMPLE_NAME"
  echo "--------------------------"
done
