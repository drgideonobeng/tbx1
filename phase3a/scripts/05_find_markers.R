#!/usr/bin/env Rscript

# =============================================================================
# Step 05: Find cluster markers
# Runs FindAllMarkers on the integrated, clustered Seurat object to identify
# positive markers for each cluster. Saves the full table and a top-5-per-
# cluster summary for use in Phase 3b annotation.
#
# Usage:
#   Rscript 05_find_markers.R <input_rds> <all_markers_csv> <top5_csv>
#
# Arguments:
#   input_rds       — 03_integrated_clustered.rds (Phase 3a step 03 output)
#   all_markers_csv — full FindAllMarkers output (05_all_markers.csv)
#   top5_csv        — top 5 markers per cluster by avg_log2FC
#                     (05_top5_markers_per_cluster.csv)
# =============================================================================

library(Seurat)
library(tidyverse)
library(glue)

SEED <- 42
set.seed(SEED)

# =============================================================================
# 1. Arguments & Validation
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 3) {
    stop(glue(
        "Usage: Rscript 05_find_markers.R <input_rds> <all_markers_csv> <top5_csv>\n",
        "Got {length(args)} argument(s)."
    ))
}

input_rds       <- args[1]
all_markers_csv <- args[2]
top5_csv        <- args[3]

if (!file.exists(input_rds)) stop(glue("Input file not found: {input_rds}"))

message(glue("
========== Find Cluster Markers (Step 05) ==========
  Input RDS       : {input_rds}
  All markers CSV : {all_markers_csv}
  Top-5 CSV       : {top5_csv}
  Seed            : {SEED}
=====================================================
"))

# =============================================================================
# 2. Load Object
# =============================================================================

message("Loading clustered Seurat object...")
seu <- readRDS(input_rds)
n_clusters <- length(levels(seu$seurat_clusters))
message(glue("Loaded: {ncol(seu)} cells, {n_clusters} clusters."))

# Seurat v5: layers are split after merge — join before DE testing
seu <- JoinLayers(seu)
message("Layers joined.")

# =============================================================================
# 3. FindAllMarkers
# =============================================================================

message("Running FindAllMarkers (only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)...")
markers <- FindAllMarkers(
    seu,
    only.pos        = TRUE,
    min.pct         = 0.25,
    logfc.threshold = 0.25,
    verbose         = FALSE
)
message(glue("Found {nrow(markers)} marker genes across {n_clusters} clusters."))

# =============================================================================
# 4. Save Full Marker Table
# =============================================================================

write_csv(markers, all_markers_csv)
message(glue("Full marker table saved to: {all_markers_csv}"))

# =============================================================================
# 5. Top 5 per Cluster
# =============================================================================

top5 <- markers |>
    group_by(cluster) |>
    slice_max(order_by = avg_log2FC, n = 5) |>
    ungroup() |>
    arrange(cluster, desc(avg_log2FC))

message("Top 5 markers per cluster:")
print(top5)

write_csv(top5, top5_csv)
message(glue("Top-5 marker table saved to: {top5_csv}"))
message("Step 05 complete.")
