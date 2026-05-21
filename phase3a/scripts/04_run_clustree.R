#!/usr/bin/env Rscript

# =============================================================================
# Step 04: Clustering resolution stability (Clustree)
# Migrated from phase1/scripts/07_run_clustree.R; adapted for Phase 3a
# integrated object. Input is the merged Seurat object produced by step 03,
# which uses the RNA assay and scVI-based neighbour graph.
#
# Usage:
#   Rscript 04_run_clustree.R <input_rds> <output_plot>
#
# Arguments:
#   input_rds    — path to 03_integrated_clustered.rds (Phase 3a step 03 output)
#   output_plot  — path for clustree PDF
#
# Outputs:
#   <output_plot> — clustree PDF showing cluster splitting/merging across
#                   resolutions 0.2, 0.4, 0.6, 0.8, 1.0, 1.2 on the scVI
#                   neighbour graph built in step 03
# =============================================================================

library(Seurat)
library(clustree)
library(tidyverse)
library(glue)

# =============================================================================
# 1. Arguments
# =============================================================================

args        <- commandArgs(trailingOnly = TRUE)
input_rds   <- args[1]
output_plot <- args[2]

SEED        <- 42
resolutions <- seq(0.2, 1.2, by = 0.2)

message(glue("
========== Resolution Stability (Step 04) ==========
  Input       : {input_rds}
  Output      : {output_plot}
  Resolutions : {paste(resolutions, collapse = ', ')}
  Seed        : {SEED}
=====================================================
"))

# =============================================================================
# 2. Load Object
# =============================================================================

seurat_obj <- readRDS(input_rds)
message(glue("Loaded: {ncol(seurat_obj)} cells x {nrow(seurat_obj)} genes"))

# =============================================================================
# 3. Cluster at Multiple Resolutions
# =============================================================================

message("Clustering across resolutions...")
seurat_obj <- FindClusters(
    object      = seurat_obj,
    resolution  = resolutions,
    random.seed = SEED,
    verbose     = FALSE
)

# =============================================================================
# 4. Summarize Cluster Counts per Resolution
# =============================================================================

cluster_cols <- glue("RNA_snn_res.{resolutions}")

cluster_summary <- tibble(
    resolution = resolutions,
    n_clusters = map_int(cluster_cols, ~ length(unique(seurat_obj@meta.data[[.x]])))
)

message("Cluster count per resolution:")
print(cluster_summary)

# =============================================================================
# 5. Generate Clustree
# =============================================================================

message("Building cluster stability tree...")
tree_plot <- clustree(seurat_obj, prefix = "RNA_snn_res.") +
    ggtitle("Clustering Resolution Stability Tree") +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))

# =============================================================================
# 6. Save
# =============================================================================

ggsave(output_plot, plot = tree_plot, width = 10, height = 12)
message(glue("Clustree saved to: {output_plot}"))
message("Step 04 complete.")
