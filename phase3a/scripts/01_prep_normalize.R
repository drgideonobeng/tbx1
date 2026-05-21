#!/usr/bin/env Rscript

# =============================================================================
# Step 01: Per-sample prep & normalization
# Adds sample metadata, log-normalizes, and selects highly variable genes.
# Prepares each Phase 1 filtered object for cross-sample scVI integration;
# raw counts are preserved in RNA@counts for VAE input downstream.
#
# Usage:
#   Rscript 01_prep_normalize.R <filtered_rds> <sample_id> <genotype>
#                               <timepoint> <n_hvgs> <output_rds>
#
# Arguments:
#   filtered_rds  — path to 04_seurat_filtered.rds (Phase 1 output)
#   sample_id     — string identifier, e.g. "wt_e85"
#   genotype      — "WT" or "KO"
#   timepoint     — "E85" or "E95"
#   n_hvgs        — number of highly variable genes to select
#   output_rds    — path for prepped Seurat object
#
# Outputs:
#   <output_rds>  — Seurat object with:
#                   * meta.data columns: sample_id, genotype, timepoint
#                   * log-normalised data (NormalizeData; raw counts intact)
#                   * n_hvgs HVGs identified with VST (FindVariableFeatures)
# =============================================================================

library(Seurat)
library(tidyverse)
library(glue)

# =============================================================================
# 1. Arguments & Validation
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 6) {
    stop(glue(
        "Usage: Rscript 01_prep_normalize.R ",
        "<filtered_rds> <sample_id> <genotype> <timepoint> <n_hvgs> <output_rds>\n",
        "Got {length(args)} argument(s)."
    ))
}

filtered_rds <- args[1]
sample_id    <- args[2]
genotype     <- args[3]
timepoint    <- args[4]
n_hvgs       <- as.integer(args[5])
output_rds   <- args[6]

if (!file.exists(filtered_rds)) stop(glue("Input file not found: {filtered_rds}"))
if (is.na(n_hvgs) || n_hvgs < 1) stop(glue("n_hvgs must be a positive integer; got: {args[5]}"))

message(glue("
========== Per-sample Prep & Normalization (Step 01) ==========
  Input      : {filtered_rds}
  Output     : {output_rds}
  sample_id  : {sample_id}
  genotype   : {genotype}
  timepoint  : {timepoint}
  n_hvgs     : {n_hvgs}
===============================================================
"))

# =============================================================================
# 2. Load Object
# =============================================================================

seurat_obj <- readRDS(filtered_rds)
message(glue("Loaded: {ncol(seurat_obj)} cells x {nrow(seurat_obj)} genes"))

# =============================================================================
# 3. Add Sample Metadata
# =============================================================================

seurat_obj$sample_id <- sample_id
seurat_obj$genotype  <- genotype
seurat_obj$timepoint <- timepoint

required_cols <- c("sample_id", "genotype", "timepoint")
missing_cols  <- required_cols[!required_cols %in% colnames(seurat_obj@meta.data)]

if (length(missing_cols) > 0) {
    stop(glue("Assertion failed — missing meta.data columns: {paste(missing_cols, collapse = ', ')}"))
}

message("Metadata assertion passed: sample_id, genotype, timepoint present.")

# =============================================================================
# 4. Log-normalize
# =============================================================================

message("Running LogNormalize (scale.factor = 10,000)...")
seurat_obj <- NormalizeData(
    seurat_obj,
    normalization.method = "LogNormalize",
    scale.factor         = 10000,
    verbose              = FALSE
)

# =============================================================================
# 5. Identify Highly Variable Genes
# =============================================================================

message(glue("Selecting top {n_hvgs} HVGs with VST..."))
seurat_obj <- FindVariableFeatures(
    seurat_obj,
    selection.method = "vst",
    nfeatures        = n_hvgs,
    verbose          = FALSE
)

top10 <- VariableFeatures(seurat_obj) |> head(10)
message(glue("Top 10 HVGs: {paste(top10, collapse = ', ')}"))

# =============================================================================
# 6. Summary & Save
# =============================================================================

n_hvgs_found <- length(VariableFeatures(seurat_obj))

tibble(
    field = c("sample_id", "genotype", "timepoint", "n_cells", "n_genes", "n_hvgs_selected"),
    value = c(sample_id, genotype, timepoint,
              as.character(ncol(seurat_obj)),
              as.character(nrow(seurat_obj)),
              as.character(n_hvgs_found))
) |> print()

saveRDS(seurat_obj, file = output_rds)
message(glue("Prepped object saved to: {output_rds}"))
message("Step 01 complete.")
