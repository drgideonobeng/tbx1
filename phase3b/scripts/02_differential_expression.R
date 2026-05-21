#!/usr/bin/env Rscript

# =============================================================================
# Step 02: Differential expression — WT vs KO, pseudobulk per cell type × timepoint
# For each cell type × timepoint combination, aggregates raw counts per
# sample_id and runs DESeq2 (design = ~ genotype). Saves one CSV per
# comparison and a summary table across all combinations.
#
# Usage:
#   Rscript 02_differential_expression.R <annotated_rds> <de_results_dir>
#       <summary_csv>
#
# Arguments:
#   annotated_rds  — 01_annotated_atlas.rds (Phase 3b step 01 output)
#   de_results_dir — directory for per-comparison CSVs
#                    (DE_{cell_type}_{timepoint}_KOvsWT.csv)
#   summary_csv    — one row per combination: cell_type, timepoint, n_wt,
#                    n_ko, n_sig_up, n_sig_down, status (run/skipped)
#
# Environment variables (set in Nextflow script block from params):
#   MIN_CELLS      — minimum cells per genotype per timepoint (default: 10)
#   LFC_THRESHOLD  — log2 fold-change threshold for sig counting (default: 0.25)
#   PADJ_THRESHOLD — adjusted p-value threshold                 (default: 0.05)
# =============================================================================

library(Seurat)
library(DESeq2)
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
        "Usage: Rscript 02_differential_expression.R <annotated_rds> ",
        "<de_results_dir> <summary_csv>\n",
        "Got {length(args)} argument(s)."
    ))
}

annotated_rds  <- args[1]
de_results_dir <- args[2]
summary_csv    <- args[3]

if (!file.exists(annotated_rds)) stop(glue("Input file not found: {annotated_rds}"))

min_cells      <- as.integer(Sys.getenv("MIN_CELLS",       unset = "10"))
lfc_threshold  <- as.numeric(Sys.getenv("LFC_THRESHOLD",   unset = "0.25"))
padj_threshold <- as.numeric(Sys.getenv("PADJ_THRESHOLD",  unset = "0.05"))

message(glue("
========== Differential Expression (Step 02) ==========
  Annotated RDS    : {annotated_rds}
  DE results dir   : {de_results_dir}
  Summary CSV      : {summary_csv}
  min_cells        : {min_cells}
  lfc_threshold    : {lfc_threshold}
  padj_threshold   : {padj_threshold}
  Seed             : {SEED}
=======================================================
"))

# =============================================================================
# 2. Load Object & Create Output Directory
# =============================================================================

message("Loading annotated Seurat object...")
seu <- readRDS(annotated_rds)
seu <- JoinLayers(seu)
message(glue("Loaded: {ncol(seu)} cells, {length(unique(seu$cell_type))} cell types."))

if (!"counts" %in% Layers(seu[["RNA"]])) {
    stop("Assertion failed: raw counts layer not found in RNA assay. Cannot run pseudobulk DE.")
}
message("Raw counts assertion passed.")

dir.create(de_results_dir, showWarnings = FALSE, recursive = TRUE)
message(glue("Output directory ready: {de_results_dir}"))

# =============================================================================
# 3. DE Loop: cell type × timepoint
# =============================================================================

cell_types <- sort(unique(as.character(seu$cell_type)))

message(glue(
    "Running DE for {length(cell_types)} cell types ",
    "(pooled timepoints, design = ~ timepoint + genotype)."
))

summary_rows <- list()

for (ct in cell_types) {

    safe_name <- gsub("[/ ]", "_", ct)
    out_csv   <- file.path(de_results_dir, glue("DE_{safe_name}_KOvsWT.csv"))

    sub  <- seu[, seu$cell_type == ct]
    n_wt <- sum(sub$genotype == "WT")
    n_ko <- sum(sub$genotype == "KO")

    if (n_wt < min_cells || n_ko < min_cells) {
        message(glue(
            "  SKIP [{ct}] — WT: {n_wt}, KO: {n_ko} (min_cells = {min_cells})"
        ))
        summary_rows[[length(summary_rows) + 1]] <- tibble(
            cell_type  = ct,  n_wt = n_wt, n_ko = n_ko,
            n_sig_up   = NA_integer_, n_sig_down = NA_integer_,
            status     = "skipped"
        )
        next
    }

    message(glue("  RUN  [{ct}] — WT: {n_wt}, KO: {n_ko}"))

    # Pseudobulk: aggregate raw counts per sample_id (4 samples: wt_e85/e95, ko_e85/e95)
    agg       <- AggregateExpression(sub, group.by = "sample_id",
                                     assays = "RNA", slot = "counts")
    count_mat <- agg$RNA

    # colData: derive genotype and timepoint per pseudobulk sample.
    # AggregateExpression replaces underscores with dashes — normalise for lookup.
    lookup_keys   <- gsub("-", "_", colnames(count_mat))
    genotype_map  <- sub@meta.data |> distinct(sample_id, genotype)  |> deframe()
    timepoint_map <- sub@meta.data |> distinct(sample_id, timepoint) |> deframe()

    col_data <- data.frame(
        row.names = colnames(count_mat),
        genotype  = factor(genotype_map[lookup_keys],  levels = c("WT", "KO")),
        timepoint = factor(timepoint_map[lookup_keys])
    )

    n_params <- ncol(model.matrix(~ timepoint + genotype, col_data))
    if (ncol(count_mat) <= n_params) {
        message(glue(
            "  SKIP [{ct}] — {ncol(count_mat)} pseudobulk sample(s), ",
            "need > {n_params} to fit design ~ timepoint + genotype"
        ))
        summary_rows[[length(summary_rows) + 1]] <- tibble(
            cell_type  = ct,  n_wt = n_wt, n_ko = n_ko,
            n_sig_up   = NA_integer_, n_sig_down = NA_integer_,
            status     = "skipped"
        )
        next
    }

    # DESeq2: block on timepoint, test genotype effect
    set.seed(SEED)
    dds <- DESeqDataSetFromMatrix(
        countData = count_mat,
        colData   = col_data,
        design    = ~ timepoint + genotype
    )
    dds <- DESeq(dds, quiet = TRUE)

    res <- results(dds,
                   contrast = c("genotype", "KO", "WT"),
                   alpha    = padj_threshold)

    res_df <- as.data.frame(res) |>
        rownames_to_column("gene") |>
        arrange(padj)

    write_csv(res_df, out_csv)

    n_sig_up   <- sum(res_df$padj < padj_threshold &
                      res_df$log2FoldChange >  lfc_threshold, na.rm = TRUE)
    n_sig_down <- sum(res_df$padj < padj_threshold &
                      res_df$log2FoldChange < -lfc_threshold, na.rm = TRUE)

    message(glue(
        "    → {nrow(res_df)} genes tested | ",
        "sig up: {n_sig_up} | sig down: {n_sig_down} | saved: {basename(out_csv)}"
    ))

    summary_rows[[length(summary_rows) + 1]] <- tibble(
        cell_type  = ct,       n_wt       = n_wt,       n_ko       = n_ko,
        n_sig_up   = n_sig_up, n_sig_down = n_sig_down, status     = "run"
    )
}

# =============================================================================
# 4. Hard Assertion: at least one comparison ran
# =============================================================================

summary_tbl <- bind_rows(summary_rows)
n_run <- sum(summary_tbl$status == "run")

if (n_run == 0) {
    stop(glue(
        "Assertion failed: all {nrow(summary_tbl)} cell types were skipped ",
        "(min_cells = {min_cells}). No DE results produced."
    ))
}
message(glue("{n_run} cell type(s) ran, {nrow(summary_tbl) - n_run} skipped."))

# =============================================================================
# 5. Save & Print Summary
# =============================================================================

write_csv(summary_tbl, summary_csv)
message(glue("Summary saved to: {summary_csv}"))

message("DE summary:")
print(summary_tbl)

message("Step 02 complete.")
