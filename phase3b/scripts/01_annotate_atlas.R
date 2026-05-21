#!/usr/bin/env Rscript

# =============================================================================
# Step 01: Annotate integrated atlas
# Applies human-provided cluster labels to the integrated Seurat object,
# validates the label mapping, produces three diagnostic UMAPs, and saves
# the annotated object for downstream DE analysis.
#
# Usage:
#   Rscript 01_annotate_atlas.R <clustered_rds> <cluster_labels_csv>
#       <output_rds> <umap_annotated_pdf> <umap_by_genotype_pdf>
#       <umap_by_timepoint_pdf>
#
# Arguments:
#   clustered_rds         — 03_integrated_clustered.rds (Phase 3a output)
#   cluster_labels_csv    — two-column CSV: cluster_id, label
#   output_rds            — annotated Seurat object with cell_type in meta.data
#   umap_annotated_pdf    — UMAP coloured by cell type (labelled)
#   umap_by_genotype_pdf  — UMAP split by genotype (WT | KO), coloured by cell type
#   umap_by_timepoint_pdf — UMAP split by timepoint (E85 | E95), coloured by cell type
# =============================================================================

library(Seurat)
library(tidyverse)
library(glue)

SEED <- 42

# =============================================================================
# 1. Arguments & Validation
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 6) {
    stop(glue(
        "Usage: Rscript 01_annotate_atlas.R <clustered_rds> <cluster_labels_csv> ",
        "<output_rds> <umap_annotated_pdf> <umap_by_genotype_pdf> ",
        "<umap_by_timepoint_pdf>\n",
        "Got {length(args)} argument(s)."
    ))
}

clustered_rds         <- args[1]
cluster_labels_csv    <- args[2]
output_rds            <- args[3]
umap_annotated_pdf    <- args[4]
umap_by_genotype_pdf  <- args[5]
umap_by_timepoint_pdf <- args[6]

for (f in c(clustered_rds, cluster_labels_csv)) {
    if (!file.exists(f)) stop(glue("Input file not found: {f}"))
}

message(glue("
========== Annotate Integrated Atlas (Step 01) ==========
  Clustered RDS         : {clustered_rds}
  Cluster labels CSV    : {cluster_labels_csv}
  Output RDS            : {output_rds}
  UMAP annotated        : {umap_annotated_pdf}
  UMAP by genotype      : {umap_by_genotype_pdf}
  UMAP by timepoint     : {umap_by_timepoint_pdf}
  Seed                  : {SEED}
=========================================================
"))

# =============================================================================
# 2. Load Object
# =============================================================================

message("Loading clustered Seurat object...")
seu <- readRDS(clustered_rds)
n_clusters <- length(levels(seu$seurat_clusters))
message(glue("Loaded: {ncol(seu)} cells, {n_clusters} clusters."))

# =============================================================================
# 3. Load & Validate Cluster Labels
# =============================================================================

message("Loading cluster labels...")
labels_df <- read_csv(cluster_labels_csv, show_col_types = FALSE)

required_cols <- c("cluster_id", "label")
missing_cols  <- required_cols[!required_cols %in% colnames(labels_df)]
if (length(missing_cols) > 0) {
    stop(glue(
        "Assertion failed — cluster_labels_csv missing columns: ",
        "{paste(missing_cols, collapse = ', ')}"
    ))
}

bad_labels <- labels_df |> filter(is.na(label) | str_trim(label) == "")
if (nrow(bad_labels) > 0) {
    stop(glue(
        "Assertion failed — empty or NA label for cluster(s): ",
        "{paste(bad_labels$cluster_id, collapse = ', ')}"
    ))
}

object_clusters <- levels(seu$seurat_clusters)
csv_clusters    <- as.character(labels_df$cluster_id)
missing_in_csv  <- object_clusters[!object_clusters %in% csv_clusters]
if (length(missing_in_csv) > 0) {
    stop(glue(
        "Assertion failed — cluster(s) in object missing from CSV: ",
        "{paste(missing_in_csv, collapse = ', ')}"
    ))
}

labels_df <- labels_df |>
    mutate(cluster_id = as.character(cluster_id)) |>
    arrange(match(cluster_id, object_clusters))

label_map <- setNames(labels_df$label, labels_df$cluster_id)

message("Label mapping (cluster → cell type):")
walk2(names(label_map), label_map, \(k, v) message(glue("  {k} → {v}")))

# =============================================================================
# 4. Apply Labels
# =============================================================================

message("Applying labels via RenameIdents...")
Idents(seu) <- "seurat_clusters"
seu <- RenameIdents(seu, label_map)
seu$cell_type <- Idents(seu)
message(glue("Labels applied. {length(unique(seu$cell_type))} unique cell types."))

# =============================================================================
# 5. UMAP Plots
# =============================================================================

message("Saving UMAP plots...")

p_annotated <- DimPlot(
    seu, reduction = "umap", group.by = "cell_type",
    label = TRUE, repel = TRUE
) +
    ggtitle("Integrated Atlas — Cell Types") +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))

p_genotype <- DimPlot(
    seu, reduction = "umap", group.by = "cell_type",
    split.by = "genotype", label = TRUE, repel = TRUE
) +
    ggtitle("Cell Types by Genotype (WT vs KO)") +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))

p_timepoint <- DimPlot(
    seu, reduction = "umap", group.by = "cell_type",
    split.by = "timepoint", label = TRUE, repel = TRUE
) +
    ggtitle("Cell Types by Timepoint (E8.5 vs E9.5)") +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))

ggsave(umap_annotated_pdf,    plot = p_annotated, width = 9,  height = 7)
ggsave(umap_by_genotype_pdf,  plot = p_genotype,  width = 14, height = 7)
ggsave(umap_by_timepoint_pdf, plot = p_timepoint, width = 14, height = 7)
message(glue(
    "UMAP plots saved: {umap_annotated_pdf}, ",
    "{umap_by_genotype_pdf}, {umap_by_timepoint_pdf}"
))

# =============================================================================
# 6. Composition Table
# =============================================================================

message("Cluster × genotype × timepoint composition (with cell type labels):")
composition <- seu@meta.data |>
    count(cell_type, seurat_clusters, genotype, timepoint, name = "n_cells") |>
    arrange(cell_type, seurat_clusters, genotype, timepoint)
print(composition)

# =============================================================================
# 7. Save Annotated Object
# =============================================================================

saveRDS(seu, file = output_rds)
message(glue("Annotated object saved to: {output_rds}"))
message("Step 01 complete.")
