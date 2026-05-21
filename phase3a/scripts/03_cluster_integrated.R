#!/usr/bin/env Rscript

# =============================================================================
# Step 03: Cluster integrated embedding
# Loads scVI embedding and metadata CSVs from Python, rebuilds a merged
# Seurat object across all four samples, inserts the scVI latent space as a
# DimReduc, clusters, runs UMAP, and saves three diagnostic plots.
#
# Usage:
#   Rscript 03_cluster_integrated.R <embedding_csv> <metadata_csv>
#       <rds_wt_e85> <rds_wt_e95> <rds_ko_e85> <rds_ko_e95>
#       <output_rds> <umap_cluster_pdf> <umap_sample_pdf>
#       <umap_genotype_pdf> <composition_csv>
#
# Arguments:
#   embedding_csv     — 02_scvi_embedding.csv (cell_barcode + scVI_1..30)
#   metadata_csv      — 02_metadata.csv (cell_barcode + sample_id, genotype,
#                       timepoint)
#   rds_wt_e85/e95,
#   rds_ko_e85/e95    — 01_*_prepped.rds files (Phase 3a step 01 outputs)
#   output_rds        — clustered merged Seurat object
#   umap_cluster_pdf  — UMAP coloured by seurat_clusters
#   umap_sample_pdf   — UMAP coloured by sample_id
#   umap_genotype_pdf — UMAP coloured by genotype (primary QC plot)
#   composition_csv   — cluster × genotype × timepoint cell counts
#
# Environment variables (set in Nextflow script block from params):
#   CLUSTER_RES  — Louvain resolution  (default: 0.4)
#   SCVI_DIMS    — scVI latent dims    (default: 30)
# =============================================================================

library(Seurat)
library(tidyverse)
library(glue)

SEED <- 42

# =============================================================================
# 1. Arguments & Validation
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 11) {
    stop(glue(
        "Usage: Rscript 03_cluster_integrated.R <embedding_csv> <metadata_csv> ",
        "<rds_wt_e85> <rds_wt_e95> <rds_ko_e85> <rds_ko_e95> ",
        "<output_rds> <umap_cluster_pdf> <umap_sample_pdf> ",
        "<umap_genotype_pdf> <composition_csv>\n",
        "Got {length(args)} argument(s)."
    ))
}

embedding_csv     <- args[1]
metadata_csv      <- args[2]
rds_wt_e85        <- args[3]
rds_wt_e95        <- args[4]
rds_ko_e85        <- args[5]
rds_ko_e95        <- args[6]
output_rds        <- args[7]
umap_cluster_pdf  <- args[8]
umap_sample_pdf   <- args[9]
umap_genotype_pdf <- args[10]
composition_csv   <- args[11]

for (f in c(embedding_csv, metadata_csv, rds_wt_e85, rds_wt_e95, rds_ko_e85, rds_ko_e95)) {
    if (!file.exists(f)) stop(glue("Input file not found: {f}"))
}

cluster_res <- as.numeric(Sys.getenv("CLUSTER_RES", unset = "0.4"))
scvi_dims   <- as.integer(Sys.getenv("SCVI_DIMS",   unset = "30"))

message(glue("
========== Cluster Integrated Embedding (Step 03) ==========
  Embedding CSV  : {embedding_csv}
  Metadata CSV   : {metadata_csv}
  cluster_res    : {cluster_res}
  scvi_dims      : {scvi_dims}
  Seed           : {SEED}
  Output RDS     : {output_rds}
=============================================================
"))

# =============================================================================
# 2. Load Embedding & Metadata CSVs
# =============================================================================

message("Loading embedding and metadata CSVs...")
emb_df  <- read_csv(embedding_csv, show_col_types = FALSE)
meta_df <- read_csv(metadata_csv,  show_col_types = FALSE)

if (nrow(emb_df) != nrow(meta_df)) {
    stop(glue(
        "Assertion failed: embedding has {nrow(emb_df)} rows ",
        "but metadata has {nrow(meta_df)} rows."
    ))
}
message(glue("CSVs loaded: {nrow(emb_df)} cells, {ncol(emb_df) - 1} latent dims."))

# =============================================================================
# 3. Load & Merge Seurat Objects
# =============================================================================

sample_ids <- c("wt_e85", "wt_e95", "ko_e85", "ko_e95")
rds_paths  <- c(rds_wt_e85, rds_wt_e95, rds_ko_e85, rds_ko_e95)

message("Loading per-sample RDS files...")
objects <- map(rds_paths, \(p) {
    message(glue("  {p}"))
    readRDS(p)
})

message("Merging Seurat objects...")
merged <- merge(
    objects[[1]],
    y            = objects[2:4],
    add.cell.ids = sample_ids
)
message(glue("Merged: {ncol(merged)} cells x {nrow(merged)} genes"))

required_cols <- c("sample_id", "genotype", "timepoint")
missing_cols  <- required_cols[!required_cols %in% colnames(merged@meta.data)]
if (length(missing_cols) > 0) {
    stop(glue(
        "Assertion failed — missing meta.data columns after merge: ",
        "{paste(missing_cols, collapse = ', ')}"
    ))
}
message("Metadata assertion passed: sample_id, genotype, timepoint present.")

# =============================================================================
# 4. Hard Assertion — Barcode Match
# =============================================================================

seurat_barcodes    <- colnames(merged)
embedding_barcodes <- emb_df$cell_barcode

if (!setequal(seurat_barcodes, embedding_barcodes)) {
    n_only_seurat    <- sum(!seurat_barcodes    %in% embedding_barcodes)
    n_only_embedding <- sum(!embedding_barcodes %in% seurat_barcodes)
    stop(glue(
        "Assertion failed — barcode mismatch:\n",
        "  {n_only_seurat} barcode(s) in Seurat not in embedding CSV\n",
        "  {n_only_embedding} barcode(s) in embedding CSV not in Seurat"
    ))
}
message(glue("Barcode assertion passed: {length(seurat_barcodes)} cells matched."))

# =============================================================================
# 5. Insert scVI Embedding as DimReduc
# =============================================================================

message("Inserting scVI embedding as DimReduc...")
emb_matrix <- emb_df |>
    column_to_rownames("cell_barcode") |>
    as.matrix()

emb_matrix <- emb_matrix[seurat_barcodes, ]

if (!all(rownames(emb_matrix) == seurat_barcodes)) {
    stop("Assertion failed: cell order mismatch after reordering embedding matrix.")
}
message("Cell order assertion passed.")

merged[["scvi"]] <- CreateDimReducObject(
    embeddings = emb_matrix,
    key        = "scVI_",
    assay      = "RNA"
)

# =============================================================================
# 6. FindNeighbors
# =============================================================================

message(glue("Running FindNeighbors (reduction = scvi, dims = 1:{scvi_dims})..."))
merged <- FindNeighbors(
    merged,
    reduction = "scvi",
    dims      = 1:scvi_dims,
    verbose   = FALSE
)

# =============================================================================
# 7. FindClusters
# =============================================================================

message(glue("Running FindClusters (resolution = {cluster_res})..."))
merged <- FindClusters(
    merged,
    resolution  = cluster_res,
    random.seed = SEED,
    verbose     = FALSE
)

n_clusters <- length(levels(merged$seurat_clusters))
message(glue("Found {n_clusters} clusters at resolution {cluster_res}."))

cluster_summary <- merged@meta.data |>
    count(seurat_clusters, name = "n_cells") |>
    arrange(seurat_clusters)
print(cluster_summary)

# =============================================================================
# 8. RunUMAP
# =============================================================================

message(glue("Running UMAP (reduction = scvi, dims = 1:{scvi_dims})..."))
merged <- RunUMAP(
    merged,
    reduction = "scvi",
    dims      = 1:scvi_dims,
    seed.use  = SEED,
    verbose   = FALSE
)

# =============================================================================
# 9. UMAP Plots
# =============================================================================

message("Saving UMAP plots...")

p_clusters <- DimPlot(
    merged, reduction = "umap", group.by = "seurat_clusters",
    label = TRUE, repel = TRUE
) +
    ggtitle("scVI Integration — Clusters") +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))

p_sample <- DimPlot(
    merged, reduction = "umap", group.by = "sample_id"
) +
    ggtitle("scVI Integration — Sample") +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))

p_genotype <- DimPlot(
    merged, reduction = "umap", group.by = "genotype"
) +
    ggtitle("scVI Integration — Genotype (WT vs KO)") +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))

ggsave(umap_cluster_pdf,  plot = p_clusters, width = 8, height = 7)
ggsave(umap_sample_pdf,   plot = p_sample,   width = 8, height = 7)
ggsave(umap_genotype_pdf, plot = p_genotype, width = 8, height = 7)
message(glue("UMAP plots saved: {umap_cluster_pdf}, {umap_sample_pdf}, {umap_genotype_pdf}"))

# =============================================================================
# 10. Composition Table
# =============================================================================

message("Computing cluster × genotype × timepoint composition...")
composition <- merged@meta.data |>
    count(seurat_clusters, genotype, timepoint, name = "n_cells") |>
    arrange(seurat_clusters, genotype, timepoint)

message("Composition table:")
print(composition)

write_csv(composition, composition_csv)
message(glue("Composition saved to: {composition_csv}"))

# =============================================================================
# 11. Save RDS
# =============================================================================

saveRDS(merged, file = output_rds)
message(glue("Clustered object saved to: {output_rds}"))
message("Step 03 complete.")
