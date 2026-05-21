#!/usr/bin/env Rscript

# =============================================================================
# Step 03: GO enrichment
# Runs clusterProfiler::enrichGO on significant DE genes (up and down
# separately) for each cell type comparison. Saves a combined results table
# and a dotplot faceted by direction (up/down).
#
# Usage:
#   Rscript 03_enrichment.R <de_results_dir> <de_summary_csv>
#       <enrichment_results_csv> <enrichment_dotplot_pdf>
#
# Arguments:
#   de_results_dir         — directory containing DE_*.csv files (step 02 output)
#   de_summary_csv         — 02_DE_summary.csv (step 02 output)
#   enrichment_results_csv — combined enrichGO results across all comparisons
#   enrichment_dotplot_pdf — dotplot: top 5 GO terms per cell type,
#                            faceted by direction (up/down)
#
# Environment variables (set in Nextflow script block from params):
#   GO_PADJ      — adjusted p-value cutoff for enrichGO (default: 0.05)
#   GO_MIN_GENES — minimum gene set size                (default: 5)
# =============================================================================

library(clusterProfiler)
library(org.Mm.eg.db)
library(tidyverse)
library(glue)

SEED       <- 42
LFC_FILTER <- 0.25   # fixed DE significance LFC threshold

# =============================================================================
# 1. Arguments & Validation
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 4) {
    stop(glue(
        "Usage: Rscript 03_enrichment.R <de_results_dir> <de_summary_csv> ",
        "<enrichment_results_csv> <enrichment_dotplot_pdf>\n",
        "Got {length(args)} argument(s)."
    ))
}

de_results_dir         <- args[1]
de_summary_csv         <- args[2]
enrichment_results_csv <- args[3]
enrichment_dotplot_pdf <- args[4]

for (f in c(de_results_dir, de_summary_csv)) {
    if (!file.exists(f)) stop(glue("Input not found: {f}"))
}

go_padj      <- as.numeric(Sys.getenv("GO_PADJ",      unset = "0.05"))
go_min_genes <- as.integer(Sys.getenv("GO_MIN_GENES",  unset = "5"))

message(glue("
========== GO Enrichment (Step 03) ==========
  DE results dir         : {de_results_dir}
  DE summary CSV         : {de_summary_csv}
  Enrichment results CSV : {enrichment_results_csv}
  Enrichment dotplot PDF : {enrichment_dotplot_pdf}
  go_padj                : {go_padj}
  go_min_genes           : {go_min_genes}
  LFC filter             : {LFC_FILTER}
=============================================
"))

# =============================================================================
# 2. Load DE Summary and Locate DE Files
# =============================================================================

message("Loading DE summary...")
summary_df <- read_csv(de_summary_csv, show_col_types = FALSE) |>
    filter(status == "run") |>
    mutate(
        safe_name = gsub("[/ ]", "_", cell_type),
        filename  = glue("DE_{safe_name}_KOvsWT.csv"),
        filepath  = file.path(de_results_dir, filename)
    )

missing_files <- summary_df$filepath[!file.exists(summary_df$filepath)]
if (length(missing_files) > 0) {
    stop(glue(
        "Assertion failed — DE result file(s) missing:\n",
        "{paste(missing_files, collapse = '\n')}"
    ))
}

message(glue("Found {nrow(summary_df)} DE result file(s) to process."))

# =============================================================================
# 3. Enrichment Loop
# =============================================================================

results_rows <- list()

for (i in seq_len(nrow(summary_df))) {

    ct       <- summary_df$cell_type[i]
    filepath <- summary_df$filepath[i]

    res_df <- read_csv(filepath, show_col_types = FALSE)

    if (nrow(res_df) == 0) {
        message(glue("  SKIP [{ct}] — empty DE result file"))
        next
    }

    universe <- res_df |> filter(!is.na(pvalue)) |> pull(gene)

    sig_df <- res_df |>
        filter(!is.na(padj), padj < go_padj, abs(log2FoldChange) > LFC_FILTER)

    genes_up   <- sig_df |> filter(log2FoldChange > 0) |> pull(gene)
    genes_down <- sig_df |> filter(log2FoldChange < 0) |> pull(gene)

    for (direction in c("up", "down")) {

        genes <- if (direction == "up") genes_up else genes_down

        if (length(genes) < go_min_genes) {
            message(glue(
                "  SKIP [{ct}] {direction} — ",
                "{length(genes)} sig genes < go_min_genes ({go_min_genes})"
            ))
            next
        }

        message(glue(
            "  RUN  [{ct}] {direction} — {length(genes)} sig genes, ",
            "{length(universe)} in universe"
        ))

        ego <- enrichGO(
            gene          = genes,
            universe      = universe,
            OrgDb         = org.Mm.eg.db,
            keyType       = "SYMBOL",
            ont           = "BP",
            pAdjustMethod = "BH",
            pvalueCutoff  = go_padj,
            minGSSize     = go_min_genes,
            readable      = FALSE
        )

        ego_df <- as.data.frame(ego)

        if (nrow(ego_df) == 0) {
            message(glue("    → no significant GO terms"))
            next
        }

        message(glue("    → {nrow(ego_df)} significant GO terms"))

        results_rows[[length(results_rows) + 1]] <- ego_df |>
            select(Description, GeneRatio, BgRatio, pvalue, p.adjust, Count) |>
            mutate(cell_type = ct, direction = direction) |>
            relocate(cell_type, direction)
    }
}

# =============================================================================
# 4. Save Combined Results
# =============================================================================

if (length(results_rows) == 0) {
    message("WARNING: no significant GO enrichment found in any comparison.")
    enrichment_tbl <- tibble(
        cell_type   = character(), direction   = character(),
        Description = character(), GeneRatio   = character(),
        BgRatio     = character(), pvalue      = double(),
        p.adjust    = double(),    Count       = integer()
    )
} else {
    enrichment_tbl <- bind_rows(results_rows)
    message(glue("{nrow(enrichment_tbl)} total significant GO terms across all comparisons."))
}

write_csv(enrichment_tbl, enrichment_results_csv)
message(glue("Enrichment results saved to: {enrichment_results_csv}"))

# =============================================================================
# 5. Dotplot — top 5 GO terms per cell type, faceted by direction
# =============================================================================

plot_data <- enrichment_tbl |>
    mutate(GeneRatio_num = map_dbl(GeneRatio, \(x) {
        p <- str_split_1(x, "/")
        as.numeric(p[1]) / as.numeric(p[2])
    })) |>
    group_by(cell_type, direction) |>
    slice_min(order_by = p.adjust, n = 5, with_ties = FALSE) |>
    ungroup()

if (nrow(plot_data) == 0) {
    message("No significant GO terms to plot — saving blank page.")
    p_blank <- ggplot() +
        annotate("text", x = 0.5, y = 0.5,
                 label = "No significant GO terms found.",
                 size = 5) +
        theme_void()
    ggsave(enrichment_dotplot_pdf, plot = p_blank, width = 8, height = 4)
} else {
    n_cell_types <- length(unique(plot_data$cell_type))

    p_dot <- ggplot(plot_data,
                    aes(x = GeneRatio_num,
                        y = reorder(Description, GeneRatio_num),
                        size = Count, color = p.adjust)) +
        geom_point() +
        scale_color_gradient(low = "#d73027", high = "#4575b4",
                             name = "p.adjust") +
        scale_size_continuous(name = "Gene count", range = c(2, 8)) +
        facet_wrap(~ direction, ncol = 2, scales = "free_y") +
        labs(
            title = "GO Biological Process Enrichment (KO vs WT)",
            x     = "Gene Ratio",
            y     = NULL
        ) +
        theme_bw() +
        theme(
            plot.title  = element_text(hjust = 0.5, face = "bold"),
            axis.text.y = element_text(size = 7),
            strip.text  = element_text(size = 9, face = "bold")
        )

    ggsave(enrichment_dotplot_pdf, plot = p_dot,
           width = 14, height = max(8, n_cell_types * 1.5))
    message(glue("Dotplot saved to: {enrichment_dotplot_pdf}"))
}

# =============================================================================
# 6. Summary
# =============================================================================

if (nrow(enrichment_tbl) > 0) {
    n_enriched <- enrichment_tbl |>
        filter(direction == "up") |>
        distinct(cell_type) |>
        nrow()
    message(glue("{n_enriched} cell type(s) had significant upregulated GO enrichment."))

    message("Top finding per cell type (upregulated, by p.adjust):")
    top_per_ct <- enrichment_tbl |>
        filter(direction == "up") |>
        group_by(cell_type) |>
        slice_min(order_by = p.adjust, n = 1, with_ties = FALSE) |>
        ungroup() |>
        select(cell_type, Description, p.adjust, Count)
    print(top_per_ct)
} else {
    message("No significant GO enrichment to summarise.")
}

message("Step 03 complete.")
