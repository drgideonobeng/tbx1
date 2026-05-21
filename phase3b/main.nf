nextflow.enable.dsl=2

// =============================================================================
// main.nf — Phase 3b: Annotation, differential expression & enrichment
// Linear DAG: ANNOTATE_ATLAS → DIFFERENTIAL_EXPRESSION → ENRICHMENT
// =============================================================================

// --- PROCESS 1: Annotate integrated atlas with cluster labels ---
process ANNOTATE_ATLAS {
    publishDir "${params.resultsdir}/objects", mode: 'copy', pattern: "*.rds"
    publishDir "${params.resultsdir}/plots",   mode: 'copy', pattern: "*.pdf"
    input:
        path clustered_rds
        path cluster_labels
    output:
        path "01_annotated_atlas.rds", emit: rds
        path "01_umap_annotated.pdf"
        path "01_umap_by_genotype.pdf"
        path "01_umap_by_timepoint.pdf"
    script:
        """
        Rscript ${projectDir}/scripts/01_annotate_atlas.R \
            ${clustered_rds} \
            ${cluster_labels} \
            01_annotated_atlas.rds \
            01_umap_annotated.pdf \
            01_umap_by_genotype.pdf \
            01_umap_by_timepoint.pdf
        """
}

// --- PROCESS 2: Differential expression (WT vs KO, pseudobulk per cluster) ---
process DIFFERENTIAL_EXPRESSION {
    publishDir "${params.resultsdir}/tables",             mode: 'copy', pattern: "de_results/**"
    publishDir "${params.resultsdir}/tables",             mode: 'copy', pattern: "*.csv"
    input:
        path annotated_rds
    output:
        path "de_results/",        emit: de_dir        
        path "02_DE_summary.csv",  emit: summary
    script:
        """
        export LFC_THRESHOLD=${params.lfc_threshold}
        export PADJ_THRESHOLD=${params.padj_threshold}
        export MIN_CELLS=${params.min_cells}
        Rscript ${projectDir}/scripts/02_differential_expression.R \
            ${annotated_rds} \
            de_results/ \
            02_DE_summary.csv
        """
}

// --- PROCESS 3: GO enrichment ---
process ENRICHMENT {
    publishDir "${params.resultsdir}/tables", mode: 'copy', pattern: "*.csv"
    publishDir "${params.resultsdir}/plots",  mode: 'copy', pattern: "*.pdf"
    input:
        path de_dir
        path de_summary
    output:
        path "03_enrichment_results.csv"
        path "03_enrichment_dotplot.pdf"
    script:
        """
        export GO_PADJ=${params.go_padj}
        export GO_MIN_GENES=${params.go_min_genes}
        Rscript ${projectDir}/scripts/03_enrichment.R \
            ${de_dir} \
            ${de_summary} \
            03_enrichment_results.csv \
            03_enrichment_dotplot.pdf
        """
}

// =============================================================================
// WORKFLOW
// =============================================================================

workflow {

    clustered_rds  = file(params.clustered_rds)
    cluster_labels = file(params.cluster_labels)

    if (!clustered_rds.exists())
        error "Input not found — clustered_rds: ${params.clustered_rds}"
    if (!cluster_labels.exists())
        error "Input not found — cluster_labels: ${params.cluster_labels}"

    ANNOTATE_ATLAS(clustered_rds, cluster_labels)
    DIFFERENTIAL_EXPRESSION(ANNOTATE_ATLAS.out.rds)
    ENRICHMENT(
        DIFFERENTIAL_EXPRESSION.out.de_dir,
        DIFFERENTIAL_EXPRESSION.out.summary
    )
}
