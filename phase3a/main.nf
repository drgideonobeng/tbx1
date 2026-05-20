nextflow.enable.dsl=2

// =============================================================================
// main.nf — Phase 3a: Cross-sample integration (Tbx1 WT vs KO, E8.5 + E9.5)
// 4-step workflow: per-sample prep/normalize → scVI integration →
// cluster integrated → clustree resolution sweep
// =============================================================================

// --- PROCESS 1: Per-sample prep & normalization (log-normalize + HVG selection) ---
process PREP_NORMALIZE {
    tag "${sample_id}"
    publishDir "${params.resultsdir}/objects", mode: 'copy', pattern: "*.rds"
    input:
        tuple val(sample_id), val(genotype), val(timepoint), path(rds)
    output:
        tuple val(sample_id), val(genotype), val(timepoint),
              path("01_${sample_id}_prepped.rds"), emit: rds
    script:
        """
        Rscript ${projectDir}/scripts/01_prep_normalize.R \
            ${rds} ${sample_id} ${genotype} ${timepoint} \
            ${params.n_hvgs} \
            01_${sample_id}_prepped.rds
        """
}

// --- PROCESS 2: scVI integration (Python / scvi-tools) ---
process INTEGRATE_SCVI {
    publishDir "${params.resultsdir}/tables",  mode: 'copy', pattern: "*.csv"
    input:
        path rds_files
        val  sample_meta   // list of "sample_id,genotype,timepoint" strings
    output:
        path "02_scvi_embedding.csv", emit: embedding
        path "02_metadata.csv",       emit: metadata
    script:
        """
        export SCVI_DIMS=${params.scvi_dims}
        export N_HVGS=${params.n_hvgs}
        python3 ${projectDir}/scripts/02_integrate_scvi.py \
            02_scvi_embedding.csv \
            02_metadata.csv \
            ${sample_meta.join(' ')} \
            ${rds_files}
        """
}

// --- PROCESS 3: Cluster integrated embedding (R / Seurat) ---
process CLUSTER_INTEGRATED {
    publishDir "${params.resultsdir}/objects", mode: 'copy', pattern: "*.rds"
    publishDir "${params.resultsdir}/tables",  mode: 'copy', pattern: "*.csv"
    publishDir "${params.resultsdir}/plots",   mode: 'copy', pattern: "*.pdf"
    
    input:
        path embedding_csv
        path metadata_csv
        path prepped_rds   // four 01_*_prepped.rds files, ordered: wt_e85 wt_e95 ko_e85 ko_e95
    output:
        path "03_integrated_clustered.rds", emit: rds
        path "03_umap_clusters.pdf"
        path "03_umap_samples.pdf"
        path "03_umap_genotype.pdf"
        path "03_composition.csv"
    script:
        """
        export CLUSTER_RES=${params.cluster_res}
        export SCVI_DIMS=${params.scvi_dims}
        Rscript ${projectDir}/scripts/03_cluster_integrated.R \
            ${embedding_csv} ${metadata_csv} \
            ${prepped_rds} \
            03_integrated_clustered.rds \
            03_umap_clusters.pdf \
            03_umap_samples.pdf \
            03_umap_genotype.pdf \
            03_composition.csv
        """
}

// --- PROCESS 4: Clustree resolution stability sweep ---
process RUN_CLUSTREE {
    publishDir "${params.resultsdir}/plots", mode: 'copy'
    input:
        path rds
    output:
        path "04_clustree_resolutions.pdf"
    script:
        """
        Rscript ${projectDir}/scripts/04_run_clustree.R \
            ${rds} \
            04_clustree_resolutions.pdf
        """
}

// =============================================================================
// WORKFLOW
// =============================================================================

workflow {

    samples_ch = Channel.of(
        ["wt_e85", "WT", "E85", file(params.wt_e85_rds)],
        ["wt_e95", "WT", "E95", file(params.wt_e95_rds)],
        ["ko_e85", "KO", "E85", file(params.ko_e85_rds)],
        ["ko_e95", "KO", "E95", file(params.ko_e95_rds)]
    )

    PREP_NORMALIZE(samples_ch)

    // Collect prepped RDS files in deterministic order for CLUSTER_INTEGRATED
    prepped_rds_ch = PREP_NORMALIZE.out.rds
        .map { id, _geno, _time, rds ->
            def order = ["wt_e85": 0, "wt_e95": 1, "ko_e85": 2, "ko_e95": 3]
            tuple(order[id], rds)
        }
        .toSortedList { a, b -> a[0] <=> b[0] }
        .map { sorted -> sorted.collect { it[1] } }

    INTEGRATE_SCVI(
        PREP_NORMALIZE.out.rds.map { _id, _geno, _time, rds -> rds }.collect(),
        PREP_NORMALIZE.out.rds.map { id, geno, time, _rds ->
            "${id},${geno},${time}" }.collect()
    )

    CLUSTER_INTEGRATED(
        INTEGRATE_SCVI.out.embedding,
        INTEGRATE_SCVI.out.metadata,
        prepped_rds_ch
    )

    RUN_CLUSTREE(CLUSTER_INTEGRATED.out.rds)
}
