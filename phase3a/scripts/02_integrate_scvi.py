#!/usr/bin/env python3

# =============================================================================
# Step 02: scVI cross-sample integration
# Converts four per-sample Seurat RDS objects to AnnData, concatenates them,
# trains an scVI model with sample_id as batch key, and exports the latent
# embedding and cell metadata as CSVs for downstream clustering in R.
#
# Usage:
#   python3 02_integrate_scvi.py <output_embedding_csv> <output_metadata_csv>
#                                <sample_meta_str>... <rds_file>...
#
# Arguments:
#   output_embedding_csv  — path for 02_scvi_embedding.csv
#   output_metadata_csv   — path for 02_metadata.csv
#   sample_meta_str...    — "sample_id,genotype,timepoint" strings, one per sample
#                           (order must match rds_file list)
#   rds_file...           — paths to 01_*_prepped.rds files, one per sample
#
# Outputs:
#   02_scvi_embedding.csv — columns: cell_barcode, scVI_1 … scVI_{n_latent}
#   02_metadata.csv       — columns: cell_barcode, sample_id, genotype, timepoint
#
# Environment variables (set in Nextflow script block from params):
#   SCVI_DIMS   — n_latent for scVI model          (default: 30)
#   N_HVGS      — HVGs selected on combined object (default: 3000)
#
# Notes:
#   - RDS → AnnData conversion via rpy2; uses Seurat v5 GetAssayData API.
#   - Raw counts extracted from RNA layer; log-normalised slot ignored.
#   - HVG selection runs on the concatenated object (sc.pp.highly_variable_genes
#     with batch_key="sample_id") per scVI best practice, not per-sample.
#   - Cell barcodes are prefixed with sample_id to guarantee uniqueness.
#   - KMP_DUPLICATE_LIB_OK=TRUE set via conda activate.d.
# =============================================================================

import os
import sys
import numpy as np
import pandas as pd
import scipy.sparse as sp
import anndata as ad
import scanpy as sc
import scvi
import rpy2.robjects as ro
from rpy2.robjects.packages import importr

SEED = 42
scvi.settings.seed = SEED

# =============================================================================
# 1. Arguments & Validation
# =============================================================================

args = sys.argv[1:]

if len(args) < 4:
    sys.exit(
        "Usage: python3 02_integrate_scvi.py <output_embedding_csv> "
        "<output_metadata_csv> <sample_meta_str>... <rds_file>...\n"
        f"Got {len(args)} argument(s)."
    )

output_embedding_csv = args[0]
output_metadata_csv  = args[1]
remaining            = args[2:]

sample_metas = [a for a in remaining if "," in a and not a.endswith(".rds")]
rds_paths    = [a for a in remaining if a.endswith(".rds")]

if len(sample_metas) != len(rds_paths):
    sys.exit(
        f"Mismatch: {len(sample_metas)} sample_meta string(s) vs "
        f"{len(rds_paths)} RDS file(s)."
    )

parsed_meta = []
for m in sample_metas:
    parts = m.split(",")
    if len(parts) != 3:
        sys.exit(f"sample_meta must be 'sample_id,genotype,timepoint'; got: '{m}'")
    parsed_meta.append({"sample_id": parts[0], "genotype": parts[1], "timepoint": parts[2]})

n_latent = int(os.environ.get("SCVI_DIMS", 30))
n_hvgs   = int(os.environ.get("N_HVGS", 3000))

print(
    f"\n========== scVI Integration (Step 02) ==========\n"
    f"  Samples     : {len(rds_paths)}\n"
    f"  n_latent    : {n_latent}\n"
    f"  n_hvgs      : {n_hvgs}\n"
    f"  Seed        : {SEED}\n"
    f"  Output emb  : {output_embedding_csv}\n"
    f"  Output meta : {output_metadata_csv}\n"
    f"=================================================\n"
)

# =============================================================================
# 2. Load Each RDS and Convert to AnnData
# =============================================================================

try:
    importr("Seurat")
    print("Seurat loaded successfully in R session.")
except Exception as e:
    sys.exit(
        f"Failed to load Seurat in R session: {e}\n"
        f"Ensure Seurat is installed in the conda environment's R library."
    )


def rds_to_anndata(rds_path, meta):
    """Load a Seurat v5 RDS and return an AnnData with raw counts."""
    R = ro.r
    print(f"  Loading {rds_path} ...")
    R(f'obj_ <- readRDS("{rds_path}")')
    R('mat_ <- GetAssayData(obj_, assay = "RNA", layer = "counts")')

    n_genes = int(R("nrow(mat_)")[0])
    n_cells = int(R("ncol(mat_)")[0])

    i_vec = np.array(R("mat_@i"), dtype=np.int32)
    p_vec = np.array(R("mat_@p"), dtype=np.int32)
    x_vec = np.array(R("mat_@x"), dtype=np.float32)
    genes = list(R("rownames(mat_)"))
    cells = [f"{meta['sample_id']}_{bc}" for bc in list(R("colnames(mat_)"))]

    # dgCMatrix is genes × cells; transpose to cells × genes for AnnData
    counts = sp.csc_matrix(
        (x_vec, i_vec, p_vec), shape=(n_genes, n_cells)
    ).T.tocsr()

    obs = pd.DataFrame(
        {
            "sample_id": meta["sample_id"],
            "genotype":  meta["genotype"],
            "timepoint": meta["timepoint"],
        },
        index=cells,
    )

    adata = ad.AnnData(X=counts, obs=obs, var=pd.DataFrame(index=genes))
    adata.var_names_make_unique()
    print(f"    → {adata.n_obs} cells × {adata.n_vars} genes")
    return adata


print("Loading and converting RDS files...")
adatas = [
    rds_to_anndata(path, meta)
    for path, meta in zip(rds_paths, parsed_meta)
]

# =============================================================================
# 3. Concatenate & Batch-aware HVG Selection
# =============================================================================

print("\nConcatenating samples...")
combined = ad.concat(adatas, join="inner")
combined.obs_names_make_unique()
print(f"Combined: {combined.n_obs} cells × {combined.n_vars} genes")

print("\nCells per sample:")
print(combined.obs["sample_id"].value_counts().to_string())

print(f"\nSelecting {n_hvgs} HVGs (seurat_v3, batch_key=sample_id)...")
sc.pp.highly_variable_genes(
    combined,
    flavor      = "seurat_v3",
    n_top_genes = n_hvgs,
    batch_key   = "sample_id",
    subset      = True,
)
print(f"Retained {combined.n_vars} HVGs after batch-aware selection.")

# =============================================================================
# 4. Set Up & Train scVI
# =============================================================================

print("\nSetting up scVI model...")
scvi.model.SCVI.setup_anndata(
    combined,
    batch_key                  = "sample_id",
    categorical_covariate_keys = ["genotype"],
)

model = scvi.model.SCVI(combined, n_latent=n_latent, n_layers=2)
print(model)

print("\nTraining (max_epochs=400, early_stopping=True)...")
model.train(max_epochs=400, early_stopping=True)

elbo = model.history["elbo_train"]
print(
    f"\nTraining ELBO — "
    f"first: {elbo.iloc[0, 0]:.2f}  "
    f"last:  {elbo.iloc[-1, 0]:.2f}  "
    f"epochs run: {len(elbo)}"
)

# =============================================================================
# 5. Extract Latent Representation
# =============================================================================

print(f"\nExtracting latent representation ({n_latent} dims)...")
latent = model.get_latent_representation()
print(f"Latent array shape: {latent.shape}")

# =============================================================================
# 6. Hard Assertion
# =============================================================================

assert latent.shape[0] == combined.n_obs, (
    f"Assertion failed: embedding has {latent.shape[0]} rows "
    f"but metadata has {combined.n_obs} rows."
)
print(f"Assertion passed: {latent.shape[0]} cells in embedding and metadata.")

# =============================================================================
# 7. Save Outputs
# =============================================================================

cell_barcodes = combined.obs_names.tolist()

dim_cols = [f"scVI_{i + 1}" for i in range(n_latent)]
emb_df   = pd.DataFrame(latent, index=cell_barcodes, columns=dim_cols)
emb_df.index.name = "cell_barcode"
emb_df.to_csv(output_embedding_csv)
print(
    f"\nEmbedding saved : {output_embedding_csv}  "
    f"({emb_df.shape[0]} cells × {emb_df.shape[1]} dims)"
)

meta_df = combined.obs[["sample_id", "genotype", "timepoint"]].copy()
meta_df.index.name = "cell_barcode"
meta_df.to_csv(output_metadata_csv)
print(
    f"Metadata saved  : {output_metadata_csv}  "
    f"({meta_df.shape[0]} cells × {meta_df.shape[1]} columns)"
)

print("\nStep 02 complete.")
