###############################################################################
# STAGATE: spatial domain identification
# Example: DLPFC section 151507
#
# Official DLPFC tutorial:
# https://stagate.readthedocs.io/en/latest/T1_DLPFC.html
#
# Output:
#   STAGATE_151507.csv
#
# This script stops after spot-level spatial-domain assignments.
# ARI/NMI evaluation and downstream analyses are not included.
###############################################################################

from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
import STAGATE

# ------------------------------ User settings -------------------------------
SLICE_ID = "151507"
N_CLUSTERS = 7  # DLPFC: 7 for 151507-151510 and 151673-151676;
                #       5 for 151669-151672

SLICE_DIR = Path(r"D:\your_DLPFC\151507")  # Change to this section's directory
H5_FILE = SLICE_DIR / f"{SLICE_ID}_filtered_feature_bc_matrix.h5"
POSITION_FILE = SLICE_DIR / "tissue_positions.csv"

OUTPUT_DIR = SLICE_DIR / "results"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# ------------------------------ Load data -----------------------------------
adata = sc.read_10x_h5(str(H5_FILE))
adata.var_names_make_unique()

positions = pd.read_csv(POSITION_FILE, header=None)

if str(positions.iloc[0, 0]).strip().lower() == "barcode":
    positions = positions.iloc[1:].reset_index(drop=True)

positions = positions.iloc[:, :6].copy()
positions.columns = [
    "barcode", "in_tissue", "array_row", "array_col",
    "pxl_row", "pxl_col"
]
positions["barcode"] = positions["barcode"].astype(str).str.strip()
positions = positions.set_index("barcode")
positions["in_tissue"] = pd.to_numeric(positions["in_tissue"])

# Keep in-tissue spots and align coordinates to expression barcodes
positions = positions.loc[positions["in_tissue"] == 1]
barcodes = [barcode for barcode in adata.obs_names if barcode in positions.index]

if not barcodes:
    raise ValueError("No matching barcodes between expression data and spatial coordinates.")

adata = adata[barcodes].copy()
positions = positions.loc[barcodes]

# Use full-resolution image coordinates, as in the Visium spatial input
adata.obsm["spatial"] = np.column_stack(
    [
        pd.to_numeric(positions["pxl_col"]).to_numpy(),
        pd.to_numeric(positions["pxl_row"]).to_numpy()
    ]
)

# -------------------------- Expression preprocessing ------------------------
sc.pp.highly_variable_genes(
    adata,
    flavor="seurat_v3",
    n_top_genes=3000
)
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# -------------------------- STAGATE workflow --------------------------------
STAGATE.Cal_Spatial_Net(
    adata,
    rad_cutoff=150
)

adata = STAGATE.train_STAGATE(
    adata,
    alpha=0
)

# Cluster the STAGATE representation
adata = STAGATE.mclust_R(
    adata,
    used_obsm="STAGATE",
    num_cluster=N_CLUSTERS
)

# ------------------------------ Save labels ---------------------------------
result = pd.DataFrame(
    {
        "barcode": adata.obs_names,
        "spatial_domain": adata.obs["mclust"].astype(str).to_numpy()
    }
)

result.to_csv(
    OUTPUT_DIR / f"STAGATE_{SLICE_ID}.csv",
    index=False
)

print("STAGATE spatial-domain identification completed.")
