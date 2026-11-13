###############################################################################
# SpaGCN: spatial domain identification
# Example: DLPFC section 151507
#
# Official tutorial:
# https://github.com/jianhuupenn/SpaGCN/blob/master/tutorial/tutorial.md
#
# Output:
#   SpaGCN_151507.csv
#
# The workflow stops after the initial spatial-domain assignments.
# Optional domain refinement and downstream analyses are not included.
###############################################################################

from pathlib import Path
import json
import random

import cv2
import numpy as np
import pandas as pd
import scanpy as sc
import torch
import SpaGCN as spg

# ------------------------------ User settings -------------------------------
SLICE_ID = "151507"
N_CLUSTERS = 7  # DLPFC: 7 for 151507-151510 and 151673-151676;
                #       5 for 151669-151672

SLICE_DIR = Path(r"D:\your_DLPFC\151507")  # Change to this section's directory

H5_FILE = SLICE_DIR / f"{SLICE_ID}_filtered_feature_bc_matrix.h5"
POSITION_FILE = SLICE_DIR / "tissue_positions.csv"
IMAGE_FILE = SLICE_DIR / "spatial" / "tissue_hires_image.png"
SCALE_FILE = SLICE_DIR / "spatial" / "scalefactors_json.json"

OUTPUT_DIR = SLICE_DIR / "results"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# Official tutorial parameters
S = 1
B = 49
P = 0.5
R_SEED = 100
T_SEED = 100
N_SEED = 100

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

# Keep in-tissue spots and align all inputs by barcode
positions["in_tissue"] = pd.to_numeric(positions["in_tissue"])
positions = positions.loc[positions["in_tissue"] == 1]

barcodes = [barcode for barcode in adata.obs_names if barcode in positions.index]
if not barcodes:
    raise ValueError("No matching barcodes between expression data and spatial coordinates.")

adata = adata[barcodes].copy()
positions = positions.loc[barcodes]

# Array coordinates and image pixel coordinates
x_array = pd.to_numeric(positions["array_row"]).to_numpy()
y_array = pd.to_numeric(positions["array_col"]).to_numpy()

with open(SCALE_FILE, "r", encoding="utf-8") as f:
    hires_scale = json.load(f)["tissue_hires_scalef"]

# Space Ranger full-resolution coordinates are scaled to the hires image
x_pixel = pd.to_numeric(positions["pxl_row"]).to_numpy() * hires_scale
y_pixel = pd.to_numeric(positions["pxl_col"]).to_numpy() * hires_scale

img = cv2.imread(str(IMAGE_FILE))
if img is None:
    raise FileNotFoundError(f"Histology image could not be read: {IMAGE_FILE}")

# ---------------------- Construct spatial adjacency -------------------------
adj = spg.calculate_adj_matrix(
    x=x_pixel,
    y=y_pixel,
    x_pixel=x_pixel,
    y_pixel=y_pixel,
    image=img,
    beta=B,
    alpha=S,
    histology=True
)

# -------------------------- Expression preprocessing ------------------------
spg.prefilter_genes(adata, min_cells=3)
spg.prefilter_specialgenes(adata)

sc.pp.normalize_per_cell(adata)
sc.pp.log1p(adata)

# -------------------------- SpaGCN clustering -------------------------------
l = spg.search_l(
    P,
    adj,
    start=0.01,
    end=1000,
    tol=0.01,
    max_run=100
)

random.seed(R_SEED)
torch.manual_seed(T_SEED)
np.random.seed(N_SEED)

res = spg.search_res(
    adata,
    adj,
    l,
    N_CLUSTERS,
    start=0.7,
    step=0.1,
    tol=5e-3,
    lr=0.05,
    max_epochs=20,
    r_seed=R_SEED,
    t_seed=T_SEED,
    n_seed=N_SEED
)

model = spg.SpaGCN()
model.set_l(l)

model.train(
    adata,
    adj,
    init_spa=True,
    init="louvain",
    res=res,
    tol=5e-3,
    lr=0.05,
    max_epochs=200
)

y_pred, _ = model.predict()
adata.obs["spatial_domain"] = pd.Categorical(y_pred)

# ------------------------------ Save labels ---------------------------------
result = pd.DataFrame(
    {
        "barcode": adata.obs_names,
        "spatial_domain": adata.obs["spatial_domain"].astype(str).to_numpy()
    }
)

result.to_csv(
    OUTPUT_DIR / f"SpaGCN_{SLICE_ID}.csv",
    index=False
)

print("SpaGCN spatial-domain identification completed.")
