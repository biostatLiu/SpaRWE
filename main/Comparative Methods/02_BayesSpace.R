###############################################################################
# BayesSpace: spot-level spatial domain identification
# Example: DLPFC section 151507
#
# Official documentation:
# https://www.bioconductor.org/packages/release/bioc/vignettes/BayesSpace/inst/doc/BayesSpace.html
#
# Output:
#   BayesSpace_151507.csv
#
# This script runs spatialCluster() only.
# It does NOT run spatialEnhance() or subspot-level resolution enhancement.
###############################################################################

rm(list = ls())

library(Seurat)
library(SingleCellExperiment)
library(S4Vectors)
library(BayesSpace)

# ------------------------------ User settings -------------------------------
slice_id <- "151507"
K <- 7                         # DLPFC: 7 for 151507-151510 and 151673-151676;
                               #       5 for 151669-151672

slice_dir <- "D:/your_DLPFC/151507"  # Change to the directory of this section
h5_file <- file.path(slice_dir, paste0(slice_id, "_filtered_feature_bc_matrix.h5"))
position_file <- file.path(slice_dir, "tissue_positions.csv")
output_dir <- file.path(slice_dir, "results")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# ------------------------------ Load data -----------------------------------
counts <- Seurat::Read10X_h5(h5_file)

if (is.list(counts)) {
  counts <- counts[["Gene Expression"]]
}
if (is.null(counts)) {
  stop("Gene Expression matrix was not found in the input H5 file.")
}

positions <- read.csv(
  position_file,
  header = FALSE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

if (tolower(trimws(as.character(positions[1, 1]))) == "barcode") {
  positions <- positions[-1, , drop = FALSE]
}

positions <- positions[, 1:6, drop = FALSE]
colnames(positions) <- c(
  "barcode", "in_tissue", "array_row", "array_col",
  "pxl_row", "pxl_col"
)

positions$barcode <- trimws(as.character(positions$barcode))
positions$in_tissue <- as.numeric(positions$in_tissue)
positions$array_row <- as.numeric(positions$array_row)
positions$array_col <- as.numeric(positions$array_col)

positions <- positions[
  !is.na(positions$in_tissue) & positions$in_tissue == 1,
  ,
  drop = FALSE
]
rownames(positions) <- positions$barcode

# Match expression data and spatial coordinates by barcode
barcodes <- colnames(counts)[colnames(counts) %in% rownames(positions)]
if (length(barcodes) == 0) {
  stop("No matching barcodes between expression data and spatial coordinates.")
}

counts <- counts[, barcodes, drop = FALSE]
positions <- positions[barcodes, , drop = FALSE]

# BayesSpace requires array coordinates in colData
col_data <- S4Vectors::DataFrame(
  array_row = positions$array_row,
  array_col = positions$array_col,
  row.names = barcodes
)

sce <- SingleCellExperiment(
  assays = list(counts = counts),
  colData = col_data
)

# -------------------------- BayesSpace workflow -----------------------------
set.seed(104)

sce <- spatialPreprocess(
  sce,
  platform = "Visium",
  n.PCs = 15,
  n.HVGs = 2000,
  log.normalize = TRUE
)

# Spot-level spatial clustering
sce <- spatialCluster(
  sce,
  q = K,
  d = 15,
  platform = "Visium",
  init.method = "mclust",
  model = "t",
  nrep = 10000,
  burn.in = 1000,
  gamma = 3
)

# ------------------------------ Save labels ---------------------------------
result <- data.frame(
  barcode = colnames(sce),
  spatial_domain = as.character(colData(sce)$spatial.cluster)
)

write.csv(
  result,
  file.path(output_dir, paste0("BayesSpace_", slice_id, ".csv")),
  row.names = FALSE
)

cat("BayesSpace spot-level spatial-domain identification completed.\n")
