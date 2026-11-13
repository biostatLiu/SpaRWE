###############################################################################
# SpatialPCA: spatial domain identification
# Example: DLPFC section 151507
#
# Official tutorial:
# https://lulushang.org/SpatialPCA_Tutorial/DLPFC.html
#
# Output:
#   SpatialPCA_151507.csv
#
# This script stops after spot-level spatial domain identification.
###############################################################################

rm(list = ls())

library(Seurat)
library(SpatialPCA)

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

# Match the count matrix and spatial coordinates by barcode
barcodes <- colnames(counts)[colnames(counts) %in% rownames(positions)]
if (length(barcodes) == 0) {
  stop("No matching barcodes between expression data and spatial coordinates.")
}

counts <- counts[, barcodes, drop = FALSE]
positions <- positions[barcodes, , drop = FALSE]

location <- as.matrix(
  positions[, c("array_row", "array_col"), drop = FALSE]
)
storage.mode(location) <- "numeric"
rownames(location) <- barcodes

# ------------------------- SpatialPCA workflow ------------------------------
set.seed(123)

SpatialPCA_obj <- CreateSpatialPCAObject(
  counts = counts,
  location = location,
  project = "SpatialPCA",
  gene.type = "spatial",
  sparkversion = "spark",
  numCores_spark = 5,
  gene.number = 3000,
  customGenelist = NULL
)

SpatialPCA_obj <- SpatialPCA_buildKernel(
  SpatialPCA_obj,
  kerneltype = "gaussian",
  bandwidthtype = "SJ",
  bandwidth.set.by.user = NULL
)

SpatialPCA_obj <- SpatialPCA_EstimateLoading(
  SpatialPCA_obj,
  fast = FALSE,
  SpatialPCnum = 20
)

SpatialPCA_obj <- SpatialPCA_SpatialPCs(
  SpatialPCA_obj,
  fast = FALSE
)

# Spatial-domain identification using Walktrap clustering
spatial_domain <- walktrap_clustering(
  clusternum = K,
  latent_dat = SpatialPCA_obj@SpatialPCs,
  knearest = 70
)

# ------------------------------ Save labels ---------------------------------

result <- data.frame(
  barcode = barcodes,
  spatial_domain = as.character(spatial_domain)
)

write.csv(
  result,
  file.path(output_dir, paste0("SpatialPCA_", slice_id, ".csv")),
  row.names = FALSE
)

cat("SpatialPCA spatial-domain identification completed.\n")
