###############################################################################
# EnSDD: ensemble spatial-domain identification
# Example: one DLPFC section (slice_id = "507")
#
# Input:
#   Four base clustering label files for the same spots, in identical row order.
#
# Output:
#   EnSDD_507.csv
#
# The ensemble procedure and parameter values follow the supplied analysis code:
#   1. Convert each base clustering into a binary similarity matrix.
#   2. Estimate the consensus similarity matrix with solve_ensemble().
#   3. Obtain final clusters with Louvain_clustering_binary().
#
# No ARI/NMI evaluation or downstream analysis is included.
###############################################################################

rm(list = ls())

library(EnSDD)

# ------------------------------ User settings -------------------------------
slice_id <- "507"
setting_k <- 7

input_dir <- "D:/your_project/ComparativeMethods/input"
output_dir <- "D:/your_project/ComparativeMethods/results"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# One-column CSV files; all four files must contain the same spots
# in exactly the same row order.
method_files <- c(
  pca     = file.path(input_dir, paste0("pca", slice_id, ".csv")),
  gcn     = file.path(input_dir, paste0("gcnpred", slice_id, ".csv")),
  stagate = file.path(input_dir, paste0("stagate", slice_id, ".csv")),
  bayes   = file.path(input_dir, paste0("bayes", slice_id, ".csv"))
)

# ------------------------------ Read labels ---------------------------------
read_cluster_label <- function(file) {
  if (!file.exists(file)) {
    stop("Input file not found: ", file)
  }

  dat <- read.csv(file, check.names = FALSE)

  if (ncol(dat) < 1L) {
    stop("No clustering-label column found in: ", file)
  }

  label <- dat[[1L]]

  if (anyNA(label)) {
    stop("Missing cluster labels found in: ", file)
  }

  as.vector(label)
}

labels <- lapply(method_files, read_cluster_label)

spot_numbers <- lengths(labels)

if (length(unique(spot_numbers)) != 1L) {
  stop(
    "The four methods have different numbers of spots. ",
    "Please align the same spots and row order before EnSDD."
  )
}

# ---------------------- Convert labels to similarity ------------------------
# Retained from the supplied EnSDD analysis code.
clust2Mat <- function(memb) {
  memb <- as.vector(memb)
  N <- length(memb)

  mat_cn <- as.numeric(outer(memb, memb, FUN = "==")) -
    outer(1:N, 1:N, FUN = "==")

  mat_cn <- matrix(mat_cn, nrow = N)
  return(mat_cn)
}

clustering_list <- lapply(labels, clust2Mat)

# --------------------------- EnSDD integration ------------------------------
# Original parameter settings are retained.
res.ensemble <- solve_ensemble(
  clustering_list,
  lambda = NULL,
  prob.quantile = 0.5,
  niter = 100,
  epsilon = 1e-5
)

# ------------------------ Louvain domain assignment -------------------------
ensdd_res <- Louvain_clustering_binary(
  similarity_mat = res.ensemble$H,
  setting_k = setting_k
)

if (length(ensdd_res) != spot_numbers[1L]) {
  stop("The number of EnSDD labels does not match the input spot number.")
}

# ------------------------------ Save labels ---------------------------------
result <- data.frame(
  spot_index = seq_len(length(ensdd_res)),
  EnSDD = as.vector(ensdd_res)
)

write.csv(
  result,
  file.path(output_dir, paste0("EnSDD_", slice_id, ".csv")),
  row.names = FALSE
)

cat(
  "EnSDD completed: ", length(ensdd_res),
  " spots; ", length(unique(ensdd_res)), " domains.\n",
  sep = ""
)

###############################################################################
# Package installation (run once if EnSDD is not installed):
#
# install.packages("devtools")
# devtools::install_github("Zhangxf-ccnu/EnSDD")
###############################################################################
