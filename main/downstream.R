###############################################################################
# SpaRWE: Downstream analysis of spatial transcriptomics data
#
# Analysis modules:
# 1. Pseudotime and trajectory analysis (Monocle3)
# 2. Pseudotime-associated genes and expression heatmap
# 3. Representative gene trends along pseudotime
# 4. Spatial gene expression visualization
# 5. Spatial-domain differential expression analysis
# 6. Differential-expression heatmap
# 7. GO and KEGG over-representation analysis (ORA)
# 8. GO biological-process chord diagram
# 9. GO and KEGG GSEA
# 10. ORA-GSEA pathway intersections
# Each analysis module retains no more than one representative visualization.
#
# Before running:
#   Set the input paths in Section 1.
#   The input expression matrix, spatial coordinates, domain labels, and
#   cell-type annotation must refer to the same spots in the same order.
###############################################################################

rm(list = ls())

###############################################################################
# 1. Packages and paths
###############################################################################

library(Seurat)
library(monocle3)
library(dplyr)
library(ggplot2)
library(viridis)
library(ClusterGVis)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(stringr)
library(GOplot)
library(DOSE)

# Project directory: edit this path for your computer
project_dir <- "D:/your_project/SpaRWE_downstream"

# Input files
count_dir   <- file.path(project_dir, "data", "filtered_feature_bc_matrix")
spatial_file <- file.path(project_dir, "data", "spatial", "tissue_positions1.csv")
label_file  <- file.path(project_dir, "data", "WEAC.csv")
cortex_file <- file.path(project_dir, "data", "cortex.RData")

# Output directory
output_dir <- file.path(project_dir, "results")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Check input paths
input_files <- c(count_dir, spatial_file, label_file, cortex_file)
if (!dir.exists(count_dir)) {
  stop("Expression-matrix directory not found: ", count_dir)
}
if (!file.exists(spatial_file)) {
  stop("Spatial-coordinate file not found: ", spatial_file)
}
if (!file.exists(label_file)) {
  stop("Spatial-domain label file not found: ", label_file)
}
if (!file.exists(cortex_file)) {
  stop("Cell-type annotation file not found: ", cortex_file)
}

###############################################################################
# 2. Load input data
###############################################################################

# Gene-expression matrix
count_sub <- Read10X(data.dir = count_dir)
count_sub1 <- as.matrix(count_sub)

# Spatial-domain labels
label <- read.csv(label_file, header = FALSE)
spatial_domain <- label$V1

# Spatial coordinates
xy_coords <- read.csv(
  spatial_file,
  header = FALSE,
  row.names = 1
)

# Cell-type annotation
load(cortex_file)  # Expected object: cortex
if (!exists("cortex")) {
  stop("The loaded RData file does not contain an object named 'cortex'.")
}
celltype <- cortex@meta.data[["predicted_celltype"]]

# Check sample numbers
stopifnot(length(spatial_domain) == ncol(count_sub1))
stopifnot(length(celltype) == ncol(count_sub1))

###############################################################################
# 3. Pseudotime and trajectory analysis (Monocle3)
###############################################################################

# Cell metadata
cell_metadata <- data.frame(
  cell = colnames(count_sub1),
  cell_type = celltype,
  spatial_domain = spatial_domain,
  row.names = colnames(count_sub1)
)

# Gene metadata
gene_metadata <- data.frame(
  gene_short_name = rownames(count_sub1),
  row.names = rownames(count_sub1)
)

# Create Monocle3 cell_data_set
cds <- new_cell_data_set(
  expression_data = count_sub1,
  cell_metadata = cell_metadata,
  gene_metadata = gene_metadata
)

# Preprocessing, dimensionality reduction, and clustering
cds <- preprocess_cds(cds, num_dim = 50)
cds <- reduce_dimension(cds, reduction_method = "UMAP")
cds <- cluster_cells(cds)

# Learn trajectory graph
cds <- learn_graph(cds, use_partition = FALSE)

# Prepare spatial-domain metadata for visualization
colData(cds)$spatial_domain <- factor(colData(cds)$spatial_domain)

# Order cells along pseudotime
# Root selection follows the original analysis setting.
cds <- order_cells(cds)

# Pseudotime visualization
plot_cells(
  cds,
  color_cells_by = "pseudotime",
  label_cell_groups = FALSE,
  label_leaves = FALSE,
  label_branch_points = FALSE
)

###############################################################################
# 4. Pseudotime-associated genes and heatmap
###############################################################################

# Test genes along the principal graph
genes <- graph_test(
  cds,
  neighbor_graph = "principal_graph",
  reduction_method = "UMAP"
)

# Select the top 50 significant genes ranked by Moran's I
top50 <- genes %>%
  filter(q_value < 0.05) %>%
  top_n(n = 50, morans_I) %>%
  pull(gene_short_name) %>%
  as.character()

# Extract and normalize expression values
raw_mat <- log1p(
  assay(cds, "counts")[top50, ] /
    Matrix::colSums(assay(cds, "counts")) * 10000
)

# Order samples by pseudotime
pseudotime_order <- order(pseudotime(cds))
raw_mat <- raw_mat[, pseudotime_order]

# Convert to a standard matrix
mat <- matrix(
  data = as.numeric(raw_mat),
  nrow = nrow(raw_mat),
  ncol = ncol(raw_mat),
  dimnames = list(rownames(raw_mat), colnames(raw_mat))
)
attr(mat, "class") <- "matrix"

# K-means clustering of pseudotime-associated genes
ck <- clusterData(
  mat,
  cluster.method = "kmeans",
  cluster.num = 5
)

pdf(
  file = file.path(output_dir, "monocle3_spatial_pseudotime.pdf"),
  height = 10,
  width = 8,
  onefile = FALSE
)

visCluster(
  object = ck,
  plot.type = "both",
  add.sampleanno = FALSE,
  markGenes = sample(rownames(mat), 30, replace = FALSE)
)

dev.off()

###############################################################################
# 5. Representative gene expression trends along pseudotime
###############################################################################

# Representative genes used in the original analysis
target_genes <- c("SLC45A3", "CXCL13", "TGM4", "SPDEF", "VEGFA")

AFD_lineage_cds <- cds[target_genes, ]

plot_genes_in_pseudotime(
  AFD_lineage_cds,
  color_cells_by = "spatial_domain",
  min_expr = 0.5
)


###############################################################################
# 6. Spatial gene expression visualization
###############################################################################

# Create Seurat object
brain <- CreateSeuratObject(
  counts = count_sub1,
  project = "Prostate",
  min.cells = 1,
  min.features = 1
)

# Align spatial coordinates with expression-matrix barcodes
xy_coords <- as.matrix(xy_coords)
xy_coords <- xy_coords[colnames(brain), , drop = FALSE]

# Add spatial coordinates to metadata
brain@meta.data$tissue_x <- xy_coords[, 1]
brain@meta.data$tissue_y <- xy_coords[, 2]

# Create Spatial assay and image object
brain[["Spatial"]] <- CreateAssayObject(counts = count_sub1)

brain@images[["slice1"]] <- new(
  Class = "SlideSeq",
  assay = "Spatial",
  key = "slice1_",
  coordinates = as.data.frame(xy_coords)
)

# SCTransform normalization
brain <- SCTransform(
  brain,
  assay = "Spatial",
  verbose = FALSE
)

DefaultAssay(brain) <- "SCT"

# Save processed Seurat object
save(
  brain,
  file = file.path(output_dir, "brain.RData")
)

# Spatial expression visualization
gene_to_plot <- c("ORM2", "SV2B", "NELL2", "MUC12", "LTF", "PTGS1")

SpatialFeaturePlot(
  brain,
  features = gene_to_plot,
  images = "slice1",
  alpha = 0.8
) +
  scale_fill_viridis_c(option = "A") +
  ggtitle(paste0("Spatial expression of ", gene_to_plot)) +
  theme(plot.title = element_text(hjust = 0.5, size = 16))

###############################################################################
# 7. Spatial-domain differential expression analysis
###############################################################################

seurat_obj <- CreateSeuratObject(
  counts = count_sub,
  project = "Prostate",
  min.cells = 1,
  min.features = 1
)

# Assign spatial-domain labels
seurat_obj$cluster <- as.factor(spatial_domain)
Idents(seurat_obj) <- seurat_obj$cluster

# Normalization and scaling
seurat_obj <- NormalizeData(seurat_obj)
seurat_obj <- FindVariableFeatures(
  seurat_obj,
  selection.method = "vst",
  nfeatures = 2000
)
seurat_obj <- ScaleData(
  seurat_obj,
  features = rownames(seurat_obj)
)

# Identify markers for all spatial domains
markers <- FindAllMarkers(
  seurat_obj,
  assay = "RNA",
  slot = "data",
  test.use = "MAST",
  min.pct = 0.1,
  logfc.threshold = 0.25
)

# Filter significant differentially expressed genes
count676R <- markers %>%
  filter(p_val_adj < 0.05 & abs(avg_log2FC) > 0.25) %>%
  select(
    SYMBOL = gene,
    Log2FoldChange = avg_log2FC,
    pvalue = p_val,
    padj = p_val_adj
  )

write.csv(
  count676R,
  file.path(output_dir, "count676R.csv"),
  row.names = TRUE
)

write.csv(
  markers,
  file.path(output_dir, "markers.csv"),
  row.names = TRUE
)

###############################################################################
# 8. Differential-expression heatmap
###############################################################################

top10 <- markers %>%
  group_by(cluster) %>%
  top_n(10, avg_log2FC)

DoHeatmap(
  seurat_obj,
  features = top10$gene
) +
  NoLegend()

###############################################################################
# 9. GO and KEGG over-representation analysis (ORA)
###############################################################################

# Convert gene symbols to Entrez IDs
gene <- bitr(
  count676R$SYMBOL,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

# Enrichment database settings
GO_database <- org.Hs.eg.db
KEGG_database <- "hsa"

# GO enrichment
GO <- enrichGO(
  gene = gene$ENTREZID,
  OrgDb = GO_database,
  keyType = "ENTREZID",
  ont = "ALL",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05,
  readable = TRUE
)

# KEGG enrichment
KEGG <- enrichKEGG(
  gene = gene$ENTREZID,
  organism = KEGG_database,
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05
)

write.csv(
  as.data.frame(GO),
  file.path(output_dir, "ORA_GO.csv"),
  row.names = FALSE
)

write.csv(
  as.data.frame(KEGG),
  file.path(output_dir, "ORA_KEGG.csv"),
  row.names = FALSE
)

# Representative ORA visualization: GO dot plot
dotplot(
  GO,
  split = "ONTOLOGY",
  showCategory = 5
) +
  facet_grid(ONTOLOGY ~ ., scale = "free")

###############################################################################
# 10. GO chord diagrams
###############################################################################

# Extract the top 10 enriched BP terms
GO_BP_result <- GO@result[GO@result$ONTOLOGY == "BP", ]

GOplotIn_BP <- GO_BP_result[
  order(GO_BP_result$p.adjust),
  c("ID", "Description", "p.adjust", "geneID")
][1:10, ]

# Format GOplot input
GOplotIn_BP$geneID <- str_replace_all(GOplotIn_BP$geneID, "/", ",")
names(GOplotIn_BP) <- c("ID", "Term", "adj_pval", "Genes")
GOplotIn_BP$Category <- "BP"

# Prepare gene-level log fold changes
deg <- count676R %>%
  filter(padj < 0.05 & abs(Log2FoldChange) > 0.25)

deg$SYMBOL2 <- mapIds(
  org.Hs.eg.db,
  keys = deg$SYMBOL,
  column = "SYMBOL",
  keytype = "SYMBOL",
  multiVals = "first"
)

deg <- deg[!is.na(deg$SYMBOL2), ]

deg_unique <- deg %>%
  group_by(SYMBOL2) %>%
  slice_max(abs(Log2FoldChange), n = 1) %>%
  ungroup()

genedata <- data.frame(
  ID = deg_unique$SYMBOL2,
  logFC = deg_unique$Log2FoldChange
)

# Generate BP chord-diagram input
circ_BP <- GOplot::circle_dat(GOplotIn_BP, genedata)
chord_BP <- chord_dat(data = circ_BP, genes = genedata)

# Representative GO visualization: biological-process chord diagram
GOChord(
  data = chord_BP,
  title = "GO-Biological Process",
  space = 0.01,
  limit = c(1, 1),
  gene.order = "logFC",
  gene.space = 0.25,
  gene.size = 5,
  lfc.col = c("red", "white", "blue"),
  process.label = 10
)

###############################################################################
# 11. GSEA analysis
###############################################################################

# Differential expression for spatial domain 1
cluster1_markers <- FindMarkers(
  seurat_obj,
  ident.1 = "1",
  test.use = "wilcox",
  logfc.threshold = 0,
  min.pct = 0
)

write.csv(
  cluster1_markers,
  file.path(output_dir, "cluster1_markers.csv"),
  row.names = TRUE
)

# Construct ranked gene list
geneList <- cluster1_markers$avg_log2FC
names(geneList) <- rownames(cluster1_markers)
geneList <- sort(geneList, decreasing = TRUE)

# Convert gene symbols to Entrez IDs
gene_df <- bitr(
  names(geneList),
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

gene_df <- gene_df[!duplicated(gene_df$ENTREZID), ]

geneList <- geneList[gene_df$SYMBOL]
names(geneList) <- gene_df$ENTREZID
geneList <- sort(geneList, decreasing = TRUE)

# GO GSEA
gsea_GO <- gseGO(
  geneList = geneList,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "ALL",
  pvalueCutoff = 0.05
)

write.csv(
  as.data.frame(gsea_GO),
  file.path(output_dir, "GSEA_GO.csv"),
  row.names = FALSE
)

# Representative GO GSEA visualization: dot plot
dotplot(
  gsea_GO,
  showCategory = 20,
  font.size = 10,
  label_format = 100
) +
  scale_size(range = c(2, 6)) +
  scale_color_gradient(low = "blue", high = "red")

# KEGG GSEA
gsea_KEGG <- gseKEGG(
  geneList = geneList,
  organism = "hsa",
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH"
)

write.csv(
  as.data.frame(gsea_KEGG),
  file.path(output_dir, "GSEA_KEGG.csv"),
  row.names = FALSE
)

###############################################################################
# 12. Intersection of ORA and GSEA pathways
###############################################################################

# GO intersection
ora_go_sig <- GO@result %>%
  filter(p.adjust < 0.05)

gsea_go_sig <- gsea_GO@result %>%
  filter(p.adjust < 0.05)

go_intersect <- intersect(
  ora_go_sig$Description,
  gsea_go_sig$Description
)

go_intersect_df <- ora_go_sig %>%
  filter(Description %in% go_intersect)

write.csv(
  go_intersect_df,
  file.path(output_dir, "GO_ORA_GSEA_intersection.csv"),
  row.names = TRUE
)

# KEGG intersection
ora_kegg_sig <- KEGG@result %>%
  filter(p.adjust < 0.05)

gsea_kegg_sig <- gsea_KEGG@result %>%
  filter(p.adjust < 0.05)

kegg_intersect <- intersect(
  ora_kegg_sig$Description,
  gsea_kegg_sig$Description
)

kegg_intersect_df <- ora_kegg_sig %>%
  filter(Description %in% kegg_intersect)

write.csv(
  kegg_intersect_df,
  file.path(output_dir, "KEGG_ORA_GSEA_intersection.csv"),
  row.names = TRUE
)

# GO intersection dot plot
if (nrow(go_intersect_df) > 0) {
  ggplot(
    go_intersect_df[1:min(10, nrow(go_intersect_df)), ],
    aes(x = GeneRatio, y = reorder(Description, GeneRatio))
  ) +
    geom_point(aes(size = Count, color = p.adjust)) +
    scale_color_gradient(low = "red", high = "blue") +
    labs(
      title = "GO Intersection (ORA vs GSEA)",
      x = "Gene Ratio",
      y = ""
    ) +
    theme_bw()
}

###############################################################################
# End of downstream analysis
###############################################################################
