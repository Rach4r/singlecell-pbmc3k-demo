library(Seurat)
library(tidyverse)
library(patchwork)

# Make sure you're in the project root
dir.create("data", showWarnings = FALSE)

download.file(
  url      = "https://cf.10xgenomics.com/samples/cell/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz",
  destfile = "data/pbmc3k_filtered_gene_bc_matrices.tar.gz"
)

untar("data/pbmc3k_filtered_gene_bc_matrices.tar.gz", exdir = "data")

# 1. Load data ------------------------------------------------------------

data_dir <- "data/filtered_gene_bc_matrices/hg19/"
pbmc.data <- Read10X(data.dir = data_dir)

pbmc <- CreateSeuratObject(
  counts = pbmc.data,
  project = "pbmc3k",
  min.cells = 3,
  min.features = 200
)
pbmc

# 2. QC metrics -----------------------------------------------------------

pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

# 1. Check where you are
getwd()

# 2. Create the 'figures' folder if it doesn't exist
if (!dir.exists("figures")) {
  dir.create("figures", recursive = TRUE)
}


# Basic QC plots (save to figures/)
png("figures/qc_vlnplot.png", width = 1600, height = 800, res = 150)
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()

png("figures/feature_scatter.png", width = 1600, height = 800, res = 150)
FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
  FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
dev.off()

# Filter cells (tweak thresholds if needed)
pbmc <- subset(
  pbmc,
  subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5
)

# 3. Normalization, variable features, scaling ----------------------------

pbmc <- NormalizeData(pbmc)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
pbmc <- ScaleData(pbmc, features = VariableFeatures(pbmc))

# 4. PCA and clustering ---------------------------------------------------

pbmc <- RunPCA(pbmc, features = VariableFeatures(pbmc))
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)
pbmc <- RunUMAP(pbmc, dims = 1:10)

# 5. Basic visualizations -------------------------------------------------

png("figures/umap_clusters.png", width = 1600, height = 1200, res = 150)
DimPlot(pbmc, reduction = "umap", label = TRUE) + NoLegend()
dev.off()

# 6. Marker genes ---------------------------------------------------------

pbmc.markers <- FindAllMarkers(
  pbmc,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)

# Save markers as CSV
write.csv(pbmc.markers, file = "pbmc3k_markers.csv", row.names = FALSE)

# Plot a few canonical markers (example genes)
genes_to_plot <- c("IL7R", "CCR7", "LYZ", "MS4A1", "CD8A", "GNLY")

png("figures/featureplots_markers.png", width = 2000, height = 1600, res = 150)
FeaturePlot(pbmc, features = genes_to_plot)
dev.off()


