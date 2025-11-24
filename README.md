# singlecell-pbmc3k-demo

This repository contains a small single-cell RNA-seq analysis project using the
10x Genomics PBMC 3k dataset and the Seurat R package.

## Data

- **Dataset:** PBMC 3k (peripheral blood mononuclear cells from a healthy donor)
- **Source:** Public 10x Genomics dataset
- **Files used:** Filtered gene-barcode matrices downloaded from the 10x support site
- **License:** The PBMC 3k dataset is released by 10x Genomics under a
  Creative Commons Attribution 4.0 International (CC BY 4.0) license.


The raw data are **not** stored in this repository. To reproduce the analysis,
download the PBMC 3k filtered matrix from 10x Genomics and place it under
`data/` as described below.

---

## Requirements

- **R:** ≥ 4.2 (tested on R 4.4.0)
- **R packages:**
  - `Seurat`
  - `dplyr`
  - `ggplot2`
  - `patchwork`
  - `Matrix`

You can install the required packages with:


install.packages(c("Seurat", "dplyr", "ggplot2", "patchwork"))

## Analysis overview

All analysis code is in `R/pbmc3k_analysis.R` and follows a standard
Seurat-based scRNA-seq workflow:

1. **Data import**
   - Read 10x filtered gene-barcode matrices with `Read10X()`.
   - Create a `Seurat` object with `CreateSeuratObject()`.

2. **Quality control**
   - Compute basic QC metrics: `nFeature_RNA`, `nCount_RNA`, and `percent.mt`.
   - Generate violin plots and scatter plots to inspect QC metrics.
   - Filter low-quality cells based on feature counts and mitochondrial percentage.

3. **Preprocessing**
   - Normalize counts (`NormalizeData()`).
   - Identify highly variable genes (`FindVariableFeatures()`).
   - Scale the data (`ScaleData()`).

4. **Dimensionality reduction and clustering**
   - Run PCA on variable features (`RunPCA()`).
   - Construct a nearest-neighbour graph (`FindNeighbors()`).
   - Perform graph-based clustering (`FindClusters()`).
   - Compute a 2D embedding using UMAP (`RunUMAP()`).

5. **Marker discovery and interpretation**
   - Identify cluster-specific marker genes with `FindAllMarkers()`.
   - Export the marker table as `pbmc3k_markers.csv`.
   - Visualize expression of selected immune marker genes with `FeaturePlot()`.

---

## Repository structure

- `R/pbmc3k_analysis.R`  
  Main analysis script implementing the workflow described above.

- `figures/`  
  Exported PNG images from the analysis, for example:
  - `qc_vlnplot.png`: QC violin plots for nFeature_RNA, nCount_RNA, percent.mt.
  - `feature_scatter.png`: Scatter plots of QC metrics.
  - `umap_clusters.png`: UMAP of cells coloured by cluster.
  - `featureplots_markers.png`: Feature plots for selected marker genes.

- `pbmc3k_markers.csv`  
  Table of marker genes for each cluster (log fold change, p-values, etc.).

- `data/` *(not tracked in this repository)*  
  Directory where the 10x PBMC 3k data should be placed locally.

---

## How to run the analysis

1. **Clone this repository** and open it as an RStudio project (optional).

2. **Download the PBMC 3k filtered matrix** from 10x Genomics and extract it under:

   data/filtered_gene_bc_matrices/hg19/

3. Install the required R packages (see the Requirements section).

4. Run the analysis from the project root in R or RStudio:

  source("R/pbmc3k_analysis.R")

This will:

Load the PBMC 3k data.

Perform QC, preprocessing, clustering, and marker analysis.

Save figures under figures/ and the marker table as pbmc3k_markers.csv.
   
## Acknowledgements
The analysis steps in this repository are inspired by the official
Seurat PBMC 3k guided clustering tutorial and adapted for my own
learning and use in job applications.
