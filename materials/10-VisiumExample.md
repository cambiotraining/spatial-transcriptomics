---
title: Visium Example
---

::: {.callout-tip}
#### Learning Objectives

-   Load a Visium data set into Seurat
-   Preprocess the data (normalization, scaling, quality control)
-   Perform clustering and visualize results
-   Identify marker genes for clusters
-   Run BANKSY analysis 
-   Visualize BANKSY results
-   Identify spatially variable features
-   Visualize spatially variable features on spatial maps
:::

## Loading Visium Data into Seurat
To load Visium data into Seurat, you can use the `Load10X_Spatial` function. This function reads the Visium data from a specified directory and creates a Seurat object.
In this example we are loading a Visium HD 3' dataset. Seurat supports loading 3' datasets through a developmental feature. Please note that this is still experimental and may not work for all datasets. We need to load additional libraries to handle this type of data. 
We are going to use a zebrafish dataset from 10x Genomics, which can be found [here](https://www.10xgenomics.com/resources/datasets/zebrafish-head-3-1-standard-3-0-0). It comes with pre-segmented data, which we will use int this case to avoid using binned data and having to deconvolute it. Please make sure the relevant directories have been unzipped. Make sure that both the `binned_outputs` and `spatial` directories are uncompressed as well as the `segmented_outputs` directory for pre-segmented high-resolution data. To keep our dataset small enough to run on a standard laptop, we will only loadthe segmented polygons data.

```r
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(sf)
library(Matrix)
library(SeuratWrappers)
library(Banksy)
v3d <- Load10X_Spatial(data.dir = "data/zebrafish_head", slice = "slice1", bin.size = c("polygons"))
```

## Filtering and Quality Control
After loading the data, we can perform some basic quality control (QC) to filter out low-quality spots. We will calculate the percentage of mitochondrial genes and filter out spots with high mitochondrial content, as well as spots with very low total counts. because we are working with polygons, we will use the `nCount_Spatial.Polygons` metadata field instead of `nCount_Spatial` and we are adjusting the filtering thresholds accordingly as well.
```r
# Calculate percentage of mitochondrial genes
v3d[["percent.mt"]] <- PercentageFeatureSet(v3d, pattern = "^mt-")
# Visualize QC metrics
VlnPlot(v3d, features = c("nCount_Spatial.Polygons", "percent.mt"), ncol = 2)
# Remove low-quality cells and check the difference in metrics
v3d <- subset(v3d, subset = nCount_Spatial.Polygons > 5 & percent.mt < 10)
# Visualize QC metrics again to check improvement
VlnPlot(v3d, features = c("nCount_Spatial.Polygons", "percent.mt"), ncol = 2)
```

## Normalization and Scaling
Next, we will normalize and scale the data using the `SCTransform` function. This function performs normalization and variance stabilization. We are also reducing the number of cells used for normalization to speed up the process and reduce memory requirements and adding a parameter to conserve memory. 

```r
# Perform SCTransform normalization
v3d <- SCTransform(v3d,assay = "Spatial.Polygons", new.assay.name = "Polygon", conserve.memory = TRUE, variable.features.n = 1000, ncells = 2000)
# Set the default assay to the newly created Polygon assay
DefaultAssay(v3d) <- "Polygon"
```

## Dimensionality Reduction and Clustering
After normalization, we can perform dimensionality reduction using PCA and then cluster the spots using the Leiden algorithm. We will visualize the clusters using UMAP.

```r
# Perform PCA
v3d <- RunPCA(v3d, assay = "Polygon", verbose = FALSE)
v3d <- RunUMAP(v3d, dims = 1:30)
v3d <- FindNeighbors(v3d, dims = 1:30)
v3d <- FindClusters(v3d, algorithm = 4, resolution = 0.3, cluster.name = "leidenClusters_03")
# Visualize UMAP with clusters
DimPlot(v3d, reduction = "umap", group.by = "leidenClusters_03", label = TRUE) + ggtitle("Leiden Clusters (res=0.3)")
SpatialDimPlot(v3d, group.by = "leidenClusters_03", label = TRUE) + ggtitle("Leiden Clusters (res=0.3)")
``` 

## Identifying Marker Genes
We can identify marker genes for each cluster using the `FindAllMarkers` function. 

```r
# Identify marker genes for each cluster
markers <- FindAllMarkers(v3d, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, assay = "Polygon")
# View top markers for each cluster
top10 <- markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
# We will look at cluster 0 as an example
top10[top10$cluster==0,]
```

## Running BANKSY Analysis
We will now run BANKSY on our dataset and visualize the results. BANKSY objects need to be pre-processed again, but without the scaling step, as BANKSY output is already scaled. 

```r
banksy <- RunBanksy(v3d,
                    lambda = 0.8, verbose = TRUE,
                    assay = "Polygon", slot = "counts",  features = "variable",
                    k_geom = 50
)
#We won't scale/normalise here, because that removes the effect of Banksy being already scaled, so we directly run PCA
banksy <- RunPCA(banksy,  assay = "BANKSY", reduction.name = "pca.banksy", npcs = 30, features = rownames(banksy))
banksy <- RunUMAP(banksy, dims = 1:30, reduction = "pca.banksy")
banksy <- FindNeighbors(banksy, dims = 1:30, assay = "BANKSY", reduction = "pca.banksy")
banksy <- FindClusters(banksy, algorithm = 4, resolution = 0.3, assay = "BANKSY", cluster.name = "banksy_cluster")

# Visualize BANKSY clusters
DimPlot(banksy, reduction = "umap", label = TRUE, label.size = 5, group.by = "banksy_cluster", raster = FALSE)
SpatialDimPlot(banksy, images = "slice1.polygons", plot_segmentations = TRUE, group.by = "banksy_cluster", label = TRUE, label.size = 5)

# Identify marker genes for BANKSY clusters
banksy_markers <- FindAllMarkers(banksy, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, assay = "BANKSY")
# View top markers for each BANKSY cluster
top10_banksy <- banksy_markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
# We will look at cluster 1 as an example
top10_banksy[top10_banksy$cluster==1,]
```

## Identifying Spatially Variable Features
We can identify spatially variable features using the `FindSpatiallyVariableFeatures` function. This function identifies genes that show spatial patterns in their expression.  

```r
# Identify spatially variable features
spatial_features <- FindSpatiallyVariableFeatures(v3d, assay = "Polygon", selection.method = "markvariogram", nfeatures = 100)
# View top spatially variable features
head(spatial_features, 10)
```

