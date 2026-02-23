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
We are going to use a dataset from 10x Genomics, which can be found [here](https://www.10xgenomics.com/resources/datasets/zebrafish-head-3-1-standard-3-0-0). This is a Visium HD 3' dataset of a zebrafish head section.
It comes with pre-segmented data, which we will use int this case to avoid using binned data and having to deconvolute it. Please make sure the relevant directories have been unzipped. Make sure that both the `binned_outputs` and `spatial` directories are uncompressed as well as the `segmented_outputs` directory for pre-segmented high-resolution data. To keep our dataset small enough to run on a standard laptop, we will only load the segmented polygons data.

```r
library(Seurat)
library(ggplot2)
library(patchwork)
libra(dplyr)                                                            
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
# Perform SCTransform normalization. This is going to take a few minutes, feel free to get a coffee.
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
## Subsetting the Object
Since this is a rather large dataset, we will subset it to only include roughly the eye of the zebrafish head for the rest of the analysis. We can do this by subsetting using an interactive plot. 
Please try selecting the innermost part of the eye - including the lens and one of the surrounding circular tissue layers. Unfortunately a larger selection might make it impossible for you to run the detection of spatially variable features later on in a reasonable timeframe.

```r
subset <- InteractiveSpatialPlot(v3d)

#Subset the object to only include the eye region
eye <- subset(v3d, cells = subset)
SpatialFeaturePlot(eye, images = "slice1.polygons", features = "nCount_Spatial.Polygons", plot_segmentations = TRUE, crop = TRUE)
```

This data now has to be preprocessed again after subsetting. 

```r
#Preprocess the subsetted data again
eye <- SCTransform(eye,assay = "Spatial.Polygons", new.assay.name = "Polygon", conserve.memory = TRUE, variable.features.n = 1000, ncells = 2000)
DefaultAssay(eye) <- "Polygon"
eye <- RunPCA(eye, assay = "Polygon", verbose = FALSE)
eye <- RunUMAP(eye, dims = 1:30)
eye <- FindNeighbors(eye, dims = 1:30)
eye <- FindClusters(eye, algorithm = 4, resolution = 0.5, cluster.name = "leidenClusters_05")

#Visualize the clusters on the subsetted data
DimPlot(eye, reduction = "umap", label = TRUE, label.size = 5, group.by = "leidenClusters_05")
SpatialDimPlot(eye, images = "slice1.polygons", plot_segmentations = TRUE, crop = TRUE, group.by = "leidenClusters_05", label = TRUE, label.size = 5)
```

## Identifying Marker Genes
We can identify marker genes for each cluster using the `FindAllMarkers` function. This function identifies genes that are differentially expressed in each cluster compared to all other clusters.

```r
# Find marker genes for each leiden cluster
eye_markers <- FindAllMarkers(eye, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, pvalue.cutoff = 0.05)
```

Additionally we want to see if a BANKSY analysis can help us identify spatial patterns in the data. 

```r
# Run BANKSY analysis
eye_banksy <- RunBanksy(eye,
                    lambda = 0.8, verbose = TRUE,
                    assay = "Polygon", slot = "counts",  features = "variable",
                    k_geom = 50
)
eye_banksy <- RunPCA(eye_banksy,  assay = "BANKSY", reduction.name = "pca.banksy", npcs = 30, features = rownames(eye_banksy))
eye_banksy <- RunUMAP(eye_banksy, dims = 1:30, reduction = "pca.banksy")
eye_banksy <- FindNeighbors(eye_banksy, dims = 1:30, assay = "BANKSY", reduction = "pca.banksy")
eye_banksy <- FindClusters(eye_banksy, algorithm = 4, resolution = 0.3, cluster.name = "banksy_cluster")

# Visualize BANKSY clusters
DimPlot(eye_banksy, reduction = "umap", label = TRUE, label.size = 5, group.by = "banksy_cluster", raster = FALSE)
SpatialDimPlot(eye_banksy, images = "slice1.polygons", plot_segmentations = TRUE, group.by = "banksy_cluster", label = TRUE, label.size = 5)

# Identify marker genes for BANKSY clusters
eye_banksy_markers <- FindAllMarkers(eye_banksy, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, assay = "BANKSY", pvalue.cutoff = 0.05)
```

## Identifying Spatially Variable Features
We can also identify spatially variable features using the `FindSpatiallyVariableFeatures` function. This function identifies genes that show spatially variable expression patterns. We are going toreduce out input to the variables features only and only request the top 20 spatially variable features to be found. We will use the Moran's I method for this analysis. This will still take a while to compute and is a good time to get a coffee or potentially do this step overnight and finish the rest of the analysis the next day.

```r
# Identify spatially variable features using the Moran's I method
eye <- FindSpatiallyVariableFeatures(eye, assay = "Polygon", method = "moransi", features = VariableFeatures(eye), nfeatures = 20)
# Visualize the top spatially variable features
top_spatial_features <- head(SpatiallyVariableFeatures(eye), 5)
SpatialFeaturePlot(eye, images = "slice1.polygons", features = top_spatial_features, plot_segmentations = TRUE, crop = TRUE)
```

You can now explore the identified marker genes and spatially variable features to gain insights into the spatial organization of cell types and gene expression patterns in the zebrafish eye tissue.

## Conclusion
In this example, we have demonstrated how to load, preprocess, and analyze Visium spatial transcriptomics data using Seurat. We have performed subsetting, normalization, dimensionality reduction, clustering, and identified marker genes and spatially variable features. Additionally, we have explored the use of BANKSY for identifying spatial patterns in the data. This workflow can be adapted and extended for other Visium datasets as well.

## Summary

::: {.callout-tip}
#### Key Points
-   Visium data can be loaded into Seurat using the `Load10X_Spatial` function.
-   Quality control is essential to filter out low-quality spots based on metrics like mitochondrial gene percentage.
-   Normalization and scaling can be performed using the `SCTransform` function.
-   Dimensionality reduction (PCA, UMAP) and clustering (Leiden algorithm) help identify distinct cell populations.
-   Marker genes for clusters can be identified using the `FindAllMarkers` function.
-   BANKSY analysis can reveal spatial patterns in gene expression.
-   Spatially variable features can be identified using the `FindSpatiallyVariableFeatures` function.
::: 