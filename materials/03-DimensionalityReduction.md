---
title: Dimensionality Reduction
---

::: {.callout-tip}
#### Learning Objectives

- Perform PCA on spatial transcriptomics data
- Visualize PCA results using UMAP
- Understand the interpretation of PCA and UMAP results
:::

```r
# Load the Seurat library
library(Seurat)
#load the precomputed preprocessed Seurat object
visium <- LoadSeuratRds("precomputed/preprocessed_human_lung_visium.rds")

```

## Principal Component Analysis (PCA)
Principal Component Analysis (PCA) is a dimensionality reduction technique that helps to reduce the complexity of high-dimensional data while retaining most of the variance. In Seurat, you can perform PCA using the `RunPCA` function. This function computes the principal components of the data and stores them in the Seurat object.

```r
# Perform PCA on the preprocessed data
visium <- RunPCA(visium, assay = "SCT", npcs = 30, verbose = FALSE)
#Check elbow plot to determine the number of significant PCs
ElbowPlot(visium, ndims = 30, reduction = 'pca')
#For good measure, we will check if it is better to use 50 PCs
visium <- RunPCA(visium, assay = "SCT", npcs = 50, verbose = FALSE,  reduction.name = 'pca50')
ElbowPlot(visium, ndims = 50, reduction = 'pca50')
``` 

If we have selected a good number of PCs, we should see a clear elbow in the plot, indicating that the first few PCs capture most of the variance in the data. To be sure we can also check the variance explained by the selected PCs:

```r  
mat <- Seurat::GetAssayData(visium, assay = "RNA", slot = "scale.data")
pca <- visium[["pca"]]

# Get the total variance:
total_variance <- sum(matrixStats::rowVars(mat))

eigValues = (pca@stdev)^2  ## EigenValues
varExplained = eigValues / total_variance
plot(varExplained, xlab = "PCs", ylab = "Variance Explained", main = "Variance Explained by PCs")
``` 

## Uniform Manifold Approximation and Projection (UMAP)
UMAP is another dimensionality reduction technique that is particularly well-suited for visualizing high-dimensional data in a low-dimensional space. In Seurat, you can perform UMAP using the `RunUMAP` function. This function computes the UMAP embedding of the data and stores it in the Seurat object. 

```r
# Perform UMAP on the PCA results
visium <- RunUMAP(visium, reduction = "pca", dims = 1:30)
#Visualize the UMAP results
DimPlot(visium, reduction = "umap", label = TRUE) + ggtitle("UMAP of Spatial Transcriptomics Data")
```