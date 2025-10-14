---
title: Dimensionality Reduction
---

::: {.callout-tip}
#### Learning Objectives

- Perform PCA on spatial transcriptomics data
- Visualize PCA results using UMAP
- Understand the interpretation of PCA and UMAP results
:::


## Introduction
We will keep working on the sagittal mouse brain dataset we have been using in previous chapters. If you have not finished the preprocessing, you can load the preprocessed data directly:

```r
# Load the precomputed preprocessed Seurat object
visium <- LoadSeuratRds(visium, file = "precomputed/preprocessed_mouse_sagittal.rds")
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
mat <- Seurat::GetAssayData(visium, assay = "SCT", slot = "scale.data")
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
umap <- DimPlot(visium, reduction = "umap", label = TRUE) + ggtitle("UMAP of Spatial Transcriptomics Data")
```

## tSNE
tSNE is another popular dimensionality reduction technique that is often used for visualizing high-dimensional data. In Seurat, you can perform tSNE using the `RunTSNE` function. This function computes the tSNE embedding of the data and stores it in the Seurat object. It's commonly accepted that UMAP performs better than tSNE, but it can still be useful to compare the results of both methods.

```r
# Perform tSNE on the PCA results
visium <- RunTSNE(visium, reduction = "pca", dims = 1:30)
#Visualize the tSNE results
tsne <- DimPlot(visium, reduction = "tsne", label = TRUE) + ggtitle("tSNE of Spatial Transcriptomics Data")

umap + tsne
``` 

These plots show the UMAP and tSNE embeddings of the spatial transcriptomics data, with each point representing a spot in the tissue. They are not particularly informative in this case, as we have not performed any clustering or identified any cell types yet. However, they can be useful for visualizing the overall structure of the data and identifying potential clusters or patterns.

## Conclusion
In this chapter, we have learned how to perform PCA and UMAP on spatial transcriptomics data using Seurat. We have also visualized the results using UMAP and tSNE plots. Dimensionality reduction techniques like PCA and UMAP are essential for analyzing high-dimensional data, as they help to reduce complexity while retaining important information. 
You can use these techniques to explore and visualize your spatial transcriptomics data, identify patterns, and gain insights into the underlying biology.

## Summary
::: {.callout-tip}
#### Key Points
- PCA is a dimensionality reduction technique that helps to reduce the complexity of high-dimensional data while retaining most of the variance.
- UMAP is another dimensionality reduction technique that is particularly well-suited for visualizing high-dimensional data in a low-dimensional space.
- Seurat provides functions for performing PCA, UMAP and tSNA on spatial transcriptomics data.
- Dimensionality reduction techniques like PCA and UMAP are essential for analyzing high-dimensional data, as they help to reduce complexity while retaining important information.
:::