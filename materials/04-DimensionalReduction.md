---
title: Dimensional Reduction
---

::: {.callout-tip}
#### Learning Objectives

- Perform PCA on spatial transcriptomics data
- Visualize PCA results using UMAP
- Understand the interpretation of PCA and UMAP results
:::

## Principal Component Analysis (PCA)
Principal Component Analysis (PCA) is a dimensionality reduction technique that helps to reduce the complexity of high-dimensional data while retaining most of the variance. In Seurat, you can perform PCA using the `RunPCA` function. This function computes the principal components of the data and stores them in the Seurat object.