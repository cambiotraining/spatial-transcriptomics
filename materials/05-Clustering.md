---
title: Clustering of Spatial Transcriptomics Data
---

::: {.callout-tip}
#### Learning Objectives

- Cluster spatial transcriptomics data using Seurat
- Visualize clusters using UMAP and spatial plots
- Identify marker genes for clusters
- Understand the interpretation of clustering results
:::

## Clustering Spatial Transcriptomics Data
Clustering is a crucial step in spatial transcriptomics analysis, allowing us to group cells based on their gene expression profiles. In Seurat, clustering can be performed using the `FindClusters` function after dimensionality reduction techniques like PCA or UMAP.
