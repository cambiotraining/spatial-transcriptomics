---
title: Loading Data into Seurat
---

::: {.callout-tip}
#### Learning Objectives

- Load a Visum data set into Seurat
- Understand the structure of Visium data in Seurat
- Familiarize with the Seurat object and its components
- Explore the data using Seurat functions
- Load an existing Seurat object (MERFISH data)
- Save a Seurat object to disk
:::

## Loading Visium Data into Seurat
To load Visium data into Seurat, you can use the `Load10X_Spatial` function. This function reads the Visium data from a specified directory and creates a Seurat object.
```r
# Load the Seurat library
library(Seurat)
# Load the spatial transcriptomics data
visium <- Load10X_Spatial("path/to/your/visium/data/directory")
```
This will create a Seurat object containing the spatial transcriptomics data, which includes gene expression data, spatial coordinates, an image of the tissue slide and other metadata.

There are other functions available in Seurat for loading different types of spatial transcriptomics data, such as `LoadXenium` for 10X Genomics Xenium data or `LoadVizgen` for loading ViZgen MerFISH data.