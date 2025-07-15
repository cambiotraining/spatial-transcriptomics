---
title: Preprocessing Spatial Transcriptomics Data
---

::: {.callout-tip}
#### Learning Objectives

- Scale and normalize spatial transcriptomics data
- Perform quality control on spatial transcriptomics data
- Identify and remove low-quality cells
:::

## Scaling and Normalization  
Scaling and normalization are essential steps in preprocessing spatial transcriptomics data to ensure that the data is comparable across different cells and conditions. In Seurat, this can be done using the `SCTransform` function, which performs normalization and variance stabilization.

```r
# Load the Seurat library
library(Seurat)
# Load the spatial transcriptomics data
visium <- LoadVisium("path/to/your/seurat_object.h5seurat")   
# Perform SCTransform normalization
visium <- SCTransform(visium, assay = "Spatial", verbose = FALSE)
```