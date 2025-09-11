---
title: Deconvolution of Spatial Transcriptomics Data
---

::: {.callout-tip}
#### Learning Objectives

- Understand the concept of deconvolution in spatial transcriptomics
- Apply deconvolution methods to spatial transcriptomics data
- Use the `RCTD` package for deconvolution analysis
- Visualize and interpret deconvolution results
:::

## Deconvolution of Spatial Transcriptomics Data
Deconvolution is a computational technique used to estimate the cellular composition of complex tissues from bulk gene expression data. In the context of spatial transcriptomics, deconvolution can help identify the proportions of different cell types within spatially resolved samples.
We will be using the `RCTD` package, which provides functions for deconvolution of spatial transcriptomics data.

```r
# Load necessary libraries
library(spacexR)
library(Seurat)
library(ggplot2)
library(dplyr)
library(pheatmap)
library(SingleCellExperiment)
library(SummarizedExperiment)
library(SpatialExperiment)
library(BiocParallel) # for parallel processing 
```

### Load Data
First, we need to load the spatial transcriptomics data and the single-cell reference data. For this example, we will use a preprocessed Seurat object for the spatial data and a SingleCellExperiment object for the single-cell reference data. 

```r
# Load the spatial transcriptomics data
visium <- LoadSeuratRds("precomputed/preprocessed_mouse_sagittal.rds")   
# Convert Seurat object to SpatialExperiment object
spe <- as(visium, "SpatialExperiment")
# Load the single-cell reference data
sce <- readRDS("precomputed/human_lung_sc_reference.rds")
```


