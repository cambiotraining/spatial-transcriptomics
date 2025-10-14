---
title: Preprocessing Spatial Transcriptomics Data
---

::: {.callout-tip}
#### Learning Objectives

- Scale and normalize spatial transcriptomics data
- Perform quality control on spatial transcriptomics data
- Identify and remove low-quality cells
:::

## Introduction
Preprocessing spatial transcriptomics data is a crucial step in ensuring the quality and reliability of the analysis. Today we will be working on the Visium mouse dataset we loaded from raw data and saved as a Seurat object in the previous session. The preprocessing steps include scaling, normalization, quality control, and filtering of low-quality cells.
We will keep working with the Seurat object `visium` that we created in the previous chapter.


## Quality Control
Quality control (QC) is crucial to ensure the reliability of the data. In Seurat, you can perform QC by calculating metrics such as the number of detected genes, percentage of mitochondrial genes, and total counts per cell. You can then filter out low-quality spots based on these metrics.

```r
# Calculate QC metrics
visium[["percent.mt"]] <- PercentageFeatureSet(visium, pattern = "^mt-")
visium[["nCount_Spatial"]] <- Matrix::colSums(visium@assays$Spatial) 
``` 
## Identifying and Removing Low-Quality Cells
After calculating the QC metrics, you can visualize them using scatter plots or histograms to identify low-quality cells. You can then remove these cells from the dataset.
```r
# Visualize QC metrics
VlnPlot(visium, features = c("nCount_Spatial", "percent.mt"), ncol = 2)
# Remove low-quality cells and check the difference in metrics
visium <- subset(visium, subset = nCount_Spatial > 1000 & percent.mt < 25)
VlnPlot(visium, features = c("nCount_Spatial", "percent.mt"), ncol = 2)
```

## Scaling and Normalization  
Scaling and normalization are essential steps in preprocessing spatial transcriptomics data to ensure that the data is comparable across different cells and conditions. In Seurat, this can be done using the `SCTransform` function, which performs normalization and variance stabilization. If you are working with a large dataset, consider using the `vars.to.regress` parameter to regress out unwanted sources of variation, such as the percentage of mitochondrial genes. 
We are also reducing the number of cells used for normalization to speed up the process and reduce memory requirements, but you can adjust this based on your dataset size and computational resources. As this step still requires a lot of memory, it is recommended to run it on a machine with sufficient RAM. 

```r
# Perform SCTransform normalization
visium <- SCTransform(visium, assay = "Spatial", verbose = FALSE, ncells = 5000, vars.to.regress = "percent.mt")
```


## Saving the Preprocessed Data
After preprocessing, you can save the Seurat object to disk for future use.
```r
# Save the preprocessed Seurat object
SaveSeuratRds(visium, file = "data/my_preprocessed_mouse_sagittal.rds")
``` 

## Loading the Preprocessed Data
We can load the preprocessed object to avoid waiting for the normalization step to complete, as it can take a significant amount of time and memory. 

```r
# Load the precomputed preprocessed Seurat object
visium <- LoadSeuratRds(visium, file = "precomputed/preprocessed_mouse_sagittal.rds")
``` 

This code provides a basic structure for preprocessing spatial transcriptomics data in Seurat, including scaling, normalization, quality control, and filtering of low-quality cells. You can adapt the parameters based on your specific dataset and analysis requirements.

## Conclusion
Preprocessing spatial transcriptomics data is a critical step in ensuring the quality and reliability of the analysis. By scaling, normalizing, and performing quality control, you can prepare your data for downstream analyses such as clustering, differential expression analysis, and spatial visualization.

## Summary

::: {.callout-tip}
#### Key Points

- Quality control helps ensure the reliability of the data by filtering out low-quality cells.
- Quality control metrics include the number of detected genes, percentage of mitochondrial genes, and total counts per cell.
- Scaling and normalization are essential for comparing spatial transcriptomics data across cells and conditions.
- Seurat provides functions for scaling, normalization, and quality control of spatial transcriptomics data.  
:::
