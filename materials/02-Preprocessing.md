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

## Quality Control
Quality control (QC) is crucial to ensure the reliability of the data. In Seurat, you can perform QC by calculating metrics such as the number of detected genes, percentage of mitochondrial genes, and total counts per cell. You can then filter out low-quality spots based on these metrics.

```r
# Calculate QC metrics
visium[["percent.mt"]] <- PercentageFeatureSet(visium, pattern = "^MT-")
visium[["nCount_Spatial"]] <- Matrix::colSums(visium@assay$Spatial) 
``` 
## Identifying and Removing Low-Quality Cells
After calculating the QC metrics, you can visualize them using scatter plots or histograms to identify low-quality cells. You can then remove these cells from the dataset.
```r
# Visualize QC metrics
VlnPlot(visium, features = c("nCount_Spatial", "percent.mt"), ncol = 2)
# Remove low-quality cells
visium <- subset(visium, subset = nCount_Spatial > 500 & percent.mt < 10)
```

## Saving the Preprocessed Data
After preprocessing, you can save the Seurat object to disk for future use.
```r
# Save the preprocessed Seurat object
saveRDS(visium, file = "path/to/your/preprocessed_visium.rds")
``` 
This code provides a basic structure for preprocessing spatial transcriptomics data in Seurat, including scaling, normalization, quality control, and filtering of low-quality cells. You can adapt the parameters based on your specific dataset and analysis requirements.

## Conclusion
Preprocessing spatial transcriptomics data is a critical step in ensuring the quality and reliability of the analysis. By scaling, normalizing, and performing quality control, you can prepare your data for downstream analyses such as clustering, differential expression analysis, and spatial visualization.

## Summary

::: {.callout-tip}
#### Key Points

- Scaling and normalization are essential for comparing spatial transcriptomics data across cells and conditions.
- Quality control helps ensure the reliability of the data by filtering out low-quality cells.
- Seurat provides functions for scaling, normalization, and quality control of spatial transcriptomics data.  
:::
