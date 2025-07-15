---
title: Spatially Variable Features
---

::: {.callout-tip}
#### Learning Objectives

-  Identify spatially variable features in spatial transcriptomics data
-  Visualize spatially variable features on spatial maps
-  Understand the biological significance of spatially variable features
:::

## Identifying Spatially Variable Features
Spatially variable features are genes or transcripts that show significant variation in expression across different spatial locations within a tissue sample. Identifying these features can provide insights into the spatial organization of gene expression and the underlying biological processes.
In Seurat, you can use the `FindSpatiallyVariableFeatures` function to identify spatially variable features in your spatial transcriptomics data. This function performs statistical tests to identify genes that exhibit significant spatial variability.