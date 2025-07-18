---
title: Visium Example
---

::: {.callout-tip}
#### Learning Objectives

-   Load a Visium data set into Seurat
-   Preprocess the data (normalization, scaling, quality control)
-   Perform clustering and visualize results
-   Identify marker genes for clusters
-   Run BANKSY analysis 
-   Visualize BANKSY results
-   Identify spatially variable features
-   Visualize spatially variable features on spatial maps
:::

## Loading Visium Data into Seurat
To load Visium data into Seurat, you can use the `Load10X_Spatial` function. This function reads the Visium data from a specified directory and creates a Seurat object.
