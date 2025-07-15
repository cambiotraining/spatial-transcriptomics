---
title: MERFISH Data Analysis
---

::: {.callout-tip}
#### Learning Objectives

- Load Seurat object containing MERFISH data
- Preprocess the data (normalization, scaling, quality control)
- Perform clustering and visualize results
- Identify marker genes for clusters
:::

## Loading MERFISH Data into Seurat
To load MERFISH data into Seurat, you can use the `LoadVizgen` function. This function reads the MERFISH data from a specified directory and creates a Seurat object. However, in this case, we will load an existing Seurat object that has already been processed and saved.
