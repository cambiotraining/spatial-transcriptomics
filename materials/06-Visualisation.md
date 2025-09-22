---
title: Visualization of spatial transcriptomics data
---

::: {.callout-tip}
#### Learning Objectives

-   Visualize spatial transcriptomics data using Seurat
:::

We have already seen some basic visualizations of spatial transcriptomics data using Seurat in previous chapters. In this chapter, we will explore their parameters and additional visualization techniques to better understand the spatial organization of gene expression in tissues.

## Basic Visualization with Seurat
Seurat provides several functions for visualizing spatial transcriptomics data. The most commonly used function is `SpatialFeaturePlot`, which allows you to visualize the expression of specific genes across the spatial coordinates of the tissue. Another useful function is `SpatialDimPlot`, which can be used to visualize clusters or other metadata across the spatial coordinates.

```r
# Visualize the expression of specific genes
SpatialFeaturePlot(visium, features = c("Gng4", "Ttr"), ncol  = 2) + ggtitle("Expression of Gng4 and Ttr")
# Visualize clusters across the spatial coordinates
SpatialDimPlot(visium, group.by = "seurat_clusters", label = TRUE) + ggtitle("Spatial Distribution of Clusters")
``` 

You can customize the appearance of these plots using various parameters, such as changing the color scale, adjusting point size, and modifying titles and labels.
```r
# Customize color scale and point size
SpatialFeaturePlot(visium, features = "Gng4", image.alpha = 0.5, cols = paletteer::paletteer_d("khroma::lapaz"))) + ggtitle("Customized Expression of Gng4")

# Visualize co-expression of two features simultaneously
FeaturePlot(visium, features = c("Gng4", "Ttr")), blend = TRUE)
``` 

## Additional Visualization Techniques
Beyond the basic visualization functions provided by Seurat, there are several additional packages that can enhance your ability to visualize and interpret spatial transcriptomics data. Some of these packages include `ggplot2`, `cowplot`, and `patchwork` for advanced plotting capabilities.

```r
# Example of using ggplot2 and patchwork for custom visualizations
library(ggplot2)
library(patchwork)  
# Create a ridge plot of gene expression across spatial coordinates
ridge_plot <- RidgePlot(visium, features = c("Gng4", "Ttr"), ncol = 1, group.by = 'seurat_clusters') + ggtitle("Ridge Plot of Gene Expression")
# Create a violin plot of gene expression across clusters
violin_plot <- VlnPlot(visium, features = c("Gng4", "Ttr"), ncol = 1, group.by = 'seurat_clusters') + ggtitle("Violin Plot of Gene Expression")
# Combine plots using patchwork
ridge_plot + violin_plot

#Create a dot plot of gene expression across clusters
DotPlot(visium, features = c("Gng4", "Ttr"), group.by = 'seurat_clusters') + ggtitle("Dot Plot of Gene Expression")

Create a heatmap of gene expression across clusters
DoHeatmap(visium, features = c("Gng4", "Ttr"), group.by = 'seurat_clusters') + ggtitle("Heatmap of Gene Expression")  
```
