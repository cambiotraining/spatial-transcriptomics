---
title: Loading Data into Seurat
---

::: {.callout-tip}
#### Learning Objectives

- Load multiple data sets into Seurat
- Understand the structure of spatial data in Seurat
- Familiarize with the Seurat object and its components
- Explore the data using Seurat functions
- Load an existing Seurat object
- Save a Seurat object to disk
:::

## Loading Visium Data into Seurat
To load Visium data into Seurat, you can use the `Load10X_Spatial` function. This function reads the Visium data from a specified directory and creates a Seurat object.
```r
# Load the Seurat library
library(Seurat)
# We will also use the `paletteer` package for color palettes
library(paletteer)
# Load the spatial transcriptomics data
visium <- Load10X_Spatial("data/mouse_visium")
```
This will create a Seurat object containing the spatial transcriptomics data, which includes gene expression data, spatial coordinates, an image of the tissue slide and other metadata.

There are other functions available in Seurat for loading different types of spatial transcriptomics data, such as `LoadXenium` for 10X Genomics Xenium data or `LoadVizgen` for loading ViZgen MerFISH data.

## Understanding the Seurat Object
The Seurat object is a complex data structure that contains multiple components, including:

- `assays`: Contains the gene expression data and other assays.
- `meta.data`: Contains metadata for each cell or spot, such as quality control metrics and experimental conditions.
- `reductions`: Contains dimensionality reduction results, such as PCA or UMAP.
- `images`: Contains images associated with the spatial transcriptomics data.
- `graphs`: Contains graph-based representations of the data, such as nearest neighbor graphs.
- `tools`: Contains various tools and methods used for analysis.

You can explore the Seurat object using various functions, such as `str()` and the RStudio-specific function `View()`, to understand its structure and contents.

```r
# Explore the Seurat object
str(visium)
View(visium)
```

## Exploring the Data
You can use `ImageDimPlot` to visualize the spatial distribution of cells or spots in the tissue slide. This function allows you to plot the spatial coordinates of the cells or spots on the tissue image. Using `SpatialFeaturePlot`, you can visualize specific features, such as the number of transcripts per spot or the expression of specific genes.

```r
# Visualize the transcript locations on the tissue slide
ImageDimPlot(visium)
# Visualize the number of transcripts per spot
SpatialFeaturePlot(visium, features = "nCount_Spatial.016um")
```   

## Loading an Existing Seurat Object
If you have an existing Seurat object saved to disk, you can load it using the `readRDS()` function. This is useful for continuing analysis from a previous session or sharing results with collaborators. We are loading a MERFISH dataset as an example. It already includes some annotations like tissue types and cluster identities.
```r
# Load an existing Seurat object
merfish <- readRDS("data/human_heart_merfish/overall_merfish.rds")
#Visualize the structure of the heart in the MERFISH data
DimPlot(merfish, reduction = "spatial", group.by = "communities", cols = paletteer_d("ggthemes::Tableau_20")
``` 

## Saving a Seurat Object
You can save a Seurat object to disk using the `saveRDS()` function. This allows you to preserve your analysis results and share them with others.
```r
# Save the Seurat object to disk
saveRDS(visium, file = "data/mouse_visium.rds")
```

## Loading Xenium Data into Seurat
If you are working with 10X Genomics Xenium data, you can use the `LoadXenium` function to load the data into a Seurat object. This function is specifically designed to handle the unique structure of Xenium data, which includes spatial transcriptomics data with high-resolution spatial coordinates and associated images.

We will use another dataset from the 10X examples. The data is available on the 10X website and has been downloaded for you in the course materials at `data/human_melanoma_xenium`.

```r
path <- "data/human_melanoma_xenium"
xenium <- LoadXenium(path, fov = "fov")
```
The FOV (Field of View) parameter allows you to specify which field of view to load if the dataset contains multiple fields. In this case, we are loading the default field of view.
 
To visualize the Xenium data we will use the `ImageDimPlot` function to plot the spatial distribution of cells or spots in the tissue slide. We can also visualize the number of features measuresd in each location using `ImageFeaturePlot`.

```r
ImageDimPlot(xenium)
# Visualize the transcript density on the tissue slide  
ImageFeaturePlot(xenium, features = "nFeature_Xenium", cols = paletteer_c("grDevices::Blue-Red 3", 150))
```


## Conclusion
Loading data into Seurat is a straightforward process that allows you to work with spatial transcriptomics data effectively. By understanding the structure of the Seurat object and using its functions, you can explore and analyze your data efficiently. It's always a good practice to save your Seurat objects after significant analysis steps to avoid losing your work and to facilitate collaboration with others.

## Summary  
::: {.callout-tip}
#### Key Points   
- Use `Load10X_Spatial` to load Visium data into Seurat.
- Understand the structure of the Seurat object, including assays, metadata, reductions, images, graphs, and tools.
- Explore the data visually using functions like `ImageDimPlot` and `SpatialFeaturePlot`.
- Load existing Seurat objects with `readRDS` and save them with `saveRDS`.
- Use `LoadXenium` to load 10X Genomics Xenium data into Seurat.
- Use `ImageDimPlot` and `ImageFeaturePlot` to visualize Xenium data.
- Familiarize yourself with the Seurat object to effectively analyze spatial transcriptomics data.
:::
