---
title: Human Melanoma Xenium Analysis
---

::: {.callout-tip}
#### Learning Objectives

- Go through standard analysis steps for Xenium data
- Understand how to visualize and interpret Xenium data
- Basic interpretation of markers and spatial patterns
:::

## Libraries
To run the code in this section, you will need to load the following packages.

```r
library(Seurat)
library(SeuratWrappers)
library(Banksy)
library(dplyr)
library(ggplot2)
library(patchwork)
library(pheatmap)
library(clusterProfiler)
library(org.Hs.eg.db) 
```

## Loading the data
We will use a dataset from the 10X example datasets. The data is available on the 10X website and has been downloaded for you in the course materials at `data/human_melanoma_xenium`.

```r
path <- "data/human_melanoma_xenium"
xenium <- LoadXenium(path, fov = "fov")
```

## Preprocessing and Dimensionality Reduction
Perform SCTransform on the Xenium data to normalize and scale the data. Then, run PCA for dimensionality reduction to 30 principal components using all measured features. Finally create a UMAP reduction for future visualization.

```r
xenium <- SCTransform(xenium)
xenium <- RunPCA(xenium, npcs = 30, features = rownames(xenium))
xenium <- RunUMAP(xenium, dims = 1:30)
```

## Clustering
Before we run clustering, we need to find the nearest neighbors of each cell in the dataset. This is done using the `FindNeighbors` function. We will use the first 30 principal components from the PCA step for this.
We will use the Louvain algorithm to cluster the cells in the Xenium dataset. The resolution parameter can be adjusted to control the granularity of the clustering. Here, we will use a resolution of 0.5, which is a common starting point.

```r
xenium <- FindNeighbors(xenium, dims = 1:30, reduction = "pca")
xenium <- FindClusters(xenium, resolution = 0.5)
```


## Visualizing Clusters
We can visualize the clusters using the UMAP reduction we created earlier. The `DimPlot` function allows us to plot the UMAP with the clusters colored by their cluster identity. For spatial data, we can also visualize the clusters on the spatial coordinates of the cells, by using the `ImageDimPlot` function.
You can also show a single cluster on the spatial image by using the `cells` argument in the `ImageDimPlot` function and selecting the cells from a specific cluster using the `WhichCells` function to select only cells from cluster 0.

```r
DimPlot(xenium, reduction = "umap", group.by = "seurat_clusters", label = TRUE) +
  labs(title = "UMAP of spatial clusters")
ImageDimPlot(xenium, group.by = "seurat_clusters") +
  labs(title = "Spatial clusters on image data")

# Show only cluster 0 on the spatial image
ImageDimPlot(xenium, cells = WhichCells(xenium, idents = 0)) +
  labs(title = "Spatial cluster 0 on image data")
```

## Finding Marker Genes
To find marker genes for each cluster, we can use the `FindAllMarkers` function.
This function will return a data frame with the marker genes for each cluster, along with their average expression, p-value, and adjusted p-value. We can then visualize the top marker genes for each cluster using the `FeaturePlot` function. This will create a scatter plot of the expression of the marker genes on the UMAP reduction, allowing us to see how the marker genes are distributed across the clusters.

```r
markers <- FindAllMarkers(xenium, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, pvalue.cutoff = 0.05)
top_markers <- markers %>%
  group_by(cluster) %>%
  top_n(1, avg_log2FC) %>%
  ungroup()
FeaturePlot(xenium, features = top_markers$gene, ncol = 5) +
  plot_annotation(title = "Top marker genes for each cluster")
ImageFeaturePlot(xenium, features = top_markers$gene) +
  plot_annotation(title = "Top marker genes for each cluster on image data")
```

## BANKSY analysis
We can use the Banksy package to perform a more detailed analysis of the spatial patterns in the Xenium data. The Banksy package provides a set of functions for spatial analysis, including spatial clustering, spatial enrichment, and spatial correlation analysis. Here, we will use the `RunBanksy` function provided by the `SeuratWrappers` library to perform spatially informed clustering on the Xenium data. This function will identify spatial clusters in the data based on the spatial coordinates of the cells and the expression of the marker genes. We will use the `lambda` parameter to control the strength of the spatial clustering, and the `k_geom` parameter to control the number of nearest neighbors used in the clustering. 

```r
banksy <- RunBanksy(xenium,
                    lambda = 0.8, verbose = TRUE,
                    assay = "Xenium", slot = "counts",  features = "variable",
                    k_geom = 50
)
``` 

The `RunBanksy` function will return a new `Seurat` object with the Banksy results stored in the `BANKSY` assay. We then need to create a new PCA reduction for the Banksy results, which will be used for visualization and further analysis. We can use the `RunPCA` function to create a PCA reduction for the Banksy results, using the `BANKSY` assay and the `counts` slot.
We will also run clustering on the Banksy results using the `FindNeighbors` and `FindClusters` functions, similar to how we did it for the original data. 

```r
banksy <- RunPCA(banksy, assay = "BANKSY", reduction.name = "pca.banksy", npcs = 30, features = rownames(banksy))
banksy <- RunUMAP(banksy, dims = 1:30, reduction = "pca.banksy")
banksy <- FindNeighbors(banksy, dims = 1:30, assay = "BANKSY", reduction = "pca.banksy")
banksy <- FindClusters(banksy, resolution = 0.3, assay = "BANKSY", cluster.name = "banksy_cluster")
```


After running the Banksy analysis, we can visualize the spatial clusters identified by Banksy on the spatial image. The `ImageDimPlot` function allows us to plot the spatial clusters on the image data.

```r
ImageDimPlot(banksy, group.by = "banksy_cluster") +
  labs(title = "Spatial clusters identified by Banksy on image data")
``` 

## Investigating BANKSY results
We can investigate the BANKSY results by looking at the marker genes for each Banksy cluster. We can use the `FindAllMarkers` function to find the marker genes for each Banksy cluster, similar to how we did it for the original clusters. We can then visualize the top marker genes for each Banksy cluster using the `FeaturePlot` function.

```r
banksy_markers <- FindAllMarkers(banksy, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, assay = "BANKSY", pvalue.cutoff = 0.05)
banksy_top_markers <- banksy_markers %>%
  group_by(cluster) %>%
  top_n(1, avg_log2FC) %>%
  ungroup()
FeaturePlot(banksy, features = banksy_top_markers$gene, ncol = 5) +
  plot_annotation(title = "Top marker genes for each Banksy cluster")
ImageFeaturePlot(banksy, features = banksy_top_markers$gene) +
  plot_annotation(title = "Top marker genes for each Banksy cluster on image data")
``` 

BANKSY adds m0 or m1 to the gene names to indicate either mean neighbor expression (m0) or AGF (Azimuthal Gabor filter) based expression (m1).This is why some of the marker genes have the suffix "m0" added here.

## Comparing BANKSY and Seurat Clusters
We can compare the clusters identified by Banksy with the original clusters identified by Seurat. We can use the `table` function to create a contingency table showing the overlap between the Banksy clusters and the Seurat clusters. This will allow us to see how well the Banksy clusters align with the original clusters by counting the number of cells in each Banksy cluster that belong to each Seurat cluster.
We can also visualize the overlap between the Banksy clusters and the Seurat clusters using a heatmap. The `pheatmap` function can be used to create a heatmap of the contingency table, showing the number of cells in each Banksy cluster that belong to each Seurat cluster.

```r
contingency_table <- table(banksy$banksy_cluster, xenium$seurat_clusters)
pheatmap(contingency_table,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         display_numbers = TRUE,
         fontsize_number = 10,
         main = "Overlap between Banksy clusters and Seurat clusters",
         color = colorRampPalette(c("white", "blue"))(100))
```

Seurat clusters have been calculated using the first 30 principal components using Louvain clustering on gene expression, while Banksy clusters are spatially informed including the mean gene expression of spatial cell neighborhoods. 
This means that the Banksy clusters are more likely to capture spatial patterns in the data, while the Seurat clusters are more likely to capture gene expression patterns. The overlap between the two clustering methods can provide insights into how well the spatial patterns align with the gene expression patterns.

## Gene Set Enrichment Analysis
We can also perform gene set enrichment analysis on the Banksy and Seurat clusters to identify enriched pathways or biological processes. To do that we will use the library `clusterProfiler` to perform the enrichment analysis. We will use the `enrichGO` function to perform Gene Ontology (GO) enrichment analysis on the marker genes for each cluster.
We will use the `org.Hs.eg.db` package to map the gene symbols to Entrez IDs, which is required for the enrichment analysis.

```r
# Convert gene symbols for one cluster to Entrez IDs (let's do cluster 0 for Seurat and Banksy as an example)
seurat_cluster0 <- markers[markers$cluster == 0,]$gene
# Convert gene symbols to Entrez IDs
entrez_seurat_cluster0 <- bitr(seurat_cluster0, fromType = "SYMBOL", 
                   toType = "ENTREZID", 
                   OrgDb = org.Hs.eg.db)
# Perform GO enrichment analysis
go_enrichment_seurat_cluster0 <- enrichGO(gene = entrez_seurat_cluster0$ENTREZID,
                          OrgDb = org.Hs.eg.db,
                          keyType = "ENTREZID",
                          ont = "BP", # Biological Process
                          pAdjustMethod = "BH",
                          qvalueCutoff = 0.05,
                          readable = TRUE)
# View the results
head(go_enrichment_seurat_cluster0)
# Plotting the GO enrichment results for Seurat cluster 0
barplot(go_enrichment_seurat_cluster0, showCategory = 10) +
  labs(title = "GO Enrichment for Seurat Cluster 0")

# For Banksy cluster 0, we will do the same
banksy_cluster0 <- banksy_markers[banksy_markers$cluster == 0,]$gene
#for the Banksy cluster, we need to remove the ".m0" suffix from the gene names
banksy_cluster0 <- gsub(".m0$", "", banksy_cluster0)
entrez_banksy_cluster0 <- bitr(banksy_cluster0, fromType = "SYMBOL", 
                   toType = "ENTREZID", 
                   OrgDb = org.Hs.eg.db)
# Perform GO enrichment analysis
go_enrichment_banksy_cluster0 <- enrichGO(gene = entrez_banksy_cluster0$ENTREZID,
                          OrgDb = org.Hs.eg.db,
                          keyType = "ENTREZID",
                          ont = "BP", # Biological Process
                          pAdjustMethod = "BH",
                          qvalueCutoff = 0.05,
                          readable = TRUE)
# View the results
head(go_enrichment_banksy_cluster0)
# Plotting the GO enrichment results for Banksy cluster 0
barplot(go_enrichment_banksy_cluster0, showCategory = 10) +
  labs(title = "GO Enrichment for Banksy Cluster 0")
``` 

## Conclusion
In this section, we have gone through the standard analysis steps for Xenium data, including loading the data, preprocessing, dimensionality reduction, clustering, and finding marker genes. We have also used the Banksy package to perform spatially informed clustering on the Xenium data, and compared the Banksy clusters with the original Seurat clusters. 
This analysis provides a comprehensive overview of how to analyze Xenium data and how to interpret the results. The use of spatially informed clustering allows us to capture spatial patterns in the data, which can provide valuable insights into the biological processes underlying the data. Gene set enrichment analysis allows us to identify enriched pathways or biological processes in the clusters, providing further insights into the biological significance of the clusters.
This analysis can be extended to other Xenium datasets or spatial transcriptomics datasets, allowing for a comprehensive understanding of spatial patterns and gene expression in various biological contexts.

## Summary

::: {.callout-tip}
#### Key Points

- Xenium data can be analyzed using standard Seurat workflows
- Spatial clustering can be performed using the Banksy package
- Marker genes can be identified for both Seurat and Banksy clusters
- Comparison of clusters can provide insights into spatial patterns and gene expression patterns
- Gene set enrichment analysis can identify enriched pathways or biological processes in clusters
:::
