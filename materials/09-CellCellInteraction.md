---
title: Cell-Cell Interaction Analysis
---

::: {.callout-tip}
#### Learning Objectives

- Understand the concept of cell-cell interactions in spatial transcriptomics
- Learn how to analyze cell-cell interactions using the `CellChat` package
- Interpret the results of cell-cell interaction analysis
- Explore different methods for visualizing cell-cell interactions
:::

## Cell-Cell Interaction Analysis
Cell-cell interactions play a crucial role in various biological processes, including development, immune response, and disease progression. In spatial transcriptomics, analyzing cell-cell interactions can provide insights into how different cell types communicate and influence each other within the tissue microenvironment.
We will be using the `CellChat` package, which provides a comprehensive framework for analyzing and visualizing cell-cell communication networks from single-cell and spatial transcriptomics data.

```r
# Load necessary libraries
library(CellChat)
```

## Prepare Seurat Data for CellChat
First, we need to prepare our Seurat object for CellChat analysis. This involves extracting the expression data and metadata from the Seurat object.

```r
data.input <- visium[["SCT"]]$data
Idents(visium) <- "first_type"
labels <- Idents(visium)  #expects the cell labels to be in the identities of the Seurat object
meta <- data.frame(labels = labels, row.names = names(labels)) 
```

To use spatial information in the analysis, we can add the spatial coordinates of each spot to the metadata. This step is optional but can provide additional context for interpreting cell-cell interactions. 
For 10X Visium, this information is in the  'tissue_positions.csv' file found in the 'spatial' directory, the data stored in the Seurat object is not appropriate for CellChat. We will therefore read in the file directly save the coordinates in a dataframe.

```r
spatial_locs = GetTissueCoordinates(visium, scale = NULL)[ , c("x", "y")] 
```
Additionally we need the conversion factor from pixel coordinates to microns. This can be computed as the ratio of the spot diameter in microns (55 microns for 10X Visium) to the spot diameter in pixels (20 pixels for 10X Visium). 
We need to read the spot diameter in pixels from the 'scalefactors_json.json' file found in the 'spatial' directory.

```r
scalefactors = jsonlite::fromJSON(txt = file.path("data/mouse_sagittal/spatial", 'scalefactors_json.json'))
spotsize = 65 # the theoretical spot size (um) in 10X Visium
conversion_factor = spotsize/scalefactors$spot_diameter_fullres
spatial_factors = data.frame(ratio = conversion_factor, tol = spotsize/2)
```

Now we can check the distances of the spots to each other, which will be used in the CellChat analysis. For 10X Visium v1 data it should be roughly 100 microns between the centers of two spots.

```r
#check distances of the spots to each other
dist <- computeCellDistance(coordinates = spatial_locs, ratio = spatial_factors$ratio, tol = spatial_factors$tol)
min(dist[dist!=0]) 
```

## Create CellChat Object & Set the ligand-receptor database
Next, we create a CellChat object using the expression data and metadata as well as the spatial coordinates we have just prepared.

```r
cellchat <- createCellChat(object = data.input, meta = meta, group.by = "labels", datatype = "spatial", coordinates = spatial_locs, spatial.factors = spatial_factors)
```

Then we load and set the ligand-receptor interaction database. CellChat provides several built-in databases, and we will use the mouse database for this example. We are excluding the "Non-protein Signaling" category for our analysis, which is a standard subset of interactions to consider. We will also remove unused expression data from the CellChat object to speed up the analysis, matching the selected database subset.

```r
CellChatDB <- CellChatDB.mouse
showDatabaseCategory(CellChatDB)

# use all CellChatDB except for "Non-protein Signaling" for cell-cell communication analysis
CellChatDB.use <- subsetDB(CellChatDB)

# set the used database in the object
cellchat@DB <- CellChatDB.use

# subset the expression data of signaling genes for saving computation cost
cellchat <- subsetData(cellchat) 
```

## Identify Overexpressed Genes and Interactions
We will identify overexpressed genes and interactions in the dataset. This step helps to focus the analysis on the most relevant genes and interactions. CellChat offers parallel processing, but this often actually causes more issues than speedup, so we will not be using it. 
To compute the communication probability at the signaling gene level, we will use the "triMean" method, which is a robust method for estimating communication probabilities. We need to set the 'contact.range' parameter to define the maximum distance for considering cell-cell interactions. Here, we set it to 100 microns, which is the expected distance between the centers of two spots in 10X Visium data. We will also filter the interactions to only consider those that are supported by at least 10 cells.

```r

# Identify overexpressed genes and interactions
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)

#The following function computes the communication probability at the signaling pathway level
#If it takes a long time to compute, so you can skip it and load precomputed results instead
cellchat <- computeCommunProb(cellchat, type = "triMean", contact.range = 100, trim = 0.1)
cellchat <- computeCommunProbPathway(cellchat)

cellchat <- filterCommunication(cellchat, min.cells = 10) #only consider interactions that are supported by at least 10 cells
cellchat <- aggregateNet(cellchat)

# cellchat <- loadRDS("precomputed/mouse_brain_cellchat.rds") #load precomputed results
```

## Visualize Cell-Cell Communication Network
Finally, we can visualize the cell-cell communication network using various plotting functions provided by the CellChat package. We will create a circle plot to visualize the overall communication network and a heatmap to show the interaction strength between different cell types as well as a bar plot ranking the signaling pathways by their acttivity scores. 

```r
#find all significant signalling pathways
pathways.show.all <- cellchat@netP$pathways
pathways.selected <- c("COLLAGEN", "LAMININ")
laminin <- c("LAMININ")

# Circle plos of the overall communication network
groupSize <- as.numeric(table(cellchat@idents))
 netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T  , label.edge= F, title.name = "Number of interactions", top = 0.3)
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T  , label.edge= F, title.name = "Interaction weights/strength", top = 0.3)  


# Heatmap of interaction strength between cell types
netVisual_heatmap(cellchat, measure = "weight") 

#Bar plot of the contribution of each ligand-receptor interaction pair to the selected pathways
netAnalysis_contribution(cellchat, signaling = pathways.selected)
netAnalysis_contribution(cellchat, signaling = laminin)

# Bubble Plot of ligand-receptor interactions from one celltype (the 1st, L2/3 IT) to other celltypes (the 3rd to 5th)
netVisual_bubble(cellchat, sources.use = 1, targets.use = c(3:5), remove.isolate = FALSE)


#Centrality scores of the selected pathways
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
netAnalysis_signalingRole_network(cellchat, signaling = pathways.selected, width =10, height = 4)

# Incoming vs outgoing signalling strength for all signalling roles in the network per celltype for all pathways and separately for the laminin pathway
netAnalysis_signalingRole_scatter(cellchat)
netAnalysis_signalingRole_scatter(cellchat, signaling = c("LAMININ"))
```

Of course, we can also visualize specific signaling pathways of interest in spatial context.

```r
netVisual_aggregate(cellchat, signaling = pathways.selected, layout = "spatial", edge.width.max = 2, vertex.size.max = 1, alpha.image = 0.2, vertex.label.cex = 3.5)
netVisual_aggregate(cellchat, signaling = laminin, layout = "spatial", edge.width.max = 2, vertex.size.max = 1, alpha.image = 0.2, vertex.label.cex = 3.5)
```

There is a plethora of other visualization and investigation options available in the CellChat package. You can explore these options in the CellChat documentation to find the best way to visualize your specific data and research questions.

## Conclusion
In this section, we have explored how to analyze cell-cell interactions using the CellChat package in R. We prepared the Seurat object for CellChat analysis, identified overexpressed genes and interactions, and visualized the cell-cell communication network. Understanding cell-cell interactions can provide valuable insights into tissue architecture and function, helping to unravel complex biological processes. 

## Summary 
::: {.callout-tip}
#### Key Points   
- Cell-cell interactions are crucial for understanding tissue microenvironments and biological processes.
- The `CellChat` package provides a comprehensive framework for analyzing and visualizing cell-cell communication networks.
- Preparing the Seurat object and identifying overexpressed genes and interactions are essential steps in the analysis.
- For spatial transcriptomics data, incorporating spatial coordinates can enhance the analysis of cell-cell interactions.
- Visualization techniques such as circle plots, heatmaps, and bar plots help interpret the results of cell-cell interaction analysis.
:::


