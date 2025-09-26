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
library(Seurat)
library(CellChat) 
library(patchwork)
library(ggplot2)
library(dplyr)
```

## Prepare Seurat Data for CellChat
First, we need to prepare our Seurat object for CellChat analysis. This involves extracting the expression data and metadata from the Seurat object.

```r
data.input <- visium[["SCT"]]$data
Idents(visium) <- "first_type"
labels <- Idents(visium)  #expects the cell labels to be in the identities of the Seurat object
meta <- data.frame(labels = labels, row.names = names(labels)) 
```

## Create CellChat Object & Set the ligand-receptor database
Next, we create a CellChat object using the expression data and metadata. 

```r
cellchat <- createCellChat(object = data.input, meta = meta, group.by = "labels")
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

```r

# Identify overexpressed genes and interactions
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)

#The following function computes the communication probability at the signaling pathway level
#If it takes a long time to compute, so you can skip it and load precomputed results instead
cellchat <- computeCommunProb(cellchat, type = "triMean")
cellchat <- computeCommunProbPathway(cellchat)

cellchat <- filterCommunication(cellchat, min.cells = 10) #only consider interactions that are supported by at least 10 cells
cellchat <- aggregateNet(cellchat)

# cellchat <- loadRDS("precomputed/cmouse_brain_cellchat.rds") #load precomputed results
```

## Visualize Cell-Cell Communication Network
Finally, we can visualize the cell-cell communication network using various plotting functions provided by the CellChat package. We will create a circle plot to visualize the overall communication network and a heatmap to show the interaction strength between different cell types as well as a bar plot ranking the signaling pathways by their acttivity scores. 

```r
#find all significant signalling pathways
pathways.show.all <- cellchat@netP$pathways
pathways.selected <- c("COLLAGEN", "LAMININ")
laminin <- c("LAMININ")

# Circle plot of the overall communication network
groupSize <- as.numeric(table(cellchat@idents))
c1 <- netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T  , label.edge= F, title.name = "Number of interactions", top = 0.3)
c2 <- netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T  , label.edge= F, title.name = "Interaction weights/strength", top = 0.3)  
c1 + c2

# Heatmap of interaction strength between cell types
netVisual_heatmap(cellchat, measure = "weight") 

#Bar plot of the contribution of each ligand-receptor interaction pair to the selected pathways
netAnalysis_contribution(cellchat, signaling = pathways.selected)
netAnalysis_contribution(cellchat, signaling = laminin)

# Bubble Plot of ligand-receptor interactions from one celltype (the 1st, L2/3 IT) to other celltypes (the 2nd to 5th)
netVisual_bubble(cellchat, sources.use = 1, targets.use = c(2:4), remove.isolate = FALSE)


#Centrality scores of the selected pathways
netAnalysis_signalingRole_network(cellchat, signaling = pathways.selected, width =10, height = 4)

# Incoming vs outgoing signalling strength for all signalling roles in the network per celltype for all pathways and separately for the laminin pathway
netAnalysis_signalingRole_scatter(cellchat)
netAnalysis_signalingRole_scatter(cellchat, signaling = c("LAMININ"))

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
- Visualization techniques such as circle plots, heatmaps, and bar plots help interpret the results of cell-cell interaction analysis.
:::


