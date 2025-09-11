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
data.input <- visium[["RNA"]]$data
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
We will identify overexpressed genes and interactions in the dataset. This step helps to focus the analysis on the most relevant genes and interactions. To speed up the process, we will use parallel processing. At the end we will project the results onto a protein-protein interaction (PPI) network.

```r
# Initialize parallel processing
future::plan("multisession", workers = 4) # do parallel

# Identify overexpressed genes and interactions
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)

# Project gene expression data onto PPI network 
cellchat <- projectData(cellchat, PPI.mouse)

##TODO: go on with the rest of the analysis
