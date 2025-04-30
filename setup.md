---
title: "Data & Setup"
number-sections: false
---

<!-- 
Note for Training Developers:
We provide instructions for commonly-used software as commented sections below.
Uncomment the sections relevant for your materials, and add additional instructions where needed (e.g. specific packages used).
Note that we use tabsets to provide instructions for all three major operating systems.
-->

::: {.callout-tip level=2}
## Workshop Attendees

If you are attending one of our workshops, we will provide a training environment with all of the required software and data.  
If you want to setup your own computer to run the analysis demonstrated on this course, you can follow the instructions below.
:::

## Data

The data used in these materials is provided as a zip file. 
Download and unzip the folder to your Desktop to follow along with the materials.

<!-- Note for Training Developers: add the link to 'href' -->
<a href="">
  <button class="btn"><i class="fa fa-download"></i> Download</button>
</a>

## Software

- Install R following <a href="https://cambiotraining.github.io/software-installation/materials/r-base.html" target="_blank">these instructions</a>
- Open RStudio and run the following command to install the packages used in this course: 

    ```r
    if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

    BiocManager::install("Seurat")
    ```