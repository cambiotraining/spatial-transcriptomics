---
title: "Spatial transcriptomics analysis"
date: today
number-sections: false
---

## Overview

This course provides a practical, end-to-end introduction to spatial transcriptomics analysis in R using Seurat and related tools.
You will learn how to assess data quality, load and manage datasets from major platforms, and preprocess data for robust downstream analysis.
The course then covers dimensionality reduction, clustering, deconvolution, and visualisation strategies for spatial expression patterns.
You will also explore spatially variable features, tissue architecture, and cell-cell communication to interpret biological structure in context.
Hands-on platform examples (Visium, Xenium, and MERFISH) are used to compare analytical choices and reproducible workflows.

::: {.callout-tip}
### Learning Objectives

- Assess spatial transcriptomics data quality using platform reports and quantitative QC metrics.
- Load, inspect, and manage Visium, Xenium, and MERFISH datasets in Seurat.
- Preprocess spatial data through filtering, normalisation, and feature preparation.
- Apply and tune dimensionality reduction and clustering for spatial datasets.
- Perform spot-level deconvolution with single-cell references to estimate cell-type composition.
- Identify and interpret marker genes, spatially variable features, and tissue architecture patterns.
- Analyse cell-cell communication networks and relate signalling to spatial tissue context.
:::

### Target Audience

This course is aimed at researchers who are new to **spatial transcriptomics data analysis** and want practical experience analysing, visualising, and interpreting spatial omics datasets in R.

### Prerequisites

- **Required**
  - Basic understanding of high-throughput sequencing technologies.
    - Watch [this iBiology video](https://youtu.be/mI0Fo9kaWqo) for an excellent overview.
  - A working knowledge of the UNIX command line ([course registration page](https://training.csx.cam.ac.uk/bioinformatics/course/bioinfo-unix2)).
    - If you are not able to attend this prerequisite course, please work through our [Unix command line materials](https://cambiotraining.github.io/unix-shell/) ahead of the course (up to section 7).
  - A working knowledge of R ([course registration page](https://training.csx.cam.ac.uk/bioinformatics/course/bioinfo-introRbio)).
    - If you are not able to attend this prerequisite course, please work through [our R materials](https://cambiotraining.github.io/intro-r/) ahead of the course.
- **Suggested**
  - Some experience with single-cell RNA-seq analysis in R.
    - For example, attend our [single-cell RNA-seq course](https://training.csx.cam.ac.uk/bioinformatics/course/bioinfo-singlecell-rna-seq) or work through [our single-cell RNA-seq materials](https://cambiotraining.github.io/single-cell-rnaseq/).

<!-- Training Developer note: comment the following section out if you did not assign levels to your exercises -->
<!--
### Exercises

Exercises in these materials are labelled according to their level of difficulty:

| Level | Description |
| ----: | :---------- |
| {{< fa solid star >}} {{< fa regular star >}} {{< fa regular star >}} | Exercises in level 1 are simpler and designed to get you familiar with the concepts and syntax covered in the course. |
| {{< fa solid star >}} {{< fa solid star >}} {{< fa regular star >}} | Exercises in level 2 combine different concepts together and apply it to a given task. |
| {{< fa solid star >}} {{< fa solid star >}} {{< fa solid star >}} | Exercises in level 3 require going beyond the concepts and syntax introduced to solve new problems. |
-->

## Citation & Authors

Please cite these materials if:

- You adapted or used any of them in your own teaching.
- These materials were useful for your research work.
  For example, you can cite us in the methods section of your paper: "We carried our analyses based on the recommendations in *YourReferenceHere*".

<!--
This is generated automatically from the CITATION.cff file.
If you think you should be added as an author, please get in touch with us.
-->

{{< citation CITATION.cff >}}

## Acknowledgements

<!-- if there are no acknowledgements we can delete this section -->

- Thank you to 10X Genomics for providing most of the datasets used in these materials.
- Thank you to the Satija lab for developing and maintaining the Seurat package and for their extensive documentation and tutorials.
- Thank you to the CambiO Training team for their support in developing and maintaining these materials.
