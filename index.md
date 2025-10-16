---
title: "Spatial transcriptomics analysis"
date: today
number-sections: false
---

## Overview 

This course provides a practical introduction to analyzing spatial transcriptomics data using the Seurat package in R and related tools. Participants will learn how to process spatial transcriptomics datasets, perform quality control, normalization, and integration, and apply clustering and dimensionality reduction techniques. The course teaches visualization of spatial gene expression patterns, identification of spatially variable features, and inference of cell-cell interactions. Through hands-on exercises, students will gain proficiency in Seurat’s spatial analysis tools and develop reproducible workflows for high-throughput spatial transcriptomics projects.

::: {.callout-tip}
### Learning Objectives

- Load and preprocess spatial transcriptomics data using Seurat
- Perform quality control and normalization of spatial transcriptomics datasets
- Apply clustering and dimensionality reduction techniques
- Visualize spatial gene expression patterns and clusters
- Identify spatially variable features and marker genes
- Infer cell-cell interactions in spatial contexts
- Develop reproducible workflows for spatial transcriptomics analysis
:::


### Target Audience

This course is aimed at researchers with **no prior experience in the analysis of spatial transcriptomics data**, who would like to learn how to analyse, visualise and extract insights from spatial omics datasets.


### Prerequisites

- Basic understanding of high-throughput sequencing technologies.
  - Watch [this iBiology video](https://youtu.be/mI0Fo9kaWqo) for an excellent overview. 
- A working knowledge of the UNIX command line ([course registration page](https://training.csx.cam.ac.uk/bioinformatics/course/bioinfo-unix2)).
  - If you are not able to attend this prerequisite course, please work through our [Unix command line materials](https://cambiotraining.github.io/unix-shell/) ahead of the course (up to section 7). 
- A working knowledge of R ([course registration page](https://training.csx.cam.ac.uk/bioinformatics/course/bioinfo-introRbio)).
  - If you are not able to attend this prerequisite course, please work through [our R materials](https://cambiotraining.github.io/intro-r/) ahead of the course.
  - It would be beneficial to have some experience with single-cell RNA-seq analysis in R, for example by attending our [single-cell RNA-seq course](https://training.csx.cam.ac.uk/bioinformatics/course/bioinfo-singlecell-rna-seq) or working through [our single-cell RNA-seq materials](https://cambiotraining.github.io/single-cell-rnaseq/).


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
- These materials were useful for your research work. For example, you can cite us in the methods section of your paper: "We carried our analyses based on the recommendations in _YourReferenceHere_".

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

