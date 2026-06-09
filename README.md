# Spatial transcriptomics analysis

This repository contains the materials for the course.

**Course Developers**: see our [guidelines page](https://cambiotraining.github.io/quarto-course-template/materials.html) if contributing materials.

These materials are released under a [CC BY 4.0](LICENSE.md) license.

--------------------------------------------------------------------------------

## Contributor guide: build and execution

This Quarto book contains the teaching materials for spatial transcriptomics (Seurat-based).

The `.qmd` files are the single source of truth for both analysis and teaching narrative.

There is no separate preprocessing pipeline.
Instead, expensive computations are selectively cached within Quarto.

At the end of each chapter, we save the intermediate `.rds` file for use in the following chapter.
This way learners can always start from each chapter, independently of having successfully completed the previous ones.

### First-time setup and rendering

The first full render of the book may take some time due to some expensive computations.
For example:

- Normalisation and scaling
- PCA / UMAP embedding
- Neighbour graph construction
- Clustering and differential expression

However, subsequent renders are significantly faster due to caching and frozen outputs.
The project uses two Quarto mechanisms:

1. Freeze (project-level)

  - Rendered outputs are stored in `_freeze/`
  - Used for GitHub Actions builds
  - Prevents unnecessary re-execution during publishing
  - If it changes, **always push the `_freeze` directory to GitHub**

2. Selective chunk caching

  - Computationally expensive R chunks use `cache: true`
  - Cached chunks store intermediate results on disk
  - This speeds up local development of the materials
  - The **cache is never pushed to GitHub**, it is local to your machine

Re-rendering only a chapter is often sufficient during development.
For example:

```bash
quarto render materials/01-introduction.qmd
```

### Cache behaviour

Caching is used only for expensive computational steps, and is marked in code chunks as:

```r
#| cache: true
```

- The chunk result is stored after first execution
- Future renders reuse the result automatically
- Code is NOT re-executed unless inputs or code change

### Unsure about caching new code chunks?

If you are unsure whether a chunk should be cached or a large file submitted to dropbox:

- leave it uncached
- open an issue or get in touch with us
