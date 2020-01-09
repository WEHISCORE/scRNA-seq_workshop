# Single-Cell RNA-seq Analysis with Bioconductor Workshop

This 2-day workshop is based on the free e-book, ['Orchestrating Single-Cell Analysis with Bioconductor'](https://osca.bioconductor.org), and the analysis strategies used by WEHI's Single Cell Open Research Endeavour (SCORE).
The workshop will focus on software tools available through [Bioconductor](https://bioconductor.org/).
Bioconductor is an open source, open development software project to provide tools for the analysis and comprehension of high-throughput genomic data.
It is based primarily on the R programming language.

## Who this workshop is for

This workshop is designed with the interested experimental biologist in mind, and we will do our best to make few assumptions on previous programming or statistical experience.
Likewise, we also welcome more seasoned bioinformaticians who are looking for a starting point from which to dive into single-cell RNA-seq analysis.

## What you will learn

The goal of this workshop is to provide a solid foundation in the usage of Bioconductor tools for single-cell RNA-seq analysis by walking through various steps of typical workflows using example datasets. 
Our examples are taken from a range of commonly used single-cell RNA-seq technologies, such as the 10X Genomics Chromium Single Cell Gene Expression and the CEL-Seq2 platforms, applied to a variety of human and mouse experimental systems.

All workflows begin with *data import* and subsequent *quality control* and *normalization*, going from a raw (count) expression matrix to a 'clean' one.
This includes adjusting for experimental factors and complications such as batch effects.
Using the clean expression matrix, *feature selection* strategies can be applied to select the features (genes) driving heterogeneity.
Furthermore, these features can then be used to perform *dimensionality reduction*, which enables downstream analysis that would not otherwise be possible and visualization in 2 or 3 dimensions.

From there, the workflows largely focus on differing downstream analyses.
*Clustering* details how to segment a single-cell RNA-seq dataset, and *differential expression* provides a means to determine what drives the differences between different groups of cells.
*Integrating datasets* walks through merging scRNA-seq datasets, an area of need as the number of scRNA-seq datasets continues to grow and comparisons between datasets must be done.

At each step in the workflow, we will pay special attention to data visualization, including interactive visualization, to help guide us in our analysis decisions and to interpret the results of our analysis.

## What you won’t learn

This workshop focuses on the downstream analysis of scRNA-seq data and our workflows will start from a raw (count) expression matrix, such as that created by 10X Genomics' CellRanger software.
We will therefore not be learning to run CellRanger (or similar software) that construct the count matrix from the raw sequencing files, although we will introduce these tools, explain when to use which, and link to tutorials.

The field of bioinformatic analysis is large and filled with many potential trajectories depending on the biological system being studied and technology being deployed.
Here, we only briefly survey some of the many tools available for the analysis of scRNA-seq, focusing on Bioconductor packages.
It is impossible to thoroughly review the plethora of tools available through R and Bioconductor for biological analysis in one workshop, but we aim to provide the means for further exploration on your own.

Thus, it goes without saying that you may not learn the optimal workflow for your own data from our examples; while we strive to provide high quality templates, they should be treated as just that, a template from which to extend upon for your own analyses.

## Expectations

This workshop assumes some familiarity with R.
We assume that you have attended WEHI's intro to R workshops or otherwise consider yourself to be able to use R.

## Course structure

### Day 1

1. Introduction
2. Data infrastructure and import
3. Quality control
4. Normalization
5. Feature selection
6. Dimensionality reduction

### Day 2

7. Clustering
8. Marker gene detection
9. Cell type annotation
10. Integrating datasets
11. Multi-sample comparisons
12. Additional topics
13. Recap
