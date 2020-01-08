# scRNA-seq workshop

The workshop will be based on https://osca.bioconductor.org.
Some of the course structure may be adapted from https://combine-australia.github.io/2018-09-26-RNAseq-Melbourne/.

## Date(s)

This workshop requires BioC 3.10, which will be released on 2019-10-30, so I propose the workshop is held in mid-November.

**TODO**: ~Propose running the workshop between Nov 11-28.~ Now pushed to February/March 2020.

I would like to teach as a 2-day workshop.
This is motivated by the Hemberg lab's scRNA-seq workshop and the COMBINE RNA-seq workshop, both of which run over two days.

## What you will learn

> The goal of this ~~book~~ workshop is to provide a solid foundation in the usage of Bioconductor tools for single-cell RNA-seq analysis by walking through various steps of typical workflows using example datasets. We strive to tackle key concepts covered in the manuscript, [“Orchestrating Single-Cell Analysis with Bioconductor”](https://www.biorxiv.org/content/10.1101/590562v1), with each workflow covering these in varying detail, as well as essential preliminaries that are important for following along with the workflows on your own.

## What you won’t learn

> The field of bioinformatic analysis is large and filled with many potential trajectories depending on the biological system being studied and technology being deployed. Here, we only briefly survey some of the many tools available for the analysis of scRNA-seq, focusing on Bioconductor packages. It is impossible to thoroughly review the plethora of tools available through R and Bioconductor for biological analysis in one ~~book~~ workshop, but we hope to provide the means for further exploration on your own.
>
> Thus, it goes without saying that you may not learn the optimal workflow for your own data from our examples - while we strive to provide high quality templates, they should be treated as just that - a template from which to extend upon for your own analyses.

- I assume that you have attended WEHI's intro R workshops or otherwise consider yourself to be able to use R
    - **TODO**: Check what is taught in WEHI's intro to R workshops
- We won't run CellRanger, scPipe, alevin, STARsolo, kalisto | bustools, etc. but I will introduce them, explain when to use which, and link to tutorials.

## Content

### Session 1

1. Overview
2. Data Infrastructure
    - Abbreviate from OSCA
3. Quality control
4. Normalization
5. Feature selection
6. Dimensionality reduction
7. Interactive interfaces

### Session 2

8. Clustering
9. Marker gene detection
10. Cell type annotation
11. Integrating datasets
12. Multi-sample comparisons
13. Doublet detection
14. Other resources

## Format

I'd like to present from within RStudio, writing code as we go.
This might be supplemented by Google Slides.

Everyone will be running RStudio.
If this a WEHI-only workshop, we could use http://rstudio.hpc.wehi.edu.au.
However, I would need to ensure sufficient resources are provisioned and that all the necessary packages are installed system-wide.
If this is open to non-WEHI people, I will use AWS.
