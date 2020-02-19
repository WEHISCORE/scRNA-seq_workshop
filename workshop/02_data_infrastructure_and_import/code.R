# Illustrative dataset: 416B ---------------------------------------------------

library(scRNAseq)
sce.416b <- LunSpikeInData(which="416b")

# Filling the assays slot ------------------------------------------------------

# Load the SingleCellExperiment package
library(SingleCellExperiment)
# Extract the count matrix from the 416b dataset
counts.416b <- counts(sce.416b)
# Construct a new SCE from the counts matrix
sce <- SingleCellExperiment(
  assays = list(counts = counts.416b))

# Inspect the object we just created
sce

# Access the counts matrix from the assays slot

# 1. The general method
assay(sce, "counts")
# The special method for the assay named "counts"
counts(sce)

# Tip: Limit the output to just a few samples & genes
counts(sce)[1:30, 1:2]

# Adding to the assays slot ----------------------------------------------------

sce <- scater::logNormCounts(sce)
# Inspect the object we just created
sce

# Access the logcounts matrix from the assays slot
# WARNING: This will flood RStudio with output!

# 1. The general method
assay(sce, "logcounts")
# 2. The special method for the assay named "logcounts"
logcounts(sce)

# Tip: Limit the output to just a few samples & genes
logcounts(sce)[1:30, 1:2]

# assign a new entry to assays slot
assay(sce, "counts_100") <- assay(sce, "counts") + 100
# List the assays in the object
assays(sce)

# Filling the colData slot -----------------------------------------------------

# Extract the sample metadata from the 416b dataset
colData.416b <- colData(sce.416b)
# Add some of the sample metadata to our SCE
colData(sce) <- colData.416b[, c("phenotype", "block")]
# Inspect the object we just updated
sce

# Adding to the colData slot ---------------------------------------------------

# Example of function that adds extra fields to colData
sce <- scater::addPerCellQC(sce)
# Access the sample metadata from our updated SCE
colData(sce)

# Using colData for subsetting -------------------------------------------------

# Subset data to just wild type cells
sce[, sce$phenotype == "wild type phenotype"]

# Adding to the rowData slot ---------------------------------------------------

# Access the feature metadata from our SCE
# It's currently empty!
rowData(sce)

# Example of function that adds extra fields to colData
sce <- scater::addPerFeatureQC(sce)
# Access the feature metadata from our updated SCE
rowData(sce)

# Download the relevant Ensembl annotation database
# using AnnotationHub resources
library(AnnotationHub)
ah <- AnnotationHub()
ah_id <- query(ah, c("Mus musculus", "Ensembl", "v97"))

# Annotate each gene with its chromosome location
ensdb <- ah[["AH73905"]]
chromosome <- mapIds(ensdb, keys=rownames(sce),
                     keytype="GENEID", column="SEQNAME")
rowData(sce)$chromosome <- chromosome

# Access the feature metadata from our updated SCE
rowData(sce)

# Adding to the metadata slot --------------------------------------------------

# Access the metadata from our SCE
# It's currently empty!
metadata(sce)

# The metadata slot is Vegas - anything goes
metadata(sce) <- list(
  favourite_genes=c("Shh", "Nck1", "Diablo"),
  analyst=c("Pete"))

# Access the metadata from our updated SCE
metadata(sce)

# Adding to the reducedDims slot -----------------------------------------------

sce <- scater::runPCA(sce)
# Inspect the object we just updated
sce
# Access the PCA matrix from the reducedDims slot
reducedDim(sce, "PCA")

sce <- scater::runTSNE(sce)
# Inspect the object we just updated
sce
# Access the t-SNE matrix from the reducedDims slot
reducedDim(sce, "TSNE")

# 'Manually' compute UMAP
u <- uwot::umap(t(logcounts(sce)), n_component=2)
# Add the UMAP matrix to the reducedDims slot
# Access the t-SNE matrix from the reducedDims slot
reducedDim(sce, "UMAP") <- u

# List the dimensionality reduction results stored in # the object
reducedDims(sce)

# Adding an alternative experiment ---------------------------------------------

# Extract the ERCC SCE from the 416b dataset
ercc.sce.416b <- altExp(sce.416b, "ERCC")
# Inspect the ERCC SCE
ercc.sce.416b

# Add the ERCC SCE as an alternative experiment to our SCE
altExp(sce, "ERCC") <- ercc.sce.416b
# Inspect the object we just updated
sce
# List the alternative experiments stored in the object
altExps(sce)

# Why use alternative experiments? ---------------------------------------------

# Subsetting the SCE by sample also subsets the
# alternative experiments
sce.subset <- sce[, 1:10]
ncol(sce.subset)
ncol(altExp(sce.subset))

# Adding size factors ----------------------------------------------------------

# Extract existing size factors (these were added
# when we ran scater::logNormCounts(sce))
sizeFactors(sce)

# 'Automatically' replace size factors
sce <- scran::computeSumFactors(sce)
sizeFactors(sce)

# 'Manually' replace size factors
sizeFactors(sce) <- scater::librarySizeFactors(sce)
sizeFactors(sce)

# I ran CellRanger -------------------------------------------------------------

# Download example data processed with CellRanger
# Aside: Using BiocFileCache means we only download the
#        data once
library(BiocFileCache)
bfc <- BiocFileCache()
pbmc.url <- "http://cf.10xgenomics.com/samples/cell-vdj/3.1.0/vdj_v1_hs_pbmc3/vdj_v1_hs_pbmc3_filtered_feature_bc_matrix.tar.gz"
pbmc.data <- bfcrpath(bfc, pbmc.url)

# Extract the files to a temporary location
untar(pbmc.data, exdir=tempdir())

# List the files we downloaded and extracted
# These files are typically CellRanger outputs
pbmc.dir <- file.path(tempdir(),
                      "filtered_feature_bc_matrix")
list.files(pbmc.dir)

# Import the data as a SingleCellExperiment
library(DropletUtils)
sce.pbmc <- read10xCounts(pbmc.dir)
# Inspect the object we just constructed
sce.pbmc

# Store the CITE-seq data in an alternative experiment
sce.pbmc <- splitAltExps(sce.pbmc, rowData(sce.pbmc)$Type)
# Inspect the object we just updated
sce.pbmc

# I ran scPipe -----------------------------------------------------------------

# Download example data processed with CellRanger
# Aside: Using BiocFileCache means we only download the
#        data once
library(BiocFileCache)
bfc <- BiocFileCache()
sis_seq.url <- "https://github.com/LuyiTian/SIS-seq_script/archive/master.zip"
sis_seq.data <- bfcrpath(bfc, sis_seq.url)

# Extract the files to a temporary location
unzip(sis_seq.data, exdir=tempdir())

# List (some of) the files we downloaded and extracted
# These files are typical scPipe outputs
sis_seq.dir <- file.path(tempdir(),
                         "SIS-seq_script-master", "data", "BcorKO_scRNAseq",
                         "RPI10")
list.files(sis_seq.dir)

# Import the data as a SingleCellExperiment
library(scPipe)
sce.sis_seq <- create_sce_by_dir(sis_seq.dir)
# Inspect the object we just constructed
sce.sis_seq

# I got a bunch of files -------------------------------------------------------

# Download example bunch o' files dataset
library(BiocFileCache)
bfc <- BiocFileCache()
lun_counts.url <- "https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-5522/E-MTAB-5522.processed.1.zip"
lun_counts.data <- bfcrpath(bfc, lun_counts.url)
lun_coldata.url <- "https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-5522/E-MTAB-5522.sdrf.txt"
lun_coldata.data <- bfcrpath(bfc, lun_coldata.url)

# Extract the counts files to a temporary location
lun_counts.dir <- tempfile("lun_counts.")
unzip(lun_counts.data, exdir=lun_counts.dir)

# List the files we downloaded and extracted
list.files(lun_counts.dir)

# Import the count matrix (for 1 plate)
lun.counts <- read.delim(
  file.path(lun_counts.dir, "counts_Calero_20160113.tsv"),
  header=TRUE,
  row.names=1,
  check.names=FALSE)
# Store the gene lengths for later
gene.lengths <- lun.counts$Length
# Convert the gene counts to a matrix
lun.counts <- as.matrix(lun.counts[, -1])

# Import the sample metadata
lun.coldata <- read.delim(lun_coldata.data,
                          check.names=FALSE, stringsAsFactors=FALSE)
library(S4Vectors)
lun.coldata <- as(lun.coldata, "DataFrame")

# Match up the sample metadata to the counts matrix
m <- match(
  colnames(lun.counts),
  lun.coldata$`Source Name`)
lun.coldata <- lun.coldata[m, ]

# Construct the feature metadata
lun.rowdata <- DataFrame(Length = gene.lengths)

# Construct the SingleCellExperiment
lun.sce <- SingleCellExperiment(
  assays = list(assays = lun.counts),
  colData = lun.coldata,
  rowData = lun.rowdata)
# Inspect the object we just constructed
lun.sce
