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

