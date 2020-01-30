# Illustrative dataset: Zeisel -------------------------------------------------

library(scRNAseq)
sce.zeisel <- ZeiselBrainData(ensembl=TRUE)

# Quality control
library(scater)
is.mito <- which(rowData(sce.zeisel)$featureType=="mito")
stats <- perCellQCMetrics(sce.zeisel, subsets=list(Mt=is.mito))
qc <- quickPerCellQC(stats, percent_subsets=c("altexps_ERCC_percent", "subsets_Mt_percent"))
sce.zeisel <- sce.zeisel[,!qc$discard]

# Library size factors ---------------------------------------------------------

lib.sf.zeisel <- librarySizeFactors(sce.zeisel)

# Examine distribution of size factors
summary(lib.sf.zeisel)
hist(log10(lib.sf.zeisel), xlab="Log10[Size factor]", col="grey80")
ls.zeisel <- colSums(counts(sce.zeisel))
plot(ls.zeisel, lib.sf.zeisel, log="xy", xlab="Library size", ylab="Size factor")

# Normalization by convolution -------------------------------------------------

library(scran)
# Pre-clustering
set.seed(100)
clust.zeisel <- quickCluster(sce.zeisel)
# Compute deconvolution size factors
deconv.sf.zeisel <- calculateSumFactors(sce.zeisel, cluster=clust.zeisel, min.mean=0.1)

# Examine distribution of size factors
summary(deconv.sf.zeisel)
hist(log10(deconv.sf.zeisel), xlab="Log10[Size factor]",
     col="grey80")
plot(ls.zeisel, deconv.sf.zeisel, log="xy", xlab="Library size",
     ylab="Size factor")

# Library size factors vs. convolution size factors ----------------------------

# Colouring points using the supplied cell-types
plot(lib.sf.zeisel, deconv.sf.zeisel, xlab="Library size factor",
     ylab="Deconvolution size factor", log='xy', pch=16,
     col=as.integer(factor(sce.zeisel$level1class)))
abline(a=0, b=1, col="red")


