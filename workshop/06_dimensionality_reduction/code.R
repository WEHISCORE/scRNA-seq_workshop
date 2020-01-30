# Illustrative dataset: Zeisel -------------------------------------------------

library(scRNAseq)
sce.zeisel <- ZeiselBrainData(ensembl=TRUE)

# Quality control
library(scater)
is.mito <- which(rowData(sce.zeisel)$featureType=="mito")
stats <- perCellQCMetrics(sce.zeisel, subsets=list(Mt=is.mito))
qc <- quickPerCellQC(stats, percent_subsets=c("altexps_ERCC_percent", "subsets_Mt_percent"))
sce.zeisel <- sce.zeisel[,!qc$discard]

# normalization
library(scran)
set.seed(1000)
clusters <- quickCluster(sce.zeisel)
sce.zeisel <- computeSumFactors(sce.zeisel, cluster=clusters)
sce.zeisel <- logNormCounts(sce.zeisel)

# variance-modelling
dec.zeisel <- modelGeneVarWithSpikes(sce.zeisel, "ERCC")

