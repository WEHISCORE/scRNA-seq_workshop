# Illustrative datasets: Filtered 10X PBMC3k and PBMC4k ------------------------

# loading data
library(TENxPBMCData)
all.sce <- list(
  pbmc3k=TENxPBMCData("pbmc3k"),
  pbmc4k=TENxPBMCData("pbmc4k"))

# quality-control
library(scater)
stats <- high.mito <- list()
for (n in names(all.sce)) {
  current <- all.sce[[n]]
  is.mito <- grep("MT", rowData(current)$Symbol_TENx)
  stats[[n]] <- perCellQCMetrics(current,
                                 subsets=list(Mito=is.mito))
  high.mito[[n]] <- isOutlier(stats[[n]]$subsets_Mito_percent,
                              type="higher")
  all.sce[[n]] <- current[,!high.mito[[n]]]
}

# normalization
all.sce <- lapply(all.sce, logNormCounts)

# variance-modelling
library(scran)
all.dec <- lapply(all.sce, modelGeneVar)
all.hvgs <- lapply(all.dec, getTopHVGs, prop=0.1)

# dimensionality-reduction
library(BiocSingular)
set.seed(10000)
all.sce <- mapply(FUN=runPCA, x=all.sce,
                  subset_row=all.hvgs,
                  MoreArgs=list(ncomponents=25, BSPARAM=RandomParam()),
                  SIMPLIFY=FALSE)

set.seed(100000)
all.sce <- lapply(all.sce, runTSNE, dimred="PCA")
set.seed(1000000)
all.sce <- lapply(all.sce, runUMAP, dimred="PCA")

# clustering
for (n in names(all.sce)) {
  g <- buildSNNGraph(all.sce[[n]], k=10, use.dimred="PCA")
  clust <- igraph::cluster_walktrap(g)$membership
  all.sce[[n]]$cluster <- factor(clust)
}

pbmc3k <- all.sce$pbmc3k
dec3k <- all.dec$pbmc3k

pbmc4k <- all.sce$pbmc4k
dec4k <- all.dec$pbmc4k

# Preparing for batch correction -----------------------------------------------

# subset all data to a common 'universe' of features
universe <- intersect(rownames(pbmc3k), rownames(pbmc4k))

# Subsetting the SingleCellExperiment object
pbmc3k <- pbmc3k[universe,]
pbmc4k <- pbmc4k[universe,]

# Also subsetting the variance modelling results
dec3k <- dec3k[universe,]
dec4k <- dec4k[universe,]

# Rescale to adjust for differences in
# sequencing depth in batches
library(batchelor)
rescaled <- multiBatchNorm(pbmc3k, pbmc4k)
pbmc3k <- rescaled[[1]]
pbmc4k <- rescaled[[2]]

# Average the variance components across batches
library(scran)
combined.dec <- combineVar(dec3k, dec4k)
chosen.hvgs <- combined.dec$bio > 0

# Diagnosing batch effects -----------------------------------------------------

# Synchronizing the metadata before combining
rowData(pbmc3k) <- rowData(pbmc4k)
pbmc3k$batch <- "3k"
pbmc4k$batch <- "4k"
uncorrected <- cbind(pbmc3k, pbmc4k)

# Perform PCA on the log-expression values for all genes
# with a positive (average) biological component
library(scater)
set.seed(0010101010)
uncorrected <- runPCA(uncorrected, subset_row=chosen.hvgs)

library(scran)
snn.gr <- buildSNNGraph(uncorrected, use.dimred="PCA")
clusters <- igraph::cluster_walktrap(snn.gr)$membership
tab <- table(Cluster=clusters, Batch=uncorrected$batch)

set.seed(1111001)
uncorrected <- runTSNE(uncorrected, dimred="PCA")
plotTSNE(uncorrected, colour_by="batch")

# Linear regression ------------------------------------------------------------

library(batchelor)
rescaled <- rescaleBatches(pbmc3k, pbmc4k)

set.seed(1010101010)
rescaled <- runPCA(rescaled, subset_row=chosen.hvgs,
                   exprs_values="corrected")

snn.gr <- buildSNNGraph(rescaled, use.dimred="PCA")
clusters.resc <- igraph::cluster_walktrap(snn.gr)$membership
tab.resc <- table(Cluster=clusters.resc,
                  Batch=rescaled$batch)

rescaled <- runTSNE(rescaled, dimred="PCA")
rescaled$batch <- factor(rescaled$batch)
plotTSNE(rescaled, colour_by="batch")

# MNN correction ---------------------------------------------------------------

set.seed(1000101001)
mnn.out <- fastMNN(pbmc3k, pbmc4k, d=50, k=20,
                   subset.row=chosen.hvgs,
                   BSPARAM=BiocSingular::RandomParam(deferred=TRUE))

# Correction diagnostics -------------------------------------------------------

library(scran)
snn.gr <- buildSNNGraph(mnn.out, use.dimred="corrected")
clusters.mnn <- igraph::cluster_walktrap(snn.gr)$membership
tab.mnn <- table(Cluster=clusters.mnn, Batch=mnn.out$batch)

library(scater)
set.seed(0010101010)
mnn.out <- runTSNE(mnn.out, dimred="corrected")
mnn.out$batch <- factor(mnn.out$batch)
plotTSNE(mnn.out, colour_by="batch")

# Encouraging consistency with marker genes ------------------------------------

stats3 <- pairwiseWilcox(logcounts(pbmc3k), pbmc3k$cluster,
                         direction="up")
markers3 <- getTopMarkers(stats3[[1]], stats3[[2]], n=10)

stats4 <- pairwiseWilcox(logcounts(pbmc4k), pbmc4k$cluster,
                         direction="up")
markers4 <- getTopMarkers(stats4[[1]], stats4[[2]], n=10)

marker.set <- unique(unlist(c(unlist(markers3), unlist(markers4))))

set.seed(1000110)
mnn.out2 <- fastMNN(pbmc3k, pbmc4k, subset.row=marker.set,
                    BSPARAM=BiocSingular::RandomParam(deferred=TRUE))
mnn.out2 <- runTSNE(mnn.out2, dimred="corrected")
mnn.out2$cluster <- factor(
  c(pbmc3k$cluster, pbmc4k$cluster))
plotTSNE(mnn.out2, colour_by="cluster",
         other_fields="batch") +
  facet_grid(~batch)

# Using the corrected values ---------------------------------------------------

m.out <- findMarkers(uncorrected, clusters.mnn,
                     block=uncorrected$batch, direction="up", lfc=1,
                     row.data=rowData(uncorrected)[,3,drop=FALSE])
# A (probably activated?) T cell subtype of some sort:
demo <- m.out[["10"]]

plotExpression(uncorrected, x=I(factor(clusters.mnn)),
               features="ENSG00000177954", colour_by="batch") +
  facet_wrap(~colour_by)
