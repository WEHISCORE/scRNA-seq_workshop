# Illustrative dataset: Unfiltered 10X PBMC4k ----------------------------------

library(BiocFileCache)
bfc <- BiocFileCache()
raw.path <- bfcrpath(bfc, file.path("http://cf.10xgenomics.com/samples",
                                    "cell-exp/2.1.0/pbmc4k/pbmc4k_raw_gene_bc_matrices.tar.gz"))
untar(raw.path, exdir=file.path(tempdir(), "pbmc4k"))

library(DropletUtils)
library(Matrix)
fname <- file.path(tempdir(), "pbmc4k/raw_gene_bc_matrices/GRCh38")
sce.pbmc <- read10xCounts(fname, col.names=TRUE)

# gene-annotation
library(scater)
rownames(sce.pbmc) <- uniquifyFeatureNames(
  rowData(sce.pbmc)$ID, rowData(sce.pbmc)$Symbol)
library(EnsDb.Hsapiens.v86)
location <- mapIds(EnsDb.Hsapiens.v86, keys=rowData(sce.pbmc)$ID,
                   column="SEQNAME", keytype="GENEID")

# cell-detection
set.seed(100)
e.out <- emptyDrops(counts(sce.pbmc))
sce.pbmc <- sce.pbmc[,which(e.out$FDR <= 0.001)]

# quality-control
stats <- perCellQCMetrics(sce.pbmc,
                          subsets=list(Mito=which(location=="MT")))
high.mito <- isOutlier(stats$subsets_Mito_percent,
                       type="higher")
sce.pbmc <- sce.pbmc[,!high.mito]

# normalization
library(scran)
set.seed(1000)
clusters <- quickCluster(sce.pbmc)
sce.pbmc <- computeSumFactors(sce.pbmc, cluster=clusters)
sce.pbmc <- logNormCounts(sce.pbmc)

# variance modelling
set.seed(1001)
dec.pbmc <- modelGeneVarByPoisson(sce.pbmc)
top.pbmc <- getTopHVGs(dec.pbmc, prop=0.1)

# dimensionality-reduction
set.seed(10000)
sce.pbmc <- denoisePCA(sce.pbmc, subset.row=top.pbmc,
                       technical=dec.pbmc)

set.seed(100000)
sce.pbmc <- runTSNE(sce.pbmc, dimred="PCA")

set.seed(1000000)
sce.pbmc <- runUMAP(sce.pbmc, dimred="PCA")

# clustering
g <- buildSNNGraph(sce.pbmc, k=10, use.dimred = "PCA")
clust <- igraph::cluster_walktrap(g)$membership
sce.pbmc$cluster <- factor(clust)

# Motivation -------------------------------------------------------------------

# Is the first gene associated with the clustering?
plotExpression(sce.pbmc, features=rownames(sce.pbmc)[1],
               x="cluster", colour_by="cluster")

# Is the second gene associated with the clustering?
plotExpression(sce.pbmc, features=rownames(sce.pbmc)[2],
               x="cluster", colour_by="cluster")

# Is gene 16863 associated with the clustering?
plotExpression(sce.pbmc, features=rownames(sce.pbmc)[16863],
               x="cluster", colour_by="cluster")

# Standard application ---------------------------------------------------------

library(scran)
markers.pbmc <- findMarkers(sce.pbmc, groups=sce.pbmc$cluster,
                            test.type="t", pval.type="any")

chosen <- "9"
interesting <- markers.pbmc[[chosen]]

plotExpression(sce.pbmc, rownames(interesting)[1:4],
               x="cluster", colour_by="cluster")

best.set <- interesting[interesting$Top <= 6,]
logFCs <- as.matrix(best.set[,-(1:3)])
colnames(logFCs) <- sub("logFC.", "", colnames(logFCs))

library(pheatmap)
pheatmap(logFCs, breaks=seq(-5, 5, length.out=101))

# Using the log-fold change ----------------------------------------------------

markers.pbmc.up <- findMarkers(sce.pbmc, groups=sce.pbmc$cluster,
                               test.type="t", direction="up", pval.type="any")
interesting.up <- markers.pbmc.up[[chosen]]


markers.pbmc.up2 <- findMarkers(sce.pbmc, groups=sce.pbmc$cluster,
                                test.type="t", direction="up", lfc=1, pval.type="any")
interesting.up2 <- markers.pbmc.up2[[chosen]]


best.set <- interesting.up2[interesting.up2$Top <= 5,]
logFCs <- as.matrix(best.set[,-(1:3)])
colnames(logFCs) <- sub("logFC.", "", colnames(logFCs))

pheatmap(logFCs, breaks=seq(-5, 5, length.out=101))

# Finding cluster-specific marker genes ----------------------------------------

markers.pbmc.up3 <- findMarkers(sce.pbmc, groups=sce.pbmc$cluster,
                                direction="up", pval.type="all")
interesting.up3 <- markers.pbmc.up3[[chosen]]


# Finding cluster-specific-ish marker genes ------------------------------------

markers.pbmc.up4 <- findMarkers(sce.pbmc, groups=sce.pbmc$cluster,
                                direction="up", pval.type="some")
interesting.up4 <- markers.pbmc.up4[[chosen]]

# Wilcoxon rank sum test -------------------------------------------------------

markers.pbmc.wmw <- findMarkers(sce.pbmc,
                                groups=sce.pbmc$cluster, test.type="wilcox",
                                direction="up", pval.type="any")
interesting.wmw <- markers.pbmc.wmw[[chosen]]


best.set <- interesting.wmw[interesting.wmw$Top <= 5,]
AUCs <- as.matrix(best.set[,-(1:3)])
colnames(AUCs) <- sub("AUC.", "", colnames(AUCs))

pheatmap(AUCs, breaks=seq(0, 1, length.out=21),
         color=viridis::viridis(21))

# Binomial test ----------------------------------------------------------------

markers.pbmc.binom <- findMarkers(sce.pbmc,
                                  groups=sce.pbmc$cluster, test.type="binom",
                                  direction="up", pval.type="any")
interesting.binom <- markers.pbmc.binom[[chosen]]


top.genes <- head(rownames(interesting.binom))
plotExpression(sce.pbmc, x="cluster", features=top.genes)


