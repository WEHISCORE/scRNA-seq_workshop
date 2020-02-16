
# Application to the PBMC data -------------------------------------------------

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
