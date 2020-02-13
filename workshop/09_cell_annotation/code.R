# Illustrative dataset ---------------------------------------------------------

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
g <- buildSNNGraph(sce.pbmc, k=10, use.dimred="PCA")
clust <- igraph::cluster_walktrap(g)$membership
sce.pbmc$cluster <- factor(clust)

# Using the in-build references ------------------------------------------------

library(SingleR)
ref <- BlueprintEncodeData()

pred <- SingleR(test=sce.pbmc, ref=ref,
                labels=ref$label.main)
table(pred$labels)

plotScoreHeatmap(pred)

# Label pruning ----------------------------------------------------------------

sum(is.na(pred$pruned.labels))
plotScoreHeatmap(pred, show.pruned=TRUE)

plotScoreDistribution(pred)

# Identifying genes driving annotation -----------------------------------------

sce.pbmc$labels <- pred$labels
all.markers <- metadata(pred)$de.genes
lab <- "B-cells"
# Get top-10 marker genes for B-cells compared to each other cell
# type
top.markers <- Reduce(union, sapply(all.markers[[lab]], head, 10))

plotHeatmap(sce.pbmc, order_columns_by="labels",
            features=top.markers, center=TRUE, zlim=c(-3, 3), main=lab)

# Comparing to clustering ------------------------------------------------------

plotScoreHeatmap(pred, clusters = sce.pbmc$cluster, order.by.clusters = TRUE)


tab <- table(Assigned=pred$pruned.labels, Cluster=sce.pbmc$cluster)

library(pheatmap)
# Proportion of cells in each cluster assigned to each label
pheatmap(prop.table(tab, margin=2),
         color=colorRampPalette(c("white", "blue"))(101))
# (log-)number of cells in each cluster assigned to each label
# Adding a pseudo-count of 10 to avoid strong color jumps with just
# 1 cell.
pheatmap(log2(tab+10),
         color=colorRampPalette(c("white", "blue"))(101))

# Additional material ----------------------------------------------------------

# NOTE: I could use this for 'Assigning cluster labels from markers' and
#       'Computing gene set activities'. This uses the PBMC 4k dataset rather
#       than introducing additional datasets. However, these examples are
#       somewhat more involved (i.e. fewer simple function calls), so may not
#       be great in a workshop setting.

markers.pbmc <- findMarkers(sce.pbmc, sce.pbmc$cluster, direction="up", lfc=1)

chosen <- "2"
cur.markers <- markers.pbmc[[chosen]]
is.de <- cur.markers$FDR <= 0.05
summary(is.de)

# goana() requires Entrez IDs, some of which map to multiple
# symbols - hence the unique() in the call below.
library(org.Hs.eg.db)
entrez.ids <- mapIds(org.Hs.eg.db, keys=rownames(cur.markers),
                     column="ENTREZID", keytype="SYMBOL")

library(limma)
go.out <- goana(unique(entrez.ids[is.de]), species="Hs",
                universe=unique(entrez.ids))

# Only keeping biological process terms that are not overly general.
go.out <- go.out[order(go.out$P.DE),]
go.useful <- go.out[go.out$Ont=="BP" & go.out$N <= 200,]
head(go.useful, 20)

plotExpression(sce.pbmc, features=c("CD79A", "CD79B"),
               x="cluster", colour_by="cluster")

# Extract symbols for each GO term; done once.
tab <- select(org.Hs.eg.db, keytype="SYMBOL",
              keys=rownames(sce.pbmc), columns="GOALL")
by.go <- split(tab[,1], tab[,2])

# Identify genes associated with an interesting term.
b_cell_differentiation <- unique(by.go[["GO:0030183"]])
head(cur.markers[rownames(cur.markers) %in% b_cell_differentiation,1:4], 10)

aggregated <- sumCountsAcrossFeatures(sce.pbmc, by.go,
                                      exprs_values="logcounts", average=TRUE)
dim(aggregated) # rows are gene sets, columns are cells


plotColData(sce.pbmc, y=I(aggregated["GO:0030183",]), x="cluster")
