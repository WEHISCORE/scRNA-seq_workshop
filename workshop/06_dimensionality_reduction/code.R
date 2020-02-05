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

# Illustrative dataset: Ulfiltered 10X PBMC4k ----------------------------------

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

# Principal component analysis -------------------------------------------------

library(scran)
top.zeisel <- getTopHVGs(dec.zeisel, n=2000)

library(scater)
set.seed(100)
sce.zeisel <- runPCA(sce.zeisel, subset_row=top.zeisel)

# Using the elbow point --------------------------------------------------------

percent.var <- attr(reducedDim(sce.zeisel), "percentVar")
chosen.elbow <- PCAtools::findElbowPoint(percent.var)
plot(percent.var, xlab="PC",
     ylab="Variance explained (%)")
abline(v=chosen.elbow, col="red")

# Based on population structure ------------------------------------------------

choices <- getClusteredPCs(reducedDim(sce.zeisel))
choices
chosen.clusters <- metadata(choices)$chosen

# Putting it all together ------------------------------------------------------

set.seed(100)
sce.zeisel <- runPCA(sce.zeisel, subset_row=top.zeisel)
# Select d (e.g., using the elbow point or based on population
# structure
reducedDim(sce.zeisel, "PCA") <- reducedDim(
  sce.zeisel, "PCA")[,1:chosen.elbow]

# If you also want to keep the 'full' set of PCs
set.seed(100)
sce.zeisel <- runPCA(sce.zeisel, subset_row=top.zeisel)
reducedDim(sce.zeisel, "PCA_elbow") <- reducedDim(
  sce.zeisel, "PCA")[,1:chosen.elbow]
reducedDim(sce.zeisel, "PCA_clusters") <- reducedDim(
  sce.zeisel, "PCA")[,1:chosen.clusters]

# Using the technical noise ----------------------------------------------------

library(scran)
set.seed(111001001)
denoised.pbmc <- denoisePCA(sce.pbmc, technical=dec.pbmc,
                            subset.row=top.pbmc)

dim(reducedDim(denoised.pbmc, "PCA"))

set.seed(001001001)
denoised.zeisel <- denoisePCA(sce.zeisel, technical=dec.zeisel, subset.row=top.zeisel)
dim(reducedDim(denoised.zeisel))

dec.pbmc2 <- modelGeneVar(sce.pbmc)
denoised.pbmc2 <- denoisePCA(sce.pbmc,
                             technical=dec.pbmc2, subset.row=top.pbmc)
dim(reducedDim(denoised.pbmc2))

# Visualizing with PCA ---------------------------------------------------------

plotReducedDim(sce.zeisel, dimred="PCA")
plotReducedDim(sce.zeisel, dimred="PCA", colour_by="level1class")

# Interactive visualization with PCA -------------------------------------------

library(iSEE)
iSEE(sce.zeisel,
     initialPanels=DataFrame(
       Name="Reduced dimension plot 1"))

# Problems with visualizing by PCA ---------------------------------------------

plotReducedDim(sce.zeisel, dimred="PCA",
               ncomponents=4, colour_by="level1class")

# Visualizing with t-SNE -------------------------------------------------------

set.seed(00101001101)
sce.zeisel <- runTSNE(sce.zeisel, dimred="PCA")
plotReducedDim(sce.zeisel, dimred="TSNE",
               colour_by="level1class")

# Visualizing with UMAP --------------------------------------------------------

set.seed(1100101001)
sce.zeisel <- runUMAP(sce.zeisel, dimred="PCA")
plotReducedDim(sce.zeisel, dimred="UMAP",
               colour_by="level1class")
