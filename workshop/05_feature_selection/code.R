# Illustrative dataset: Unfiltered 10X PBMCs -----------------------------------

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
stats <- perCellQCMetrics(sce.pbmc, subsets=list(Mito=which(location=="MT")))
high.mito <- isOutlier(stats$subsets_Mito_percent,
                       type="higher")
sce.pbmc <- sce.pbmc[,!high.mito]

# normalization
library(scran)
set.seed(1000)
clusters <- quickCluster(sce.pbmc)
sce.pbmc <- computeSumFactors(sce.pbmc, cluster=clusters)
sce.pbmc <- logNormCounts(sce.pbmc)

# Illustrative dataset: 416B ---------------------------------------------------

library(scRNAseq)
sce.416b <- LunSpikeInData(which="416b")
sce.416b$block <- factor(sce.416b$block)

# gene-annotation
library(AnnotationHub)
ens.mm.v97 <- AnnotationHub()[["AH73905"]]
rowData(sce.416b)$ENSEMBL <- rownames(sce.416b)
rowData(sce.416b)$SYMBOL <- mapIds(ens.mm.v97,
                                   keys=rownames(sce.416b),
                                   keytype="GENEID", column="SYMBOL")
rowData(sce.416b)$SEQNAME <- mapIds(ens.mm.v97,
                                    keys=rownames(sce.416b),
                                    keytype="GENEID", column="SEQNAME")
library(scater)
rownames(sce.416b) <- uniquifyFeatureNames(
  rowData(sce.416b)$ENSEMBL,
  rowData(sce.416b)$SYMBOL)

# quality-control
mito <- which(rowData(sce.416b)$SEQNAME=="MT")
stats <- perCellQCMetrics(sce.416b, subsets=list(Mt=mito))
qc <- quickPerCellQC(
  stats,
  percent_subsets=c("subsets_Mt_percent", "altexps_ERCC_percent"),
  batch=sce.416b$block)
sce.416b <- sce.416b[,!qc$discard]

# normalization
library(scran)
sce.416b <- computeSumFactors(sce.416b)
sce.416b <- logNormCounts(sce.416b)

# Variance of the log-counts ---------------------------------------------------

library(scran)
dec.pbmc <- modelGeneVar(sce.pbmc)

# Visualizing the fit:
fit.pbmc <- metadata(dec.pbmc)
plot(fit.pbmc$mean, fit.pbmc$var,
     xlab="Mean of log-expression",
     ylab="Variance of log-expression")
curve(fit.pbmc$trend(x), col="dodgerblue", add=TRUE,
      lwd=2)

# Ordering by most interesting genes for inspection.
dec.pbmc[order(dec.pbmc$bio, decreasing=TRUE),]

# Coefficient of variation -----------------------------------------------------

dec.cv2.pbmc <- modelGeneCV2(sce.pbmc)
# Visualizing the fit:
fit.cv2.pbmc <- metadata(dec.cv2.pbmc)
plot(fit.cv2.pbmc$mean, fit.cv2.pbmc$cv2, log="xy")
curve(fit.cv2.pbmc$trend(x), col="dodgerblue", add=TRUE, lwd=2)

# Ordering by most interesting genes for inspection.
dec.cv2.pbmc[order(dec.cv2.pbmc$ratio, decreasing=TRUE),]

# In the presence of spike-ins -------------------------------------------------

dec.spike.416b <- modelGeneVarWithSpikes(sce.416b, "ERCC")
# Ordering by most interesting genes for inspection.
dec.spike.416b[order(dec.spike.416b$bio, decreasing=TRUE),]

# In the absence of spike-ins --------------------------------------------------

set.seed(0010101)
dec.pois.pbmc <- modelGeneVarByPoisson(sce.pbmc)
# Ordering by most interesting genes for inspection.
dec.pois.pbmc[order(dec.pois.pbmc$bio,decreasing=TRUE),]

# Considering experimental factors ---------------------------------------------

dec.block.416b <- modelGeneVarWithSpikes(sce.416b, "ERCC",
                                         block=sce.416b$block)
dec.block.416b[order(dec.block.416b$bio, decreasing=TRUE),]
dec.block.416b$per.block
dec.block.416b$per.block$X20160113

# Selecting HVGs on the largest metrics ----------------------------------------

# Works with modelGeneVar() output
hvg.pbmc.var <- getTopHVGs(dec.pbmc, n=1000)
str(hvg.pbmc.var)

# Works with modelGeneVarWithSpikes() output
hvg.416b.var <- getTopHVGs(dec.spike.416b, n=1000)
str(hvg.416b.var)

# Also works with modelGeneCV2() but note `var.field`
hvg.pbmc.cv2 <- getTopHVGs(dec.cv2.pbmc,
                           var.field="ratio", n=1000)
str(hvg.pbmc.cv2)

# Selecting HVGs on statistical significance -----------------------------------

# Works with modelGeneVar() output
hvg.pbmc.var.2 <- getTopHVGs(dec.pbmc, fdr.threshold=0.05)
str(hvg.pbmc.var.2)
# Works with modelGeneVarWithSpikes() output
hvg.416b.var.2 <- getTopHVGs(dec.spike.416b,
                             fdr.threshold=0.05)
str(hvg.416b.var.2)

# Also works with modelGeneCV2() but note `var.field`
hvg.pbmc.cv2.2 <- getTopHVGs(dec.cv2.pbmc,
                             var.field="ratio", fdr.threshold=0.05)
str(hvg.pbmc.cv2.2)

# Selecting genes above the trend as HVGs --------------------------------------

# Works with modelGeneVar() output
hvg.pbmc.var.3 <- getTopHVGs(dec.pbmc, var.threshold=0)
str(hvg.pbmc.var.3)
# Works with modelGeneVarWithSpikes() output
hvg.416b.var.3 <- getTopHVGs(dec.spike.416b,
                             var.threshold=0)
str(hvg.416b.var.3)

# Also works with modelGeneCV2() but note `var.field` and
# value of `var.threshold`
hvg.pbmc.cv2.3 <- getTopHVGs(dec.cv2.pbmc,
                             var.field="ratio", var.threshold=1)
str(hvg.pbmc.cv2.2)

# Putting it all together ------------------------------------------------------

dec.pbmc <- modelGeneVar(sce.pbmc)
chosen <- getTopHVGs(dec.pbmc, prop=0.1)
str(chosen)

# Subsetting to just the HVGs --------------------------------------------------

sce.pbmc.hvg <- sce.pbmc[chosen,]
sce.pbmc.hvg

# Specifying HVGs in downstream functions --------------------------------------

# Example of specifying HVGs in a downstream function
# Performing PCA only on the chosen HVGs.
library(scater)
sce.pbmc <- runPCA(sce.pbmc, subset_row=chosen)
sce.pbmc

# Witchcraft -------------------------------------------------------------------

# Add the full SCE to the subsetted data SCE
altExp(sce.pbmc.hvg, "original") <- sce.pbmc
sce.pbmc.hvg
altExp(sce.pbmc.hvg, "original")
