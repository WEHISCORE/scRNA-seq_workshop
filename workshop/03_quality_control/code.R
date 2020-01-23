# Illustrative dataset: 416B ---------------------------------------------------

library(scRNAseq)
sce.416b <- LunSpikeInData(which="416b")
sce.416b$block <- factor(sce.416b$block)

# Identifying the mitochondrial transcripts ------------------------------------

library(AnnotationHub)
ens.mm.v97 <- AnnotationHub()[["AH73905"]]
location <- mapIds(ens.mm.v97, keys=rownames(sce.416b),
                   keytype="GENEID", column="SEQNAME")
is.mito <- which(location=="MT")

# Computing the QC metrics -----------------------------------------------------

library(scater)
sce.416b <- addPerCellQC(sce.416b, subsets=list(Mito=is.mito))

# Visualizing the QC metrics ---------------------------------------------------

plotColData(sce.416b, x="block", y="detected")

plotColData(sce.416b, x="block", y="detected") +
  scale_y_log10()

plotColData(sce.416b, x="block", y="detected", other_fields="phenotype") +
  scale_y_log10() +
  facet_wrap(~phenotype)

# Using fixed thresholds -------------------------------------------------------

# Example thresholds
qc.lib <- sce.416b$sum < 100000
qc.nexprs <- sce.416b$detected < 5000
qc.spike <- sce.416b$altexps_ERCC_percent > 10
qc.mito <- sce.416b$subsets_Mito_percent > 10
discard <- qc.lib | qc.nexprs | qc.spike | qc.mito

# Summarize the number of cells removed for each reason
DataFrame(LibSize=sum(qc.lib), NExprs=sum(qc.nexprs),
          SpikeProp=sum(qc.spike), MitoProp=sum(qc.mito), Total=sum(discard))

# Using adaptive thresholds ----------------------------------------------------

# Adaptive thresholds
qc.lib2 <- isOutlier(sce.416b$sum, log=TRUE, type="lower")
qc.nexprs2 <- isOutlier(sce.416b$detected, log=TRUE, type="lower")
qc.spike2 <- isOutlier(sce.416b$altexps_ERCC_percent, type="higher")
qc.mito2 <- isOutlier(sce.416b$subsets_Mito_percent, type="higher")
discard2 <- qc.lib2 | qc.nexprs2 | qc.spike2 | qc.mito2

# Extract the thresholds
attr(qc.lib2, "thresholds")
attr(qc.nexprs2, "thresholds")

# Summarize the number of cells removed for each reason.
DataFrame(LibSize=sum(qc.lib2), NExprs=sum(qc.nexprs2),
          SpikeProp=sum(qc.spike2), MitoProp=sum(qc.mito2), Total=sum(discard2))

# Considering experimental factors ---------------------------------------------

plotColData(sce.416b, x="block", y="detected", other_fields="phenotype") +
  scale_y_log10() +
  facet_wrap(~phenotype)

# Adaptive thresholds accounting for batch
batch <- paste0(sce.416b$phenotype, "-", sce.416b$block)
qc.lib3 <- isOutlier(sce.416b$sum, log=TRUE, type="lower", batch=batch)
qc.nexprs3 <- isOutlier(sce.416b$detected, log=TRUE, type="lower", batch=batch)
qc.spike3 <- isOutlier(sce.416b$altexps_ERCC_percent, type="higher", batch=batch)
qc.mito3 <- isOutlier(sce.416b$subsets_Mito_percent, type="higher", batch=batch)
discard3 <- qc.lib3 | qc.nexprs3 | qc.spike3 | qc.mito3

# Extract the thresholds
attr(qc.lib3, "thresholds")
attr(qc.nexprs3, "thresholds")

# Illustratative dataset: Grun -------------------------------------------------

sce.grun <- GrunPancreasData()
sce.grun <- addPerCellQC(sce.grun)

# What to do if an entire batch failed -----------------------------------------

plotColData(sce.grun, x="donor", y="altexps_ERCC_percent")

# Adaptive outlier detection using a subset ------------------------------------

discard.ercc <- isOutlier(sce.grun$altexps_ERCC_percent,
                          type="higher", batch=sce.grun$donor)
discard.ercc2 <- isOutlier(sce.grun$altexps_ERCC_percent,
                           type="higher", batch=sce.grun$donor,
                           subset=sce.grun$donor %in% c("D17", "D2", "D7"))

plotColData(sce.grun, x="donor", y="altexps_ERCC_percent",
            colour_by=data.frame(discard=discard.ercc))
plotColData(sce.grun, x="donor", y="altexps_ERCC_percent",
            colour_by=data.frame(discard=discard.ercc2))

# Checking diagnostic plots ----------------------------------------------------

# Add info about which cells are outliers
sce.416b$discard <- discard2

# Make these plots for each QC metric
plotColData(sce.416b, x="block", y="sum", colour_by="discard", other_fields="phenotype") +
  facet_wrap(~phenotype) +
  scale_y_log10()

# Another useful diagnostic plot
plotColData(sce.416b, x="sum", y="subsets_Mito_percent",
            colour_by="discard", other_fields=c("block", "phenotype")) +
  facet_grid(block~phenotype)

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

# Barcode rank plot ------------------------------------------------------------

bcrank <- barcodeRanks(counts(sce.pbmc))

# Only showing unique points for plotting speed.
uniq <- !duplicated(bcrank$rank)
plot(bcrank$rank[uniq], bcrank$total[uniq], log="xy",
     xlab="Rank", ylab="Total UMI count", cex.lab=1.2)
abline(h=metadata(bcrank)$inflection, col="darkgreen", lty=2)
abline(h=metadata(bcrank)$knee, col="dodgerblue", lty=2)
legend("bottomleft", legend=c("Inflection", "Knee"),
       col=c("darkgreen", "dodgerblue"), lty=2, cex=1.2)

# Testing for empty droplets ---------------------------------------------------

set.seed(100)
e.out <- emptyDrops(counts(sce.pbmc))

# See ?emptyDrops for an explanation of why there are NA values.
summary(e.out$FDR <= 0.001)

# Removing empty droplets ------------------------------------------------------

sce.pbmc <- sce.pbmc[,which(e.out$FDR <= 0.001)]

# Relationship to other QC metrics ---------------------------------------------

is.mito <- grep("^MT-", rowData(sce.pbmc)$Symbol)
sce.pmbc <- addPerCellQC(sce.pbmc, subsets=list(MT=is.mito))
discard.mito <- isOutlier(sce.pmbc$subsets_MT_percent, type="higher")
plot(sce.pmbc$sum, sce.pmbc$subsets_MT_percent, log="x",
     xlab="Total count", ylab='Mitochondrial %')
abline(h=attr(discard.mito, "thresholds")["higher"], col="red")

# Handling low-quality cells in the 416B dataset -------------------------------

# Removing low-quality cells
# Keeping the columns we DON'T want to discard
filtered <- sce.416b[,!discard2]
# Marking low-quality cells
marked <- sce.416b
marked$discard <- discard2
