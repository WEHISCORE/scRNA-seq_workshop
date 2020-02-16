# Grun -------------------------------------------------------------------------

#--- loading ---#
library(scRNAseq)
sce.grun <- GrunPancreasData()

#--- gene-annotation ---#
library(org.Hs.eg.db)
gene.ids <- mapIds(org.Hs.eg.db, keys=rowData(sce.grun)$symbol,
                   keytype="SYMBOL", column="ENSEMBL")

keep <- !is.na(gene.ids) & !duplicated(gene.ids)
sce.grun <- sce.grun[keep,]
rownames(sce.grun) <- gene.ids[keep]

#--- quality-control ---#
library(scater)
stats <- perCellQCMetrics(sce.grun)

qc <- quickPerCellQC(stats, percent_subsets="altexps_ERCC_percent",
                     batch=sce.grun$donor,
                     subset=sce.grun$donor %in% c("D17", "D7", "D2"))

sce.grun <- sce.grun[,!qc$discard]

#--- normalization ---#
library(scran)
set.seed(1000) # for irlba.
clusters <- quickCluster(sce.grun)
sce.grun <- computeSumFactors(sce.grun, clusters=clusters)
sce.grun <- logNormCounts(sce.grun)

#--- variance-modelling ---#
block <- paste0(sce.grun$sample, "_", sce.grun$donor)
dec.grun <- modelGeneVarWithSpikes(sce.grun, spikes="ERCC", block=block)
top.grun <- getTopHVGs(dec.grun, prop=0.1)

# Muraro -----------------------------------------------------------------------

#--- loading ---#
library(scRNAseq)
sce.muraro <- MuraroPancreasData()

#--- gene-annotation ---#
library(AnnotationHub)
edb <- AnnotationHub()[["AH73881"]]
gene.symb <- sub("__chr.*$", "", rownames(sce.muraro))
gene.ids <- mapIds(edb, keys=gene.symb,
                   keytype="SYMBOL", column="GENEID")

# Removing duplicated genes or genes without Ensembl IDs.
keep <- !is.na(gene.ids) & !duplicated(gene.ids)
sce.muraro <- sce.muraro[keep,]
rownames(sce.muraro) <- gene.ids[keep]

#--- quality-control ---#
library(scater)
stats <- perCellQCMetrics(sce.muraro)
qc <- quickPerCellQC(stats, percent_subsets="altexps_ERCC_percent",
                     batch=sce.muraro$donor, subset=sce.muraro$donor!="D28")
sce.muraro <- sce.muraro[,!qc$discard]

#--- normalization ---#
library(scran)
set.seed(1000)
clusters <- quickCluster(sce.muraro)
sce.muraro <- computeSumFactors(sce.muraro, clusters=clusters)
sce.muraro <- logNormCounts(sce.muraro)

#--- variance-modelling ---#
block <- paste0(sce.muraro$plate, "_", sce.muraro$donor)
dec.muraro <- modelGeneVarWithSpikes(sce.muraro, "ERCC", block=block)
top.muraro <- getTopHVGs(dec.muraro, prop=0.1)

universe <- intersect(rownames(sce.grun), rownames(sce.muraro))
sce.grun2 <- sce.grun[universe,]
dec.grun2 <- dec.grun[universe,]
sce.muraro2 <- sce.muraro[universe,]
dec.muraro2 <- dec.muraro[universe,]

# Lawlor -----------------------------------------------------------------------

#--- loading ---#
library(scRNAseq)
sce.lawlor <- LawlorPancreasData()

#--- gene-annotation ---#
library(AnnotationHub)
edb <- AnnotationHub()[["AH73881"]]
anno <- select(edb, keys=rownames(sce.lawlor), keytype="GENEID",
               columns=c("SYMBOL", "SEQNAME"))
rowData(sce.lawlor) <- anno[match(rownames(sce.lawlor), anno[,1]),-1]

#--- quality-control ---#
library(scater)
stats <- perCellQCMetrics(sce.lawlor,
                          subsets=list(Mito=which(rowData(sce.lawlor)$SEQNAME=="MT")))
qc <- quickPerCellQC(stats, percent_subsets="subsets_Mito_percent",
                     batch=sce.lawlor$`islet unos id`)
sce.lawlor <- sce.lawlor[,!qc$discard]

#--- normalization ---#
library(scran)
set.seed(1000)
clusters <- quickCluster(sce.lawlor)
sce.lawlor <- computeSumFactors(sce.lawlor, clusters=clusters)
sce.lawlor <- logNormCounts(sce.lawlor)

#--- variance-modelling ---#
dec.lawlor <- modelGeneVar(sce.lawlor, block=sce.lawlor$`islet unos id`)
chosen.genes <- getTopHVGs(dec.lawlor, n=2000)

# Seger ------------------------------------------------------------------------

#--- loading ---#
library(scRNAseq)
sce.seger <- SegerstolpePancreasData()

#--- gene-annotation ---#
library(AnnotationHub)
edb <- AnnotationHub()[["AH73881"]]
symbols <- rowData(sce.seger)$symbol
ens.id <- mapIds(edb, keys=symbols, keytype="SYMBOL", column="GENEID")
ens.id <- ifelse(is.na(ens.id), symbols, ens.id)

# Removing duplicated rows.
keep <- !duplicated(ens.id)
sce.seger <- sce.seger[keep,]
rownames(sce.seger) <- ens.id[keep]

#--- sample-annotation ---#
emtab.meta <- colData(sce.seger)[,c("cell type",
                                    "individual", "single cell well quality")]
colnames(emtab.meta) <- c("CellType", "Donor", "Quality")
colData(sce.seger) <- emtab.meta

sce.seger$CellType <- gsub(" cell", "", sce.seger$CellType)
sce.seger$CellType <- paste0(
  toupper(substr(sce.seger$CellType, 1, 1)),
  substring(sce.seger$CellType, 2))

#--- quality-control ---#
low.qual <- sce.seger$Quality == "low quality cell"

library(scater)
stats <- perCellQCMetrics(sce.seger)
qc <- quickPerCellQC(stats, percent_subsets="altexps_ERCC_percent",
                     batch=sce.seger$Donor,
                     subset=!sce.seger$Donor %in% c("HP1504901", "HP1509101"))

sce.seger <- sce.seger[,!(qc$discard | low.qual)]

#--- normalization ---#
library(scran)
clusters <- quickCluster(sce.seger)
sce.seger <- computeSumFactors(sce.seger, clusters=clusters)
sce.seger <- logNormCounts(sce.seger)

#--- variance-modelling ---#
for.hvg <- sce.seger[,librarySizeFactors(altExp(sce.seger)) > 0
                     & sce.seger$Donor!="AZ"]
dec.seger <- modelGeneVarWithSpikes(for.hvg, "ERCC", block=for.hvg$Donor)
chosen.hvgs <- getTopHVGs(dec.seger, n=2000)

# Application to pancreas data: The Good ---------------------------------------

library(batchelor)
normed.pancreas <- multiBatchNorm(sce.grun2, sce.muraro2)
sce.grun2 <- normed.pancreas[[1]]
sce.muraro2 <- normed.pancreas[[2]]

library(scran)
combined.pan <- combineVar(dec.grun2, dec.muraro2)
chosen.genes <- rownames(combined.pan)[combined.pan$bio > 0]

library(scater)
rescaled.pancreas <- rescaleBatches(sce.grun2, sce.muraro2)
# NOTE: Make batch variable a factor
rescaled.pancreas$batch <- factor(
  ifelse(rescaled.pancreas$batch == 1, "Grun", "Muraro"))

set.seed(100101)
rescaled.pancreas <- runPCA(rescaled.pancreas, subset_row=chosen.genes,
                            exprs_values="corrected")

rescaled.pancreas <- runTSNE(rescaled.pancreas, dimred="PCA")
plotTSNE(rescaled.pancreas, colour_by="batch")

set.seed(1011011)
mnn.pancreas <- fastMNN(sce.grun2, sce.muraro2, subset.row=chosen.genes)
mnn.pancreas$batch <- factor(
  ifelse(mnn.pancreas$batch == 1, "Grun", "Muraro"))

snn.gr <- buildSNNGraph(mnn.pancreas, use.dimred="corrected")
clusters <- igraph::cluster_walktrap(snn.gr)$membership
tab <- table(Cluster=clusters, Batch=mnn.pancreas$batch)
tab

mnn.pancreas <- runTSNE(mnn.pancreas, dimred="corrected")
plotTSNE(mnn.pancreas, colour_by="batch")

# My summary plot
uncorrected.pancreas <- noCorrect(sce.grun2, sce.muraro2)
set.seed(1001011)
uncorrected.pancreas <- runPCA(uncorrected.pancreas, subset_row=chosen.genes,
                            exprs_values="merged")
uncorrected.pancreas <- runTSNE(uncorrected.pancreas, dimred="PCA")
uncorrected.pancreas$batch <- factor(
  ifelse(uncorrected.pancreas$batch == 1, "Grun", "Muraro"))
cowplot::plot_grid(
  plotTSNE(uncorrected.pancreas, colour_by="batch") +
    ggtitle("Uncorrected") +
    cowplot::theme_cowplot(font_size = 20),
  plotTSNE(rescaled.pancreas, colour_by="batch") +
    ggtitle("Linear regression") +
    cowplot::theme_cowplot(font_size = 20),
  plotTSNE(mnn.pancreas, colour_by="batch") +
    ggtitle("MNN") +
    cowplot::theme_cowplot(font_size = 20),
  ncol = 3)

# Application to pancreas data: The Bad ----------------------------------------

# NOTE: These results look different to OSCA.

all.sce <- list(Grun=sce.grun, Muraro=sce.muraro,
                Lawlor=sce.lawlor, Seger=sce.seger)
all.dec <- list(Grun=dec.grun, Muraro=dec.muraro,
                Lawlor=dec.lawlor, Seger=dec.seger)

universe <- Reduce(intersect, lapply(all.sce, rownames))
all.sce <- lapply(all.sce, "[", i=universe,)
all.dec <- lapply(all.dec, "[", i=universe,)

normed.pancreas <- do.call(multiBatchNorm, all.sce)
combined.pan <- do.call(combineVar, all.dec)
chosen.genes <- rownames(combined.pan)[combined.pan$bio > 0]

set.seed(1011110)
# mnn.pancreas <- fastMNN(normed.pancreas)
mnn.pancreas <- fastMNN(
  normed.pancreas$Grun, normed.pancreas$Muraro, normed.pancreas$Lawlor, normed.pancreas$Seger)
# NOTE: Make batch variable a factor
mnn.pancreas$batch <- factor(
  dplyr::case_when(
    mnn.pancreas$batch == 1 ~ "Grun",
    mnn.pancreas$batch == 2 ~ "Muraro",
    mnn.pancreas$batch == 3 ~ "Lawlor",
    mnn.pancreas$batch == 4 ~ "Seger"))

# Bumping up 'k' to get broader clusters for this demonstration.
snn.gr <- buildSNNGraph(mnn.pancreas, use.dimred="corrected", k=20)
clusters <- igraph::cluster_walktrap(snn.gr)$membership
clusters <- factor(clusters)
tab <- table(Cluster=clusters, Batch=mnn.pancreas$batch)
tab

mnn.pancreas <- runTSNE(mnn.pancreas, dimred="corrected")
gridExtra::grid.arrange(
  plotTSNE(mnn.pancreas, colour_by="batch", text_by=I(clusters)),
  plotTSNE(mnn.pancreas, colour_by=I(clusters), text_by=I(clusters)),
  ncol=2
)

donors <- c(
  normed.pancreas$Grun$donor,
  normed.pancreas$Muraro$donor,
  normed.pancreas$Lawlor$`islet unos id`,
  normed.pancreas$Seger$Donor
)

seger.donors <- donors
seger.donors[mnn.pancreas$batch!="Seger"] <- NA
plotTSNE(mnn.pancreas, colour_by=I(seger.donors))

# Application to pancreas data: The Ugly ---------------------------------------

# combined <- noCorrect(normed.pancreas)
combined <- do.call(noCorrect, normed.pancreas)
assayNames(combined) <- "logcounts"
combined$donor <- donors

donors.by.batches <- lapply(split(combined$donor, combined$batch), unique)
donors.by.batches

ndonors <- rep(lengths(donors.by.batches), lengths(donors.by.batches))
names(ndonors) <- unlist(donors.by.batches)
ndonors

set.seed(1010100)
multiout <- fastMNN(combined, batch=donors, subset.row=chosen.genes,
                    weights=1/ndonors)

# Renaming metadata fields for easier communication later.
multiout$dataset <- combined$batch
multiout$donor <- multiout$batch
multiout$batch <- NULL

library(scater)
g <- buildSNNGraph(multiout, use.dimred=1, k=20)
clusters <- igraph::cluster_walktrap(g)$membership
tab <- table(clusters, multiout$dataset)
tab

multiout <- runTSNE(multiout, dimred="corrected")
gridExtra::grid.arrange(
  plotTSNE(multiout, colour_by="dataset", text_by=I(clusters)),
  plotTSNE(multiout, colour_by=data.frame(donor = seger.donors)),
  ncol=2
)
