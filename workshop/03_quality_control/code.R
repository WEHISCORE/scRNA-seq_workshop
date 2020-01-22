# Load the data ----------------------------------------------------------------

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

library(iSEE)
iSEE(sce.416b)

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
qc.nexprs2 <- isOutlier(sce$detected, log=TRUE, type="lower")
qc.spike2 <- isOutlier(sce$altexps_ERCC_percent, type="higher")
qc.mito2 <- isOutlier(sce$subsets_Mito_percent, type="higher")
discard2 <- qc.lib2 | qc.nexprs2 | qc.spike2 | qc.mito2

# Extract the thresholds
attr(qc.lib2, "thresholds")
attr(qc.nexprs2, "thresholds")

# Summarize the number of cells removed for each reason.
DataFrame(LibSize=sum(qc.lib2), NExprs=sum(qc.nexprs2),
          SpikeProp=sum(qc.spike2), MitoProp=sum(qc.mito2), Total=sum(discard2))














