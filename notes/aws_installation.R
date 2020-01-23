# Packages ---------------------------------------------------------------------

install.packages("BiocManager")

# TODO: Fully use renv (i.e. not just for finding dependencies)?
pkgs <- renv::dependencies()[, "Package"]
BiocManager::install(pkgs)

# Datasets ---------------------------------------------------------------------

# QC
library(scRNAseq)
LunSpikeInData(which="416b")
GrunPancreasData()
library(AnnotationHub)
AnnotationHub()[["AH73905"]]
