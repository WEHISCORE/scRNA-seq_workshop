# Interactive visualization of Day1 content ------------------------------------

# TODO: Would be good to have gene symbols rather than Ensembl IDs.
sce.zeisel <- addPerCellQC(sce.zeisel, subsets = list(Mito = rowData(sce.zeisel)$featureType == "mito"))

# TODO: Would be good to have some interesting genes to highlight
library(iSEE)
iSEE(sce.zeisel,
     initialPanels=DataFrame(
       Name=c("Reduced dimension plot 1", "Feature assay plot 1", "Column data plot 1")))
