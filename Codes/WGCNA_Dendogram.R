#### This file is a script that is responsible for the Analysis of  WGCNA clustering 
#### of the whole genome and to test for significant clusters

#### Prerequisites: 1. WGCNA.R

#### Follow up files: 1. Fig1B_Prep.R and 2.Clustering_Reactome.R

## Setup env

library(here)
dev.off()
source(here("R/ThirdSubmission/Final/startup.R"))

## Load clustering results
load(str_c(here("Res/ThirdSubmission/"), "WGCNA_Res.RData"))
mergedColors = labels2colors(net$colors)
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)


