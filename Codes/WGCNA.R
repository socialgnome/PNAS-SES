#### This file is a script that is responsible for WGCNA clustering of the whole genome 

#### Prerequisites: 1. ExpressionSet_Creation.R

#### Follow up files: 1. WGCNA_Analysis.R

## Setup env

library(here)
dev.off()
source(here("R/ThirdSubmission/Final/startup.R"))
set.seed(1234)

filt = readRDS(str_c(here("data/"), "Filtered_ExpressionSet_Clustering_AllSubData.rds"))
signatures = readRDS(str_c(here("data/"), "Signatures.rds"))
exprs = exprs(filt)

## WGCNA
enableWGCNAThreads()
# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 11, to=30, by=1))
# Call the network topology analysis function
sft = pickSoftThreshold(t(exprs), powerVector = powers, verbose = 5)
# Plot the results:
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

pow = sft$powerEstimate
net = blockwiseModules(t(exprs), power = pow,networkType = "signed hybrid",
                       TOMType = "unsigned", minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.2,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,
                       saveTOMFileBase = "AllGenes", 
                       verbose = 3)

moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
moduleLabels = data.frame(Gene = names(moduleLabels), Module = as.numeric(moduleLabels))
geneTree = net$dendrograms[[1]]
MEs0 = moduleEigengenes(t(exprs), moduleColors)$eigengenes
MEs = orderMEs(MEs0)
save.image(str_c(here("Res/ThirdSubmission/"), "WGCNA_Res.RData"))