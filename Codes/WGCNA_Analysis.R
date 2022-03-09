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
rm(list= ls()[!(ls() %in% c('net'))])

## Load expression data
dat = readRDS(str_c(here("data/"), "Filtered_ExpressionSet_Clustering_SelSubData.rds"))
exprs = exprs(dat)

## Module characteristics and numbers
moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
moduleInfo = data.frame(Gene = names(moduleLabels), Module = as.numeric(moduleLabels), Colors = moduleColors)
ColInfo = moduleInfo[!duplicated(moduleInfo[-c(1)]),-c(1)]
MEs = net$MEs;
# MEs0 = moduleEigengenes(t(exprs), moduleColors)$eigengenes
# MEs1 = orderMEs(MEs0)
mod_nos = as.data.frame(table(moduleLabels))

## Subject data to test for significance
subData = pData(dat)
controls = c(
  "sex_interv", "re","age_w5"
  ,"pregnant_biow5", "FastHrs","Plate"
  ,"H5INFECT", "H5SUBCLN", "H5CRP8"
)

## 1KI analysis
signatures = readRDS("~/Projects/PNAS/data/Signatures.rds")

treatment = c("ses_sss_composite","sss_5","SEI_ff5",
              "edu_max","income_hh_ff5")

## Remove subjects with diseases
sub_rem = list()
for (i in 1:length(signatures)) {
  sub_rem[[i]] = NA
}
sub_rem[[1]] = which(subData$heartatk==1) # CVD
sub_rem[[2]] = which(subData$diabetes==1) # Diabetes
sub_rem[[6]] = which(subData$H5ID6Q == 1 | subData$H5ID6QM == 1) # Aortic Aneurysm
sub_rem[[8]] = which(subData$H5ID6F == 1 | subData$H5ID6FM == 1) # Asthma
sub_rem[[9]] = which(subData$H5ID6C == 1 | subData$H5ID6CM == 1) # Hypertension


res = data.frame(Signature = rep(names(signatures), each = length(treatment)), 
                 Treatment = rep(treatment, length(names(signatures))),
                 Total_Genes = 0, Up_Genes = 0, Down_Genes = 0,Total_Clusters = 0, Up_Clusters = 0,
                 Down_Clusters = 0, Up_Sig_Genes = 0, Down_Sig_Genes = 0, Up_Sig_Clusters = 0,
                 Down_Sig_Clusters = 0, Up_Min_FDR = 0, Up_Av_FDR = 0, Down_Min_FDR = 0, Down_Av_FDR = 0)

allres = list()
allres[[1]] = list() #Save up module information
allres[[2]] = list() #Save down module information
allres[[3]] = list() #Save up module average expression
allres[[4]] = list() #Save down module average expression
allres[[5]] = list() #Save up module average expression DE results
allres[[6]] = list() #Save up module average expression DE results
allres[[7]] = list() #Treatment
allres[[8]] = list() #Signature
naming = c() #List names
names = c("Up Gene Modules", "Down Gene Modules", "Up Av Expression", "Down Av Expression", "Up Av Expression DE", "Down Av Expression DE", "Treatment", "Signature")
allres = setNames(allres, names)

for (i in 1:length(signatures)) {
  
  for (j in 1:length(treatment)) {
    
    genes = signatures[[i]]
    res$Total_Genes[((i-1)*length(treatment))+j] = length(genes)
    
    sel = moduleInfo[moduleInfo$Gene %in% genes,]
    df = readRDS(str_c(here("Res/ThirdSubmission/DE/"), treatment[[j]],".rds"))
    df = df[df$gene %in% genes,]
    all = sel[sel$Gene %in% df$gene,]
    all_split = split(all$Gene, all$Module)
    up = sel[sel$Gene %in% df$gene[df$logFC>0],]
    down = sel[sel$Gene %in% df$gene[df$logFC<0],]
    up_split = split(up$Gene, up$Module)
    allres[[1]][[((i-1)*length(treatment))+j]] = up_split
    down_split = split(down$Gene, down$Module)
    allres[[2]][[((i-1)*length(treatment))+j]] = down_split
    
    res$Up_Genes[((i-1)*length(treatment))+j] = length(up$Gene)
    res$Down_Genes[((i-1)*length(treatment))+j] = length(down$Gene)
    res$Total_Clusters[((i-1)*length(treatment))+j] = length(all_split)
    res$Up_Clusters[((i-1)*length(treatment))+j] = length(up_split)
    res$Down_Clusters[((i-1)*length(treatment))+j] =length(down_split)
    
    av = c()
    for (k in 1:length(up_split)) {
      avg = exprs[up_split[[k]],,drop=F]
      av = rbind(av,colMeans(avg))
    }
    
    rownames(av) = names(up_split)
    rhs = str_c(c(treatment[[j]],controls), collapse = " + ")
    model_formula = str_c(" ~ ",rhs) %>% as.formula()
    sub = subData %>% dplyr::select(treatment[[j]], all_of(controls))
    sub = droplevels(sub)
    
    missing = sub_rem[[i]]
    keep = setdiff(1:length(subData$AID), missing)
    sub = sub[keep, ]
    sub = droplevels(sub)
    
    av = av[,keep]
    allres[[3]][[((i-1)*length(treatment))+j]] = av
    
    if (class(sub[,c(1)])=="factor") {
      sub[,c(1)] = as.numeric(sub[,c(1)])
    }
    
    design  = model.matrix(model_formula, data = sub)
    fit = lmFit(av, design) 
    fit = eBayes(fit, trend = T)
    tab = topTable(fit, coef = treatment[[j]], n= Inf) %>%
      rownames_to_column(var = "gene")
    allres[[5]][[((i-1)*length(treatment))+j]] = tab
    res$Up_Sig_Clusters[((i-1)*length(treatment))+j] = length(tab$gene[tab$adj.P.Val<0.05 & tab$logFC>0])
    res$Up_Sig_Genes[((i-1)*length(treatment))+j] = length(up$Gene[up$Module %in% as.numeric(tab$gene[which(tab$adj.P.Val<0.05)])])
    res$Up_Min_FDR[((i-1)*length(treatment))+j] = min(tab$adj.P.Val)
    res$Up_Av_FDR[((i-1)*length(treatment))+j] = combine.test(tab$adj.P.Val)
    
    
    av = c()
    for (k in 1:length(down_split)) {
      avg = exprs[down_split[[k]],,drop=F]
      av = rbind(av,colMeans(avg))
    }
    
    rownames(av) = names(down_split)
    rhs = str_c(c(treatment[[j]],controls), collapse = " + ")
    model_formula = str_c(" ~ ",rhs) %>% as.formula()
    sub = subData %>% dplyr::select(treatment[[j]], all_of(controls))
    sub = droplevels(sub)
    
    missing = sub_rem[[i]]
    keep = setdiff(1:length(subData$AID), missing)
    sub = sub[keep, ]
    sub = droplevels(sub)
    
    av = av[,keep]
    allres[[4]][[((i-1)*length(treatment))+j]] = av
    
    if (class(sub[,c(1)])=="factor") {
      sub[,c(1)] = as.numeric(sub[,c(1)])
    }
    
    design  = model.matrix(model_formula, data = sub)
    fit = lmFit(av, design) 
    fit = eBayes(fit, trend = T)
    tab = topTable(fit, coef = treatment[[j]], n= Inf) %>%
      rownames_to_column(var = "gene")
    allres[[6]][[((i-1)*length(treatment))+j]] = tab
    res$Down_Sig_Clusters[((i-1)*length(treatment))+j] = length(tab$gene[tab$adj.P.Val<0.05 & tab$logFC<0])
    res$Down_Sig_Genes[((i-1)*length(treatment))+j] = length(down$Gene[down$Module %in% as.numeric(tab$gene[which(tab$adj.P.Val<0.05)])])
    res$Down_Min_FDR[((i-1)*length(treatment))+j] = min(tab$adj.P.Val)
    res$Down_Av_FDR[((i-1)*length(treatment))+j] = combine.test(tab$adj.P.Val)
    
    allres[[7]][[((i-1)*length(treatment))+j]] = treatment[[j]]
    allres[[8]][[((i-1)*length(treatment))+j]] = names(signatures)[[i]]
  }
}

res$Class = "1KI Genes"
res_1k = res
saveRDS(allres, str_c(here("Res/ThirdSubmission/Clustering/"), "1KI.rds"))


## Without 1KI analysis
inflam_genes = signatures[[which(names(signatures)=="inflam1k_mRNA")]]
sub_rem = sub_rem[-c(which(names(signatures)=="inflam1k_mRNA"))]
signatures = signatures[-c(which(names(signatures)=="inflam1k_mRNA"))]

res = data.frame(Signature = rep(names(signatures), each = length(treatment)), 
                 Treatment = rep(treatment, length(names(signatures))),
                 Total_Genes = 0, Up_Genes = 0, Down_Genes = 0,Total_Clusters = 0, Up_Clusters = 0,
                 Down_Clusters = 0, Up_Sig_Genes = 0, Down_Sig_Genes = 0, Up_Sig_Clusters = 0,
                 Down_Sig_Clusters = 0, Up_Min_FDR = 0, Up_Av_FDR = 0, Down_Min_FDR = 0, Down_Av_FDR = 0)

allres = list()
allres[[1]] = list() #Save up module information
allres[[2]] = list() #Save down module information
allres[[3]] = list() #Save up module average expression
allres[[4]] = list() #Save down module average expression
allres[[5]] = list() #Save up module average expression DE results
allres[[6]] = list() #Save up module average expression DE results
allres[[7]] = list() #Treatment
allres[[8]] = list() #Signature
naming = c() #List names
names = c("Up Gene Modules", "Down Gene Modules", "Up Av Expression", "Down Av Expression", "Up Av Expression DE", "Down Av Expression DE", "Treatment", "Signature")
allres = setNames(allres, names)

for (i in 1:length(signatures)) {
  
  for (j in 1:length(treatment)) {
    
    genes = setdiff(signatures[[i]], inflam_genes)
    res$Total_Genes[((i-1)*length(treatment))+j] = length(genes)
    
    sel = moduleInfo[moduleInfo$Gene %in% genes,]
    df = readRDS(str_c(here("Res/ThirdSubmission/DE/"), treatment[[j]],".rds"))
    df = df[df$gene %in% genes,]
    all = sel[sel$Gene %in% df$gene,]
    all_split = split(all$Gene, all$Module)
    up = sel[sel$Gene %in% df$gene[df$logFC>0],]
    down = sel[sel$Gene %in% df$gene[df$logFC<0],]
    up_split = split(up$Gene, up$Module)
    allres[[1]][[((i-1)*length(treatment))+j]] = up_split
    down_split = split(down$Gene, down$Module)
    allres[[2]][[((i-1)*length(treatment))+j]] = down_split
    
    res$Up_Genes[((i-1)*length(treatment))+j] = length(up$Gene)
    res$Down_Genes[((i-1)*length(treatment))+j] = length(down$Gene)
    res$Total_Clusters[((i-1)*length(treatment))+j] = length(all_split)
    res$Up_Clusters[((i-1)*length(treatment))+j] = length(up_split)
    res$Down_Clusters[((i-1)*length(treatment))+j] =length(down_split)
    
    av = c()
    for (k in 1:length(up_split)) {
      avg = exprs[up_split[[k]],,drop=F]
      av = rbind(av,colMeans(avg))
    }
    
    rownames(av) = names(up_split)
    rhs = str_c(c(treatment[[j]],controls), collapse = " + ")
    model_formula = str_c(" ~ ",rhs) %>% as.formula()
    sub = subData %>% dplyr::select(treatment[[j]], all_of(controls))
    sub = droplevels(sub)
    
    missing = sub_rem[[i]]
    keep = setdiff(1:length(subData$AID), missing)
    sub = sub[keep, ]
    sub = droplevels(sub)
    
    av = av[,keep]
    allres[[3]][[((i-1)*length(treatment))+j]] = av
    
    if (class(sub[,c(1)])=="factor") {
      sub[,c(1)] = as.numeric(sub[,c(1)])
    }
    
    design  = model.matrix(model_formula, data = sub)
    fit = lmFit(av, design) 
    fit = eBayes(fit, trend = T)
    tab = topTable(fit, coef = treatment[[j]], n= Inf) %>%
      rownames_to_column(var = "gene")
    allres[[5]][[((i-1)*length(treatment))+j]] = tab
    res$Up_Sig_Clusters[((i-1)*length(treatment))+j] = length(tab$gene[tab$adj.P.Val<0.05 & tab$logFC>0])
    res$Up_Sig_Genes[((i-1)*length(treatment))+j] = length(up$Gene[up$Module %in% as.numeric(tab$gene[which(tab$adj.P.Val<0.05)])])
    res$Up_Min_FDR[((i-1)*length(treatment))+j] = min(tab$adj.P.Val)
    res$Up_Av_FDR[((i-1)*length(treatment))+j] = combine.test(tab$adj.P.Val)
    
    
    av = c()
    for (k in 1:length(down_split)) {
      avg = exprs[down_split[[k]],,drop=F]
      av = rbind(av,colMeans(avg))
    }
    
    rownames(av) = names(down_split)
    rhs = str_c(c(treatment[[j]],controls), collapse = " + ")
    model_formula = str_c(" ~ ",rhs) %>% as.formula()
    sub = subData %>% dplyr::select(treatment[[j]], all_of(controls))
    sub = droplevels(sub)
    
    missing = sub_rem[[i]]
    keep = setdiff(1:length(subData$AID), missing)
    sub = sub[keep, ]
    sub = droplevels(sub)
    
    av = av[,keep]
    allres[[4]][[((i-1)*length(treatment))+j]] = av
    
    if (class(sub[,c(1)])=="factor") {
      sub[,c(1)] = as.numeric(sub[,c(1)])
    }
    
    design  = model.matrix(model_formula, data = sub)
    fit = lmFit(av, design) 
    fit = eBayes(fit, trend = T)
    tab = topTable(fit, coef = treatment[[j]], n= Inf) %>%
      rownames_to_column(var = "gene")
    allres[[6]][[((i-1)*length(treatment))+j]] = tab
    res$Down_Sig_Clusters[((i-1)*length(treatment))+j] = length(tab$gene[tab$adj.P.Val<0.05 & tab$logFC<0])
    res$Down_Sig_Genes[((i-1)*length(treatment))+j] = length(down$Gene[down$Module %in% as.numeric(tab$gene[which(tab$adj.P.Val<0.05)])])
    res$Down_Min_FDR[((i-1)*length(treatment))+j] = min(tab$adj.P.Val)
    res$Down_Av_FDR[((i-1)*length(treatment))+j] = combine.test(tab$adj.P.Val)
    
    allres[[7]][[((i-1)*length(treatment))+j]] = treatment[[j]]
    allres[[8]][[((i-1)*length(treatment))+j]] = names(signatures)[[i]]
  }
}
res$Class = "Without 1KI Genes"
res_wo_1k = res
saveRDS(allres, str_c(here("Res/ThirdSubmission/Clustering/"), "Without_1KI.rds"))

## All results
res = rbind(res_1k, res_wo_1k)
saveRDS(res, str_c(here("Res/ThirdSubmission/"), "Clustering_Fig1B.rds"))