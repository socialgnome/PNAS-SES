#### This file is a script that is responsible for enrichment of the  
#### results in the significant clusters

#### Prerequisites: 1. WGCNA_Analysis.R

#### Follow up files: None

## Setup env
library(here)
dev.off()
source(here("R/ThirdSubmission/Final/startup.R"))


## Read clustering results
res = readRDS(str_c(here("Res/ThirdSubmission/Clustering"), "/1KI.rds"))

treatment = unlist(res[[7]])
signature = unlist(res[[8]])
up_clus = res[[1]]
up_de = res[[5]]
down_clus = res[[2]]
down_de = res[[6]]
up_gene_list = list()
down_gene_list = list()

treat_up = c()
geneset_up = c()
cluster_up = c()
treat_down = c()
geneset_down = c()
cluster_down = c()

for (i in 1:length(treatment)) {
  sel = which(up_de[[i]]$adj.P.Val<0.05)
  if (length(sel)>0) {
    for (j in 1:length(sel)) {
      up_gene_list[[length(up_gene_list)+1]] = up_clus[[i]][[sel[j]]]
      treat_up = c(treat_up, treatment[[i]])
      geneset_up = c(geneset_up, signature[[i]])
      cluster_up = c(cluster_up, up_de[[i]]$gene[sel[j]])
    }
  }
  sel = which(down_de[[i]]$adj.P.Val<0.05)
  if (length(sel)>0) {
    for (j in 1:length(sel)) {
      down_gene_list[[length(down_gene_list)+1]] = down_clus[[i]][[sel[j]]]
      treat_down = c(treat_down, treatment[[i]])
      geneset_down = c(geneset_down, signature[[i]])
      cluster_down = c(cluster_down, up_de[[i]]$gene[sel[j]])
    }
  }
}


names(up_gene_list) = paste(geneset_up, treat_up, cluster_up,"Upregulated", sep = "--")
names(down_gene_list) = paste(geneset_down, treat_down, cluster_down,"Downregulated", sep = "--")

glists = c(up_gene_list, down_gene_list)
## Biomart (Downloaded 18.01.22)
sapiens_ensembl <- readRDS(str_c(here("data/"), "sapiens_ensembl.rds"))

glists_ent = list()
conv = list()
for (i in 1:length(glists)) {
  df = biomaRt::getBM(attributes = c("hgnc_symbol", "entrezgene_id"), filters = "hgnc_symbol", values =glists[[i]], mart = sapiens_ensembl, useCache = F)
  conv[[i]] = df
  df = as.character(df$entrezgene_id[!is.na(df$entrezgene_id)])
  glists_ent[[i]] = df
}

clusternames = names(glists)

glists_ent = setNames(glists_ent, names(glists))
eres = compareCluster(glists_ent, fun = "enrichPathway", organism = "human", pvalueCutoff = 0.05, pAdjustMethod = "BH")
eres = eres@compareClusterResult

eres$genes = 0
for (i in 1:length(eres$Cluster)) {
  eres$genes[[i]] = paste(conv[[which(clusternames==eres$Cluster[[i]])]]$hgnc_symbol[match(unlist(strsplit(eres$geneID[[i]], "/")),conv[[which(clusternames==eres$Cluster[[i]])]]$entrezgene_id)],collapse = ",")
}

eres$Cramer = 0
for (i in 1:length(eres$Cluster)) {
  n1 = as.numeric(strsplit(eres$GeneRatio[[i]], "/")[[1]][[1]])
  n2 = as.numeric(strsplit(eres$GeneRatio[[i]], "/")[[1]][[2]])
  n3 = as.numeric(strsplit(eres$BgRatio[[i]], "/")[[1]][[1]])
  n4 = as.numeric(strsplit(eres$BgRatio[[i]], "/")[[1]][[2]])
  dd = matrix(c(n1,n3-n1,n2-n1,n4-n2-n3+n1),nrow = 2,ncol = 2)
  eres$Cramer[i] = cramerV(dd)
}


relation = read.table("~/Projects/Helper files/Reactome/ReactomePathwaysRelation.txt", sep = "\t", stringsAsFactors = F)
names = read.table("~/Projects/Helper files/Reactome/ReactomePathways.txt", sep = "\t", stringsAsFactors = F, quote = "")
names = names[names$V3=="Homo sapiens",]
test = eres
Des = data.frame(Result = test$Description, stringsAsFactors = F)
Des$Level1 = 0
Des$Level2 = 0
Des$Level3 = 0

for (i in 1:length(test$ID)) {
  
  k3 = tryCatch(which(relation$V2 == test$ID[i]),error = function(e) NULL)
  if (length(k3)>0) {
    
    k2 = tryCatch(relation$V1[k3], error = function(e) NULL)
    k3 = tryCatch(relation$V2[k3], error = function(e) NULL)
    
    k1 = tryCatch(which(relation$V2 == k2[1]), error = function(e) NULL)
    k1 = tryCatch(relation$V1[k1], error = function(e) NULL)
    
    k0 = tryCatch(which(relation$V2 == k1[1]), error = function(e) NULL)
    while (length(k0)>0) {
      k3= k2
      k2 = k1
      k1 = tryCatch(relation$V1[k0], error = function(e) NULL)
      k0 = tryCatch(which(relation$V2==k1[1]), error = function(e) NULL)
    }  
    k = c(k1,k2,k3)
    k = k[complete.cases(k)]
    
    
    if (length(names$V2[which(names$V1==k[1])])>0) {
      Des$Level1[i] = names$V2[which(names$V1==k[1])]
    }
    if (length(names$V2[which(names$V1==k[2])])>0) {
      Des$Level2[i] = names$V2[which(names$V1==k[2])]
    }
    if (length(names$V2[which(names$V1==k[3])])>0) {
      Des$Level3[i] = names$V2[which(names$V1==k[3])]
    }
  } else {
    kp = which(names$V1 == test$ID[i])
    if (length(kp)>0) {
      Des$Level1[i] = names$V2[kp]
    }
  }
  
}

Des[Des==0] <- NA
#Des = Des[,-c(4)]
Des$Show = 0
Des$ID = test$ID
test = test[complete.cases(Des[,c(2)]),]
Des = Des[complete.cases(Des[,c(2)]),]
for (i in 1:length(Des$Result)) {
  k = as.character(Des[i,-c(1,4,5,6)])
  k = k[complete.cases(k)]
  Des$Show[i] = k[length(k)]
}


Final = test
Final$Show = Des$Show
Final$Level1 = Des$Level1
Final$Level2 = Des$Level2
Final$Result = Des$Result
Final = Final[order(Final$Level1, Final$p.adjust),]
Final$FinalCount = 0
Final$Terms = 0
Final$Genes = 0
Final$MedC = 0
Final$MaxC = 0
Final$GenesSymbol = 0

for (i in 1:length(Final$ID)) {
  m = which(Final$Show==Final$Show[i] & Final$Cluster==Final$Cluster[i])
  k = strsplit(Final$geneID[m],"/")
  l = unique(unlist(k))
  Final$FinalCount[i] = length(l)
  Final$Terms[i] = length(m)
  Final$Genes[i] = paste(l, collapse = "/")
  Final$MedC[i] = median(Final$Cramer[m])
  Final$MaxC[i] = max(Final$Cramer[m])
  k = strsplit(Final$genes[m],",")
  l = unique(unlist(k))
  Final$GenesSymbol[i] = paste(l, collapse = ",")
}

write.table(Final, str_c(here("Res/ThirdSubmission/Clustering/"), "Reactome_Sig_All.txt"),row.names = F, col.names = T, sep = "\t", quote = F)
