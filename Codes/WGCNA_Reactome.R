#### This file is a script that is responsible for the analysing the functional 
#### roles of each cluster that is identified by WGCNA

#### Prerequisites: 1. WGCNA.R

#### Follow up files: None

## Setup env

library(here)
dev.off()
source(here("R/ThirdSubmission/Final/startup.R"))

## Load clustering results
load(str_c(here("Res/ThirdSubmission/"), "WGCNA_Res.RData"))
rm(list= ls()[!(ls() %in% c('net'))])

## Module characteristics and numbers
moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
moduleInfo = data.frame(Gene = names(moduleLabels), Module = as.numeric(moduleLabels), Colors = moduleColors)
ColInfo = moduleInfo[!duplicated(moduleInfo[-c(1)]),-c(1)]
MEs = net$MEs;
# MEs0 = moduleEigengenes(t(exprs), moduleColors)$eigengenes
# MEs1 = orderMEs(MEs0)
mod_nos = as.data.frame(table(moduleLabels))

## Find enriched pathways
glists = split(moduleInfo$Gene, moduleInfo$Module)

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

lists = split(eres$ID,eres$Cluster)

# tiff(str_c(here("Res/ThirdSubmission/Figures/"), "WGCNA_Pathway.tiff"), units="px", width=(3*1150), height=(3*675), res=300)
# t1 = upset(fromList(lists),
#            mainbar.y.label = "Intersecting enriched pathways",
#            sets.x.label = "Total enriched pathways",
#            keep.order = T,
#            main.bar.color = "#555555",
#            matrix.color  = "#326a97",
#            set_size.show = T,
#            point.size = 2.8,
#            shade.alpha = 0.2,
#            matrix.dot.alpha = 0.7,
#            line.size = 0.6,
#            mb.ratio = c(0.6,0.4), text.scale = 1.7, 
#            set_size.numbers_size = 6, set_size.scale_max = 100)
# t1
# dev.off()


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
  k = as.character(Des[i,-c(1,3,4,5,6)])
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

for (i in 1:length(Final$ID)) {
  m = which(Final$Show==Final$Show[i] & Final$Cluster==Final$Cluster[i])
  k = strsplit(Final$geneID[m],"/")
  l = unique(unlist(k))
  Final$FinalCount[i] = length(l)
  Final$Terms[i] = length(m)
  Final$Genes[i] = paste(l, collapse = "/")
  Final$MedC[i] = median(Final$Cramer[m])
  Final$MaxC[i] = max(Final$Cramer[m])
}
tab = data.frame(Treatment = Final$Cluster, ReactomeID = Final$ID, Pathway = Final$Description,
                 ParentNode = Final$Level1, ChildNode =  Final$Level2, GeneRatio = Final$GeneRatio,
                 BgRatio = Final$BgRatio, Adj.P.Val = Final$p.adjust, CramerV = Final$Cramer,
                 CollapsedCramerV = Final$MaxC,CollapsedCount = Final$Count, Genes = Final$genes)
tab = tab[order(tab$Treatment),]
write.xlsx(tab, str_c(here("Res/ThirdSubmission/Clustering/"), "Reactome_All.xlsx"), overwrite = T)
tab = tab[,c(1,2,3,8,9,12)]
tab$Adj.P.Val  = sprintf("%.3f x 10^(%d)", tab$Adj.P.Val/10^floor(log10(abs(tab$Adj.P.Val))), floor(log10(abs(tab$Adj.P.Val))))
tab$CramerV  = sprintf("%.3f x 10^(%d)", tab$CramerV/10^floor(log10(abs(tab$CramerV))), floor(log10(abs(tab$CramerV))))
write.xlsx(tab, str_c(here("Res/ThirdSubmission/Clustering/"), "Reactome_Fig.xlsx"), overwrite = T)

Final = Final[!duplicated(Final[,c(1,13)]),] #Cluster, Group, Show
Final$FinalSize = 7

FinalT = Final
FinalT = droplevels(FinalT)
FinalT$Show[which(FinalT$Show=="The citric acid (TCA) cycle and respiratory electron transport")] = "TCA cycle"
FinalT$Show[which(FinalT$Show=="Signaling by Rho GTPases, Miro GTPases and RHOBTB3")] = "Signaling by Rho GTPases"
FinalT$Show[which(FinalT$Show=="Metabolism of amino acids and derivatives")] = "Metabolism of amino acids"

FinalT$Neg = -log(FinalT$p.adjust)
try = as.character(unique(FinalT$Show))
FinalT$Show = factor(FinalT$Show, levels = unique(FinalT$Show))
FinalT$Show = str_wrap(FinalT$Show, width = 60)
FinalT$Level1 = str_wrap(FinalT$Level1, width = 25)

tiff(str_c(here("Res/ThirdSubmission/Figures/"), "WGCNA_Reactome.tiff"), units="px", width=(6*900), height=(6*700), res=600)

p <- ggplot(FinalT, aes(x = Cluster, y = Show, size = FinalCount, color = MaxC, fill = MaxC)) +
  geom_point(shape = 21, stroke = 1) +
  theme_bw(base_size = 12) +
  scale_color_viridis(alpha = 0.9,
                      breaks = c(min(FinalT$MaxC),(min(FinalT$MaxC)+ max(FinalT$MaxC))/2, max(FinalT$MaxC)),
                      labels = c("0.03", "0.15", "0.30"),
                      discrete = F, name = "Cramer's V\n")+
  scale_fill_viridis(alpha = 0.9,
                     breaks = c(min(FinalT$MaxC),(min(FinalT$MaxC)+ max(FinalT$MaxC))/2, max(FinalT$MaxC)),
                     labels = c("0.03", "0.15", "0.30"),
                     discrete = F, name = "Cramer's V\n")+
  scale_size_continuous(range = c(2,6), name = "Gene Count", breaks = c(min(FinalT$FinalCount),150,max(FinalT$FinalCount)), labels = c("3","150", "300")) + 
  #facet_grid(Level1~.,scales="free",space="free")+
  theme(strip.text.y = element_text(angle = 0, family = "Calibri", face = "bold", size = 12))+
  theme(panel.spacing =unit(0.05, "lines"),
        panel.border = element_rect(color = "#476b6b", fill = NA, size = 1), 
        strip.background = element_rect(color = "#476b6b", size = 1, fill = "#d6dce4"))+
  theme(axis.ticks = element_line(size=1,color='#476b6b')) +
  scale_y_discrete(limits = rev) +
  xlab("Cluster Labels") + ylab("Pathways") + ggtitle("") + theme(panel.background = element_rect(colour = "grey70")) +
  theme(axis.title = element_text(size = 13, face = "bold", family = "Calibri")) +
  theme(axis.text.y = element_text(size=11, face = "bold", family = "Calibri")) +
  theme(axis.text.x = element_text(size=11,face = "bold" ,family = "Calibri")) +
  theme(legend.text=element_text(size=11, family = "Calibri")) +
  theme(legend.title = element_text(size = 12, face = "bold", family = "Calibri")) +
  theme(legend.position="right") 

p
dev.off()
