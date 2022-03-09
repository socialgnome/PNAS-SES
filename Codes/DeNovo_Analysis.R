#### This file is a script that is responsible for creating 
#### Fig2 in the PNAS manuscript - DeNovo analysis

#### Prerequisites: 1. DE.R

#### Follow up files: None

## Setup env
library(here)
dev.off()
source(here("R/ThirdSubmission/Final/startup.R"))
library(openxlsx)
signatures = readRDS("~/Projects/PNAS/data/Signatures.rds")

## DeNovo analysis
disgenes = unique(unlist(signatures))
treatment = c("ses_sss_composite","sss_5","SEI_ff5",
              "edu_max","income_hh_ff5")

all = list()
up = list()
down = list()
for (i in 1:length(treatment)) {
  df = readRDS(str_c(here("Res/ThirdSubmission/DE/"), treatment[[i]],".rds"))
  df = df[-which(df$gene %in% disgenes),]
  df$FDR = p.adjust(df$P.Value)
  all[[i]] = df$gene[df$FDR<0.05]
  up[[i]] = df$gene[df$FDR<0.05 & df$logFC>0]
  down[[i]] = df$gene[df$FDR<0.05 & df$logFC<0]
}
treatment = c("SES Composite","Subjective Social Status","Occupation","Education","Income")

all = setNames(all, treatment)
down = setNames(down, treatment)
up = setNames(up, treatment)

## UpsetR Plots
tiff(str_c(here("Res/ThirdSubmission/Figures/"), "DeNovo_all.tiff"), units="px", width=(3*1150), height=(3*575), res=300)
t1 = upset(fromList(all),
           mainbar.y.label = "Intersecting Genes",
           sets.x.label = "Total DE genes",
           keep.order = T,
           main.bar.color = "#555555",
           matrix.color  = "#326a97",
           set_size.show = T,
           point.size = 2.8,
           shade.alpha = 0.2,
           matrix.dot.alpha = 0.7,
           line.size = 0.6,
           mb.ratio = c(0.65,0.35), text.scale = 1.7, 
           set_size.numbers_size = 6, set_size.scale_max = 200)
t1
dev.off()

tiff(str_c(here("Res/ThirdSubmission/Figures/"), "DeNovo_up.tiff"), units="px", width=(3*1150), height=(3*575), res=300)
t1 = upset(fromList(up),
           mainbar.y.label = "Intersecting Genes",
           sets.x.label = "Total Upregulated genes",
           keep.order = T,
           main.bar.color = "#555555",
           matrix.color  = brewer.pal(8, "RdYlGn")[1],
           set_size.show = T,
           point.size = 2.8,
           shade.alpha = 0.2,
           matrix.dot.alpha = 0.7,
           line.size = 0.6,
           mb.ratio = c(0.65,0.35), text.scale = 1.7, 
           set_size.numbers_size = 6, set_size.scale_max = 80)
t1
dev.off()

tiff(str_c(here("Res/ThirdSubmission/Figures/"), "DeNovo_down.tiff"), units="px", width=(3*1150), height=(3*575), res=300)
t1 = upset(fromList(down),
           mainbar.y.label = "Intersecting Genes",
           sets.x.label = "Total Downregulated genes",
           keep.order = T,
           main.bar.color = "#555555",
           matrix.color  = brewer.pal(8, "RdYlGn")[8],
           set_size.show = T,
           point.size = 2.8,
           shade.alpha = 0.2,
           matrix.dot.alpha = 0.7,
           line.size = 0.6,
           mb.ratio = c(0.65,0.35), text.scale = 1.7, 
           set_size.numbers_size = 6, set_size.scale_max = 80)
t1
dev.off()

## Biomart (Downloaded 18.01.22)
sapiens_ensembl <- readRDS(str_c(here("data/"), "sapiens_ensembl.rds"))

all_ent = list()
conv = list()
for (i in 1:length(all)) {
  df = biomaRt::getBM(attributes = c("hgnc_symbol", "entrezgene_id"), filters = "hgnc_symbol", values =all[[i]], mart = sapiens_ensembl, useCache = F)
  conv[[i]] = df
  df = as.character(df$entrezgene_id[!is.na(df$entrezgene_id)])
  all_ent[[i]] = df
}

all_ent = setNames(all_ent, treatment)
eres = compareCluster(all_ent, fun = "enrichPathway", organism = "human", pvalueCutoff = 0.05, pAdjustMethod = "BH")
eres = eres@compareClusterResult
eres = eres[eres$Count>2,]

eres$genes = 0
for (i in 1:length(eres$Cluster)) {
  eres$genes[[i]] = paste(conv[[which(treatment==eres$Cluster[[i]])]]$hgnc_symbol[match(unlist(strsplit(eres$geneID[[i]], "/")),conv[[which(treatment==eres$Cluster[[i]])]]$entrezgene_id)],collapse = ",")
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

tiff(str_c(here("Res/ThirdSubmission/Figures/"), "DeNovo_Pathway.tiff"), units="px", width=(3*1150), height=(3*675), res=300)
t1 = upset(fromList(lists[1:5]),
           mainbar.y.label = "Intersecting enriched pathways",
           sets.x.label = "Total enriched pathways",
           keep.order = T,
           main.bar.color = "#555555",
           matrix.color  = "#326a97",
           set_size.show = T,
           point.size = 2.8,
           shade.alpha = 0.2,
           matrix.dot.alpha = 0.7,
           line.size = 0.6,
           mb.ratio = c(0.6,0.4), text.scale = 1.7, 
           set_size.numbers_size = 6, set_size.scale_max = 75)
t1
dev.off()


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
Final$FinalUp = 0
Final$FinalDown = 0

for (i in 1:length(Final$ID)) {
  m = which(Final$Show==Final$Show[i] & Final$Cluster==Final$Cluster[i])
  k = strsplit(Final$geneID[m],"/")
  l = unique(unlist(k))
  Final$FinalCount[i] = length(l)
  Final$Terms[i] = length(m)
  Final$Genes[i] = paste(l, collapse = "/")
  Final$MedC[i] = median(Final$Cramer[m])
  Final$MaxC[i] = max(Final$Cramer[m])
  Final$FinalUp[i] = length(intersect(conv[[which(treatment==Final$Cluster[[i]])]]$hgnc_symbol[match(unlist(strsplit(Final$geneID[m], "/")),conv[[which(treatment==Final$Cluster[[i]])]]$entrezgene_id)],up[[which(treatment==Final$Cluster[[i]])]]))
  Final$FinalDown[i] = length(intersect(conv[[which(treatment==Final$Cluster[[i]])]]$hgnc_symbol[match(unlist(strsplit(Final$geneID[m], "/")),conv[[which(treatment==Final$Cluster[[i]])]]$entrezgene_id)],down[[which(treatment==Final$Cluster[[i]])]]))
}
tab = data.frame(Treatment = Final$Cluster, ReactomeID = Final$ID, Pathway = Final$Description,
                 ParentNode = Final$Level1, ChildNode =  Final$Level2, GeneRatio = Final$GeneRatio,
                 BgRatio = Final$BgRatio, Adj.P.Val = Final$p.adjust, CramerV = Final$Cramer,
                 CollapsedCramerV = Final$MaxC,CollapsedCount = Final$FinalCount,UpregulatedCount=Final$FinalUp, DownregulatedCount = Final$FinalDown, Genes = Final$genes)
tab = tab[order(tab$Treatment),]
write.xlsx(tab, str_c(here("Res/ThirdSubmission/DeNovo/"), "Reactome_All.xlsx"), overwrite = T)
tab = tab[,c(1,2,3,8,9,12)]
tab$Adj.P.Val  = sprintf("%.3f x 10^(%d)", tab$Adj.P.Val/10^floor(log10(abs(tab$Adj.P.Val))), floor(log10(abs(tab$Adj.P.Val))))
tab$CramerV  = sprintf("%.3f x 10^(%d)", tab$CramerV/10^floor(log10(abs(tab$CramerV))), floor(log10(abs(tab$CramerV))))
write.xlsx(tab, str_c(here("Res/ThirdSubmission/DeNovo/"), "Reactome_Fig.xlsx"), overwrite = T)

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

FinalT$Cluster  = factor(FinalT$Cluster, levels =  c("SES Composite", "Income", "Subjective Social Status"))

tiff(str_c(here("Res/ThirdSubmission/Figures/"), "Fig2.tiff"), units="px", width=(6*780), height=(6*725), res=600)

p <- ggplot(FinalT, aes(x = Cluster, y = Show, size = FinalCount, color = MaxC, fill = MaxC)) +
  geom_point(shape = 21, stroke = 1) +
  theme_bw(base_size = 12) +
  scale_color_viridis(alpha = 0.9,
                      breaks = c(min(FinalT$MaxC),(min(FinalT$MaxC)+ max(FinalT$MaxC))/2, max(FinalT$MaxC)),
                      labels = c("0.04", "0.20", "0.40"),
                      discrete = F, name = "Cramer's V\n")+
  scale_fill_viridis(alpha = 0.9,
                     breaks = c(min(FinalT$MaxC),(min(FinalT$MaxC)+ max(FinalT$MaxC))/2, max(FinalT$MaxC)),
                     labels = c("0.04", "0.20", "0.40"),
                     discrete = F, name = "Cramer's V\n")+
  scale_size_continuous(range = c(4,12), name = "Gene Count", breaks = c(min(FinalT$FinalCount),20,max(FinalT$FinalCount)), labels = c("3","20", "40")) + 
  facet_grid(Level1~.,scales="free",space="free")+
  theme(strip.text.y = element_text(angle = 0, family = "Calibri", face = "bold", size = 12))+
  theme(panel.spacing =unit(0.1, "lines"),
        panel.border = element_rect(color = "#476b6b", fill = NA, size = 1), 
        strip.background = element_rect(color = "#476b6b", size = 1, fill = "#d6dce4"))+
  theme(axis.ticks = element_line(size=1,color='#476b6b')) +
  xlab("SES Indicators") + ylab("Pathways") + ggtitle("") + theme(panel.background = element_rect(colour = "grey70")) +
  theme(axis.title = element_text(size = 13, face = "bold", family = "Calibri")) +
  theme(axis.text.y = element_text(size=12, face = "bold", family = "Calibri")) +
  theme(axis.text.x = element_text(size=12,face = "bold" ,family = "Calibri", angle = 30,hjust = 1)) +
  theme(legend.text=element_text(size=11, family = "Calibri")) +
  theme(legend.title = element_text(size = 12, face = "bold", family = "Calibri")) +
  theme(legend.position="right") 

p
dev.off()


### Regulation
FinalT$Reg = 1 -  FinalT$FinalDown/(FinalT$FinalCount)

tiff(str_c(here("Res/ThirdSubmission/Figures/"), "Fig2_v2.tiff"), units="px", width=(6*780), height=(6*725), res=600)

p <- ggplot(FinalT, aes(x = Cluster, y = Show, size = FinalCount, color = Reg, fill = Reg)) +
  geom_point(shape = 21, stroke = 1) +
  theme_bw(base_size = 12) +
  scale_color_gradient(low = "#1A9850", high = "#D73027", 
                       limits = c(0,1), 
                       breaks = c(0,1), 
                       labels = c("Downregulated", "Upregulated"), 
                       name = "Direction\nof change") +
  scale_fill_gradient(low = "#1A9850", high = "#D73027", 
                       limits = c(0,1), 
                       breaks = c(0,1), 
                       labels = c("Downregulated", "Upregulated"), 
                       name = "Direction\nof change") +
  scale_size_continuous(range = c(4,12), name = "Gene Count", breaks = c(min(FinalT$FinalCount),20,max(FinalT$FinalCount)), labels = c("3","20", "40")) + 
  facet_grid(Level1~.,scales="free",space="free")+
  theme(strip.text.y = element_text(angle = 0, family = "Calibri", face = "bold", size = 12))+
  theme(panel.spacing =unit(0.1, "lines"),
        panel.border = element_rect(color = "#476b6b", fill = NA, size = 1), 
        strip.background = element_rect(color = "#476b6b", size = 1, fill = "#d6dce4"))+
  theme(axis.ticks = element_line(size=1,color='#476b6b')) +
  xlab("SES Indicators") + ylab("Pathways") + ggtitle("") + theme(panel.background = element_rect(colour = "grey70")) +
  theme(axis.title = element_text(size = 13, face = "bold", family = "Calibri")) +
  theme(axis.text.y = element_text(size=12, face = "bold", family = "Calibri")) +
  theme(axis.text.x = element_text(size=12,face = "bold" ,family = "Calibri", angle = 30,hjust = 1)) +
  theme(legend.text=element_text(size=11, family = "Calibri")) +
  theme(legend.title = element_text(size = 12, face = "bold", family = "Calibri")) +
  theme(legend.position="right") 

p
dev.off()