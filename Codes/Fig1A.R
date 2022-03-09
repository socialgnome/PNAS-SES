#### This file is a script that is responsible for creating 
#### Fig1A in the PNAS manuscript

#### Prerequisites: 1. OmnibusTest_Fig1A.R

#### Follow up files: None

## Setup env
library(here)
dev.off()
rm(list=ls(all=TRUE))
library(ggplot2)
library(extrafont)
library(here)
library(stringr)
loadfonts()
res = readRDS(str_c(here("Res/ThirdSubmission/"), "Omnibus_Fig1A.rds"))

res$pval2 = -log10(res$P)
res$pval2[which(res$pval2>-log10(0.001))] = -log10(0.001)
res = res[order(res$Class),]

res$Instance = 0
for (i in 1:length(res$Treatment)) {
  res$Instance[i] = length(which(res$Treatment==res$Treatment[i] & res$Signature==res$Signature[i]))
}

## Combined figure

#tiff(str_c(here("Res/ThirdSubmission/Figures/"), "Fig1A_Combined.tiff"), units="px", width=(6*600), height=(6*800), res=600)
p = ggplot(res[res$Instance==1,], aes(x = Treatment,y  = Signature, color = Class, fill = Class,  size = pval2)) +
  geom_point(alpha  = 0.5,shape = 21, stroke = 1) +
  
  geom_point(data = res[res$Instance==2 & !duplicated(res[,c(1,2)]),], inherit.aes = T, position = position_nudge(y = +0.15), alpha = 0.5, stroke = 1,show.legend = F) + 
  geom_point(data = res[res$Instance==2 & duplicated(res[,c(1,2)]),], inherit.aes = T, position = position_nudge(y = -0.15), alpha = 0.5, stroke = 1,show.legend = F) +
  
  theme_bw(base_size = 12) +
  scale_color_manual(values = c("#01665e","#ff7f0e"), name = "1KI Genes") +
  scale_fill_manual(values = c("#01665e","#ff7f0e"), name = "1KI Genes") +
  scale_size_continuous(range = c(6,12),
                        name = "Adjusted\n p-value", 
                        limits = c(-log10(0.05), -log10(0.001)),
                        breaks = c(-log10(0.05),-log10(0.01),-log10(0.001)), 
                        labels = c(expression(italic("p")~'<'~0.05), expression(italic("p")~'<'~0.01),expression(italic("p")~'<'~0.001))) + 
  xlab("SES Indicators") + ylab("mRNA Signatures") + ggtitle("") + theme(panel.background = element_rect(colour = "grey70", size = 1.5)) +
  theme(axis.title = element_text(size = 13, face = "bold", family = "Calibri")) +
  theme(axis.text.y = element_text(size=12, face = "bold", family = "Calibri")) +
  theme(axis.text.x = element_text(size=12,face = "bold" ,family = "Calibri", angle = 30,hjust = 1)) +
  theme(legend.text=element_text(size=11, family = "Calibri")) +
  theme(legend.text.align = 0) +
  theme(legend.title = element_text(size = 12, face = "bold", family = "Calibri")) +
  theme(legend.position="right") +
  scale_y_discrete(limits = rev)  +
  guides(color = guide_legend(override.aes = list(size = 7))) +
  guides(fill = guide_legend(override.aes = list(size = 7))) +
  guides(shape = guide_legend(override.aes = list(size = 7))) 

#dev.off()
allres = res

## 1KI Figure
res = allres[allres$Class == "1KI Genes",]
tiff(str_c(here("Res/ThirdSubmission/Figures/"), "Fig1A.tiff"), units="px", width=(6*600), height=(6*720), res=600)
p = ggplot(res, aes(x = Treatment,y  = Signature, color = Class, fill = Class,  size = pval2)) +
  geom_point(alpha  = 0.5,shape = 21, stroke = 1) +
  
  theme_bw(base_size = 12) +
  scale_color_manual(values = c("#01665e"), name = "1KI Genes", guide = "none") +
  scale_fill_manual(values = c("#01665e"), name = "1KI Genes", guide = "none") +
  scale_size_continuous(range = c(6,12),
                        name = "Adjusted\n p-value", 
                        limits = c(-log10(0.05), -log10(0.001)),
                        breaks = c(-log10(0.05),-log10(0.01),-log10(0.001)), 
                        labels = c(expression(italic("p")~'<'~0.05), expression(italic("p")~'<'~0.01),expression(italic("p")~'<'~0.001)),
                        guide = guide_legend(override.aes = list(shape = c(21), fill = "NA", colour = "grey30"))) + 
  xlab("SES Indicators") + ylab("mRNA Signatures") + ggtitle("") + theme(panel.background = element_rect(colour = "grey70", size = 1.5)) +
  theme(axis.title = element_text(size = 13, face = "bold", family = "Calibri")) +
  theme(axis.text.y = element_text(size=12, face = "bold", family = "Calibri")) +
  theme(axis.text.x = element_text(size=12,face = "bold" ,family = "Calibri", angle = 30,hjust = 1)) +
  theme(legend.text=element_text(size=11, family = "Calibri")) +
  theme(legend.text.align = 0) +
  theme(legend.title = element_text(size = 12, face = "bold", family = "Calibri")) +
  theme(legend.position="right") 
p
dev.off()

## Without 1KI Figure
res = allres[allres$Class == "Without 1KI Genes",]
tiff(str_c(here("Res/ThirdSubmission/Figures/"), "Fig1A_S.tiff"), units="px", width=(6*600), height=(6*700), res=600)
p = ggplot(res, aes(x = Treatment,y  = Signature, color = Class, fill = Class,  size = pval2)) +
  geom_point(alpha  = 0.5,shape = 21, stroke = 1) +
  
  theme_bw(base_size = 12) +
  scale_color_manual(values = c("#ff7f0e"), name = "1KI Genes", guide = "none") +
  scale_fill_manual(values = c("#ff7f0e"), name = "1KI Genes", guide = "none") +
  scale_size_continuous(range = c(6,12),
                        name = "Adjusted\n p-value", 
                        limits = c(-log10(0.05), -log10(0.001)),
                        breaks = c(-log10(0.05),-log10(0.01),-log10(0.001)), 
                        labels = c(expression(italic("p")~'<'~0.05), expression(italic("p")~'<'~0.01),expression(italic("p")~'<'~0.001)),
                        guide = guide_legend(override.aes = list(shape = c(21), fill = "NA", colour = "grey30"))) + 
  xlab("SES Indicators") + ylab("mRNA Signatures") + ggtitle("") + theme(panel.background = element_rect(colour = "grey70", size = 1.5)) +
  theme(axis.title = element_text(size = 13, face = "bold", family = "Calibri")) +
  theme(axis.text.y = element_text(size=12, face = "bold", family = "Calibri")) +
  theme(axis.text.x = element_text(size=12,face = "bold" ,family = "Calibri", angle = 30,hjust = 1)) +
  theme(legend.text=element_text(size=11, family = "Calibri")) +
  theme(legend.text.align = 0) +
  theme(legend.title = element_text(size = 12, face = "bold", family = "Calibri")) +
  theme(legend.position="right") +
  scale_y_discrete(drop = T) 
p
dev.off()