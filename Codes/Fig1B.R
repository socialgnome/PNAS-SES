#### This file is a script that is responsible for creating 
#### Fig1B in the PNAS manuscript

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
library(ggtext) 
loadfonts()

total_res = readRDS(str_c(here("Res/ThirdSubmission/"), "Plot_Fig1B.rds"))
total_res$FDR_p = -log10(total_res$FDR)
total_res$FDR_p[which(total_res$FDR_p>-log10(0.001))] = -log10(0.001)
total_res = total_res[total_res$Signature!="Aortic Aneurysm",]

## Combined figure

#tiff(str_c(here("Res/ThirdSubmission/Figures/"), "Fig1B_Combined.tiff"), units="px", width=(6*800), height=(6*800), res=600)

p = ggplot(total_res[total_res$colInstance==0 & total_res$Reg=="Downregulated",], aes(x = Treatment,y  = YAxisLabel, color = Class, fill = Class,  size = FDR_p)) +
  geom_point(alpha  = 0.5, shape="\u25D6", stroke = 1,position = position_nudge(x = -0.15)) +
  geom_point(data = total_res[total_res$colInstance==0 & total_res$Reg=="Upregulated",], shape="\u25D7",inherit.aes = T, position = position_nudge(x = +0.15), alpha = 0.5, stroke = 1,show.legend = F) +
  geom_point(data = total_res[total_res$colInstance!=0 & total_res$Reg=="Downregulated" & total_res$Class=="1KI Genes",], shape="\u25D6",inherit.aes = T, position = position_nudge(x = -0.15, y = +0.15), alpha = 0.5, stroke = 1,show.legend = F) +
  geom_point(data = total_res[total_res$colInstance!=0 & total_res$Reg=="Upregulated" & total_res$Class=="1KI Genes",], shape="\u25D7",inherit.aes = T, position = position_nudge(x = +0.15, y = +0.15), alpha = 0.5, stroke = 1,show.legend = F) +
  geom_point(data = total_res[total_res$colInstance!=0 & total_res$Reg=="Downregulated" & total_res$Class=="Without 1KI Genes",], shape="\u25D6",inherit.aes = T, position = position_nudge(x = -0.15, y = -0.15), alpha = 0.5, stroke = 1,show.legend = F) +
  geom_point(data = total_res[total_res$colInstance!=0 & total_res$Reg=="Upregulated" & total_res$Class=="Without 1KI Genes",], shape="\u25D7",inherit.aes = T, position = position_nudge(x = +0.15, y = -0.15), alpha = 0.5, stroke = 1,show.legend = F) +
  
  geom_text(data = total_res[total_res$colInstance==0 & total_res$Reg=="Downregulated",], aes(x = Treatment,y  = YAxisLabel, label = as.character(SelSigClusters)), color = "black", family = "Calibri", fontface="bold",alpha=1, size=3, inherit.aes = FALSE,position = position_nudge(x = -0.15)) +
  geom_text(data = total_res[total_res$colInstance==0 & total_res$Reg=="Upregulated",], aes(x = Treatment,y  = YAxisLabel, label = as.character(SelSigClusters)), color = "black", family = "Calibri", fontface="bold",alpha=1, size=3, inherit.aes = FALSE,position = position_nudge(x = +0.15)) +
  geom_text(data = total_res[total_res$colInstance!=0 & total_res$Reg=="Downregulated" & total_res$Class=="1KI Genes",], aes(x = Treatment,y  = YAxisLabel, label = as.character(SelSigClusters)), color = "black", family = "Calibri", fontface="bold",alpha=1, size=3, inherit.aes = FALSE,position = position_nudge(x = -0.15, y=+0.15)) +
  geom_text(data = total_res[total_res$colInstance!=0 & total_res$Reg=="Upregulated" & total_res$Class=="1KI Genes",], aes(x = Treatment,y  = YAxisLabel, label = as.character(SelSigClusters)), color = "black", family = "Calibri", fontface="bold",alpha=1, size=3, inherit.aes = FALSE,position = position_nudge(x = +0.15, y = +0.15)) +
  geom_text(data = total_res[total_res$colInstance!=0 & total_res$Reg=="Downregulated" & total_res$Class=="Without 1KI Genes",], aes(x = Treatment,y  = YAxisLabel, label = as.character(SelSigClusters)), color = "black", family = "Calibri", fontface="bold",alpha=1, size=3, inherit.aes = FALSE,position = position_nudge(x = -0.15, y = -0.15)) +
  geom_text(data = total_res[total_res$colInstance!=0 & total_res$Reg=="Upregulated" & total_res$Class=="Without 1KI Genes",], aes(x = Treatment,y  = YAxisLabel, label = as.character(SelSigClusters)), color = "black", family = "Calibri", fontface="bold",alpha=1, size=3, inherit.aes = FALSE,position = position_nudge(x = +0.15, y = -0.15)) +
  
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
  theme(axis.text.y = element_markdown(size=12, face = "bold", family = "Calibri")) +
  theme(axis.text.x = element_text(size=12,face = "bold" ,family = "Calibri", angle = 30,hjust = 1)) +
  theme(legend.text=element_text(size=11, family = "Calibri")) +
  theme(legend.text.align = 0) +
  theme(legend.title = element_text(size = 12, face = "bold", family = "Calibri")) +
  theme(legend.position="right") +
  scale_x_discrete(drop = FALSE) +
  scale_y_discrete(limits = rev)  +
  guides(color = guide_legend(override.aes = list(size = 7))) +
  guides(fill = guide_legend(override.aes = list(size = 7))) +
  guides(shape = guide_legend(override.aes = list(size = 7))) 


#p

#dev.off()


## 1KI Figure

res = total_res[total_res$Class=="1KI Genes",]
tiff(str_c(here("Res/ThirdSubmission/Figures/"), "Fig1B.tiff"), units="px", width=(6*625), height=(6*720), res=600)

p = ggplot(res[res$Reg=="Downregulated",], aes(x = Treatment,y  = YAxisLabel_A, color = Reg, fill = Reg,  size = FDR_p)) +
  geom_point(alpha  = 0.90, shape="\u25D6", stroke = 1,position = position_nudge(x = -0.18)) +
  geom_point(data = res[res$Reg=="Upregulated",], shape="\u25D7",inherit.aes = T, position = position_nudge(x = +0.18), alpha = 0.90, stroke = 1,show.legend = F) +
  
  geom_text(data = res[res$Reg=="Downregulated",], aes(x = Treatment,y  = YAxisLabel_A, label = as.character(SelSigClusters)), color = "black", family = "Calibri", fontface="bold",alpha=1, size=3.25, inherit.aes = FALSE,position = position_nudge(x = -0.18)) +
  geom_text(data = res[res$Reg=="Upregulated",], aes(x = Treatment,y  = YAxisLabel_A, label = as.character(SelSigClusters)), color = "black", family = "Calibri", fontface="bold",alpha=1, size=3.25, inherit.aes = FALSE,position = position_nudge(x = +0.18)) +
  
  theme_bw(base_size = 12) +
  scale_color_manual(values = c("#4379a7","#c93311"), name = "Direction\nof change",
                     #guide = guide_legend(override.aes = list(shape = c("\u25D0","\u25D1"), size =7))) +
                    guide = guide_legend(override.aes = list(shape = c("\u25D6","\u25D7"), size =7))) +
  scale_fill_manual(values = c("#4379a7","#c93311"), name = "Direction\nof change",
                    #guide = guide_legend(override.aes = list(shape = c("\u25D0","\u25D1"), size =7))) +
                    guide = guide_legend(override.aes = list(shape = c("\u25D6","\u25D7"), size =7))) +
  scale_size_continuous(range = c(6,12),
                        name = "Adjusted\n p-value", 
                        limits = c(-log10(0.05), -log10(0.001)),
                        breaks = c(-log10(0.05),-log10(0.01),-log10(0.001)), 
                        labels = c(expression(italic("p")~'<'~0.05), expression(italic("p")~'<'~0.01),expression(italic("p")~'<'~0.001)),
                        guide = guide_legend(override.aes = list(shape = c(21), fill = "NA", colour = "grey30"))) + 
  xlab("SES Indicators") + ylab("mRNA Signatures") + ggtitle("") + theme(panel.background = element_rect(colour = "grey70", size = 1.5)) +
  theme(axis.title = element_text(size = 13, face = "bold", family = "Calibri")) +
  theme(axis.text.y = element_markdown(size=12, face = "bold", family = "Calibri")) +
  theme(axis.text.x = element_text(size=12,face = "bold" ,family = "Calibri", angle = 30,hjust = 1)) +
  theme(legend.text=element_text(size=11, family = "Calibri")) +
  theme(legend.text.align = 0) +
  theme(legend.title = element_text(size = 12, face = "bold", family = "Calibri")) +
  theme(legend.position=c(1.325,0.55)) +
  theme(plot.margin = unit(c(5, 120 , 5, 5), "pt")) +
  scale_x_discrete(drop = FALSE)  
p

dev.off()


## Without 1KI Figure

res = total_res[total_res$Class=="Without 1KI Genes",]
tiff(str_c(here("Res/ThirdSubmission/Figures/"), "Fig1B_S.tiff"), units="px", width=(6*625), height=(6*700), res=600)

p = ggplot(res[res$Reg=="Downregulated",], aes(x = Treatment,y  = YAxisLabel_A, color = Reg, fill = Reg,  size = FDR_p)) +
  geom_point(alpha  = 0.90, shape="\u25D6", stroke = 1,position = position_nudge(x = -0.18)) +
  geom_point(data = res[res$Reg=="Upregulated",], shape="\u25D7",inherit.aes = T, position = position_nudge(x = +0.18), alpha = 0.9, stroke = 1,show.legend = F) +
  
  geom_text(data = res[res$Reg=="Downregulated",], aes(x = Treatment,y  = YAxisLabel_A, label = as.character(SelSigClusters)), color = "black", family = "Calibri", fontface="bold",alpha=1, size=3.25, inherit.aes = FALSE,position = position_nudge(x = -0.18)) +
  geom_text(data = res[res$Reg=="Upregulated",], aes(x = Treatment,y  = YAxisLabel_A, label = as.character(SelSigClusters)), color = "black", family = "Calibri", fontface="bold",alpha=1, size=3.25, inherit.aes = FALSE,position = position_nudge(x = +0.18)) +
  
  theme_bw(base_size = 12) +
  scale_color_manual(values = c("#4379a7","#c93311"), name = "Direction\nof change",
                     #guide = guide_legend(override.aes = list(shape = c("\u25D0","\u25D1"), size =7))) +
                     guide = guide_legend(override.aes = list(shape = c("\u25D6","\u25D7"), size =7))) +
  scale_fill_manual(values = c("#4379a7","#c93311"), name = "Direction\nof change",
                    #guide = guide_legend(override.aes = list(shape = c("\u25D0","\u25D1"), size =7))) +
                    guide = guide_legend(override.aes = list(shape = c("\u25D6","\u25D7"), size =7))) +
  scale_size_continuous(range = c(6,12),
                        name = "Adjusted\n p-value", 
                        limits = c(-log10(0.05), -log10(0.001)),
                        breaks = c(-log10(0.05),-log10(0.01),-log10(0.001)), 
                        labels = c(expression(italic("p")~'<'~0.05), expression(italic("p")~'<'~0.01),expression(italic("p")~'<'~0.001)),
                        guide = guide_legend(override.aes = list(shape = c(21), fill = "NA", colour = "grey30"))) + 
  xlab("SES Indicators") + ylab("mRNA Signatures") + ggtitle("") + theme(panel.background = element_rect(colour = "grey70", size = 1.5)) +
  theme(axis.title = element_text(size = 13, face = "bold", family = "Calibri")) +
  theme(axis.text.y = element_markdown(size=12, face = "bold", family = "Calibri")) +
  theme(axis.text.x = element_text(size=12,face = "bold" ,family = "Calibri", angle = 30,hjust = 1)) +
  theme(legend.text=element_text(size=11, family = "Calibri")) +
  theme(legend.text.align = 0) +
  theme(legend.title = element_text(size = 12, face = "bold", family = "Calibri")) +
  theme(legend.position=c(1.325,0.55)) +
  theme(plot.margin = unit(c(5, 120 , 5, 5), "pt")) +
  scale_x_discrete(drop = FALSE)  
p

dev.off()