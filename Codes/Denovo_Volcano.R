#### This file is a script that is responsible for creating 
#### volcano plots in the PNAS manuscript - DeNovo analysis

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

for (i in 1:length(treatment)) {
  df = readRDS(str_c(here("Res/ThirdSubmission/DE/"), treatment[[i]],".rds"))
  df = df[-which(df$gene %in% disgenes),]
  df$FDR = p.adjust(df$P.Value)
  df$Neg = -log10(df$FDR)
  df$class = "Non significant"
  df$class[df$FDR<0.05 & df$logFC>0] = "Upregulated"
  df$class[df$FDR<0.05 & df$logFC<0] = "Downregulated"
  
  all[[i]] = df
}
treatment = c("SES Composite","Subjective Social Status","Occupation","Education","Income")
all = setNames(all, treatment)

####Volcano plots
data = all[[1]]
tiff(str_c(here("Res/ThirdSubmission/Figures/"), "Volcano_SES.tiff"), units="px", width=(6*750), height=(6*400), res=600)
p = ggplot(data, aes(x = logFC, y = Neg, color = class, fill = class)) +
  geom_point(shape = 21, size = 3, alpha = 0.7) +
  theme_bw(base_size = 12) +
  scale_color_manual(values = c("#4379a7","#767676","#c93311"), name = "Direction\nof change") +
  scale_fill_manual(values = c("#4379a7","#767676","#c93311"), name = "Direction\nof change") +
  xlab(expression(Log[2]~fold~change)) + ylab(expression(-Log[10]~italic(p))) + 
  ggtitle(treatment[[1]]) +
  theme(plot.title = element_text(size = 14, face = "bold", family = "calibri",hjust = 0.5)) +
  theme(panel.background = element_rect(colour = "grey70", size = 1.5)) +
  xlim(c(-0.08, 0.08)) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  theme(axis.title = element_text(size = 13, face = "bold", family = "Calibri")) +
  theme(axis.text.y = element_text(size=12, face = "bold", family = "Calibri")) +
  theme(axis.text.x = element_text(size=12,face = "bold" ,family = "Calibri")) +
  theme(legend.text=element_text(size=11, family = "Calibri")) +
  theme(legend.text.align = 0) + 
  theme(legend.title = element_text(size = 12, face = "bold", family = "Calibri")) 
p 
dev.off()
  

data = all[[2]]
tiff(str_c(here("Res/ThirdSubmission/Figures/"), "Volcano_SSS.tiff"), units="px", width=(6*750), height=(6*400), res=600)
p = ggplot(data, aes(x = logFC, y = Neg, color = class, fill = class)) +
  geom_point(shape = 21, size = 3, alpha = 0.7) +
  theme_bw(base_size = 12) +
  scale_color_manual(values = c("#4379a7","#767676","#c93311"), name = "Direction\nof change") +
  scale_fill_manual(values = c("#4379a7","#767676","#c93311"), name = "Direction\nof change") +
  xlab(expression(Log[2]~fold~change)) + ylab(expression(-Log[10]~italic(p))) + 
  ggtitle(treatment[[2]]) +
  theme(plot.title = element_text(size = 14, face = "bold", family = "calibri",hjust = 0.5)) +
  theme(panel.background = element_rect(colour = "grey70", size = 1.5)) +
  xlim(c(-0.1, 0.1)) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  theme(axis.title = element_text(size = 13, face = "bold", family = "Calibri")) +
  theme(axis.text.y = element_text(size=12, face = "bold", family = "Calibri")) +
  theme(axis.text.x = element_text(size=12,face = "bold" ,family = "Calibri")) +
  theme(legend.text=element_text(size=11, family = "Calibri")) +
  theme(legend.text.align = 0) + 
  theme(legend.title = element_text(size = 12, face = "bold", family = "Calibri")) 
p 
dev.off()


data = all[[3]]
tiff(str_c(here("Res/ThirdSubmission/Figures/"), "Volcano_Occupation.tiff"), units="px", width=(6*750), height=(6*400), res=600)
p = ggplot(data, aes(x = logFC, y = Neg, color = class, fill = class)) +
  geom_point(shape = 21, size = 3, alpha = 0.7) +
  theme_bw(base_size = 12) +
  scale_color_manual(values = c("#4379a7","#767676","#c93311"), name = "Direction\nof change") +
  scale_fill_manual(values = c("#4379a7","#767676","#c93311"), name = "Direction\nof change") +
  xlab(expression(Log[2]~fold~change)) + ylab(expression(-Log[10]~italic(p))) + 
  ggtitle(treatment[[3]]) +
  theme(plot.title = element_text(size = 14, face = "bold", family = "calibri",hjust = 0.5)) +
  theme(panel.background = element_rect(colour = "grey70", size = 1.5)) +
  #xlim(c(-0.1, 0.1)) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  theme(axis.title = element_text(size = 13, face = "bold", family = "Calibri")) +
  theme(axis.text.y = element_text(size=12, face = "bold", family = "Calibri")) +
  theme(axis.text.x = element_text(size=12,face = "bold" ,family = "Calibri")) +
  theme(legend.text=element_text(size=11, family = "Calibri")) +
  theme(legend.text.align = 0) + 
  theme(legend.title = element_text(size = 12, face = "bold", family = "Calibri")) 
p 
dev.off()

data = all[[4]]
tiff(str_c(here("Res/ThirdSubmission/Figures/"), "Volcano_Education.tiff"), units="px", width=(6*750), height=(6*400), res=600)
p = ggplot(data, aes(x = logFC, y = Neg, color = class, fill = class)) +
  geom_point(shape = 21, size = 3, alpha = 0.7) +
  theme_bw(base_size = 12) +
  scale_color_manual(values = c("#4379a7","#767676","#c93311"), name = "Direction\nof change") +
  scale_fill_manual(values = c("#4379a7","#767676","#c93311"), name = "Direction\nof change") +
  xlab(expression(Log[2]~fold~change)) + ylab(expression(-Log[10]~italic(p))) + 
  ggtitle(treatment[[4]]) +
  theme(plot.title = element_text(size = 14, face = "bold", family = "calibri",hjust = 0.5)) +
  theme(panel.background = element_rect(colour = "grey70", size = 1.5)) +
  xlim(c(-0.6, 0.6)) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  theme(axis.title = element_text(size = 13, face = "bold", family = "Calibri")) +
  theme(axis.text.y = element_text(size=12, face = "bold", family = "Calibri")) +
  theme(axis.text.x = element_text(size=12,face = "bold" ,family = "Calibri")) +
  theme(legend.text=element_text(size=11, family = "Calibri")) +
  theme(legend.text.align = 0) + 
  theme(legend.title = element_text(size = 12, face = "bold", family = "Calibri")) 
p 
dev.off()

data = all[[5]]
tiff(str_c(here("Res/ThirdSubmission/Figures/"), "Volcano_Income.tiff"), units="px", width=(6*750), height=(6*400), res=600)
p = ggplot(data, aes(x = logFC, y = Neg, color = class, fill = class)) +
  geom_point(shape = 21, size = 3, alpha = 0.7) +
  theme_bw(base_size = 12) +
  scale_color_manual(values = c("#4379a7","#767676","#c93311"), name = "Direction\nof change") +
  scale_fill_manual(values = c("#4379a7","#767676","#c93311"), name = "Direction\nof change") +
  xlab(expression(Log[2]~fold~change)) + ylab(expression(-Log[10]~italic(p))) +
  ggtitle(treatment[[5]]) +
  theme(plot.title = element_text(size = 14, face = "bold", family = "calibri",hjust = 0.5)) +
  theme(panel.background = element_rect(colour = "grey70", size = 1.5)) +
  xlim(c(-0.2, 0.2)) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  theme(axis.title = element_text(size = 13, face = "bold", family = "Calibri")) +
  theme(axis.text.y = element_text(size=12, face = "bold", family = "Calibri")) +
  theme(axis.text.x = element_text(size=12,face = "bold" ,family = "Calibri")) +
  theme(legend.text=element_text(size=11, family = "Calibri")) +
  theme(legend.text.align = 0) + 
  theme(legend.title = element_text(size = 12, face = "bold", family = "Calibri")) 
p 
dev.off()