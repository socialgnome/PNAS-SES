#### This file is a script that is responsible for the preparation for plotting
#### DeNovo clustering 

#### Prerequisites: 1. DeNovo_Clustering.R

#### Follow up files: 1. Denovo_Clustering_Fig.R 

## Setup env

library(here)
dev.off()
source(here("R/ThirdSubmission/Final/startup.R"))

## Read Clustering results
res = readRDS(str_c(here("Res/ThirdSubmission/DeNovo/"), "Clustering_Tab.rds"))
fig = rbind(data.frame(Signature = res$Signature, Treatment = res$Treatment, 
                       TotalGenes = res$Total_Genes, TotalClusters = res$Total_Clusters,
                       SelGenes = res$Up_Genes, SelClusters = res$Up_Clusters, SelSigGenes = res$Up_Sig_Genes,
                       SelSigClusters = res$Up_Sig_Clusters, FDR = res$Up_Min_FDR, Reg = "Upregulated"),
            data.frame(Signature = res$Signature, Treatment = res$Treatment, 
                       TotalGenes = res$Total_Genes, TotalClusters = res$Total_Clusters,
                       SelGenes = res$Down_Genes, SelClusters = res$Down_Clusters, SelSigGenes = res$Down_Sig_Genes,
                       SelSigClusters = res$Down_Sig_Clusters, FDR = res$Down_Min_FDR, Reg = "Downregulated")
)

fig = fig[fig$FDR<0.05,]

## Mutate data frame for Plotting
fig = fig %>%
  mutate(Treatment  = case_when(Treatment == "ses_sss_composite" ~ "SES Composite",
                                Treatment == "sss_5" ~ "Subjective Social Status", 
                                Treatment == "SEI_ff5" ~ "Occupation",
                                Treatment == "edu_max" ~ "Education",
                                Treatment == "income_hh_ff5" ~ "Income"),
         Signature = case_when(Signature == "ses_sss_composite" ~ "SES Composite",
                               Signature == "sss_5" ~ "Subjective Social Status", 
                               Signature == "SEI_ff5" ~ "Occupation",
                               Signature == "edu_max" ~ "Education",
                               Signature == "income_hh_ff5" ~ "Income")
  )

fig$YAxisLabel_A = 0
for (i in 1:length(fig$Signature)) {
  fig$YAxisLabel_A[i] = paste("(",head(fig$TotalClusters[which(fig$Signature==fig$Signature[i])],n=1), ")", sep = "")
}
fig$YAxisLabel_A = paste(fig$Signature,"<br/>",fig$YAxisLabel_A, sep = "")

fig$Treatment = factor(fig$Treatment, levels = c("SES Composite","Education","Income","Occupation","Subjective Social Status"))

fig$YAxisLabel_A = factor(fig$YAxisLabel_A, levels = c("SES Composite<br/>(21)",
                                                       "Education<br/>(5)",
                                                       "Income<br/>(14)",
                                                       "Occupation<br/>(4)",
                                                       "Subjective Social Status<br/>(15)"))

library(ggplot2)
library(extrafont)
library(here)
library(stringr)
library(ggtext)
loadfonts()

total_res = fig
total_res$FDR_p = -log10(total_res$FDR)
total_res$FDR_p[which(total_res$FDR_p>-log10(0.001))] = -log10(0.001)

res = total_res
tiff(str_c(here("Res/ThirdSubmission/DeNovo//"), "DeNovo_Clustering.tiff"), units="px", width=(6*625), height=(6*520), res=600)

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
  xlab("SES Indicators") + ylab(expression(italic("De Novo")~"Genes")) + ggtitle("") + theme(panel.background = element_rect(colour = "grey70", size = 1.5)) +
  theme(axis.title = element_text(size = 13, face = "bold", family = "Calibri")) +
  theme(axis.text.y = element_markdown(size=12, face = "bold", family = "Calibri")) +
  theme(axis.text.x = element_text(size=12,face = "bold" ,family = "Calibri", angle = 30,hjust = 1)) +
  theme(legend.text=element_text(size=11, family = "Calibri")) +
  theme(legend.text.align = 0) +
  theme(legend.title = element_text(size = 12, face = "bold", family = "Calibri")) +
  theme(legend.position=c(1.325,0.55)) +
  theme(plot.margin = unit(c(5, 120 , 5, 5), "pt")) +
  scale_x_discrete(drop = FALSE)  + 
  scale_y_discrete(limits = rev) 
p

dev.off()
