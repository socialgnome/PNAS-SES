library(here)
dev.off()
source(here("R/ThirdSubmission/Final/startup.R"))

## Define treatment, signatures
signatures = readRDS("~/Projects/PNAS/data/Signatures.rds")

treatment = c("ses_sss_composite","sss_5","SEI_ff5",
              "edu_max","income_hh_ff5")


df = readRDS(str_c(here("Res/ThirdSubmission/Omnibus/"), names(signatures)[[5]], "_", treatment[[1]],".rds"))

#### 14163 genes = 25 Modules
#### 368 down and 482 upregulated
#### 14163 Genes = 850 Significant
#### 126 genes in alzeheimers/ 101 have expression data = 18 clusters
#### 51 genes upregulated and 50 genes downregulated (From the DE analysis)
#### 15 clusters are up and 11 clusters are down
#### 32 and 35 up/down sig genes in 3/4 sig clusters


df$Neg = -log10(df$adj.P.Val)
df$class = "Non significant"
df$class[df$adj.P.Val<0.05 & df$logFC>0] = "Upregulated"
df$class[df$adj.P.Val<0.05 & df$logFC<0] = "Downregulated"
data = df
tiff(str_c(here("Res/ThirdSubmission/Figures/"), "Volcano_SES_Eg.tiff"), units="px", width=(6*750), height=(6*400), res=600)
p = ggplot(data, aes(x = logFC, y = Neg, color = class, fill = class)) +
  geom_point(shape = 21, size = 3, alpha = 0.7) +
  theme_bw(base_size = 12) +
  scale_color_manual(values = c("#4379a7","#767676","#c93311"), name = "Direction\nof change") +
  scale_fill_manual(values = c("#4379a7","#767676","#c93311"), name = "Direction\nof change") +
  xlab(expression(Log[2]~fold~change)) + ylab(expression(-Log[10]~italic(p))) + 
  ggtitle("SES Composite") +
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