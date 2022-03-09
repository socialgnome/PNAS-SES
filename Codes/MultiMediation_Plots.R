#### This file is a script that is responsible for plotting of the mediational 
#### results in the significant clusters

#### Prerequisites: 1. Mediation.R

#### Follow up files: None

## Setup env
library(here)
dev.off()
source(here("R/ThirdSubmission/Final/startup.R"))

## Load clustering results - Up 
res = readRDS(str_c(here("Res/ThirdSubmission/Multimediation/"), "Up_1KI_parallel.rds"))
treatment = res[[1]]
signature = res[[2]]
mediator = res[[3]]
cluster = res[[4]]
medres = res[[5]]

allres = c()
for (i in 1:length(treatment)) {
  df = data.frame(Treatment = treatment[[i]], Signature = signature[[i]], Cluster= cluster[[i]],
                  Regulation = "Upregulated", PropMed = medres[[i]][[2,2]], p = medres[[i]][[2,5]])
  allres = rbind(allres,df)
}


## Load clustering results - Down
res = readRDS(str_c(here("Res/ThirdSubmission/Multimediation/"), "Down_1KI_parallel.rds"))
treatment = res[[1]]
signature = res[[2]]
mediator = res[[3]]
cluster = res[[4]]
medres = res[[5]]

for (i in 1:length(treatment)) {
  df = data.frame(Treatment = treatment[[i]], Signature = signature[[i]], Cluster= cluster[[i]],
                  Regulation = "Downregulated", PropMed = medres[[i]][[2,2]], p = medres[[i]][[2,5]])
  allres = rbind(allres,df)
}

## Computing adjusted P-values
allres$adjP = 0
allres$FDR = 0
allres$AvMed = 0

for (i in 1:length(allres$Treatment)) {
  sel = which(allres$Treatment==allres$Treatment[i]  & allres$Regulation==allres$Regulation[i])
  allres$adjP[sel] = p.adjust(allres$p[sel])
}
allres = allres[allres$adjP<0.05,]


## For significantly mediated clusters, combine results
for (i in 1:length(allres$Treatment)) {
  sel = which(allres$Treatment==allres$Treatment[i] & allres$Signature==allres$Signature[i] & allres$Regulation==allres$Regulation[i])
  allres$FDR[i] = combine.test(allres$adjP[sel], method = c("fisher"))
  allres$AvMed[i] = median(allres$PropMed[sel])
}

tab = allres
allres = allres[!duplicated(allres[,c(1,2,4)]),]

## Prepare data for plotting
allres = allres %>%
  mutate(Treatment  = case_when(Treatment == "ses_sss_composite" ~ "SES Composite",
                                Treatment == "sss_5" ~ "Subjective Social Status", 
                                Treatment == "SEI_ff5" ~ "Occupation",
                                Treatment == "edu_max" ~ "Education",
                                Treatment == "income_hh_ff5" ~ "Income"),
         Signature = case_when(Signature == "CVD_mRNA" ~ "CVD",
                               Signature == "diabetes_mRNA" ~ "Diabetes",
                               Signature == "Rheumatoid_Arthritis_mRNA" ~ "Rheumatoid Arthritis",
                               Signature == "Alzheimers_mRNA" ~ "Alzheimers",
                               Signature == "COPD_mRNA" ~ "COPD",
                               Signature == "Asthma_mRNA" ~ "Asthma",
                               Signature == "Hypertension_mRNA" ~ "Hypertension",
                               Signature == "Depression_mRNA" ~ "Depression",
                               Signature == "CKD_mRNA" ~ "CKD",
                               Signature == "Aortic_Aneurysm_mRNA" ~ "Aortic Aneurysm",
                               Signature == "inflam1k_mRNA" ~ "1KI")
  )
allres$Treatment = factor(allres$Treatment, levels = c("SES Composite","Education","Income","Occupation","Subjective Social Status"))
allres$Signature = factor(allres$Signature, levels = c("Rheumatoid Arthritis","Hypertension","Diabetes","Depression","CVD",
                                                       "COPD","CKD","Asthma","Aortic Aneurysm","Alzheimers",
                                                       "1KI"))

## Text results
tab = tab %>%
  mutate(Treatment  = case_when(Treatment == "ses_sss_composite" ~ "SES Composite",
                                Treatment == "sss_5" ~ "Subjective Social Status", 
                                Treatment == "SEI_ff5" ~ "Occupation",
                                Treatment == "edu_max" ~ "Education",
                                Treatment == "income_hh_ff5" ~ "Income"),
         Signature = case_when(Signature == "CVD_mRNA" ~ "CVD",
                               Signature == "diabetes_mRNA" ~ "Diabetes",
                               Signature == "Rheumatoid_Arthritis_mRNA" ~ "Rheumatoid Arthritis",
                               Signature == "Alzheimers_mRNA" ~ "Alzheimers",
                               Signature == "COPD_mRNA" ~ "COPD",
                               Signature == "Asthma_mRNA" ~ "Asthma",
                               Signature == "Hypertension_mRNA" ~ "Hypertension",
                               Signature == "Depression_mRNA" ~ "Depression",
                               Signature == "CKD_mRNA" ~ "CKD",
                               Signature == "Aortic_Aneurysm_mRNA" ~ "Aortic Aneurysm",
                               Signature == "inflam1k_mRNA" ~ "1KI")
  )

tab$Treatment = factor(tab$Treatment, levels = c("SES Composite","Education","Income","Occupation","Subjective Social Status"))
tab$Signature = factor(tab$Signature, levels = c("Rheumatoid Arthritis","Hypertension","Diabetes","Depression","CVD",
                                                 "COPD","CKD","Asthma","Aortic Aneurysm","Alzheimers",
                                                 "1KI"))

tab_signature = data.frame(Signature = levels(tab$Signature), Min = sapply(split(tab$PropMed, tab$Signature), min),
                           Max = sapply(split(tab$PropMed, tab$Signature), max), Median = sapply(split(tab$PropMed, tab$Signature), median))

tab_treatment = data.frame(Signature = levels(tab$Treatment), Min = sapply(split(tab$PropMed, tab$Treatment), min),
                           Max = sapply(split(tab$PropMed, tab$Treatment), max), Median = sapply(split(tab$PropMed, tab$Treatment), median))

allres$pval2 = -log10(allres$FDR)
allres$pval2[which(allres$pval2>-log10(0.001))] = -log10(0.001)

allres$AvMed  = round(allres$AvMed*100, digits = 1)
allres$AvMed = format(allres$AvMed,nsmall = 1)

res = allres

tiff(str_c(here("Res/ThirdSubmission/Figures/"), "MultiMed.tiff"), units="px", width=(6*750), height=(6*670), res=600)
p = ggplot(res[res$Regulation=="Downregulated",], aes(x = Treatment,y  = Signature, color = Regulation, fill = Regulation,  size = pval2)) +
  geom_point(alpha  = 0.90, shape="\u25D6", stroke = 1,position = position_nudge(x = -0.13)) +
  geom_point(data = res[res$Regulation=="Upregulated",], shape="\u25D7",inherit.aes = T, position = position_nudge(x = +0.13), alpha = 0.90, stroke = 1,show.legend = F) +
  
  geom_text(data = res[res$Regulation=="Downregulated",], aes(x = Treatment,y  = Signature, label = AvMed), color = "black", family = "Calibri", fontface="bold",alpha=1, size=3.5, inherit.aes = FALSE,position = position_nudge(x = -0.23)) +
  geom_text(data = res[res$Regulation=="Upregulated",], aes(x = Treatment,y  = Signature, label = AvMed), color = "black", family = "Calibri", fontface="bold",alpha=1, size=3.5, inherit.aes = FALSE,position = position_nudge(x = +0.18)) +
  
  theme_bw(base_size = 12) +
  scale_color_manual(values = c("#4379a7","#c93311"), name = "Direction\nof change",
                     #guide = guide_legend(override.aes = list(shape = c("\u25D0","\u25D1"), size =7))) +
                     guide = guide_legend(override.aes = list(shape = c("\u25D6","\u25D7"), size =7))) +
  scale_fill_manual(values = c("#4379a7","#c93311"), name = "Direction\nof change",
                    #guide = guide_legend(override.aes = list(shape = c("\u25D0","\u25D1"), size =7))) +
                    guide = guide_legend(override.aes = list(shape = c("\u25D6","\u25D7"), size =7))) +
  scale_size_continuous(range = c(6,12),
                        name = "Aggregated\nadjusted\np-value", 
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
  theme(legend.position=c(1.2,0.55)) +
  theme(plot.margin = unit(c(5, 120 , 5, 5), "pt")) +
  scale_x_discrete(drop = FALSE)  +
  scale_y_discrete(limits = rev)
p
dev.off()

