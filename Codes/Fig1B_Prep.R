#### This file is a script that is responsible for the preparation for plotting
#### Fig1B 

#### Prerequisites: 1. WGCNA_Analysis.R

#### Follow up files: 1. Fig1B.R 

## Setup env

library(here)
dev.off()
source(here("R/ThirdSubmission/Final/startup.R"))

## Read Clustering results
res = readRDS(str_c(here("Res/ThirdSubmission/"), "Clustering_Fig1B.rds"))
fig = rbind(data.frame(Signature = res$Signature, Treatment = res$Treatment, 
                       TotalGenes = res$Total_Genes, TotalClusters = res$Total_Clusters,
                       SelGenes = res$Up_Genes, SelClusters = res$Up_Clusters, SelSigGenes = res$Up_Sig_Genes,
                       SelSigClusters = res$Up_Sig_Clusters, FDR = res$Up_Min_FDR, Class = res$Class, Reg = "Upregulated"),
            data.frame(Signature = res$Signature, Treatment = res$Treatment, 
                       TotalGenes = res$Total_Genes, TotalClusters = res$Total_Clusters,
                       SelGenes = res$Down_Genes, SelClusters = res$Down_Clusters, SelSigGenes = res$Down_Sig_Genes,
                       SelSigClusters = res$Down_Sig_Clusters, FDR = res$Down_Min_FDR, Class = res$Class, Reg = "Downregulated")
)

fig = fig[fig$FDR<0.05,]

## Mutate data frame for Plotting
fig$colInstance = 0
for (i in 1:length(fig$Signature)) {
  fig$colInstance[i] = length(which(fig$Signature==fig$Signature[i] & fig$Treatment==fig$Treatment[i] & fig$Class!=fig$Class[i]))
}

fig = fig %>%
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

fig = fig[order(fig$Class),]
fig$YAxisLabel = 0
fig$YAxisLabel_A = 0
fig$YAxisLabel_B = 0
for (i in 1:length(fig$Signature)) {
  fig$YAxisLabel[i] = paste("( <span style='color:#01665e'>",head(fig$TotalClusters[which(fig$Signature==fig$Signature[i])],n=1), "</span> | <span style='color:#ff7f0e'>", tail(fig$TotalClusters[which(fig$Signature==fig$Signature[i])],n=1), "</span> )", sep = "")
  fig$YAxisLabel_A[i] = paste("(",head(fig$TotalClusters[which(fig$Signature==fig$Signature[i])],n=1), ")", sep = "")
  fig$YAxisLabel_B[i] = paste("(",tail(fig$TotalClusters[which(fig$Signature==fig$Signature[i])],n=1), ")", sep = "")
}
fig$YAxisLabel[which(fig$Signature=="1KI")]="( <span style='color:#01665e'>25</span> )" 
fig$YAxisLabel = paste(fig$Signature,"<br/>",fig$YAxisLabel, sep = "")
fig$YAxisLabel_A = paste(fig$Signature,"<br/>",fig$YAxisLabel_A, sep = "")
fig$YAxisLabel_B = paste(fig$Signature,"<br/>",fig$YAxisLabel_B, sep = "")

fig$Treatment = factor(fig$Treatment, levels = c("SES Composite","Education","Income","Occupation","Subjective Social Status"))
fig$YAxisLabel = factor(fig$YAxisLabel, levels = c("Rheumatoid Arthritis<br/>( <span style='color:#01665e'>14</span> | <span style='color:#ff7f0e'>12</span> )",
                                                   "Hypertension<br/>( <span style='color:#01665e'>18</span> | <span style='color:#ff7f0e'>14</span> )",
                                                   "Diabetes<br/>( <span style='color:#01665e'>18</span> | <span style='color:#ff7f0e'>17</span> )",
                                                   "Depression<br/>( <span style='color:#01665e'>21</span> | <span style='color:#ff7f0e'>19</span> )",
                                                   "CVD<br/>( <span style='color:#01665e'>12</span> | <span style='color:#ff7f0e'>11</span> )",
                                                   "COPD<br/>( <span style='color:#01665e'>15</span> | <span style='color:#ff7f0e'>14</span> )",
                                                   "CKD<br/>( <span style='color:#01665e'>21</span> | <span style='color:#ff7f0e'>21</span> )",
                                                   "Asthma<br/>( <span style='color:#01665e'>14</span> | <span style='color:#ff7f0e'>13</span> )",
                                                   "Alzheimers<br/>( <span style='color:#01665e'>18</span> | <span style='color:#ff7f0e'>17</span> )",
                                                   "Aortic Aneurysm<br/>( <span style='color:#01665e'>7</span> | <span style='color:#ff7f0e'>7</span> )",
                                                   "1KI<br/>( <span style='color:#01665e'>25</span> )"))

fig$YAxisLabel_A = factor(fig$YAxisLabel_A, levels = c("Rheumatoid Arthritis<br/>(14)",
                                                       "Hypertension<br/>(18)",
                                                       "Diabetes<br/>(18)",
                                                       "Depression<br/>(21)",
                                                       "CVD<br/>(12)",
                                                       "COPD<br/>(15)",
                                                       "CKD<br/>(21)",
                                                       "Asthma<br/>(14)",
                                                       "Alzheimers<br/>(18)",
                                                       "Aortic Aneurysm<br/>(7)",
                                                       "1KI<br/>(25)"))

fig$YAxisLabel_B = factor(fig$YAxisLabel_B, levels = c("Rheumatoid Arthritis<br/>(12)",
                                                       "Hypertension<br/>(14)",
                                                       "Diabetes<br/>(17)",
                                                       "Depression<br/>(19)",
                                                       "CVD<br/>(11)",
                                                       "COPD<br/>(13)",
                                                       "CKD<br/>(21)",
                                                       "Asthma<br/>(13)",
                                                       "Alzheimers<br/>(17)",
                                                       "Aortic Aneurysm<br/>(7)",
                                                       "1KI<br/>(25)"))

saveRDS(fig, str_c(here("Res/ThirdSubmission/"), "Plot_Fig1B.rds"))