#### This file is a script that is responsible for enrichment of the  
#### results in the significant clusters

#### Prerequisites: 1. WGCNA_Analysis.R

#### Follow up files: None

## Setup env
library(here)
dev.off()
source(here("R/ThirdSubmission/Final/startup.R"))


## Read clustering results
Final = read.table(str_c(here("Res/ThirdSubmission/Clustering/"), "Reactome_Sig_All.txt"), header = T, stringsAsFactors = F, sep = "\t", quote = "")
Final$Treatment = 0
Final$Signature = 0
Final$ClusterID = 0
Final$Regulation = 0

for (i in 1:length(Final$Cluster)) {
  Final$Treatment[[i]] = strsplit(Final$Cluster[[i]], "--")[[1]][[2]]
  Final$Signature[[i]] = strsplit(Final$Cluster[[i]], "--")[[1]][[1]]
  Final$ClusterID[[i]] = strsplit(Final$Cluster[[i]], "--")[[1]][[3]]
  Final$Regulation[[i]] = strsplit(Final$Cluster[[i]], "--")[[1]][[4]]
}
#Final = Final[which(Final$FinalCount>2),]
Final = Final[order(Final$Level1, -Final$FinalCount, Final$p.adjust),]
Final = Final[!duplicated(Final[,c(1,13)]),] #Cluster, Group, Show

tab = data.frame(Signature = Final$Signature, Treatment = Final$Treatment, Regulation = Final$Regulation,
                 Pathway = Final$Show,
                 Adj.P.Val = Final$p.adjust,
                 CollapsedCramerV = Final$MaxC,
                 CollapsedCount = Final$FinalCount, Genes = Final$GenesSymbol)


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
tabw = tab[tab$Regulation=="Upregulated",]
openxlsx::write.xlsx(tabw, str_c(here("Res/ThirdSubmission/Clustering/"), "Reactome_Sig_Up.xlsx"), overwrite = T, rowNames = F, colNames = T)

tabw = tab[tab$Regulation=="Downregulated",]
openxlsx::write.xlsx(tabw, str_c(here("Res/ThirdSubmission/Clustering/"), "Reactome_Sig_Down.xlsx"), overwrite = T, rowNames = F, colNames = T)
