#### This file is a script that is responsible for creating the object for the 
#### creation of Fig1A in the PNAS manuscript

#### Prerequisites: 1. Omnibus_DE.R

#### Follow up files: 1. Fig1A.R

## Setup env
library(here)
dev.off()
source(here("R/ThirdSubmission/Final/startup.R"))

## Define treatment, signatures
signatures = readRDS("~/Projects/PNAS/data/Signatures.rds")

treatment = c("ses_sss_composite","sss_5","SEI_ff5",
              "edu_max","income_hh_ff5")


## With 1K - Analysis
res_1k = data.frame(Treatment = rep(treatment, each = length(signatures)), 
                    Signature = rep(names(signatures), length(treatment)),
                    P = 0, Class = "1KI Genes")

p = c()
t = c()
for (i in 1:length(treatment)) {
  
  for (j in 1:length(signatures)) {
    
    df = readRDS(str_c(here("Res/ThirdSubmission/Omnibus/"), names(signatures)[[j]], "_", treatment[[i]],".rds"))
    
    genes = signatures[[j]]
    sel = df[df$gene %in% genes,]
    p = c(p, min(sel$adj.P.Val))
    t = c(t, sel$t[which(sel$adj.P.Val==min(sel$adj.P.Val))[1]])
  }
  
}
res_1k$P = p
res_1k$t = t

## Without 1K - Analysis
inflam_genes = signatures[[which(names(signatures)=="inflam1k_mRNA")]]
signatures = signatures[-c(which(names(signatures)=="inflam1k_mRNA"))]
res_wo_1k = data.frame(Treatment = rep(treatment, each = length(signatures)), 
                       Signature = rep(names(signatures), length(treatment)),
                       P = 0, Class = "Without 1KI Genes")

p = c()
t = c()
for (i in 1:length(treatment)) {
  
  for (j in 1:length(signatures)) {
    
    df = readRDS(str_c(here("Res/ThirdSubmission/Omnibus/"), names(signatures)[[j]], "_", treatment[[i]],".rds"))
    
    genes = setdiff(signatures[[j]], inflam_genes)
    sel = df[df$gene %in% genes,]
    p = c(p, min(sel$adj.P.Val))
    t = c(t, sel$t[which(sel$adj.P.Val==min(sel$adj.P.Val))[1]])
  }
  
}
res_wo_1k$P = p
res_wo_1k$t = t

## All res
res = rbind(res_1k, res_wo_1k)
res = res %>%
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
                               Signature == "Aortic_Aneurysm_mRNA" ~ "Aortic_Aneurysm",
                               Signature == "inflam1k_mRNA" ~ "1KI")
  )
res$Treatment = factor(res$Treatment, levels = c("SES Composite","Education","Income","Occupation","Subjective Social Status"))
res$Signature = factor(res$Signature, levels = c("Rheumatoid Arthritis","Hypertension","Diabetes","Depression","CVD",
                                                 "COPD","CKD","Asthma","Aortic Aneurysm","Alzheimers",
                                                 "1KI"))
res = res[res$P<0.05,]

saveRDS(res, str_c(here("Res/ThirdSubmission/"), "Omnibus_Fig1A.rds"))

openxlsx::write.xlsx(res, str_c(here("Res/ThirdSubmission/Tables/"), "Omnibus_Results_for_Fig1A.xlsx"),overwrite = T, colNames = T, colWidths =25)
