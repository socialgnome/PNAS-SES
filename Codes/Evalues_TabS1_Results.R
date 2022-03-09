#### This file is a script that is responsible for creating the table of e-values for
#### results in the significant clusters

#### Prerequisites: 1. Evalues_TabS1.R

#### Follow up files: None

## Setup env
library(here)
dev.off()
rm(list=ls(all=TRUE))
library(stringr)
library(EValue)

## Load clustering results - Up 
res = readRDS(str_c(here("Res/ThirdSubmission/Clustering/"), "Evalues_Up_With_1KI.rds"))
treatment = res[[1]]
signature = res[[2]]
mediator = res[[3]]
cluster = res[[4]]
eres = res[[5]]

allres = c()
for (i in 1:length(treatment)) {
  df = data.frame(Treatment = treatment[[i]], Signature = signature[[i]], Mediator = mediator[[i]], Cluster= cluster[[i]],
                  Regulation = "Upregulated", Estimate = eres[[i]][1], SE = eres[[i]][2], SD = eres[[i]][[3]])
  allres = rbind(allres,df)
}

## Load clustering results - Down 
res = readRDS(str_c(here("Res/ThirdSubmission/Clustering/"), "Evalues_Down_With_1KI.rds"))
treatment = res[[1]]
signature = res[[2]]
mediator = res[[3]]
cluster = res[[4]]
eres = res[[5]]

for (i in 1:length(treatment)) {
  df = data.frame(Treatment = treatment[[i]], Signature = signature[[i]], Mediator = mediator[[i]], Cluster= cluster[[i]],
                  Regulation = "Downregulated", Estimate = eres[[i]][1], SE = eres[[i]][2], SD = eres[[i]][[3]])
  allres = rbind(allres,df)
}


## Evalue calculations
allres$Evalue = 0
for (i in 1:length(allres$Treatment)) {
  allres$Evalue[i] = evalues.OLS(est = allres$Estimate[i], se = allres$SE[i], sd = allres$SD[i])[2,1]
}

allres$AvE = 0
for (i in 1:length(allres$Treatment)) {
  sel = which(allres$Treatment==allres$Treatment[i] & allres$Signature==allres$Signature[i]  & allres$Regulation==allres$Regulation[i])
  allres$AvE[i] = mean(allres$Evalue[sel])
}

allres = allres[!duplicated(allres[,c(1,2,5)]),]

library(dplyr)
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
allres$Signature = factor(allres$Signature, levels = rev(c("Rheumatoid Arthritis","Hypertension","Diabetes","Depression","CVD",
                                                       "COPD","CKD","Asthma","Aortic Aneurysm","Alzheimers",
                                                       "1KI")))

res = allres[allres$Regulation=="Upregulated",]
tabup = matrix(NA, nrow = length(levels(res$Signature)), ncol = length(levels(res$Treatment)))
for (i in 1:length(res$Treatment)) {
  tabup[which(levels(res$Signature)==res$Signature[i]), which(levels(res$Treatment)==res$Treatment[i])] = res$AvE[i]
}
rownames(tabup) = levels(res$Signature)
colnames(tabup) = levels(res$Treatment)
tabup = as.data.frame(tabup)

res = allres[allres$Regulation=="Downregulated",]
tabdown = matrix(NA, nrow = length(levels(res$Signature)), ncol = length(levels(res$Treatment)))
for (i in 1:length(res$Treatment)) {
  tabdown[which(levels(res$Signature)==res$Signature[i]), which(levels(res$Treatment)==res$Treatment[i])] = res$AvE[i]
}
rownames(tabdown) = levels(res$Signature)
colnames(tabdown) = levels(res$Treatment)
tabdown = as.data.frame(tabdown)

openxlsx::write.xlsx(tabup, str_c(here("Res/ThirdSubmission/Tables/"),"TabS1_1KI_Up.xlsx"),overwrite = T, colNames = T, rowNames = T)
openxlsx::write.xlsx(tabdown, str_c(here("Res/ThirdSubmission/Tables/"),"TabS1_1KI_Down.xlsx"),overwrite = T, colNames = T, rowNames = T)



############# Without 1KI
#### This file is a script that is responsible for creating the table of e-values for
#### results in the significant clusters

#### Prerequisites: 1. Evalues_TabS1.R

#### Follow up files: None

## Setup env
library(here)
dev.off()
rm(list=ls(all=TRUE))
library(stringr)
library(EValue)

## Load clustering results - Up 
res = readRDS(str_c(here("Res/ThirdSubmission/Clustering/"), "Evalues_Up_Without_1KI.rds"))
treatment = res[[1]]
signature = res[[2]]
mediator = res[[3]]
cluster = res[[4]]
eres = res[[5]]

allres = c()
for (i in 1:length(treatment)) {
  df = data.frame(Treatment = treatment[[i]], Signature = signature[[i]], Mediator = mediator[[i]], Cluster= cluster[[i]],
                  Regulation = "Upregulated", Estimate = eres[[i]][1], SE = eres[[i]][2], SD = eres[[i]][[3]])
  allres = rbind(allres,df)
}

## Load clustering results - Down 
res = readRDS(str_c(here("Res/ThirdSubmission/Clustering/"), "Evalues_Down_Without_1KI.rds"))
treatment = res[[1]]
signature = res[[2]]
mediator = res[[3]]
cluster = res[[4]]
eres = res[[5]]

for (i in 1:length(treatment)) {
  df = data.frame(Treatment = treatment[[i]], Signature = signature[[i]], Mediator = mediator[[i]], Cluster= cluster[[i]],
                  Regulation = "Downregulated", Estimate = eres[[i]][1], SE = eres[[i]][2], SD = eres[[i]][[3]])
  allres = rbind(allres,df)
}


## Evalue calculations
allres$Evalue = 0
for (i in 1:length(allres$Treatment)) {
  allres$Evalue[i] = evalues.OLS(est = allres$Estimate[i], se = allres$SE[i], sd = allres$SD[i])[2,1]
}

allres$AvE = 0
for (i in 1:length(allres$Treatment)) {
  sel = which(allres$Treatment==allres$Treatment[i] & allres$Signature==allres$Signature[i]  & allres$Regulation==allres$Regulation[i])
  allres$AvE[i] = mean(allres$Evalue[sel])
}

allres = allres[!duplicated(allres[,c(1,2,5)]),]

library(dplyr)
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
allres$Signature = factor(allres$Signature, levels = rev(c("Rheumatoid Arthritis","Hypertension","Diabetes","Depression","CVD",
                                                           "COPD","CKD","Asthma","Aortic Aneurysm","Alzheimers",
                                                           "1KI")))
allres = droplevels(allres)
res = allres[allres$Regulation=="Upregulated",]
tabup = matrix(NA, nrow = length(levels(res$Signature)), ncol = length(levels(res$Treatment)))
for (i in 1:length(res$Treatment)) {
  tabup[which(levels(res$Signature)==res$Signature[i]), which(levels(res$Treatment)==res$Treatment[i])] = res$AvE[i]
}
rownames(tabup) = levels(res$Signature)
colnames(tabup) = levels(res$Treatment)
tabup = as.data.frame(tabup)

res = allres[allres$Regulation=="Downregulated",]
tabdown = matrix(NA, nrow = length(levels(res$Signature)), ncol = length(levels(res$Treatment)))
for (i in 1:length(res$Treatment)) {
  tabdown[which(levels(res$Signature)==res$Signature[i]), which(levels(res$Treatment)==res$Treatment[i])] = res$AvE[i]
}
rownames(tabdown) = levels(res$Signature)
colnames(tabdown) = levels(res$Treatment)
tabdown = as.data.frame(tabdown)

openxlsx::write.xlsx(tabup, str_c(here("Res/ThirdSubmission/Tables/"),"TabS1_Without_1KI_Up.xlsx"),overwrite = T, colNames = T, rowNames = T)
openxlsx::write.xlsx(tabdown, str_c(here("Res/ThirdSubmission/Tables/"),"TabS1_Without_1KI_Down.xlsx"),overwrite = T, colNames = T, rowNames = T)
