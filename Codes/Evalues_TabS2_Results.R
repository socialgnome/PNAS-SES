#### This file is a script that is responsible for creating the table of e-values for
#### results in the significant clusters

#### Prerequisites: 1. Mediation.R

#### Follow up files: None

## Setup env
library(here)
dev.off()
rm(list=ls(all=TRUE))
library(stringr)
library(EValue)

## Load clustering results - Up 
res = readRDS(str_c(here("Res/ThirdSubmission/Mediation/"), "Up_1KI_parallel.rds"))
treatment = res[[1]]
signature = res[[2]]
mediator = res[[3]]
cluster = res[[4]]
medres = res[[5]]

allres = c()
for (i in 1:length(treatment)) {
  df = data.frame(Treatment = treatment[[i]], Signature = signature[[i]], Mediator = mediator[[i]], Cluster= cluster[[i]],
                  Regulation = "Upregulated", Estimate = medres[[i]][3], SE = medres[[i]][8], SD = medres[[i]][[11]],p = medres[[i]][2])
  allres = rbind(allres,df)
}


## Load clustering results - Down 
res = readRDS(str_c(here("Res/ThirdSubmission/Mediation/"), "Down_1KI_parallel.rds"))
treatment = res[[1]]
signature = res[[2]]
mediator = res[[3]]
cluster = res[[4]]
medres = res[[5]]

for (i in 1:length(treatment)) {
  df = data.frame(Treatment = treatment[[i]], Signature = signature[[i]], Mediator = mediator[[i]], Cluster= cluster[[i]],
                  Regulation = "Downregulated", Estimate = medres[[i]][3], SE = medres[[i]][8], SD = medres[[i]][[11]],p = medres[[i]][2])
  allres = rbind(allres,df)
}

## Computing adjusted P-values
allres$adjP = 0

for (i in 1:length(allres$Treatment)) {
  sel = which(allres$Treatment==allres$Treatment[i] & allres$Mediator==allres$Mediator[i] & allres$Regulation==allres$Regulation[i])
  allres$adjP[sel] = p.adjust(allres$p[sel])
}
allres = allres[allres$adjP<0.05,]

## Evalue calculations
allres$Evalue = 0
for (i in 1:length(allres$Treatment)) {
  allres$Evalue[i] = evalues.OLS(est = allres$Estimate[i], se = allres$SE[i], sd = allres$SD[i])[2,1]
}

allres$AvE = 0
for (i in 1:length(allres$Treatment)) {
  sel = which(allres$Treatment==allres$Treatment[i] & allres$Signature==allres$Signature[i] & allres$Mediator==allres$Mediator[i]  & allres$Regulation==allres$Regulation[i])
  allres$AvE[i] = mean(allres$Evalue[sel])
}

allres = allres[!duplicated(allres[,c(1,2,3,5)]),]

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
                               Signature == "inflam1k_mRNA" ~ "1KI"),
         Mediator = case_when(Mediator == "stress_perceived" ~ "Stress Perceived",
                              Mediator == "w5bmi" ~ "BMI",
                              Mediator == "bills" ~ "Bills",
                              Mediator == "currentsmoke" ~ "Current Smoke",
                              Mediator == "insurance_lack" ~ "Lack of Insurance")
  )

allres$Treatment = factor(allres$Treatment, levels = c("SES Composite","Education","Income","Occupation","Subjective Social Status"))
allres$Signature = factor(allres$Signature, levels = rev(c("Rheumatoid Arthritis","Hypertension","Diabetes","Depression","CVD",
                                                           "COPD","CKD","Asthma","Aortic Aneurysm","Alzheimers",
                                                           "1KI")))

for (j in 1:length(unique(allres$Mediator))) {
  
  res = allres[allres$Regulation=="Upregulated" & allres$Mediator==unique(allres$Mediator)[j],]
  if (length(res$Treatment)>0) {
    tabup = matrix(NA, nrow = length(levels(res$Signature)), ncol = length(levels(res$Treatment)))
    for (i in 1:length(res$Treatment)) {
      tabup[which(levels(res$Signature)==res$Signature[i]), which(levels(res$Treatment)==res$Treatment[i])] = res$AvE[i]
    }
    rownames(tabup) = levels(res$Signature)
    colnames(tabup) = levels(res$Treatment)
    tabup = as.data.frame(tabup)
    openxlsx::write.xlsx(tabup, str_c(here("Res/ThirdSubmission/Tables/"),"TabS2_1KI_Up_",unique(allres$Mediator)[j],".xlsx"),overwrite = T, colNames = T, rowNames = T)
  }
  if (length(res$Treatment)>0) {
    res = allres[allres$Regulation=="Downregulated" & allres$Mediator==unique(allres$Mediator)[j],]
    tabdown = matrix(NA, nrow = length(levels(res$Signature)), ncol = length(levels(res$Treatment)))
    for (i in 1:length(res$Treatment)) {
      tabdown[which(levels(res$Signature)==res$Signature[i]), which(levels(res$Treatment)==res$Treatment[i])] = res$AvE[i]
    }
    rownames(tabdown) = levels(res$Signature)
    colnames(tabdown) = levels(res$Treatment)
    tabdown = as.data.frame(tabdown)
    openxlsx::write.xlsx(tabdown, str_c(here("Res/ThirdSubmission/Tables/"),"TabS2_1KI_Down_",unique(allres$Mediator)[j],".xlsx"),overwrite = T, colNames = T, rowNames = T)
  }
}




##################### Without 1KI

#### This file is a script that is responsible for creating the table of e-values for
#### results in the significant clusters

#### Prerequisites: 1. Mediation.R

#### Follow up files: None

## Setup env
library(here)
dev.off()
rm(list=ls(all=TRUE))
library(stringr)
library(EValue)

## Load clustering results - Up 
res = readRDS(str_c(here("Res/ThirdSubmission/Mediation/"), "Up_Without_1KI_parallel.rds"))
treatment = res[[1]]
signature = res[[2]]
mediator = res[[3]]
cluster = res[[4]]
medres = res[[5]]

allres = c()
for (i in 1:length(treatment)) {
  df = data.frame(Treatment = treatment[[i]], Signature = signature[[i]], Mediator = mediator[[i]], Cluster= cluster[[i]],
                  Regulation = "Upregulated", Estimate = medres[[i]][3], SE = medres[[i]][8], SD = medres[[i]][[11]],p = medres[[i]][2])
  allres = rbind(allres,df)
}


## Load clustering results - Down 
res = readRDS(str_c(here("Res/ThirdSubmission/Mediation/"), "Down_Without_1KI_parallel.rds"))
treatment = res[[1]]
signature = res[[2]]
mediator = res[[3]]
cluster = res[[4]]
medres = res[[5]]

for (i in 1:length(treatment)) {
  df = data.frame(Treatment = treatment[[i]], Signature = signature[[i]], Mediator = mediator[[i]], Cluster= cluster[[i]],
                  Regulation = "Downregulated", Estimate = medres[[i]][3], SE = medres[[i]][8], SD = medres[[i]][[11]],p = medres[[i]][2])
  allres = rbind(allres,df)
}

## Computing adjusted P-values
allres$adjP = 0

for (i in 1:length(allres$Treatment)) {
  sel = which(allres$Treatment==allres$Treatment[i] & allres$Mediator==allres$Mediator[i] & allres$Regulation==allres$Regulation[i])
  allres$adjP[sel] = p.adjust(allres$p[sel])
}
allres = allres[allres$adjP<0.05,]

## Evalue calculations
allres$Evalue = 0
for (i in 1:length(allres$Treatment)) {
  allres$Evalue[i] = evalues.OLS(est = allres$Estimate[i], se = allres$SE[i], sd = allres$SD[i])[2,1]
}

allres$AvE = 0
for (i in 1:length(allres$Treatment)) {
  sel = which(allres$Treatment==allres$Treatment[i] & allres$Signature==allres$Signature[i] & allres$Mediator==allres$Mediator[i]  & allres$Regulation==allres$Regulation[i])
  allres$AvE[i] = mean(allres$Evalue[sel])
}

allres = allres[!duplicated(allres[,c(1,2,3,5)]),]

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
                               Signature == "inflam1k_mRNA" ~ "1KI"),
         Mediator = case_when(Mediator == "stress_perceived" ~ "Stress Perceived",
                              Mediator == "w5bmi" ~ "BMI",
                              Mediator == "bills" ~ "Bills",
                              Mediator == "currentsmoke" ~ "Current Smoke",
                              Mediator == "insurance_lack" ~ "Lack of Insurance")
  )

allres$Treatment = factor(allres$Treatment, levels = c("SES Composite","Education","Income","Occupation","Subjective Social Status"))
allres$Signature = factor(allres$Signature, levels = rev(c("Rheumatoid Arthritis","Hypertension","Diabetes","Depression","CVD",
                                                           "COPD","CKD","Asthma","Aortic Aneurysm","Alzheimers",
                                                           "1KI")))
allres = droplevels(allres)

for (j in 1:length(unique(allres$Mediator))) {
  
  res = allres[allres$Regulation=="Upregulated" & allres$Mediator==unique(allres$Mediator)[j],]
  if (length(res$Treatment)>0) {
    tabup = matrix(NA, nrow = length(levels(res$Signature)), ncol = length(levels(res$Treatment)))
    for (i in 1:length(res$Treatment)) {
      tabup[which(levels(res$Signature)==res$Signature[i]), which(levels(res$Treatment)==res$Treatment[i])] = res$AvE[i]
    }
    rownames(tabup) = levels(res$Signature)
    colnames(tabup) = levels(res$Treatment)
    tabup = as.data.frame(tabup)
    openxlsx::write.xlsx(tabup, str_c(here("Res/ThirdSubmission/Tables/"),"TabS2_Without_1KI_Up_",unique(allres$Mediator)[j],".xlsx"),overwrite = T, colNames = T, rowNames = T)
  }
  if (length(res$Treatment)>0) {
    res = allres[allres$Regulation=="Downregulated" & allres$Mediator==unique(allres$Mediator)[j],]
    tabdown = matrix(NA, nrow = length(levels(res$Signature)), ncol = length(levels(res$Treatment)))
    for (i in 1:length(res$Treatment)) {
      tabdown[which(levels(res$Signature)==res$Signature[i]), which(levels(res$Treatment)==res$Treatment[i])] = res$AvE[i]
    }
    rownames(tabdown) = levels(res$Signature)
    colnames(tabdown) = levels(res$Treatment)
    tabdown = as.data.frame(tabdown)
    openxlsx::write.xlsx(tabdown, str_c(here("Res/ThirdSubmission/Tables/"),"TabS2_Without_1KI_Down_",unique(allres$Mediator)[j],".xlsx"),overwrite = T, colNames = T, rowNames = T)
  }
}
