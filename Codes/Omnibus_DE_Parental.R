#### This file is a script that is responsible for creating the DE objects for the 
#### omnibus test for subjects without the corresponding diseases in the signatures.

#### Prerequisites: 1. BioMart Download of gene symbols

#### Follow up files: 1. OmnibusTest_Fig1A.R


## Setup env
library(here)
dev.off()
source(here("R/ThirdSubmission/Final/startup.R"))


## Biomart (Downloaded 18.01.22)
sapiens_ensembl <- readRDS(str_c(here("data/"), "sapiens_ensembl.rds"))


## Signatures
signatures = readRDS("~/Projects/PNAS/data/Signatures.rds")


## Read raw counts from all batches for all subjects <- Brandts data
raw = readRDS("~/Data/wave5/all.batches.expression.set.070121.Rds")
temp = pData(raw)
raw = exprs(raw)
genes = rownames(raw)


## Get pheno data for all subjects from the RDS file created using preprocess <- Wenjia's recoded phenotypes
df = readRDS("/home/share/preprocessing/preprocessed_two_batches/all.batches.expression.set.rawcount_waves_05.08.2021.rds")
pheno = pData(df)
rownames(pheno) = pheno$AID
pheno = pheno[rownames(temp),]
all.equal(rownames(temp), pheno$AID)
pheno$AnyFlag = temp$AnyFlag
pheno$batch = temp$batch
all.equal(colnames(raw), pheno$AID)


# Inflammation report from AddHealth
data_input = "~/Data/wave5/"
df = read.xport(str_c(data_input, "bcrp5.xpt"))
pheno$H5INFECT = as.factor(df$H5INFECT[match(pheno$AID, df$AID)])
pheno$H5SUBCLN = as.factor(df$H5SUBCLN[match(pheno$AID, df$AID)])
pheno$H5CRP8 = df$H5CRP8[match(pheno$AID, df$AID)]


# Alcohol use report from AddHealth
df = read_xpt(str_c(data_input, "WAVE5.xpt"))
pheno$alcohol_use = df$H5TO14[match(pheno$AID, df$AID)]


## Wave 5 Biocovariates for flagging of subjects
df = read_sas(str_c(data_input, "w5biocovars.sas7bdat"))
pheno$KITCOND = df$KITCOND[match(pheno$AID, df$AID)]
pheno$TUBECOND = df$TUBECOND[match(pheno$AID, df$AID)]
pheno$BloodQual = df$Q086[match(pheno$AID, df$AID)]
pheno$QualCode = df$QualCode[match(pheno$AID, df$AID)]
pheno$AlqQual = df$AlqQuality[match(pheno$AID, df$AID)]


## Convert gene names to HGNC
genedat = getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"), filters ="ensembl_gene_id", values = genes, mart = sapiens_ensembl)
rawgenedata = data.frame(Ensembl = genes, stringsAsFactors = F)
rawgenedata$HGNC = NA
for (i in 1:length(rawgenedata$Ensembl)) {
  match = which(genedat$ensembl_gene_id==rawgenedata$Ensembl[i])
  if  (!is_empty(match)) {
    rawgenedata$HGNC[i] = genedat$hgnc_symbol[match[1]]
  }
}
keepgenes = which(rawgenedata$HGNC!="" | rawgenedata$HGNC!=NA)
raw = raw[keepgenes,]
genes = genes[keepgenes]
rawgenedata = rawgenedata[keepgenes,]


## Count exploration and data filtering
keepgenes = filterByExpr(raw, min.count = 10, min.pro = 0.05)
raw = raw[keepgenes,]
genes = genes[keepgenes]
rawgenedata = rawgenedata[keepgenes,]
rownames(raw) = rawgenedata$HGNC


## Filter subjects based on QC metrics (AnyFlag for subjects)
sel = which(pheno$AnyFlag==0)
pheno = pheno[sel,]
raw = raw[,sel]


## Recode variables in Pheno data
pheno = pheno %>%
  mutate_at(
    .vars = vars(matches("^edu_p$|^edu_max$")),
    .funs = list(~ .x %>%
                   factor() %>%
                   fct_collapse(
                     "high or less" = "high",
                     "more than high" = c("votec","college","post")
                   ))
  ) %>% 
  mutate_at(vars(c("H5SUBCLN",
                   "H5CRP8",
                   "H5INFECT")), .funs = list(~ ifelse(.x %in% c(0,1,2,3), .x, NA))) %>%
  
  mutate(alcohol_use = case_when(alcohol_use %in% c(0,997) ~ "Never",
                                 alcohol_use %in% c(1,2) ~ "Mild",
                                 alcohol_use %in% c(3,4,5) ~ "Moderate",
                                 alcohol_use > 5 ~ "Severe"                                 
  ))


## Define treatment, control, mediation and disease vars
treatment = c("ses_composite_pp1","SEI_max_p_w12","edu_p","income_pp1_log")

controls = c(
  "sex_interv", "re","age_w5"
  ,"pregnant_biow5", "FastHrs","Plate"
  ,"H5INFECT", "H5SUBCLN", "H5CRP8"
)

imp_vars = c("batch")

mediators = c(
  "stress_perceived",
  "w5bmi",
  "bills",
  "currentsmoke",
  "insurance_lack",
  "alcohol_use"
)


disease = c("diabetes", "heartatk", "H5ID6F", "H5ID6FM", "H5ID6C", "H5ID6CM",
            "H5ID6Q","H5ID6QM", "H5ID6A", "H5ID6AM")

subData = pheno %>% dplyr::select(AID,all_of(treatment), all_of(controls), all_of(imp_vars), all_of(mediators), all_of(disease))
selsubData = pheno %>% dplyr::select(AID,all_of(treatment), all_of(controls), all_of(imp_vars))


## Filter subjects and genes based on missing values
non_missing = complete.cases(selsubData)
subData = subData[non_missing, ]
subData = droplevels(subData)
selsubData = selsubData[non_missing, ]
selsubData = droplevels(selsubData)
pheno = pheno[non_missing, ]

keepgenes = filterByExpr(raw[,non_missing], min.count = 10, min.pro = 0.05)
counts = raw[keepgenes,non_missing]


## Remove subjects with diseases
sub_rem = list()
for (i in 1:length(signatures)) {
  sub_rem[[i]] = NA
}
sub_rem[[1]] = which(subData$heartatk==1) # CVD
sub_rem[[2]] = which(subData$diabetes==1) # Diabetes
sub_rem[[6]] = which(subData$H5ID6Q == 1 | subData$H5ID6QM == 1) # Aortic Aneurysm
sub_rem[[8]] = which(subData$H5ID6F == 1 | subData$H5ID6FM == 1) # Asthma
sub_rem[[9]] = which(subData$H5ID6C == 1 | subData$H5ID6CM == 1) # Hypertension


## Run DE for every treatment - signature combo
for (j in 1:length(signatures)) {
  missing = sub_rem[[j]]
  keep = setdiff(1:length(subData$AID), missing)
  selsub = selsubData[keep, ]
  selsub = droplevels(selsub)
  
  keepgenes = filterByExpr(counts[,keep], min.count = 10, min.pro = 0.05)
  countssel = counts[keepgenes,keep]
  
  dge = DGEList(counts = countssel)
  
  ## Normalization - TMM
  dge = calcNormFactors(dge, method = "TMM")
  
  
  ## Linear modeling
  for (i in 1:length(treatment)) {
    rhs = str_c(c(treatment[[i]],controls), collapse = " + ")
    model_formula = str_c(" ~ ",rhs) %>% as.formula()
    sub = selsub %>% dplyr::select(treatment[[i]], all_of(controls))
    sub = droplevels(sub)
    
    if (class(sub[,c(1)])=="factor") {
      sub[,c(1)] = as.numeric(sub[,c(1)])
    }
    
    design  = model.matrix(model_formula, data = sub)
    
    ## voom based linear modeling
    v = voom(dge, design, plot = T)
    
    ## Batch correct using sva::ComBat
    v$E = ComBat(v$E, selsub$batch)
    
    ## DE
    fit = lmFit(v$E, design) 
    fit = eBayes(fit, trend = T)
    
    res = topTable(fit, coef = treatment[[i]], n= Inf) %>%
      rownames_to_column(var = "gene")
    res$treatment = treatment[[i]]
    saveRDS(res, str_c(here("Res/ThirdSubmission/Omnibus/Parental/"), names(signatures)[[j]], "_", treatment[[i]],".rds"))
  }
  
}
