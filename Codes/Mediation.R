rm(list=ls(all=TRUE))
dev.off()
library(here)
library(mediation)
library(stringr)
library(purrr)
library(foreach)
library(parallel)
library(doParallel)
library(doSNOW)
library(progress)
library(pbmcapply)
closeAllConnections()


## Read clustering results
res = readRDS(str_c(here("Res/ThirdSubmission/Clustering"), "/1KI.rds"))

## Load expression data
dat = readRDS(str_c(here("data"), "/Filtered_ExpressionSet_Clustering_SelSubData.rds"))

## Define treatment, control, mediation and disease vars
subData = Biobase::pData(dat)
treatment = c("ses_sss_composite","sss_5","SEI_ff5",
              "edu_max","income_hh_ff5")

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

subData$alcohol_use = factor(subData$alcohol_use, levels = c("Never", "Mild", "Moderate", "Severe"))
subData[,which(colnames(subData) %in% treatment | colnames(subData) %in% mediators)] <- sapply(subData[,which(colnames(subData) %in% treatment | colnames(subData) %in% mediators)], as.numeric)


## Mediation function and extraction
med_fun = function(treatment, controls, mediator, subjectData, AID, expressionData) {
  ## Data
  subData = subjectData[AID,]
  data = data.frame(Y = expressionData, subData)
  ## fit.x
  mod_formula = as.formula(str_c(mediator," ~ ",str_c(c(treatment,controls), collapse = " + ")))
  fit_x = lm(mod_formula, data)
  ## fit.y
  mod_formula = as.formula(str_c("Y ~ ",str_c(c(treatment,controls, mediator), collapse = " + ")))
  fit_y = lm(mod_formula,data)
  ## Mediate
  m = mediate(fit_x, fit_y, treat = treatment, mediator = mediator)
  ## Extract results
  ACME_sd = m$d0.sims %>% sd
  ADE_sd = m$z0.sims %>% sd
  Total_sd = m$tau.sims %>% sd 
  y_sd = summary(m$model.y)$sigma
  p = m$d1.p
  med_prop = m$n1
  med_ACME = m$d1
  med_ADE = m$z1
  med_ACME_p = m$d1.p
  med_ADE_p = m$z1.p
  med_Total = m$tau.coef
  allres = c(med_prop,p,med_ACME, med_ADE,med_Total, med_ACME_p, med_ADE_p,ACME_sd, ADE_sd,Total_sd, y_sd)
  return(allres)
}


### For upregulated average expression
up_av_expr = res[[3]] ## 55 different combo of treatment - signature
up_av_expr_de = res[[5]]
treat = res[[7]]
signature = res[[8]]

expr = list()
subjectdata = list()
treat_m = list()
signature_m = list()
mediator_m = list()
clus_m = list()
for (i in 1:length(up_av_expr)) {
  sel = which(up_av_expr_de[[i]]$adj.P.Val<0.05)
  if (length(sel)>0) {
    for (j in 1:length(sel)) {
      expr[[length(expr)+1]] = up_av_expr[[i]][sel[j],]
      subjectdata[[length(subjectdata)+1]] = colnames(up_av_expr[[i]])
      treat_m[[length(treat_m)+1]] = treat[[i]]
      signature_m[[length(signature_m)+1]] = signature[[i]]
      clus_m[[length(clus_m)+1]] = rownames(up_av_expr[[i]])[sel[j]]
    }
  }
}
mediator_m = rep(mediators, each = length(expr))
expr = rep(expr, length(mediators))
subjectdata = rep(subjectdata, length(mediators))
treat_m = rep(treat_m, length(mediators))
signature_m = rep(signature_m, length(mediators))
clus_m = rep(clus_m, length(mediators))


### parallel - SNOW 
strt<-Sys.time()
cl <- makeSOCKcluster(16)
registerDoSNOW(cl)
clusterEvalQ(cl, library(stringr))
clusterEvalQ(cl, library(mediation))
pb <- txtProgressBar(min = 1, max = length(expr), style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress=progress)
fullout_up = foreach(i=1:length(expr) , .options.snow = opts) %dopar% {  
  out = med_fun(treat_m[[i]], controls, mediator_m[[i]], subData, subjectdata[[i]], expr[[i]])
  out
}
close(pb)
parallel::stopCluster(cl)
print(Sys.time()-strt)

### Results
allresup = list()
allresup[[1]] = treat_m
allresup[[2]] = signature_m
allresup[[3]] = mediator_m
allresup[[4]] = clus_m
allresup[[5]] = fullout_up

### For Downregulated average expression
down_av_expr = res[[4]] ## 55 different combo of treatment - signature
down_av_expr_de = res[[6]]
treat = res[[7]]
signature = res[[8]]

expr = list()
subjectdata = list()
treat_m = list()
signature_m = list()
mediator_m = list()
clus_m = list()
for (i in 1:length(down_av_expr)) {
  sel = which(down_av_expr_de[[i]]$adj.P.Val<0.05)
  if (length(sel)>0) {
    for (j in 1:length(sel)) {
      expr[[length(expr)+1]] = down_av_expr[[i]][sel[j],]
      subjectdata[[length(subjectdata)+1]] = colnames(down_av_expr[[i]])
      treat_m[[length(treat_m)+1]] = treat[[i]]
      signature_m[[length(signature_m)+1]] = signature[[i]]
      clus_m[[length(clus_m)+1]] = rownames(down_av_expr[[i]])[sel[j]]
    }
  }
}
mediator_m = rep(mediators, each = length(expr))
expr = rep(expr, length(mediators))
subjectdata = rep(subjectdata, length(mediators))
treat_m = rep(treat_m, length(mediators))
signature_m = rep(signature_m, length(mediators))
clus_m = rep(clus_m, length(mediators))


### parallel - SNOW 
strt<-Sys.time()
cl <- makeSOCKcluster(16)
registerDoSNOW(cl)
clusterEvalQ(cl, library(stringr))
clusterEvalQ(cl, library(mediation))
pb <- txtProgressBar(min = 1, max = length(expr), style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress=progress)
fullout_down = foreach(i=1:length(expr) , .options.snow = opts) %dopar% {  
  out = med_fun(treat_m[[i]], controls, mediator_m[[i]], subData, subjectdata[[i]], expr[[i]])
  out
}
close(pb)
parallel::stopCluster(cl)
print(Sys.time()-strt)


### Results
allresdown = list()
allresdown[[1]] = treat_m
allresdown[[2]] = signature_m
allresdown[[3]] = mediator_m
allresdown[[4]] = clus_m
allresdown[[5]] = fullout_down

### Save results
saveRDS(allresup, str_c(here("Res/ThirdSubmission/Mediation"),"/Up_1KI_parallel.rds"))
saveRDS(allresdown, str_c(here("Res/ThirdSubmission/Mediation"),"/Down_1KI_parallel.rds"))



######################
######################
######################

rm(list=ls(all=TRUE))
dev.off()
library(here)
library(mediation)
library(stringr)
library(purrr)
library(foreach)
library(parallel)
library(doParallel)
library(doSNOW)
library(progress)
library(pbmcapply)
closeAllConnections()


## Read clustering results
res = readRDS(str_c(here("Res/ThirdSubmission/Clustering"), "/Without_1KI.rds"))
  
## Load expression data
dat = readRDS(str_c(here("data"), "/Filtered_ExpressionSet_Clustering_SelSubData.rds"))

## Define treatment, control, mediation and disease vars
subData = Biobase::pData(dat)
treatment = c("ses_sss_composite","sss_5","SEI_ff5",
              "edu_max","income_hh_ff5")

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

subData$alcohol_use = factor(subData$alcohol_use, levels = c("Never", "Mild", "Moderate", "Severe"))

subData[,which(colnames(subData) %in% treatment | colnames(subData) %in% mediators)] <- sapply(subData[,which(colnames(subData) %in% treatment | colnames(subData) %in% mediators)], as.numeric)


## Mediation function and extraction
med_fun = function(treatment, controls, mediator, subjectData, AID, expressionData) {
  ## Data
  subData = subjectData[AID,]
  data = data.frame(Y = expressionData, subData)
  ## fit.x
  mod_formula = as.formula(str_c(mediator," ~ ",str_c(c(treatment,controls), collapse = " + ")))
  fit_x = lm(mod_formula, data)
  ## fit.y
  mod_formula = as.formula(str_c("Y ~ ",str_c(c(treatment,controls, mediator), collapse = " + ")))
  fit_y = lm(mod_formula,data)
  ## Mediate
  m = mediate(fit_x, fit_y, treat = treatment, mediator = mediator)
  ## Extract results
  ACME_sd = m$d0.sims %>% sd
  ADE_sd = m$z0.sims %>% sd
  Total_sd = m$tau.sims %>% sd 
  y_sd = summary(m$model.y)$sigma
  p = m$d1.p
  med_prop = m$n1
  med_ACME = m$d1
  med_ADE = m$z1
  med_ACME_p = m$d1.p
  med_ADE_p = m$z1.p
  med_Total = m$tau.coef
  allres = c(med_prop,p,med_ACME, med_ADE,med_Total, med_ACME_p, med_ADE_p,ACME_sd, ADE_sd,Total_sd, y_sd)
  return(allres)
}


### For upregulated average expression
up_av_expr = res[[3]] ## 50 different combo of treatment - signature
up_av_expr_de = res[[5]]
treat = res[[7]]
signature = res[[8]]

expr = list()
subjectdata = list()
treat_m = list()
signature_m = list()
mediator_m = list()
clus_m = list()
for (i in 1:length(up_av_expr)) {
  sel = which(up_av_expr_de[[i]]$adj.P.Val<0.05)
  if (length(sel)>0) {
    for (j in 1:length(sel)) {
      expr[[length(expr)+1]] = up_av_expr[[i]][sel[j],]
      subjectdata[[length(subjectdata)+1]] = colnames(up_av_expr[[i]])
      treat_m[[length(treat_m)+1]] = treat[[i]]
      signature_m[[length(signature_m)+1]] = signature[[i]]
      clus_m[[length(clus_m)+1]] = rownames(up_av_expr[[i]])[sel[j]]
    }
  }
}
mediator_m = rep(mediators, each = length(expr))
expr = rep(expr, length(mediators))
subjectdata = rep(subjectdata, length(mediators))
treat_m = rep(treat_m, length(mediators))
signature_m = rep(signature_m, length(mediators))
clus_m = rep(clus_m, length(mediators))


### parallel - SNOW 
strt<-Sys.time()
cl <- makeSOCKcluster(16)
registerDoSNOW(cl)
clusterEvalQ(cl, library(stringr))
clusterEvalQ(cl, library(mediation))
pb <- txtProgressBar(min = 1, max = length(expr), style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress=progress)
fullout_up = foreach(i=1:length(expr) , .options.snow = opts) %dopar% {  
  out = med_fun(treat_m[[i]], controls, mediator_m[[i]], subData, subjectdata[[i]], expr[[i]])
  out
}
close(pb)
parallel::stopCluster(cl)
print(Sys.time()-strt)

### Results
allresup = list()
allresup[[1]] = treat_m
allresup[[2]] = signature_m
allresup[[3]] = mediator_m
allresup[[4]] = clus_m
allresup[[5]] = fullout_up

### For Downregulated average expression
down_av_expr = res[[4]] ## 50 different combo of treatment - signature
down_av_expr_de = res[[6]]
treat = res[[7]]
signature = res[[8]]

expr = list()
subjectdata = list()
treat_m = list()
signature_m = list()
mediator_m = list()
clus_m = list()
for (i in 1:length(down_av_expr)) {
  sel = which(down_av_expr_de[[i]]$adj.P.Val<0.05)
  if (length(sel)>0) {
    for (j in 1:length(sel)) {
      expr[[length(expr)+1]] = down_av_expr[[i]][sel[j],]
      subjectdata[[length(subjectdata)+1]] = colnames(down_av_expr[[i]])
      treat_m[[length(treat_m)+1]] = treat[[i]]
      signature_m[[length(signature_m)+1]] = signature[[i]]
      clus_m[[length(clus_m)+1]] = rownames(down_av_expr[[i]])[sel[j]]
    }
  }
}
mediator_m = rep(mediators, each = length(expr))
expr = rep(expr, length(mediators))
subjectdata = rep(subjectdata, length(mediators))
treat_m = rep(treat_m, length(mediators))
signature_m = rep(signature_m, length(mediators))
clus_m = rep(clus_m, length(mediators))


### parallel - SNOW 
strt<-Sys.time()
cl <- makeSOCKcluster(16)
registerDoSNOW(cl)
clusterEvalQ(cl, library(stringr))
clusterEvalQ(cl, library(mediation))
pb <- txtProgressBar(min = 1, max = length(expr), style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress=progress)
fullout_down = foreach(i=1:length(expr) , .options.snow = opts) %dopar% {  
  out = med_fun(treat_m[[i]], controls, mediator_m[[i]], subData, subjectdata[[i]], expr[[i]])
  out
}
close(pb)
parallel::stopCluster(cl)
print(Sys.time()-strt)


### Results
allresdown = list()
allresdown[[1]] = treat_m
allresdown[[2]] = signature_m
allresdown[[3]] = mediator_m
allresdown[[4]] = clus_m
allresdown[[5]] = fullout_down

### Save results
saveRDS(allresup, str_c(here("Res/ThirdSubmission/Mediation"),"/Up_Without_1KI_parallel.rds"))
saveRDS(allresdown, str_c(here("Res/ThirdSubmission/Mediation"),"/Down_Without_1KI_parallel.rds"))
