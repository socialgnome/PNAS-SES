rm(list=ls(all=TRUE))
dev.off()
library(here)
library(multimediate)
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
  data = data %>% dplyr::select(Y, all_of(treatment),all_of(controls), all_of(mediator))
  data = data[complete.cases(data),]
  ## M.reg
  M1 = as.formula(str_c(mediator[[1]]," ~ ",str_c(c(treatment,controls), collapse = " + ")))
  fit_M1 = lm(M1, data)
  M2 = as.formula(str_c(mediator[[2]]," ~ ",str_c(c(treatment,controls), collapse = " + ")))
  fit_M2 = lm(M2, data)
  M3 = as.formula(str_c(mediator[[3]]," ~ ",str_c(c(treatment,controls), collapse = " + ")))
  fit_M3 = lm(M3, data)
  M4 = as.formula(str_c(mediator[[4]]," ~ ",str_c(c(treatment,controls), collapse = " + ")))
  fit_M4 = lm(M4, data)
  M5 = as.formula(str_c(mediator[[5]]," ~ ",str_c(c(treatment,controls), collapse = " + ")))
  fit_M5 = lm(M5, data)
  ## Y.reg
  mod_formula = as.formula(str_c("Y ~ ",str_c(c(treatment, mediator, controls), collapse = " + ")))
  fit_y = lm(mod_formula,data)
  ## Mediate
  m = multimediate(lmodel.m = list(fit_M1,fit_M2,fit_M3, fit_M4, fit_M5),
                   model.y = fit_y,
                   treat = treatment, 
                   treat.value = 1,
                   control.value = 0,
                   J = 1000,
                   conf.level = 0.95)
  ## Extract results
  allres = list(m) %>% map(~.x %>% summary(opt = "avg"))
  allres = allres[[1]]
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



### parallel - SNOW 
strt<-Sys.time()
cl <- makeSOCKcluster(16)
registerDoSNOW(cl)
clusterEvalQ(cl, library(stringr))
clusterEvalQ(cl, library(purrr))
clusterEvalQ(cl, library(multimediate))
pb <- txtProgressBar(min = 1, max = length(expr), style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress=progress)
fullout_up = foreach(i=1:length(expr) , .options.snow = opts) %dopar% {  
  out = med_fun(treat_m[[i]], controls, mediators, subData, subjectdata[[i]], expr[[i]])
  out
}
close(pb)
parallel::stopCluster(cl)
print(Sys.time()-strt)

### Results
allresup = list()
allresup[[1]] = treat_m
allresup[[2]] = signature_m
allresup[[3]] = mediators
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


### parallel - SNOW 
strt<-Sys.time()
cl <- makeSOCKcluster(16)
registerDoSNOW(cl)
clusterEvalQ(cl, library(stringr))
clusterEvalQ(cl, library(purrr))
clusterEvalQ(cl, library(multimediate))
pb <- txtProgressBar(min = 1, max = length(expr), style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress=progress)
fullout_down = foreach(i=1:length(expr) , .options.snow = opts) %dopar% {  
  out = med_fun(treat_m[[i]], controls, mediators, subData, subjectdata[[i]], expr[[i]])
  out
}
close(pb)
parallel::stopCluster(cl)
print(Sys.time()-strt)


### Results
allresdown = list()
allresdown[[1]] = treat_m
allresdown[[2]] = signature_m
allresdown[[3]] = mediators
allresdown[[4]] = clus_m
allresdown[[5]] = fullout_down

### Save results
saveRDS(allresup, str_c(here("Res/ThirdSubmission/Multimediation"),"/Up_1KI_parallel.rds"))
saveRDS(allresdown, str_c(here("Res/ThirdSubmission/Multimediation"),"/Down_1KI_parallel.rds"))



######################
######################
######################

rm(list=ls(all=TRUE))
dev.off()
library(here)
library(multimediate)
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
  data = data %>% dplyr::select(Y, all_of(treatment),all_of(controls), all_of(mediator))
  data = data[complete.cases(data),]
  ## M.reg
  M1 = as.formula(str_c(mediator[[1]]," ~ ",str_c(c(treatment,controls), collapse = " + ")))
  fit_M1 = lm(M1, data)
  M2 = as.formula(str_c(mediator[[2]]," ~ ",str_c(c(treatment,controls), collapse = " + ")))
  fit_M2 = lm(M2, data)
  M3 = as.formula(str_c(mediator[[3]]," ~ ",str_c(c(treatment,controls), collapse = " + ")))
  fit_M3 = lm(M3, data)
  M4 = as.formula(str_c(mediator[[4]]," ~ ",str_c(c(treatment,controls), collapse = " + ")))
  fit_M4 = lm(M4, data)
  M5 = as.formula(str_c(mediator[[5]]," ~ ",str_c(c(treatment,controls), collapse = " + ")))
  fit_M5 = lm(M5, data)
  ## Y.reg
  mod_formula = as.formula(str_c("Y ~ ",str_c(c(treatment, mediator, controls), collapse = " + ")))
  fit_y = lm(mod_formula,data)
  ## Mediate
  m = multimediate(lmodel.m = list(fit_M1,fit_M2,fit_M3, fit_M4, fit_M5),
                   model.y = fit_y,
                   treat = treatment, 
                   treat.value = 1,
                   control.value = 0,
                   J = 1000,
                   conf.level = 0.95)
  ## Extract results
  allres = list(m) %>% map(~.x %>% summary(opt = "avg"))
  allres = allres[[1]]
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



### parallel - SNOW 
strt<-Sys.time()
cl <- makeSOCKcluster(16)
registerDoSNOW(cl)
clusterEvalQ(cl, library(stringr))
clusterEvalQ(cl, library(purrr))
clusterEvalQ(cl, library(multimediate))
pb <- txtProgressBar(min = 1, max = length(expr), style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress=progress)
fullout_up = foreach(i=1:length(expr) , .options.snow = opts) %dopar% {  
  out = med_fun(treat_m[[i]], controls, mediators, subData, subjectdata[[i]], expr[[i]])
  out
}
close(pb)
parallel::stopCluster(cl)
print(Sys.time()-strt)

### Results
allresup = list()
allresup[[1]] = treat_m
allresup[[2]] = signature_m
allresup[[3]] = mediators
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



### parallel - SNOW 
strt<-Sys.time()
cl <- makeSOCKcluster(16)
registerDoSNOW(cl)
clusterEvalQ(cl, library(stringr))
clusterEvalQ(cl, library(purrr))
clusterEvalQ(cl, library(multimediate))
pb <- txtProgressBar(min = 1, max = length(expr), style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress=progress)
fullout_down = foreach(i=1:length(expr) , .options.snow = opts) %dopar% {  
  out = med_fun(treat_m[[i]], controls, mediators, subData, subjectdata[[i]], expr[[i]])
  out
}
close(pb)
parallel::stopCluster(cl)
print(Sys.time()-strt)


### Results
allresdown = list()
allresdown[[1]] = treat_m
allresdown[[2]] = signature_m
allresdown[[3]] = mediators
allresdown[[4]] = clus_m
allresdown[[5]] = fullout_down

### Save results
saveRDS(allresup, str_c(here("Res/ThirdSubmission/Multimediation"),"/Up_Without_1KI_parallel.rds"))
saveRDS(allresdown, str_c(here("Res/ThirdSubmission/Multimediation"),"/Down_Without_1KI_parallel.rds"))
