rm(list=ls(all=TRUE))
dev.off()
library(here)
library(mediation)
library(stringr)
library(purrr)
library(foreach)
library(parallel)
library(doParallel)

## Read clustering results
res = readRDS(str_c(here("Res/ThirdSubmission/Clustering/"), "1KI.rds"))

## Load expression data
dat = readRDS(str_c(here("data/"), "Filtered_ExpressionSet_Clustering_SelSubData.rds"))

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
  "insurance_lack"
)

subData[,which(colnames(subData) %in% treatment | colnames(subData) %in% mediators)] <- sapply(subData[,which(colnames(subData) %in% treatment | colnames(subData) %in% mediators)], as.numeric)


## Mediation function and extraction
med_fun = function(treatment, controls, mediator, subjectData, expressionData) {
  ## Data
  data = data.frame(Y = expressionData, subjectData)
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
treat = res[[7]]
signature = res[[8]]


### Mediation results
strt<-Sys.time()
cl <- parallel::makeCluster(24)
doParallel::registerDoParallel(cl)
clusterEvalQ(cl, library(stringr))
clusterEvalQ(cl, library(mediation))
clusterExport(cl=cl, c('up_av_expr', 'subData','treat','controls','mediators'))
fullout_up = foreach(i=1:length(treat) ,.packages = "foreach") %dopar% {  
  foreach(j=1:length(mediators),.packages = "foreach") %dopar% {
    clusexpr = up_av_expr[[i]]
    subjectdata = subData[colnames(clusexpr),]
    foreach(k=1:nrow(clusexpr),.combine ='cbind',.packages = "foreach") %dopar% {
      out = med_fun(treat[[i]], controls, mediators[[j]], subjectdata, clusexpr[k,])
      out
    }
  }
}
parallel::stopCluster(cl)
print(Sys.time()-strt)



### For downregulated average expression
down_av_expr = res[[4]] ## 55 different combo of treatment - signature
treat = res[[7]]
signature = res[[8]]


### Mediation results
strt<-Sys.time()
cl <- parallel::makeCluster(24)
doParallel::registerDoParallel(cl)
clusterEvalQ(cl, library(stringr))
clusterEvalQ(cl, library(mediation))
clusterExport(cl=cl, c('down_av_expr', 'subData','treat','controls','mediators'))
fullout_down = foreach(i=1:length(treat) ,.packages = "foreach") %dopar% {  
  foreach(j=1:length(mediators),.packages = "foreach") %dopar% {
    clusexpr = down_av_expr[[i]]
    subjectdata = subData[colnames(clusexpr),]
    foreach(k=1:nrow(clusexpr),.combine ='cbind',.packages = "foreach") %dopar% {
      out = med_fun(treat[[i]], controls, mediators[[j]], subjectdata, clusexpr[k,])
      out
    }
  }
}
parallel::stopCluster(cl)
print(Sys.time()-strt)


### Save results
saveRDS(fullout_up, str_c(here("Res/ThirdSubmission/Mediation"), "Up_1KI_parallel.rds"))
saveRDS(fullout_down, str_c(here("Res/ThirdSubmission/Mediation"), "Down_1KI_parallel.rds"))

