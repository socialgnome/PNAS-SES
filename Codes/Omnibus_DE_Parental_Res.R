## Setup env
library(here)
dev.off()
source(here("R/ThirdSubmission/Final/startup.R"))

## Signatures
signatures = readRDS("~/Projects/PNAS/data/Signatures.rds")

## Define treatment
treatment = c("ses_composite_pp1","SEI_max_p_w12","edu_p","income_pp1_log")

## With 1K - Analysis
res = data.frame(Treatment = rep(treatment, each = length(signatures)), 
                 Signature = rep(names(signatures), length(treatment)),
                 P = 0)

p = c()
t = c()
for (i in 1:length(treatment)) {
  
  for (j in 1:length(signatures)) {
    
    df = readRDS(str_c(here("Res/ThirdSubmission/Omnibus/Parental/"), names(signatures)[[j]], "_", treatment[[i]],".rds"))
    
    genes = signatures[[j]]
    sel = df[df$gene %in% genes,]
    p = c(p, min(sel$adj.P.Val))
    t = c(t, sel$t[which(sel$adj.P.Val==min(sel$adj.P.Val))[1]])
  }
  
}
res$P = p
res$t = t

res = res[res$P<0.05,]

## DeNovo analysis
disgenes = unique(unlist(signatures))
df = readRDS(str_c(here("Res/ThirdSubmission/Omnibus/Parental/"), names(signatures)[[3]], "_", treatment[[1]],".rds"))
length(df$gene[df$adj.P.Val<0.05])
df = df[-which(df$gene %in% disgenes),]
length(df$gene[df$adj.P.Val<0.05])
