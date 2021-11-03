rm(list = ls())
library(data.table);library(plyr);library(igraph)
setwd("d:/data/NAFLD/")

#
metabolite = "lipid_bile"
disease = "non_obese_nafl"

#checked!!!
dat1 = read.table("SNU_DEG_nafld.txt")
gene_up = rownames(dat1)[dat1$logFC > 0 & dat1$adj.P.Val < 0.05]
gene_do = rownames(dat1)[dat1$logFC < 0 & dat1$adj.P.Val < 0.05]

#
storage = paste0(metabolite, "_", disease, ".txt")
dat1 = fread(storage)
dat1 = data.frame(dat1)

#
idx1 = which(dat1$p < 0.2 & dat1$fc > 1)
idx2 = which(dat1$p < 0.2 & dat1$fc < 1)
lipid_up = dat1$SNU_ID[idx1]
lipid_do = dat1$SNU_ID[idx2]

storage = paste0("diff_features_", metabolite, "_",disease, ".rdata")
save(file = storage, gene_up, gene_do, lipid_up, lipid_do)
