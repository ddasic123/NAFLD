rm(list = ls())
library(data.table);library(plyr);library(igraph)
setwd("d:/data/NAFLD/")

#
load("diff_features_non_obese_nash.rdata")
dat1 = fread("lipid_non_obese_nafl.txt")
dat1 = data.frame(dat1)

#
idx1 = which(dat1$p < 0.2 & dat1$fc > 1)
idx2 = which(dat1$p < 0.2 & dat1$fc < 1)
lipid_up = dat1$SNU_ID[idx1]
lipid_do = dat1$SNU_ID[idx2]

save(file = "diff_features_lipid_non_obese_nafl.rdata", gene_up, gene_do, lipid_up, lipid_do)
