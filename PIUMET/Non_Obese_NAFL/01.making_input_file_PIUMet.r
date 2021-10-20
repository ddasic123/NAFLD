rm(list = ls())
options(stringsAsFactors = F)
setwd("d:/data/NAFLD/")
library(data.table);library(matrixStats);library(limma);library(plyr)

#
dat1 = fread("apt16066-sup-0002-tables2.txt")
dat1 = data.frame(dat1)
#dat1$
dat.p = dat1[, c("Lipids", "non_obese_nl_nafld")]
colnames(dat.p) = c("SNU_ID", "p")

#
dat1 = fread("apt16066-sup-0002-tables_FC.csv")
dat1 = data.frame(dat1)
#dat1$SNU_ID
dat.fc = dat1[, c("SNU_ID", "non_obese_nafld")]
colnames(dat.fc) = c("SNU_ID", "fc")

#
dat1 = join(dat.fc, dat.p, by = "SNU_ID")
dat1 = dat1[dat1$p < 0.2, ]

#
annotation = read.table("annotation_final_210721.txt", sep = "\t", header = T)
annotation$paper = gsub("TG", "TAG", annotation$SNU)
annotation$paper = gsub("DG", "DAG", annotation$paper)
annotation$paper = gsub("Lyso", "L", annotation$paper)
idx1 = match(dat1$SNU_ID, annotation$paper)
dat1$HMDB = annotation$HMDB[idx1]
dat1 = dat1[!is.na(dat1$HMDB), ]
colnames(dat1)[4] = "HMDB_ID"
write.table(dat1, "lipid_non_obese_nafl.txt", sep = "\t", row.names = F)

#
df_piumet = NULL
df = dat1
for(i in 1:nrow(df)){
  df.t = df[i, ]
  hmdb = unlist(strsplit(df.t$HMDB_ID, split = ";"))
  prize = -log(df.t$p)
  df.t = data.frame(id = hmdb, prize = prize, type = "Metabolite")
  df_piumet = rbind(df_piumet, df.t)
}
#11#32
write.table(df_piumet, "piumet_lip_nobefl.txt", sep = "\t", row.names = F, col.names = F, quote = F)

