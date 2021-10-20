rm(list = ls())
#
library(plyr)
options(stringsAsFactors = F)
setwd("e:/data/NAFLD/")

#
anno = NULL
for(i in 1:20){
  dir = paste0("./lipid_non_obese_nafl/piumet_output (", i, ")/")
  #
  file1 = "result_node_frequency_w10.0_b2.0_mu0.0005_R3.txt"
  file1 = paste0(dir, file1)
  #node_info
  res1 = fread(file1)
  res1 = data.frame(res1)
  res1 = res1[, -2]
  colnames(res1) = c("name", "type", "hmdb")
  anno = rbind(anno, res1)
}  
anno = unique(anno)

#
dat1 = read.table("lipid_non_obese_nafl.txt", header = T)
anno2 = NULL
for(i in 1:nrow(dat1)){
  df.t = dat1[i, ]
  hmdb = unlist(strsplit(df.t$HMDB_ID, split = ";"))
  df.t = data.frame(hmdb = hmdb, snu = df.t$SNU_ID)
  anno2 = rbind(anno2, df.t)
}
anno2 = unique(anno2)
anno = join(anno, anno2, by = "hmdb")
anno[is.na(anno)] = ""

write.table(anno, "anno_lipid_non_obese_nafl.txt", sep = "\t", row.names = F)