rm(list = ls())
#
library(plyr);library(data.table)
options(stringsAsFactors = F)
setwd("d:/data/NAFLD/")

metabolite = "lipid_bile"
disease = "obese_nash"

#
anno = NULL
for(i in 1:20){
  dir = paste0("./",metabolite, "_", disease, "/piumet_output (", i, ")/")
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
storage = paste0(metabolite, "_", disease, ".txt")
dat1 = read.table(storage, header = T)
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

storage = paste0("anno_", metabolite, "_", disease, ".txt")
write.table(anno, storage, sep = "\t", row.names = F)
