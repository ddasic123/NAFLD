rm(list = ls());options(stringsAsFactors = F)
library(data.table);library(plyr)
setwd("e:/data/NAFLD/")

#
edge_list = list()
for(j in 1:20){
  dir = paste0("./lipid_non_obese_nafl/piumet_output (", j, ")/")
  #
  file1 = "result_node_frequency_w10.0_b2.0_mu0.0005_R3.txt"
  file2 = "result_edge_frequency_w10.0_b2.0_mu0.0005_R3.txt"
  file3 = "result_union_net_w10.0_b2.0_mu0.0005_R3.txt"
  file1 = paste0(dir, file1)
  file2 = paste0(dir, file2)
  file3 = paste0(dir, file3)
  #node_info
  res1 = fread(file1)
  #edge_score
  res2 = fread(file2)
  #only_edge
  res3 = fread(file3, header = F)
  #
  res1 = data.frame(res1)
  res2 = data.frame(res2)
  res3 = data.frame(res3)
  
  #
  res3 = cbind(res3, res2)
  colnames(res3) = c("id1", "id2", "edge", "robust", "score")
  edge = res3$edge
  edge_list[[j]] = edge
}

#
edge_all = NULL
for(i in 1:length(edge_list)){
  edge_all = c(edge_all, edge_list[[i]])
}
edge_all = unique(edge_all)

#
edge_res = data.frame(edge = edge_all)
edge_res$num = 0
for(i in 1:length(edge_list)){
  edge.t = edge_list[[i]]
  idx1 = match(edge.t, edge_res$edge)
  edge_res[idx1, "num"] = edge_res[idx1, "num"] + 1
}
pdf("./lipid_non_obese_nafl/lipid_non_obese_nalf.pdf")
hist(edge_res$num, 20, col = "grey90",
     xlab = "# of replicated results [Gene-lipd pair] among 20 iterations")
abline(v = 10, col = "red", lwd = 2)
dev.off()

#
edge = edge_res$edge
edge = strsplit(edge, split = "pp")
for(i in 1:length(edge)){
  edge.t = edge[[i]]
  n1 = gsub(" \\(", "", edge.t[1])
  n2 = gsub("\\) ", "", edge.t[2])
  edge_res[i, "id1"] = n1
  edge_res[i, "id2"] = n2
}

#
write.table(edge_res, "edge_res_all.txt", sep = "\t", row.names = F)

#
edge_robust = edge_res$edge[edge_res$num >= 10]
save(file = "edge_robust_10_lipid_non_obese_nafl.rdata", edge_robust, edge_res)
