rm(list = ls())
library(data.table);library(plyr);library(org.Hs.eg.db)
setwd("d:/data/NAFLD/")

#
metabolite = "lipid"
disease = "obese_nash"
#Check!!!
deg_nash = read.table("SNU_DEG_NASH.txt")

#
storage = paste0("edge_robust_10_", metabolite, "_", disease, ".rdata") 
load(storage)
edge = edge_res[edge_res$num >= 10, ]
edge = unique(edge[, c("id1", "id2", "num")])

#
features = unlist(edge[, c(1, 2)])
features = unique(features)

#
storage = paste0("anno_", metabolite, "_", disease, ".txt")
anno = read.table(storage, header = T, stringsAsFactors = F)
anno_meta = anno[anno$type != "Protein", ]

#
idx1 = which(anno_meta$snu == "")
anno_meta$snu[idx1] = anno_meta$name[idx1]

#
idx1 = match(edge$id1, anno_meta$name)
idx_edge = which(!is.na(idx1))
idx_anno = idx1[!is.na(idx1)]
edge2 = edge
edge[idx_edge, "id1"] = anno_meta$snu[idx_anno]

#
idx1 = match(edge$id2, anno_meta$name)
idx_edge = which(!is.na(idx1))
idx_anno = idx1[!is.na(idx1)]
edge[idx_edge, "id2"] = anno_meta$snu[idx_anno]
edge = unique(edge)
idx1 = which(edge$id1 == "")
idx2 = which(edge$id2 == "")
idx3 = union(idx1, idx2)
if(length(idx3) > 0){
  edge = edge[-idx3, ]
}
storage = paste0(metabolite, "_", disease, "_network.rdata")
save(file = storage, edge)

#
edge1 = edge[, c(1, 2)]
edge2 = edge[, c(2, 1)]
colnames(edge1) = colnames(edge2) = c("id1", "id2")
edge = rbind(edge1, edge2)
edge = unique(edge)
edge_non_obese_nash = data.frame(table(edge$id1))
edge_non_obese_nash = as.character(edge_non_obese_nash$Var1[edge_non_obese_nash$Freq >= 3])
gene_all = anno$name[anno$type == "Protein"]
hub_gene = intersect(edge_non_obese_nash, gene_all)

#
HUB = select(org.Hs.eg.db, 
             keys = hub_gene, 
             keytype = "SYMBOL", columns = "GENENAME")

#
idx1 = match(HUB$SYMBOL, rownames(deg_nash))
idx_hub = which(!is.na(idx1))
idx_deg = idx1[!is.na(idx1)]

#
HUB$logFC = ""
HUB$logFC[idx_hub] = deg_nash$logFC[idx_deg]
HUB = HUB$SYMBOL[HUB$logFC != ""]
storage = paste0(metabolite, "_", disease, "_DEG_HUB.rdata")
save(file = storage, HUB)
