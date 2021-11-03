rm(list = ls())
library(data.table);library(plyr);library(igraph)
setwd("d:/data/NAFLD/")

#
metabolite = "lipid_bile"
disease = "non_obese_nafl"

#
storage = paste0(metabolite, "_", disease, "_DEG_HUB.rdata")
load(storage)
candidate.gene = HUB

#
storage = paste0("diff_features_", metabolite, "_",disease, ".rdata")
load(storage)
f_up = union(gene_up, lipid_up)
f_down = union(gene_do, lipid_do)

#step4
storage = paste0(metabolite, "_", disease, "_network.rdata")
load(storage)
res = edge

#
i = 3
for(i in 1:length(candidate.gene)){
  gene = candidate.gene[i]
  idx1 = which(res$id1 == gene | res$id2 == gene)
  res.t = res[idx1, ]
  idx1 = which(res.t$id2 == gene)
  res.t[idx1, ] = res.t[idx1, c(2, 1, 3)]
  df = as.matrix(res.t[, c(1, 2)])
  df = unique(df)
  g = graph.edgelist(df, directed = T)
  g = simplify(g)
  
  #node.shape
  V(g)$shape = "circle"
  
  #nodel.color
  V(g)$color = adjustcolor("black", alpha.f = .1)
  
  #node.size
  V(g)$size = 5
  
  #lable.size
  V(g)$label.cex = 0.7
  idx1 = which(V(g)$name %in% f_up)
  V(g)[idx1]$label.cex = 1.2
  idx1 = which(V(g)$name %in% f_down)
  V(g)[idx1]$label.cex = 1.2
  idx1 = which(V(g)$name %in% gene)
  V(g)[idx1]$label.cex = 1.5
  
  #label.color
  V(g)$label.color = "black"
  idx1 = which(V(g)$name %in% f_up)
  V(g)[idx1]$label.color = "red"
  idx1 = which(V(g)$name %in% f_down)
  V(g)[idx1]$label.color = "blue"
  
  #
  edge = get.edgelist(g)
  pair = paste0(edge[, 1], "_", edge[, 2])
  E(g)$width = 2
  E(g)$label = res.t$num
  
  #edge_color
  E(g)$color = adjustcolor("grey", alpha.f = 0.2)
  # idx1 = which(pair %in% sig.pair.up)
  # E(g)[idx1]$color = adjustcolor("red", alpha.f = 0.2)
  # idx1 = which(pair %in% sig.pair.down)
  # E(g)[idx1]$color = adjustcolor("blue", alpha.f = 0.2)
  
  if(nrow(edge) > 30){
    size = 12
  } else  if(nrow(edge) > 10){
    size = 9
  } else{
    size = 6
  }
  
  #
  set.seed(2)
  storage = paste0("./", metabolite, "_", disease, "_gene/", gene,"_", i, "_", metabolite, "_", disease, "_nafl.pdf")
  pdf(storage, width = size, height = size)
  plot(g,
       vertex.frame.color = NA,
       vertex.label.degree = 0,
       edge.arrow.size = 0,
       edge.label = E(g)$label
  )
  dev.off()
  
}


#
df = as.matrix(res[, c(1, 2)])
df = unique(df)
g = graph.edgelist(df, directed = T)
g = simplify(g)

#node.shape
V(g)$shape = "circle"

#nodel.color
V(g)$color = adjustcolor("black", alpha.f = .1)

#node.size
V(g)$size = 5

#lable.size
V(g)$label.cex = 0.3
idx1 = which(V(g)$name %in% f_up)
V(g)[idx1]$label.cex = 1
idx1 = which(V(g)$name %in% f_down)
V(g)[idx1]$label.cex = 1
idx1 = which(V(g)$name %in% candidate.gene)
V(g)[idx1]$label.cex = 2

#label.color
V(g)$label.color = "black"
idx1 = which(V(g)$name %in% f_up)
V(g)[idx1]$label.color = "red"
idx1 = which(V(g)$name %in% f_down)
V(g)[idx1]$label.color = "blue"

#
edge = get.edgelist(g)
pair = paste0(edge[, 1], "_", edge[, 2])
E(g)$width = 2
idx1 = which(pair %in% sig.pair.up)
E(g)[idx1]$width = 3
idx1 = which(pair %in% sig.pair.down)
E(g)[idx1]$width = 3

#edge_color
E(g)$color = adjustcolor("grey", alpha.f = 0.2)
idx1 = which(pair %in% sig.pair.up)
E(g)[idx1]$color = adjustcolor("red", alpha.f = 0.2)
idx1 = which(pair %in% sig.pair.down)
E(g)[idx1]$color = adjustcolor("blue", alpha.f = 0.2)

if(nrow(edge) > 30){
  size = 10
} else  if(nrow(edge) > 10){
  size = 6
} else{
  size = 4
}

#
set.seed(2)
storage = paste0("network_non_obese_nash.pdf")
pdf(storage, width = 30, height = 30)
plot(g,
     vertex.frame.color = NA,
     vertex.label.degree = 0,
     edge.arrow.size = 0
)
dev.off()
