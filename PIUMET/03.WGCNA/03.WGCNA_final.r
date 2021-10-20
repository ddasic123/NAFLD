rm(list=ls()); options(stringsAsFactors = F)

library(data.table);library(WGCNA);library(ggplot2);library(gridExtra);library(flashClust)
setwd("d:/data/NAFLD/")

#DS=3;MMS=5;DCOR=0.1;PAM=FALSE
#DS=1;MMS=5;DCOR=0.1;PAM=TRUE
#DS=4;MMS=5;DCOR=0.2;PAM=FALSE
#DS=3;MMS=5;DCOR=0.1;PAM=FALSE
DS=4;MMS=15;DCOR=0.1;PAM=TRUE


load("Coexpression_metabolite_211020.RData")
#final
tree = cutreeHybrid(
  dendro = geneTreeA1, 
  minClusterSize = MMS, ############
  pamStage = PAM, 
  cutHeight = 0.999, 
  deepSplit = DS, ################
  distM=dissTOMA1)

#
merged = mergeCloseModules(exprData= t(dat), 
                           colors = tree$labels, 
                           cutHeight=DCOR   ############
)

#
colors = labels2colors(merged$colors)
metabolite = rownames(dat)

#
plotDendroAndColors(dendro = geneTreeA1,
                    colors = colors,
                    groupLabels = "MOD",
                    cex.colorLabels = 1,
                    addGuide=T,
                    dendroLabels=F)

MEs = moduleEigengenes(expr = t(dat), colors, softPower = 15)
kMEtable = signedKME(t(dat), MEs$eigengenes)
save(file = "Coexpression_metabolite_211020.RData", dat, dissTOMA1, geneTreeA1, colors, metabolite, MEs)

#
pdf("F.module_metabolite.pdf",width=8,height=4)
plotDendroAndColors(geneTreeA1,
                    #cbind(colors, beta),
                    colors,
                    #groupLabels = c("Module", "NAFLD"),
                    groupLabels = c("Module"),
                    cex.colorLabels = 1,
                    addGuide=T,
                    dendroLabels=F)
dev.off()

module = data.frame(module = colors, metabolite = metabolite)
bile = read.table("bile_cn_nafl_nash.txt")
idx1 = which(module$metabolite %in% rownames(bile))
#
module$type = "lipid"
module$type[idx1] = "bile"
write.table(module, "lipid_bile_module.txt", sep = "\t", row.names = F)
