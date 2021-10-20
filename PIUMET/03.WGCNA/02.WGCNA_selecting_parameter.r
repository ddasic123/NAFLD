rm(list=ls()); options(stringsAsFactors = F)

library(data.table);library(WGCNA);library(ggplot2);library(gridExtra);library(flashClust)
setwd("d:/data/NAFLD/")

#48#245
dat = read.table("lipid_bile_batch.txt")

#supp
powers = c(seq(1,10,by=1),seq(11,20,by=2))
sft = pickSoftThreshold(data= t(dat), 
                        networkType = "signed", 
                        verbose=5,
                        corFnc = "cor",
                        powerVector=powers)

#
pdf("FS.soft.power.pdf",width=12,height=6)
par(mfrow=c(1,2))
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], 
     xlab="Soft Thresh Power", ylab="Scale free R^2",type="n")
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], 
     labels = powers, cex = 0.7, col="red",  xlab="Soft Thresh Power", ylab="Scale free R^2")
abline(h=0.8, col="black")
plot(sft$fitIndices[,1], sft$fitIndices[,5], 
     xlab = "Soft threshold power", 
     ylab = "Mean connectivity", type = "n")
text(sft$fitIndices[,1], sft$fitIndices[,5], labels = powers, cex = 0.7, col="black")
dev.off()

#9
adjacencyA1 = adjacency(t(dat),
                         power=9,
                         type="signed", 
                         corFnc = "cor")
diag(adjacencyA1) = 0
dissTOMA1 = 1-TOMsimilarity(adjacencyA1, TOMType="signed")
geneTreeA1 = flashClust(as.dist(dissTOMA1), method="average")

#Iterate WGCNA parameters for robustness -- this takes a while
colors = vector(mode="list")
labels = vector(mode="list")
for(pam in c(TRUE))
  for (minModSize in c(10, 15, 20, 25)) {
    for (dthresh in c(0.1, 0.2)) {
      for(ds in c(0:4)) { 
        print(paste("DS=", ds, ",MMS=", minModSize, ",DCOR=",dthresh,",PAM=",pam,sep=""))
        
        tree = cutreeHybrid(dendro = geneTreeA1, 
                            minClusterSize= minModSize, 
                            pamStage = pam, 
                            cutHeight = 0.999, 
                            deepSplit = ds, 
                            distM = dissTOMA1)
        merged = mergeCloseModules(exprData= t(dat), colors = tree$labels, cutHeight=dthresh)
        colors = cbind(colors, labels2colors(merged$colors))
        
        labels = c(labels, paste("DS=", ds, ",MMS=", minModSize, ",DCOR=",dthresh,",PAM=",pam,sep=""))
      }
    }
  }

#
pdf("FS.module_diffParams.pdf",width=12,height=14)
plotDendroAndColors(geneTreeA1, colors, groupLabels = labels,
                    addGuide=T, dendroLabels=F, cex.colorLabels=0.4)
dev.off()

#DS=3,MMS=5,DCOR=0.1,PAM=FALSE
#DS=1,MMS=5,DCOR=0.1,PAM=TRUE
#DS=4,MMS=5,DCOR=0.2,PAM=FALSE
#DS=3,MMS=5,DCOR=0.1,PAM=FALSE
#DS=4,MMS=15,DCOR=0.1,PAM=TRUE
save(file="Coexpression_metabolite_211020.RData", dat, dissTOMA1, geneTreeA1, colors, labels)







