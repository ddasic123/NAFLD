rm(list = ls())
#################
library(sva)
library(data.table)
library(plyr)
library(limma)
setwd("d:/data/NAFLD/")

#
dat = fread("lipid_bile.txt")
dat = data.frame(dat)
rownames(dat) = dat$V1
dat = dat[, -1]
dat = log2(dat+1)
boxplot(temp)

#
dat = normalizeQuantiles(t(dat))
a = rep(1, 11)
b = rep(2, 239-11)
batch = c(a, b)
batch = as.factor(batch)
#batch_removal
dat = ComBat(dat = dat, batch = batch, par.prior=TRUE, prior.plots=FALSE)


#
label = read.table("lipid_bile_label.txt", header = T, stringsAsFactors = F)
label_num = table(label$label)


#batch
a = rep(1, 55)
b = rep(2, 92)
c = rep(3, 76)
batch = c(a, b, c)
batch = as.factor(batch)
#batch_removal
dat = t(dat)
combat_edata = ComBat(dat = dat, batch = batch, par.prior=TRUE, prior.plots=FALSE)
write.table(combat_edata, "lipid_bile_batch.txt", sep = '\t')
