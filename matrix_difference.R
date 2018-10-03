#libraries
library(tidyverse)
library(gplots)

#load data 
basal = readRDS(file = "data/basalCleanMatrix.RDS")
cntrl = readRDS(file = "data/cntrlCleanMatrix.RDS")
annot = readRDS(file = "data/gene_annotation.RDS")

#check similarity in gene order 
all(rownames(basal)==rownames(cntrl))
all(colnames(basal)==colnames(cntrl))

#add lower matrix triangle

basal[lower.tri(basal)] <- basal[upper.tri(basal)]
cntrl[lower.tri(cntrl)] <- cntrl[upper.tri(cntrl)]

#difference between basal and cntrl
vs_basal_cntrl = basal - cntrl

#extract genes from chromosome 1 
head(annot)
chromo_1 = sort(filter(.data = annot, Chr == 1)$symbol)

