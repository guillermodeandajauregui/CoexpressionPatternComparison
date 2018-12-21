#libraries
library(tidyverse)
library(gplots)
library(reshape2)
library(textmineR)
library(fpc)

#load data 
basal = readRDS(file = "data/basalCleanMatrix.RDS")
cntrl = readRDS(file = "data/cntrlCleanMatrix.RDS")
annot = readRDS(file = "data/gene_annotation.RDS")

cromos = c(1:22, "X") #define chromosomes, skip chromosome Y 
names(cromos) = cromos

#sort genes by chromosome, then alphabetically 
annot2 = arrange(annot, Chr)




ks_cases_cntrl_ij = data.frame(i     = character(),
                         j     = character(),
                         ks.d = numeric(),
                         ks.p = numeric()
)


for(i in cromos){
  for(j in cromos){
    print(i)
    print(j)
  #extract ij matrices
    
  chromo_i = sort(filter(.data = annot2, Chr == i)$symbol)
  chromo_j = sort(filter(.data = annot2, Chr == j)$symbol)
  
  xx = which(rownames(basal)%in%chromo_i)
  yy = which(rownames(basal)%in%chromo_j)
  cases_ij = basal[xx,yy]
  cntrl_ij = cntrl[xx,yy]
  print(cases_ij[1:4,1:4])
  print(cntrl_ij[1:4,1:4])
  #calculate KS divergence 
  
  kissy.d = ks.test(cntrl_ij, cases_ij)$statistic
  kissy.p = ks.test(cntrl_ij, cases_ij)$p.value
  
  #return
  regreso = data.frame(i = i,
                       j = j, 
                       ks.d = kissy.d,
                       ks.p = kissy.p
                       )
  
  
  ks_cases_cntrl_ij = rbind(ks_cases_cntrl_ij, regreso)
  }
}

#plotting
max(ks_cases_cntrl_ij$ks.d)

p = ggplot(data = ks_cases_cntrl_ij, mapping = aes(i,j, fill=-log(ks.d)))
#p = ggplot(data = ks_basal_ij, mapping = aes(i,j, fill=(ks.d)))
p = p + geom_raster()
#p = p + scale_fill_gradient(low = "white", high = "red", limits = c(0,7))
#p = p + scale_fill_gradient(low = "white", high = "pink", limits = c(0,3))
p = p + scale_fill_distiller(palette = "BuPu", direction = 1, limits = c(1,3))
p = p + ggtitle(label = "KS distance (Basal breast cancer  vs Healthy breast tissue)")
p = p + ylab(label = "Chromosome") + xlab(label = "Chromosome")
p
ggsave(filename = "results/figure02_alt.pdf", 
       height = 8.5, 
       width = 14, 
       units = "in")


arrange(ks_cases_cntrl_ij, ks.d)%>%head
ks_cases_cntrl_ij%>%filter(i == j)%>%arrange(ks.d)%>%tail
ks_cases_cntrl_ij%>%filter(i != j)%>%arrange(ks.d)%>%head
