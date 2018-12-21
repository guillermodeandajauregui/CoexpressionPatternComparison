#compare chromosome ii vs all chromosomes ij

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


#compare cases

df_basal_ij = data.frame(i     = character(),
                         j     = character(),
                         wil.v = numeric(),
                         wil.p = numeric()
                         )

for(i in cromos){
  print(i)
  #extract ii matrix 
  chromo_i = sort(filter(.data = annot2, Chr == i)$symbol)
  xx = which(rownames(basal)%in%chromo_i)
  yy = which(rownames(basal)%in%chromo_i)
  matrix_ii = basal[xx,yy]
  
  #extract the ij matrix, calculate JSD
  for(j in cromos){
    print(j)
    chromo_j = sort(filter(.data = annot2, Chr == j)$symbol)
    yy = which(rownames(basal)%in%chromo_j)
    matrix_ij = basal[xx,yy]
    
    willy = wilcox.test(x = as.vector(matrix_ii), 
                        y = as.vector(matrix_ij) 
                        )
    
    wil.p  = willy$p.value
    wil.v  = willy$statistic
    r      = data.frame(i   = i,
                        j   = j,
                        wil.v = wil.v,
                        wil.p = wil.p
                        )
    
    df_basal_ij = rbind(df_basal_ij, r)
    
  }
}

head(df_basal_ij)
min(df_basal_ij$wil.v)
max(df_basal_ij$wil.v)

which(df_basal_ij$wil.v==min(df_basal_ij$wil.v))
df_basal_ij[480:482,]

p = ggplot(data = df_basal_ij, mapping = aes(i,j, fill=log(wil.v)))
p = p + geom_raster()
p

#KS


ks_basal_ij = data.frame(i     = character(),
                         j     = character(),
                         ks.d = numeric(),
                         ks.p = numeric()
)

for(i in cromos){
  print(i)
  #extract ii matrix 
  chromo_i = sort(filter(.data = annot2, Chr == i)$symbol)
  xx = which(rownames(basal)%in%chromo_i)
  yy = which(rownames(basal)%in%chromo_i)
  matrix_ii = basal[xx,yy]
  
  #extract the ij matrix, calculate JSD
  for(j in cromos){
    print(j)
    chromo_j = sort(filter(.data = annot2, Chr == j)$symbol)
    yy = which(rownames(basal)%in%chromo_j)
    matrix_ij = basal[xx,yy]
    
    kissy = ks.test(x = as.vector(matrix_ii), 
                        y = as.vector(matrix_ij) 
    )
    
    ks.d  = kissy$statistic
    ks.p  = kissy$p.value
    r      = data.frame(i   = i,
                        j   = j,
                        ks.d = ks.d,
                        ks.p = ks.p
    )
    
    ks_basal_ij = rbind(ks_basal_ij, r)
    
  }
}



ks_cntrl_ij = data.frame(i     = character(),
                         j     = character(),
                         ks.d = numeric(),
                         ks.p = numeric()
)

for(i in cromos){
  print(i)
  #extract ii matrix 
  chromo_i = sort(filter(.data = annot2, Chr == i)$symbol)
  xx = which(rownames(cntrl)%in%chromo_i)
  yy = which(rownames(cntrl)%in%chromo_i)
  matrix_ii = cntrl[xx,yy]
  
  #extract the ij matrix, calculate JSD
  for(j in cromos){
    print(j)
    chromo_j = sort(filter(.data = annot2, Chr == j)$symbol)
    yy = which(rownames(cntrl)%in%chromo_j)
    matrix_ij = cntrl[xx,yy]
    
    kissy = ks.test(x = as.vector(matrix_ii), 
                    y = as.vector(matrix_ij) 
    )
    
    ks.d  = kissy$statistic
    ks.p  = kissy$p.value
    r      = data.frame(i   = i,
                        j   = j,
                        ks.d = ks.d,
                        ks.p = ks.p
    )
    
    ks_cntrl_ij = rbind(ks_cntrl_ij, r)
    
  }
}


max(ks_cntrl_ij$ks.d)
min(ks_cntrl_ij$ks.d)

max(ks_basal_ij$ks.d)
min(ks_basal_ij$ks.d)


p1 = ggplot(data = ks_basal_ij, mapping = aes(i,j, fill=-log(ks.d)))
#p = ggplot(data = ks_basal_ij, mapping = aes(i,j, fill=(ks.d)))
p1 = p1 + geom_raster()
#p1 = p1 + scale_fill_gradient(low = "white", high = "steelblue", limits = c(0,7))
#p1 = p1 + scale_fill_distiller(palette = "Purples", direction = 1, limits = c(0,10))
#p1 = p1 + scale_fill_distiller(palette = "Spectral", direction = -1, limits = c(0,10))
#p1 = p1 + scale_fill_distiller(palette = "Spectral", direction = -1, limits = c(1,7))
p1 = p1 + scale_fill_distiller(palette = "RdYlBu", direction = -1, limits = c(1,7))
p1 = p1 + ggtitle(label = "Basal Breast Cancer")
p1 = p1 + ylab(label = "chromosome l") + xlab(label = "chromosome k")
p1


p2 = ggplot(data = ks_cntrl_ij, mapping = aes(i,j, fill=-log(ks.d)))
p2 = p2 + geom_raster()
#p2 = p2 + scale_fill_gradient(low = "white", high = "steelblue", limits = c(0,7))
#p2 = p2 + scale_fill_distiller(palette = "Purples", direction = 1, limits = c(0,10))
#p2 = p2 + scale_fill_distiller(palette = "Spectral", direction = -1, limits = c(0,10))
#p2 = p2 + scale_fill_distiller(palette = "Spectral", direction = -1, limits = c(1,7))
p2 = p2 + scale_fill_distiller(palette = "RdYlBu", direction = -1, limits = c(1,7))
#p2 = p2 + ggtitle(label = "control")
p2 = p2 + ggtitle(label = "Healthy Breast Tissue")
p2 = p2 + ylab(label = "chromosome l") + xlab(label = "chromosome k")
p2

pdf(file = "results/Figure03.pdf", 
    width = 15,
    height = 7 
    #paper = "US"
    
    )
cowplot::plot_grid(p1, p2, align = "v", labels = "AUTO")
dev.off()

filter(ks_cntrl_ij, i=="20")

arrange(ks_basal_ij, ks.d)[24:30,]
arrange(ks_cntrl_ij, ks.d)[24:30,]
