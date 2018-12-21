#Comparison of density functions

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



basal_densities = 
lapply(X = cromos, FUN = function(i){
  lapply(X = cromos, FUN = function(j){
    print(paste0(i, ";", j))
    
    chromo_i = sort(filter(.data = annot2, Chr == i)$symbol)
    xx = which(rownames(basal)%in%chromo_i)
    
    chromo_j = sort(filter(.data = annot2, Chr == j)$symbol)
    yy = which(rownames(basal)%in%chromo_j)
    
    matrix_ij = basal[xx,yy]
    
    #make diagonal NA if i == j 
    if(i == j){
      diag(matrix_ij)<-NA
    }
    
    #make density 
    resultado = density(matrix_ij, na.rm = TRUE)
    plot(resultado)
    return(resultado)
  })
})

#result: a list of 23 lists, each with 23 probability density functions 

#plotting
lapply(X = basal_densities, FUN = function(j){
  lapply(X = j, FUN = function(i){
  x2 = data.frame(x = i$x, y = i$y)
  ggplot(x2, aes(x,y)) + geom_line(colour = "red")
  })
})

#

testing_ii = basal_densities$`3`$`3`
testing_ij = basal_densities$`3`$`5`
testing_ik = basal_densities$`3`$`7`

textmineR::CalcHellingerDist(testing_ij$y, 
                              testing_ik$y)




density_comparisons_basal = data.frame(cromo_i = character(),
                                 cromo_j = character(),
                                 hellinger_D = numeric(),
                                 JSD_D       = numeric()
                                 )


for(i in seq_along(cromos)){
  cromo_i = cromos[i]
  pdf_ii = basal_densities[[i]][[i]]$y
  
  for(j in seq_along(cromos)){
    cromo_j = cromos[j]
    pdf_ij = basal_densities[[i]][[j]]$y
    
    my_result        = data.frame(cromo_i = cromo_i,
                         cromo_j = cromo_j,
                         hellinger_D = textmineR::CalcHellingerDist(pdf_ii, 
                                                                    pdf_ij),
                         JSD_D       = textmineR::CalcJSDivergence(pdf_ii, 
                                                                    pdf_ij)
                          )
    
    density_comparisons_basal = rbind(density_comparisons_basal,
                                      my_result)
    
  }
}

p = ggplot(data = density_comparisons_basal, 
           mapping = aes(cromo_i,cromo_j, fill=hellinger_D)
           )


#p = ggplot(data = ks_basal_ij, mapping = aes(i,j, fill=(ks.d)))
p = p + geom_raster()
p = p + scale_fill_gradient(low = "white", high = "steelblue", limits = c(0,1))
p = p + ggtitle(label = "basal")
p

p = ggplot(data = density_comparisons_basal, 
           mapping = aes(cromo_i,cromo_j, fill=JSD_D)
)


#p = ggplot(data = ks_basal_ij, mapping = aes(i,j, fill=(ks.d)))
p = p + geom_raster()
p = p + scale_fill_gradient(low = "white", high = "steelblue", limits = c(0,1))
p = p + ggtitle(label = "basal")
p

head(density_comparisons_basal)

basal_densities[[1]][[2]]
basal_densities[[2]][[1]]

density_comparisons_basal%>%
  filter(cromo_i == 1 & cromo_j == 2)
density_comparisons_basal%>%
  filter(cromo_i == 2 & cromo_j == 1)

########################################
basal_ijkl =              data.frame(i = character(),
                             j = character(),
                             k = character(),
                             l = character(),
                             hellinger_D = numeric(),
                             JSD_D       = numeric()
                  )

for(i in seq_along(cromos)){
  for(j in seq_along(cromos)){
    pdf_ij = pdf_ij = basal_densities[[i]][[j]]$y
    for(k in seq_along(cromos)){
      for(l in seq_along(cromos)){
        pdf_kl = basal_densities[[k]][[l]]$y
        
        my_result = data.frame(i = cromos[i],
                               j = cromos[j],
                               k = cromos[k],
                               l = cromos[l],
                               hellinger_D = textmineR::CalcHellingerDist(pdf_ij, 
                                                                          pdf_kl),
                               JSD_D       = textmineR::CalcJSDivergence(pdf_ij, 
                                                                         pdf_kl)
        )
        
        basal_ijkl = rbind(basal_ijkl, my_result)
      }
    }
  }
}

saveRDS(basal_ijkl, "results/basal_ijkl.RDS")
basal_ijkl%>%arrange(desc(hellinger_D))%>%head(10)


basal_ijkl%>%
  filter(i==j & j==k)%>%
  mutate("i:j" = paste0(i,":",j),
         "k:l" = paste0(k,":",l)
         )%>%
  select(i, hellinger_D, l)%>%
  filter(i!=l)%>%
  mutate(inv_D = 1/hellinger_D)%>%
  write_delim(path = "results/basal_cromoDistance.txt")

