#libraries
library(tidyverse)
library(gplots)
library(reshape2)

#load data 
basal = readRDS(file = "data/basalCleanMatrix.RDS")
cntrl = readRDS(file = "data/cntrlCleanMatrix.RDS")
annot = readRDS(file = "data/gene_annotation.RDS")

#check similarity in gene order 
all(rownames(basal)==rownames(cntrl))
all(colnames(basal)==colnames(cntrl))

#check if they are annotated
all(rownames(basal)%in%annot$symbol)
all(rownames(cntrl)%in%annot$symbol)

all(colnames(basal)%in%annot$symbol)
all(colnames(cntrl)%in%annot$symbol)

##not all annotated, let's remove those that aren't in the annotation
basal = basal[-which(!(rownames(basal)%in%annot$symbol)),
             -which(!(colnames(basal)%in%annot$symbol)),
             ]

cntrl = cntrl[-which(!(rownames(cntrl)%in%annot$symbol)),
              -which(!(colnames(cntrl)%in%annot$symbol)),
              ]

#add lower matrix triangle

#basal[lower.tri(basal)] <- basal[upper.tri(basal)]
#cntrl[lower.tri(cntrl)] <- cntrl[upper.tri(cntrl)]

basal = as.matrix(Matrix::forceSymmetric(as.matrix(basal), uplo = "U"))
cntrl = as.matrix(Matrix::forceSymmetric(as.matrix(cntrl), uplo = "U"))

#check for NA

head(which(is.na(basal), arr.ind = TRUE))
head(which(is.na(cntrl), arr.ind = TRUE))


#convert NAs to zeros

basal<-as.matrix(basal)
basal[which(is.na(basal))]<-0

cntrl<-as.matrix(cntrl)
cntrl[which(is.na(cntrl))]<-0
which(is.na(cntrl))
#difference between basal and cntrl
vs_basal_cntrl = basal - cntrl

#sort genes by chromosome, then alphabetically 
annot2 = arrange(annot, Chr)

#extract genes from chromosome 1 
#head(annot2)
chromo_1 = sort(filter(.data = annot2, Chr == 1)$symbol)

#heatmap of genes from chromosome vs themselves 
vs_basal_cntrl_1 = vs_basal_cntrl[rownames(vs_basal_cntrl)%in%chromo_1,
                                  colnames(vs_basal_cntrl)%in%chromo_1]

# gplots::heatmap.2(as.matrix(vs_basal_cntrl_1), 
#                   trace = "none", 
#                   col = bluered, 
#                   key = TRUE, 
#                   dendrogram = "none", Rowv = NULL, Colv = NULL, 
#                   main = "diff basal healthy chr1 ")


m = melt(as.matrix(vs_basal_cntrl_1))
p = ggplot(m, aes(Var1,Var2, fill=value)) + 
  geom_raster() + 
  scale_fill_gradient2(low = "blue", high = "red") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 1)) +
  theme(axis.text.y = element_text(size = 1))
p  
#chromosome 1 vs all others 
xx = which(rownames(vs_basal_cntrl)%in%chromo_1)
yy = match(table = colnames(vs_basal_cntrl), x = annot2$symbol)
yy = yy[-which(is.na(yy))]

vs_basal_cntrl_1vALL = vs_basal_cntrl[xx,yy]

m = melt(as.matrix(vs_basal_cntrl_1vALL))
p = ggplot(m, aes(Var1,Var2, fill=value)) + 
  geom_raster() + 
  scale_fill_gradient2(low = "blue", high = "red") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 1)) +
  theme(axis.text.y = element_text(size = 1))
p

cromos = as.character(unique(annot2$Chr))

# 
# 
# for(i in cromos){
#   chromo_col = sort(filter(.data = annot2, Chr == i)$symbol)
#   xx = which(rownames(vs_basal_cntrl)%in%chromo_1)
#   yy = which(colnames(vs_basal_cntrl)%in%chromo_col)
#   
#   my_matrix = vs_basal_cntrl[xx,yy]
#   
#   print(i)
#   print("positives")
#   print(length(which(my_matrix>0)))
#   print("negatives")
#   print(length(which(my_matrix<0)))
#   print("mas positivos que negativos?")
#   print(length(which(my_matrix>0))>length(which(my_matrix<0)))
#   
#   m = melt(as.matrix(my_matrix))
#   p = ggplot(m, aes(Var1,Var2, fill=value)) + 
#     geom_raster() + 
#     scale_fill_gradient2(low = "blue", high = "red") +
#     theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 1)) +
#     theme(axis.text.y = element_text(size = 1))+
#     ggtitle(paste0("chromosome ", i))
# #  plot(p)
# }

# my_frame = data.frame(cromo_i = NA,
#                       cromo_j = NA,
#                       pluses  = NA,
#                       minus   = NA,
#                       zeros   = NA,
#                       mostly  = NA,
#                       sum     = NA)

my_frame = data.frame(cromo_i     = factor(),
                      cromo_j     = factor(),
                      pluses      = numeric(),
                      minus       = numeric(),
                      zeros       = numeric(),
                      mostly      = numeric(),
                      dif.pos.neg = numeric(),
                      sum         = numeric(),
                      avg         = numeric(),
                      avg.pos     = numeric(),
                      avg.neg     = numeric(),
                      median      = numeric(),
                      max         = numeric(),
                      min         = numeric()
                      )

for(i in cromos){
  chromo_row = sort(filter(.data = annot2, Chr == i)$symbol)
  for(j in cromos){
    chromo_col = sort(filter(.data = annot2, Chr == j)$symbol)
    xx = which(rownames(vs_basal_cntrl)%in%chromo_row)
    yy = which(colnames(vs_basal_cntrl)%in%chromo_col)
    
    my_matrix = vs_basal_cntrl[xx,yy]
    
    positivos = length(which(my_matrix>0))
    negativos = length(which(my_matrix<0))
    #print(positivos>negativos)
    dif.pos.neg = positivos-negativos
    mostly      = ifelse(positivos>negativos, "pos", "neg")
    suma        = sum(my_matrix)
    zeros       = length(which(my_matrix==0))
    avg         = mean(my_matrix)
    avg.pos     = mean(my_matrix[which(my_matrix>0)])
    avg.neg     = mean(my_matrix[which(my_matrix<0)])
    maxy        = max(my_matrix)
    miny        = min(my_matrix)
    med         = median(my_matrix)
    #r = c(i,j,positivos, negativos, zeros, NA, suma)
    r = data.frame(cromo_i     = i,
                   cromo_j     = j,
                   pluses      = positivos,
                   minus       = negativos,
                   zeros       = zeros, 
                   mostly      = mostly,
                   dif.pos.neg = dif.pos.neg,
                   sum         = suma, 
                   avg         = avg,
                   avg.pos     = avg.pos,
                   avg.neg     = avg.neg,
                   median      = med,
                   max         = maxy,
                   min         = miny
                   )
    my_frame = rbind(my_frame, r)
    
  }}

# my_frame<-my_frame[-1,]
# my_frame$pluses <- as.numeric(my_frame$pluses)
# my_frame$minus <- as.numeric(my_frame$minus)
# my_frame$sum <- as.numeric(my_frame$sum)
# my_frame$zeros <- as.numeric(my_frame$zeros)
# 
# my_frame$mostly<-ifelse(my_frame$minus<my_frame$pluses, "pos", "neg")

my_better_frame = my_frame[!(my_frame$cromo_i=="Y"),]
my_better_frame = my_better_frame[!(my_better_frame$cromo_j=="Y"),]

#change factor: abs(avg.pos)/abs(avg.neg) 
#Larger than 1 if avg of gains is larger than avg of losses.
#Less than 1 if avg of losses is larger than avg of gains

my_better_frame$ChangeFactor = abs(my_better_frame$avg.pos)/abs(my_better_frame$avg.neg)
my_better_frame[which(log(my_better_frame$ChangeFactor)>0),]

my_better_frame$signFactor = (my_better_frame$pluses - my_better_frame$minus)/(my_better_frame$pluses + my_better_frame$minus + my_better_frame$zeros)

my_better_frame$intra.inter = ifelse(my_better_frame$cromo_i==my_better_frame$cromo_j, "intra", "inter")


for(i in cromos){
  chromo_row = sort(filter(.data = annot2, Chr == i)$symbol)
  for(j in cromos){
    chromo_col = sort(filter(.data = annot2, Chr == j)$symbol)
    xx = which(rownames(vs_basal_cntrl)%in%chromo_row)
    yy = which(colnames(vs_basal_cntrl)%in%chromo_col)
    
    my_matrix = vs_basal_cntrl[xx,yy]
    my_reduced = my_matrix[,abs(colMeans(my_matrix))<abs(mean(my_matrix))]
    
    
    
    m = melt(as.matrix(my_reduced))
    p = ggplot(m, aes(Var1,Var2, fill=value)) + 
      geom_raster() + 
      scale_fill_gradient2(low = "blue", high = "red") +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 1)) +
      theme(axis.text.y = element_text(size = 1))+
      ggtitle(paste0("chromosome ", i, "versus ", j))
    plot(p)
    my_file = paste0("results/chromos",i,"_",j,".pdf")
    ggsave(filename = my_file)
  }}

xx = which(rownames(vs_basal_cntrl)%in%chromo_1)
yy = which(colnames(vs_basal_cntrl)%in%chromo_1)
my_test_matrix = vs_basal_cntrl[xx,yy]

max(my_test_matrix)
min(my_test_matrix)
mean(my_test_matrix)
colSums(my_test_matrix)
colMeans(my_test_matrix)

table(abs(colMeans(my_test_matrix))<abs(mean(my_test_matrix)))
ncol(my_test_matrix)


