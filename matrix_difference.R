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

#difference between basal and cntrl
vs_basal_cntrl = basal - cntrl

#sort genes by chromosome, then alphabetically 
annot2 = arrange(annot, Chr)

#generate dataframe of results

cromos = c(1:22, "X") #define chromosomes, skip chromosome Y 

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
    dif.pos.neg   = positivos-negativos
    mostly        = ifelse(positivos>negativos, "pos", "neg")
    suma          = sum(my_matrix)
    zeros         = length(which(my_matrix==0))
    avg           = mean(my_matrix)
    avg.pos       = mean(my_matrix[which(my_matrix>0)])
    avg.neg       = mean(my_matrix[which(my_matrix<0)])
    maxy          = max(my_matrix)
    miny          = min(my_matrix)
    med           = median(my_matrix)
    #signFactor: describes if there are more MI gains or losses
    signFactor    = (positivos-negativos)/(positivos+negativos+zeros)
    #change factor: abs(avg.pos)/abs(avg.neg) 
    #Larger than 1 if avg of gains is larger than avg of losses.
    #Less than 1 if avg of losses is larger than avg of gains
    changeFactor  = abs(avg.pos)/abs(avg.neg)
    intra.inter   = ifelse(i==j, "intra", "inter")
    #r = c(i,j,positivos, negativos, zeros, NA, suma)
    r = data.frame(cromo_i        = i,
                   cromo_j        = j,
                   pluses         = positivos,
                   minus          = negativos,
                   zeros          = zeros, 
                   mostly         = mostly,
                   dif.pos.neg    = dif.pos.neg,
                   sum            = suma, 
                   avg            = avg,
                   avg.pos        = avg.pos,
                   avg.neg        = avg.neg,
                   median         = med,
                   max            = maxy,
                   min            = miny,
                   signFactor     = signFactor,
                   changeFactor   = changeFactor,
                   intra.inter    = intra.inter
                   )
    my_frame = rbind(my_frame, r)
    
  }}

write.csv(x = my_frame, file = "results/basal_v_ctrl.csv", 
          quote = FALSE, 
          row.names = FALSE)

#make a scatter plot 
p = ggplot(data = my_frame, mapping = aes(signFactor, 
                                                 log(changeFactor), 
                                                 colour = factor(intra.inter),
                                                 label  = paste0(cromo_i, ";", cromo_j)
)
)
p = p + geom_point()
p = p + geom_text(check_overlap=TRUE, size = 5)
p = p + geom_hline(yintercept = 0) + geom_vline(xintercept = 0)
plot(p)
ggsave(filename = "results/scatter_signchangefactors_basalVhealth.pdf", 
       height = 8.5, 
       width = 14, 
       units = "in")


# for(i in cromos){
#   chromo_row = sort(filter(.data = annot2, Chr == i)$symbol)
#   for(j in cromos){
#     chromo_col = sort(filter(.data = annot2, Chr == j)$symbol)
#     xx = which(rownames(vs_basal_cntrl)%in%chromo_row)
#     yy = which(colnames(vs_basal_cntrl)%in%chromo_col)
#     
#     my_matrix = vs_basal_cntrl[xx,yy]
  #   my_reduced = my_matrix[,abs(colMeans(my_matrix))<abs(mean(my_matrix))]
  # 
  # 
  # 
  #   m = melt(as.matrix(my_reduced))
  #   m = melt(as.matrix(my_matrix))
  #   p = ggplot(m, aes(Var1,Var2, fill=value)) +
  #     geom_raster() +
  #     scale_fill_gradient2(low = "blue", high = "red") +
  #     theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 1)) +
  #     theme(axis.text.y = element_text(size = 1))+
  #     ggtitle(paste0("chromosome ", i, "versus ", j))
  #   plot(p)
  #   my_file = paste0("results/chromos",i,"_",j,".pdf")
  #   ggsave(filename = my_file)
  # }}

# xx = which(rownames(vs_basal_cntrl)%in%chromo_1)
# yy = which(colnames(vs_basal_cntrl)%in%chromo_1)
# my_test_matrix = vs_basal_cntrl[xx,yy]
# 
# max(my_test_matrix)
# min(my_test_matrix)
# mean(my_test_matrix)
# colSums(my_test_matrix)
# colMeans(my_test_matrix)
# 
# table(abs(colMeans(my_test_matrix))<abs(mean(my_test_matrix)))
# ncol(my_test_matrix)


