#load data 
basal = readRDS(file = "data/basal_adjmatrix.RDS")
cntrl = readRDS(file = "data/healthy_adjmatrix.RDS")
annot = readRDS(file = "data/gene_annotation.RDS")

#check similarity in gene order 

all(rownames(basal)==rownames(cntrl))

#Error!

which(rownames(basal)!=rownames(cntrl))
which(rownames(cntrl)!=rownames(basal))
#check names
rownames(basal)[which(rownames(basal)!=rownames(cntrl))]
rownames(cntrl)[which(rownames(basal)!=rownames(cntrl))]

rownames(basal)[13624:13627]
rownames(cntrl)[13624:13627]

#extra problem: some long non coding and some 

x = grep(pattern = "LINC", x = rownames(basal))
z = grep(pattern = "LINC", x = colnames(cntrl))

y = grep(pattern = "LINC", x = rownames(cntrl))
w = grep(pattern = "LINC", x = colnames(basal))

basal = basal[-x, -z]
cntrl = cntrl[-y, -w]

#does this fix the problem
all(sort(rownames(basal))==sort(rownames(cntrl)))
all(sort(colnames(basal))==sort(colnames(cntrl)))
#yes it does

basal = basal[order(rownames(basal)), order(colnames(basal))]
cntrl = cntrl[order(rownames(cntrl)), order(colnames(cntrl))]

basal[13624:13627, 13624:13627]
cntrl[13624:13627, 13624:13627]

##############
##not all annotated, let's remove those that aren't in the annotation
basal = basal[-which(!(rownames(basal)%in%annot$symbol)),
              -which(!(colnames(basal)%in%annot$symbol)),
              ]

cntrl = cntrl[-which(!(rownames(cntrl)%in%annot$symbol)),
              -which(!(colnames(cntrl)%in%annot$symbol)),
              ]

#add lower matrix triangle

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

#write matrices 
saveRDS(object = basal, file = "data/basalCleanMatrix.RDS")
saveRDS(object = cntrl, file = "data/cntrlCleanMatrix.RDS")
