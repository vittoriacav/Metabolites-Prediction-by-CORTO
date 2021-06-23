# Metabolic profiling of cancer cells reveals genomewide
# crosstalk between transcriptional regulators
# and metabolism
# Nature Communications 2019, Ortmayr, Dubuis and Zampieri
setwd("/home/student/Desktop/vittoria/code")

library(xlsx)


# Possible NCI-60 data (from Cellminer)
fnames<-read.delim("data/zampieri/U133.Plus2.filenames.txt")

# Relative metabolite abundance (from Ortmayr/Zampieri paper)
metabs<-read.xlsx2("data/zampieri/new/SuppData1_metabolites_onlymetabs.xlsx",1)
dim(metabs)
unique(sort(fnames[,2]))
unique(sort(colnames(metabs)[6:ncol(metabs)]))

# Make (by hand) a conversion table based on these
write.table(unique(sort(fnames[,2])),quote=FALSE,row=FALSE,col=FALSE)
write.table(unique(sort(colnames(metabs)[6:ncol(metabs)])),quote=FALSE,row=FALSE,col=FALSE)

# Reload id conversion table
idconv<-read.delim("data/zampieri/idconversion.txt",header=TRUE)



### Load and normalize microarray data with CustomCDF ----
# Get Custom CDF here http://brainarray.mbni.med.umich.edu/Brainarray/Database/CustomCDF/24.0.0/entrezg.asp
# Install it
install.packages("data/zampieri/hgu133plus2hsentrezgcdf_24.0.0.tar.gz",repos=NULL)


### THIS IS THE CUSTOM CDF package, it will have the name of the chip in it
source("../shared/functions/geneids.R")
source("../shared/functions/expmerger.R")
library("affy")
library("hgu133plus2hsentrezgcdf")
abatch<-ReadAffy(cdfname="hgu133plus2hsentrezgcdf",filenames=dir("data/zampieri/CEL-files/",pattern="cel",full.names=TRUE))
eset<-rma(abatch)
entrezmat<-as.matrix(exprs(eset))
# Remove "_at": they are now magically entrez ids
rownames(entrezmat)<-gsub("_at","",rownames(entrezmat))

# Now, conversion to symbols (multiple entrezs may match to the same currently valid gene symbol)
convlist<-eg2sym(rownames(entrezmat))
rawexpmat<-squish(entrezmat,convlist=convlist,method="average",verbose=TRUE)
save(rawexpmat,file="data/zampieri/rawexpmat.rda")


rm(list=ls())


################### RESTART HERE
library(xlsx)
### Now, save both datasets, but with the same column names ---
idconv<-read.delim("data/zampieri/idconversion.txt",header=TRUE)
idconv<-setNames(idconv[,2],idconv[,1])

# Format metabolite matrix
rawmetab<-read.xlsx2("data/zampieri/SuppData1_metabolites_onlymetabs.xlsx",1)
metmat<-rawmetab[,6:ncol(rawmetab)]
rownames(metmat)<-rawmetab[,3]
rm(rawmetab)

# Format transcript matrix
source("../shared/functions/expmerger.R")
load("data/zampieri/rawexpmat.rda")
dim(rawexpmat) # 20419    174
# We must merge by average samples mapping to the same cell type
fnames<-read.delim("data/zampieri/U133.Plus2.filenames.txt")
colnames(rawexpmat)<-gsub("\\.cel","",colnames(rawexpmat))
transpost<-t(rawexpmat)
convlist<-setNames(fnames$Cellminer.name,fnames$Sample.name)
transpost2<-squish(transpost,convlist=convlist,method="average",verbose=TRUE)
rawexpmat<-t(transpost2)
dim(rawexpmat) # 20419    59
colnames(rawexpmat)<-idconv[colnames(rawexpmat)]
rawexpmat<-rawexpmat[,!is.na(colnames(rawexpmat))]
dim(rawexpmat) # 20419    53
# 
common<-intersect(colnames(metmat),colnames(rawexpmat))
expmat<-rawexpmat[,common]
metmat<-metmat[,common]
dim(expmat) # 20419    53
dim(metmat) # 2181   53
#
metmat2<-apply(metmat,2,as.numeric)
dimnames(metmat2)<-dimnames(metmat)
metmat<-metmat2
save(metmat,file="data/zampieri/metmat.rda")
write.table(sort(rownames(metmat)),quote=FALSE,row=FALSE,col=FALSE,file="data/zampieri/metabolites.txt")
save(expmat,file="data/zampieri/expmat.rda")
mat<-rbind(metmat,expmat)
dim(mat) # 22600    53
save(mat,file="data/zampieri/mat.rda")



rm(list=ls())











