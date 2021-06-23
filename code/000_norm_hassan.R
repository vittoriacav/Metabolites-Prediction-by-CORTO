# NORMALIZE HASSAN RNAseq DATASET
#setwd("C:/Users/vitto/OneDrive - Alma Mater Studiorum Universit? di Bologna/UNIVERSITA'/Thesis/code")
setwd("/home/student/Desktop/vittoria")
source("code/vst.R")

raw <- read.delim("data/datasets/hassan/raw/GSE148892_a047.all.readcnt.known.genes.csv", sep=",", skip=9)
#sort(raw$GENE_NAME)[1:30]

# converting geneids
source("code/geneids.R")

ensID <- raw[,4] 
egID <- ens2eg(ensID)
symID <- eg2sym(egID)
raw$GENE_NAME <- symID

# Delete columns not needed
raw[1:4] <- NULL
raw$SOURCE <- NULL

#Checking for duplicates
genes <- raw[,1]
duplicates <- genes[duplicated(genes)]

# Deleting duplicates
raw_no_duplicates <- aggregate(raw[-1], by = list(raw$GENE_NAME), FUN = sum)
rownames(raw_no_duplicates) = raw_no_duplicates[, 1]
raw_no_duplicates = raw_no_duplicates[, -1] 

dim(raw) # 63772    27
dim(raw_no_duplicates) #33086    26

# conversion df to matrix
rawcounts <- as.matrix(raw_no_duplicates)

# running normalization
expmat <- vst(rawcounts)
write.table(expmat, file="data/datasets/hassan/raw/genes_hassan.txt", row.names=TRUE, col.names=TRUE, quote = FALSE, sep = "\t")


