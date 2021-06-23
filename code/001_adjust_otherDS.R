setwd("/home/student/Desktop/vittoria")

library("readxl")
library("corto")
library("VennDiagram")

source("code/clean_df.R")

set_rownames <- function(df){
  rownames(df) <- df[,1]
  return (df[,-1])
}

# 1. CLLE ####

# 1.1 Load and clean data ####
load("data/datasets/ccle/raw/ccle-metab.rda")
load("data/datasets/ccle/raw/ccle-expmat.rda")
metab_ccle <- metab
genes_ccle <- expmat

# dimension before removal of non-common cell lines
dim(genes_ccle) # 24274 genes  1019 cell lines
dim(metab_ccle) # 225 metabolites 928 cell lines

# Removing non-common cell lines 
clean_list_ccle <- clean_dfs(expmat, metab)
genes_ccle <- clean_list_ccle[[1]]
metab_ccle <- clean_list_ccle[[2]]
# dimension after removal of non-common cell lines
dim(genes_ccle) # 24274 genes  898 cell lines
dim(metab_ccle) # 225 metabolites 898 cell line

# Save cleaned data for later use
save(genes_ccle,file="data/final_ds/genes_ccle.rda")
save(metab_ccle,file="data/final_ds/metab_ccle.rda")

# 1.2 Plotting to check distributions of these data sets ####
png(paste0("results/figures/1_dist/ccle_genes.png"),w=2500,h=2000,res=600)
plot(density(genes_ccle), main = "Gene Expression\nDistribution CCLE", xlab = "mRNA Expression")
dev.off()
png(paste0("results/figures/1_dist/ccle_metab.png"),w=2500,h=2000,res=600)
plot(density(metab_ccle), main = "Metabolites Distribution CCLE", xlab = "metabolite measurement AU")
dev.off()


# 2. NCI60 SIDDAQUI ####

# 2.1 Load and clean data ####
genes_NCI60 <- read.delim("data/datasets/NCI60/raw/geneData.csv", sep=",")
genes_NCI60 <- as.matrix(set_rownames(genes_NCI60))
metab_NCI60 <- read.delim("data/datasets/NCI60/raw/metabData.csv", sep=",")
metab_NCI60 <- as.matrix(set_rownames(metab_NCI60))
# dimension before removal of non-common cell lines
dim(genes_NCI60) # 17987 genes   57 cell lines
dim(metab_NCI60) # 280 metabolites  58 cell lines

# Removing non-common cell lines 
clean_list_NCI60 <- clean_dfs(genes_NCI60, metab_NCI60)
genes_NCI60 <- clean_list_NCI60[[1]]
metab_NCI60 <- clean_list_NCI60[[2]]
# dimension after removal of non-common cell lines
dim(genes_NCI60) # 17987 genes   57 cell lines
dim(metab_NCI60) # 280 metabolites  57 cell lines

# Save cleaned data for later use 
save(genes_NCI60,file="data/final_ds/genes_NCI60.rda")
save(metab_NCI60,file="data/final_ds/metab_NCI60.rda")

# 2.2 Plotting to check distributions ####
png(paste0("results/figures/1_dist/NCI60_genes.png"),w=2500,h=2000,res=600)
plot(density(genes_NCI60), main = "Gene Expression\nDistribution NCI60 (Siddiqui)", xlab= "mRNA Expression")
dev.off()
png(paste0("results/figures/1_dist/NCI60_metab.png"),w=2500,h=2000,res=600)
plot(density(metab_NCI60), main = "Metabolites Distribution\nNCI60 (Siddiqui)", xlab = "metabolite measurement AU")
dev.off()


# 3.BREAST SIDDAQUI ####
# 3.1 Load and clean data ####
genes_breast <- read.delim("data/datasets/breast/raw/geneData.csv", sep=",")
genes_breast <- as.matrix(set_rownames(genes_breast))
metab_breast <- read.delim("data/datasets/breast/raw/metabData.csv", sep=",")
metab_breast <- as.matrix(set_rownames(metab_breast))
# dimension before removal of non-common cell lines
dim(genes_breast) # 20254 genes 108 cell lines
dim(metab_breast) # 536 metabolites 132 cell lines

# Removing non-common cell lines 
clean_list_breast <- clean_dfs(genes_breast, metab_breast)
genes_breast <- clean_list_breast[[1]]
metab_breast <- clean_list_breast[[2]]
# dimension after removal of non-common cell lines
dim(genes_breast) # 20254 genes 108 cell lines
dim(metab_breast) # 536 metabolites 108 cell lines

# Standardizing metabolite names to CCLE ones
real_names_breast <- as.matrix(read.delim("data/datasets/breast/changes/changed_names.txt", sep="\t", header = F))
rownames(metab_breast) <- real_names_breast

# Save cleaned data for later use 
save(genes_breast,file="data/final_ds/genes_breast.rda")
save(metab_breast,file="data/final_ds/metab_breast.rda")

# 3.2 Plotting to check distributions ####
png(paste0("results/figures/1_dist/breast_genes.png"),w=2500,h=2000,res=600)
plot(density(genes_breast), main = "Gene Expression\nDistribution GSE37751", xlab = "mRNA Expression")
dev.off()
png(paste0("results/figures/1_dist/breast_metab.png"),w=2500,h=2000,res=600)
plot(density(metab_breast), main = "Metabolites Distribution GSE37751", xlab = "metabolite measurement AU")
dev.off()

# 4. HASSAN ####
# 4.1 Load and clean data ####
genes_hassan <- as.matrix(read.delim("data/datasets/hassan/raw/genes_hassan.txt", sep="\t"))

# Hassan stored metabolite raw data in a strange format: "tbl_df"     "tbl"        "data.frame"
# The adjustment of this dataset in fact is a little different from the others and more laborious 
colspec <- c("skip", "skip", "guess", rep("skip", 20), rep("numeric", 36))
metab_hassan <- read_excel("data/datasets/hassan/raw/Metabolomics raw data.xlsx", col_types = colspec)
new_colnames <- c("X52_S26_A044", "X10_S5_A044", "X99_S23_A044", "X91_S30_A044", "X6_S4_A044", "X9_S14_A044", "X35_S18_A044", "X75_S20_A044",
                  "X36_S31_A044", "X103_S27_A044", "X84_S19_A044", "X82_S7_A044")

# Putting metabolites names in temporary variable
temp_metab_hassan <- metab_hassan[,1]
# Removing metabolites names from original metab df and make them numeric
metab_hassan <- metab_hassan[,-1]
metab_hassan <- mapply(metab_hassan, FUN=as.numeric)
# "FUN = as.numeric" causes the matrix to become one single consecutive vector so we have to reshape the matrix
metab_hassan <- matrix(data=metab_hassan, ncol=36, nrow=131)

# For each sample we have 3 measurements for each metabolite, with this for loop we calculate the mean of each triplicate
start <- 1
for(i in 1:12){
  end <- i*3
  new_col <- matrix(rowMeans(metab_hassan[,start:end]), ncol=1, nrow=131)
  # With the temporary variable created before (that initially contains metab names only)
  #  we add, column by column, the mean of each triplicate
  temp_metab_hassan <- cbind(temp_metab_hassan, new_col)
  start <- end
}


temp_metab_hassan <- set_rownames(temp_metab_hassan)
colnames(temp_metab_hassan) <- new_colnames
metab_hassan <- as.matrix(temp_metab_hassan)

# dimension before removal of non-common cell lines
dim(genes_hassan) # 33086 genes 26 cell lines
dim(metab_hassan) #131 metabolites 12 cell lines

# Removing non-common cell lines 
clean_list_hassan <- clean_dfs(genes_hassan, metab_hassan)
genes_hassan <- clean_list_hassan[[1]]
metab_hassan <- clean_list_hassan[[2]]
# dimension after removal of non-common cell lines
dim(genes_hassan) # 33086 genes 12 cell lines
dim(metab_hassan) #131 metabolites 12 cell lines

#Changing names
real_names_hassan <- as.matrix(read.delim("data/datasets/hassan/changes/changed_names.txt", sep="\t", header = F))
rownames(metab_hassan) <- real_names_hassan

# Save cleaned data for later use
save(genes_hassan,file="data/final_ds/genes_hassan.rda")
save(metab_hassan,file="data/final_ds/metab_hassan.rda")


# 4.2 Plotting to check distributions ####
png(paste0("results/figures/1_dist/hassan_genes"),w=2500,h=2000,res=600)
plot(density(genes_hassan), main = "Gene Expression\nDistribution GSE148892", xlab = "mRNA Expression")
dev.off()
png(paste0("results/figures/1_dist/hassan_metab"),w=2500,h=2000,res=600)
plot(density(log10(metab_hassan)), main = "Metabolites Distribution GSE148892", xlab = "metabolite measurement AU") # LOG10 !!!
dev.off()



# 5. ZAMPIERI ####
# 5.1 Load and clean data ####
load("data/datasets/zampieri/raw/metmat.rda")
load("data/datasets/zampieri/raw/expmat.rda")
genes_zamp <- expmat
metab_zamp <- metmat

# Standardizing metabolite names to CCLE ones
real_names_zamp <- as.matrix(read.delim("data/datasets/zampieri/changes/changed_names.txt", sep="\t"))
rownames(metab_zamp) <- real_names_zamp

# cell lines are already the same
setdiff(colnames(genes_zamp), colnames(metab_zamp))

# Checking dimension of dataset
dim(genes_zamp) #20419 genes 53 cell lines
dim(metab_zamp) # 2181 genes 53 cell lines

# Save cleaned data for later use 
save(genes_zamp,file="data/final_ds/genes_zamp.rda")
save(metab_zamp,file="data/final_ds/metab_zamp.rda")

# 5.2 Plotting to check distributions ####
png(paste0("results/figures/1_dist/zamp_genes.png"),w=2500,h=2000,res=600)
plot(density(genes_zamp), main = "Gene Expression\nDistribution GSE32474",  xlab = "mRNA Expression")
dev.off()
png(paste0("results/figures/1_dist/zamp_metab.png"),w=2500,h=2000,res=600)
plot(density(log10(metab_zamp), na.rm=T), main = "Metabolites Distribution GSE32474", xlab = "metabolite measurement AU")# LOG10 !!!
dev.off()

# VENN DIAGRAM ####

CCLE <- rownames(metab_ccle);
NCI60 <- rownames(metab_NCI60);
GSE37751 <- rownames(metab_breast);
GSE148892 <- rownames(metab_hassan);
GSE32474 <- rownames(metab_zamp);
x <- list(CCLE = CCLE, NCI60 = NCI60, GSE37751 = GSE37751, GSE148892 = GSE148892, GSE32474=GSE32474)
length_ds <- c(length(CCLE), length(NCI60), length(GSE37751), length(GSE148892), length(GSE32474))
names(length_ds) <- c("CCLE", "NCI60", "GSE37751 (br)", "GSE148892 (ha)", "GSE32474 (za)")

venn.diagram(x, filename = "results/figures/0_Venn/venn.png",
             col = "transparent", fill = c("brown1","cornflowerblue","yellow","darkorchid1","darkorange"),
             alpha = 0.50, label.col = "black", cex = 1, fontfamily = "serif", fontface = "bold",
             main.fontface = "bold", main.fontfamily = "sans",
             cat.col = c("brown3", "darkblue", "orange", "darkorchid4","darkorange3"), cat.cex = 1.1,
             cat.pos = 0, cat.dist = 0.07, cat.fontfamily = "sans", rotation.degree = 270,
             margin = 0.2, main = "Measured metabolites intersection\nbetween all datasets", 
             main.cex = 2, )

# just 2 metabolites are common to all data sets: serine and cytidine
intersect(GSE32474, intersect(GSE148892, intersect(intersect(CCLE, NCI60), GSE37751)))



