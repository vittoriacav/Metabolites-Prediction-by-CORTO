# To check changes of metabolite names
setwd("/home/student/Desktop/vittoria")
set_rownames <- function(df){
  rownames(df) <- df[,1]
  return (df[,-1])
}


load("data/final_ds/metab_NCI60.rda")
load("data/final_ds/metab_breast.rda")
load("data/final_ds/metab_hassan.rda")
load("data/final_ds/metab_zamp.rda")

# NCI60 ####
# no changes in NCI60
raw_NCI60 <- read.delim("data/datasets/NCI60/raw/metabData.csv", sep=",")
raw_NCI60 <- as.matrix(set_rownames(raw_NCI60))

setdiff(rownames(raw_NCI60),rownames(metab_NCI60)) # No changes in NCI60


# BREAST ####
# 6 names changed
raw_breast <- read.delim("data/datasets/breast/raw/metabData.csv", sep=",")
raw_breast <- as.matrix(set_rownames(raw_breast))

old <- setdiff(rownames(raw_breast), rownames(metab_breast)) 
new <- setdiff(rownames(metab_breast),rownames(raw_breast))

breast_metab_hist_names <- matrix(c(old, new), ncol = 2, nrow = length(old))
colnames(breast_metab_hist_names) <- c("ORIGINAL", "CHANGED")

save(breast_metab_hist_names, file="data/datasets/breast/changes/breast_hist.rda") 

# HASSAN ####
# 5 names changed
colspec <- c("skip", "skip", "guess", rep("skip", 20), rep("numeric", 36))
raw_hassan <- read_excel("data/datasets/hassan/raw/Metabolomics raw data.xlsx", col_types = colspec)
new_colnames <- c("X52_S26_A044", "X10_S5_A044", "X99_S23_A044", "X91_S30_A044", "X6_S4_A044", "X9_S14_A044", "X35_S18_A044", "X75_S20_A044",
                  "X36_S31_A044", "X103_S27_A044", "X84_S19_A044", "X82_S7_A044")

# Putting metabolites names in temporary variable
temp_metab_hassan <- raw_hassan[,1]
# Removing metabolites names from original metab df and make them numeric
raw_hassan <- raw_hassan[,-1]
raw_hassan <- mapply(raw_hassan, FUN=as.numeric)
# FUN = as.numeric cause the matrix to become one single vector so we have to rashape the matrix
raw_hassan <- matrix(data=raw_hassan, ncol=36, nrow=131)

# For each sample we have 3 measurements, with this for loop we calculate the mean of each triplicate
start <- 1
for(i in 1:12){
  end <- i*3
  new_col <- matrix(rowMeans(raw_hassan[,start:end]), ncol=1, nrow=131)
  temp_metab_hassan <- cbind(temp_metab_hassan, new_col)
  start <- end
}

temp_metab_hassan <- set_rownames(temp_metab_hassan)
colnames(temp_metab_hassan) <- new_colnames
raw_hassan <- as.matrix(temp_metab_hassan)

old <- setdiff(rownames(raw_hassan), rownames(metab_hassan))
new <- setdiff(rownames(metab_hassan),rownames(raw_hassan)) 

hassan_metab_hist_names <- matrix(c(old, new), ncol = 2, nrow = length(old))
colnames(hassan_metab_hist_names) <- c("ORIGINAL", "CHANGED")

save(hassan_metab_hist_names, file="data/datasets/hassan/changes/hassan_hist.rda")



# ZAMPIERI ####
# 51 names changed
load("data/datasets/zampieri/raw/metmat.rda")
raw_zamp <- metmat

old <- rownames(raw_zamp)
new <- rownames(metab_zamp)

zamp_metab_hist_names <- c("ORIGINAL", "CHANGED")

for(i in 1:length(old)){
  if(old[i] != new[i]){
    zamp_metab_hist_names <- rbind(zamp_metab_hist_names, c(old[i], new[i]))
  }
}

colnames(zamp_metab_hist_names) <- zamp_metab_hist_names[1,]
zamp_metab_hist_names <- zamp_metab_hist_names[-1,]

save(zamp_metab_hist_names, file="data/datasets/zampieri/changes/zamp_hist.rda")




