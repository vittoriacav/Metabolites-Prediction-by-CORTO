setwd("C:/Users/vitto/OneDrive - Alma Mater Studiorum Università di Bologna/UNIVERSITA'/Thesis/vittoria")
library(gplots)
library(corto)

quads <- function(ur, ul, br, bl){
  colours = c(ur, ul, br, bl)
  limits = par()$usr
  rect(0,0,limits[2],limits[4],col=colours[1])
  rect(0,0,limits[1],limits[4],col=colours[2])
  rect(0,0,limits[1],limits[3],col=colours[3])
  rect(0,0,limits[2],limits[3],col=colours[4])
}

pvalues_col <- function(cc, ps){
  #ps = pvalues
  col_vector <- c()
  for(i in 1:length(ps)){
    if(cc[i] < 0){ #neg ones
      if(ps[i] <= 0.001){
        col_vector <- c(col_vector, "cadetblue4")  
      }else if(ps[i] <= 0.01){
        col_vector <- c(col_vector, "cadetblue3")  
      }else if(ps[i] <= 0.05){
        col_vector <- c(col_vector, "cadetblue2")  
      }else if(ps[i] <= 0.1){
        col_vector <- c(col_vector, "cadetblue1")  
      }else{
        col_vector <- c(col_vector, "azure1")  
      }
    }else{ #pos ones
      if(ps[i]<= 0.001){
        col_vector <- c(col_vector, "brown4")  
      }else if(ps[i] <= 0.01){
        col_vector <- c(col_vector, "brown3")  
      }else if(ps[i] <= 0.05){
        col_vector <- c(col_vector, "brown2")  
      }else if(ps[i] <= 0.1){
        col_vector <- c(col_vector, "brown1")  
      }else{
        col_vector <- c(col_vector, "coral")  
      }
    }
  }
  return(col_vector)
}

correlation <- function(real, pred, thr){
  cors_spearman <- c()
  cors_pearson <- c()
  significant <- c()
  ps <- c()
  for(i in rownames(real)){
    cors_spearman <- c(cors_spearman,cor(real[i,],pred[i,],method="spearman"))
    cors_pearson <- c(cors_pearson,cor(real[i,],pred[i,],method="pearson"))
    ps <- c(ps, cor.test(real[i,],pred[i,], method = "spearman")$p.value)
    # Do you want to apply the threshold on spearman, pearson or both ?
    if(cors_spearman[length(cors_spearman)] > thr){
      significant <- c(significant, i)
    }
    #if(cors_pearson[length(cors_pearson)] > thr){
    #  significant <- c(significant, i)
    #}
    #if(cors_spearman[length(cors_spearman)] > thr && cors_pearson[length(cors_pearson)] > thr){
    #  significant <- c(significant, i)
    #}
  }
  names(ps) <- rownames(real)
  return(list(cors_spearman, cors_pearson, significant, ps))
}

# CCLE ####
# In this section I'll use CCLE data set to train and construct the predicting network with metabolites 
#  as centroids

# Load CCLE metabolite data, data used to create the network
load("data/final_ds/metab_ccle.rda")
names_metab_ccle <- rownames(metab_ccle)
#metab["alanine",]

# Load expression matrix
load("data/final_ds/genes_ccle.rda")
names_gene_ccle <- rownames(genes_ccle)

# Load expression + metabolite data (expression data is VST normalized)
load("data/datasets/ccle/raw/ccle-metabexp_mat.rda")
dim(mat) # 24499 molecular species (225 metabolites, 24274 transcripts) in 898 common samples
data_for_network_ccle <- mat

# Generate an overall network (it will take a while, adjust the nthreads according to your pc)
#network <- corto(data_for_network_ccle,centroids=names_metab_ccle, nbootstraps=100, p=1e-10, nthreads=7, verbose=FALSE)
#save(network,file="code/network.rda")
load("code/network.rda")

# EXPLORE NETWORK ####
#match("C16:0_CE", names(network))
#network[["xanthosine"]]
#network[["xanthosine"]]["tfmode"] #network[["xanthosine"]][1]
#network[["xanthosine"]]["likelihood"] #network[["xanthosine"]][2]

# Single gene e.g. xanthosine --> PKD1
#network[["xanthosine"]][["likelihood"]][104]
#network[["xanthosine"]][["tfmode"]][104]

# get all targets genes of one met
#names(network[["C22:1 SM"]][["tfmode"]])

#names(network[["xanthosine"]][["tfmode"]])

# FINDING BEST CORRELATION IN THE NETWORK
best <- list(metab = "", gene = "", tfmode = 0, likelihood = 0)
nodes <- names(network)
for(metab in nodes){
  for(i in 1:length(network[[metab]][["tfmode"]]))
  if(network[[metab]][["tfmode"]][[i]] >= best$tfmode && 
     network[[metab]][["likelihood"]][[i]] >= best$likelihood){
    best$metab <- metab
    best$gene <- names(network[[metab]][["tfmode"]])[i]
    best$tfmode <- network[[metab]][["tfmode"]][[i]]
    best$likelihood <- network[[metab]][["likelihood"]][[i]]
  }
}


# EDGE MEANING explained with correlation plot ####
# 1-methylnicotinamide and NNMT

png(paste0("results/figures/2_networks/edge_explanation.png"),w=3200,h=2500,res=600)
plot(genes_ccle[best$gene,], metab_ccle[best$metab,], xlab = paste0(best$gene," expression"), 
     ylab = paste0(best$metab," measurement AU"),
     pch = 18, main = paste0(best$gene," - ", toupper(best$metab)," CORRELATION"), col = "darkgreen", panel.first = grid(8,8))
fit <- lm(metab_ccle[best$metab,]~genes_ccle[best$gene,])
mtext(paste0("CC = ",round(cor(genes_ccle[best$gene,], metab_ccle[best$metab,], method="s"), 3)))
abline(fit$coefficients)
dev.off()

# METABOLITE LEVEL PREDICTION ####

# Set correlation coefficient threshold to later filter the non-significant predictions
cc_thr <- 0.3

# 1. NCI60 ####
load("data/final_ds/genes_NCI60.rda")
load("data/final_ds/metab_NCI60.rda")

real_metab_names_NCI60 <- rownames(metab_NCI60)

# 1.1 CCLE-based prediction####
pred_NCI60 <- mra(genes_NCI60, regulon = network)
pred_metab_names_NCI60 <- rownames(pred_NCI60)
common_NCI60 <- intersect(pred_metab_names_NCI60, real_metab_names_NCI60)

length(common_NCI60) # 21 correctly predicted metabolites, however...

# 1.2 Correlation checks ####
# Correlation between real and predicted
corr_outputs_NCI60 <- correlation(metab_NCI60[common_NCI60,], pred_NCI60[common_NCI60,], cc_thr)
cors_NCI60_spearman <- corr_outputs_NCI60[[1]]
cors_NCI60_pearson <- corr_outputs_NCI60[[2]]

names(cors_NCI60_spearman) <- common_NCI60
names(cors_NCI60_pearson) <- common_NCI60

significant_NCI60 <- corr_outputs_NCI60[[3]]

ps_NCI60 <- corr_outputs_NCI60[[4]]

length(significant_NCI60) # ...however just just 10/21 are >0.3 (Spearman Correlation)

# Saving with different names to better manage the step in which I will compare the correlations performances of 
#  different tools
NCI60_sp_corto <- cors_NCI60_spearman
NCI60_pe_corto <- cors_NCI60_pearson

save(NCI60_sp_corto, file = "data/final_ds/corr/NCI60_sp_corto.rda")
save(NCI60_pe_corto, file = "data/final_ds/corr/NCI60_pe_corto.rda")

# 1.3 Testing distribution ####
# Testing how much dist is diff from zero
shapiro.test(cors_NCI60_spearman) #p-value = 0.1316 --> can use t-test
shapiro.test(cors_NCI60_pearson) #p-value = 0.352 --> can use t-test
t1_NCI60 <- t.test(cors_NCI60_spearman, alternative = "greater") #p-value = 0.0003278619
t2_NCI60 <- t.test(cors_NCI60_pearson, alternative = "greater") # p-value = 0.0001724943

# 1.4 Plotting correlation distributions ####

png(paste0("results/corr_distributions/correlations_NCI60.png"),w=3000,h=2500,res=600)
plot(density(cors_NCI60_spearman),col="orange red",lwd=3, xlim=c(-1,1), ylim = c(0,1.7), 
     main = "NCI60 correlations distribution", xlab = "metabolite CC")
lines(density(cors_NCI60_pearson), lty = 2, lwd = 2)
abline(v=0,lty=2, col="gray")
mtext(paste0("Spearman p-value = ",round(t1_NCI60$p.value,3),"   Pearson p-value = ",round(t2_NCI60$p.value,3)))
legend("topleft", legend = c("Spearman", "Pearson"), lwd = c(3,2),
       col = c("orange red", "black"), lty = c(1,2))
dev.off()

# Same image but little different format (for thesis)
png(paste0("results/figures/3_prediction_performace/correlations_NCI60.png"),w=2900,h=2500,res=600)
plot(density(cors_NCI60_spearman),col="darkgoldenrod3",lwd=3, xlim=c(-1,1), ylim = c(0,1.7), 
     main = "NCI60 correlations ditribution", xlab = "metabolite CC")
lines(density(cors_NCI60_pearson), col = "brown3", lty = 3, lwd = 2)
abline(v=0,lty=2, col="gray")
mtext(paste0("Spearman p-value = ",round(t1_NCI60$p.value,3),"     Pearson p-value = ",round(t2_NCI60$p.value,3)))
legend("topleft", legend = c("Spearman", "Pearson"), lwd = 2,
       col = c("darkgoldenrod3", "brown3"), lty = c(1,2))
dev.off()

# 1.5 Barplot of metabolite correlations ####

sorted_spear <- sort(cors_NCI60_spearman, decreasing = T)
sorted_pear <- cors_NCI60_pearson[names(sorted_spear)]
sorted_ps <- ps_NCI60[names(sorted_spear)]

all_sorted <- rbind(sorted_spear, sorted_pear)

png(paste0("results/figures/3_prediction_performace/barplot.png"),w=3000,h=2800,res=600)
op <- par(no.readonly = TRUE)
par(mar=c(6.5, 4, 2, 2) + 0.1)
bp <- barplot(all_sorted,beside=T, 
              las = 2, 
              ylim = c(min(all_sorted)-0.15, max(all_sorted)+0.2), 
              col = c("darkgoldenrod3", "brown3"),
              main="Metabolite CCs in NCI60",
              ylab = "Correlation Coefficient"
              )
grid(nx=NA, ny=NULL)

bp <- barplot(all_sorted,beside=T, 
              las = 2, 
              ylim = c(0, max(all_sorted)+0.1), 
              col = c("darkgoldenrod3", "brown3"),
              main="Metabolite CCs in NCI60",
              add = T,
              ylab = "Correlation Coefficient")
legend("topright", legend = c("Spearman", "Pearson"),
       fill = c("darkgoldenrod3", "brown3"))

x_c <- 1.5
for(i in 1:length(colnames(all_sorted))){
  if(sorted_ps[i] < 0.05){
    text(x_c, sorted_spear[i] + 0.03, "*")
  }
  x_c <- x_c + 3
}

par(op)

dev.off()



# 1.6 Scatterplot of each single correlation ####
c <- 1
for(i in common_NCI60){
  png(paste0("results/test_NCI60/test_",c,".png"),w=2900,h=2500,res=600)
  scatter(pred_NCI60[i,],metab_NCI60[i,], xlab = paste0("measured_",i), ylab = paste0("pred_",i), method = "spearman", 
          main = paste0(toupper(i), "\nmeasured vs predicted in NCI60"))
  dev.off()
  c <- c+1
}

# Chosen scatter of Sorbitol for thesis
png(paste0("results/figures/3_prediction_performace/scatter_sorbitol.png"),w=2500,h=2500,res=700)
scatter(pred_NCI60["sorbitol",],metab_NCI60["sorbitol",], xlab = "Measured Sorbitol in NCI60", 
        ylab = "Sorbitol CCLE-based prediction in NCI60", 
        method = "spearman", main = "Sorbitol prediction in NCI60")
dev.off()


# 2. BREAST ####
load("data/final_ds/genes_breast.rda")
load("data/final_ds/metab_breast.rda")
real_names_breast <- rownames(metab_breast)

# 2.1 Checking phenotype ####
pheno_breast_df <- read.delim("data/datasets/breast/raw/pData.csv", sep=",")
pheno_breast_df <- set_rownames(pheno_breast_df)
pheno_breast <- c()
for(i in colnames(genes_breast)){
  if(pheno_breast_df[i, "DIAG"] == "NORMAL"){
    pheno_breast <- c(pheno_breast, "n") # n = non-tumoral
  } else if(pheno_breast_df[i, "DIAG"] == "TUMOR"){
    pheno_breast <- c(pheno_breast, "t") # t = tumoral
  }
}

#     1:47 non-tumoral phenotype
#     48:108 tumoral phenotype

# 2.2 Master regulator analysis ####
# TUMORAL vs NON-TUMORAL
#differential_breast <- mra(genes_breast[,48:108], genes_breast[,1:47], regulon = network, nthreads = 6)
#save(differential_breast, file = "code/differential_breast.rda")
load("code/differential_breast.rda")

# 2.3 MRA plot ####
png(paste0("results/figures/5_contrast_analysis/MRA_breast.png"),w=4000,h=2000,res=300)
mraplot(differential_breast, title = "MRA on tumoral vs non-tumoral cell lines")
dev.off()

mra_names_breast <- names(differential_breast$regulon)
common_breast_mra <- intersect(mra_names_breast, real_names_breast)

# 2.4 Boxplot tumoral vs non-tumoral ####
# I calculate the average measurement for each metabolite in both groups (tumor, non-tumor)
tumor <- rowMeans(metab_breast[,48:108])
non.tumor <- rowMeans(metab_breast[,1:47])

png(paste0("results/figures/5_contrast_analysis/boxplot_breast.png"),w=2500,h=2500,res=700)
boxplot(non.tumor, tumor, col = c("#E1FAD3","#FADAD3"), names = c("Normal", "Breast Cancer"),
        main = "Tumoral vs Non-Tumoral\nboxplot in GSE37751", ylab= "average metabolite measurement")
dev.off()

# 2.4.2 Scatter replacing boxplot to show that all metab are measured as up-regulated in  this dataset ####

tumor_common <- tumor[common_breast_mra]
non.tumor_common <- non.tumor[common_breast_mra]

png(paste0("results/figures/5_contrast_analysis/scatter_breast.png"),w=2500,h=2500,res=500)
plot(non.tumor, tumor,  xlab = "non-tumoral metabolite measurement", 
     ylab = "tumoral metabolite measurement", type = "n", 
     xlim = c(-max(abs(tumor)),max(abs(tumor))),
     ylim = c(-max(abs(non.tumor)),max(abs(non.tumor))),
     main = "Non-tumoral vs Tumoral\nmetabolite measurement in GSE37751")
quads("white","#FADAD3","white","#C9CEF2")
points(non.tumor, tumor,pch = 19)
points(non.tumor_common, tumor_common,pch = 19, col="chartreuse3")
abline(coef = c(0,1), lty = 2)
legend("top", fill = "chartreuse3", legend = "metabolites involved in prediction", cex = 0.6 )
dev.off()

# 2.5 Plotting pred vs real in MRA ####

# NES: Normalized Enrichment Score -> measure of up/down-regulation
nes_y_breast <- differential_breast[["nes"]]
nes_y_breast <- nes_y_breast[common_breast_mra]

# Using t-test to calculate the level of up/down-regulation with the real measured data
ttest_x_breast <- c()
for(i in common_breast_mra){
  ttest_x_breast <- c(ttest_x_breast, t.test(metab_breast[i,48:108], metab_breast[i,1:47])$stat)
}

png(paste0("results/figures/5_contrast_analysis/breast_scatter_Ttest.png"),w=2500,h=2500,res=600)
plot(ttest_x_breast, nes_y_breast,  xlab = "t-test(tumoral vs non-tumoral)", 
     ylab = "Normalized Enrichment Score (NES)", type = "n", 
     xlim = c(-max(abs(ttest_x_breast)),max(abs(ttest_x_breast))),
     ylim = c(-max(abs(nes_y_breast)),max(abs(nes_y_breast))),
     main = "Measured vs Predicted\nup/down-regulation of metabolites\nin tumoral samples of GSE37751")
quads("#E1FAD3","#FADAD3","#E1FAD3","#FADAD3")
points(ttest_x_breast, nes_y_breast,pch = 19)
dev.off()

# Using log2(fold change) to calculate the level of up/down-regulation with the real measured data
# Shifting the distribution to avoid negative values
metab_breast_shifted <- metab_breast + abs(min(metab_breast))
log2fc_x <- c()
for(i in common_breast_mra){
  log2fc_x <- c(log2fc_x, log2((median(metab_breast_shifted[i,48:108], na.rm = T)/median(metab_breast_shifted[i,1:47], na.rm = T))))
} 

png(paste0("results/figures/5_contrast_analysis/breast_scatter_log2FC.png"),w=2500,h=2500,res=600)
plot(log2fc_x, nes_y_breast,  xlab = "log2(FC)", 
     ylab = "Normalized Enrichment Score (NES)", type = "n", xlim = c(-max(abs(log2fc_x)),max(abs(log2fc_x))),
     ylim = c(-max(abs(nes_y_breast)),max(abs(nes_y_breast))),
     main = "Measured vs Predicted\nup/down-regulation of metabolites\nin tumoral samples of GSE37751")
quads("#E1FAD3","#FADAD3","#E1FAD3","#FADAD3")
points(log2fc_x, nes_y_breast,pch = 19)
dev.off()



# 2.6 CCLE-based prediction ####
# Prediction in tumoral samples only
#metab_breast_tumoral <- metab_breast[,48:108]
pred_breast <- mra(genes_breast, regulon = network)
pred_metab_names_breast <- rownames(pred_breast)

#common in both normal pred and MRA
common_breast_pred <- intersect(pred_metab_names_breast, real_names_breast)

length(common_breast_pred) # 43 correctly predicted metabolites, however...

# 2.7 Correlation checks ####

corr_outputs_breast <- correlation(metab_breast[common_breast_pred,], pred_breast[common_breast_pred,], cc_thr)
cors_breast_spearman <- corr_outputs_breast[[1]]
cors_breast_pearson <- corr_outputs_breast[[2]]

names(cors_breast_spearman) <- common_breast_pred
names(cors_breast_pearson) <- common_breast_pred

significant_breast<- corr_outputs_breast[[3]]

length(significant_breast) # ...however just just 3/43 are >0.3

# Saving with different names to better manage the step in which I will compare the correlations performances of 
#  different tools
breast_sp_corto <- cors_breast_spearman
breast_pe_corto <- cors_breast_pearson

save(breast_sp_corto, file = "data/final_ds/corr/breast_sp_corto.rda")
save(breast_pe_corto, file = "data/final_ds/corr/breast_pe_corto.rda")

# 2.8 Testing distribution ####
# Testing how much dist is diff from zero
shapiro.test(cors_breast_spearman) #p-value = 0.423 --> can t-test
shapiro.test(cors_breast_pearson) #p-value = 0.4316 --> can t-test
t1_breast<- t.test(cors_breast_spearman, alternative = "greater") #p-value =  0.995302
t2_breast <- t.test(cors_breast_pearson, alternative = "greater") # p-value =  0.9937332

# 2.9 Plotting correlation distributions ####
png(paste0("results/corr_distributions/correlations_breast.png"),w=3000,h=2500,res=600)
plot(density(cors_breast_spearman),col="orange red",lwd=3, xlim=c(-1,1), ylim = c(0,2), 
     main = "GSE37751 correlations distribution", xlab = "metabolite CC")
lines(density(cors_breast_pearson), lty = 2, lwd = 2)
abline(v=0,lty=2, col="gray")
mtext(paste0("Spearman p-value = ",round(t1_breast$p.value,3)," Pearson p-value = ",round(t2_breast$p.value,3)))
legend("topleft", legend = c("Spearman", "Pearson"), lwd = c(3,2),
       col = c("orange red", "black"), lty = c(1,2))
dev.off()

# 2.10 Scatterplot of each single correlation ####
c <- 1
for(i in common_breast_pred){
  png(paste0("results/test_breast/test_",c,".png"),w=2900,h=2500,res=600)
  scatter(pred_breast[i,],metab_breast[i,], xlab = paste0("measured_",i), ylab = paste0("pred_",i), 
          method = "spearman", main = paste0(toupper(i), "\npredicted vs measured in GSE37751")) 
  dev.off()
  c <- c+1
}



# 3. HASSAN ####

load("data/final_ds/genes_hassan.rda")
load("data/final_ds/metab_hassan.rda")
real_names_hassan <- rownames(metab_hassan)

# 3.1 Checking phenotype (manually) ####
# non obese: 52_S26, 10_S5, 99_S23, 91_S30, 6_S4, 9_S14 [,1:6]
# obese: 35_S18, 75_S20, 36_S31, 103_S27, 84_S19, 82_S7 [,7:12]

# 3.2 Master regulator analysis ####
#differential_hassan <- mra(genes_hassan[,7:12], genes_hassan[,1:6], regulon = network)
#save(differential_hassan,file="differential_hassan.rda")
load("code/differential_hassan.rda")

# 3.3 MRA plot ####
png(paste0("results/figures/5_contrast_analysis/MRA_hassan.png"),w=4000,h=2000,res=300)
mraplot(differential_hassan, pthr = 0.05, title = "MRA on obese vs non-obese")
dev.off()

# 3.4 Boxplot obese vs non-obese ####
obese <- rowMeans(metab_hassan[,7:12])
non.obese <- rowMeans(metab_hassan[,1:6])

png(paste0("results/figures/5_contrast_analysis/boxplot_hassan.png"),w=2500,h=2500,res=700)
boxplot(log10(non.obese), log10(obese), col = c("#E1FAD3","#FADAD3"), names = c("Non-obese", "Obese"),
        main = "Obese vs Non-Obese\nboxplot in GSE148892", ylab = "log10(average metabolite measurements)")
dev.off()

# 3.5 Plotting pred vs real in MRA ####
mra_names_hassan <- names(differential_hassan$regulon)
common_hassan_mra <- intersect(mra_names_hassan, real_names_hassan)

nes_y_hassan <- differential_hassan[["nes"]]
nes_y_hassan <- nes_y_hassan[common_hassan_mra]

nes_y_hassan_pvalues <- differential_hassan[["pvalue"]]
nes_y_hassan_pvalues <- nes_y_hassan_pvalues[common_hassan_mra]

# Using t-test
ttest_x_hassan <- c()
ttest_x_hassan_pvalues <- c()
for(i in common_hassan_mra){
  t <- t.test(metab_hassan[i,7:12], metab_hassan[i,1:6])
  ttest_x_hassan <- c(ttest_x_hassan, t$stat)
  ttest_x_hassan_pvalues <- c(ttest_x_hassan_pvalues, t$p.value)
}
names(ttest_x_hassan) <-common_hassan_mra
names(ttest_x_hassan_pvalues) <-common_hassan_mra

png(paste0("results/figures/5_contrast_analysis/hassan_scatter_Ttest.png"),w=2500,h=2500,res=600)
plot(ttest_x_hassan, nes_y_hassan,  xlab = "t-test(obese vs non-obese)", 
     ylab = "Normalized Enrichment Score (NES)", type = "n", xlim = c(-max(abs(ttest_x_hassan)),max(abs(ttest_x_hassan))),
     ylim = c(-max(abs(nes_y_hassan)),max(abs(nes_y_hassan))),
     main = "Measured vs Predicted\nup/down-regulation of metabolites\nin obese samples of GSE148892")
quads("#E1FAD3","#FADAD3","#E1FAD3","#FADAD3")
points(ttest_x_hassan, nes_y_hassan,pch = 19)
dev.off()

# Using log2(fold change)
metab_hassan_shifted <- metab_hassan + abs(min(metab_hassan))
log2fc_x <- c()
for(i in common_hassan_mra){
  log2fc_x <- c(log2fc_x, log2((median(metab_hassan_shifted[i,7:12], na.rm = T)/median(metab_hassan_shifted[i,1:6], na.rm = T))))
} 

png(paste0("results/figures/5_contrast_analysis/hassan_scatter_log2FC.png"),w=2500,h=2500,res=600)
plot(log2fc_x, nes_y_hassan,  xlab = "log2(FC)", 
     ylab = "Normalized Enrichment Score (NES)", type = "n", xlim = c(-max(abs(log2fc_x)),max(abs(log2fc_x))),
     ylim = c(-max(abs(nes_y_hassan)),max(abs(nes_y_hassan))),
     main = "Measured vs Predicted\nup/down-regulation of metabolites\nin obese samples of GSE148892")
quads("#E1FAD3","#FADAD3","#E1FAD3","#FADAD3")
points(log2fc_x, nes_y_hassan,pch = 19)
dev.off()

# 3.5.2 Barplot pred vs real in MRA ####

ttest_x_hassan <- sort(ttest_x_hassan)
nes_y_hassan <- nes_y_hassan[names(ttest_x_hassan)]
png(paste0("results/figures/5_contrast_analysis/hassan_barplot.png"),w=7000,h=4000,res=600)
lim <- max(abs(min(ttest_x_hassan)), abs(max(ttest_x_hassan)), abs(min(nes_y_hassan)), abs(max(nes_y_hassan)))
par(mfrow=c(1,2))
barplot(ttest_x_hassan, 
        ylim = c(-lim-1, lim+1), 
        col = pvalues_col(ttest_x_hassan, ttest_x_hassan_pvalues[names(ttest_x_hassan)]),
        main="Measured up/down-regulation in GSE148892",
        ylab = "t-test(obese vs non-obese)",
        names.arg = names(ttest_x_hassan)
)
abline(h=0)
legend("top", legend = c(rep(c("<0.001", "<0.01", "<0.05", "<0.1", "<1"), 2)), ncol = 2, 
       fill = c("cadetblue4","cadetblue3", "cadetblue2", "cadetblue1", "azure1", "brown4",
                "brown3","brown2","brown1","coral"),
       title = "t-test p-value"
)
barplot(nes_y_hassan, 
        ylim = c(-lim-1, lim+1), 
        col = pvalues_col(nes_y_hassan, nes_y_hassan_pvalues[names(ttest_x_hassan)]),
        main="Predicted up/down-regulation in GSE148892",
        ylab = "NES",
        names.arg = names(nes_y_hassan)
)
abline(h=0)
legend("top", legend = c(rep(c("<0.001", "<0.01", "<0.05", "<0.1", "<1"), 2)), ncol = 2, 
       fill = c("cadetblue4","cadetblue3", "cadetblue2", "cadetblue1", "azure1", "brown4",
                "brown3","brown2","brown1","coral"),
       title = "NES p-value"
)
dev.off()

# 3.6 CCLE-based prediction ####

pred_hassan <- mra(genes_hassan, regulon = network)
pred_metab_names_hassan <- rownames(pred_hassan)
common_hassan_pred <- intersect(pred_metab_names_hassan, real_names_hassan)

length(common_hassan_pred) # 5 correctly predicted metabolites, however...

# 3.7 Correlation checks #### 

corr_outputs_hassan <- correlation(metab_hassan[common_hassan_pred,], pred_hassan[common_hassan_pred,], cc_thr)
cors_hassan_spearman <- corr_outputs_hassan[[1]]
cors_hassan_pearson <- corr_outputs_hassan[[2]]

names(cors_hassan_spearman) <- common_hassan_pred
names(cors_hassan_pearson) <- common_hassan_pred

significant_hassan<- corr_outputs_hassan[[3]]

length(significant_hassan) # ...however just 3/5 are > 0.3

# Saving with different names to better manage the step in which I will compare the correlations performances of 
#  different tools
hassan_sp_corto <- cors_hassan_spearman
hassan_pe_corto <- cors_hassan_pearson

save(hassan_sp_corto, file = "data/final_ds/corr/hassan_sp_corto.rda")
save(hassan_pe_corto, file = "data/final_ds/corr/hassan_pe_corto.rda")

# 3.8 Testing distribution ####
# Testing how much dist is diff from zero
shapiro.test(cors_hassan_spearman) #p-value = 0.09383 --> can use t-test
shapiro.test(cors_hassan_pearson) #p-value = 0.76 --> can t-test
t1_hassan <- t.test(cors_hassan_spearman, alternative = "greater") #p-value = 0.2864553
t2_hassan <- t.test(cors_hassan_pearson, alternative = "greater") # p-value = 0.2273071

# 3.9 Plotting correlation distributions ####
png(paste0("results/corr_distributions/correlations_hassan.png"),w=2900,h=2500,res=600)
plot(density(cors_hassan_spearman),col="orange red",lwd=3, xlim=c(-1,1), ylim = c(0,1.5), 
     main = "GSE148892 correlations distribution", xlab = "metabolite CC")
lines(density(cors_hassan_pearson), lty = 2, lwd = 2)
abline(v=0,lty=2, col="gray")
mtext(paste0("Spearman p-value = ",round(t1_hassan$p.value,3)," Pearson p-value = ",round(t2_hassan$p.value,3)))
legend("topleft", legend = c("Spearman", "Pearson"), lwd = c(3,2),
       col = c("orange red", "black"), lty = c(1,2))
dev.off()


# 3.10 Scatterplot of each single correlation ####
c <- 1
for(i in common_hassan_pred){
  png(paste0("results/test_hassan/test_",c,".png"),w=2900,h=2500,res=600)
  scatter(pred_hassan[i,],metab_hassan[i,], xlab = paste0("measured_",i), ylab = paste0("pred_",i), 
          method = "spearman", main = paste0(toupper(i),"\npredicted vs measured in GSE148892")) 
  dev.off()
  c <- c + 1
}




# 4. ZAMPIERI ####
load("data/final_ds/genes_zamp.rda")
load("data/final_ds/metab_zamp.rda")

real_metab_names_zamp <- rownames(metab_zamp)

# 4.1 CCLE-based prediction####
pred_zamp <- mra(genes_zamp, regulon = network)
pred_metab_names_zamp <- rownames(pred_zamp)
common_zamp <- intersect(pred_metab_names_zamp, real_metab_names_zamp)

length(common_zamp) # 48 correctly predicted metabolites, however...

# 4.2 Correlation checks ####
# Need to manage this dataset in a different way to calculate the correlation due to missing values 
cors_zamp_spearman<-c()
cors_zamp_pearson<-c()
significant_zamp <- c()
for(j in common_zamp){
  temp_real <- c()
  temp_pred <- c()
  for(v in 1:length(metab_zamp[j,])){
    if(metab_zamp[j,v] != "NaN"){
      temp_real <- c(temp_real, metab_zamp[j,v])
      temp_pred <- c(temp_pred, pred_zamp[j,v])
    }
  }
  cors_zamp_spearman<-c(cors_zamp_spearman,cor(temp_real,temp_pred,method="spearman"))
  cors_zamp_pearson<-c(cors_zamp_pearson,cor(temp_real,temp_pred,method="pearson"))
  if(cors_zamp_spearman[length(cors_zamp_spearman)]>cc_thr && cors_zamp_pearson[length(cors_zamp_pearson)]>cc_thr){
    significant_zamp <- c(significant_zamp, j)
  }
}
names(cors_zamp_spearman) <- common_zamp
names(cors_zamp_pearson) <- common_zamp

length(significant_zamp) # ...however just 6/48 are > 0.3

# Saving with different names to better manage the step in which I will compare the correlations performances of 
#  different tools
zamp_sp_corto <- cors_zamp_spearman
zamp_pe_corto <- cors_zamp_pearson

save(zamp_sp_corto, file = "data/final_ds/corr/zamp_sp_corto.rda")
save(zamp_pe_corto, file = "data/final_ds/corr/zamp_pe_corto.rda")



# 4.3 Testing distribution ####
# Testing how much dist is diff from zero
shapiro.test(cors_zamp_spearman) #p-value = 0.9398 --> can t-test
shapiro.test(cors_zamp_pearson) #p-value = 0.832 --> can t-test 
t1_zamp <- t.test(cors_zamp_spearman, alternative = "greater") # p-value = 0.01126951
t2_zamp <- t.test(cors_zamp_pearson, alternative = "greater") # p-value = 0.01919159 

# 4.4 Plotting correlation distributions ####
png(paste0("results/corr_distributions/correlations_zamp.png"),w=2900,h=2500,res=600)
plot(density(cors_zamp_spearman),col="orange red",lwd=3, xlim=c(-1,1), ylim = c(0,2), 
     main = "GSE32474 correlations distribution", xlab = "metabolite CC")
lines(density(cors_zamp_pearson), lty = 2, lwd = 2)
abline(v=0,lty=2, col="gray")
mtext(paste0("Spearman p-value = ",round(t1_zamp$p.value,3),"   Pearson p-value = ",round(t2_zamp$p.value,3)))
legend("topleft", legend = c("Spearman", "Pearson"), lwd = c(3,2),
       col = c("orange red", "black"), lty = c(1,2))
dev.off()

# 4.5 Scatterplot of each single correlation ####
c <- 1
for(i in common_zamp){
  png(paste0("results/test_zamp/test_",c,".png"),w=1500,h=1500,res=300)
  scatter(pred_zamp[i,],metab_zamp[i,], xlab = paste0("measured_",i), ylab = paste0("pred_",i), 
          method = "spearman",main = paste0(toupper(i), "\npredicted vs measured in GSE32474")) 
  dev.off()
  c <- c+1
}

# URACIL CANDIDATE ####
# Hassan not included, otherwise no metabolite will be in all datasets with CC > 0.3
common_predictions <-intersect(common_zamp,intersect(common_breast_pred, common_NCI60))

# Isolating good candidates
candidates <- c()
for(i in common_predictions){
  if(cors_NCI60_spearman[i] >= 0.3 && cors_breast_spearman[i] >= 0.3 && cors_zamp_spearman[i] >= 0.3){
    candidates <- c(candidates, i) # URACIL only (excluding Hassan, not measured there)
  }
}

#For thesis
png(paste0("results/figures/4_other_datasets/NCI60_uracil.png"),w=1500,h=1500,res=400)
scatter(pred_NCI60["uracil",],metab_NCI60["uracil",], xlab = "Measured Uracil", ylab = "Predicted Uracil", 
        method = "spearman", main = "Uracil Prediction in NCI60") 
dev.off()

png(paste0("results/figures/4_other_datasets/breast_uracil.png"),w=1500,h=1500,res=400)
scatter(pred_breast["uracil",],metab_breast["uracil",], xlab = "Measured Uracil", ylab = "Predicted Uracil", 
        method = "spearman", main = "Uracil Prediction in GSE37751") 
dev.off()

png(paste0("results/figures/4_other_datasets/zamp_uracil.png"),w=1500,h=1500,res=400)
scatter(pred_zamp["uracil",],metab_zamp["uracil",], xlab = "Measured Uracil", ylab = "Predicted Uracil", 
        method = "spearman", main = "Uracil Prediction in GSE32474") 
dev.off()


# CYTOSCAPE ####
# Get names of the metabolites nodes present in the network 
nodes <- names(network) # metabolite nodes
nodes <- gsub(" ", "_", nodes) #substituting the space with underscore to remove problems in Cytoscape
names(network) <- nodes

# Creating txt matrix for the whole network
# each row contains: metabolite name; gene name; tfmode; likelihood

net_matrix<-matrix(nrow=0,ncol=4)

pb<-txtProgressBar(0,length(nodes),style=3)
for(i in 1:length(nodes)){
  met_genes <- names(network[[nodes[i]]][["tfmode"]])
  met_tfmode <- network[[nodes[i]]][["tfmode"]]
  met_lik <- network[[nodes[i]]][["likelihood"]]
  met<-rep(nodes[i],length(met_genes))
  new_row<-cbind(met,met_genes, met_tfmode, met_lik)
  net_matrix<-rbind(net_matrix,new_row)
  setTxtProgressBar(pb,i)
}

# Removing quotes due to problems in Cytoscape
net_matrix <- noquote(net_matrix)
write.table(net_matrix, file="tables_cyt/net_matrix.txt", row.names=FALSE, col.names=TRUE, quote = FALSE, sep = "\t")

# Creating mapping table used to differentiate metabolites from genes on Cytoscape
type <- rep("node", length(nodes))
nodes_matrix <- cbind(nodes,type)
write.table(nodes_matrix, file="tables_cyt/nodes_matrix.txt", row.names=FALSE, col.names=TRUE, quote = FALSE, sep = "\t")

# Creating network for the candidate uracil
met_lik_uracil <- network[["uracil"]][["likelihood"]]
met_genes_uracil <- names(network[["uracil"]][["tfmode"]])
met_tfmode_uracil <- network[["uracil"]][["tfmode"]]
met_uracil <- rep("uracil", length(met_genes_uracil))
uracil_matrix <- cbind(met_uracil, met_genes_uracil, met_tfmode_uracil, met_lik_uracil)
write.table(uracil_matrix, file=paste0("tables_cyt/uracil_matrix.txt"), row.names=FALSE, col.names=TRUE, quote = FALSE, sep = "\t")

# mapping tfmode value to customize size of nodes
uracil_nodes_size <- cbind(names(network[["uracil"]][["tfmode"]]), abs(network[["uracil"]][["tfmode"]]))
write.table(uracil_nodes_size, file=paste0("tables_cyt/uracil_nodes_size.txt"), row.names=FALSE, col.names=TRUE, quote = FALSE, sep = "\t")


# extracting negative and positive tfmode of uracil for gene ontology
pos_uracil <- c()
neg_uracil <- c()
for(i in names(network[["uracil"]][["tfmode"]])){
  if(network[["uracil"]][["tfmode"]][[i]] > 0){
    pos_uracil <- c(pos_uracil, i)
  }
  else{
    neg_uracil <- c(neg_uracil, i)
  }
}
uracil_ordered <- c(pos_uracil, neg_uracil)
save(uracil_ordered,file="results/GO/uracil_ordered.rda")

write.table(pos_uracil, file=paste0("results/GO/pos_uracil.txt"), row.names=FALSE, col.names=FALSE, quote = FALSE, sep = "\t")
write.table(neg_uracil, file=paste0("results/GO/neg_uracil.txt"), row.names=FALSE, col.names=FALSE, quote = FALSE, sep = "\t")
write.table(genes, file=paste0("results/GO/all_net_genes.txt"), row.names=FALSE, col.names=FALSE, quote = FALSE, sep = "\t")


# GO: EXTRACTING THE SETS OF GENES CONNECTED TO THE METABOLITES ####
#  (IN THE CCLE-BASED NETWORK) THAT WERE PREDICTED IN NCI60 

for(i in common_NCI60){
  plus <- 0
  minus <- 0
  go_table <- c("GENE", "SIGN")
  for(j in names(network[[i]][["tfmode"]])){
    if(network[[i]][["tfmode"]][[j]] > 0){
      go_table <- rbind(go_table, c(j, "+"))
      plus <- plus +1
    }
    else{
      go_table <- rbind(go_table, c(j, "-"))
      minus <- minus+1
    }
  }
  colnames(go_table) <-go_table[1,]
  go_table <- go_table[-1,]
  print(i)
  print(plus)
  print(minus)
  
  save(go_table, file =  paste0("data/final_ds/GO/",i,"_go.rda"))
}



