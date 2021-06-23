setwd("/home/student/Desktop/vittoria")
library(corto)

correlation <- function(real, pred, thr){
  cors_spearman <- c()
  cors_pearson <- c()
  significant <- c()
  for(i in rownames(real)){
    cors_spearman <- c(cors_spearman,cor(real[i,],pred[i,],method="spearman"))
    cors_pearson <- c(cors_pearson,cor(real[i,],pred[i,],method="pearson"))
    
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
  
  return(list(cors_spearman, cors_pearson, significant))
}

# Load CCLE metabolite data
load("data/datasets/ccle/raw/ccle-metab.rda")
metabolites <- rownames(metab)
ncol(metab) # 928 samples (cell lines)
rm(metab)
length(metabolites) # 225 measured metabolites


# Load expression + metabolite data (expression data is VST normalized)
load("data/datasets/ccle/raw/ccle-metabexp_mat.rda")
dim(mat) # 24499 molecular species (225 metabolites, 24274 transcripts) in 898 common samples


# Check if one half can predict the other half 
i <- 1 
set.seed(i) # Reproducible analysis
s1 <- sample(1:ncol(mat), floor(ncol(mat)/2))
s2 <- setdiff(1:ncol(mat),s1)
half1 <- mat[,s1]
half2 <- mat[,s2]
# Generate a metabolite/transcript network
network <- corto(half1,centroids=metabolites, nbootstraps=100, p=1e-10, nthreads=7, verbose=FALSE)

# Use the network to predict metabolite levels in the second half
# We need to give the function a matrix with genes only --> "setdiff" will remove all metabolites rows 
predscore <- mra(half2[setdiff(rownames(half2),metabolites),],regulon=network,verbose=FALSE)
realscore <- half2[metabolites,]
common<-intersect(rownames(predscore),rownames(realscore))

# Check how well each metabolite is predicted
cc_thr <- 0.3
corr_outputs <- correlation(realscore[common,], predscore[common,], cc_thr)
cors_spearman <- corr_outputs[[1]]
cors_pearson <- corr_outputs[[2]]

names(cors_spearman) <- common
names(cors_pearson) <- common

significant <- corr_outputs[[3]]

length(common) # among the 98 correctly predicted metabolites...
length(significant) # 75 are also significant

# Testing distribution ####
# Testing how much dist is diff from zero
shapiro.test(cors_spearman) #p-value = 0.6024 --> can use t-test
shapiro.test(cors_pearson) #p-value = 0.5932 --> can use t-test
t1 <- t.test(cors_spearman, alternative = "greater") #p-value = 2.2e-16
t2 <- t.test(cors_pearson, alternative = "greater") # p-value = 2.2e-16

png(paste0("results/corr_distributions/CCLE_train_test.png"),w=3000,h=2500,res=600)
plot(density(cors_spearman),col="darkgoldenrod3",lwd=3, xlim=c(-1,1), ylim = c(0,max(density(cors_spearman)$y)), 
     main = "CCLE prediction performance", xlab = "metabolite CC")
lines(density(cors_pearson), col = "brown3", lty = 3, lwd = 2)
abline(v=0,lty=2, col="gray")
mtext(paste0("Spearman p-value = ",round(t1$p.value,3),"   Pearson p-value = ",round(t2$p.value,3)))
legend("topleft", legend = c("Spearman", "Pearson"), lwd = 3,
       col = c("darkgoldenrod3", "brown3"), lty = c(1,2))
dev.off()

# 1.6 Scatterplot of each single correlation ####
c <- 1
for(i in common){
  png(paste0("results/test_CCLE/test_",c,".png"),w=2900,h=2500,res=600)
  scatter(predscore[i,],realscore[i,], xlab = paste0("measured_",i), ylab = paste0("pred_",i), method = "spearman", 
          main = paste0(toupper(i), "\nmeasured vs predicted in CCLE"))
  dev.off()
  c <- c+1
}

# Plot overall performance on different halves of the CCLE dataset

# results <- list(ps_s = c(), ps_p = c(), sig = c())
# for(i in 1:100){
#   set.seed(i)
#   print(i)
#   # generation of the 2 halves
#   s1 <- sample(1:ncol(mat), floor(ncol(mat)/2))
#   s2 <- setdiff(1:ncol(mat),s1)
#   half1 <- mat[,s1]
#   half2 <- mat[,s2]
# 
#   # train network
#   #network <- corto(half1,centroids=metabolites, nbootstraps=100, p=1e-10, nthreads=7, verbose=FALSE)
#   #save(network, file=paste0("data/ccle_halves/ccle_half_",i,".rda"))
#   load(paste0("data/ccle_halves/ccle_half_",i,".rda"))
# 
#   # prediction
#   predscore <- mra(half2[setdiff(rownames(half2),metabolites),],regulon=network,verbose=FALSE)
#   realscore <- half2[metabolites,]
#   common<-intersect(rownames(predscore),rownames(realscore))
# 
#   # evaluate performance
#   cc_thr <- 0.3
#   corr_outputs <- correlation(realscore[common,], predscore[common,], cc_thr)
#   current_spearman <- corr_outputs[[1]]
#   current_pearson <- corr_outputs[[2]]
# 
#   names(current_spearman) <- common
#   names(current_pearson) <- common
# 
#   current_significant <- length(corr_outputs[[3]])/length(common)
# 
#   # Testing how much dist is diff from zero
#   if(shapiro.test(current_spearman)$p.value > 0.05){
#     t1 <- t.test(current_spearman, alternative = "greater")
#   }else{
#     t1 <- wilcox.test(current_spearman, alternative = "greater")
#   }
# 
#   if(shapiro.test(current_pearson)$p.value > 0.05){
#     t2 <- t.test(current_pearson, alternative = "greater")
#   }else{
#     t2 <- wilcox.test(current_pearson, alternative = "greater")
#   }
# 
#   results$ps_s <- c(results$ps_s, t1$p.value)
#   results$ps_p <- c(results$ps_p, t2$p.value)
#   results$sig <- c(results$sig, current_significant)
# 
# }
# save(results, file= "data/ccle_halves/results.rda")
load("data/ccle_halves/results.rda")
ps  <- list(Spearman = log10(results$ps_s), Pearson = log10(results$ps_p))

png(paste0("results/figures/3_prediction_performace/ccle_halves_boxplot.png"),w=2900,h=2500,res=600)
boxplot(ps,xlab= c("Correlation Method"), ylab="log10(pvalue)",main="CCLE performance")
mtext("Metabolite Prediction, Train: half1 CCLE, Test: half2 CCLE",cex=0.9,font=2)
dev.off()

png(paste0("results/figures/3_prediction_performace/significants_ccle_halves.png"),w=2900,h=2500,res=600)
plot(density(results$sig*100), col="coral2",lwd=3, 
     main = "CCLE predictions CC > 0.3", xlab = "% metabolites with CC > 0.3")
abline(v = mean(results$sig*100), lty = 2)
dev.off()



