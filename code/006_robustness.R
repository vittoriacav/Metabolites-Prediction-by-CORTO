# Testing for robustness of the CCLE predicting network
# I'll compare the performance of the prediction on NCI60 between a the real CCLE-based network a
#  and a random-genes CCLE-network
setwd("C:/Users/vitto/OneDrive - Alma Mater Studiorum Università di Bologna/UNIVERSITA'/Thesis/vittoria")
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

# Loading CCLE data
load("data/final_ds/metab_ccle.rda")
names_metab_ccle <- rownames(metab_ccle)
load("data/final_ds/genes_ccle.rda")
names_gene_ccle <- rownames(genes_ccle)

# Loading the real network, previously trained in "003_main.R"
load("code/network.rda")

# Loading NCI60 dataset in order to test of the prediction perform using the 2 different networks
load("data/final_ds/genes_NCI60.rda")
load("data/final_ds/metab_NCI60.rda")

# Prediction using real network (genes correctly "placed")
real_metab_names_NCI60 <- rownames(metab_NCI60)
pred <- mra(genes_NCI60, regulon = network)
pred_metab_names_NCI60<- rownames(pred)
common <- intersect(pred_metab_names_NCI60, real_metab_names_NCI60)
length(common) # 21 

cc_thr <- 0.3
corr_outputs <- correlation(metab_NCI60[common,], pred[common,], cc_thr)
cors_spearman <- corr_outputs[[1]]
cors_pearson <- corr_outputs[[2]]

names(cors_spearman) <- common
names(cors_pearson) <- common

significant <- corr_outputs[[3]]
length(significant) #10/21 significant

png(paste0("results/figures/3_prediction_performace/ROBUSTNESS_NCI60.png"),w=3000,h=2500,res=600)
plot(density(cors_perm_spearman),col="darkgoldenrod3",lwd=3, xlim=c(-1,1), ylim = c(0,3), 
     main = "NCI60 correlations ditribution\n(with random network)", xlab = "metabolite CC")
lines(density(cors_perm_pearson), col = "brown3", lty = 3, lwd = 2)
abline(v=0,lty=2, col="gray")
mtext(paste0("Spearman p-value = ",round(t1_perm$p.value,3),"     Pearson p-value = ",round(t2_perm$p.value,3)))
legend("topleft", legend = c("Spearman", "Pearson"), lwd = 2,
       col = c("darkgoldenrod3", "brown3"), lty = c(1,2))
dev.off()

# T-test
t1 <- t.test(cors_spearman, cors_perm_spearman, alternative = "greater") # p-value = 3.848124e-05
t2 <- t.test(cors_pearson, cors_perm_pearson, alternative = "greater") # p-value = 1.013e-05

png(paste0("results/figures/3_prediction_performace/ROBUSTNESS_boxplot.png"),w=4000,h=2700,res=600)
data <- cbind(cors_spearman,cors_perm_spearman)
ylim=c(min(data), max(data))

layout(matrix(1:2, nc=2), widths=c(9,9))
par(las=1, mar=c(4,4,5,2))
boxplot(cors_spearman, cors_perm_spearman, col=c("darkgoldenrod4", "darkgoldenrod1"),
        names = c("REAL", "PERMUTED"), ylab = "measured vs predicted metabolite CC")
title(main = "Spearman Correlation", line = 1.4)
mtext(paste0("p-value = ",round(t1$p.value,3)))

par(mar=c(4,4,5,2))
boxplot(cors_pearson, cors_perm_pearson, col=c("brown4", "brown1"), names = c("REAL", "PERMUTED"),
        ylab = "measured vs predicted metabolite CC")
title(main = "Pearson Correlation",  line = 1.4)
mtext(paste0("p-value = ",round(t2$p.value,3)))

title(main="Correlations boxplots: real vs permuted network", line = -1.5, outer=TRUE, cex.main=1.5)
dev.off()


# Shuffle gene names completely random
results <- list(c_sp = list(), c_pe=list())
for(i in 1:500){
  print(i)
  tester_genes_ccle <- genes_ccle
  set.seed(i)
  rownames(tester_genes_ccle) <- sample(rownames(genes_ccle))
  tester_mat <- rbind(tester_genes_ccle, metab_ccle)
  tester_network <- corto(tester_mat,centroids=names_metab_ccle, nbootstraps=100, p=1e-10, nthreads=7, verbose=FALSE)
  save(tester_network, file=paste0("data/robustness_ds/tester_network_",i,".rda"))
  load(paste0("data/robustness_ds/tester_network_",i,".rda"))
  #predict on NCI60 using this random network
  pred_perm <- mra(genes_NCI60, regulon = tester_network)
  pred_metab_names_perm <- rownames(pred_perm)
  common_perm <- intersect(pred_metab_names_perm, real_metab_names_NCI60) 
  corr_outputs_perm <- correlation(metab_NCI60[common_perm,], pred_perm[common_perm,], cc_thr)
  cors_perm_spearman <- corr_outputs_perm[[1]]
  cors_perm_pearson <- corr_outputs_perm[[2]]
  
  names(cors_perm_spearman) <- common_perm
  results$c_sp[[as.character(i)]] <- cors_perm_spearman
  names(cors_perm_pearson) <- common_perm
  results$c_pe[[as.character(i)]] <- cors_perm_pearson
  setTxtProgressBar(pb,i)
  
}
save(results, file ="data/robustness_ds/results.rda")

ps <- list(Spearman = c(), Pearson = c())

for(i in 1:500){
  # Testing how much dist is diff from zero
  if(shapiro.test(results$c_sp[[i]])$p.value > 0.05){
    t1 <- t.test(results$c_sp[[i]], alternative = "greater")
  }else{
    t1 <- wilcox.test(results$c_sp[[i]], alternative = "greater")
  }
  
  if(shapiro.test(results$c_pe[[i]])$p.value > 0.05){
    t2 <- t.test(results$c_pe[[i]], alternative = "greater")
  }else{
    t2 <- wilcox.test(results$c_pe[[i]], alternative = "greater")
  }
  ps$Spearman <- c(ps$Spearman, t1$p.value)
  ps$Pearson <- c(ps$Pearson, t2$p.value)
}


png(paste0("results/figures/3_prediction_performace/robustness_ps_dist.png"),w=2900,h=2500,res=600)
boxplot(ps,xlab="method",ylab="pvalue",main="Permuted Networks Performance")
mtext("Metabolite Prediction, Train: permuted CCLE, Test: NCI-60",cex=0.9,font=2)
dev.off()











