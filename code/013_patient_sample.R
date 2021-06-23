# Last chance for Siddiqui_breast and Hassan to show some significant result: failed

setwd("/home/student/Desktop/vittoria")

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
    # if(cors_pearson[length(cors_pearson)] > thr){
    #   significant <- c(significant, i)
    # }
    #if(cors_spearman[length(cors_spearman)] > thr && cors_pearson[length(cors_pearson)] > thr){
    #  significant <- c(significant, i)
    #}
  }
  
  return(list(cors_spearman, cors_pearson, significant))
}

library(corto)

load("data/final_ds/genes_breast.rda")
load("data/final_ds/metab_breast.rda")
names_metab_breast <- rownames(metab_breast)
load("data/final_ds/genes_hassan.rda")
load("data/final_ds/metab_hassan.rda")


real_names_hassan2 <- as.matrix(read.delim("data/datasets/hassan/changes/changed_names2.txt", sep="\t", header = F))
rownames(metab_hassan) <- real_names_hassan2

length(intersect(rownames(metab_breast), rownames(metab_hassan)))

data_for_network_breast <- rbind(genes_breast, metab_breast)
network_breast <- corto(data_for_network_breast,centroids=names_metab_breast, nbootstraps=100, p=1e-10, nthreads=7, verbose=FALSE)

real_names_hassan <- rownames(metab_hassan)
pred_hassan_breast <- mra(genes_hassan, regulon = network_breast)
pred_metab_names_hassan <- rownames(pred_hassan_breast)
common_hassan_pred <- intersect(pred_metab_names_hassan, real_names_hassan)

length(common_hassan_pred) # 11 correctly predicted metabolites, however...

# Correlation checks 

cc_thr <- 0.3
corr_outputs_hassan <- correlation(metab_hassan[common_hassan_pred,], pred_hassan_breast[common_hassan_pred,], cc_thr)
cors_hassan_spearman <- corr_outputs_hassan[[1]]
cors_hassan_pearson <- corr_outputs_hassan[[2]]

names(cors_hassan_spearman) <- common_hassan_pred
names(cors_hassan_pearson) <- common_hassan_pred

significant_hassan<- corr_outputs_hassan[[3]]

length(significant_hassan) # just 1/11 is significant (Sp), 4/11 (Pe)

# Testing how much dist is diff from zero
shapiro.test(cors_hassan_spearman) #p-value = 0.687 --> can use t-test
shapiro.test(cors_hassan_pearson) #p-value = 0.6912 --> can t-test
t1_hassan <- t.test(cors_hassan_spearman, alternative = "greater") #p-value = 0.07298849
t2_hassan <- t.test(cors_hassan_pearson, alternative = "greater") # p-value = 0.007549814

# Plotting correlation distributions
png(paste0("results/corr_distributions/correlations_trBreast_teHassan.png"),w=2900,h=2500,res=600)
plot(density(cors_hassan_spearman),col="orange red",lwd=3, xlim=c(-1,1), ylim = c(0,max(density(cors_hassan_spearman)$y)), 
     main = "GSE148892 correlations distribution", xlab = "metabolite CC")
lines(density(cors_hassan_pearson), lty = 2, lwd = 2)
abline(v=0,lty=2, col="gray")
mtext(paste0("Spearman p-value = ",round(t1_hassan$p.value,3)," Pearson p-value = ",round(t2_hassan$p.value,3)))
legend("topleft", legend = c("Spearman", "Pearson"), lwd = c(3,2),
       col = c("orange red", "black"), lty = c(1,2))
dev.off()

#scatter(metab_hassan["serine",], pred_hassan_breast["serine",])



