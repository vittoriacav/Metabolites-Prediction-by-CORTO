setwd("C:/Users/vitto/OneDrive - Alma Mater Studiorum Università di Bologna/UNIVERSITA'/Thesis/vittoria")

library(randomForest)
library(glmnet)
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

settings <- function(genes1, genes2, metab1, metab2){
  settings <- list()
  settings$com_met <- intersect(rownames(metab1), rownames(metab2))
  com_genes <- intersect(rownames(genes1), rownames(genes2))
  predictors <- apply(genes1[com_genes,],1,var)
  predictors <- names(sort(predictors,dec=TRUE))[1:1000]
  
  settings$preds <- predictors
  settings$new_genes1 <- genes1[predictors,]
  
  temp_new_genes1 <- rownames(settings$new_genes1) 
  temp_new_genes1 <- gsub("-", "_", temp_new_genes1) 
  rownames(settings$new_genes1) <- temp_new_genes1
  
  settings$new_genes2 <- genes2[predictors,]
  temp_new_genes2 <- rownames(settings$new_genes2) 
  temp_new_genes2 <- gsub("-", "_", temp_new_genes2) 
  rownames(settings$new_genes2) <- temp_new_genes2
  
  return(settings)
}

train_pred <- function(common_met, met1, genes1, genes2){
  pred_mat <- c()
  pb = txtProgressBar(min = 0, max = length(common_met), initial = 0, style = 3)
  for(met in common_met){
    Tmet1 <- t(met1)
    Tgenes1 <- t(genes1)
    metabolite <- Tmet1[,met]
    
    # create a temporary variable which stores the response and the significant predictors only
    current <- cbind(metabolite, Tgenes1)
    
    # Train the random forest
    rf <- randomForest(metabolite ~ ., data = current, ntree = 1000)
    #print(rf)
    #plot(rf)
    
    # Predict
    pred <- predict(rf, t(genes2))
    
    # Fill the prediction matrix
    if(length(pred_mat) == 0){
      pred_mat <- pred
    }else{
      pred_mat <- rbind(pred_mat, pred)
    }
    setTxtProgressBar(pb,match(met,common_met))
  }
  close(pb)
  rownames(pred_mat) <- common_met
  return(pred_mat)
}



#Loading CCLE data: I'll always use CCLE data set to train my predictive tool
load("data/final_ds/genes_ccle.rda")
load("data/final_ds/metab_ccle.rda")




# NCI60 ####

load("data/final_ds/genes_NCI60.rda")
load("data/final_ds/metab_NCI60.rda")


settings <- settings(genes_ccle, genes_NCI60, metab_ccle, metab_NCI60)
common_met_NCI60_ccle <- settings$com_met
rf_genes_ccle <- settings$new_genes1
rf_genes_NCI60 <- settings$new_genes2
#predictors <- settings$preds

pred_rf_NCI60 <- train_pred(common_met_NCI60_ccle, metab_ccle, rf_genes_ccle, rf_genes_NCI60)

save(pred_rf_NCI60, file="data/final_ds/rf/pred_rf_NCI60.rda")
load("data/final_ds/rf/pred_rf_NCI60.rda")

cc_thr <- 0.3
corr_outputs_NCI60 <- correlation(metab_NCI60[common_met_NCI60_ccle,], pred_rf_NCI60[common_met_NCI60_ccle,], cc_thr)
cors_NCI60_spearman <- corr_outputs_NCI60[[1]]
cors_NCI60_pearson <- corr_outputs_NCI60[[2]]

names(cors_NCI60_spearman) <- common_met_NCI60_ccle
names(cors_NCI60_pearson) <- common_met_NCI60_ccle

significant_NCI60 <- corr_outputs_NCI60[[3]]

length(common_met_NCI60_ccle) # Among 35 correctly predicted
length(significant_NCI60) #  just 5/35 are significant

# Saving with different names to better manage the step in which I will compare the correlations performances of 
#  different tools
NCI60_sp_rf <- cors_NCI60_spearman
NCI60_pe_rf <- cors_NCI60_pearson

save(NCI60_sp_rf, file = "data/final_ds/corr/NCI60_sp_rf.rda")
save(NCI60_pe_rf, file = "data/final_ds/corr/NCI60_pe_rf.rda")


# Testing how much dist is diff from zero
shapiro.test(cors_NCI60_spearman) #p-value = 0.08919 --> cannot use t-test
shapiro.test(cors_NCI60_pearson) #p-value = 0.1629 --> can use t-test
t1_NCI60 <- t.test(cors_NCI60_spearman, alternative = "greater") #p-value = 0.1238
t2_NCI60 <- t.test(cors_NCI60_pearson, alternative = "greater") # p-value = 0.08029


# For thesis
png(paste0("results/figures/8_rf_pred/correlations_NCI60.png"),w=2900,h=2500,res=600)
plot(density(cors_NCI60_spearman),col="darkgoldenrod3",lwd=3, xlim=c(-1,1), ylim = c(0,max(density(cors_NCI60_spearman)$y)), 
     main = "NCI60 correlations ditribution", xlab = "metabolite CC")
lines(density(cors_NCI60_pearson), col = "brown3", lty = 3, lwd = 2)
abline(v=0,lty=2, col="gray")
mtext(paste0("Spearman p-value = ",round(t1_NCI60$p.value,3),"     Pearson p-value = ",round(t2_NCI60$p.value,3)))
legend("topleft", legend = c("Spearman", "Pearson"), lwd = 2,
       col = c("darkgoldenrod3", "brown3"), lty = c(1,2))
dev.off()




# BREAST ####

load("data/final_ds/genes_breast.rda")
load("data/final_ds/metab_breast.rda")


settings <- settings(genes_ccle, genes_breast, metab_ccle, metab_breast)
common_met_breast_ccle <- settings$com_met
rf_genes_ccle <- settings$new_genes1
rf_genes_breast <- settings$new_genes2
#predictors <- settings$preds

pred_rf_breast <- train_pred(common_met_breast_ccle, metab_ccle, rf_genes_ccle, rf_genes_breast)

save(pred_rf_breast, file="data/final_ds/rf/pred_rf_breast.rda")
load("data/final_ds/rf/pred_rf_breast.rda")

cc_thr <- 0.3
corr_outputs_breast <- correlation(metab_breast[common_met_breast_ccle,], pred_rf_breast[common_met_breast_ccle,], cc_thr)
cors_breast_spearman <- corr_outputs_breast[[1]]
cors_breast_pearson <- corr_outputs_breast[[2]]

names(cors_breast_spearman) <- common_met_breast_ccle
names(cors_breast_pearson) <- common_met_breast_ccle

significant_breast <- corr_outputs_breast[[3]]

length(common_met_breast_ccle) # Among 70 correctly predicted metabolites...
length(significant_breast) # just 3/70 are significant

# Saving with different names to better manage the step in which I will compare the correlations performances of 
#  different tools
breast_sp_rf <- cors_breast_spearman
breast_pe_rf <- cors_breast_pearson

save(breast_sp_rf, file = "data/final_ds/corr/breast_sp_rf.rda")
save(breast_pe_rf, file = "data/final_ds/corr/breast_pe_rf.rda")


# Testing how much dist is diff from zero
shapiro.test(cors_breast_spearman) #p-value = 0.661 --> cannot use t-test
shapiro.test(cors_breast_pearson) #p-value = 0.809 --> can use t-test
t1_breast <- t.test(cors_breast_spearman, alternative = "greater") #p-value = 0.9792
t2_breast <- t.test(cors_breast_pearson, alternative = "greater") # p-value = 0.9838


# For thesis
png(paste0("results/figures/8_rf_pred/correlations_breast.png"),w=2900,h=2500,res=600)
plot(density(cors_breast_spearman),col="darkgoldenrod3",lwd=3, xlim=c(-1,1), ylim = c(0,max(density(cors_breast_spearman)$y)), 
     main = "GSE37751 correlations ditribution", xlab = "metabolite CC")
lines(density(cors_breast_pearson), col = "brown3", lty = 3, lwd = 2)
abline(v=0,lty=2, col="gray")
mtext(paste0("Spearman p-value = ",round(t1_breast$p.value,3),"     Pearson p-value = ",round(t2_breast$p.value,3)))
legend("topright", legend = c("Spearman", "Pearson"), lwd = 2,
       col = c("darkgoldenrod3", "brown3"), lty = c(1,2))
dev.off()





# HASSAN ####

load("data/final_ds/genes_hassan.rda")
load("data/final_ds/metab_hassan.rda")


settings <- settings(genes_ccle, genes_hassan, metab_ccle, metab_hassan)
common_met_hassan_ccle <- settings$com_met
rf_genes_ccle <- settings$new_genes1
rf_genes_hassan <- settings$new_genes2
#predictors <- settings$preds

pred_rf_hassan <- train_pred(common_met_hassan_ccle, metab_ccle, rf_genes_ccle, rf_genes_hassan)

save(pred_rf_hassan, file="data/final_ds/rf/pred_rf_hassan.rda")
load("data/final_ds/rf/pred_rf_hassan.rda")

cc_thr <- 0.3
corr_outputs_hassan <- correlation(metab_hassan[common_met_hassan_ccle,], pred_rf_hassan[common_met_hassan_ccle,], cc_thr)
cors_hassan_spearman <- corr_outputs_hassan[[1]]
cors_hassan_pearson <- corr_outputs_hassan[[2]]

names(cors_hassan_spearman) <- common_met_hassan_ccle
names(cors_hassan_pearson) <- common_met_hassan_ccle

significant_hassan <- corr_outputs_hassan[[3]]

length(common_met_hassan_ccle) # Among 5 correctly predicted metabolites...
length(significant_hassan) # ... just 3/35 are significant

# Saving with different names to better manage the step in which I will compare the correlations performances of 
#  different tools
hassan_sp_rf <- cors_hassan_spearman
hassan_pe_rf <- cors_hassan_pearson

save(hassan_sp_rf, file = "data/final_ds/corr/hassan_sp_rf.rda")
save(hassan_pe_rf, file = "data/final_ds/corr/hassan_pe_rf.rda")


# Testing how much dist is diff from zero
shapiro.test(cors_hassan_spearman) #p-value = 0.1418 --> cannot use t-test
shapiro.test(cors_hassan_pearson) #p-value = 0.1021 --> can use t-test
t1_hassan <- t.test(cors_hassan_spearman, alternative = "greater") #p-value = 0.1545
t2_hassan <- t.test(cors_hassan_pearson, alternative = "greater") # p-value = 0.1974


# For thesis
png(paste0("results/figures/8_rf_pred/correlations_hassan.png"),w=2900,h=2500,res=600)
plot(density(cors_hassan_spearman),col="darkgoldenrod3",lwd=3, xlim=c(-1,1), ylim = c(0,max(density(cors_hassan_pearson)$y)), 
     main = "GSE148892 correlations ditribution", xlab = "metabolite CC")
lines(density(cors_hassan_pearson), col = "brown3", lty = 3, lwd = 2)
abline(v=0,lty=2, col="gray")
mtext(paste0("Spearman p-value = ",round(t1_hassan$p.value,3),"     Pearson p-value = ",round(t2_hassan$p.value,3)))
legend("topleft", legend = c("Spearman", "Pearson"), lwd = 2,
       col = c("darkgoldenrod3", "brown3"), lty = c(1,2))
dev.off()




# ZAMPIERI ####

load("data/final_ds/genes_zamp.rda")
load("data/final_ds/metab_zamp.rda")


settings <- settings(genes_ccle, genes_zamp, metab_ccle, metab_zamp)
common_met_zamp_ccle <- settings$com_met
rf_genes_ccle <- settings$new_genes1
rf_genes_zamp <- settings$new_genes2
#predictors <- settings$preds

pred_rf_zamp <- train_pred(common_met_zamp_ccle, metab_ccle, rf_genes_ccle, rf_genes_zamp)

save(pred_rf_zamp, file="data/final_ds/rf/pred_rf_zamp.rda")
load("data/final_ds/rf/pred_rf_zamp.rda")

cc_thr <- 0.3
cors_zamp_spearman<-c()
cors_zamp_pearson<-c()
significant_zamp <- c()
for(j in common_met_zamp_ccle){
  temp_real <- c()
  temp_pred <- c()
  for(v in 1:length(metab_zamp[j,])){
    if(metab_zamp[j,v] != "NaN"){
      temp_real <- c(temp_real, metab_zamp[j,v])
      temp_pred <- c(temp_pred, pred_rf_zamp[j,v])
    }
  }
  cors_zamp_spearman<-c(cors_zamp_spearman,cor(temp_real,temp_pred,method="spearman"))
  cors_zamp_pearson<-c(cors_zamp_pearson,cor(temp_real,temp_pred,method="pearson"))
  if(cors_zamp_spearman[length(cors_zamp_spearman)]>cc_thr && cors_zamp_pearson[length(cors_zamp_pearson)]>cc_thr){
    significant_zamp <- c(significant_zamp, j)
  }
}
names(cors_zamp_spearman) <- common_met_zamp_ccle
names(cors_zamp_pearson) <- common_met_zamp_ccle

length(common_met_zamp_ccle) # Among 48 correctly predicted metabolites...
length(significant_zamp) # ... just 1/48 are significant (> 0.3)

# Saving with different names to better manage the step in which I will compare the correlations performances of 
#  different tools
zamp_sp_rf <- cors_zamp_spearman
zamp_pe_rf <- cors_zamp_pearson

save(zamp_sp_rf, file = "data/final_ds/corr/zamp_sp_rf.rda")
save(zamp_pe_rf, file = "data/final_ds/corr/zamp_pe_rf.rda")


# Testing how much dist is diff from zero
shapiro.test(cors_zamp_spearman) #p-value = 0.05493 --> can use t-test
shapiro.test(cors_zamp_pearson) #p-value = 0.02533 --> cannot use t-test
t1_zamp <- t.test(cors_zamp_spearman, alternative = "greater") #p-value = 0.7674
t2_zamp <- wilcox.test(cors_zamp_pearson, alternative = "greater") # p-value = 0.4656


# For thesis
png(paste0("results/figures/8_rf_pred/correlations_zamp.png"),w=2900,h=2500,res=600)
plot(density(cors_zamp_spearman),col="darkgoldenrod3",lwd=3, xlim=c(-1,1), ylim = c(0,max(density(cors_zamp_spearman)$y)), 
     main = "GSE32474 correlations ditribution", xlab = "metabolite CC")
lines(density(cors_zamp_pearson), col = "brown3", lty = 3, lwd = 2)
abline(v=0,lty=2, col="gray")
mtext(paste0("Spearman p-value = ",round(t1_zamp$p.value,3),"     Pearson p-value = ",round(t2_zamp$p.value,3)))
legend("topleft", legend = c("Spearman", "Pearson"), lwd = 2,
       col = c("darkgoldenrod3", "brown3"), lty = c(1,2))
dev.off()

