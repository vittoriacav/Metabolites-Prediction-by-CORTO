setwd("C:/Users/vitto/OneDrive - Alma Mater Studiorum Università di Bologna/UNIVERSITA'/Thesis/vittoria")
library(glmnet)

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

# LASSO ####

#Loading CCLE data: I'll always use CCLE data set to train my predictive tool
load("data/final_ds/genes_ccle.rda")
load("data/final_ds/metab_ccle.rda")

# NCI60 ####
# trying linear model (with lasso) on NCI60
load("data/final_ds/genes_NCI60.rda")
load("data/final_ds/metab_NCI60.rda")

# First I need to select common metabolites and common genes between CCLE(train) and NCI60 (test)
#  because this tool must have exactly same predictors as input of the testing dataset
# I also select common metabolites otherwise later I would have to do unnecessary and computationally expensive work
common_met_NCI60_ccle <- intersect(rownames(metab_ccle), rownames(metab_NCI60))
common_genes_NCI60_ccle <- intersect(rownames(genes_ccle), rownames(genes_NCI60))

# Deleting non-common genes
new_genes_NCI60<-genes_NCI60[common_genes_NCI60_ccle,]
new_genes_ccle <- genes_ccle[common_genes_NCI60_ccle,]

# Creating a empty matrix that will host the metabolites prediction of NCI60 
pred_metab_mat_NCI60 <- matrix(ncol = length(colnames(genes_NCI60)))
colnames(pred_metab_mat_NCI60) <- colnames(genes_NCI60)

pb = txtProgressBar(min = 0, max = length(common_met_NCI60_ccle), initial = 0, style = 3)
for(i in common_met_NCI60_ccle){
  metabolite <- metab_ccle[i,]
  model<-cv.glmnet(t(new_genes_ccle),metabolite,alpha=1)
  #plot(model)

  # Get the best model
  bestcv<-which(model$lambda==model$lambda.min)
  betas<-model$glmnet.fit$beta[,bestcv]
  intercept<-model$glmnet.fit$a0[bestcv]
  sigvars<-names(betas)[betas!=0]
  bestmodel<-betas[sigvars]

  # R^2
  #model$glmnet.fit$dev.ratio[bestcv]

  # Predict
  predvalues<-as.numeric(predict(model,newx=t(new_genes_NCI60),s = model$lambda.min))
  #plot(values,predvalues)
  pred_metab_mat_NCI60 <- rbind(pred_metab_mat_NCI60, predvalues)
  setTxtProgressBar(pb,match(i,common_met_NCI60_ccle))

}
close(pb)

pred_metab_mat_NCI60 <- pred_metab_mat_NCI60[-1,]
rownames(pred_metab_mat_NCI60) <- common_met_NCI60_ccle
save(pred_metab_mat_NCI60, file="data/final_ds/lm/pred_metab_mat_NCI60.rda")
load("data/final_ds/lm/pred_metab_mat_NCI60.rda")


cc_thr <- 0.3
corr_outputs_NCI60 <- correlation(metab_NCI60[common_met_NCI60_ccle,], pred_metab_mat_NCI60[common_met_NCI60_ccle,], cc_thr)
cors_NCI60_spearman <- corr_outputs_NCI60[[1]]
cors_NCI60_pearson <- corr_outputs_NCI60[[2]]

names(cors_NCI60_spearman) <- common_met_NCI60_ccle
names(cors_NCI60_pearson) <- common_met_NCI60_ccle

significant_NCI60 <- corr_outputs_NCI60[[3]]
length(common_met_NCI60_ccle) # Among 35 common metabolites predicted...
length(significant_NCI60) # just 6/35 predictions are significant

# Saving with different names to better manage the step in which I will compare the correlations performances of 
#  different tools
NCI60_sp_lm <- cors_NCI60_spearman
NCI60_pe_lm <- cors_NCI60_pearson

save(NCI60_sp_lm, file = "data/final_ds/corr/NCI60_sp_lm.rda")
save(NCI60_pe_lm, file = "data/final_ds/corr/NCI60_pe_lm.rda")

# Testing how much dist is diff from zero
shapiro.test(cors_NCI60_spearman) #p-value = 0.1733 --> can use t-test
shapiro.test(cors_NCI60_pearson) #p-value = 0.0782 --> can use t-test
t1_NCI60 <- t.test(cors_NCI60_spearman, alternative = "greater") #p-value = 5.784e-07
t2_NCI60 <- t.test(cors_NCI60_pearson, alternative = "greater") # p-value = 5.404e-07


# For thesis
png(paste0("results/figures/7_lasso_pred/correlations_NCI60.png"),w=2900,h=2500,res=600)
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
common_met_breast_ccle <- intersect(rownames(metab_ccle), rownames(metab_breast))
common_genes_breast_ccle <- intersect(rownames(genes_ccle), rownames(genes_breast))

new_genes_breast<-genes_breast[common_genes_breast_ccle,]
new_genes_ccle <- genes_ccle[common_genes_breast_ccle,]

pred_metab_mat_breast <- matrix(ncol = length((colnames(genes_breast))))
colnames(pred_metab_mat_breast) <- colnames(genes_breast)

# pb = txtProgressBar(min = 0, max = length(common_met_breast_ccle), initial = 0, style = 3)
# for(i in common_met_breast_ccle){
#   metabolite <- metab_ccle[i,]
#   model<-cv.glmnet(t(new_genes_ccle),metabolite,alpha=1)
#   #plot(model)
# 
#   # Get the best model
#   bestcv<-which(model$lambda==model$lambda.min)
#   betas<-model$glmnet.fit$beta[,bestcv]
#   intercept<-model$glmnet.fit$a0[bestcv]
#   sigvars<-names(betas)[betas!=0]
#   bestmodel<-betas[sigvars]
# 
#   # R^2
#   #model$glmnet.fit$dev.ratio[bestcv]
# 
#   # Predict
#   predvalues<-as.numeric(predict(model,newx=t(new_genes_breast),s = model$lambda.min))
#   #plot(values,predvalues)
#   pred_metab_mat_breast <- rbind(pred_metab_mat_breast, predvalues)
#   setTxtProgressBar(pb,match(i,common_met_breast_ccle))
# 
# }
# close(pb)
# pred_metab_mat_breast <- pred_metab_mat_breast[-1,]
# rownames(pred_metab_mat_breast) <- common_met_breast_ccle
# save(pred_metab_mat_breast, file="data/final_ds/lm/pred_metab_mat_breast.rda")
load("data/final_ds/lm/pred_metab_mat_breast.rda")

cc_thr <- 0.3
corr_outputs_breast <- correlation(metab_breast[common_met_breast_ccle,], pred_metab_mat_breast[common_met_breast_ccle,], cc_thr)
cors_breast_spearman <- corr_outputs_breast[[1]]
cors_breast_pearson <- corr_outputs_breast[[2]]

names(cors_breast_spearman) <- common_met_breast_ccle
names(cors_breast_pearson) <- common_met_breast_ccle

significant_breast <- corr_outputs_breast[[3]]

length(common_met_breast_ccle)# Among 70 common metabolites predicted...
length(significant_breast) # just 7/70 predictions are significant

# Saving with different names to better manage the step in which I will compare the correlations performances of 
#  different tools
breast_sp_lm <- cors_breast_spearman
breast_pe_lm <- cors_breast_pearson

save(breast_sp_lm, file = "data/final_ds/corr/breast_sp_lm.rda")
save(breast_pe_lm, file = "data/final_ds/corr/breast_pe_lm.rda")

# Testing how much dist is diff from zero
shapiro.test(cors_breast_spearman) #p-value = 0.4594 --> can use t-test
shapiro.test(cors_breast_pearson) #p-value = 0.5746 --> can use t-test
t1_breast <- t.test(cors_breast_spearman, alternative = "greater") #p-value = 0.5008
t2_breast <- t.test(cors_breast_pearson, alternative = "greater") # p-value = 0.462


# For thesis
png(paste0("results/figures/7_lasso_pred/correlations_breast.png"),w=2900,h=2500,res=600)
plot(density(cors_breast_spearman),col="darkgoldenrod3",lwd=3, xlim=c(-1,1), ylim = c(0,max(density(cors_breast_spearman)$y)), 
     main = "GSE37751 correlations ditribution", xlab = "metabolite CC")
lines(density(cors_breast_pearson), col = "brown3", lty = 3, lwd = 2)
abline(v=0,lty=2, col="gray")
mtext(paste0("Spearman p-value = ",round(t1_breast$p.value,3),"     Pearson p-value = ",round(t2_breast$p.value,3)))
legend("topleft", legend = c("Spearman", "Pearson"), lwd = 2,
       col = c("darkgoldenrod3", "brown3"), lty = c(1,2))
dev.off()

# HASSAN #### 

load("data/final_ds/genes_hassan.rda")
load("data/final_ds/metab_hassan.rda")
common_met_hassan_ccle <- intersect(rownames(metab_ccle), rownames(metab_hassan))
common_genes_hassan_ccle <- intersect(rownames(genes_ccle), rownames(genes_hassan))

new_genes_hassan <- genes_hassan[common_genes_hassan_ccle,]
new_genes_ccle <- genes_ccle[common_genes_hassan_ccle,]

pred_metab_mat_hassan <- matrix(ncol = length((colnames(genes_hassan))))
colnames(pred_metab_mat_hassan) <- colnames(genes_hassan)

# pb = txtProgressBar(min = 0, max = length(common_met_hassan_ccle), initial = 0, style = 3) 
# for(i in common_met_hassan_ccle){
#   metabolite <- metab_ccle[i,]
#   model<-cv.glmnet(t(new_genes_ccle),metabolite,alpha=1)
#   #plot(model)
#   
#   # Get the best model
#   bestcv<-which(model$lambda==model$lambda.min)
#   betas<-model$glmnet.fit$beta[,bestcv]
#   intercept<-model$glmnet.fit$a0[bestcv]
#   sigvars<-names(betas)[betas!=0]
#   bestmodel<-betas[sigvars]
#   
#   # R^2
#   #model$glmnet.fit$dev.ratio[bestcv]
#   
#   # Predict
#   predvalues<-as.numeric(predict(model,newx=t(new_genes_hassan),s = model$lambda.min))
#   #plot(values,predvalues)
#   pred_metab_mat_hassan <- rbind(pred_metab_mat_hassan, predvalues)
#   setTxtProgressBar(pb,match(i,common_met_hassan_ccle))
#   
# }
# close(pb)
# pred_metab_mat_hassan <- pred_metab_mat_hassan[-1,]
# rownames(pred_metab_mat_hassan) <- common_met_hassan_ccle
# save(pred_metab_mat_hassan, file="data/final_ds/lm/pred_metab_mat_hassan.rda")
load("data/final_ds/lm/pred_metab_mat_hassan.rda")

cc_thr <- 0.3
corr_outputs_hassan <- correlation(metab_hassan[common_met_hassan_ccle,], pred_metab_mat_hassan[common_met_hassan_ccle,], cc_thr)
cors_hassan_spearman <- corr_outputs_hassan[[1]]
cors_hassan_pearson <- corr_outputs_hassan[[2]]

names(cors_hassan_spearman) <- common_met_hassan_ccle
names(cors_hassan_pearson) <- common_met_hassan_ccle

significant_hassan <- corr_outputs_hassan[[3]]

length(common_met_hassan_ccle)# Among 5 common metabolites predicted...
length(significant_hassan) # ...just 1/5 predictions are significant

# Saving with different names to better manage the step in which I will compare the correlations performances of 
#  different tools
hassan_sp_lm <- cors_hassan_spearman
hassan_pe_lm <- cors_hassan_pearson

save(hassan_sp_lm, file = "data/final_ds/corr/hassan_sp_lm.rda")
save(hassan_pe_lm, file = "data/final_ds/corr/hassan_pe_lm.rda")

# Testing how much dist is diff from zero
shapiro.test(cors_hassan_spearman) #p-value = 0.9724 --> can use t-test
shapiro.test(cors_hassan_pearson) #p-value = 0.5708 --> can use t-test
t1_hassan <- t.test(cors_hassan_spearman, alternative = "greater") #p-value = 0.4523
t2_hassan <- t.test(cors_hassan_pearson, alternative = "greater") # p-value = 0.5051


# For thesis
png(paste0("results/figures/7_lasso_pred/correlations_hassan.png"),w=2900,h=2500,res=600)
plot(density(cors_hassan_spearman),col="darkgoldenrod3",lwd=3, xlim=c(-1.5,1.5), ylim = c(0,max(density(cors_hassan_spearman)$y)), 
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
common_met_zamp_ccle <- intersect(rownames(metab_ccle), rownames(metab_zamp))
common_genes_zamp_ccle <- intersect(rownames(genes_ccle), rownames(genes_zamp))

new_genes_zamp<-genes_zamp[common_genes_zamp_ccle,]
new_genes_ccle <- genes_ccle[common_genes_zamp_ccle,]

pred_metab_mat_zamp <- matrix(ncol = length((colnames(genes_zamp))))
colnames(pred_metab_mat_zamp) <- colnames(genes_zamp)

# pb = txtProgressBar(min = 0, max = length(common_met_zamp_ccle), initial = 0, style = 3) 
# for(i in common_met_zamp_ccle){
#   metabolite <- metab_ccle[i,]
#   model<-cv.glmnet(t(new_genes_ccle),metabolite,alpha=1)
#   #plot(model)
#   
#   # Get the best model
#   bestcv<-which(model$lambda==model$lambda.min)
#   betas<-model$glmnet.fit$beta[,bestcv]
#   intercept<-model$glmnet.fit$a0[bestcv]
#   sigvars<-names(betas)[betas!=0]
#   bestmodel<-betas[sigvars]
#   
#   # R^2
#   #model$glmnet.fit$dev.ratio[bestcv]
#   
#   # Predict 
#   predvalues<-as.numeric(predict(model,newx=t(new_genes_zamp),s = model$lambda.min))
#   #plot(values,predvalues)
#   pred_metab_mat_zamp <- rbind(pred_metab_mat_zamp, predvalues)
#   setTxtProgressBar(pb,match(i,common_met_zamp_ccle))
#   
# }
# close(pb)
# pred_metab_mat_zamp <- pred_metab_mat_zamp[-1,]
# rownames(pred_metab_mat_zamp) <- common_met_zamp_ccle
# save(pred_metab_mat_zamp, file="data/final_ds/lm/pred_metab_mat_zamp.rda")
load("data/final_ds/lm/pred_metab_mat_zamp.rda")

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
      temp_pred <- c(temp_pred, pred_metab_mat_zamp[j,v])
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

length(common_met_zamp_ccle) # Among 48 common metabolites predicted...
length(significant_zamp) # ...just 8/48 predictions are significant

# Saving with different names to better manage the step in which I will compare the correlations performances of 
#  different tools
zamp_sp_lm <- cors_zamp_spearman
zamp_pe_lm <- cors_zamp_pearson

save(zamp_sp_lm, file = "data/final_ds/corr/zamp_sp_lm.rda")
save(zamp_pe_lm, file = "data/final_ds/corr/zamp_pe_lm.rda")

# Testing how much dist is diff from zero
shapiro.test(cors_zamp_spearman) #p-value = 0.5774 --> can use t-test
shapiro.test(cors_zamp_pearson) #p-value = 0.445 --> can use t-test
t1_zamp <- t.test(cors_zamp_spearman, alternative = "greater") #p-value = 0.008662
t2_zamp <- t.test(cors_zamp_pearson, alternative = "greater") # p-value = 0.0118

# For thesis
png(paste0("results/figures/7_lasso_pred/correlations_zamp.png"),w=2900,h=2500,res=600)
plot(density(cors_zamp_spearman),col="darkgoldenrod3",lwd=3, xlim=c(-1,1), ylim = c(0,max(density(cors_zamp_spearman)$y)), 
     main = "GSE32474 correlations ditribution", xlab = "metabolite CC")
lines(density(cors_zamp_pearson), col = "brown3", lty = 3, lwd = 2)
abline(v=0,lty=2, col="gray")
mtext(paste0("Spearman p-value = ",round(t1_zamp$p.value,3),"     Pearson p-value = ",round(t2_zamp$p.value,3)))
legend("topleft", legend = c("Spearman", "Pearson"), lwd = 2,
       col = c("darkgoldenrod3", "brown3"), lty = c(1,2))
dev.off()
s