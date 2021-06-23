setwd("/home/student/Desktop/vittoria")

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

train_pred <- function(ntree, common_met, met1, genes1, genes2){
  pred_mat <- c()
  #pb = txtProgressBar(min = 0, max = length(common_met), initial = 0, style = 3)
  for(met in common_met){
    Tmet1 <- t(met1)
    Tgenes1 <- t(genes1)
    metabolite <- Tmet1[,met]
    
    # create a temporary variable which stores the response and the significant predictors only
    current <- cbind(metabolite, Tgenes1)
    
    # Train the random forest
    rf <- randomForest(metabolite ~ ., data = current, ntree = ntree)
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
    #setTxtProgressBar(pb,match(met,common_met))
  }
  #close(pb)
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

cc_thr <- 0.3

ntrees <- c(10, 20, 50, 100) 
plist <- list() 
for(ntree in ntrees){
  print(ntree)
  ps <- c()
  pb = txtProgressBar(min = 0, char = "*" , max = 100, initial = 0, style = 3) 
  for(i in 1:100){ 
    set.seed(i)
    pred_rf_NCI60 <- train_pred(ntree, common_met_NCI60_ccle, metab_ccle, rf_genes_ccle, rf_genes_NCI60)
    rownames(pred_rf_NCI60) <- rnames
    colnames(pred_rf_NCI60) <- cnames
    save(pred_rf_NCI60, file=paste0("data/ntrees_selection/pred_rf_NCI60_",ntree,"_",i,".rda"))
    
    corr_outputs_NCI60 <- correlation(metab_NCI60[common_met_NCI60_ccle,], pred_rf_NCI60[common_met_NCI60_ccle,], cc_thr)
    cors_NCI60_spearman <- corr_outputs_NCI60[[1]]
    names(cors_NCI60_spearman) <- common_met_NCI60_ccle
    t1 <- wilcox.test(cors_NCI60_spearman, alternative = "greater")
    ps <- c(ps, t1$p.value)
    setTxtProgressBar(pb,i)
  }
  close(pb)
  plist[[as.character(ntree)]]<-ps
}

png(paste0("results/figures/8_rf_pred/ntree_selection.png"),w=2900,h=2500,res=600)
boxplot(plist,xlab="ntrees",ylab="pvalue",main="Random Forest Performance")
mtext("Metabolite Prediction, Train: CCLE, Test: NCI-60",cex=0.9,font=2)
dev.off()























