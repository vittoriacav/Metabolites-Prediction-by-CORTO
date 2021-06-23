setwd("/home/student/Desktop/vittoria")
library(corto)
library(beeswarm)
library(gplots)
library(vioplot)


load("data/final_ds/corr/breast_pe_corto.rda")
load("data/final_ds/corr/breast_pe_lm.rda")
load("data/final_ds/corr/breast_pe_rf.rda")

load("data/final_ds/corr/breast_sp_corto.rda")
load("data/final_ds/corr/breast_sp_lm.rda")
load("data/final_ds/corr/breast_sp_rf.rda")


load("data/final_ds/corr/hassan_pe_corto.rda")
load("data/final_ds/corr/hassan_pe_lm.rda")
load("data/final_ds/corr/hassan_pe_rf.rda")

load("data/final_ds/corr/hassan_sp_corto.rda")
load("data/final_ds/corr/hassan_sp_lm.rda")
load("data/final_ds/corr/hassan_sp_rf.rda")


load("data/final_ds/corr/NCI60_pe_corto.rda")
load("data/final_ds/corr/NCI60_pe_lm.rda")
load("data/final_ds/corr/NCI60_pe_rf.rda")

load("data/final_ds/corr/NCI60_sp_corto.rda")
load("data/final_ds/corr/NCI60_sp_lm.rda")
load("data/final_ds/corr/NCI60_sp_rf.rda")


load("data/final_ds/corr/zamp_pe_corto.rda")
load("data/final_ds/corr/zamp_pe_lm.rda")
load("data/final_ds/corr/zamp_pe_rf.rda")

load("data/final_ds/corr/zamp_sp_corto.rda")
load("data/final_ds/corr/zamp_sp_lm.rda")
load("data/final_ds/corr/zamp_sp_rf.rda")

dbs <- c("NCI60", "breast", "hassan", "zamp")
names(dbs) <- c("NCI60", "GSE37751", "GSE148892", "GSE32474")


for(db in dbs){
  png(paste0("results/figures/9_tools_comparison/",db,".png"), w=6000, h=2700, res =600)
  par(mfrow = c(1,2))
  for(cty in c("sp", "pe")){
    corr_list <- list(CORTO = get(paste0(db,"_",cty,"_corto")), 
                      LASSO = get(paste0(db,"_",cty,"_lm")), 
                      RF = get(paste0(db,"_",cty,"_rf")))
    
    vioplot(corr_list, col = c("#F6CC92", "#FAEBAA", "#FEC5C2"), ylab = "measured vs predicted metabolite CC")  
    abline(h=0, lty = 2)
    
    vioplot(corr_list, col = c("#F6CC92", "#FAEBAA", "#FEC5C2"), ylab = "measured vs predicted metabolite CC", add = T)
    boxplot(corr_list, add = T, col = "#00000000")
    #border = c("#F7AF4D", "#EECD3D", "#FB8D88")
    beeswarm(corr_list, pch = 19, col = c("#DC7F00", "#D1AA00", "#DC0A00"),
             add = T)
    if(cty == "sp"){
      title(main = "Spearman Correlation", 
            sub = paste0("p-value (CORTO - LASSO) = ", round(wilcox.test(corr_list[[1]], corr_list[[2]])$p.value,3),
                         "\np-value (CORTO - RF) = ", round(wilcox.test(corr_list[[1]], corr_list[[3]])$p.value,3)))
      title(main = paste0("Predicting tools comparison in ", names(dbs[match(db, dbs)])),
            outer = T, line = -1, cex.main = 1.5)
    }else{
      title(main = "Pearson Correlation",
            sub = paste0("p-value (CORTO - LASSO) = ", round(wilcox.test(corr_list[[1]], corr_list[[2]])$p.value,3),
                         "\np-value (CORTO - RF) = ", round(wilcox.test(corr_list[[1]], corr_list[[3]])$p.value,3)))
    }
  }
  dev.off()
} 


