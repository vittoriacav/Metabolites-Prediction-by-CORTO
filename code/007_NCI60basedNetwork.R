# Here we try to train the network on the NCI60 dataset and predict on CLLE to validate it
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

# Loading genes and metabolites of NCI60
load("data/final_ds/genes_NCI60.rda")
names_gene_NCI60 <- rownames(genes_NCI60)
load("data/final_ds/metab_NCI60.rda")
names_metab_NCI60 <- rownames(metab_NCI60)


# Loading genes and metabolites of CCLE
load("data/final_ds/metab_ccle.rda")
real_metab_names_ccle <- rownames(metab_ccle)
load("data/final_ds/genes_ccle.rda")
names_gene_ccle <- rownames(genes_ccle)

# Expression + metabolite data NCI60
NCI60_netMat <- rbind(genes_NCI60, metab_NCI60)
# Train the network
#network_NCI60 <- corto(NCI60_netMat,centroids=names_metab_NCI60, nbootstraps=100, p=1e-3, nthreads=7, verbose=FALSE)
#save(network_NCI60,file="code/network_NCI60.rda")
load("code/network_NCI60.rda")

# Prediction on CCLE using NCI60-based network
pred_ccle <- mra(genes_ccle, regulon = network_NCI60)
pred_metab_names_ccle <- rownames(pred_ccle)
common_ccle <- intersect(pred_metab_names_ccle, real_metab_names_ccle) 
length(common_ccle) # 34 in common, however....

cc_thr <- 0.3
corr_outputs <- correlation(metab_ccle[common_ccle,], pred_ccle[common_ccle,], cc_thr)
cors_spearman <- corr_outputs[[1]]
cors_pearson <- corr_outputs[[2]]

names(cors_spearman) <- common_ccle
names(cors_pearson) <- common_ccle

significant <- corr_outputs[[3]]

length(significant) # ...however just just 8/34 are >0.3

# Testing how much dist is diff from zero
shapiro.test(cors_spearman) #p-value = 0.3205 --> can use t-test
shapiro.test(cors_pearson) #p-value = 0.2168 --> can use t-test
t1<- t.test(cors_spearman, alternative = "greater") #p-value = 1.014228e-05
t2<- t.test(cors_pearson, alternative = "greater") # p-value = 1.944393e-05



png(paste0("results/figures/6_validation/correlations_CCLE.png"),w=2900,h=2500,res=600)
plot(density(cors_spearman),col="darkgoldenrod3",lwd=3, xlim=c(-1,1), ylim = c(0,2.5), 
     main = "CCLE correlations ditribution", xlab = "metabolite CC")
lines(density(cors_pearson), col = "brown3", lty = 3, lwd = 2)
abline(v=0,lty=2, col="gray")
mtext(paste0("Spearman p-value = ",round(t1$p.value,3),"     Pearson p-value = ",round(t2$p.value,3)))
legend("topleft", legend = c("Spearman", "Pearson"), lwd = 2,
       col = c("darkgoldenrod3", "brown3"), lty = c(1,2))
dev.off()

# CYTOSCAPE ####
# Get names of the metabolites nodes present in the network 
nodes <- names(network_NCI60) # metabolite nodes
nodes <- gsub(" ", "_", nodes) #substituting the space with underscore to remove problems in Cytoscape
names(network_NCI60) <- nodes

# Creating txt matrix for the whole network_NCI60
# each row contains: metabolite name; gene name; tfmode; likelihood

net_matrix_NCI60<-matrix(nrow=0,ncol=4)

pb<-txtProgressBar(0,length(nodes),style=3)
for(i in 1:length(nodes)){
  met_genes <- names(network_NCI60[[nodes[i]]][["tfmode"]])
  met_tfmode <- network_NCI60[[nodes[i]]][["tfmode"]]
  met_lik <- network_NCI60[[nodes[i]]][["likelihood"]]
  met<-rep(nodes[i],length(met_genes))
  new_row<-cbind(met,met_genes, met_tfmode, met_lik)
  net_matrix_NCI60<-rbind(net_matrix_NCI60,new_row)
  setTxtProgressBar(pb,i)
}

# Removing quotes due to problems in Cytoscape
net_matrix_NCI60 <- noquote(net_matrix_NCI60)
write.table(net_matrix_NCI60, file="tables_cyt/net_matrix_NCI60.txt", row.names=FALSE, col.names=TRUE, quote = FALSE, sep = "\t")

# Creating mapping table used to differentiate metabolites from genes on Cytoscape
type <- rep("node", length(nodes))
nodes_matrix <- cbind(nodes,type)
write.table(nodes_matrix, file="tables_cyt/nodes_matrix_NCI60.txt", row.names=FALSE, col.names=TRUE, quote = FALSE, sep = "\t")



