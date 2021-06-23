#source("../shared/functions/qol.R")
setwd("/home/student/Desktop/vittoria")

# EnrichR
library(enrichR)
library(ggplot2)

loadfile <- function(f) {
    env <- new.env()
    nm <- load(f, env)[1]
    env[[nm]]
}



pvalues_col <- function(ps, lenDn){
    #ps = pvalues
    col_vector <- c()
    for(i in 1:length(ps)){
        if(i <= lenDn){ #neg ones
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
            if(ps[i] <= 0.001){
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

dbs <- listEnrichrDbs()
dbs<-c("Reactome_2016",
       "WikiPathways_2019_Human",
       "GO_Biological_Process_2018",
       "GO_Cellular_Component_2018",
       "GO_Molecular_Function_2018")

metabolites_NCI60 <- setNames(c("data/final_ds/GO/alanine_go.rda",
                      "data/final_ds/GO/allantoin_go.rda",
                      "data/final_ds/GO/arginine_go.rda",
                      "data/final_ds/GO/creatinine_go.rda",
                      "data/final_ds/GO/cytidine_go.rda",
                      "data/final_ds/GO/guanosine_go.rda",
                      "data/final_ds/GO/histidine_go.rda",
                      "data/final_ds/GO/inosine_go.rda",
                      "data/final_ds/GO/isoleucine_go.rda",
                      "data/final_ds/GO/lactate_go.rda",
                      "data/final_ds/GO/leucine_go.rda",
                      "data/final_ds/GO/methionine_go.rda",
                      "data/final_ds/GO/phenylalanine_go.rda",
                      "data/final_ds/GO/serine_go.rda",
                      "data/final_ds/GO/sorbitol_go.rda",
                      "data/final_ds/GO/tyrosine_go.rda",
                      "data/final_ds/GO/uracil_go.rda",
                      "data/final_ds/GO/xanthine_go.rda",
                      "data/final_ds/GO/xanthosine_go.rda"),
                    c("ALANINE", "ALLANTOIN", "ARGININE", "CREATININE",
                      "CYTIDINE", "GUANOSINE", "HISTIDINE", "INOSINE",
                      "ISOLEUCINE", "LACTATE", "LEUCINE", "METHIONINE",
                      "PHENYLALANINE", "SERINE", "SORBITOL", "TYROSINE",
                      "URACIL", "XANTHINE", "XANTHOSINE"))


for(i in 1:length(metabolites_NCI60)){
    print(names(metabolites_NCI60)[i])
    file<-metabolites_NCI60[i]
    metab<-names(metabolites_NCI60[i])
    current_file<-loadfile(file)
    
    gup <- c()
    gdn <- c()
    
    for(j in 1:length(rownames(current_file))){
        if(current_file[j,2] == "+"){
            gup <- c(gup, current_file[j,1])
        }else{
            gdn <- c(gdn, current_file[j,1])
        }
    }
    
    if(length(gup) >= 10){
        eup <- enrichr(gup,dbs)
    }
    if(length(gdn) >= 10){
        edn <- enrichr(gdn,dbs)
    }
    if(length(gup) < 10 && length(gdn) < 10){
        #print("here")
        next
    }
    for(db in dbs){
        if(length(gup)>=10){ # GO has been done
            up<-eup[[db]]
            up$type<-"up"
            #up<-up[order(up$Combined.Score),]
        }else{
            up <- matrix(nrow=0, ncol=0)
        }
        if(length(gdn)>=10){ # GO has been done
            dn<-edn[[db]]
            dn$type<-"down"
            #dn<-dn[order(-dn$Combined.Score),]
            #dn$Combined.Score<-(-1)*dn$Combined.Score
        }else{
            dn <- matrix(nrow=0, ncol=0)
        }
        
        if(dim(up)[1]+dim(dn)[1] >= 30){
            #print("yes")
            if(dim(up)[1] >= 15 && dim(dn)[1] >= 15){
                #print("1")
                up <- up[1:15,]
                up<-up[order(up$Combined.Score),]
                dn <- dn[1:15,]
                dn<-dn[order(-dn$Combined.Score),]
                dn$Combined.Score<-(-1)*dn$Combined.Score
            }else if(dim(up)[1] < 15){
                #print("2")
                dn <- dn[1:(30-dim(up)[1]),]
                if(dim(up)[1]!=0){
                    up<-up[order(up$Combined.Score),]
                }
                dn<-dn[order(-dn$Combined.Score),]
                dn$Combined.Score<-(-1)*dn$Combined.Score
            }else if(dim(dn)[1] < 15){
                #print("3")
                up <- up[1:(30-dim(dn)[1]),]
                up<-up[order(up$Combined.Score),]
                if(dim(dn)[1]!=0){
                    dn<-dn[order(-dn$Combined.Score),]
                    dn$Combined.Score<-(-1)*dn$Combined.Score
                }
            }
        }else{
            #print("no")
            if(dim(up)[1] ==0){
                dn<-dn[order(-dn$Combined.Score),]
                dn$Combined.Score<-(-1)*dn$Combined.Score  
            }else if (dim(dn)[1] == 0){
                up<-up[order(up$Combined.Score),]
            }else{
                up<-up[order(up$Combined.Score),]
                dn<-dn[order(-dn$Combined.Score),]
                dn$Combined.Score<-(-1)*dn$Combined.Score
            }
    
            
        }
        
        gos<-rbind(dn,up)
        gos$Term<-gsub(" WP.+","",gos$Term)
        gos$Term<-gsub("_Homo sapiens_.+","",gos$Term)
        
        
        # Bar plot
        toplot<-setNames(gos$Combined.Score,gos$Term)
        pvalues <- gos$P.value
        png(paste0("results/GO/plots/met_NCI60/",metab,"_",db,".png"),w=8000,h=5000,res=400,family="xkcd")
        par(mar=c(4,1,3,1))
        bp<-barplot(toplot,horiz=TRUE,xlab="Combined EnrichR score",
                    xlim= 1.3*c(-max(abs(toplot)),max(abs(toplot))),
                    main=paste0("Top ",db," pathways in ",metab),
                    col=pvalues_col(pvalues, dim(dn)[1]),
                    yaxt="n",cex.main=2,
                    xpd = F
        )
        #legend(x =  max(abs(toplot))-100, y = 30, xjust =0, legend = c(rep(c("<0.001", "<0.01", "<0.05", "<0.1", "<1"), 2)), ncol = 2, 
        #       fill = c("cadetblue4","cadetblue3", "cadetblue2", "cadetblue1", "azure1", "brown4","brown3","brown2","brown1","coral")
        #)
        if(dim(up)[1] > 0 && dim(dn)[1] > 0){
            text(0,bp[1:dim(dn)[1],1],names(toplot)[1:dim(dn)[1]],pos=4)
            text(0,bp[(dim(dn)[1]+1):dim(gos)[1],1],names(toplot)[(dim(dn)[1]+1):dim(gos)[1]],pos=2)
        }else if(dim(up)[1] == 0){
            text(0,bp[1:dim(gos)[1]],names(toplot)[1:dim(gos)[1]],pos=4)
        }else if(dim(dn)[1] == 0){
            text(0,bp[1:dim(gos)[1]],names(toplot)[1:dim(gos)[1]],pos=2)
        }
        
        dev.off() 
    }
}
