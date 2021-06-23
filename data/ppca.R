library(softImpute)
library(admixr)

psvd.eig <- function(eig=test_file, ...){
    e <- eigenstrat(eig)
    data <- admixr::read_geno(e)  
    labels <- names(data)
    data <- as.matrix(data)
    data[data==9] <- NA
    colnames(data) <- labels
    return(psvd(data, ...))
}

psvd <- function(data, max_miss_snp=0.5, max_miss_ind=0.0, lambda=0, type='als', rank=4, 
                 only_pcs=TRUE, ...){

    if(is.null(colnames(data)))colnames(data) <- 1:ncol(data)

    #X is random allele
    X <- data

    Z <- X[rowMeans(!is.na(X)) > max_miss_snp, ]


    ind_completeness <- colMeans(!is.na(Z)) 

    if(min(ind_completeness) < max_miss_ind){
        print("removing individuals:")
        w <- which(ind_completeness < max_miss_ind)
        for(i in w){
            print(sprintf("%s: %s", colnames(data)[i], ind_completeness[i]))
        }
    }
    Z <- Z[,ind_completeness >= max_miss_ind] 

    print(sprintf("PCA-matrix size:  %s", dim(Z)))
    G <- t(Z - rowMeans(Z, na.rm=T)) #- colMeans(Z, na.rm=T)  + mean(Z, na.rm=T)

    fit <- softImpute(G, lambda=lambda, rank=rank, type=type, ...)

    PCS = t(t(fit$u) * fit$d)
    PCS <- data.frame(name=rownames(G), PCS)
    names(PCS)[-1] <- sprintf("PC%s", 1:ncol(fit$u))

    if(only_pcs)return(PCS)

    return(list(PCS=PCS, d=fit$d, loadings=fit$v, data=data))
}

test <- function(){
    library(tidyverse)
    test_file="~benjamin_peter/projects/100a/plink/100a2_all/archaics"
    pca <- psvd.eig(test_file, lambda=0, max_miss_snp=0.48, 
                    max_miss_ind=.2,
                        rank=6)  
    ggplot(pca, aes(x=PC1, y=PC2, label=name)) + geom_text()     

}
