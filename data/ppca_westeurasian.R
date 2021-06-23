library(softImpute)
library(tidyverse)
library(admixr)
source("ppca.R")

set.seed(1)

eig <- "subdata/westeurasian1"
e <- eigenstrat(eig)
data <- admixr::read_geno(e)  
labels <- names(data)
data <- as.matrix(data)
data[data==9] <- NA
colnames(data) <- labels
pmiss <- colMeans(!is.na(data))
d2 <- data[rowMeans(!is.na(data)) > .75,]

if(F){
fit <- psvd(d2, max_miss_ind=.5, thresh=1e-4, 
            trace.it=T, rank=10)
saveRDS(fit, "fit.rds")
}

ind = admixr::read_ind(e)
PCS = fit$PCS
PCS$pop = ind$label
pcs = PCS %>% group_by(pop) %>% summarize_at(vars(PC1:PC10), .funs=mean)

#get circle based on basque, turkish
circle = pcs %>% filter(pop %in% c("Basque", "Turkish")) %>% select(PC1, PC2)
radius = (dist(circle) %>% as.matrix)[1,2] / 2
origin = colMeans(circle)
circle_df=tibble(x0=origin[1], y0=origin[2], radius=radius)




rd2 = rowMeans(d2)
d3 = d2[!is.na(rd2),] / 2
mu = rowMeans(d3)
Y = d3 - mu
SVD = svd(Y)
lambda = SVD$d
#lambda= lambda/sum(lambda)
P = SVD$v
S = t(t(P) * lambda)
tS = t(S)


h2 = mean(mu * (1-mu)) + mean(Y[,1] * (1-Y[,1])) - 2 * mean(Y[,1] * mu)
h1 = mean(d3[,1] * (1-d3[,1]))

T1 = sum(colSums((1 - 2 * mu) * SVD$u) * pki ) / n
M = sum(mu * (1-mu)) / n
SQ = sum(S[1,]^2)/ n
h3 = T1 + M - SQ

n = nrow(Y)
npcs = ncol(Y)

T1t = colSums((1-2*mu)*SVD$u) * pki
SQt = S[1,]^2
Mt = M *n / npcs


f3test = function(idx, id1, id2){
    f3 =  (tS[,idx] - tS[,id1]) * (tS[,idx]  - tS[,id2])
    pki = (SVD$v[idx,] * SVD$d)
    T1t = colSums((1-2*mu)*SVD$u) * pki
    Mt = M *n / npcs
    SQt = S[idx,]^2

}
