require(admixtools)
require(tidyverse)
require(Matrix)
require(ggforce)


pca_from_f2s <- function(f2s, force_nonneg=F){
    f2 = f2(f2s, unique_only=F, sure=T) 
    n = dim(f2s)[1]
    f2mat = f2 %>%
        mutate(est=if(force_nonneg){pmax(est, 0)}else{est}) %>%
        select(-se) %>%
        pivot_wider(names_from=pop2, values_from=est, values_fill=0, names_sort=T) %>%
        arrange(pop1) %>%
        column_to_rownames('pop1') %>% as.matrix

    cmat = double_center(f2mat)

    #E = eigen(nearPD(-cmat/2)[[1]])
    E = eigen(-cmat/2)
    E$values[E$values < 0] = 0

    pcs = t(t(E$vectors) * sqrt(pmax(E$values, 0))) 
    rownames(pcs) = rownames(f2mat)
    return(pcs)
}

double_center = function(m){
    cmat = t(t(m - rowMeans(m) ) - colMeans(m)) + mean(m)
}

pca_from_f4 <- function(f4){
    f4mat = f4 %>% 
        select(pop1, pop2, est) %>% 
        pivot_wider(names_from=pop2, values_from=est, values_fill=0, names_sort=T) %>% 
        column_to_rownames('pop1') %>% as.matrix %>% abs

    cmat = double_center(f4mat)

    #E = eigen(nearPD(-cmat/2)[[1]])
    E = eigen(-cmat/2)
    E$values[E$values < 0] = 0

    pcs = t(t(E$vectors) * sqrt(pmax(E$values, 0))) 
    rownames(pcs) = rownames(f4mat)
    return(pcs)
}

pca_from_pcmat <- function(pcs){
    pcs = pcs %>% as.data.frame %>% rownames_to_column("pop")
    names(pcs)[-1] = sprintf("PC%d", 1:(ncol(pcs) - 1))
    return(pcs)
}


f3_from_pc_matrix <-function(pcmat, idx="Basque", idy="Turkish"){
    pcmat = t(pcmat)
    f3mat = t((pcmat - pcmat[,idx]) * (pcmat - pcmat[,idy]))
    colnames(f3mat) = sprintf("PC%d", 1:ncol(pcmat))
    return(f3mat)
}

of3_from_pc_matrix <-function(pcmat, idx="Basque", idy="Turkish"){
    pcmat = t(pcmat)
    f3mat = t((pcmat[,idx] - pcmat[,idy]) * (pcmat[,idx] - pcmat))
    colnames(f3mat) = sprintf("PC%d", 1:ncol(pcmat))
    return(f3mat)
}

f4_from_pc_matrix_single <-function(pcmat, A, B, C, D){
    pcmat = t(pcmat)
    f4mat = t((pcmat[,A] - pcmat[,B]) * (pcmat[,C] - pcmat[,D]))
    return(f4mat)
}

f4_from_pc_matrix <-function(pcmat, idx="Basque", idz="Turkish", idy="Scottish"){
    pcmat = t(pcmat)
    f4mat = t((pcmat - pcmat[,idx]) * (pcmat[,idy] - pcmat[,idz]))
    colnames(f4mat) = sprintf("PC%d", 1:ncol(f4mat))
    return(f4mat)
}
f2_from_pc_matrix <-function(pcmat, idx="Basque"){
    pcmat = t(pcmat)
    f2mat = t((pcmat - pcmat[,idx]) * (pcmat - pcmat[,idx]))
    colnames(f2mat) = sprintf("PC%d", 1:ncol(f2mat))
    return(f2mat)
}

f4_from_pc_matrix2 <-function(pcmat, px="Basque", py="Turkish"){
    n = nrow(pcmat)
    x = lapply(1:n, function(id) f4_from_pc_matrix(pcmat, idx=id, idy=px, idz=py) %>% 
               as.data.frame %>% rownames_to_column("pop2"))
    names(x) = rownames(pcmat)
    f4mat = bind_rows(x, .id="pop1") %>% arrange(pop1, pop2)
    names(f4mat)[-1:-2] = sprintf("PC%d", 1:n)
    return(f4mat)
}
f2_from_pc_matrix2 <-function(pcmat){
    n = nrow(pcmat)
    x = lapply(1:n, function(id) f2_from_pc_matrix(pcmat, idx=id) %>% 
               as.data.frame %>% rownames_to_column("pop2"))
    names(x) = rownames(pcmat)
    f2mat = bind_rows(x, .id="pop1") %>% arrange(pop1, pop2)
    names(f2mat)[-1:-2] = sprintf("PC%d", 1:n)
    return(f2mat)
}

ggsvg <- function(P, fname, ...){
    svg(fname, ...)
    print(P)
    dev.off()
}

