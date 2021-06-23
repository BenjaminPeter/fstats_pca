require(admixtools)
require(tidyverse)
require(Matrix)


pca_from_f2s <- function(f2s){
    f2 = f2(f2s)
    n = dim(f2s)[1]
    f2mat = f2 %>%
        mutate(tmp=pop2, pop2=pop1, pop1=tmp) %>% 
        select(-tmp) %>% bind_rows(f2) %>% 
        select(-se) %>% 
        pivot_wider(names_from=pop2, values_from=est, values_fill=0, names_sort=T) %>%
        arrange(pop1) %>%
        column_to_rownames('pop1') %>% as.matrix

    cmat = t(t(f2mat - rowMeans(f2mat) ) - colMeans(f2mat)) + mean(f2mat)

    E = eigen(nearPD(-cmat/2)[[1]])

    pcs = t(t(E$vectors) * sqrt(pmax(E$values, 0))) 
    rownames(pcs) = rownames(f2mat)
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
    return(f4mat)
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


#f2s = admixtools::read_f2("worldfoci2/")
f2s = admixtools::read_f2("f4ratio/")
#to_remove = c("Denisova.DG")
#pops0 = colnames(f2s)
#sample_list = pops0[!pops0 %in% to_remove]
#f2s = admixtools::read_f2("f4ratio/", pops = sample_list)
pcmat = pca_from_f2s(f2s)
pcs = pca_from_pcmat(pcmat)
idx = "Primate_Chimp"
idy = "Altai_Neanderthal.DG"

sample_list = pcs$pop 
#f3 = f3(f2s, pop3=idx, pop2=idy, pop1=sample_list) #admixture
f3 = f3(f2s, pop1=idx, pop2=idy, pop3=sample_list) #outgroup
f4 = f4(f2s, pop3=idx, pop4 = idy, pop1=sample_list, pop2=sample_list) %>% arrange(pop1, pop2)

circle = pcs %>% filter(pop %in% c(idx, idy)) %>% select(-pop)
circle_df = colMeans(circle) %>% t %>% as.data.frame
circle_df$radius = (dist(circle) %>% as.matrix)[1,2] / 2

