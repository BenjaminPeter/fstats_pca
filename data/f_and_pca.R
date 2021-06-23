require(admixtools)
require(tidyverse)
require(Matrix)

f3 = f3(f2s, pop3="Basque", pop2="Turkish", pop1=sample_list)
f4 = f4(f2s, pop2="Basque", pop3="Scottish", pop4="Turkish", pop1=sample_list, pop2) %>% arrange(pop1)
f4.2 = f4(f2s, pop4="Basque", pop3="Turkish", pop1=sample_list, pop2=sample_list) %>% arrange(pop1, pop2)

f2mat = f2 %>%
    mutate(tmp=pop2, pop2=pop1, pop1=tmp) %>% 
    select(-tmp) %>% bind_rows(f2) %>% 
    select(-se) %>% 
    pivot_wider(names_from=pop2, values_from=est, values_fill=0, names_sort=T) %>%
    arrange(pop1) %>%
    column_to_rownames('pop1') %>% as.matrix

f2mat %>% cmdscale(k=61) %>% as.data.frame %>% rename(PC1=V1, PC2=V2) %>% rownames_to_column("pop")

cmat = t(t(f2mat - rowMeans(f2mat) ) - colMeans(f2mat)) + mean(f2mat)

E = eigen(nearPD(-cmat/2)[[1]])

pcs = t(t(E$vectors) * sqrt(pmax(E$values, 0)))
rownames(pcs) = rownames(f2mat)
pcmat = pcs
pcs = pcs %>% as.data.frame %>% rownames_to_column("pop")
names(pcs)[-1] = sprintf("PC%d", 1:(ncol(pcs) - 1))

circle = pcs %>% filter(pop %in% c("Basque", "Turkish")) %>% select(PC1, PC2)
radius = (dist(circle) %>% as.matrix)[1,2] / 2
origin = colMeans(circle)
circle_df=tibble(x0=origin[1], y0=origin[2], radius=radius, pop="0")

circle34 = pcs %>% filter(pop %in% c("Basque", "Turkish")) %>% select(PC3, PC4)
radius = (dist(circle34) %>% as.matrix)[1,2] / 2
origin = colMeans(circle34)
circle_df34=tibble(x0=origin[1], y0=origin[2], radius=radius, pop="0")

f3_from_pc_matrix <-function(pcmat, idx="Basque", idy="Turkish"){
    pcmat = t(pcmat)
    f3mat = t((pcmat - pcmat[,idx]) * (pcmat - pcmat[,idy]))
    return(f3mat)
}

f4_from_pc_matrix <-function(pcmat, idx="Basque", idz="Turkish", idy="Scottish"){
    pcmat = t(pcmat)
    f4mat = t((pcmat - pcmat[,idx]) * (pcmat[,idy] - pcmat[,idz]))
    return(f4mat)
}

f4_from_pc_matrix2 <-function(pcmat, idx="Basque", idy="Turkish"){
    n = nrow(pcmat)
    x = lapply(1:n, function(id) f4_from_pc_matrix(pcmat, idx=id, idy="Basque", idz="Turkish") %>% 
               as.data.frame %>% rownames_to_column("pop2"))
    names(x) = rownames(pcmat)
    f4mat = bind_rows(x, .id="pop1") %>% arrange(pop1, pop2)
    names(f4mat)[-1:-2] = sprintf("PC%d", 1:n)
    return(f4mat)
}

