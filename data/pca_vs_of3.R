require(tidyverse)
require(admixtools)

OUTGROUP = 'Mbuti'

df2mat <- function(f2_df){
    f2_df%>% select(-se) %>% 
        arrange(pop1, pop2) %>%
        pivot_wider(names_from=pop2, values_from=est) %>%
        column_to_rownames('pop1') %>% as.matrix
}

double_center <- function(mat){
    m = rowMeans(mat)
    x = t(mat - m)-m; 
    x = x-mean(x)
}

pca_from_mat <- function(mat){
    df = tibble()
    e = eigen(mat)
    e$vectors = t(t(e$vectors) )#* e$values)
    colnames(e$vectors) = sprintf("PC%d", 1:nrow(mat))
    as_tibble(e$vectors) %>% mutate(pop=colnames(mat)) %>%
        mutate_at(vars(starts_with("PC")), function(X)X *sign(X[1]))  %>%
        select(pop, starts_with("PC"))

}

project_on_pcs <- function(pcs, evals, new_dists, old_dist_to_center){
    1/evals * t(pcs[,-1)a
}

f2s = admixtools::f2_from_precomp("pca_test1/")
f2_raw = f2(f2s, unique=F)
pops = rownames(f2s[,,1])
inpops = pops[pops != OUTGROUP]


f2mat = f2_raw %>% filter(pop1 != OUTGROUP, pop2 != OUTGROUP) %>% df2mat
f3s = f3(f2s, pop1=OUTGROUP, pop2=inpops, pop3=inpops, unique=F)
f3mat = f3s %>% mutate(pop1=pop3) %>% select(-z, -p, -pop3) %>% df2mat


cf2mat = -double_center(f2mat) /2
#diag(cf2mat) = diag(cf2mat) + 1e-3
#cf2mat = double_center(cf2mat)
pcs = cf2mat %>% pca_from_mat
e = eigen(cf2mat)$values

OUTGROUP = 'Bulgarian'
inpops = pops[!pops %in% c("Mbuti", OUTGROUP)]
xf3s = f3(f2s, pop1=OUTGROUP, pop2=inpops, pop3=inpops, unique=F)
xf3mat = xf3s %>% mutate(pop1=pop3) %>% select(-z, -p, -pop3) %>% df2mat


f3_spectrum <- function(x,y,z, ...)f4_spectrum(x,y,x,z, ...)
f4_spectrum <- function(x,y,z, w, pca_df, evals){
    v1 = pca_df[pca_df$pop == x, -1]
    v2 = pca_df[pca_df$pop == y, -1]
    v3 = pca_df[pca_df$pop == z, -1]
    v4 = pca_df[pca_df$pop == w, -1]
    return(evals * (v1 - v2) * (v3 - v4) %>% unname %>% unlist)
}


