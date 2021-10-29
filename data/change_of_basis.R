source("fscripts.R")
library(glue)
BTHEME = theme_classic() + theme(legend.position="none")

f2s = admixtools::read_f2("westeurasian1")
pcmat = pca_from_f2s(f2s)
pcs = pca_from_pcmat(pcmat)
idw = "AA"
idx = "Sardinian"
idy = "Basque"
idz = "Han"

ids = c(idx, idy)
n = nrow(pcs)

f2mat = f2(f2s, unique_only=F) %>% select(-se) %>% 
    pivot_wider(names_from=pop2, values_from=est) %>%
    column_to_rownames('pop1') %>% as.matrix
cmat = double_center(-f2mat)

a = treelet::Run_JTree(cmat, n-1, 1:(n-1))
B = lapply(a$basis, function(b){
               rownames(b) = rownames(cmat)
               colnames(b) = rownames(cmat)
               b %>% as.data.frame %>%
               rownames_to_column('pop1')
    })


M = B %>% 
    bind_rows(.id='lvl') %>% 
    mutate(lvl=as.integer(lvl)) %>% 
    pivot_longer(-1:-2, names_to='pop2', values_to='bval')

B = a$basis[[n-1]]
x = function(b)ifelse( abs(b) < 1e-5, NA, (b - min(b)) / diff(range(b)))
S = apply(B, 2,x)


