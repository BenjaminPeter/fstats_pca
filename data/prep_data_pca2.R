require(admixtools)
require(tidyverse)
require(glue)

anno = read_tsv("/mnt/ancient/ModernHuman/ReichLab/reich_public_geno_v44.3/v44.3_HO_public.anno")
anno = anno %>% select(iid=2, pub=4, age=6, lat=11, long=12, country=10, group=8,
cov = 14, n_snps=15, sex=16, lib=17, qc=18)


MISSING = c("Irish_Ulster", "Polish", "German", "Irish", "Shetlandic", "Sorb")

list1 = read_table("WestEurasian_PCA.list", col_names=F)[,1] %>% unlist %>% 
    unname
sample_list = list1[!list1 %in% MISSING]
sample_list_og = c(sample_list, "Mbuti")

l =read_tsv("global_pca.list", col_names=F)

data_loc = '/mnt/ancient/ModernHuman/ReichLab/reich_public_geno_v44.3/v44.3_HO_public'
inds = read_table2(glue("{data_loc}.ind"), col_names=c("iid","sex", "pop"))

f2loc1 = 'pca_test3'

f2s = admixtools::extract_f2(data_loc, pops=l$X1, outdir=f2loc1, maxmiss=0.03, maxmem=64000)

