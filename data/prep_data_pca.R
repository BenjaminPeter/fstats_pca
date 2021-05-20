require(admixtools)
require(tidyverse)
require(glue)


MISSING = c("Irish_Ulster", "Polish", "German", "Irish", "Shetlandic", "Sorb")

list1 = read_table("WestEurasian_PCA.list", col_names=F)[,1] %>% unlist %>% 
    unname
sample_list = list1[!list1 %in% MISSING]
sample_list_og = c(sample_list, "Mbuti")

data_loc = '/mnt/ancient/ModernHuman/ReichLab/reich_public_geno_v44.3/v44.3_HO_public'
inds = read_table2(glue("{data_loc}.ind"), col_names=c("iid","sex", "pop"))

#data_loc = '/mnt/ancient/ModernHuman/ReichLab/reich_public_geno_v44.3/v44.3_1240K_public'
#inds = read_table2(glue("{data_loc}.ind"), col_names=c("iid","sex", "pop"))


f2loc1 = 'pca_test1'
f2loc1 = 'pca_test1c'

afs = admixtools::extract_afs(data_loc, pops=sample_list_og, 
                              outdir=f2loc1,format="eigenstrat", 
                              cols_per_chunk=1000)
#f2s = admixtools::extract_f2(data_loc, pops=sample_list_og, 
#                             outdir=f2loc1, maxmem=64000, blgsize=.5,
                            # format="eigenstrat")

