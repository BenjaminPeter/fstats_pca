require(admixtools)
require(tidyverse)
require(glue)


data_loc = 'subdata/westeurasian1ind'
inds = read_table2(glue("{data_loc}.ind"), col_names=c("iid","sex", "pop"))



f2loc1 = 'westeurasian1ind'

#afs = admixtools::extract_afs(data_loc, pops=sample_list, 
#                              outdir=glue("{f2loc1}_afs"),format="eigenstrat", 
#                              cols_per_chunk=1000)
f2s = admixtools::extract_f2(data_loc, 
                             inds=inds$iid,
                             outdir=f2loc1, maxmem=500000, 
                             n_cores=140,
                             overwrite=T)

