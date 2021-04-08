require(tidyverse)

IND = 'I0104'

#ibd_all = read_tsv("./ibd.pairs_all.v1.addinfo.tsv") 
ibd = ibd_all %>%
    filter(iid1 == IND | iid2 == IND) %>% 
    mutate(iid=ifelse(iid1==IND, iid2, iid1)) %>% 
    select(-iid1, -iid2) %>% relocate(iid)

a = read_csv("~/fstats_pca/data/test_I0104.bin.xz", col_types=cols(chrom='f'))

b = a %>% select(chrom, map, YAM:YAMCOW) %>% pivot_longer(cols=YAM:YAMCOW) %>% 
    filter(value>.2, name != "YAM")
