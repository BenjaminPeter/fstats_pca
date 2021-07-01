require(tidyverse)
system("ln -s westeurasian1.geno westeurasian1ind.geno")
system("ln -s westeurasian1.snp westeurasian1ind.snp")
x = read.table("westeurasian1.ind")
x = x %>% group_by(V3) %>% mutate(V3=sprintf("%s%d", V3, order(V1)))
write.table(x, "westeurasian1ind.ind", quote=F, col.names=F, row.names=F)

