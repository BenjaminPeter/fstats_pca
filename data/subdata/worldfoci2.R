require(tidyverse)
system("ln -s worldfoci2.geno worldfoci2ind.geno")
system("ln -s worldfoci2.snp worldfoci2ind.snp")
x = read.table("worldfoci2.ind")
x = x %>% group_by(V3) %>% mutate(V3=sprintf("%s%d", V3, order(V1)))
write.table(x, "worldfoci2ind.ind", quote=F, col.names=F, row.names=F)

