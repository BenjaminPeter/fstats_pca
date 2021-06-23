source("fscripts.R")
library(ggpubr)
BTHEME = theme_classic() + theme(legend.position="none")

if(T){
f2s = admixtools::read_f2("westeurasian1/")
pcmat = pca_from_f2s(f2s)
pcs = pca_from_pcmat(pcmat)
idx = "Basque"
idy = "Turkish"

sample_list = pcs$pop 

f3 = f3(f2s, pop3=idx, pop2=idy, pop1=sample_list) #admixture
f3pc = f3_from_pc_matrix(pcmat, idx, idy)
#f3 = f3(f2s, pop1=idx, pop2=idy, pop3=sample_list) #outgroup
f4 = f4(f2s, pop3=idx, pop4 = idy, pop1=sample_list, pop2=sample_list) %>% arrange(pop1, pop2)

circle = pcs %>% filter(pop %in% c(idx, idy)) %>% select(-pop)
circle_df = colMeans(circle) %>% t %>% as.data.frame
circle_df$radius = (dist(circle) %>% as.matrix)[1,2] / 2
circle_df$r0 = (dist(circle[,1:2]) %>% as.matrix)[1,2] / 2

pcs2 = pcs %>% left_join(f3, by=c(pop="pop1")) %>% 
    mutate(col=as.factor(cut(est, c(-1,-1e-10,1e-10,2), labels=F) ))
P1 = ggplot() + 
    geom_circle(aes(x0=PC1, y0 =PC2, r=radius), data=circle_df, color='lightgray', fill='lightgray') +
#    geom_circle(aes(x0=PC1, y0 =PC2, r=r0), data=circle_df, color='darkgray', fill='darkgray') +
    geom_text(aes(x=PC1, y=PC2, label=pop, color=col), data=pcs2) +
    coord_fixed()  + BTHEME

f3pc = f3_from_pc_matrix(pcmat, idx, idy) %>% 
    as.data.frame %>% rownames_to_column("pop1") %>% 
    left_join(f3, .) %>% 
    mutate(PC1_2=rowSums(across(PC1:PC2)), PC1_10=rowSums(across(PC1:PC10)))

P2 = f3pc %>% select(pop1, est, PC1_2, PC1_10) %>% 
    pivot_longer(starts_with("PC"), names_to='PC', values_to="f3") %>% 
    ggplot(aes(x=est, y=f3, color=PC)) + geom_point() + geom_abline() +
    xlab("f3(X; Turkish, Basque)") + ylab("approx. f3")  + BTHEME +
    coord_cartesian(xlim=c(-0.001, 0.005), ylim=c(-0.001, 0.005))

P5 = f3pc %>% pivot_longer(PC1:PC10) %>% 
    mutate(PC=as.integer(substr(name, 3,10))) %>% 
    ggplot(aes(x=PC, y=value, group=PC)) + 
    geom_hline(yintercept=0, color='lightgray') + 
    geom_boxplot() + 
    ylab("f3(X; Tur., Bas.)") + 
    scale_x_continuous(breaks=1:10) + BTHEME
}
if(T){
f2s = admixtools::read_f2("worldfoci2")
pcmat = pca_from_f2s(f2s)
pcs = pca_from_pcmat(pcmat)
idx = "Mbuti"
idy = "Mozabite"

sample_list = pcs$pop 

f3 = f3(f2s, pop1=idx, pop2=idy, pop3=sample_list) #outgroup
f3pc = of3_from_pc_matrix(pcmat, idx, idy) %>% 
    as.data.frame %>% rownames_to_column("pop3") %>% 
    left_join(f3, .) %>% 
    mutate(PC1_2=rowSums(across(PC1:PC2)), PC1_10=rowSums(across(PC1:PC10)))

#stuff for projection
smat = pcmat[c(idx, idy),]
slope = diff(smat[,2]) / diff(smat[,1])
intercept = smat[2,2] - slope * smat[2,1]
bv = diff(smat[,1:2])
norm = (bv %*% t(bv))[1,1]
bvall = diff(smat)
norm2 = (bvall %*% t(bvall))[1,1]
f3pc$k = (f3pc$PC1 + f3pc$PC2) / norm
f3pc$k2 = f3pc$est / norm2
df = tibble(x=pcs$PC1,
            y=pcs$PC2,
            xend=f3pc$k2 * bv[,1] + pcs$PC1[pcs$pop==idx], 
            yend=f3pc$k2 * bv[,2] + pcs$PC2[pcs$pop==idx])


pcs2 = pcs %>% left_join(f3, by=c(pop="pop3")) 
P3 = ggplot() + 
    geom_abline(intercept=intercept, slope=slope) +
    #geom_circle(aes(x0=PC1, y0 =PC2, r=radius), data=circle_df, color='lightgray', fill='lightgray') +
    #geom_circle(aes(x0=PC1, y0 =PC2, r=r0), data=circle_df, color='darkgray', fill='darkgray') +
    geom_text(aes(x=PC1, y=PC2, label=pop, color=est), data=pcs2) +
    scale_color_viridis_c(direction=1) + 
    geom_segment(data=df, aes(x=x, y=y, xend=xend, yend=yend), alpha=.3) +
    coord_fixed() + BTHEME


P4 = f3pc %>% select(pop3, est, PC1_2, PC1_10) %>% 
    filter(!pop3 %in% c(idx, idy)) %>%
    pivot_longer(starts_with("PC"), names_to='PC', values_to="f3") %>% 
    ggplot(aes(x=est, y=f3, color=PC)) + geom_point() + geom_abline() + 
    xlab("f3(Mbuti, Mozabite, X)") + ylab("approx. f3")  + BTHEME
P6 = f3pc %>% pivot_longer(PC1:PC10) %>% 
    mutate(PC=as.integer(substr(name, 3,10))) %>% 
    ggplot(aes(x=PC, y=value, group=PC)) + 
    geom_hline(yintercept=0, color='lightgray') + 
    geom_boxplot() + 
    ylab("f3(Mbu.; Moz., X)") + 
    scale_x_continuous(breaks=1:10) + BTHEME
}


ggsvg(P1, "f3_westeurasia1.svg", width=3.5, height=3)
ggsvg(P2, "f3_westeurasia2.svg", width=3.5, height=1)
ggsvg(P3, "f3_world1.svg", width=3.5, height=3)
ggsvg(P4, "f3_world2.svg", width=3.5, height=1)
ggsave( "f3_westeurasia1.png", P1, width=3.5, height=3, scale=1.5)
ggsave( "f3_westeurasia2.png", P2, width=1.75, height=1, scale=1.5)
ggsave( "f3_westeurasia3.png", P5, width=1.75, height=1, scale=1.5)
ggsave( "f3_world1.png", P3, width=3.5, height=3, scale=1.5)
ggsave( "f3_world2.png", P4, width=1.75, height=1, scale=1.5)
ggsave( "f3_world3.png", P6, width=1.75, height=1, scale=1.5)



