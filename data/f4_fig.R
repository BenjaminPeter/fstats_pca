source("fscripts.R")
library(glue)
BTHEME = theme_classic() + theme(legend.position="none")

if(T){
f2s = admixtools::read_f2("westeurasian1/")
pcmat = pca_from_f2s(f2s)
pcs = pca_from_pcmat(pcmat)
idx = "Finnish"
idy = "Canary_Islander"
idz = "French"
ids = c(idx, idy)

sample_list = pcs$pop 

f4 = f4(f2s, pop3=idx, pop4 = idy, pop1=sample_list, pop2=sample_list) %>% arrange(pop1, pop2)
f4mat = f4 %>% 
    select(pop1, pop2, est) %>% 
    pivot_wider(names_from=pop2, values_from=est, values_fill=0, names_sort=T) %>% 
    column_to_rownames('pop1')
f4pc = f4_from_pc_matrix2(pcmat, idx, idy) %>%
    left_join(f4, .) %>% 
    filter(! pop1 %in% ids, ! pop2 %in% ids) %>%
    mutate(PC1_2=rowSums(across(PC1:PC2)), PC1_10=rowSums(across(PC1:PC10)))
v = pcmat[ids,] %>% diff
pmat = t(v) %*% v / (v %*% t(v))[1,1]
qmat = diag(length(v)) - pmat
a = pcmat %*% pmat
b = pcmat %*% qmat
A = t(t(a) - a[idz,])
B = t(t(b) - b[idz,])
E = eigen(b %*% t(b))
K = t(t(E$vectors) * sqrt(pmax(E$values, 0)))
K = K %>% as.data.frame %>% mutate(pop = rownames(B))

EA = eigen(a %*% t(a))
KA = t(t(EA$vectors) * sqrt(pmax(EA$values, 0)))
KA = KA %>% as.data.frame %>% mutate(pop = rownames(B))
X = KA %>% select(pop, A=V1) %>% left_join(K)

P3 = X %>% 
    ggplot(aes(x=A, y=V1, label=pop, color=V2)) + 
    geom_text(size=3) + 
    coord_fixed() + 
    BTHEME + 
    scale_color_viridis_c() + 
    xlab("proj (Finnish, Canary Islander)") + 
    ylab("Residual PC1")

f42 = f4(f2s, pop3=idx, pop4 = idy, pop1=sample_list, pop2=idz) %>% arrange(pop1, pop2)
f4pc2 = f4_from_pc_matrix(pcmat, idz, idx, idy) %>%
        as.data.frame %>% rownames_to_column("pop1") %>%
        left_join(f42, .) %>% 
        filter(! pop1 %in% ids) %>%
        mutate(PC1_2=rowSums(across(PC1:PC2)), 
               PC1_3=rowSums(across(PC1:PC3)),
               PC1_10=rowSums(across(PC1:PC10)))
    v = f4pc2 %>% 
        filter(pop1 %in% c("Greek", "Iranian", "Sicilian", "Estonian", "Saudi", "Georgian", "Spanish")) %>% 
        pivot_longer(PC1:PC33) %>% 
        mutate(PC=as.integer(substr(name, 3,20))) %>% 
        #group_by(pop1) %>% mutate(value=cumsum(value)) %>% 
        ungroup 
    P4= v %>%   ggplot(aes(x=PC, y=value, group=pop1, color=pop1, fill=pop1)) + 
        geom_col(position="dodge") + 
        facet_grid(pop1~.) + 
        geom_hline(aes(yintercept=-est), color="lightgray") + 
        geom_hline(yintercept=0) + BTHEME + 
        ylab(glue("f4(X, {idz}; {idx}, {idy})")) +
        coord_cartesian(xlim=c(1,10))

#pct var plot
x = c(EA$values[1], E$values[1:61]) %>% as.data.frame
names(x)[1] = "pctv"
x$lab = c("proj", sprintf("RPC%d", 1:61))
x$lab = factor(x$lab, levels=x$lab)
x$pctv = x$pctv / sum(x$pctv)
P5 = x[1:10,] %>% 
    ggplot(aes(y=pctv*100, x=lab)) + geom_col() + 
    xlab(NULL) + ylab("% variance") + BTHEME + theme(axis.text.x=element_text(angle=90, vjust=.5))

P6 = f4pc %>% 
    filter(!pop1 %in% ids, pop1 != pop2, pop1 != pop3, pop1 != pop4, pop2=="French") %>%
    select(pop3, est, PC1_2, PC1_10) %>% 
    pivot_longer(starts_with("PC"), names_to='PC', values_to="f4") %>% 
    ggplot(aes(x=-est, y=f4, color=PC)) + geom_abline() + 
    geom_hline(yintercept=0, color="lightgray") + 
    geom_vline(xintercept=0, color="lightgray") +
    geom_point() + 
    xlab("f4(X, FRA; FIN, CI)") + ylab("approx. f4")  + BTHEME

f4all = f4(f2s, pop3=idx, pop4 = idy, pop1=sample_list, pop2=sample_list) %>% arrange(pop1, pop2)
f4pcall = f4_from_pc_matrix2(pcmat, px=idx, py=idy) %>%
    #as.data.frame %>% rownames_to_column("pop1") %>%
    left_join(f4all, .) %>% 
    filter(! pop1 %in% ids) %>%
    mutate(PC1_2=rowSums(across(PC1:PC2)), PC1_4=rowSums(across(PC1:PC4)))



}
if(F){
    f2s = admixtools::read_f2("worldfoci2")
    pcmat = pca_from_f2s(f2s)
    pcs = pca_from_pcmat(pcmat)
    idx = "Sardinian"
    idy = "Mozabite"
    idz = "Yoruba"
    ids = c(idx, idy, idz)

    sample_list = pcs$pop 

    f2 = f2(f2s, unique_only=F) %>% select(pop1, pop2, f2=est)
    f4 = f4(f2s, pop3=idy, pop4 = idz, pop1=sample_list, pop2=idx) %>% arrange(pop1, pop2)
    f4mat = f4 %>% 
        select(pop1, pop2, est) %>% 
        pivot_wider(names_from=pop2, values_from=est, values_fill=0, names_sort=T) %>% 
        column_to_rownames('pop1')
    f4all = f4(f2s, pop3=idy, pop4 = idz, pop1=sample_list, pop2=sample_list) %>% arrange(pop1, pop2)
    f4pcall = f4_from_pc_matrix2(pcmat, px=idy, py=idz) %>%
        #as.data.frame %>% rownames_to_column("pop1") %>%
        left_join(f4all, .) %>% 
        filter(! pop1 %in% ids) %>%
        mutate(PC1_2=rowSums(across(PC1:PC2)), PC1_10=rowSums(across(PC1:PC20)))

    f4pc = f4_from_pc_matrix(pcmat, idx, idy, idz) %>%
        as.data.frame %>% rownames_to_column("pop1") %>%
        left_join(f4, .) %>% 
        filter(! pop1 %in% ids) %>%
        mutate(PC1_2=rowSums(across(PC1:PC2)), PC1_10=rowSums(across(PC1:PC20)))

    f2x = f2pc %>% 
        mutate(f2_2 = rowSums(across(PC1:PC2)), 
               f2_10 = rowSums(across(PC1:PC10))) %>% 
        select(pop1, pop2, f2_2, f2_10)

    f4x = f4pc %>% left_join(f2) %>% rename(f2.12=f2) %>% 
        left_join(f2, by=c(pop3="pop1", pop4="pop2")) %>% rename(f2.34=f2) %>% 
        mutate(angle=acos(-est / sqrt(f2.12 * f2.34)))  %>%
        left_join(f2x) %>% 
        left_join(f2x, by=c(pop3="pop1", pop4="pop2"), suffix =c(".12", ".34")) %>%
        mutate(angle_2 = acos(PC1_2 / sqrt (f2_2.12 * f2_2.34)),
               angle_10 = acos(pmin(pmax(PC1_10 / sqrt (f2_10.12 * f2_10.34), -1),1 )))

    P1 = f4x %>% select(pop1, starts_with("angle")) %>% 
        pivot_longer(angle:angle_10) %>% 
        ggplot(aes(x=value / pi * 180, y=pop1, color=name)) + 
        geom_vline(xintercept=c(0, 90), color='lightgray') + 
        geom_point() + BTHEME + 
        xlab("angle (degrees)") + ylab(NULL) + 
        scale_x_continuous(breaks=c(0, 45, 90, 135, 180))  

    P2 = f4pc %>% 
        filter(pop1 %in% c("Han", "Basque", "Papuan", "Mbuti", "Surui", "AA", "GujaratiD")) %>% 
        pivot_longer(PC1:PC33) %>% 
        mutate(PC=as.integer(substr(name, 3,20))) %>% 
        #group_by(pop1) %>% mutate(value=cumsum(value)) %>% 
        ggplot(aes(x=PC, y=value, group=pop1, color=pop1, fill=pop1)) + 
        geom_col(position="dodge") + 
        facet_grid(pop1~.) + 
        geom_hline(aes(yintercept=-est), color="lightgray") + 
        geom_hline(yintercept=0) + BTHEME + 
        ylab("f4(X, Sardinian; Mozabite, Yoruba)") +
        coord_cartesian(xlim=c(1,10))

}

if(F){
f2s = admixtools::read_f2("f4ratio"); sample_list= rownames(f2s)[!rownames(f2s) %in% "Denisova.DG"]
f2s = admixtools::read_f2("f4ratio", pops=sample_list )
pcmat = pca_from_f2s(f2s)
pcs = pca_from_pcmat(pcmat)
#idx = "Primate_Chimp"
#idz = "Dinka.DG"
#idy = "Altai_Neanderthal.DG"
idx = "Primate_Chimp"
idy = "Yoruba"
idz = "Altai_Neanderthal.DG"
ids = c(idx, idy)

sample_list = pcs$pop 

v = pcmat[ids,] %>% diff
pmat = t(v) %*% v / (v %*% t(v))[1,1]
qmat = diag(length(v)) - pmat
a = pcmat %*% pmat
b = pcmat %*% qmat
A = t(t(a) - a[idz,])
B = t(t(b) - b[idz,])
E = eigen(B %*% t(B))
K = t(t(E$vectors) * sqrt(pmax(E$values, 0)))
K = K %>% as.data.frame %>% mutate(pop = rownames(B))

E = eigen(A %*% t(A))
KA = t(t(E$vectors) * sqrt(pmax(E$values, 0)))
KA = KA %>% as.data.frame %>% mutate(pop = rownames(B))
X = KA %>% select(pop, A=V1) %>% left_join(K)


f4 = f4(f2s, pop3=idx, pop4 = idy, pop1=sample_list, pop2=sample_list) %>% arrange(pop1, pop2)
f4mat = f4 %>% 
    select(pop1, pop2, est) %>% 
    pivot_wider(names_from=pop2, values_from=est, values_fill=0, names_sort=T) %>% 
    column_to_rownames('pop1')
f4pc = f4_from_pc_matrix2(pcmat, idx, idy) %>%
    left_join(f4, .) %>% 
    filter(! pop1 %in% ids, ! pop2 %in% ids) %>%
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


ggsave( "f4_westeurasia1.png", P3, width=2.75, height=2.8, scale=1.5)
ggsave( "f4_westeurasia2.png", P4, width=1.25, height=4, scale=1.5)
ggsave( "f4_westeurasia3.png", P5, width=1.3, height=1, scale=1.5)
ggsave( "f4_westeurasia4.png", P6, width=1.3, height=1, scale=1.5)
ggsave( "f4_world1.png", P1, width=1.75, height=4, scale=1.5)
ggsave( "f4_world2.png", P2, width=1.25, height=4, scale=1.5)




