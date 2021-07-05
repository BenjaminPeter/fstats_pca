source("fscripts.R")
library(glue)
BTHEME = theme_classic() + theme(legend.position="none")

if(F){
#f2s = admixtools::read_f2("f4ratio"); sample_list= rownames(f2s)[!rownames(f2s) %in% "Denisova.DG"]
f2s = admixtools::read_f2("f4ratio", pops=sample_list )
pcmat = pca_from_f2s(f2s)
pcs = pca_from_pcmat(pcmat)
idx = "Primate_Chimp"
#idz = "Dinka.DG"
idy = "Altai_Neanderthal.DG"
idy = "Altai_Neanderthal.DG"

#idx = "Primate_Chimp"
#idy = "Yoruba"
#idz = "Altai_Neanderthal.DG"
ids = c(idx, idy)

sample_list = pcs$pop 

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



f4 = f4(f2s, pop3=idx, pop4 = idy, pop1=sample_list, pop2=sample_list) %>% arrange(pop1, pop2)
f4mat = f4 %>% 
    select(pop1, pop2, est) %>% 
    pivot_wider(names_from=pop2, values_from=est, values_fill=0, names_sort=T) %>% 
    column_to_rownames('pop1')
f4pc = f4_from_pc_matrix2(pcmat, idx, idy) %>%
    left_join(f4, .) %>% 
    filter(! pop1 %in% ids, ! pop2 %in% ids) %>%
    mutate(PC1_2=rowSums(across(PC1:PC2)), PC1_3=rowSums(across(PC1:PC3)))

f4a = f4_from_pc_matrix2(a, idx, idy)  %>% 
    left_join(f4, .) %>% 
    filter(! pop1 %in% ids, ! pop2 %in% ids) %>%
    mutate(PC1_2=rowSums(across(PC1:PC2)), PC1_3=rowSums(across(PC1:PC3)))
f4b = f4_from_pc_matrix2(b, idx, idy)  %>% 
    left_join(f4, .) %>% 
    filter(! pop1 %in% ids, ! pop2 %in% ids) %>%
    mutate(PC1_2=rowSums(across(PC1:PC2)), PC1_3=rowSums(across(PC1:PC3)))



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


P4 = f3pc %>% select(pop3, est, PC1_2, PC1_3) %>% 
filter(!pop3 %in% c(idx, idy)) %>%
pivot_longer(starts_with("PC"), names_to='PC', values_to="f3") %>% 
ggplot(aes(x=est, y=f3, color=PC)) + geom_point() + geom_abline() + 
xlab("f3(Mbuti, Mozabite, X)") + ylab("approx. f3")  + BTHEME
P6 = f3pc %>% pivot_longer(PC1:PC3) %>% 
mutate(PC=as.integer(substr(name, 3,10))) %>% 
ggplot(aes(x=PC, y=value, group=PC)) + 
geom_hline(yintercept=0, color='lightgray') + 
geom_boxplot() + 
ylab("f3(Mbu.; Moz., X)") + 
scale_x_continuous(breaks=1:10) + BTHEME
}

if(T){
    f2s = admixtools::read_f2("f4ratio")
    rn = !rownames(f2s) %in% "Primate_Chimp"
    f2s = f2s[rn,,]
    pcmat = pca_from_f2s(f2s)
    pcs = pca_from_pcmat(pcmat)

    idx = "Primate_Chimp"
    idx = "Denisova.DG"
    #idy = "Altai_Neanderthal.DG"
    idy = "Yoruba"
    ids = c(idx, idy)
    idz = "Dinka.DG"

    sample_list = pcs$pop 

    f3 = f3(f2s, pop1=idx, pop2=idy, pop3=sample_list) #outgroup
    f3pc = of3_from_pc_matrix(pcmat, idx, idy) %>% 
        as.data.frame %>% rownames_to_column("pop3") %>% 
        left_join(f3, .) %>% 
        mutate(PC1_2=rowSums(across(PC1:PC2)), PC1_10=rowSums(across(PC1:PC10)))


    #get (A-B)
    v = pcmat[ids,] %>% diff

    #projection matrix P on (A-B)
    pmat = t(v) %*% v / (v %*% t(v))[1,1]

    #residual of the projection P
    qmat = diag(length(v)) - pmat

    # location of samples on the projection
    a = pcmat %*% pmat

    # location of samples in null-space
    b = pcmat %*% qmat
    A = t(t(a) - a[idz,])
    B = t(t(b) - b[idz,])
    E = eigen(b %*% t(b))
    K = t(t(E$vectors) * sqrt(pmax(E$values, 0)))
    K = K %>% as.data.frame %>% mutate(pop = rownames(B))

    EA = eigen(a %*% t(a))
    KA = t(t(EA$vectors) * sqrt(pmax(EA$values, 0)))

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
        geom_text_repel(aes(x=PC1, y=PC2, label=pop, color=est), data=pcs2) +
        geom_point(aes(x=PC1, y=PC2, color=est), data=pcs2) +
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

