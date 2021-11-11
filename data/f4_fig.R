source("fscripts.R")
library(glue)
library(ggrepel)
BTHEME = theme_classic() + theme(legend.position="none")
R=1.5

if(T){
    f2s = admixtools::read_f2("worldfoci2")
    pcmat = pca_from_f2s(f2s)
    pcs = pca_from_pcmat(pcmat)
    idx = "Surui"
    idz = "Han"
    idy = "Mbuti"
ids = c(idx, idy)
n = nrow(pcs)

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
    mutate(V2 = ifelse(pop %in% c(idx, idy), NA, V2)) %>%
    ggplot(aes(x=A, y=V1, label=pop, color=V2)) + 
    #mutate(A = ifelse(pop %in% c(idx, idy), NA, A)) %>%
    #ggplot(aes(x=V1, y=V2, label=pop, color=A)) + 
    geom_hline(yintercept=X$V1[X$pop == idx], color='lightgray') +
    geom_text_repel() + 
    geom_point() +
    coord_fixed(xlim=c(-.1, .2), ylim=c(-.1, .1)) + BTHEME +
    scale_color_viridis_c(na.value='red') + 
    xlab(glue("⟨X; {idx}, {idy}⟩")) + 
    ylab("Residual PC1")

f42 = f4(f2s, pop3=idx, pop4 = idy, pop1=sample_list, pop2=idz) %>% arrange(pop1, pop2)
f4pc2 = f4_from_pc_matrix(pcmat, idz, idx, idy) %>%
        as.data.frame %>% rownames_to_column("pop1") %>%
        left_join(f42, .) %>% 
        filter(! pop1 %in% ids) %>%
        mutate(PC1_2=rowSums(across(PC1:PC2)), 
               PC1_3=rowSums(across(PC1:PC3)),
               PC1_5=rowSums(across(PC1:PC5)),
               PC1_10=rowSums(across(PC1:PC10)))
v = f4pc2 %>% 
    filter(pop1 %in% c("Masai", "Basque", "Papuan", "Mbuti", "Surui", "AA", "GujaratiD")) %>% 
    pivot_longer(PC1:PC33) %>% 
    mutate(PC=as.integer(substr(name, 3,20))) %>% 
    #group_by(pop1) %>% mutate(value=cumsum(value)) %>% 
    ungroup 
ids = c(idx, idy)
P4= v %>%   
    #group_by(pop1) %>% mutate(value=cumsum(value)) %>%
    ggplot(aes(x=PC, y=value, group=pop1, color=pop1, fill=pop1)) + 
        geom_col() + 
        facet_grid(. ~ pop1) + 
        geom_hline(aes(yintercept=-est), color="lightgray") + 
        geom_hline(yintercept=0) + BTHEME + 
        ylab(glue("F4(X, {idz}; {idx}, {idy})")) +
        coord_flip(xlim=c(.5,10.5)) +
        scale_x_continuous(breaks=1:10) +
        scale_y_continuous(breaks=c(-0.02, .02))


#pct var plot
x = c(EA$values[1], E$values[1:(n-1)]) %>% as.data.frame
names(x)[1] = "pctv"
x$lab = c("proj", sprintf("RPC%d", 1:(n-1)))
x$lab = factor(x$lab, levels=x$lab)
x$pctv = x$pctv / sum(x$pctv)
P5 = x[1:10,] %>% 
    ggplot(aes(y=pctv*100, x=lab)) + geom_col() + 
    xlab(NULL) + ylab("% variance") +
    scale_y_continuous(breaks=seq(0, 30, 5), lim=c(0,30)) + 
    BTHEME + theme(axis.text.x=element_text(angle=90, vjust=.5))

idx2=substr(idx, 1, 3)
idy2="YRI"
P6 = f4pc %>% 
    filter(!pop1 %in% ids, pop1 != pop2, pop1 != pop3, pop1 != pop4, pop2==idz) %>%
    select(pop3, est, PC1_2, PC1_10) %>% 
    pivot_longer(starts_with("PC"), names_to='PC', values_to="f4") %>% 
    ggplot(aes(x=-est, y=f4, color=PC)) + geom_abline() + 
    geom_hline(yintercept=0, color="lightgray") + 
    geom_vline(xintercept=0, color="lightgray") +
    geom_point() + 
    xlab(glue("F4(X, {idz}; {idx2}, {idy2})")) + ylab("approx. F4")  + BTHEME

f4all = f4(f2s, pop3=idx, pop4 = idy, pop1=sample_list, pop2=sample_list) %>% arrange(pop1, pop2)
f4pcall = f4_from_pc_matrix2(pcmat, px=idx, py=idy) %>%
    #as.data.frame %>% rownames_to_column("pop1") %>%
    left_join(f4all, .) %>% 
    filter(! pop1 %in% ids) %>%
    mutate(PC1_2=rowSums(across(PC1:PC2)), PC1_4=rowSums(across(PC1:PC4)))

e1 = eigen(pcmat %*% t(pcmat))$values 
#pct var plot
y = e1 %>% as.data.frame
names(y)[1] = "pctv"
y$lab = sprintf("PC%d", 1:(n))
y$lab = factor(y$lab, levels=y$lab)
y$pctv = y$pctv / sum(y$pctv)
P6 = y[1:10,] %>% 
    ggplot(aes(y=pctv*100, x=lab)) + geom_col() + 
    xlab(NULL) + ylab("% variance") +
    scale_y_continuous(breaks=seq(0, 30, 5), lim=c(0,41)) + 
    BTHEME + theme(axis.text.x=element_text(angle=90, vjust=.5))
P6 = cbind(lab=1:33, pctv=y$pctv, PCA=x[,1]) %>% 
    as_tibble %>% 
    pivot_longer(pctv:PCA) %>% group_by(name) %>% mutate(value=cumsum(value)) %>% 
    filter(lab<=10) %>% 
    ggplot(aes(x=lab, y=value, fill=name)) + 
    geom_col(position='dodge') + 
    BTHEME + ylim(0:1) + 
    scale_fill_manual(values=c("lightgray", "darkgray")) + 
    geom_hline(yintercept=c(0,1), color='lightgray') + 
    scale_x_continuous("Projection axis", breaks=1:10) + ylab("%var explained") 


ggsvg(P3,  "f4_world1.svg", width = R * 3.5, height=2.5 * R)
ggsvg(P4, "f4_world2.svg",  width = R * 4, height=1.5 *R)
ggsvg(P5, "f4_world3.svg", width = R * 1.3, height=1 *R)

ggsvg(P6,  "f4_world4.svg", width = R * 2.6, height=1.3 *R)
}
if(T){
f2s = admixtools::read_f2("westeurasian1/")
pcmat = pca_from_f2s(f2s)
pcs = pca_from_pcmat(pcmat)
idx = "Saudi"
idy = "Sardinian"
idz = "Finnish"
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
        mutate(PC1_2=rowSums(across(PC1:PC2)), 
               PC1_5=rowSums(across(PC1:PC5)),
               PC1_3=rowSums(across(PC1:PC3)))

    f4pc = f4_from_pc_matrix(pcmat, idx, idy, idz) %>%
        as.data.frame %>% rownames_to_column("pop1") %>%
        left_join(f4, .) %>% 
        filter(! pop1 %in% ids) %>%
        mutate(PC1_2=rowSums(across(PC1:PC2)), 
               PC1_5=rowSums(across(PC1:PC5)), 
               PC1_3=rowSums(across(PC1:PC3)))

    f2pc = f2_from_pc_matrix2(pcmat) %>%
        left_join(f2, .) %>% 
        mutate(PC1_2=rowSums(across(PC1:PC2)), 
               PC1_5=rowSums(across(PC1:PC5)),
               PC1_3=rowSums(across(PC1:PC3)))

    f2x = f2pc %>% 
        mutate(f2_2 = rowSums(across(PC1:PC2)), 
               f2_5 = rowSums(across(PC1:PC5)),
               f2_3 = rowSums(across(PC1:PC3))) %>% 
        select(pop1, pop2, f2_2, f2_3, f2_5)

    f4x = f4pc %>% left_join(f2) %>% rename(f2.12=f2) %>% 
        left_join(f2, by=c(pop3="pop1", pop4="pop2")) %>% rename(f2.34=f2) %>% 
        mutate(angle=acos(abs(-est / sqrt(f2.12 * f2.34))))  %>%
        left_join(f2x) %>% 
        left_join(f2x, by=c(pop3="pop1", pop4="pop2"), suffix =c(".12", ".34")) %>%
        mutate(angle_2 = acos(abs(PC1_2 / sqrt (f2_2.12 * f2_2.34))),
               #angle_5 = acos(abs(PC1_5 / sqrt (f2_5.12 * f2_5.34))),
               angle_3 = acos(abs(pmin(pmax(PC1_3 / sqrt (f2_3.12 * f2_3.34), -1),1 ))))

    P1 = f4x %>% select(pop1, starts_with("angle")) %>% 
        #filter(!startsWith(pop1, "Jew")) %>% 
        mutate(pop1 = substr(pop1, 1, 8)) %>% 
        mutate(pop1 = fct_reorder(pop1, angle)) %>% 
        pivot_longer(angle:angle_3) %>% 
        ggplot(aes(x=value / pi * 180, y=pop1, color=name)) + 
        geom_vline(xintercept=c(0, 90), color='lightgray') + 
        geom_point() + BTHEME + 
        xlab("angle (degrees)") + ylab(NULL) + 
        scale_x_continuous(breaks=seq(0, 180, 22.5), 
                           sec.axis=sec_axis(~ cos(pi * ./180), 
                                             breaks=c(.99, round(seq(-.9, .9, .3),1)), name="correlation (r)"))


    P2 = f4pc %>% 
        filter(pop1 %in% c("Greek", "Iranian", "Sicilian", "Estonian", "Canary_Islander", 
                           "French", "Spanish", 'Georgian')) %>% 
        mutate(pop1 = substr(pop1, 1, 8)) %>% 
        pivot_longer(PC1:PC10) %>% 
        mutate(PC=as.integer(substr(name, 3,20))) %>% 
        #group_by(pop1) %>% mutate(value=cumsum(value)) %>% 
        ggplot(aes(x=PC, y=value, group=pop1, color=pop1, fill=pop1)) + 
        geom_col() + 
        facet_grid(. ~ pop1) + 
        geom_hline(aes(yintercept=-est), color="lightgray") + 
        geom_hline(yintercept=0) + BTHEME + 
        ylab(glue("F4(X, {idx}; {idy}, {idz})")) +
        coord_flip(xlim=c(.5,10.5)) +
        scale_x_continuous(breaks=1:10) +
        scale_y_continuous(breaks=c(-0.001, .002))
ggsvg(P1, "f4_westeurasia1.svg", width = R * 1.5, height=4 *R)
ggsvg(P2, "f4_westeurasia2.svg", width = R * 5.5, height=1.5 *R)
}



