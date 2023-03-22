require(MASS)
set.seed(4)

N_SNPS = 100000
N_SAMPLES = 50
TRUE_RANK = 2
REP_SIZE = 8 
BLOCK_SIZE = 21
N_PCS=3

cov2dist <- function(C)outer(diag(C),diag(C),`+`) - 2 * C


points = matrix(runif(N_SAMPLES * TRUE_RANK), nrow=N_SAMPLES)
points <- rbind(points)

C = points %*% t(points)
D_pred = cov2dist(C)

#equivalent constructions
dist = as.matrix(dist((points)))^2
all.equal(c(D_pred), c(dist))

#good data X0, and fake data X with replicated SNPs
X0 = mvrnorm(N_SNPS, runif(N_SAMPLES), C)
#X0 = mvrnorm(N_SNPS, rep(0, N_SAMPLES), C)
rd = diff(range(X0))
X0= (X0 - min(X0)) / rd
X = X0[rep(1:(N_SNPS), each=REP_SIZE),]

N0 = nrow(X0)
N = nrow(X)


#same with genotype data
G0 = matrix(rbinom(prod(dim(X0)), 2, X0), nrow=dim(X0))  / 2
G = G0[rep(1:(N_SNPS), each=REP_SIZE),]
H0 =  colMeans(G0 * (1-G0))
H =  colMeans(G * (1-G))

CG0 = t(G0) %*% G0 / N0 - diag(H0)
CG = t(G) %*% G / N - diag(H)

f2_g0 = cov2dist(CG0)
f2_g = cov2dist(CG)

EG0 = eigen(CG0)
EG = eigen(CG)

#predicted data, using equation from p4 of van Waaij et al.
R0 = (G0 %*% EG0$vec[,1:N_PCS] %*% t(EG0$vec[,1:N_PCS]))
R = (G %*% EG$vec[,1:N_PCS] %*% t(EG$vec[,1:N_PCS]))

SVD_R0 = svd(R0)
SVD_R = svd(R)



#mean estimates are unaffected
f2_x = as.matrix(dist(t(X)))^2 / N 
f2_x0 = as.matrix(dist(t(X0)))^2 / N0


SVD = svd(X)
P = diag(SVD$d) %*% t(SVD$v)
L = SVD$u
D = SVD$d

SVD0 = svd(X0)
P0 = diag(SVD0$d) %*% t(SVD0$v)
L0 = SVD0$u
D0 = SVD0$d

#test for equality
all.equal(c(L %*% P), c(X))
all.equal(c(L0 %*% P0), c(X0))


#bjk standard error of a 1d-vector x.
#blocks has one of 3 formats:
# 1. a single number, giving average block size (will be rounded to fit data)
# 2. a vector of length n, giving block for each bin
# 3. a vector of sum n, giving size of each block
block_jk <- function(x, blocks){
    n = length(x)

    #goal is to create a vector blocks with assignments, and a vector
    #m having block lengths
    if(length(blocks) == 1){ 
        blocks = ceiling( 1:n / blocks)
    } 
    if(length(blocks) == length(x)){ 
        m = tabulate(blocks)
    } 
    if(sum(blocks) == length(x)){ 
        m = blocks
        nb = length(m)
        blocks = rep(1:nb, m)
    } 

    pseudo_val = c(rowsum(x, blocks))
    sqrt(mean((pseudo_val / m - mean(x))^2 / (n-m) * m ))
}

se <- function(x){
    sqrt(var(x) / length(x))
}


bj_l = apply(L, 2, function(i)block_jk(i, 100)) #normal PCA, with replicates
bj_l0 = apply(L0, 2, function(i)block_jk(i, 100)) #normal PCA, no replicates
bj_pl0 = apply(SVD_R0$u, 2, function(i)block_jk(i, 100)) #PPCA, replicates
bj_pl = apply(SVD_R$u, 2, function(i)block_jk(i, 100)) #PPCA, no replicates


df = data.frame(bj, bj0, l, l0, pc=1:length(bj))
df %>% pivot_longer(-pc) %>% filter(pc<=10) %>% ggplot(aes(x=pc, y=value, col=name)) + geom_point()

#factor is reasonably similar across PCs, a natural summary is to weight them by
#associated eigenvalues
iL = apply(L, 2, function(i)block_jk(i, 100) / sqrt(var(i) / length(i)))
iL0 = apply(L0, 2, function(i)block_jk(i, 100) / sqrt(var(i) / length(i)))
iPL0 = apply(SVD_R0$u, 2, function(i)block_jk(i, 100) / sqrt(var(i) / length(i)))
iPL = apply(SVD_R$u, 2, function(i)block_jk(i, 100) / sqrt(var(i) / length(i)))

pcwL = weighted.mean(iL, D)
pcwL0 = weighted.mean(iL0, D0)
pcwPL = weighted.mean(iPL, SVD_R$d)
pcwPL0 = weighted.mean(iPL0, SVD_R0$d)

#and again getting "effective sizes" requires squares for some reason
NeL0 = N0 / pcwL0^2
NeL = N / pcwL^2
NePL0 = N0 / pcwPL0^2
NePL = N / pcwPL^2
