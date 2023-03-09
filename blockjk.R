require(MASS)
set.seed(2)

N_SNPS = 100000
N_SAMPLES = 7
TRUE_RANK = 4
REP_SIZE = 10
BLOCK_SIZE = 21


points = matrix(runif(N_SAMPLES * TRUE_RANK), nrow=N_SAMPLES)
points <- rbind(points)

C = points %*% t(points)
D = outer(diag(C),diag(C),`+`) - 2 * C

#equivalent constructions
dist = as.matrix(dist((points)))^2
all.equal(c(D), c(dist))

#good data X0, and fake data X with replicated SNPs
X0 = mvrnorm(N_SNPS, rep(0, N_SAMPLES), C)
X = X0[rep(1:(N_SNPS), each=REP_SIZE),]


N0 = nrow(X0)
N = nrow(X)


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

ess <- function(x, blocks){
    return(1 / length(x) /  block_jk(x, blocks) ^2)
}


#let's look at the statistic F2(1, 2)
x_1 = f2_x[1,2]
x0_1 = f2_x0[1,2]
p0_1 = sum((P0[,1] - P0[,2])^2) / N0
p_1 = sum((P[,1] - P[,2])^2) / N


#we can predict SE analytically
se_true = sqrt(2*dist[1,2]^2 / N0)

#standard errors are biased for replicated data
se_1 = sqrt(var((X[,1]-X[,2])^2)/N) 
se0_1 = sqrt(var((X0[,1]-X0[,2])^2)/N0)

#jk errors are better
se_2 = block_jk((X0[,1]-X0[,2])^2,BLOCK_SIZE)
se0_2 =block_jk((X[,1]-X[,2])^2,BLOCK_SIZE)

#we can calculate inflation factor
i_1 = se_2 / se_1
i0_1 = se0_2 / se0_1

#this seems to match, (but not sure where square comes from)
Ne0_1 = N0 / i0_1^2
Ne_1 = N / (i_1^2)


#can we also obtain inflation factor from PCA?
i_2 = apply(L, 2, function(i)block_jk(i, 100) / sqrt(var(i) / length(i)))
i0_2 = apply(L0, 2, function(i)block_jk(i, 100) / sqrt(var(i) / length(i)))
Ne0_2 = N0 / i0_2^2
Ne_2 = N / (i_2^2)

#in fact, since the columns of L are normalized, we can do better:
i_3 = apply(L, 2, function(i)block_jk(i, 100) * length(i))
i0_3 = apply(L0, 2, function(i)block_jk(i, 100) * length(i))
Ne0_3 = N0 / i0_3^2
Ne_3 = N / (i_3^2)

i_4 = apply(L, 2, function(i)block_jk(i, 100) )
i0_4 = apply(L0, 2, function(i)block_jk(i, 100) )
Ne0_4 = 1 / N0 /  i0_4^2
Ne_4 = 1 / N /  i_4^2


#factor is reasonably similar across PCs, a natural summary is to weight them by
#associated eigenvalues
pcw = weighted.mean(i_2, D)
pcw_0 = weighted.mean(i0_2, D0)

#and again getting "effective sizes" requires squares for some reason
Ne0_w = N0 / pcw_0^2
Ne_w = N / pcw^2




data_test <- function(){
    G = eigenstrat('~/fstats_pca/data/subdata/worldfoci2') %>% read_geno
    G[is.na(G)] = 0
    G = as.matrix(G)
    H = colMeans(G * (2-G))
    C = (t(G) %*% G - diag(H)) / nrow(G)
    s = svd(C)
    P = s$u %*% diag(s$d)

    UV = svd(G)
}
