require(ape)

double_center <- function(dist)
    t(dist - rowMeans(dist)) - rowMeans(dist) + mean(dist)
dc <- double_center

rdist <- function(n, dist){
}

#helper function for display
csf <- function(...)cumsum(...)[c(1:10, length(...))]
f2 <- function(data, S1, S2){(data[S1,] - data[S2,])^2 }
f3 <- function(data, S1, S2, S3){(data[S1,] - data[S2,]) * (data[S1,] - data[S3,])}
f4 <- function(data, S1, S2, S3, S4){(data[S1,] - data[S2,]) * (data[S3,] - data[S4,])}



require(glue)
N_SNPS = 1000
N_SAMPLES = 10
TRUE_RANK = 2
S1 = 1; S2 = 2
S3 = 4; S4 = 3
NOISE = 0.15


require(MASS)
set.seed(12)
points = matrix(runif(N_SAMPLES * TRUE_RANK), nrow=N_SAMPLES)
points <- rbind(points)
test_cov = points %*% t(points)
data = mvrnorm(N_SNPS, mu=rep(0, N_SAMPLES), Sigma=test_cov); 
d2=data + mvrnorm(N_SNPS, mu=rep(0, N_SAMPLES), 
                  Sigma=diag(NOISE, N_SAMPLES))


dist = as.matrix(dist((points)))^2

#d1 = cophenetic.phylo(rtree(10))
#d2 = cophenetic.phylo(rtree(10))
#dist = d1 + d2


require(phangorn)

#eq 1
calc_q_ij <- function(dist, i, j){
    m <- nrow(dist)
    (m-2)*dist[i, j] - sum(dist[,i]) - sum(dist[,j])
}
q_mat = function(dist){
    M = (nrow(dist)-2) * dist - outer(rowSums(dist), colSums(dist), `+`) 
    return (M)
}

#! get next pair to merge
get_next_pair <- function(dist){
    Q = q_mat(dist)
    diag(Q) <- Inf
    pair = c(arrayInd(which.min(Q), dim(Q)))
}
    

nj_dist_update <- function(dist, ij){
    i = ij[1]; j = ij[2]
    n = nrow(dist)
    Si = sum(dist[,i]);  Sj = sum(dist[,j])

    di_new = (dist[i,j] + (Si - Sj) / (n - 2)) / 2
    dj_new = (dist[i,j] + (Sj - Si) / (n - 2)) / 2
    dist_new = dist[,i]/2 + dist[,j]/2 - dist[i, j] / 2
}

local_pca <- function(dist, ij){
    
}

mynj <- function(dist, pairs=c(), labels=NULL){
    if(is.null(labels)) labels=rownames(dist)
    if(is.null(labels)) labels=1:nrow(dist)
    if(dim(dist)[1] == 1) return(list(pairs, labels))

    pair = get_next_pair(dist)
    pairs = rbind(pairs, labels[pair])
    labels = c(labels, paste(labels[pair],collapse=''))[-pair]

    new_dists =  nj_dist_update(dist, pair)[-pair]
    reduced_dmat = dist[-pair, -pair]
    new_dmat = cbind(reduced_dmat, new_dists)
    new_dmat = rbind(new_dmat, c(new_dists, 0))
    mynj(new_dmat, pairs=pairs, labels=labels)
}

