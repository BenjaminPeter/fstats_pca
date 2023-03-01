require(ape)

double_center <- function(dist)
    t(dist - rowMeans(dist)) - rowMeans(dist) + mean(dist)
dc <- double_center

u_center <- function(dist){
    n = ncol(dist)
    r = rowSums(dist) / (n-2)
    dist = t(dist - r ) - r + sum(dist) / (n-1) / (n-2) 
    diag(dist) = 0
    return(dist)
}
uc = u_center

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
#dist = dc(dist)

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
q_mat1 = function(dist){
    n = nrow(dist)
    mv = rowMeans(dist)
    #M = n * (dist - outer(mv, mv, `+`) ) - 2 * dist
    M =  dc(dist) - 2/n * dist - sum(dist) / n / n
    return (M * n)
}

#! get next pair to merge
get_next_pair <- function(dist){
    Q = q_mat(dist)
    diag(Q) <- Inf
    pair = c(arrayInd(which.min(Q), dim(Q)))
    return(pair)
}
    

nj_dist_update <- function(dist, ij){
    i = ij[1]; j = ij[2]
    n = nrow(dist)
    #Si = sum(dist[,i]);  Sj = sum(dist[,j])

    #di_new = (dist[i,j] + (Si - Sj) / (n - 2)) / 2
    #dj_new = (dist[i,j] + (Sj - Si) / (n - 2)) / 2
    dist[,i]/2 + dist[,j]/2 - dist[i, j] / 2
}

local_pca <- function(dist, ij){
    
}

mynj <- function(dist, labels=NULL){
    if(is.null(labels)) labels=rownames(dist)
    if(is.null(labels)) labels=1:nrow(dist)
    tree <- list()
    tree$edge = matrix(NA, nrow=2 * length(labels) - 3 , ncol=2)
    tree$edge.length = rep(NA, 2 * length(labels) - 3)
    tree$tip.label = labels
    tree$Nnode = length(labels) - 2
    tree$dist = list()
    n_tips = length(tree$tip.label)

    active_nodes = 1:n_tips


    for (i in 1:tree$Nnode){
        pair = get_next_pair(dist)
        new_node = n_tips + i
        tree$edge[2 * i - 1,] = c(new_node, active_nodes[pair[1]])
        tree$edge[2 * i + 0,] = c(new_node, active_nodes[pair[2]])
        tree$edge.length[ 2 * i - 1 ] = 1 #placeholder
        tree$edge.length[ 2 * i + 0 ] = 1 #placeholder

        active_nodes = c( active_nodes[-pair], new_node)
        new_dists =  nj_dist_update(dist, pair)[-pair]
        new_dmat = cbind(dist[-pair, -pair], new_dists)
        dist = rbind(new_dmat, c(new_dists, 0))
        rownames(dist) = active_nodes
        colnames(dist) = active_nodes
        dist = dc(dist)
        tree$dist[[i]] = dist
    }
    tree$edge[2 * i + 1,] = active_nodes[2:1]
    tree$edge.length[ 2 * i + 1 ] = 1 #placeholder
    
    tree$edge[tree$edge > n_tips] = 3 * n_tips - tree$edge[tree$edge > n_tips] - 1
    class(tree) = 'phylo'

    
    return(tree)

}

