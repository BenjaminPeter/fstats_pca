require(ape)

tplot <- function(..., label){
    plot(..., type='n')
    text(..., label=label)
}

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


#helper function for display
csf <- function(...)cumsum(...)[c(1:10, length(...))]
f2 <- function(data, S1, S2){(data[S1,] - data[S2,])^2 }
f3 <- function(data, S1, S2, S3){(data[S1,] - data[S2,]) * (data[S1,] - data[S3,])}
f4 <- function(data, S1, S2, S3, S4){(data[S1,] - data[S2,]) * (data[S3,] - data[S4,])}



require(glue)
N_SNPS = 1000
N_SAMPLES = 12
TRUE_RANK = 1
S1 = 1; S2 = 2
S3 = 4; S4 = 3
NOISE = 0.15


require(MASS)
set.seed(1)
#points = matrix(runif(N_SAMPLES * TRUE_RANK), nrow=N_SAMPLES)
#points <- rbind(points)
#dist = as.matrix(dist((points)))^2
#test_cov = points %*% t(points)
#data = mvrnorm(N_SNPS, mu=rep(0, N_SAMPLES), Sigma=test_cov); 
#d2=data + mvrnorm(N_SNPS, mu=rep(0, N_SAMPLES), 
#                  Sigma=diag(NOISE, N_SAMPLES))


random_tree <- function(N_SAMPLES=10, seed=1){
    TREE = rtree(N_SAMPLES, rooted=F)
    TREE$tip.label = 1:N_SAMPLES
    return(TREE)
}
TREE = random_tree(N_SAMPLES, seed=1)

dist = cophenetic(TREE)
U = -uc(dist)
D = -dc(dist)
eU = eigen(U)
vU = eU$vectors
eU = eU$values


require(phangorn)

#eq 1
calc_q_ij <- function(dist, i, j){
    m <- nrow(dist)
    (m-2)*dist[i, j] - sum(dist[,i]) - sum(dist[,j])
}
q_mat = function(dist){
    #minimum Q is best fit
    M = (nrow(dist)-2) * dist - outer(rowSums(dist), colSums(dist), `+`) 
    diag(M) = 0
    return (M)
}


umerge <- function(U, x, y){
    U_new = (U[x,] + U[y,] - U[x,y]) / 2
}

unext <- function(U, x, y){
    xy = c(x, y)
    u_new = umerge(U, x, y)
    u2 = rbind(cbind(U[-xy, -xy], u_new[-xy]), c(u_new[-xy], 0))
    return(u2)
}


q_mat1 = function(dist){
    n = nrow(dist)
    mv = rowMeans(dist)
    #M = n * (dist - outer(mv, mv, `+`) ) - 2 * dist
    M =  dc(dist) - 2/n * dist - sum(dist) / n / n
    return (M * n)
}

#! get next pair to merge
get_next_pair <- function(dist, tol=1e-8){
    Q = q_mat(dist)
    diag(Q) <- NA
    Q[abs(Q) < tol] <- NA

    pair = c(arrayInd(which.min(Q), dim(Q)))
    return(pair)
}

get_next_pair_u <- function(dist, tol=1e-8){
    Q = dist 
    diag(Q) <- NA
    Q[abs(Q) < tol] <- NA

    pair = c(arrayInd(which.max(abs(Q)), dim(Q)))
    return(pair)
}


givens <- function(M, i, j){
    theta = atan( 2 * M[i,j] / (M[i,i] - M[j,j])) / 2
    ij = c(i,j)
    G = diag(nrow(M))
    G[i,i] = cos(theta) 
    G[j,j] = cos(theta) 
    G[i, j]= sin(theta) 
    G[j, i]= -sin(theta) 
    return(G)

}
    

nj_dist_update <- function(dist, ij){
    i = ij[1]; j = ij[2]
    n = nrow(dist)
    #Si = sum(dist[,i]);  Sj = sum(dist[,j])

    #di_new = (dist[i,j] + (Si - Sj) / (n - 2)) / 2
    #dj_new = (dist[i,j] + (Sj - Si) / (n - 2)) / 2
    dist[,i]/2 + dist[,j]/2 - dist[i, j] / 2
}


test_jacobi <- function(U){
    n = nrow(dist)
    V = diag(n)
    for(i in 1:100){
        pair = get_next_pair_u(U)
        print(sprintf('[%d]: %d -> %d: %e', i, pair[1], pair[2], U[pair[1], pair[2]]))
        G = givens(U, pair[1], pair[2])
        U = G %*% U %*% t(G)
        V = V %*% G
    }
    return(list(U=U, V=V))

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

