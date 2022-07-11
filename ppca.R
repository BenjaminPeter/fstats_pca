#dummy inits
init_sigma <- function() {1}
init_cov = function(n_dims,n_pcs)diag(1, nrow=n_dims, ncol=n_pcs)


ppca <- function(data, n_pcs=4, n_iter=1000, tol=1e-4, verbose=T){
    n_dims = ncol(data)
    #eq 31
    centered_data = data - rowMeans(data)
    S = t(centered_data) %*% centered_data / nrow(data)

    sigma = init_sigma()
    W = init_cov(n_dims, n_pcs)
    
    score = Inf

    for(i in 1:n_iter){
        #M is defined after eq 6
        M = t(W) %*% W + diag(sigma, n_pcs)
        M_inv = solve(M)

        #eq 29
        aux = diag(sigma, n_pcs) + M_inv %*% t(W) %*% S %*% W
        W_new = S %*% W %*% solve(aux)

        #eq 30
        sigma_aux = S - S %*% W %*% M_inv %*% t(W_new)
        sigma = 1 / n_dims * sum(diag(sigma_aux))

        W = W_new


        #simple convergence test
        new_score = sum(abs(W)) + sigma
        if(verbose)print(sprintf("iter: %d, score : |%.5f|",i, score-new_score))
        if(abs(score - new_score) < tol){
            break
        } else {
            score = new_score
        }


    }

    e = eigen(W %*% t(W))
    pcs = t(sqrt(abs(e$values)) * t(e$vectors))
    return(list(W=W, sigma=sigma, pcs=pcs, eval=e$values, evec=e$vectors))
}


test_ppca <- function(){
    N_SNPS = 5000
    N_SAMPLES = 100
    TRUE_RANK = 2 
    S1 = 1; S2 = 19
    NOISE = 0.1


    require(MASS)
    set.seed(10)
    points = matrix(runif(N_SAMPLES * TRUE_RANK), nrow=N_SAMPLES)
    test_cov = points %*% t(points)
    data = mvrnorm(N_SNPS, mu=rep(0, N_SAMPLES), Sigma=test_cov); 
    d2=data + mvrnorm(N_SNPS, mu=rep(0, N_SAMPLES), Sigma=diag(NOISE, N_SAMPLES))

    x = ppca(d2, verbose=F, n_pcs=2)
    y = ppca(d2, verbose=F, n_pcs=10)
    z = prcomp(d2)
    z$pcs = unname(t(t(z$rotation) * z$sdev))

    print("PPCA 2 PCs")
    print(x$eval[1:10]) 
    print("PPCA 20 PCs")
    print(y$eval[1:10]) 
    print("PCA first 10 PCs")
    print(z$sdev[1:10]) 
    

    f2_true = as.matrix(dist(points))[S1,S2]^2
    f2_data = sum((d2[,S1] - d2[,S2])^2) / N_SNPS

        
    print("-------------------")
    print(sprintf("f2 from true points: %f", f2_true))
    print(sprintf("f2 from data: %f", f2_data))
    print("f2 from PPCA 2 PCs")
    print(cumsum((x$pcs[S1,] - x$pcs[S2,])^2)[c(1:10, N_SAMPLES)])
    print("f2 from PPCA 10 PCs")
    print(cumsum((y$pcs[S1,] - y$pcs[S2,])^2)[c(1:10, N_SAMPLES)])
    print("f2 from PCA first 10 PCs")
    print(cumsum((z$pcs[S1,] - z$pcs[S2,])^2)[c(1:10, N_SAMPLES)])
    
}

ppca2 = function(Y){
    #formula from agrawal
    C = diag(1, ncol=4, nrow=100)
    X = solve(t(C) %*% C) %*% t(C) %*% t(Y)
    C = t(Y) %*% t(X) %*% solve(X %*% t(X))
}
