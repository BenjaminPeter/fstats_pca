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
    return(list(W=W, sigma=sigma, eval=e$values, evec=e$vectors))
}


test_ppca <- function(){
    require(MASS)
    set.seed(10)
    v1 = runif(100)
    v2 = runif(100)
    test_cov = v1 %*% t(v1) + v2 %*% t(v2) #rank 2
    data = mvrnorm(1000, mu=rep(0, nrow(test_cov)), Sigma=test_cov); 
    d2=data + mvrnorm(1000, mu=rep(0, nrow(test_cov)), Sigma=diag(nrow(test_cov))/10)

    x = ppca(d2, verbose=F, n_pcs=2)
    y = ppca(d2, verbose=F, n_pcs=10)
    z = prcomp(d2)

    print("PPCA 2 PCs")
    print(x$eval[1:10]) 
    print("PPCA 10 PCs")
    print(y$eval[1:10]) 
    print("PCA first 10 PCs")
    print(z$sdev[1:10]) 
    
        
    
}

ppca2 = function(Y){
    #formula from agrawal
    C = diag(1, ncol=4, nrow=100)
    X = solve(t(C) %*% C) %*% t(C) %*% t(Y)
    C = t(Y) %*% t(X) %*% solve(X %*% t(X))
}
