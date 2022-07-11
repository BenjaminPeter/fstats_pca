#dummy inits
init_sigma <- function() {1}
init_cov = function(n_dims,n_pcs)diag(1, nrow=n_dims, ncol=n_pcs)

#helper function for display
csf <- function(...)cumsum(...)[c(1:10, N_SAMPLES)]
f2 <- function(data, S1, S2){(data[S1,] - data[S2,])^2 }
f3 <- function(data, S1, S2, S3){(data[S1,] - data[S2,]) * (data[S1,] - data[S3,])}
f4 <- function(data, S1, S2, S3, S4){(data[S1,] - data[S2,]) * (data[S3,] - data[S4,])}


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
    S3 = 3; S4 = 6
    NOISE = 0.05


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

    print("PPCA 2 PCs eigenvalues")
    print(sqrt(x$eval[1:10])) 
    print("PPCA 20 PCs eigenvalues")
    print(sqrt(y$eval[1:10])) 
    print("PCA first 10 PCs eigenvalues")
    print(z$sdev[1:10]) 
    

    #test f2s
    f2_true = as.matrix(dist(points))[S1,S2]^2
    f2_data = sum((d2[,S1] - d2[,S2])^2) / N_SNPS


        
    print("-------------------")
    print(sprintf("f2 from true points: %f", f2_true))
    print(sprintf("f2 from data: %f", f2_data))
    print("f2 from PPCA 2 PCs")
    print(csf(f2(x$pcs, S1, S2)))
    print("f2 from PPCA 10 PCs")
    print(csf(f2(y$pcs, S1, S2)))
    print("f2 from PCA first 10 PCs")
    print(csf(f2(z$pcs, S1, S2)))

    par(mfrow=c(1,3))

    plot(NA, xlim=c(0, N_SAMPLES), ylim=c(0, max(f2_data, f2_true)))
    abline(h=c(f2_true, f2_data), col=2:3, lty=2)
    lines(cumsum(f2(x$pcs, S1, S2)),  col=4,lwd=2)
    lines(cumsum(f2(y$pcs, S1, S2)),  col=5)
    lines(cumsum(f2(z$pcs, S1, S2)),  col=1)

    #test f3s
    dmat = as.matrix(dist(points))^2
    f3_true = (dmat[S1, S2] + dmat[S1, S3] - dmat[S2, S3]) / 2
    f3_data = sum(f3(t(d2), S1, S2, S3)) / N_SNPS

        
    print("-------------------")
    print(sprintf("f3 from true points: %f", f3_true))
    print(sprintf("f3 from data: %f", f3_data))
    print("f3 from PPCA 2 PCs")
    print(csf(f3(x$pcs, S1, S2, S3)))
    print("f3 from PPCA 10 PCs")
    print(csf(f3(y$pcs, S1, S2, S3)))
    print("f3 from PCA first 10 PCs")
    print(csf(f3(z$pcs, S1, S2, S3)))

    plot(NA, xlim=c(0, N_SAMPLES), ylim=c(min(0, f3_data, f3_true), max(f3_data, f3_true)))
    abline(h=c(f3_true, f3_data), col=2:3, lty=2)
    lines(cumsum(f3(x$pcs, S1, S2, S3)), col=4,lwd=2)
    lines(cumsum(f3(y$pcs, S1, S2, S3)), col=5)
    lines(cumsum(f3(z$pcs, S1, S2, S3)), col=1)

    #test f4s
    f4_true = -(dmat[S1, S3] + dmat[S2, S4] - dmat[S1, S4] - dmat[S2, S3]) / 2
    f4_data = sum(f4(t(d2), S1, S2, S3, S4)) / N_SNPS

        
    print("-------------------")
    print(sprintf("f4 from true points: %f", f4_true))
    print(sprintf("f4 from data: %f", f4_data))
    print("f4 from PPCA 2 PCs")
    print(csf(f4(x$pcs, S1, S2, S3)))
    print("f4 from PPCA 10 PCs")
    print(csf(f4(y$pcs, S1, S2, S3)))
    print("f4 from PCA first 10 PCs")
    print(csf(f4(z$pcs, S1, S2, S3)))

    plot(NA, xlim=c(0, N_SAMPLES), ylim=c(min(0, f4_data, f4_true), max(f4_data, f4_true)))
    abline(h=c(f4_true, f4_data), col=2:3, lty=2)
    lines(cumsum(f4(x$pcs, S1, S2, S3, S4)), col=4,lwd=2)
    lines(cumsum(f4(y$pcs, S1, S2, S3, S4)), col=5)
    lines(cumsum(f4(z$pcs, S1, S2, S3, S4)), col=1)

}

ppca2 = function(Y){
    #formula from agrawal
    C = diag(1, ncol=4, nrow=100)
    X = solve(t(C) %*% C) %*% t(C) %*% t(Y)
    C = t(Y) %*% t(X) %*% solve(X %*% t(X))
}
