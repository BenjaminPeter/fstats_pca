#dummy inits
init_sigma <- function() {1}
init_cov = function(n_dims,n_pcs)diag(1, nrow=n_dims, ncol=n_pcs)

#helper function for display
csf <- function(...)cumsum(...)[c(1:10, length(...))]
f2 <- function(data, S1, S2){(data[S1,] - data[S2,])^2 }
f3 <- function(data, S1, S2, S3){(data[S1,] - data[S2,]) * (data[S1,] - data[S3,])}
f4 <- function(data, S1, S2, S3, S4){(data[S1,] - data[S3,]) * (data[S4,] - data[S2,])}

double_center <- function(dist)
    t(dist - rowMeans(dist)) - rowMeans(dist) + mean(dist)


ppca_iter <- function(data, n_pcs=4, n_iter=1000, tol=1e-4, verbose=T){
    n_dims = ncol(data)
    #eq 31
    mu = rowMeans(data)
    centered_data = data - mu
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

    #posterior sanity checks, eq 9
    M = t(W) %*% W + diag(sigma, n_pcs)
    M_inv = solve(M)

    posterior_mean = M_inv %*% t(W) %*% t(centered_data )
    approx_data = W %*% solve(t(W) %*% W) %*% M %*% posterior_mean + mu

    return(list(W=W, sigma=sigma, pcs=pcs, eval=e$values, evec=e$vectors,
                data=approx_data, est_data=posterior_mean + mu))
}

ppca_ml <- function(data, n_pcs=4, n_iter=1000, tol=1e-4, verbose=T){
    n_dims = ncol(data)
    #eq 31
    mu = rowMeans(data)
    centered_data = data - mu

    S = t(centered_data) %*% centered_data / nrow(data)
    SVD = svd(S)

    #eq 8
    sigma = sum(SVD$d[(n_pcs+1):n_dims]) / (n_dims - n_pcs)

    #eq 7
    W = SVD$u[,1:n_pcs] %*% sqrt(diag(SVD$d[1:n_pcs] - sigma))
    e = svd(W)
    pcs = e$u %*% diag(e$d)
    pcs = cbind(pcs, matrix(0, ncol=n_dims-n_pcs, nrow=n_dims))


    return(list(W=W, sigma=sigma, pcs=pcs, eval=e$d, evec=e$u))
}

ppca = ppca_ml

ppca_cor <- function(data, n_pcs=4, normalized=F){
    n_dims = ncol(data)
    n_snps = nrow(data)

    if(normalized){ #such that data is between 0 and 1
        dhat = colMeans( data * (1-data))  
        gram_mat <- t(data) %*% data / n_snps - diag(dhat)
        cm = double_center(gram_mat) 
    } else{ #unnormalized version
        dhat = colMeans( data * (2-data))  
        gram_mat <- t(data) %*% data / n_snps - diag(dhat)
        cm = double_center(gram_mat) 
    }

    SVD = svd(cm)

    #eq 8
    sigma = sum(SVD$d[(n_pcs+1):n_dims]) / (n_dims - n_pcs)

    #eq 7
    W = SVD$u[,1:n_pcs] %*% sqrt(diag(SVD$d[1:n_pcs] - sigma))
    e = svd(W)
    pcs = e$u %*% diag(e$d)
    pcs = cbind(pcs, matrix(0, ncol=n_dims-n_pcs, nrow=n_dims))


    return(list(W=W, sigma=sigma, pcs=pcs, eval=e$d, evec=e$u))
}
pca_cor <- function(data, n_pcs=4, tol=1e-4, verbose=T, normalized=T){
    n_dims = ncol(data)
    n_snps = nrow(data)

    if(normalized){ #such that data is between 0 and 1
        dhat = colMeans( data * (1-data))   
        gram_mat <- t(data) %*% data / n_snps - diag(dhat)
        cm = double_center(gram_mat) 
    } else{ #unnormalized version
        dhat = colMeans( data * (2-data))  
        gram_mat <- t(data) %*% data / n_snps - diag(dhat)
        cm = double_center(gram_mat) 
    }

    e = svd(cm)
    pcs = e$u %*% sqrt(diag(e$d)  )


    return(list(pcs=pcs, eval=e$d))
}

#pca_cor = ppca_cor

test_ppca <- function(){
    require(glue)
    require(MASS)
    N_SNPS = 100000
    N_SAMPLES = 100
    TRUE_RANK = 4
    S1 = 5; S2 = 39
    S3 = 3; S4 = 30

    N_PCS1 = 5
    N_PCS2 = 20
    set.seed(10)


    points = matrix(runif(N_SAMPLES * TRUE_RANK), nrow=N_SAMPLES) / 10
    test_cov = points %*% t(points)
    data = mvrnorm(N_SNPS, mu=rep(0.5, N_SAMPLES), Sigma=test_cov); 
    dceil = pmax(pmin(data, 1), 0)
    #d2=data + mvrnorm(N_SNPS, mu=rep(0, N_SAMPLES), Sigma=diag(NOISE, N_SAMPLES))
    d2 = t(apply(dceil, 1, function(p)rbinom(ncol(dceil), 2, p))) / 2
    het = colMeans(d2 * (1-d2))

    x = ppca(d2, verbose=F, n_pcs=N_PCS1)
    y = ppca(d2, verbose=F, n_pcs=N_PCS2)
    z = prcomp(d2)
    z$pcs = unname(t(t(z$rotation) * z$sdev))
    w = pca_cor(d2 , normalized=T)
    v1 = ppca_cor(d2 , n_pcs=N_PCS1, normalized=T)
    v2 = ppca_cor(d2 , n_pcs=N_PCS2, normalized=T)

    #test f2s
    f2_true = as.matrix(dist(points))[S1,S2]^2
    f2_data = sum((d2[,S1] - d2[,S2])^2) / N_SNPS
    f2_cor = f2_data - het[S1] - het[S2]

        
    par(mfrow=c(1,3))

    plot(NA, 
         xlab="# PCS", ylab="stat",
         xlim=c(0, 30), ylim=c(0, 0.03))
    title(sprintf("f2(%d, %d)", S1, S2))
    abline(h=c(f2_true, f2_cor), col=2:3, lty=2)
    lines(cumsum(f2(x$pcs, S1, S2)),  col=4,lwd=2)
    lines(cumsum(f2(y$pcs, S1, S2)),  col=5)
    lines(cumsum(f2(z$pcs, S1, S2)),  col=1)
    lines(cumsum(f2(w$pcs, S1, S2)),  col='purple', lwd=2)
    lines(cumsum(f2(v1$pcs, S1, S2)),  col='orange', lwd=2)
    lines(cumsum(f2(v2$pcs, S1, S2)),  col='brown', lwd=2)

    #test f3s
    dmat = as.matrix(dist(points))^2
    f3_true = (dmat[S1, S2] + dmat[S1, S3] - dmat[S2, S3]) / 2
    f3_data = sum(f3(t(d2), S1, S2, S3)) / N_SNPS
    f3_cor = f3_data - het[S1] 

        
    plot(NA,
         xlab="# PCS", ylab="stat",
         xlim=c(0, 30),
         #ylim=c(min(0, f3_data, f3_true), max(f3_data, f3_true)))
         ylim=c(0, 0.04))
    title(sprintf("f3(%d; %d, %d)", S1, S2, S3))
    abline(h=c(f3_true, f3_cor), col=2:3, lty=2)
    lines(cumsum(f3(x$pcs, S1, S2, S3)), col=4,lwd=2)
    lines(cumsum(f3(y$pcs, S1, S2, S3)), col=5)
    lines(cumsum(f3(z$pcs, S1, S2, S3)), col=1)
    lines(cumsum(f3(w$pcs, S1, S2, S3)), col='purple', lwd=2)
    lines(cumsum(f3(v1$pcs, S1, S2, S3)),  col='orange', lwd=2)
    lines(cumsum(f3(v2$pcs, S1, S2, S3)),  col='brown', lwd=2)

    #test f4s
    f4_true = sum(f4(points, S1, S2, S3, S4))
    f4_data = sum(f4(t(d2), S1, S2, S3, S4)) / N_SNPS

        
    plot(NA, 
         xlab="# PCS", ylab="stat",
         xlim=c(0, 30), 
         #ylim=c(min(0, f4_data, f4_true), max(f4_data, f4_true)))
         ylim=c(0, 0.015))
    title(sprintf("f4(%d, %d; %d, %d)", S1, S2, S3, S4))
    abline(h=c(f4_true, f4_true), col=2:3, lty=2)
    lines(cumsum(f4(x$pcs, S1, S2, S3, S4)), col=4,lwd=2)
    lines(cumsum(f4(y$pcs, S1, S2, S3, S4)), col=5)
    lines(cumsum(f4(z$pcs, S1, S2, S3, S4)), col=1)
    lines(cumsum(f4(w$pcs, S1, S2, S3, S4)), col='purple', lwd=2)
    lines(cumsum(f4(v1$pcs, S1, S2, S3, S4)),  col='orange', lwd=2)
    lines(cumsum(f4(v2$pcs, S1, S2, S3, S4)),  col='brown', lwd=2)

    legend("bottomleft", col=c(4,5,1, 3,2, 'purple', 'orange', 'brown'), 
           legend=c(glue("PPCA rank {N_PCS1}"), 
                    glue("PPCA rank {N_PCS2}"), 
                    "PCA", "admixtools estimate", "truth",
                    'corrected PCA', 
                    glue('corrected PPCA rank {N_PCS1})'), 
                    glue('corrected PPCA rank {N_PCS2})') ), 
           lty=c(1,1,1,2,2, 1, 1, 1))

}

ppca2 = function(Y){
    #formula from agrawal
    C = diag(1, ncol=4, nrow=100)
    X = solve(t(C) %*% C) %*% t(C) %*% t(Y)
    C = t(Y) %*% t(X) %*% solve(X %*% t(X))
}


