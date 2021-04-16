require(tidyverse)
require(mvtnorm)

set.seed(3)

n_snps = 100000
outgroup_branch_length <- .1

#covariance matrix of samples
n_samples <- 100
posX <- sort(runif(n_samples))
posY <- sort(runif(n_samples))
local_cov = sqrt(2)-as.matrix(dist(cbind(posX,posY)))

#subpop
n_subpop <- 100
posX2 <- c(sort(runif(n_subpop, 0, .1)), posX)
posY2 <- c(sort(runif(n_subpop, 0, .1)), posY)
subpop_cov = sqrt(2)-as.matrix(dist(cbind(posX2,posY2)))

X = rmvnorm(n=n_snps, sigma=local_cov) %>% t
Y = t(t(X) - colMeans(X))

X2 = rmvnorm(n=n_snps, sigma=subpop_cov) %>% t
Y2 = t(t(X2) - colMeans(X2))

Y3 = t(t(X2) - colMeans(X))

pcs1 <- svd(Y)$u[,1:2]
pcs2 <- svd(Y2)$u[,1:2]
pcs3 <- svd(Y3)$u[,1:2]
pcs3b <- svd(Y3)$u[,2:3]
