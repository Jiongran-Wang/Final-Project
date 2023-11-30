# Final-Project
Logistic-tree normal mixture model with diagonal covariance matrix



## Installation

```{r eval=FALSE}
devtools::install_github("Jiongran-Wang/LTNMD")
```

## Toy Example
```{r eval=FALSE}
library(LTNMD)
set.seed(1111)
K = 3 # number of clusters
M = 6 # number of OTU
N = 300 # number of samples

true_alpha <- c(0.4, 0.8, 0.8, -0.3, -0.5)
true_beta1 <- c(0.5, 0.3, 1, 1, 1)
true_beta2 <- c(-0.5, -0.3, -1, -1, -1)
true_beta3 <- c(0.1, 0.6, -0.1, 0.1, -0.1)
true_beta <- cbind(true_beta1, true_beta2, true_beta3)
true_gamma <- c(1, 0, 1, 1, 1)
mu1 <- true_alpha * (1-true_gamma) + true_beta1 * true_gamma
mu2 <- true_alpha * (1-true_gamma) + true_beta2 * true_gamma
mu3 <- true_alpha * (1-true_gamma) + true_beta3 * true_gamma
mu <- rbind(mu1, mu2, mu3)

sigma1 <- diag(c(0.001, 0.001, 0.003, 0.003, 0.003))
sigma2 <- diag(c(0.001, 0.001, 0.003, 0.003, 0.003))
sigma3 <- diag(c(0.001, 0.001, 0.003, 0.003, 0.003))

pi1 <- 0.3
pi2 <- 0.3
pi3 <- 0.4

true_labels <- c()
psi <- matrix(0, nrow = M - 1, ncol = N)
for (i in 1:N){
  u <- runif(1)
  if (u <= pi1){
    true_labels[i] <- 1
    psi[, i] <- rmvnorm(1, mean = mu1, sigma = sigma1)
  } else if (u > pi3){
    true_labels[i] <- 3
    psi[, i] <- rmvnorm(1, mean = mu3, sigma = sigma3)
  }else{
    true_labels[i] <- 2
    psi[, i] <- rmvnorm(1, mean = mu2, sigma = sigma2)
  }
}


Ya <- matrix(0, nrow = M - 1, ncol = N)
Yal <- matrix(0, nrow = M - 1, ncol = N)
Ya[1, ] <- sample(seq(9000, 15000), size = N, replace = TRUE)

for (j in 1:N){
  Ya[2, j] <- rbinom(1, Ya[1, j], prob = exp(psi[1, j]) / (1 + exp(psi[1, j])))
}
Yal[1, ] <- Ya[2, ]
Ya[3, ] <- Ya[1, ] - Ya[2, ]

for (j in 1:N){
  Ya[4, j] <- rbinom(1, Ya[2, j], prob = exp(psi[2, j]) / (1 + exp(psi[2, j])))
}
Yal[2, ] <- Ya[4, ]

for (j in 1:N){
  Ya[5, j] <- rbinom(1, Ya[3, j], prob = exp(psi[3, j]) / (1 + exp(psi[3, j])))
}

Yal[3, ] <- Ya[5, ]

for (j in 1:N){
  Yal[4, j] <- rbinom(1, Ya[4, j], prob = exp(psi[4, j]) / (1 + exp(psi[4, j])))
  Yal[5, j] <- rbinom(1, Ya[5, j], prob = exp(psi[5, j]) / (1 + exp(psi[5, j])))
}

res <- gibbs_sampling(S = 1000, Ya = Ya, Yal = Yal, K = 20)
predicted_label <- representative_cluster(res, burnin = 500, K = 20)

```
