#' Title gibbs sampling function
#'
#' @param S The number of iterations (default value is 1000)
#' @param Ya (M * N) matrix containing count data on all internal nodes of binary tree (M) and N is the number of samples (tree)
#' @param Yal (M * N) matrix containing count data on all left child of internal nodes of binary tree (M)
#' @param K Predefined number of clusters. This value is usually pretty large and will automatically decrease to the correct number of clusters (default value is 20)
#' @param nu_0 The prior distribution of diagonal variance term is gamma and nu_0 / 2 is the shape parameter of gamma distribution (default value is 10)
#' @param sigma2_0 (nu_0 * sigma2_0) / 2 is the scale parameter of gamma distribution (default value is 6)
#'
#' @return gibbs_sampling function returns a list containing one matrix (N * S) which contains posterior samples of label
#' @export
#'
#' @examples gibbs_sampling(S = 1000, Ya = matrix(c(10, 20, 25, 8, 12, 20), nrow = 2, byrow = TRUE), Yal = matrix(c(8, 12, 20, 4, 6, 10), nrow = 2, byrow = TRUE), K = 20, nu_0 = 10, sigma2_0 = 6)
gibbs_sampling <- function(S = 1000, Ya, Yal, K = 20, nu_0 = 10, sigma2_0 = 6){
  # d = L - 1, where L is # of OTUs
  d <- nrow(Ya)
  L <- d + 1
  # number of samples
  N <- ncol(Ya)
  # kappa
  kappa <- Yal - Ya/2

  # alpha
  set.seed(1)
  ALPHA <- matrix(0, nrow = d, ncol = S)
  ALPHA[,1] <- current_alpha <- rmvnorm(1, mean = rep(0, d), sigma = diag(rep(5, d)))

  # beta
  BETA <- list()
  current_beta <- t(rmvnorm(K, mean = rep(0, d), sigma = diag(rep(5, d)))) # d x K
  BETA[[1]] <- current_beta

  # p
  P <- rep(0, S)
  a <- 1
  b <- 1
  P[1] <- current_p <- rbeta(1, a, b)

  # gamma
  GAMMA <- matrix(0, nrow = d, ncol = S)
  GAMMA[,1] <- current_gamma <- rcat(d, prob = c(.5, .5)) - 1

  # pi
  zeta <- rep(1 / K, K)
  PI <- matrix(0, nrow = K, ncol = S)
  PI[,1] <- current_pi <- rdirichlet(1, zeta)

  # c
  C <- matrix(0, nrow = N, ncol = S)
  C[,1] <- current_c <- rcat(N, p = current_pi)

  # Omega
  OMEGA <- list()
  current_omega <- list()
  for (k in 1:K){
    current_omega[[k]] <- diag(rep(1, d))
  }
  OMEGA[[1]] <- current_omega

  # Psi
  # this prior setting is based on the meaning of Psi
  PSI <- list()
  current_psi <- matrix(0, nrow = d, ncol = N)
  for (i in 1:d){
    for (j in 1:N){
      if (Ya[i, j] == 0){
        current_psi[i, j] <- 0
      } else{
        tmp <- Yal[i, j] / Ya[i, j]
        if (tmp == 0){
          current_psi[i, j] <- -5
        } else if (tmp == 1){
          current_psi[i, j] <- 5
        } else{
          current_psi[i, j] <- qlogis(tmp)
        }
      }
    }
  }

  PSI[[1]] <- current_psi

  # W: auxiliary variable (Polya Gamma data augmentation)
  W <- list()
  current_w <- matrix(0, nrow = d, ncol = N)
  for (i in 1:d){
    for (j in 1:N){
      if (Ya[i, j] != 0){
        current_w[i, j] <- rpg.sp(1, Ya[i, j] + 0.0001, 0)
        # first check if it is NA
        # then check if it is positive
        while (is.na(current_w[i, j]) || current_w[i, j] < 0){
          current_w[i, j] <- rpg.sp(1, Ya[i, j] + 0.0001, 0)
          warning('NA in PG variable')
        }
      }
    }
  }
  W[[1]] <- current_w

  # Gibbs Sampling
  for (s in 2:S){
    # obtain current values
    current_alpha <- ALPHA[,s-1]
    current_beta <- BETA[[s-1]]
    current_p <- P[s-1]
    current_gamma <- GAMMA[,s-1]
    current_pi <- PI[,s-1]
    current_c <- C[,s-1]
    current_psi <- PSI[[s-1]]
    current_w <- W[[s-1]]
    current_omega <- OMEGA[[s-1]] #list

    # update Psi
    tmp_psi <- current_psi

    for (i in 1:N){
      current_label <- current_c[i]
      current_mu <- current_alpha * (1 - current_gamma) + current_beta[,current_label] * current_gamma

      cov_mat <- chol2inv(chol(current_omega[[current_label]] + diag(current_w[,i])))

      mean_vector <- cov_mat %*% (current_omega[[current_label]] %*% current_mu + kappa[,i])
      tmp_psi[,i] <- rmvnorm(1, mean_vector, cov_mat)
    }
    PSI[[s]] <- tmp_psi

    # update Alpha
    tmp_cov <- matrix(0, nrow = d, ncol = d)
    tmp_mean <- rep(0, d)
    for (i in 1:N){
      current_label <- current_c[i]
      tmp_cov <- tmp_cov + current_omega[[current_label]]
      tmp_mean <- tmp_mean + (tmp_psi[,i] - current_beta[,current_label] * current_gamma) %*% current_omega[[current_label]]
    }

    cov_mat <- chol2inv(chol(tmp_cov * t(t(1 - current_gamma)) %*%  (1 - current_gamma) + diag(rep(1/5, d))))

    mean_vector <- cov_mat %*% t(tmp_mean * (1 - current_gamma))
    tmp_alpha <- rmvnorm(1, mean_vector, cov_mat)
    ALPHA[, s] <- tmp_alpha

    # update Beta
    tmp_beta <- current_beta
    for (k in 1:K){
      tmp_cardinality <- sum(current_c == k)

      cov_mat <- chol2inv(chol((tmp_cardinality * current_omega[[k]]) * t(t(current_gamma)) %*% current_gamma  + diag(rep(1/5, d))))

      if (tmp_cardinality == 0){
        tmp_beta[,k] <- rmvnorm(1, rep(0, d), diag(rep(5, d)))
        next
      }

      if (tmp_cardinality == 1){
        mean_vector <- cov_mat %*% t(((tmp_psi[, which(current_c == k)] - tmp_cardinality * tmp_alpha * (1 - current_gamma)) %*% current_omega[[k]]) * current_gamma)
        tmp_beta[,k] <- rmvnorm(1, mean_vector, cov_mat)
      }else{
        mean_vector <- cov_mat %*% t(((rowSums(tmp_psi[, which(current_c == k)]) - tmp_cardinality * tmp_alpha * (1 - current_gamma)) %*% current_omega[[k]]) * current_gamma)
        tmp_beta[,k] <- rmvnorm(1, mean_vector, cov_mat)
      }
    }
    BETA[[s]] <- tmp_beta

    # update W
    tmp_w <- current_w
    for (i in 1:N){
      for (j in 1:d){
        # we only update w when Ya[j, i] > 0
        if (Ya[j, i] != 0){
          tmp_w[j, i] <- rpg.sp(1, as.numeric(Ya[j, i]), tmp_psi[j, i])
          while (is.na(tmp_w[j, i]) || tmp_w[j, i] < 0){
            tmp_w[j, i] <- rpg.sp(1, as.numeric(Ya[j, i]), tmp_psi[j, i])
            warning('NA in PG variable')
          }
        }
      }
    }
    W[[s]] <- tmp_w

    # update OMEGA
    for (k in 1:K){
      num <- sum(current_c == k)
      if (num == 0){
        for (i in 1:d){
          current_omega[[k]][i,i] <- rgamma(1,
                                            shape = nu_0/2,
                                            (nu_0 * sigma2_0) / 2)
        }
      }else{
        for (i in 1:d){
          current_omega[[k]][i, i] <- rgamma(1,
                                             shape = (nu_0 + num)/2,
                                             rate = (nu_0 * sigma2_0 + sum((tmp_psi[i, current_c == k] - tmp_alpha[i] * (1 - current_gamma[i]) - tmp_beta[i, k] * current_gamma[i])^2))/2)
        }
      }
    }
    OMEGA[[s]] <- current_omega

    # update Pi
    new_zeta <- zeta
    for (k in 1:K){
      new_zeta[k] <- new_zeta[k] + sum(current_c == k)
    }
    tmp_pi <- rdirichlet(1, alpha = new_zeta)
    PI[, s] <- tmp_pi

    # update c
    for (i in 1:N){
      probs <- c()
      for (k in 1:K){
        probs[k] <- log(tmp_pi[k]) + dmvnorm(x = tmp_psi[, i],
                                             mean = tmp_alpha * (1 - current_gamma) + tmp_beta[, k] * current_gamma,
                                             sigma = chol2inv(chol(current_omega[[k]])),
                                             log = TRUE)
      }
      current_c[i] <- rcat(1, p = exp(probs) / sum(exp(probs)))
    }

    C[, s] <- current_c

    # update Gamma
    tmp_gamma <- current_gamma
    for (i in 1:d){
      gamma_0 <- gamma_1 <- tmp_gamma
      gamma_0[i] <- 0
      gamma_1[i] <- 1
      tmp_0 <- log(1 - current_p) + sum(apply(rbind(tmp_psi, current_c),
                                              MARGIN = 2,
                                              FUN = function(x){dmvnorm(x[1:d],
                                                                        mean = tmp_alpha * (1-gamma_0) + tmp_beta[, x[d+1]] * gamma_0,
                                                                        sigma = diag(1 / diag(current_omega[[x[d+1]]])),
                                                                        log = TRUE)}))

      tmp_1 <- log(current_p) + sum(apply(rbind(tmp_psi, current_c),
                                          MARGIN = 2,
                                          FUN = function(x){dmvnorm(x[1:d],
                                                                    mean = tmp_alpha * (1-gamma_1) + tmp_beta[, x[d+1]] * gamma_1,
                                                                    sigma = diag(1 / diag(current_omega[[x[d+1]]])),
                                                                    log = TRUE)}))
      tmp_o <- exp(tmp_0 - tmp_1)
      tmp_gamma[i] <- rbinom(1, 1, prob = 1 / (1 + tmp_o))
    }
    GAMMA[, s] <- tmp_gamma

    # update p
    tmp_p <- rbeta(1, a + sum(tmp_gamma), b + d - sum(tmp_gamma))
    P[s] <- tmp_p

  }
  return(list(C = C))
}
