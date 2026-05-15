## Repro Samples: Gaussian Mixture Model with Unknown K
## Functions for Fisher inversion and nuclear mapping
## Reference: Xie and Wang (2022), ReproGeneralWithNames

library(ClusterR)   # GMM()
library(flexmix)    # stepFlexmix(), FLXMRglmfix(), getModel()
library(dplyr)      # pipe %>%, tibble operations


# Generate repro samples Y* given membership vector M
generate_Y_star <- function(Y, M, repro_size) {

  n <- length(Y)
  u_star_matrix <- matrix(rnorm(n * repro_size), ncol = repro_size)
  Y_star_matrix <- matrix(NA, nrow = n, ncol = repro_size)

  if (is.vector(M) | NCOL(M) == 1)
    lapply(M %>% unique() %>% sort %>% as.list, function(x) which(M == x)) -> cluster_index_list
  else
    M %>% as_tibble %>% lapply(function(x) which(x == 1)) -> cluster_index_list

  for (i in 1:length(cluster_index_list)) {
    cluster_index_i <- cluster_index_list[[i]]
    if (length(cluster_index_i) == 1)
      Y_star_matrix[cluster_index_i, ] <- Y[cluster_index_i]
    else
      Y_star_matrix[cluster_index_i, ] <- sapply(1:repro_size, function(col)
        mean(Y[cluster_index_i]) +
          sd(Y[cluster_index_i]) / sd(u_star_matrix[cluster_index_i, col]) *
          (u_star_matrix[cluster_index_i, col] - mean(u_star_matrix[cluster_index_i, col])))
  }
  return(Y_star_matrix)
}


# Estimate K via BIC over a grid of candidate values
Estimate_K <- function(Y, K_grid, seed_local) {

  n <- length(Y)
  Y <- matrix(Y, ncol = 1)
  log.lik.vec <- BIC.vec <- rep(0, length(K_grid))

  for (kk in K_grid) {
    fit.tmp       <- GMM(Y, kk, seed = seed_local)
    log.lik.tmp   <- sum(log(exp(fit.tmp$Log_likelihood) %*% fit.tmp$weights))
    log.lik.vec[kk] <- log.lik.tmp
    BIC.vec[kk]     <- -2 * log.lik.tmp + (3 * kk - 1) * log(n)
  }
  return(list(BIC = BIC.vec, K_hat = K_grid[which.min(BIC.vec)]))
}


# Select best K from a stepFlexmix run using a modified BIC
getModel_ModifiedBIC <- function(temp_flexmix, n, K_h_grid) {
  LogLik       <- sapply(temp_flexmix@models, function(m) m@logLik)
  RSS_d_n      <- exp(-2 * LogLik / n)
  modified_BIC <- n * log(RSS_d_n + 1/n) + log(n) * (2 * K_h_grid + 1)
  getModel(temp_flexmix, which.min(modified_BIC))
}


# Fisher inversion: build candidate membership matrix M via repro draws
# Each draw u* ~ N(0,I_n) is used as regressor in stepFlexmix(y ~ u*)
Get_M_u_BIC_modified <- function(ys.obs, K_h_grid, repro_M_size) {

  n <- length(ys.obs)
  M_u_matrix <- sapply(1:repro_M_size, function(i) {
    u_star <- rnorm(n)
    stepFlexmix(ys.obs ~ u_star, k = K_h_grid,
                model = FLXMRglmfix(varFix = TRUE)) %>%
      getModel_ModifiedBIC(n, K_h_grid) %>%
      `@`(cluster)
  })
  return(unique(M_u_matrix, MARGIN = 2))
}


# Nuclear mapping: compute nuclear statistic for a candidate membership M
# Returns Nuclear_values_2 = P(K_hat* less extreme than K_hat_obs | M)
Get_nuclear_M <- function(ys.obs, M, K_hat_obs, repro_size, K_h_grid, seed_local_repro) {

  Y_star_matrix  <- generate_Y_star(ys.obs, M, repro_size)
  solve_K_list   <- apply(Y_star_matrix, 2, function(x)
    Estimate_K(x, K_grid = K_h_grid, seed_local_repro))
  K_hat_vector   <- sapply(solve_K_list, function(x) x$K_hat)
  K_hat_vector %>% table %>% sort -> K_hat_rank_table
  K_hat_rank_table %>% names %>% as.numeric() -> K_rank_vector

  if (is.vector(M) | NCOL(M) == 1)
    K_h <- length(unique(M))
  else
    K_h <- NCOL(M)

  Nuclear_values_1 <- mean(abs(K_hat_vector - K_h) < abs(K_hat_obs - K_h))

  if (!K_hat_obs %in% K_rank_vector)
    Nuclear_values_2 <- 1
  else
    Nuclear_values_2 <- 1 - sum(K_hat_rank_table[1:which(K_rank_vector == K_hat_obs)]) /
      sum(K_hat_rank_table)

  return(list(Nuclear_values_1 = Nuclear_values_1,
              Nuclear_values_2 = Nuclear_values_2))
}
