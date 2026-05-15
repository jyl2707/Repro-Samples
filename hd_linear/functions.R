## Repro Samples: High-Dimensional Linear Regression
## Functions for Fisher inversion, model CS, and coefficient CIs
## Reference: Wang et al. (2022), FINITE AND LARGE-SAMPLE INFERENCE FOR MODEL
##            AND COEFFICIENTS IN HIGH-DIMENSIONAL LINEAR REGRESSION WITH REPRO SAMPLES

library(glmnet)
library(MASS)       # ginv()
library(tidyverse)
library(pbapply)
library(rlang)      # parse_expr()
library(intervals)
library(Matrix)     # sparse matrix in tuning_solution_range()


# Projection matrix onto column space of X (handles rank-deficient X via ginv)
opm <- function(X) {
  if (nrow(X) > ncol(X))
    return(X %*% solve(t(X) %*% X) %*% t(X))
  else
    return(X %*% ginv(t(X) %*% X) %*% t(X))
}


# Adaptive lasso with the u* column unpenalized (sigma estimation via u* coef)
ada_lasso_glmnet_sigma <- function(X = X, Y = Y, lambda.min.ratio = 0.01, ...) {
  fl1     <- glmnet(X, Y, family = "gaussian", alpha = 1,
                    lambda.min.ratio = lambda.min.ratio, ...)
  p       <- dim(X)[2] - 1
  cv_fl1  <- cv.glmnet(x = X, y = Y, penalty.factor = c(rep(1, p), 0), ...)
  lambda_tem   <- cv_fl1$lambda.min
  first_coef   <- fl1$beta[, which.min(abs(fl1$lambda - lambda_tem))]
  penalty_factor <- abs(first_coef + 1/sqrt(nrow(X)))^(-1)
  penalty_factor <- penalty_factor[1:p]
  fl1 <- glmnet(x = X, y = Y, penalty.factor = c(penalty_factor, 0),
                family = "gaussian", alpha = 1,
                lambda.min.ratio = lambda.min.ratio, ...)
  return(fl1)
}


# Adaptive lasso (for model estimation inside cal_p_tau)
ada_lasso <- function(X = X, Y = Y, intercept = FALSE, cv = TRUE,
                      standardize = FALSE, ...) {
  cv_fl1       <- cv.glmnet(x = X, y = Y, intercept = intercept,
                             standardize = standardize, ...)
  lambda_tem   <- cv_fl1$lambda.min
  first_coef   <- cv_fl1$glmnet.fit$beta[, which.min(abs(cv_fl1$lambda - lambda_tem))]
  penalty_factor <- abs(first_coef + 1/sqrt(nrow(X)))^(-1)
  if (cv == TRUE)
    fl1 <- cv.glmnet(x = X, y = Y, penalty.factor = penalty_factor,
                     family = "gaussian", intercept = intercept,
                     standardize = standardize, ...)
  else
    fl1 <- glmnet(x = X, y = Y, penalty.factor = penalty_factor,
                  family = "gaussian", alpha = 1, intercept = FALSE,
                  standardize = standardize, ...)
  return(fl1)
}


# BIC / EBIC for a single coefficient vector (used inside tuning_solution_range)
cal_tuning_betahat <- function(beta.hat.long, Y, X, u_star, intercept = intercept) {
  gamma <- 1
  n     <- length(Y)
  p     <- dim(X)[2]
  k     <- sum(beta.hat.long[1:(p+1)] != 0)
  SSE   <- sum((Y - cbind(1, X, u_star) %*% beta.hat.long)^2)
  AIC   <- n * log(SSE/n) + 2 * k
  BIC   <- n * log(SSE/n) + log(n) * k
  EBIC  <- BIC + 2 * gamma * log(choose(p, k))
  return(list(AIC = AIC, BIC = BIC, EBIC = EBIC))
}


# Select models from the BIC-to-EBIC range of the solution path
tuning_solution_range <- function(RSM_fl_1, Y, X, u_star, intercept = intercept,
                                  tuning_method) {
  nlambda <- length(RSM_fl_1$lambda)
  p       <- dim(X)[2]
  beta.hat.long.list <- rbind(RSM_fl_1$a0, RSM_fl_1$beta) %>%
    split(rep(1:nlambda, each = p + 2))
  beta.hat.tuning.list <- lapply(beta.hat.long.list, cal_tuning_betahat,
                                  Y, X, u_star, intercept)

  if (tuning_method == "BIC range") {
    BIC_tuning_index  <- which.min(sapply(beta.hat.tuning.list, function(x) x$BIC))
    EBIC_tuning_index <- which.min(sapply(beta.hat.tuning.list, function(x) x$EBIC))
    beta.hat.long.list[EBIC_tuning_index:BIC_tuning_index] %>%
      do.call(cbind, .) %>% Matrix(sparse = TRUE) %>% return
  } else if (tuning_method == "BIC") {
    BIC_tuning_index <- which.min(sapply(beta.hat.tuning.list, function(x) x$BIC))
    beta.hat.long.list[[BIC_tuning_index]] %>% Matrix(sparse = TRUE) %>% return
  } else if (tuning_method == "AIC") {
    AIC_tuning_index <- which.min(sapply(beta.hat.tuning.list, function(x) x$AIC))
    beta.hat.long.list[[AIC_tuning_index]] %>% Matrix(sparse = TRUE) %>% return
  } else
    stop("incorrect tuning criterion")
}


# Estimate model tau from data using adaptive lasso with k-variable constraint
est.tau.hat <- function(Y, X, k, intercept = FALSE, method) {
  if (method == "cv")
    est.tau.hat.cv(Y, X, intercept)
  else if (method == "constraint")
    est.tau.hat.k(Y, X, k, intercept)
  else
    stop("incorrect method: use 'cv' or 'constraint'")
}

est.tau.hat.cv <- function(Y, X, intercept = FALSE) {
  beta.hat <- ada_lasso(X, Y, nfolds = 3, intercept = intercept) %>%
    coef("lambda.1se")
}

est.tau.hat.k <- function(Y, X, k, intercept = FALSE) {
  lasso.fit <- ada_lasso(X, Y, intercept = intercept, cv = FALSE)
  nzero     <- coef(lasso.fit, lasso.fit$lambda) %>%
    apply(2, function(x) sum(x != 0))
  tau.index <- which(nzero <= k) %>% max
  return(coef(lasso.fit, s = lasso.fit$lambda[tau.index]))
}


# Fisher-Dempster p-value for a candidate model tau
cal_p_tau <- function(tau, Y, X, repro_sample_size = 100,
                      intercept = FALSE, method) {
  print(tau)
  k <- sum(tau != 0)
  n <- length(Y)

  if (k == 0) {
    if (intercept == FALSE) {
      X.tau <- 0
      H.tau <- matrix(0, nrow = n, ncol = n)
    } else {
      X.tau <- matrix(rep(1, n), ncol = 1)
      H.tau <- opm(X.tau)
    }
  } else {
    if (intercept == FALSE)
      X.tau <- X[, tau, drop = FALSE]
    else
      X.tau <- cbind(1, X[, tau, drop = FALSE])
    H.tau <- opm(X.tau)
  }

  tau.hat.long    <- est.tau.hat(Y, X, k, intercept, method)
  epsilon.matrix  <- matrix(rnorm(repro_sample_size * n), ncol = repro_sample_size)
  R2              <- ((diag(n) - H.tau) %*% Y)^2 %>% sum
  Ystar.matrix    <- apply(epsilon.matrix, MARGIN = 2, function(e)
    H.tau %*% Y +
      sqrt(sum(((diag(n) - H.tau) %*% e)^2)) / sqrt(R2) *
      ((diag(n) - H.tau) %*% e))

  if (method == "cv")
    tau.star.list <- apply(Ystar.matrix, 2, function(y_star)
      est.tau.hat.cv(y_star, X, intercept))
  else if (method == "constraint")
    tau.star.list <- apply(Ystar.matrix, 2, function(y_star)
      est.tau.hat.k(y_star, X, k, intercept))
  else
    stop("incorrect method")

  tau.star.list %>% lapply(function(x) which(x[-1] != 0)) -> tau.star.var.list
  tau.star.var.list %>% lapply(function(x) paste(x, collapse = " ")) -> tau.star.char.list
  tau.star.table <- tau.star.char.list %>% unlist %>% table %>% as_tibble

  tau.hat.char <- paste(which(tau.hat.long[-1] != 0), collapse = " ")
  tau.star.table %>% filter_at(1, all_vars(. == tau.hat.char)) %>%
    pull(2) -> freq.tau.hat

  if (is_empty(freq.tau.hat))
    p.tau <- 0
  else
    p.tau <- (tau.star.table %>%
                filter_at(2, all_vars(. <= (freq.tau.hat))) %>%
                pull(2) %>% sum) / repro_sample_size

  return(p.tau)
}


# Construct model confidence set: retain candidate models with p-value > 1 - conf.level
cal_model_cs <- function(unique.model.list, Y, X, repro_sample_size,
                          intercept, conf.level, method) {
  p_tau_vector <- pbsapply(unique.model.list, cal_p_tau,
                            Y, X, repro_sample_size, intercept, method)
  return(list(model_cs      = unique.model.list[which(p_tau_vector > (1 - conf.level))],
              model_list    = unique.model.list,
              p_value_vector = p_tau_vector))
}


# Coefficient confidence sets: union of OLS intervals across models in CS
cal_beta_i_cs <- function(unique.model.list, Y, X, i_list,
                           conf_level, intercept) {
  null_model_ind    <- sapply(unique.model.list, function(mm) length(mm) == 0)
  unique.model.list <- unique.model.list[!null_model_ind]

  if (length(unique.model.list) == 0)
    return(c(0, 0) %>% Intervals(closed = c(TRUE, TRUE), type = "R") %>%
             list %>% rep(length(i_list)))

  p <- dim(X)[2]
  if (intercept == FALSE) {
    conf_table_list <- lapply(unique.model.list, function(m)
      lm(Y ~ X[, m] - 1) %>% confint(level = conf_level))
    cs_matrix_list  <- lapply(i_list, function(i)
      lapply(1:length(unique.model.list), function(im) {
        i_location <- which(unique.model.list[[im]] == i)
        if (length(i_location) == 0) c(0, 0)
        else conf_table_list[[im]][i_location, ]
      }) %>% do.call(rbind, .))
  } else {
    conf_table_list <- lapply(unique.model.list, function(m)
      lm(Y ~ X[, m]) %>% confint(level = conf_level))
    cs_matrix_list  <- sapply(i_list, function(i)
      lapply(1:length(unique.model.list), function(im) {
        i_location <- which(c(0, unique.model.list[[im]]) == i)
        if (length(i_location) == 0) c(0, 0)
        else conf_table_list[[im]][i_location, ]
      }) %>% do.call(rbind, .))
  }

  if (any(null_model_ind))
    cs_matrix_list %>% lapply(function(csm)
      Intervals(rbind(csm, c(0, 0)), closed = c(TRUE, TRUE), type = "R") %>%
        interval_union)
  else
    cs_matrix_list %>% lapply(function(csm)
      Intervals(csm, closed = c(TRUE, TRUE), type = "R") %>%
        interval_union)
}


# Main inference function: Fisher inversion + model CS + coefficient CIs
rps_hlm_inf <- function(Y, X, rps_size_candidate, rps_size_modelcs, beta_index,
                         u_star_matrix = NULL, conf_level_modelcs = 0.95,
                         conf_level_beta_cs = 0.95,
                         u_star_seed = 100 * runif(1),
                         tuning_method = "BIC range",
                         keep_solution_path = FALSE, ...) {
  n         <- length(Y)
  p         <- dim(X)[2]
  intercept <- as.list(match.call())$intercept
  if (is.null(intercept)) intercept <- TRUE

  set.seed(u_star_seed)
  if (is.null(u_star_matrix))
    u_star_matrix <- rnorm(n * rps_size_candidate) %>% matrix(ncol = rps_size_candidate)

  print("Calculating Solution Paths")
  rps_lm_solutions_list <- pbapply(u_star_matrix, MARGIN = 2, function(u_star)
    ada_lasso_glmnet_sigma(X = cbind(X, u_star), Y = Y, ...))

  print("Tuning")
  rps_tau_u_list <- pblapply(1:rps_size_candidate, function(i)
    tuning_solution_range(rps_lm_solutions_list[[i]], Y, X,
                          u_star_matrix[, i], intercept = intercept, tuning_method))

  mod_can_coef_matrix <- do.call(cbind, rps_tau_u_list)[2:(p+1), ] %>% t

  if (is.null(nrow(mod_can_coef_matrix))) {
    unique.model.list <- list(which(mod_can_coef_matrix != 0))
  } else if (nrow(mod_can_coef_matrix) == 1) {
    unique.model.list <- list(which(mod_can_coef_matrix != 0))
  } else {
    all.model.table   <- apply(mod_can_coef_matrix, 1,
                                function(x) toString(which(x != 0))) %>%
      table %>% as_tibble
    unique.model.list <- lapply(as.list(all.model.table %>% pull(1)),
                                 function(x) parse_expr(paste0("c(", x, ")")) %>% eval)
  }

  print("Obtaining model CS")
  model_inf <- cal_model_cs(unique.model.list, Y, X,
                              repro_sample_size = rps_size_modelcs,
                              intercept = intercept,
                              conf.level = conf_level_modelcs,
                              method = "constraint")

  screen_cand_model_list <- unique.model.list[which(model_inf$p_value_vector > 0)]

  print("Obtaining beta_i CS")
  beta_i_inf <- cal_beta_i_cs(screen_cand_model_list, Y, X,
                                i_list = beta_index,
                                conf_level = conf_level_beta_cs,
                                intercept = intercept)

  if (keep_solution_path)
    return(list(candidate_model_list  = unique.model.list,
                model_inf             = model_inf,
                beta_i_cs             = beta_i_inf,
                rps_lm_solutions_list = rps_lm_solutions_list,
                rps_tau_u_list        = rps_tau_u_list,
                u_star_matrix         = u_star_matrix))
  else
    return(list(candidate_model_list = unique.model.list,
                model_inf            = model_inf,
                beta_i_cs            = beta_i_inf,
                rps_tau_u_list       = rps_tau_u_list,
                u_star_matrix        = u_star_matrix))
}
