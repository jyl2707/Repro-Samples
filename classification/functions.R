library(glmnet)
library(foreach)
library(doParallel)
library(doRNG)
library(hdi)
library(intervals)
library(gcdnet)
library(e1071)


get_tau_s <- function(X, Y, epsilons, s, method){
  ########### model candidate set with sparsity upper bound s ###############
  # epsilons: all the randomly generated logistic noises
  # method: "svm adapt" or "logis adapt", they correspond to hinge loss or logistic loss with adpative lasso penalty
  # with_intercept: whether include intercept in logistic regression
  
  n <- nrow(X)
  p <- ncol(X)
  tau_set <- NULL
  tau_candidate <- NULL
  
  for(i in 1:ncol(epsilons)){
    epsilon <- epsilons[, i]
    
    X_work <- cbind(X, epsilon)
    Data <- data.frame(Y = Y, X = X_work)
    
    if(method == "svm adapt"){
      fit_init <- svm(Y~., data = Data, kernel = "linear", type = "C-classification", scale = FALSE)
      coef_fit_init <- as.vector(t(fit_init$coefs) %*% fit_init$SV)
      
      # initial weights for adaptive lasso
      w <- 1 / abs(coef_fit_init[1:p])
      fit <- gcdnet(X_work, Y, method = "hhsvm", pf = c(w, 0), lambda2 = 0)
      betas <- as.matrix(coef(fit))
      df <- colSums(betas[2:(p+1), ] != 0)
      df_index <- which(df <= s)
      hinge_loss <- 1 - diag(2*Y-1) %*% (X_work %*% betas[-1, df_index] + matrix(rep(1, n), nrow = n) %*% matrix(betas[1, df_index], nrow = 1))
      hinge_loss <- colSums(hinge_loss + abs(hinge_loss)) / 2
      
      # EBIC_gamma has a parameter gamma, we include all models between those selected by EBIC_0 and EBIC_1
      bic <- hinge_loss + df[df_index] * log(n)
      ebic <- hinge_loss + df[df_index] * log(n) + 2 * log(choose(p, df[df_index]))
      ind_l <- min(which.min(bic), which.min(ebic))
      ind_u <- max(which.min(bic), which.min(ebic))
      index <- ind_l:ind_u
    } else if(method == "logis adapt"){
      fit_init <- cv.glmnet(X_work, Y, family = "binomial", intercept = with_intercept, nfolds = 3,
                            alpha = 0)
      
      # initial weights for adaptive lasso
      w <- 1 / abs(coef.glmnet(fit_init, s = "lambda.1se")[2:(p+1)])
      fit <- glmnet(X_work, Y, family = "binomial", intercept = with_intercept, penalty.factor = c(w,0))
      
      betas <- as.matrix(coef.glmnet(fit))
      df <- colSums(betas[2:(p+1), ] != 0)
      df_index <- which(df <= s)
      
      # EBIC_gamma has a parameter gamma, we include all models between those selected by EBIC_0 and EBIC_1
      ebic <- deviance(fit)[df_index] + log(n) * df[df_index] + 2 * log(choose(p, df[df_index]))
      bic <- deviance(fit)[df_index] + log(n) * df[df_index]
      ind_l <- min(which.min(bic), which.min(ebic))
      ind_u <- max(which.min(bic), which.min(ebic))
      index <- ind_l:ind_u
      
    }
    
    
    tau_fit <- (betas[2:(p+1), df_index[index]] != 0)
    tau_set <- cbind(tau_set, tau_fit)
    
  }
  
  
  # exclude repetitive models
  tau_set <- t(unique(t(tau_set)))
  
  for(j in 1:ncol(tau_set)){
    tau_candidate <- c(tau_candidate, list(which(tau_set[, j] != 0)))
  }
  return(tau_candidate)
}



nuclear_tau_nuisance <- function(intercept_beta, X, epsilons, tau_fit, tau, with_intercept, n_mc){
  ################# model confidence set based on beta MLE #################
  # inercept_beta: concatenated vector of intercept and beta in MLE
  # epsilons: another set of randomly generated logistic noises for Monte Carlo approximation
  # tau_fit: fitted model support based on the observed data
  # tau: the candidate model we want to test
  # with_intercept: whether include intercept in logistic regression
  # n_mc: number of Monte Carlo samples
  
  tau_prob <- NULL
  tau_fit <- toString(tau_fit)
  n <- nrow(X)
  for(i in 1:n_mc){
    epsilon <- epsilons[, i]
    
    # generate fake responses based on the candidate model tau and generated noise
    if(with_intercept){
      beta <- intercept_beta[-1]
      intercept <- intercept_beta[1]
      Y_star <- as.numeric(intercept + matrix(X[, tau], nrow = n) %*% matrix(beta, ncol = 1) + epsilon >= 0)
    } else{
      beta <- intercept_beta
      Y_star <- as.numeric(matrix(X[, tau], nrow = n) %*% matrix(beta, ncol = 1) + epsilon >= 0)
    }
    
    # fit lasso on the generated fake data to get estimated model support
    fit_star <- glmnet(X, Y_star, family = "binomial", intercept = with_intercept)
    
    tau_index <- max(which(fit_star$df <= length(tau)))
    
    if(fit_star$df[tau_index] == 0){
      tau_star <- 0
    } else{
      tau_star <- which(coef.glmnet(fit_star)[-1, tau_index] != 0)
    }
    
    # empirical distribution of estimated model support
    tau_star <- toString(tau_star)
    if(tau_star %in% names(tau_prob)){
      tau_prob[tau_star] <- tau_prob[tau_star] + 1
    } else{
      tau_prob[tau_star] <- 1
    }
    
    
  }
  if(!(tau_fit %in% names(tau_prob))){
    tau_prob[tau_fit] <- 0
  }
  tau_prob <- tau_prob / n_mc
  nuclear_stat <- sum(tau_prob[tau_prob > tau_prob[tau_fit]])
  return(nuclear_stat)
}




nuclear_tau_profile <- function(X, epsilons, tau_fit, tau, beta_mle, intercept_mle = 0, with_intercept, n_mc){
  ################# model confidence set based on profile method #################
  if(with_intercept){
    opt <- optim(par = c(intercept_mle, beta_mle), fn = nuclear_tau_nuisance, X = X, epsilons = epsilons, tau_fit = tau_fit, 
                 tau = tau, with_intercept = with_intercept, n_mc = n_mc)
  } else{
    opt <- optim(par = beta_mle, fn = nuclear_tau_nuisance, X = X, epsilons = epsilons, tau_fit = tau_fit, 
                 tau = tau, with_intercept = with_intercept, n_mc = n_mc)
  }
  nuclear_stat <- opt$value
  return(nuclear_stat)
}



betaj_cs_wald <- function(X, Y, beta0, tau_set, alpha){
  ############ confidence set for betaj based on Wald test ##############
  # beta0: true beta vector
  # tau_set: model candidate set
  # alpha: confidence level
  
  p <- ncol(X)
  n <- nrow(X)
  betaj_ci_set <- vector(mode = "list", length = p)
  betaj_cs_results <- vector(mode = "list", length = p)
  
  for(tau in tau_set){
    X_work <- X[, tau]
    Data <- data.frame(Y = Y, X_work = X_work)
    beta_mle <- glm(Y ~.-1, data = Data, family = binomial(link = "logit"))$coefficients
    
    gradient <- X_work %*% beta_mle
    gradient[Y==0] <- -gradient[Y==0]
    gradient <- 1 / (1+exp(gradient))
    gradient[Y==0] <- -gradient[Y==0] 
    gradient <- diag(as.vector(gradient)) %*% X_work
    hessian <- exp(X_work%*%beta_mle) / (1+exp(X_work%*%beta_mle))^2
    hessian <- -t(X_work) %*% diag(as.vector(hessian)) %*% X_work / n
    prec_matrix <- solve(hessian)
    variance <- prec_matrix %*% cov(gradient) %*% prec_matrix
    
    for(j in 1:p){
      if(all(j != tau)){
        betaj_ci_set[[j]] <- rbind(betaj_ci_set[[j]], c(0, 0))
      } else{
        betaj_mle <- beta_mle[which(tau == j)]
        variancej <- variance[which(tau == j), which(tau == j)]
        betaj_ci <- betaj_mle + sqrt(variancej/n) * c(qnorm(1/2-alpha/2), qnorm(1/2+alpha/2))
        betaj_ci_set[[j]] <- rbind(betaj_ci_set[[j]], betaj_ci)
      }
    }
  }
  for(j in 1:p){
    betaj_ci_set[[j]] <- unique(betaj_ci_set[[j]])
    betaj_cs_results[[j]] <- sum_beta_interval_matrix(betaj_ci_set[[j]], beta0[j])
  }
  return(betaj_cs_results)
}


sum_beta_interval_matrix<-function(conf_interval_matrix, beta0_i){
  ############# union of intervals and coverage calculation ##############
  # conf_interval_matrix: matrix of confidence intervals for beta_i from different models, two columns: lower and upper bounds
  # beta0_i: true value of beta_i
  x <- Intervals(
    conf_interval_matrix,
    closed = c( TRUE, TRUE ),
    type = "R")
  
  
  conf_interval<-interval_union(x)
  
  coverage<-(distance_to_nearest(beta0_i, conf_interval)==0)
  width_of_interval<-sum(size(conf_interval))
  
  return(list(conf_interval=conf_interval,
              coverage=coverage,
              width_of_interval=width_of_interval))
  
}

