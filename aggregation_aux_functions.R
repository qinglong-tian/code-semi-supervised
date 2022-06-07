Fit_Lasso <- function(y_centered, X_centered)
{
  fit_b <- glmnet(
    X_centered,
    y_centered,
    family = "gaussian",
    lambda =
      cv.glmnet(X_centered, y_centered, family = "gaussian")$lambda.1se
  )
  beta_ll <- coef(fit_b)[-1]
  
  return(beta_ll)
}

Compute_Noise_Lasso <- function(y_centered, X_centered, beta_lasso)
{
  n <- length(y_centered)
  Y_all <- X_centered %*% beta_lasso
  noise <-
    sum((Y_all[1:n] - y_centered) ^ 2) / (n - length(which(beta_ll != 0)))
  return(noise)
}

Compute_Zeta <- function(y_centered,
                         X_centered,
                         beta_lasso_original,
                         beta_lasso_for_m,
                         X_centered_all_data)
  # Use "Zeta" to be consistent with the legacy code; which is actually $\xi$ for the the safe estimator
{
  n <- length(y_centered)
  p <- ncol(X_centered)
  
  newy <-
    as.numeric(y_centered - X_centered %*% beta_lasso_original)
  sdxinv <- 1 / sqrt(colSums(X_centered ^ 2) / (n - 1))
  newY <- diag(newy) %*% X_centered %*% diag(sdxinv)
  
  Y_all <- X_centered_all_data %*% beta_lasso_for_m
  XX2 <-
    diag(as.numeric(Y_all)) %*% X_centered_all_data %*% diag(sdxinv)
  XXs <- scale(XX2[1:n,], scale = F)
  XX_new <- scale(XX2, scale = F)[c(1:n),]
  
  zeta_new <-
    colMeans(diag(as.numeric(y_centered)) %*% X_centered %*% diag(sdxinv))
  
  BB <- matrix(0, nrow = p, ncol = p)
  for (i in 1:p)
  {
    y_new <- newY[, i]
    fit_new <-
      glmnet(XXs,
             y_new,
             lambda = cv.glmnet(XXs, y_new, family = "gaussian", nfold = 5)$lambda.1se)
    BB[, i] <- coef(fit_new)[-1]
    zeta_new[i] <-
      zeta_new[i] - colMeans(XX_new) %*% coef(fit_new)[-1]
  }
  return(list(Zeta_New = zeta_new, XX_new = XX_new))
}

Compute_Zeta_Two_Models <- function(y_centered,
                                    X_centered,
                                    beta_lasso_origial,
                                    beta_lasso_for_m_current,
                                    beta_lasso_for_m_next,
                                    X_centered_all_data)
  # Aggregate two models (in the safe estimators)
{
  n <- length(y_centered)
  p <- ncol(X_centered)
  
  # The response part of the lasso (estimating B)
  sdxinv <- 1 / sqrt(colSums(X_centered ^ 2) / (n - 1))
  
  newy <-
    as.numeric(y_centered - X_centered %*% beta_lasso_origial)
  newY <- diag(newy) %*% X_centered %*% diag(sdxinv)
  
  # The covariate part
  Y_all_cur <- X_centered_all_data %*% beta_lasso_for_m_current
  XX2_cur <-
    diag(c(Y_all_cur)) %*% X_centered_all_data %*% diag(sdxinv)
  
  Y_all_next <- X_centered_all_data %*% beta_lasso_for_m_next
  XX2_next <-
    diag(c(Y_all_next)) %*% X_centered_all_data %*% diag(sdxinv)
  
  XX2 <- cbind(XX2_cur, XX2_next)
  
  XXs <- scale(XX2[1:n, ], scale = F)
  XX_new <- scale(XX2, scale = F)[1:n, ]
  
  zeta_new <-
    colMeans(diag(as.numeric(y_centered)) %*% X_centered %*% diag(sdxinv))
  
  BB <- matrix(0, nrow = 2 * p, ncol = p)
  for (i in 1:p)
  {
    y_new <- newY[, i]
    fit_new <-
      glmnet(XXs,
             y_new,
             lambda = cv.glmnet(XXs, y_new, family = "gaussian", nfold = 5)$lambda.1se)
    BB[, i] <- coef(fit_new)[-1]
    zeta_new[i] <-
      zeta_new[i] - colMeans(XX_new) %*% coef(fit_new)[-1]
  }
  
  return(list(
    Zeta_New = zeta_new,
    XX2 = XX2,
    BB = BB
  ))
}

Compute_Omega_Hat <- function(X_original_mat, num_of_labeled)
{
  n <- num_of_labeled
  N <- nrow(X_original_mat) - n
  
  X_n1 <- scale(X_original_mat[1:n, ], center = T, scale = T)
  X_new <- scale(X, center = T, scale = T)
  sigma_n <- t(X_new) %*% X_new / (n + N)
  M <-
    InverseLinfty(
      sigma_n,
      n + N,
      resol = 1.2,
      maxiter = 50,
      threshold = 1e-2,
      verbose = T
    )
  
  return(M)
}

Compute_Debias_Safe_Estimator <- function(beta_for_debias,
                                          omega_hat,
                                          zeta_new,
                                          x_mat_original,
                                          sdxinv,
                                          n)
  # beta_for_debias could be \theta_S or lasso estimators
{
  X_n1 <- scale(x_mat_original[1:n, ], center = T, scale = T)
  Xs <- scale(x_mat_original[1:n, ], scale = F)
  
  c(beta_for_debias + omega_hat %*% (zeta_new - t(X_n1) %*% Xs %*% beta_for_debias /
                                       n) * sdxinv)
}

Compute_SE_Debias_Safe_Estimator <- function(y_centered,
                                             x_mat_original,
                                             beta_hat_for_var,
                                             omega_hat,
                                             BB_in_zeta,
                                             XX2_in_zeta,
                                             sdxinv)
{
  n <- length(y_centered)
  N <- nrow(x_mat_original) - n
  
  Xs <- scale(x_mat_original[1:n, ], scale = F)
  X_n1 <- scale(x_mat_original[1:n, ], center = T, scale = T)
  XXs <- scale(XX2_in_zeta[1:n, ], scale = F)
  
  AA <- diag(c(y_centered - Xs %*% beta_hat_for_var)) %*% X_n1
  Ar <- (AA - XXs %*% BB_in_zeta)
  var1 <-
    omega_hat %*% (N / (n + N) * cov(Ar)) %*% t(omega_hat) + n / (n + N) *
    omega_hat * cov(AA) %*% t(omega_hat)
  se1 <- sqrt(diag(var1)) / sqrt(n) * sdxinv
  
  return(se1)
}


