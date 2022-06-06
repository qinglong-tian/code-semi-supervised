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

Compute_Zeta <-
  function(y_centered,
           X_centered,
           beta_lasso_original,
           beta_lasso_for_m,
           X_centered_all_data)
    # Use "Zeta" to be consistent with the legacy code; which is actually $\xi$ for the the safe estimator
  {
    n <- length(y_centered)
    p <- ncol(X_centered)
    
    newy <-
      as.numeric(y_centered - X_centered[1:n, ] %*% beta_lasso_original)
    sdxinv <- 1 / sqrt(colSums(X_centered ^ 2) / (n - 1))
    newY <- diag(newy) %*% X_centered[1:n, ] %*% diag(sdxinv)
    
    Y_all <- X_centered_all_data %*% beta_lasso_for_m
    XX2 <-
      diag(as.numeric(Y_all)) %*% X_centered_all_data %*% diag(sdxinv)
    XXs <- scale(XX2[1:n, ], scale = F)
    XX_new <- scale(XX2, scale = F)[c(1:n), ]
    
    zeta_new <-
      colMeans(diag(as.numeric(y_centered)) %*% X_centered %*% diag(sdxinv))
    zeta_original <- zeta_new
    
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
    return(zeta_new)
  }
