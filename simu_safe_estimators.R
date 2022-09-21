library(parallel)
setwd("legacy/")
source("nec.R")
setwd("../")
source("aggregation_aux_functions.R")

B <- 3

# Factor levels
n <- 200
p <- 500
alpha <- 0.05
N <- 200
ratio <- 0.08
nlam <- 10
df <- 5
set.seed(7)

# True values of $\theta\ast$
beta_t <- rep(0, p)
beta_t[1:6] <- c(1.1, 0, 2.4, 4, 4, 2)

# Data Generation
data_mc_list <- mclapply(1:B, function(x)
{
  data_generate4(r = .3, n, N, p)
},
mc.cores = detectCores())

# Estimation and Inference
results <- mclapply(data_mc_list, function(data1)
{
  X <- data1$X
  Ytrue <- scale(data1$y, scale = F)
  y1 <- scale(data1$y[1:n], center = T, scale = F)
  Xs <- scale(X[1:n, ], scale = F)
  Xs2 <- scale(X, scale = F)
  sdxinv <- 1 / sqrt(colSums(Xs ^ 2) / (n - 1))
  # Fit Lasso Using All Columns
  beta_ll <- Fit_Lasso(y_centered = y1, X_centered = Xs)
  
  # Fit Lasso Using the Cubic of the 3rd Column
  Xs_E3 <- Xs
  Xs_E3[, 3] <- Xs_E3[, 3]^3
  beta_ll_E3 <- Fit_Lasso(y_centered = y1, X_centered = Xs_E3)
  zeta_new_E3 <- Compute_Zeta(
    y_centered = y1,
    X_centered = Xs,
    beta_lasso_original = beta_ll,
    beta_lasso_for_m = beta_ll_E3,
    X_centered_all_data = Xs2
  )
  
  # Fit Lasso Using the Square the 4th Column
  Xs_E4 <- Xs
  Xs_E4[, 1] <- Xs_E4[, 1]^2
  beta_ll_E4 <- Fit_Lasso(y_centered = y1, X_centered = Xs_E4)
  zeta_new_E4 <- Compute_Zeta(
    y_centered = y1,
    X_centered = Xs,
    beta_lasso_original = beta_ll,
    beta_lasso_for_m = beta_ll_E4,
    X_centered_all_data = Xs2
  )
  
  # Combining No 3rd With No 4th
  Compute_Zeta_Two_Models(
    y_centered = y1,
    X_centered = Xs,
    beta_lasso_origial = beta_ll,
    beta_lasso_for_m_current = beta_ll_E3,
    beta_lasso_for_m_next = beta_ll_E4,
    X_centered_all_data = Xs2
  ) -> zeta_new_two_models_list
  
  ## Prepare Omega Hat
  omegaHat <- Compute_Omega_Hat(X, n)
  
  ## Prepare for the Theta_S
  ### No 3rd or 4th
  Estimate_Inference_Safe_Estimator(
    X = X,
    y1 = y1,
    zeta_list = zeta_new_two_models_list,
    ratio = ratio,
    Xs = Xs,
    omegaHat = omegaHat,
    sdxinv = sdxinv,
    n = n,
    beta_ll = beta_ll,
    alpha = alpha,
    beta_t = beta_t
  ) -> est_inf_34
  
  ### No 3rd
  Estimate_Inference_Safe_Estimator(
    X = X,
    y1 = y1,
    zeta_list = zeta_new_E3,
    ratio = ratio,
    Xs = Xs,
    omegaHat = omegaHat,
    sdxinv = sdxinv,
    n = n,
    beta_ll = beta_ll,
    alpha = alpha,
    beta_t = beta_t
  ) -> est_inf_3
  
  ### No 4th
  Estimate_Inference_Safe_Estimator(
    X = X,
    y1 = y1,
    zeta_list = zeta_new_E4,
    ratio = ratio,
    Xs = Xs,
    omegaHat = omegaHat,
    sdxinv = sdxinv,
    n = n,
    beta_ll = beta_ll,
    alpha = alpha,
    beta_t = beta_t
  ) -> est_inf_4
  
  return(list(
    rmv3 = est_inf_3,
    rmv4 = est_inf_4,
    rmv34 = est_inf_34
  ))
},
mc.cores = detectCores())
