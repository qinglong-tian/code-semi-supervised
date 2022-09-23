library(parallel)
setwd("legacy/")
source("nec.R")
setwd("../")
source("aggregation_aux_functions.R")

B <- 5

# Factor levels
n <- 200
p <- 500
alpha <- 0.05
N <- 200
ratio <- 0.08
nlam <- 10
df <- 5
set.seed(7)

t0 <- Sys.time()

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
  Xs <- scale(X[1:n,], scale = F)
  Xs2 <- scale(X, scale = F)
  sdxinv <- 1 / sqrt(colSums(Xs ^ 2) / (n - 1))
  
  # Fit Lasso Using All Columns
  beta_ll <- Fit_Lasso(y_centered = y1, X_centered = Xs)
  
  # Debiased Lasso 1
  object <- lasso.proj(X[c(1:n), ], y1, robust = T)
  beta_ll_debias1 <- object$bhat
  beta_ll_debias1_se <- object$se
  beta_ll_debias1_cp <-
    (beta_ll_debias1 + 1.96 * beta_ll_debias1_se > beta_t) &
    (beta_ll_debias1 - 1.96 * beta_ll_debias1_se < beta_t)
  debias_lasso <- list(betaHat = beta_ll_debias1,
                       se = beta_ll_debias1_se,
                       cp = beta_ll_debias1_cp)
  
  # Additive Model
  Y_all <- c()
  bsMat <- c()
  for (i in 1:p) {
    bsmat <- bs(X[, i], df = df, degree = 2)
    bsMat <- cbind(bsMat, bsmat)
  }
  bsMat1 <- scale(bsMat, scale = F)
  index <- c()
  for (i in 1:p) {
    index <- c(index, rep(i, df))
  }
  XX <- bsMat1[1:n, ] # Regressor in group lasso
  lambda.grp <-
    lambdamax(
      XX,
      y = y1,
      index = index,
      penscale = sqrt,
      center = F,
      model = LinReg()
    ) * 0.8 ^ (0:nlam)
  fit1 <-
    grplasso(
      XX,
      y = y1,
      index = index,
      lambda = lambda.grp,
      model = LinReg(),
      penscale = sqrt,
      center = F,
      control = grpl.control(update.hess = "lambda", trace = 0)
    )
  clambda <- gBIC(fit1, y1)
  beta_cc <- fit1$coefficients[, clambda]
  Y_all_additive <- bsMat1 %*% beta_cc
  zeta_new_add <- Compute_Zeta(
    y_centered = y1,
    X_centered = Xs,
    beta_lasso_original = beta_ll,
    predY_m = Y_all_additive,
    X_centered_all_data = Xs2
  )
  
  # Fit Lasso Using the Cubic of the 3rd Column
  X_full_m1 <- X
  X_full_m1[, 3] <- X_full_m1[, 3] ^ 3
  X_full_m1[, 1] <- X_full_m1[, 1] ^ 2
  # X_full_m1[, -c(1, 3)] <- 0
  X_full_m1 <- scale(X_full_m1, scale = F)
  Xm1 <- X_full_m1[1:n, ]
  beta_ll_m1 <- Fit_Lasso(y_centered = y1, X_centered = Xm1)
  yAll_m1 <- X_full_m1 %*% beta_ll_m1
  
  zeta_new_m1 <- Compute_Zeta(
    y_centered = y1,
    X_centered = Xs,
    beta_lasso_original = beta_ll,
    predY_m = yAll_m1,
    X_centered_all_data = Xs2
  )
  
  # Combining m1 and m2
  Compute_Zeta_Two_Models(
    y_centered = y1,
    X_centered = Xs,
    beta_lasso_origial = beta_ll,
    predY_m1 = Y_all_additive,
    predY_m2 = yAll_m1,
    X_centered_all_data = Xs2
  ) -> zeta_new_two_models_list
  
  ## Prepare Omega Hat
  omegaHat <- Compute_Omega_Hat(X, n)
  
  ## Prepare for the Theta_S
  ### m1+m2
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
  ) -> est_inf_comb
  
  ### m1
  Estimate_Inference_Safe_Estimator(
    X = X,
    y1 = y1,
    zeta_list = zeta_new_m1,
    ratio = ratio,
    Xs = Xs,
    omegaHat = omegaHat,
    sdxinv = sdxinv,
    n = n,
    beta_ll = beta_ll,
    alpha = alpha,
    beta_t = beta_t
  ) -> est_inf_m1

  ### additive
  Estimate_Inference_Safe_Estimator(
    X = X,
    y1 = y1,
    zeta_list = zeta_new_add,
    ratio = ratio,
    Xs = Xs,
    omegaHat = omegaHat,
    sdxinv = sdxinv,
    n = n,
    beta_ll = beta_ll,
    alpha = alpha,
    beta_t = beta_t
  ) -> est_inf_add
  
  return(list(
    rmv3 = est_inf_m1,
    rmv34 = est_inf_comb,
    add = est_inf_add,
    debias_lasso = debias_lasso
  ))
},
mc.cores = detectCores())
saveRDS(results, file = "run120timesN200.RDS")
Sys.time()-t0