## Load functions
setwd("legacy/")
source("nec.R")
setwd("../")
source("aggregation_aux_functions.R")

t0 <- Sys.time()

## The data setting: copy from the previous code; data generating
# Factor levels
n <- 200
p <- 500
alpha <- 0.05
N <- 200
ratio <- 0.08
nlam <- 10
df <- 5
an <- 30000
set.seed(7)

# True values of $\theta\ast$
beta_t <- rep(0, p)
beta_t[1:6] <- c(1.1, 0, 2.4, 4, 4, 2)

# Generate data and run lasso

data1 <- data_generate4(r = 0.3, n, N, p) # 400 rows, 500 columns
X <- data1$X # Original X all; use X from all 400 rows

Ytrue <-
  scale(data1$y, scale = F) #Centered Y using both observed and unobserved data

y1 <-
  scale(data1$y[c(1:n)], center = T, scale = F) # Centering Y using only the observed data
Xs <-
  scale(X[c(1:n),], scale = F) # Centering X using only the observed data
Xs2 <-
  scale(X, scale = F) #  Centering the X (each column) using both observed and unobserved data
sdxinv <- 1/sqrt(colSums(Xs^2)/(n - 1))

# Three settings: two sets of lasso
## Total (1,3,4,5,6): this one is used as the "\theta_L" in (3.15)
beta_ll <- Fit_Lasso(y_centered = y1, X_centered = Xs)

## The next two settings are used to compute "\hat{m}(\cdot)" in (3.15)
## Setting 1: (1,4,5,6)
Xs_E3 <- Xs
Xs_E3[, 3] <- 0
beta_ll_E3 <- Fit_Lasso(y_centered = y1, X_centered = Xs_E3)
# Compute Zeta for the Lasso estimator with zero column 3
zeta_new_E3 <-
  Compute_Zeta(
    y_centered = y1,
    X_centered = Xs,
    beta_lasso_original = beta_ll,
    beta_lasso_for_m = beta_ll_E3,
    X_centered_all_data = Xs2
  )$Zeta_New

## Setting 2: (1,3,5,6)
Xs_E4 <- Xs
Xs_E4[, 4] <- 0
beta_ll_E4 <- Fit_Lasso(y_centered = y1, X_centered = Xs_E4)
# Compute Zeta for the Lasso estimator with zero column 3
zeta_new_E4 <-
  Compute_Zeta(
    y_centered = y1,
    X_centered = Xs,
    beta_lasso_original = beta_ll,
    beta_lasso_for_m = beta_ll_E4,
    X_centered_all_data = Xs2
  )$Zeta_New

## Setting 3: Combining settings 1 and 2
Compute_Zeta_Two_Models(
  y_centered = y1,
  X_centered = Xs,
  beta_lasso_origial = beta_ll,
  beta_lasso_for_m_current = beta_ll_E3,
  beta_lasso_for_m_next = beta_ll_E4,
  X_centered_all_data = Xs2
) -> zeta_new_two_models_list

### Estimate Omega
omegaHat <- Compute_Omega_Hat(X, n)

### The paper says that we should use (3.17) for the inference
### But the \theta_L in (3.17) can be other estimator
### Here we can use both the standard lasso estimator and the \theta_S

# Use Standard lasso
beta_for_debias <- beta_ll

# Use \theta_S in (4.1) of the manuscript
dantzig_n2 <-
  dantz(
    X[c(1:n),],
    y1,
    zeta_new_two_models_list$Zeta_New,
    lambda.min.ratio = ratio,
    nlambda = 50,
    clam = 1,
    verbose = T
  ) # Equation 4.1
beta_n2 <- dantzig_n2$beta
chooselam <-
  cv_nonadd(
    nfolds = 5,
    XX2 = zeta_new_two_models_list$XX2,
    Xs = Xs,
    y1 = y1,
    zetafull = zeta_new_two_models_list$Zeta_New,
    BB = zeta_new_two_models_list$BB,
    nlambda = 50,
    ratio = ratio
  )
clam <- chooselam$lambda.min
beta_for_debias <- beta_n2[, clam]

beta_debiased <- Compute_Debias_Safe_Estimator(
  beta_for_debias = beta_for_debias,
  omega_hat = omegaHat,
  zeta_new = zeta_new_two_models_list$Zeta_New,
  x_mat_original = X,
  sdxinv = sdxinv,
  n = n
)

se_beta_debiased <- Compute_SE_Debias_Safe_Estimator(
  y_centered = y1,
  x_mat_original = X,
  beta_hat_for_var = beta_ll,
  omega_hat = omegaHat,
  BB_in_zeta = zeta_new_two_models_list$BB,
  XX2_in_zeta = zeta_new_two_models_list$XX2,
  sdxinv = sdxinv
)

Sys.time()-t0

lwb <- beta_debiased-qnorm(1-alpha/2)*se_beta_debiased
upb <- beta_debiased+qnorm(1-alpha/2)*se_beta_debiased
cp <- mean(beta_t>lwb & beta_t < upb)
