## Load functions
setwd("legacy/")
source("nec.R")
setwd("../")

## The data setting: copy from the previous code; data generating
# Factor levels
n = 200
p = 500
alpha = 0.05
N = 200
ratio = 0.08
nlam = 10
df = 5
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
  scale(X[c(1:n), ], scale = F) # Centering X using only the observed data
Xs2 <-
  scale(X, scale = F) #  Centering the X (each column) using both observed and unobserved data

# Two settings: two sets of lasso
## Total (1,3,4,5,6)

## Setting 1: (1,4,5,6)
Xs_E3 <- Xs
Xs_E3[, 3] <- 0

fit_b_E3 <- glmnet(
  Xs_E3,
  y1,
  family = "gaussian",
  lambda =
    cv.glmnet(Xs_E3, y1, family = "gaussian")$lambda.1se
)
beta_ll_E3 <- coef(fit_b_E3)[-1] # Lasso coefficients

## Setting 2: (1,3,5,6)
Xs_E4 <- Xs
Xs_E4[, 4] <- 0

fit_b_E4 <-
  glmnet(
    Xs_E4,
    y1,
    family = "gaussian",
    lambda = cv.glmnet(Xs_E4, y1, family = "gaussian")$lambda.1se
  )
beta_ll_E4 <- coef(fit_b_E4)[-1] # Lasso coefficients

