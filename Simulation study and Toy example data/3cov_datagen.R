# 3-covariate data generation #
# initiated: 14 april 2020 #

# This function generates the dataset from a 3-covariate structure
# This data generation function is used for simulation study and the toy example

# Covariate structures:
# 1 common comvariate, 
# 1 covariate distinctive to block 1, 
# 1 covariate distinctive to block 2

ck_data <- function(I, J, d1_active, c_active_block1, c_active_block2, 
                    d2_active, comp_sd, signal_level,
                    Py, seed, strict_uncorrelated){
  
  # Inputs #
  # I: number of observations
  # J: number of predictor variables
  # d1_active: indices of non-zero loadings for the covariate distinctive to the first block
  # c_active_block1, c_active_block2: indices of non-zero loadings for the common covariate
  # d2_active: indices of non-zero loadings for the covariate distinctive to the second block
  # comp_sd: magnitude of standard deviation defined for the covariates
  # signal_level: the proportion of variance accounted for by the covariates
  # Py: regression coefficients (will be normalized)
  # seed: random seed value
  # strict_uncorrelated: if TRUE, the covariates are orthonormalized. This can only work for low-dimensionality.
  # thus, it is always FALSE in the paper. The covariates are anyhow generated from multivariate normal distribution,
  # with a diagonal covariance matrix.
  
  library(MASS)
  
  I <- I
  J1 <- J/2
  J2 <- J/2
  
  set.seed(seed)
  
  # generating the covariates and noise variables, such that they are all uncorrelated
  H_noise <- MASS::mvrnorm(n = I, mu = rep(0, 3+J+1), Sigma = diag(c(comp_sd, rep(comp_sd[1], J), comp_sd[1])^2))
  
  H <- H_noise[,1:3]
  
  # if TRUE, the covariates are orthonormalized (usually FALSE)
  if (strict_uncorrelated){
    H <- scale(H,T,F)
    
    H <- qr.Q(qr(H))
    
    H <- H * sqrt(comp_sd[1]^2 * I)
  }
  
  X_noise <- H_noise[,4:(J+3)]
  
  # defining the loadings matrix Px
  truePx <- matrix(0, nrow = J, ncol = 3)
  
  truePx[d1_active,1] <- 1
  truePx[c_active_block1,2] <- 1
  truePx[c_active_block2,2] <- 1
  truePx[d2_active,3] <- 1
  
  # loadings orthonormalized
  truePx <- qr.Q(qr(truePx))
  
  # Xtrue defined
  Xtrue <- H %*% t(truePx)
  
  # amount of SSE (sum of squared errors) needed
  SS_noise <- (1 - signal_level) / signal_level * sum(Xtrue^2)
  
  # normalize the noise vectors to unit length
  X_noise_normalized <- X_noise %*% diag(1/sqrt(diag(t(X_noise) %*% X_noise)))
  
  # correctly scaling the X_noise to achieve the desired signal_level
  X_noise_scaled <- X_noise_normalized * sqrt(SS_noise/J)
  
  sum(X_noise_scaled^2)
  
  # sum(Xtrue^2) / (sum(Xtrue^2) + sum(X_noise_scaled^2))
  
  # observed X defined
  X_train <- Xtrue + X_noise_scaled
  
  # signal level check
  signal_level_train <- sum(Xtrue^2) / sum(X_train^2)
  
  # regression coefficients scaled
  Py <- Py / sqrt(as.numeric(t(Py) %*% Py))
  
  ymodel <- H %*% Py
  
  probs <- 1 / (1 + exp(-ymodel))
  
  # sampling from binomial distribution
  y_train <- rbinom(n = I, size = 1, prob = probs)
  
  
  # test set: exactly the same procedure as the train set #
  H_noise <- MASS::mvrnorm(n = I, mu = rep(0, 3+J+1), Sigma = diag(c(comp_sd, rep(comp_sd[1], J), comp_sd[1])^2))
  
  H <- H_noise[,1:3]
  
  if (strict_uncorrelated){
    H <- scale(H,T,F)
    
    H <- qr.Q(qr(H))
    
    H <- H * sqrt(comp_sd[1]^2 * I)
  }
  
  X_noise <- H_noise[,4:(J+3)]
  
  truePx <- matrix(0, nrow = J, ncol = 3)
  
  truePx[d1_active,1] <- 1
  truePx[c_active_block1,2] <- 1
  truePx[c_active_block2,2] <- 1
  truePx[d2_active,3] <- 1
  
  truePx <- qr.Q(qr(truePx))
  
  Xtrue <- H %*% t(truePx)
  
  # amount of SSE needed
  SS_noise <- (1 - signal_level) / signal_level * sum(Xtrue^2)
  
  # normalize the noise vectors to unit length
  X_noise_normalized <- X_noise %*% diag(1/sqrt(diag(t(X_noise) %*% X_noise)))
  
  # correctly scaling the X_noise to achieve the desired signal_level
  X_noise_scaled <- X_noise_normalized * sqrt(SS_noise/J)
  
  sum(X_noise_scaled^2)
  
  # sum(Xtrue^2) / (sum(Xtrue^2) + sum(X_noise_scaled^2))
  
  X_test <- Xtrue + X_noise_scaled
  
  signal_level_test <- sum(Xtrue^2) / sum(X_test^2)
  
  Py <- Py / sqrt(as.numeric(t(Py) %*% Py))
  
  ymodel_test <- H %*% Py
  
  probs_test <- 1 / (1 + exp(-ymodel_test))
  
  y_test <- rbinom(n = I, size = 1, prob = probs_test)
  
  
  result <- list(X_train = X_train, y_train = y_train, X_test = X_test, y_test = y_test, truePx = truePx, 
                 Py = Py, signal_level_test = signal_level_test, signal_level_train = signal_level_train)
  
  return (result)
  
}

# 
# dat <- ck_data(I = I, J = J, d1_active = d1_active, 
#                               c_active_block1 = c_active_block1, c_active_block2 = c_active_block2,
#                               d2_active = d2_active, comp_sd = rep(50,3), signal_level = 0.9,
#                               Py = Py, seed = 1313, strict_uncorrelated = F)

