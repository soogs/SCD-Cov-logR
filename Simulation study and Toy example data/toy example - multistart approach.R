# toy example - multistart for scd-cov-logR #
# initiated: 30-april-2021 #
# last modified: 10-may-2021 #

# 0. Introduction ####
# On the toy example dataset, this R script performs model selection and analysis using the multistart approach.
# 20 different random starting values are adopted.
# For each set of starting values, the model selection and analysis are conducted.
# and in the end, the final loss value (objective criterion) is observed to see the minimal value.

# 1. Source functions, load packages ####

setwd("E:\\Users\\park\\Desktop\\Final functions for Github") # TO BE CHANGED

# data generation function (for toy example and simulation study)
source("./3cov_datagen.R") 

Rcpp::sourceCpp("./updateW_diffpen.cpp")

source("./scd_cov_logr.R")

source("./scd_cv.R")


library(nFactors)

# 2. Generate dataset ####

cond <- data.frame(dimension = "low",
                   py_pattern = 2,
                   signal_level = 0.5,
                   reps = 1)


if (cond$dimension == "low"){
  I <- 100
  J <- 30
  
  d1_active <- 1:4
  c_active_block1 <- 5:8
  c_active_block2 <- 16:19
  d2_active <- 20:23
}

if (cond$dimension == "high"){
  I <- 100
  J <- 200
  
  d1_active <- 1:25
  c_active_block1 <- 26:45
  c_active_block2 <- 101:120
  d2_active <- 121:145
  
}

if (cond$py_pattern == 1){
  Py <- c(3, 0.05, -4) # common component is irrelevant
}

if (cond$py_pattern == 2){
  Py <- c(3, -4, 0.05) # d2 is irrelevant
}

cond$signal_level <- as.numeric(cond$signal_level)
cond$py_pattern <- as.numeric(cond$py_pattern)
cond$reps <- as.numeric(cond$reps)

truePx <- matrix(0, nrow = 30, ncol = 3)

truePx[d1_active,1] <- 1
truePx[c_active_block1,2] <- 1
truePx[c_active_block2,2] <- 1
truePx[d2_active,3] <- 1

# for the toy dataset, the following inputs were used
dat <- ck_data(I = I, J = J, d1_active = d1_active, 
               c_active_block1 = c_active_block1, c_active_block2 = c_active_block2,
               d2_active = d2_active, comp_sd = rep(50,3), signal_level = 0.5,
               Py = Py, seed = 11, strict_uncorrelated = F)

X_train <- dat$X_train
y_train <- dat$y_train

# standardize the predictor variables
X_train <- scale(X_train,T,T)

X_test <- dat$X_test
y_test <- dat$y_test

X_test <- scale(X_test,T,T)

# 3. Model selection ####
svdd <- svd(X_train)

nsc <- nScree(svdd$d) # 3 component chosen
nsc2 <- nScree(svdd$d^2)

R <- 3

blockcols <- c(J/2, J/2)

blockcols2 <- cumsum(blockcols)

blockindex <- list()

blockindex[[1]] <- 1:blockcols2[1]

for (i in 2:length(blockcols)){
  blockindex[[i]] <- (blockcols2[i-1]+1):blockcols2[i]
}

svdd <- svd(X_train)


# rational start: provided in the paper #
scd1 <- scd_cov_logr(X = X_train, y = y_train, blockcols = c(J/2, J/2), 
                                 R = 3, alpha = 0.2, lasso = c(45/2,45, 45),
                                 glasso = rep(2,3), ridge_y = 1, inits = "rational",
                                 nrstart = 1, MAXITER = 5000, seed = 1, include_rational = TRUE,
                                 stop_value = 1e-5)

scd1$result$W


# 3.1. Cross validation for alpha and ridge ####
# cross validation for alpha and ridge - multistart #

# i set aside 20 starting seeds. then do the cross validation for each of the 20 seeds. 
# then, i will get 20 sets of model parameters, and 20 different models 
# i can compare all of the 20 to see if there is a big impact played by the starting values

alpha_range <- c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)

ridge_range <-  c(50, 30, 10, 5, 3, 1, 0.5, 0.1)

ranges <- expand.grid(alpha_range, ridge_range)

colnames(ranges) <- c("alpha", "ridge")

cv_multistart <- list()

for (nr in 1:20){
  
  scd_cv_results <- matrix(NA, ncol = 3+2, nrow = nrow(ranges))
  
  
  for (i in 1:nrow(ranges)){
    
    
    scdcv1 <- scd_cv(X = X_train, y = y_train, R = R, blockcols = c(J/2, J/2),
                         alpha = ranges[i,]$alpha, lasso = rep(0,3), glasso = rep(0,3), 
                         ridge_y = ranges[i,]$ridge,
                         MAXITER = 500, nrFolds = 5, seed_cv = i+1, seed_method = nr, 
                         inits = "multistart", include_rational = FALSE,
                         type_measure = "mse")
    
    scd_cv_results[i,] <- c(scdcv1$cve, scdcv1$se, sum(scdcv1$convergence),
                            ranges[i,]$alpha, ranges[i,]$ridge)
    
  }
  
  cv_multistart[[nr]] <- scd_cv_results
  
  print(nr)
}


# function to carry out the 1se rule # 
se_select <- function(x){
  
  cv_min_index <- which.min(x[,1])
  
  serule <- x[cv_min_index,1] + x[cv_min_index,2]
  
  scd_possible <- x[x[,1] < serule,]
  
  scd_possible <- scd_possible[order(scd_possible[,4], decreasing = F),]
  
  if (is.vector(scd_possible)){
    scd_chosen <- scd_possible
  } else {
    scd_chosen <- scd_possible[1,]
  }
  
  return(scd_chosen)
}

alpha_ridge_cv_chosen <- lapply(cv_multistart, se_select)

time2 <- Sys.time()




# 3.3. cross validation for lasso and glasso  ####
lasso_range <- c(100, 45, 30, 15, 10, 7, 5, 1, 0.5)

glasso_range <-  c(10, 5, 2, 1, 0.5, 0.1)

ranges <- expand.grid(lasso_range, glasso_range)

colnames(ranges) <- c("lasso", "glasso")


cv_multistart_lasso_glasso <- list()

for (nr in 1:20){
  
  scd_cv_results <- matrix(NA, ncol = 3+2, nrow = nrow(ranges))
  
  
  for (i in 1:nrow(ranges)){
    
    
    scdcv1 <- scd_cv(X = X_train, y = y_train, R = R, blockcols = c(J/2, J/2),
                         alpha = alpha_ridge_cv_chosen[[nr]][4], 
                         lasso = c(ranges[i,]$lasso/2, ranges[i,]$lasso, ranges[i,]$lasso), 
                         glasso = rep(ranges[i,]$glasso,3), 
                         ridge_y = alpha_ridge_cv_chosen[[nr]][5],
                         MAXITER = 500, nrFolds = 5, seed_cv = i+1, seed_method = nr, 
                         inits = "multistart", include_rational = FALSE,
                         type_measure = "mse")
    
    scd_cv_results[i,] <- c(scdcv1$cve, scdcv1$se, sum(scdcv1$convergence),
                            ranges[i,]$lasso, ranges[i,]$glasso)
    
  }
  
  cv_multistart_lasso_glasso[[nr]] <- scd_cv_results
  
  print(nr)
}


# function to carry out the 1se rule: lasso and glasso # 
se_select_lasso_glasso <- function(x){
  
  cv_min_index <- which.min(x[,1])
  
  serule <- x[cv_min_index,1] + x[cv_min_index,2]
  
  scd_possible <- x[x[,1] < serule,]
  
  scd_possible <- scd_possible[order(scd_possible[,4], decreasing = T),]
  
  if (is.vector(scd_possible)){
    scd_chosen <- scd_possible
  } else {
    scd_chosen <- scd_possible[1,]
  }
  
  return(scd_chosen)
}


lasso_glasso_cv_chosen <- lapply(cv_multistart_lasso_glasso, se_select_lasso_glasso)


# actual estimation: with the 20 multiple starts #

scd_multi <- list()

for (nr in 1:20){
  
  alpha_chosen <- alpha_ridge_cv_chosen[[nr]][4]
  
  ridge_chosen <- alpha_ridge_cv_chosen[[nr]][5]
  
  lasso_chosen <- lasso_glasso_cv_chosen[[nr]][4]
  
  glasso_chosen <- lasso_glasso_cv_chosen[[nr]][5]
  
  
  scd1 <- scd_cov_logr(X = X_train, y = y_train, blockcols = c(J/2, J/2), 
                                   R = R, alpha = alpha_chosen, 
                                   lasso = c(lasso_chosen/2, lasso_chosen, lasso_chosen),
                                   glasso = rep(glasso_chosen,3), 
                                   ridge_y = ridge_chosen, 
                                   inits = "multistart", include_rational = FALSE,
                                   nrstart = 1, MAXITER = 5000, seed = nr,
                                   stop_value = 1e-5)
  
  scd_multi[[nr]] <- scd1
  
  print(nr)
}

losses <- unlist(lapply(scd_multi, function(x){(x$loss)}))

which.min(losses)

lapply(scd_multi, function(x){(x$loss)})

scd_multi[[11]]$result$W # very similar weights matrix from the rational start, reported in the paper

scd_multi[[11]]$result$Py

scd_multi[[11]]$result$py0


predicted <- X_test %*% scd_multi[[11]]$result$W %*% scd_multi[[11]]$result$Py + scd_multi[[11]]$result$py0

fitted <- X_train %*% scd_multi[[11]]$result$W %*% scd_multi[[11]]$result$Py + scd_multi[[11]]$result$py0

sum(((1 / (1 + exp(-predicted))) > 0.5) == y_test)

sum(((1 / (1 + exp(-fitted))) > 0.5) == y_train)
# one or two more observations are classified more from the multistart approach actually

lapply(scd_multi, function(x){(x$result$W)})
# actually, the W matrix is pretty different across different starting values
# (there are a couple of weights matrix where the common covariate is not retrieved well)

lapply(scd_multi, function(x){(x$result$Py)})
# Py is very similar everywhere



# model found with rational start (see the other toy example dataset analysis Rscript) #
scd_rational <- scd_cov_logr(X = X_train, y = y_train, blockcols = c(J/2, J/2), 
                                 R = R, alpha = 0.2, lasso = c(45/2, 45, 45),
                                 glasso = rep(2,3), ridge_y = 1, inits = "rational",
                                 include_rational = TRUE, 
                                 nrstart = 1, MAXITER = 5000, seed = 11,
                                 stop_value = 1e-5)

scd_rational$loss 
# much smaller loss from rational start, although the final solutions lead to the same interpretation
# so the results from the rational starting values are reported.

scd_rational$result$W
scd_rational$result$Py



# save(cv_multistart_alpha_ridge, cv_multistart_lasso_glasso, scd_multi, scd_rational, 
#      file = "../toy example/toy_signal050_seed11_multistart.Rdata")

