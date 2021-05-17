# simulation 11-may-2021 #
# scd_cov_logr #
# initiaed: 11 may 2021 #

# supplement: same as simulation 11-may-2021, I'm just conducting 20 more reps per condition

# 1. functions and libraries ####
setwd("E:\\Users\\park\\Desktop\\Final functions for Github/") # TO BE CHANGED BY USER

Rcpp::sourceCpp("./updateW_diffpen.cpp") # updateW cpp

source("./scd_cov_logR.R") # scdcovlogR function

source("./scd_cv.R") # crossvalidation for scdcovlogr

source("./3cov_datagen.R") # data generation function

data_b <- ck_data

source("./simulation 11-may-2021/conditions_making_supplement.R")  # SUPPLEMENT CONDITIONS - MAKING USED
# PATH TO BE CHANGED BY USER

library(mixOmics)

permutations <- function(n){
  if(n==1){
    return(matrix(1))
  } else {
    sp <- permutations(n-1)
    p <- nrow(sp)
    A <- matrix(nrow=n*p,ncol=n)
    for(i in 1:n){
      A[(i-1)*p+1:p,] <- cbind(i,sp+(sp>=i))
    }
    return(A)
  }
}

# congruence for 2 vectors
congruence <- function(a, b){
  
  if(anyNA(cbind(a, b))){
    return (0)
  }
  
  if (sum(abs(a)) < 1e-10 | sum(abs(b)) < 1e-10){
    result <- 0
    
    return(result)
  }
  
  result <- (t(a) %*% b) / sqrt(sum(a^2) * sum(b^2))
  return(c(result))
}


# evaluation criteria functions #
corrects <- function(estimate, defined){
  
  if (is.matrix(estimate)){
    total <- prod(dim(estimate))
    
    tucker <- RegularizedSCA::TuckerCoef(defined, estimate)
    
    estimate <- estimate[,tucker$perm]
    
    ratio <- (sum((abs(estimate) > 1e-7) + (abs(defined) > 1e-7) == 2) + 
                sum((abs(estimate) < 1e-7) + (abs(defined) < 1e-7) == 2)) / (total)
    
  } else {
    total <- length(estimate)
    
    ratio <- (sum((abs(estimate) > 1e-7) + (abs(defined) > 1e-7) == 2) + 
                sum((abs(estimate) < 1e-7) + (abs(defined) < 1e-7) == 2)) / (total)
    
  }
  
  return(ratio)
}

nonzero <- function(estimate, defined){
  
  if (is.matrix(estimate)){
    nonzeros <-  sum(abs(defined) > 1e-7)
    
    tucker <- RegularizedSCA::TuckerCoef(defined, estimate)
    
    estimate <- estimate[,tucker$perm]
    
    ratio <- sum((abs(estimate) > 1e-7) + 
                   (abs(defined) > 1e-7) == 2) / (nonzeros)
    
  } else {
    nonzeros <-  sum(abs(defined) > 1e-7)
    
    ratio <- sum((abs(estimate) > 1e-7) + 
                   (abs(defined) > 1e-7) == 2) / (nonzeros)
    
  }
  
  return (ratio)
}


# balanced error rate #
ber <- function(pred, true){
  true_neg <- sum(((pred == 0) + (true == 0)) == 2)
  true_pos <- sum(((pred == 1) + (true == 1)) == 2)
  
  false_neg <- sum(((pred == 0) + (true == 1)) == 2)
  false_pos <- sum(((pred == 1) + (true == 0)) == 2)
  
  fpr <- false_pos / (true_neg + false_pos)
  fnr <- false_neg / (true_pos + false_neg)
  
  result <- (fpr + fnr) /2 
  
  return(result)
}

# common and distinctive identifier #
comdis <- function(estimate){
  J <- nrow(estimate)
  
  d1zero <- colSums(estimate[1:(J/2),]) == 0
  d2zero <- colSums(estimate[(J/2 + 1):J,]) == 0
  
  commz <- (d1zero + d2zero) == 0
  
  result <- list(d1 = sum(d2zero), d2 = sum(d1zero), common = sum(commz))
  
  return(result)
}

# data.frame to save the results #
results <- data.frame(matrix(NA, nrow = nrow(condition_df), ncol = 76))

results[,1] <- factor(results[,1], levels = c("low", "high"))

# there are two different ways of doing SCD here:
# A. model selection according to the paper. First alpha and ridge CV with fixed lasso and glasso values
# -> and then one more run of scdcovlogr (no CV) to find lasso and glasso values that return the correct weights.
# B. First alpha and ridge CV with fixed lasso and glasso values. And the fixed lasso and glasso values are 
# used for the analysis in the end.

# >>>>> in the paper, method A (objects starting with "scd_") is reported.
# method B: objects starting with "scd2_"

colnames(results) <- c("dimensions", "signal_level", "py_pattern", "reps",
                       "model_seed", "noise_seed",
                       
                       "scd_alpha", "scd_ridge",
                       "scd_fit", "scd_pred", "scd_fit_ber", "scd_pred_ber",
                       "scd_tucker",
                       "scd_correct", 
                       "scd_nonzero", 
                       "scd_d1", "scd_d2", "scd_common",
                       "scd_lasso", "scd_glasso", 
                       "scd_zeros1", "scd_zeros2", "scd_zeros3",
                       
                       "scd2_alpha", "scd2_ridge",
                       "scd2_fit", "scd2_pred", "scd2_fit_ber", "scd2_pred_ber",
                       "scd2_tucker",
                       "scd2_correct",
                       "scd2_nonzero",
                       "scd2_d1", "scd2_d2", "scd2_common",
                       "scd2_lasso1", "scd2_lasso", "scd2_glasso",
                       "scd2_zeros1", "scd2_zeros2", "scd2_zeros3",
                       
                       "diacon_fit", "diacon_pred", "diacon_fit_ber", "diacon_pred_ber",
                       "diacon_tucker",
                       "diacon_correct", 
                       "diacon_nonzero", 
                       "diacon_d1", "diacon_d2", "diacon_common",
                       
                       "diathree_fit", "diathree_pred", "diathree_fit_ber", "diathree_pred_ber",
                       "diathree_tucker",
                       "diathree_correct", 
                       "diathree_nonzero", 
                       "diathree_comcheck",
                       
                       "diaweighted_fit", "diaweighted_pred", "diaweighted_fit_ber", "diaweighted_pred_ber",
                       "diaweighted_tucker",
                       "diaweighted_correct", 
                       "diaweighted_nonzero", 
                       "diaweighted_comcheck",
                       
                       "diacor_fit", "diacor_pred", "diacor_fit_ber", "diacor_pred_ber",
                       "diacor_tucker",
                       "diacor_correct", 
                       "diacor_nonzero", 
                       "diacor_comcheck", 
                       "sumdesigncor")




# 2. replication starts ####
set.seed(6397) # for original 11-may-2021 simulation, seed(256) used
model_seed <- sample(x = 1:100000, size = nrow(condition_df))
noise_seed <- sample(x = 1:100000, size = nrow(condition_df))

dim(condition_df)

unique(condition_df$signal_level)

# dividing the task such that each run takes equal number of high dimensionality problem #
repshigh <- rownames(subset(condition_df, dimension == "high"))[1:20]
repslow <- rownames(subset(condition_df, dimension == "low"))[1:20]

reps_to_do <- c(as.numeric(repslow), as.numeric(repshigh))

length(reps_to_do)


for (rrr in reps_to_do){ 
  
  modelseeding <- model_seed[rrr]
  noiseseeding <- noise_seed[rrr]
  
  cond <- condition_df[rrr,]
  
  
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
  
  
  # data generating #
  dat <- data_b(I = I, J = J, d1_active = d1_active, 
                c_active_block1 = c_active_block1, c_active_block2 = c_active_block2,
                d2_active = d2_active, comp_sd = rep(50, 3), 
                signal_level = cond$signal_level,
                strict_uncorrelated = F,
                Py = Py, seed = modelseeding)
  
  # scaling the data #
  X_train <- dat$X_train
  y_train <- dat$y_train
  
  X_train <- scale(X_train,T,T)
  
  X_test <- dat$X_test
  y_test <- dat$y_test
  
  X_test <- scale(X_test,T,T)
  
  truePx <- matrix(0, nrow = ncol(X_test), ncol = 3)
  
  truePx[d1_active,1] <- 1
  truePx[c_active_block1,2] <- 1
  truePx[c_active_block2,2] <- 1
  truePx[d2_active,3] <- 1
  
  irr <- which(abs(Py) < 1) 
  # index for the covariate that is irrelevant
  
  # splitting the truePx matrix for use for DIABLO #
  truePx_d1 <- truePx[1:(J/2), 1]
  truePx_c <- truePx[, 2]
  truePx_d2 <- truePx[(J/2 +1):J,3]
  
  # scdcovlogr #
  
  R <- 3 # number of covariates
  
  blockcols <- c(J/2, J/2)
  
  blockcols2 <- cumsum(blockcols)
  
  blockindex <- list()
  
  blockindex[[1]] <- 1:blockcols2[1] # blockindex is necessary input for scdcovlogr
  
  for (i in 2:length(blockcols)){
    blockindex[[i]] <- (blockcols2[i-1]+1):blockcols2[i]
  }
  
  
  # initial lasso glasso values per condition #
  if (cond$dimension == "low"){
    
    if (cond$signal_level == 0.8){
      if (cond$py_pattern == 1){
        glasso <- rep(0.5,3)
        lasso <- c(20,10,20)
      } else {
        glasso <- rep(0.5,3)
        lasso <- c(30,15,30)
      }
    }
    
    if (cond$signal_level == 0.5){
      if (cond$py_pattern == 1){
        glasso <- rep(0.5,3)
        lasso <- c(30,15,30)
      } else {
        glasso <- rep(0.5,3)
        lasso <- c(30,15,30)
      }
    }
    
    if (cond$signal_level == 0.2){
      if (cond$py_pattern == 1){
        glasso <- c(0.5,0.5,0.5)
        lasso <- c(30,15,30)
      } else {
        glasso <- c(5,5,5)
        lasso <- c(30,15,30)
      }
    }
    
  } else {
    
    if (cond$signal_level == 0.8){
      if (cond$py_pattern == 1){
        glasso <- rep(3,3)
        lasso <- c(15, 7.5, 15)
      } else {
        glasso <- c(1,1,1)
        lasso <- c(10,10,10)
      }
    }
    
    
    if (cond$signal_level == 0.5){
      if (cond$py_pattern == 1){
        glasso <- rep(2,3)
        lasso <- c(30,15,30)
      } else {
        glasso <- c(1,1,1)
        lasso <- c(30,15,30)
      }
    }
    
    
    if (cond$signal_level == 0.2){
      if (cond$py_pattern == 1){
        glasso <- rep(1,3)
        lasso <- c(20,10,20)
      } else {
        glasso <- c(1,1,1)
        lasso <- c(10,10,10)
      }
    }
    
  }
  
  
  # cross validation for alpha and ridge #
  time1 <- Sys.time()
  
  # ranges of alpha and ridge values considered
  alpha_range <- c(0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)
  
  ridge_range <-  c(50, 30, 10, 5, 1, 0.5)
  
  
  ranges <- expand.grid(alpha_range, ridge_range)
  
  colnames(ranges) <- c("alpha", "ridge")
  
  ranges <- ranges[order(ranges$alpha, decreasing = F),]
  
  scd_cv_results <- c()
  scd_cv_conv <- c()
  
  scd_cv_results <- matrix(NA, ncol = 3+2, nrow = nrow(ranges))
  
  colnames(ranges) <- c("alpha", "ridge")
  
  for (i in 1:nrow(ranges)){
    
    scdcv1 <- scd_cv(X = X_train, y = y_train, R = R, blockcols = c(J/2, J/2),
                     alpha = ranges[i,]$alpha, lasso = lasso, glasso = glasso, 
                     ridge_y = ranges[i,]$ridge, seed_cv = i+1, seed_method = 111,
                     inits = "rational", include_rational = TRUE, 
                     MAXITER = 500, nrFolds = 5, type_measure = "mse")
    
    
    
    scd_cv_results[i,] <- c(scdcv1$cve, scdcv1$se, sum(scdcv1$convergence),
                            ranges[i,]$alpha, ranges[i,]$ridge)
    
    print(i)
  }
  
  time2 <- Sys.time()
  
  cv_min_index <- which.min(scd_cv_results[,1]) # smallest CVE
  
  serule <- scd_cv_results[cv_min_index,1] + scd_cv_results[cv_min_index,2] # 1SE rule
  
  scd_possible <- scd_cv_results[scd_cv_results[,1] <= serule,]
  
  if (is.vector(scd_possible)){
    scd_chosen <- scd_possible
  } else {
    scd_chosen <- scd_possible[1,]
  }
  
  scd_alpha <- scd_chosen[4]
  scd_ridge <- scd_chosen[5]
  
  # scd without the second run to find the fitting lasso / glasso  (method B) #
  scd_initial_lasso <- scd_cov_logr(X = X_train, y = y_train, blockcols = c(J/2, J/2), 
                                    R = R, alpha = scd_alpha, lasso = lasso,
                                    glasso = glasso, 
                                    ridge_y = scd_ridge, inits = "rational", include_rational = TRUE,
                                    nrstart = 1, MAXITER = 3000, seed = noiseseeding,
                                    stop_value = 1e-5)
  
  
  # if the model does not converge, try other ridge and alpha values #
  if (scd_initial_lasso$result$iter > 2000){
    
    scd_possible_try <- scd_cv_results[order(scd_cv_results[,1], decreasing = F),]
    
    for (iii in 1:nrow(scd_possible_try)){
      
      scd_chosen_try <- scd_possible_try[iii,]
      
      scd_alpha_try <- scd_chosen_try[4]
      scd_ridge_try <- scd_chosen_try[5]
      
      # scd without the second lasso / glasso fix
      scd_initial_lasso_try <- scd_cov_logr(X = X_train, y = y_train, blockcols = c(J/2, J/2), 
                                            R = R, alpha = scd_alpha_try, lasso = lasso,
                                            glasso = glasso, 
                                            ridge_y = scd_ridge_try, inits = "rational", include_rational = TRUE,
                                            nrstart = 1, MAXITER = 2000, seed = noiseseeding,
                                            stop_value = 1e-5)
      
      print(iii)
      
      
      if (scd_initial_lasso_try$result$iter < 1000){
        scd_initial_lasso <- scd_initial_lasso_try
        break
      }
      
    }
    
    scd_alpha <- scd_alpha_try
    scd_ridge <- scd_ridge_try
    
  }
  
  
  # given the ridge and alpah value selected,
  # we figure out the best lasso and glasso that fit the population weights (method A)
  lasso_range <- c(3, 5, 10, 15, 20, 30, 50, 80)
  
  glasso_range <- c(0.5, 1, 2, 3, 5, 10)
  
  ranges <- expand.grid(lasso_range, glasso_range)
  
  colnames(ranges) <- c("lasso", "glasso")
  
  sparse_conv <- c()
  
  sparse_results <- matrix(NA, ncol = 6+2, nrow = nrow(ranges))
  
  for (i in 1:nrow(ranges)){
    
    scdlasso <- scd_cov_logr(X = X_train, y = y_train, blockcols = c(J/2, J/2), 
                             R = R, alpha = scd_alpha, lasso = c(ranges[i,]$lasso, ranges[i,]$lasso/2, ranges[i,]$lasso),
                             glasso = rep(ranges[i,]$glasso,3), 
                             ridge_y = scd_ridge, inits = "rational", include_rational = TRUE,
                             nrstart = 1, MAXITER = 1000, seed = noiseseeding,
                             stop_value = 1e-5)
    
    scdlasso_fit <- X_train %*% scdlasso$result$W %*% scdlasso$result$Py + scdlasso$result$py0
    
    scdlasso_fit_ber <- ber(pred = ((1 / (1 + exp(-scdlasso_fit))) > 0.5), true = y_train)
    
    scdlasso_correct <- corrects(estimate = scdlasso$result$W, defined = truePx)
    
    scdlasso_nonzero <- nonzero(estimate = scdlasso$result$W, defined = truePx)
    
    scdlasso_tucker <- tryCatch({
      RegularizedSCA::TuckerCoef(MatrixA = scdlasso$result$W, MatrixB = truePx)$tucker_value
    }, error = function(error_condition) {
      NA
    })
    
    scdlasso_zeros <- sum(scdlasso$result$W == 0)
    
    sparse_results[i,] <- c(scdlasso_fit_ber, scdlasso_correct, scdlasso_nonzero, 
                            scdlasso_tucker, scdlasso_zeros,
                            scdlasso$result$iter <= 500,
                            ranges[i,]$lasso, ranges[i,]$glasso)
    
    print(i)
  }
  
  
  colnames(sparse_results) <- c("fit", "corrects", "nonzeros", "tucker", "zeros", "conv", "lasso", "glasso")
  
  sparse_results <- as.data.frame(sparse_results)
  
  # according to the dimensionality condition, 
  # we select the lasso and group lasso values that are close to the defined population weights
  
  if (cond$dimension == "high"){
    
    sparse_results_490 <- sparse_results[sparse_results$zeros > 490,]
    
    sparse_results_490 <- sparse_results_490[sparse_results_490$zeros < 530,]
    
    if (nrow(sparse_results_490) == 0){
      
      sparse_results_490 <- sparse_results[sparse_results$zeros > 490,]
      
      sparse_results_490 <- sparse_results_490[complete.cases(sparse_results_490),]
      
      sparse_results_490_morezero <- sparse_results_490[order(sparse_results_490$zeros, decreasing = T),]
      
      sparse_results_490_ordered <- sparse_results_490_morezero[order(sparse_results_490_morezero$corrects, decreasing = T),]
      
      new_sparse_possible <- sparse_results_490_ordered[1,]
      
    } else {
      sparse_results_490 <- sparse_results_490[complete.cases(sparse_results_490),]
      
      sparse_results_490_morezero <- sparse_results_490[order(sparse_results_490$zeros, decreasing = T),]
      
      sparse_results_490_ordered <- sparse_results_490_morezero[order(sparse_results_490_morezero$corrects, decreasing = T),]
      
      new_sparse_possible <- sparse_results_490_ordered[1,]
      
    }
    
  } else {
    
    sparse_results_64 <- sparse_results[sparse_results$zeros > 64,]
    
    sparse_results_64 <- sparse_results_64[sparse_results_64$zeros < 85,]
    
    if (nrow(sparse_results_64) == 0){
      
      sparse_results_64 <- sparse_results[sparse_results$zeros > 64,]
      
      sparse_results_64 <- sparse_results_64[complete.cases(sparse_results_64),]
      
      sparse_results_64_morezero <- sparse_results_64[order(sparse_results_64$zeros, decreasing = T),]
      
      sparse_results_64_ordered <- sparse_results_64_morezero[order(sparse_results_64_morezero$corrects, decreasing = T),]
      
      new_sparse_possible <- sparse_results_64_ordered[1,]
      
    } else {
      
      sparse_results_64 <- sparse_results_64[complete.cases(sparse_results_64),]
      
      sparse_results_64_morezero <- sparse_results_64[order(sparse_results_64$zeros, decreasing = T),]
      
      sparse_results_64_ordered <- sparse_results_64_morezero[order(sparse_results_64_morezero$corrects, decreasing = T),]
      
      new_sparse_possible <- sparse_results_64_ordered[1,]
      
    }
    
    
  }
  
  
  
  scd_lasso <- c(new_sparse_possible$lasso, new_sparse_possible$lasso/2, new_sparse_possible$lasso)
  scd_glasso <- rep(new_sparse_possible$glasso,3)
  
  # actual estimation #
  scd1 <- scd_cov_logr(X = X_train, y = y_train, blockcols = c(J/2, J/2), 
                       R = R, alpha = scd_alpha, lasso = scd_lasso,
                       glasso = scd_glasso, ridge_y = scd_ridge, 
                       inits = "rational", include_rational = TRUE,
                       nrstart = 1, MAXITER = 5000, seed = noiseseeding,
                       stop_value = 1e-5)
  
  
  # ber and classifiation rate #
  scd_fit <- X_train %*% scd1$result$W %*% scd1$result$Py + scd1$result$py0
  
  scd_fit_ber <- ber(pred = ((1 / (1 + exp(-scd_fit))) > 0.5), true = y_train)
  
  scd_fit <- sum(((1 / (1 + exp(-scd_fit))) > 0.5)  == y_train)
  
  scd_pred <- X_test %*% scd1$result$W %*% scd1$result$Py + scd1$result$py0
  
  scd_pred_ber <- ber(pred = ((1 / (1 + exp(-scd_pred))) > 0.5), true = y_test)
  
  scd_pred <- sum(((1 / (1 + exp(-scd_pred))) > 0.5)  == y_test)
  
  # tucker congruence #
  tucks <- RegularizedSCA::TuckerCoef(truePx, scd1$result$W)
  
  scd_tucker <- tucks$tucker_value
  
  # nonzero #
  scd_nonzero <- nonzero(estimate = scd1$result$W, defined = truePx)
  
  scd_irr_nonzero <- nonzero(estimate = (scd1$result$W[,tucks$perm])[,irr],
                             defined = truePx[,irr])
  
  # correct classificaiton #
  scd_correct <- corrects(estimate = scd1$result$W, defined = truePx)
  
  scd_irr_correct <- corrects(estimate = (scd1$result$W[,tucks$perm])[,irr],
                              defined = truePx[,irr])
  
  # common and distinctive component identification #
  scd_comdis <- comdis(estimate = scd1$result$W)
  
  scd_d1 <- scd_comdis$d1
  scd_d2 <- scd_comdis$d2
  scd_common <- scd_comdis$common
  
  scd_wzeros <- apply(scd1$result$W == 0,2,sum)
  
  
  
  # second estimation - using initial lasso, glasso (method B) #
  scd2 <- scd_initial_lasso
  
  
  # ber and classifiation rate #
  scd2_fit <- X_train %*% scd2$result$W %*% scd2$result$Py + scd2$result$py0
  
  scd2_fit_ber <- ber(pred = ((1 / (1 + exp(-scd2_fit))) > 0.5), true = y_train)
  
  scd2_fit <- sum(((1 / (1 + exp(-scd2_fit))) > 0.5)  == y_train)
  
  scd2_pred <- X_test %*% scd2$result$W %*% scd2$result$Py + scd2$result$py0
  
  scd2_pred_ber <- ber(pred = ((1 / (1 + exp(-scd2_pred))) > 0.5), true = y_test)
  
  scd2_pred <- sum(((1 / (1 + exp(-scd2_pred))) > 0.5)  == y_test)
  
  # tucker congruence #
  tucks <- RegularizedSCA::TuckerCoef(truePx, scd2$result$W)
  
  scd2_tucker <- tucks$tucker_value
  
  # nonzero #
  scd2_nonzero <- nonzero(estimate = scd2$result$W, defined = truePx)
  
  scd2_irr_nonzero <- nonzero(estimate = (scd2$result$W[,tucks$perm])[,irr],
                              defined = truePx[,irr])
  
  # correct classificaiton #
  scd2_correct <- corrects(estimate = scd2$result$W, defined = truePx)
  
  scd2_irr_correct <- corrects(estimate = (scd2$result$W[,tucks$perm])[,irr],
                               defined = truePx[,irr])
  
  # common and distinctive component identification #
  scd2_comdis <- comdis(estimate = scd2$result$W)
  
  scd2_d1 <- scd2_comdis$d1
  scd2_d2 <- scd2_comdis$d2
  scd2_common <- scd2_comdis$common
  
  scd2_wzeros <- apply(scd2$result$W == 0,2,sum)
  
  
  # DIABLO-concatenated block.splsda (not reported in the paper) #
  
  Y_train <- mixOmics::unmap(y_train)
  Y_test <- mixOmics::unmap(y_test)
  
  # cross validation -> I don't do it because I know the oracle information #
  
  # first arranging the objects for the diablo function #
  design <- matrix(1, ncol = 1, nrow = 1)
  
  block1 <- data.frame(X_train)
  
  plsX <- list(block1 = block1)
  
  if (cond$dimension == "low"){
    list.keepX = list(block1 = c(4,8,4))
    
  } else {
    list.keepX = list(block1 = c(25, 50, 25)) # must check #
  }
  
  # actual model fitting
  diacon1 <- block.splsda(X = plsX, Y = y_train,
                          ncomp = R, keepX = list.keepX, design = design)
  
  # fit #
  diacon_fit <- predict(object = diacon1, newdata = plsX, dist = "max.dist")
  
  diacon_fit$MajorityVote$max.dist[,R]
  
  diacon_fit_ber <- ber(pred =  diacon_fit$MajorityVote$max.dist[,R], y_train)
  
  diacon_fit <- sum(y_train == diacon_fit$MajorityVote$max.dist[,R])
  
  # ber and classification #
  block1test <- data.frame(X_test)
  
  plsXtest <- list(block1 = block1test)
  
  diacon_pred <- predict(object = diacon1,
                         newdata = plsXtest,
                         dist = "max.dist")
  
  diacon_pred$MajorityVote$max.dist[,R]
  
  diacon_pred_ber <- ber(pred = diacon_pred$MajorityVote$max.dist[,R],
                         true = y_test)
  
  diacon_pred <- sum(y_test == diacon_pred$MajorityVote$max.dist[,R])
  
  # tucker congruence #
  diacon_tucks <- RegularizedSCA::TuckerCoef(truePx,
                                             diacon1$loadings$block1)
  
  diacon_tucker <- diacon_tucks$tucker_value
  
  # correct classification #
  diacon_correct <- corrects(estimate = diacon1$loadings$block1,
                             defined = truePx)
  
  diacon_irr_correct <- corrects(
    estimate = (diacon1$loadings$block1[,diacon_tucks$perm])[,irr],
    defined = truePx[,irr])
  
  # nonzero classification #
  diacon_nonzero <- nonzero(estimate = diacon1$loadings$block1,
                            defined = truePx)
  
  diacon_irr_nonzero <- nonzero(
    estimate = (diacon1$loadings$block1[,diacon_tucks$perm])[,irr],
    defined = truePx[,irr])
  
  # common and distinctive component identification #
  diacon_comdis <- comdis(estimate = diacon1$loadings$block1)
  diacon_d1 <- diacon_comdis$d1
  diacon_d2 <- diacon_comdis$d2
  diacon_common <- diacon_comdis$common
  
  # DIABLO 3 BLOCKS block.splsda (reported in the paper) #
  
  # no cross validation as we know oracle sparsity structure #
  
  # first arranging the objects for the diablo function #
  block1_3 <- data.frame(X_train[,1:(J/2)]) # first datablock
  block2_3 <- data.frame(X_train[,(J/2 + 1):J]) # second datablock
  block3_3 <- data.frame(X_train) # superblock
  
  
  plsX_3 <- list(block1_3 = block1_3, 
                 block2_3 = block2_3,
                 block3_3 = block3_3)
  
  # design matrix specification #
  
  # 1. null design (not reported in the paper) #
  design_null <- matrix(0, ncol = 3, nrow = 3)
  diag(design_null) <- 0
  
  # 2. full weighted design (not reported in the paper)
  design_weighted <- matrix(0.1, ncol = 3, nrow = 3)
  diag(design_weighted) <- 0
  
  # 3. Through PLS and observing pairwise correlation (reported in the paper)
  design_full <- matrix(1, ncol = 3, nrow = 3)
  diag(design_full) <- 0
  
  pls_solution <- block.plsda(X = plsX_3, Y = y_train, ncomp = c(1), design = design_full)
  
  pls_design_tmat <- matrix(unlist(pls_solution$variates), ncol = 4, byrow = TRUE)
  
  pls_design_cor <-cor(pls_design_tmat[,-4])
  
  diag(pls_design_cor) <- 0
  
  pls_design_cor <- abs(pls_design_cor)
  
  design_threshold <- pls_design_cor > 0.8
  
  pls_design_cor[] <- 0 
  
  pls_design_cor[design_threshold] <- 1
  
  design_cor <- pls_design_cor
  
  sumdesigncor <- sum(design_cor)
  
  if (cond$dimension == "low"){
    list.keepX_3 = list(block1_3 = c(4), 
                        block2_3 = c(4),
                        block3_3 = c(8))
  } else {
    list.keepX_3 = list(block1_3 = c(25), 
                        block2_3 = c(25),
                        block3_3 = c(40))
    
  }
  
  # Null design (not reported in the paper) #
  diathree1 <- block.splsda(X = plsX_3, Y = y_train,
                            ncomp = c(1), keepX = list.keepX_3, 
                            design = design_null)
  
  # fit #
  diathree_fit <- predict(object = diathree1, newdata = plsX_3, dist = "max.dist")
  
  diathree_fit$WeightedVote$max.dist[,1]
  
  diathree_fit_ber <- ber(pred = diathree_fit$WeightedVote$max.dist[,1], true = y_train)
  
  diathree_fit <- sum(y_train == diathree_fit$WeightedVote$max.dist[,1])
  
  # test set classification #
  block1test_3 <- data.frame(X_test[,1:(J/2)])
  block2test_3 <- data.frame(X_test[,(J/2 + 1):J])
  block3test_3 <- data.frame(X_test)
  
  plsXtest_3 <- list(block1_3 = block1test_3,
                     block2_3 = block2test_3,
                     block3_3 = block3test_3)
  
  diathree_pred <- predict(object = diathree1, 
                           newdata = plsXtest_3, 
                           dist = "max.dist")
  
  diathree_pred$WeightedVote$max.dist[,1]
  
  diathree_pred_ber <- ber(pred = diathree_pred$WeightedVote$max.dist[,1], 
                           true = y_test)
  
  diathree_pred <- sum(y_test == diathree_pred$WeightedVote$max.dist[,1])
  
  
  # calculate the correct classification and the tucker congruence 
  # for all of the distinctive components,
  # and then select which is the one to be reported, based on the congruence 
  
  # all of the congruences for distinctive components
  congs <- c(congruence(truePx_d1, diathree1$loadings$block1_3),
             congruence(truePx_d2, diathree1$loadings$block1_3),
             congruence(truePx_d1, diathree1$loadings$block2_3),
             congruence(truePx_d2, diathree1$loadings$block2_3))
  
  diathree_corrects <- c(corrects(estimate = c(diathree1$loadings$block1_3), defined = truePx_d1),
                         corrects(estimate = c(diathree1$loadings$block1_3), defined = truePx_d2),
                         corrects(estimate = c(diathree1$loadings$block2_3), defined = truePx_d1),
                         corrects(estimate = c(diathree1$loadings$block2_3), defined = truePx_d2))
  
  diathree_nonzeros <- c(nonzero(estimate = c(diathree1$loadings$block1_3), defined = truePx_d1),
                         nonzero(estimate = c(diathree1$loadings$block1_3), defined = truePx_d2),
                         nonzero(estimate = c(diathree1$loadings$block2_3), defined = truePx_d1),
                         nonzero(estimate = c(diathree1$loadings$block2_3), defined = truePx_d2))
  
  # the order #
  # 1. dia1 vs true1
  # 2. dia1 vs true2
  # 3. dia2 vs true1
  # 4. dia2 vs true2
  
  # determine which matches are the correct matches:
  matched <- which.max(congs)
  
  # if the best match is match no. 4; dia2 vs true2,
  # then the other match should be match no. 1: dia 1 vs true 1.
  
  # this code figures that out
  if (matched > 2){
    dia1_matched <- 5 - matched
    dia2_matched <- matched
  } else {
    dia1_matched <- matched
    dia2_matched <- (matched - 5) * -1
  }
  
  # tucker congruence #
  diathree_congs_matched <- congs[c(dia1_matched, dia2_matched)]
  
  # common component
  diathree_congs_matched[3] <- congruence(truePx_c, diathree1$loadings$block3_3)
  
  diathree_tucker <- mean(abs(diathree_congs_matched))
  
  # correct classification #
  diathree_corrects_matched <- diathree_corrects[c(dia1_matched, dia2_matched)]
  
  # common coponent
  diathree_corrects_matched[3] <- corrects(estimate = truePx_c, diathree1$loadings$block3_3)
  
  diathree_correct <- mean(diathree_corrects_matched)
  
  # nonzero classification #
  diathree_nonzeros_matched <- diathree_nonzeros[c(dia1_matched, dia2_matched)]
  
  # common coponent
  diathree_nonzeros_matched[3] <- nonzero(estimate = truePx_c, diathree1$loadings$block3_3)
  
  diathree_nonzero <- mean(diathree_nonzeros_matched)
  
  # there is no need for doing comdis for diathree, as we just spoonfeed the answer for them. 
  diathree_comcheck <- T
  
  if(sum(diathree1$loadings$block3_3[1:(J/2)] != 0) == 0) {diathree_comcheck <- F}
  if(sum(diathree1$loadings$block3_3[(J/2 + 1):J] != 0) == 0) {diathree_comcheck <- F}
  
  
  
  
  ## 2. full weighted design (not reported in the paper)  ##
  diaweighted <- block.splsda(X = plsX_3, Y = y_train,
                              ncomp = c(1), keepX = list.keepX_3, 
                              design = design_weighted)
  
  
  # full weighted deign: fit #
  diaweighted_fit <- predict(object = diaweighted, newdata = plsX_3, dist = "max.dist")
  
  diaweighted_fit$WeightedVote$max.dist[,1]
  
  diaweighted_fit_ber <- ber(pred = diaweighted_fit$WeightedVote$max.dist[,1], true = y_train)
  
  diaweighted_fit <- sum(y_train == diaweighted_fit$WeightedVote$max.dist[,1])
  
  # test set classification #
  block1test_3 <- data.frame(X_test[,1:(J/2)])
  block2test_3 <- data.frame(X_test[,(J/2 + 1):J])
  block3test_3 <- data.frame(X_test)
  
  plsXtest_3 <- list(block1_3 = block1test_3,
                     block2_3 = block2test_3,
                     block3_3 = block3test_3)
  
  diaweighted_pred <- predict(object = diaweighted, 
                              newdata = plsXtest_3, 
                              dist = "max.dist")
  
  diaweighted_pred$WeightedVote$max.dist[,1]
  
  diaweighted_pred_ber <- ber(pred = diaweighted_pred$WeightedVote$max.dist[,1], 
                              true = y_test)
  
  diaweighted_pred <- sum(y_test == diaweighted_pred$WeightedVote$max.dist[,1])
  
  
  # calculate the correct classification and the tucker congruence 
  # for all of the distinctive components,
  # and then select which is the one to be reported, based on the congruence 
  
  # all of the congruences for distinctive components
  diaweighted_congs <- c(congruence(truePx_d1, diaweighted$loadings$block1_3),
                         congruence(truePx_d2, diaweighted$loadings$block1_3),
                         congruence(truePx_d1, diaweighted$loadings$block2_3),
                         congruence(truePx_d2, diaweighted$loadings$block2_3))
  
  diaweighted_corrects <- c(corrects(estimate = c(diaweighted$loadings$block1_3), defined = truePx_d1),
                            corrects(estimate = c(diaweighted$loadings$block1_3), defined = truePx_d2),
                            corrects(estimate = c(diaweighted$loadings$block2_3), defined = truePx_d1),
                            corrects(estimate = c(diaweighted$loadings$block2_3), defined = truePx_d2))
  
  diaweighted_nonzeros <- c(nonzero(estimate = c(diaweighted$loadings$block1_3), defined = truePx_d1),
                            nonzero(estimate = c(diaweighted$loadings$block1_3), defined = truePx_d2),
                            nonzero(estimate = c(diaweighted$loadings$block2_3), defined = truePx_d1),
                            nonzero(estimate = c(diaweighted$loadings$block2_3), defined = truePx_d2))
  
  # the order #
  # 1. dia1 vs true1
  # 2. dia1 vs true2
  # 3. dia2 vs true1
  # 4. dia2 vs true2
  
  # determine which matches are the correct matches:
  diaweighted_matched <- which.max(diaweighted_congs)
  
  # if the best match is match no. 4; dia2 vs true2,
  # then the other match should be match no. 1: dia 1 vs true 1.
  
  # this code figures that out
  if (diaweighted_matched > 2){
    diaweighted1_matched <- 5 - matched
    diaweighted2_matched <- matched
  } else {
    diaweighted1_matched <- matched
    diaweighted2_matched <- (matched - 5) * -1
  }
  
  # tucker congruence #
  diaweighted_congs_matched <- congs[c(dia1_matched, dia2_matched)]
  
  # common component
  diaweighted_congs_matched[3] <- congruence(truePx_c, diaweighted$loadings$block3_3)
  
  diaweighted_tucker <- mean(abs(diaweighted_congs_matched))
  
  # correct classification #
  diaweighted_corrects_matched <- diaweighted_corrects[c(dia1_matched, dia2_matched)]
  
  # common coponent
  diaweighted_corrects_matched[3] <- corrects(estimate = truePx_c, diaweighted$loadings$block3_3)
  
  diaweighted_correct <- mean(diaweighted_corrects_matched)
  
  # nonzero classification #
  diaweighted_nonzeros_matched <- diaweighted_nonzeros[c(dia1_matched, dia2_matched)]
  
  # common coponent
  diaweighted_nonzeros_matched[3] <- nonzero(estimate = truePx_c, diaweighted$loadings$block3_3)
  
  diaweighted_nonzero <- mean(diaweighted_nonzeros_matched)
  
  # there is no need for doing comdis for diathree, as we just spoonfeed the answer for them. 
  diaweighted_comcheck <- T
  
  if(sum(diaweighted$loadings$block3_3[1:(J/2)] != 0) == 0) {diaweighted_comcheck <- F}
  if(sum(diaweighted$loadings$block3_3[(J/2 + 1):J] != 0) == 0) {diaweighted_comcheck <- F}
  
  
  
  
  ## 3. design through pre-PLS (reported in the paper)  ##
  diacor <- block.splsda(X = plsX_3, Y = y_train,
                         ncomp = c(1), keepX = list.keepX_3, 
                         design = design_cor)
  
  
  # full weighted deign: fit #
  diacor_fit <- predict(object = diacor, newdata = plsX_3, dist = "max.dist")
  
  diacor_fit$WeightedVote$max.dist[,1]
  
  diacor_fit_ber <- ber(pred = diacor_fit$WeightedVote$max.dist[,1], true = y_train)
  
  diacor_fit <- sum(y_train == diacor_fit$WeightedVote$max.dist[,1])
  
  # test set classification #
  block1test_3 <- data.frame(X_test[,1:(J/2)])
  block2test_3 <- data.frame(X_test[,(J/2 + 1):J])
  block3test_3 <- data.frame(X_test)
  
  plsXtest_3 <- list(block1_3 = block1test_3,
                     block2_3 = block2test_3,
                     block3_3 = block3test_3)
  
  diacor_pred <- predict(object = diacor, 
                         newdata = plsXtest_3, 
                         dist = "max.dist")
  
  diacor_pred$WeightedVote$max.dist[,1]
  
  diacor_pred_ber <- ber(pred = diacor_pred$WeightedVote$max.dist[,1], 
                         true = y_test)
  
  diacor_pred <- sum(y_test == diacor_pred$WeightedVote$max.dist[,1])
  
  
  # calculate the correct classification and the tucker congruence 
  # for all of the distinctive components,
  # and then select which is the one to be reported, based on the congruence 
  
  # all of the congruences for distinctive components
  diacor_congs <- c(congruence(truePx_d1, diacor$loadings$block1_3),
                    congruence(truePx_d2, diacor$loadings$block1_3),
                    congruence(truePx_d1, diacor$loadings$block2_3),
                    congruence(truePx_d2, diacor$loadings$block2_3))
  
  diacor_corrects <- c(corrects(estimate = c(diacor$loadings$block1_3), defined = truePx_d1),
                       corrects(estimate = c(diacor$loadings$block1_3), defined = truePx_d2),
                       corrects(estimate = c(diacor$loadings$block2_3), defined = truePx_d1),
                       corrects(estimate = c(diacor$loadings$block2_3), defined = truePx_d2))
  
  diacor_nonzeros <- c(nonzero(estimate = c(diacor$loadings$block1_3), defined = truePx_d1),
                       nonzero(estimate = c(diacor$loadings$block1_3), defined = truePx_d2),
                       nonzero(estimate = c(diacor$loadings$block2_3), defined = truePx_d1),
                       nonzero(estimate = c(diacor$loadings$block2_3), defined = truePx_d2))
  
  # the order #
  # 1. dia1 vs true1
  # 2. dia1 vs true2
  # 3. dia2 vs true1
  # 4. dia2 vs true2
  
  # determine which matches are the correct matches:
  diacor_matched <- which.max(diacor_congs)
  
  # if the best match is match no. 4; dia2 vs true2,
  # then the other match should be match no. 1: dia 1 vs true 1.
  
  # this code figures that out
  if (diacor_matched > 2){
    diacor1_matched <- 5 - matched
    diacor2_matched <- matched
  } else {
    diacor1_matched <- matched
    diacor2_matched <- (matched - 5) * -1
  }
  
  # tucker congruence #
  diacor_congs_matched <- congs[c(dia1_matched, dia2_matched)]
  
  # common component
  diacor_congs_matched[3] <- congruence(truePx_c, diacor$loadings$block3_3)
  
  diacor_tucker <- mean(abs(diacor_congs_matched))
  
  # correct classification #
  diacor_corrects_matched <- diacor_corrects[c(dia1_matched, dia2_matched)]
  
  # common coponent
  diacor_corrects_matched[3] <- corrects(estimate = truePx_c, diacor$loadings$block3_3)
  
  diacor_correct <- mean(diacor_corrects_matched)
  
  # nonzero classification #
  diacor_nonzeros_matched <- diacor_nonzeros[c(dia1_matched, dia2_matched)]
  
  # common coponent
  diacor_nonzeros_matched[3] <- nonzero(estimate = truePx_c, diacor$loadings$block3_3)
  
  diacor_nonzero <- mean(diacor_nonzeros_matched)
  
  # there is no need for doing comdis for diathree, as we just spoonfeed the answer for them. 
  diacor_comcheck <- T
  
  if(sum(diacor$loadings$block3_3[1:(J/2)] != 0) == 0) {diacor_comcheck <- F}
  if(sum(diacor$loadings$block3_3[(J/2 + 1):J] != 0) == 0) {diacor_comcheck <- F}
  
  
  
  
  
  # recording the results #
  
  performances <- t(c(scd_alpha, scd_ridge,
                      scd_fit, scd_pred, scd_fit_ber, scd_pred_ber,
                      scd_tucker,
                      scd_correct, 
                      scd_nonzero,
                      scd_d1, scd_d2, scd_common,
                      scd_lasso[1], scd_glasso[1],
                      scd_wzeros[1], scd_wzeros[2], scd_wzeros[3],
                      
                      scd_alpha, scd_ridge,
                      scd2_fit, scd2_pred, scd2_fit_ber, scd2_pred_ber,
                      scd2_tucker,
                      scd2_correct, 
                      scd2_nonzero,
                      scd2_d1, scd2_d2, scd2_common,
                      lasso[1], lasso[2], glasso[1],
                      scd2_wzeros[1], scd2_wzeros[2], scd2_wzeros[3],
                      
                      
                      diacon_fit, diacon_pred, diacon_fit_ber, diacon_pred_ber,
                      diacon_tucker,
                      diacon_correct, 
                      diacon_nonzero, 
                      diacon_d1, diacon_d2, diacon_common,
                      
                      diathree_fit, diathree_pred, diathree_fit_ber, diathree_pred_ber,
                      diathree_tucker,
                      diathree_correct, 
                      diathree_nonzero, 
                      diathree_comcheck,
                      
                      diaweighted_fit, diaweighted_pred, diaweighted_fit_ber, diaweighted_pred_ber,
                      diaweighted_tucker,
                      diaweighted_correct, 
                      diaweighted_nonzero, 
                      diaweighted_comcheck, 
                      
                      diacor_fit, diacor_pred, diacor_fit_ber, diacor_pred_ber,
                      diacor_tucker,
                      diacor_correct, 
                      diacor_nonzero, 
                      diacor_comcheck, sumdesigncor
  ))
  
  
  resulting <- data.frame(cond$dimension, cond$signal_level, cond$py_pattern,
                          cond$reps, modelseeding, noiseseeding)
  
  resulting <- cbind(resulting, performances)
  
  # recording the current result in the data frame that I've defined above
  results[rrr,] <- resulting
  
  # if there was a problem, stop replication
  if(anyNA(results[rrr,])){
    stop("NA resulted")
  }
  
  print(rep(rrr, 10))
  
  # saving the current "results" object, in case it crashes in the middle 
  save("results", file = "./simulation 11-may-2021/results1_sup_11_may_2021.Rdata")
  write.table(x = results, file = "./simulation 11-may-2021/results1_sup_11_may_2021.txt", sep=",")
  
  flush.console()
  
}

