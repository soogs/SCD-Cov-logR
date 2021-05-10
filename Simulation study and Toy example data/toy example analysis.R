# toy example dataset analysis #
# initiated: 17-april-2021 #
# last modified (comments only): 10-may-2021

# 0. Introduction ####
# For the toy example dataset, SCD-Cov-logR and the other methods are administered
# model selection and final analysis are all included in this Rscript.


# 1. Source functions, load packages ####
setwd("E:\\Users\\park\\Desktop\\Final functions for Github") # TO BE CHANGED BY USERS

# data generation function (for toy example and simulation study)
source("./3cov_datagen.R") 

Rcpp::sourceCpp("./updateW_diffpen.cpp")

source("./scd_cov_logr.R")

source("./scd_cv.R")

# SCaDS by de Schipper and Van Deun (Zeitschrift für Psychologie, 2019)
# this function can be found in the Github repository of Niek de Schipper: https://github.com/trbKnl/SCaDS
Rcpp::sourceCpp("../toy example/niek spca/mmsca.cpp")

library(mixOmics)
library(psych)
library(data.table)
library(semPlot)
library(regsem)
library(lavaan)
library(nFactors)
library(glmnet)

# function to be used to evaluate the loadings matrices provided by RegSEM
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

X_train <- scale(X_train,T,T)

X_test <- dat$X_test
y_test <- dat$y_test

X_test <- scale(X_test,T,T)


# scree test to determine the number of covariates: three
svdd <- svd(X_train)

nsc <- nScree(svdd$d)
nsc2 <- nScree(svdd$d^2)
  

# 3. Lavaan and RegSEM ####

# I do not conduct ESEM, but I merely use the package to make variable names,
# so that I can apply Lavaan and RegSEM

# giving the data matrix column names for SEM
thirty <- 1:30

x_thirty <- paste("x", thirty, sep = "")

colnames(X_train) <- x_thirty

dat_train <- cbind(X_train, y_train)

colnames(dat_train)[ncol(dat_train)] <- c("y")

colnames(dat_train)


colnames(X_test) <- x_thirty

dat_test <- cbind(X_test, y_test)

colnames(dat_test)[ncol(dat_test)] <- "y"

colnames(dat_test)



esem_efa <- fa(X_train, nfactors =3,
               fm = 'ml')

esem_efa$loadings

esem_loadings <- data.table(matrix(round(esem_efa$loadings, 2),
                                   nrow = 30, ncol = 3))

names(esem_loadings) <- c("F1","F2", "F3")

esem_loadings$item <- paste0("x", c(1:30))

loadings1 <- paste0("p", c(1:30), "_1")
loadings2 <- paste0("p", c(1:30), "_2")
loadings3 <- paste0("p", c(1:30), "_3")

esem_loadings <- melt(esem_loadings, "item", variable.name = "latent")

esem_loadings

esem_loadings[, loadingindex := c(loadings1, loadings2, loadings3)]

anchors <- c(F1 = "x1", F2 = "x5", F3 = "x20") 
# from the truePx structure, i know that these do not have cross loadings

make_esem_model <- function (loadings_dt, anchors){
  
  # make is_anchor variable
  loadings_dt[, is_anchor := 0]
  for (l in names(anchors)) loadings_dt[latent != l & item == anchors[l], is_anchor := 1]
  
  # make syntax column per item; syntax is different depending on is_anchor
  loadings_dt[is_anchor == 0, syntax := paste0(loadingindex, "*", item)]
  loadings_dt[is_anchor == 1, syntax := paste0(value,"*", item)]
  
  #Make syntax for each latent variable
  each_syntax <- function (l){
    paste(l, "=~", paste0(loadings_dt[latent == l, syntax], collapse = "+"),"\n")
  }
  
  # Put all syntaxes together
  paste(sapply(unique(loadings_dt$latent), each_syntax), collapse = " ")
}


#make model
esem_model <- make_esem_model(esem_loadings, anchors)
#print model
writeLines(esem_model)

# following model is specified finally:
# the variance of the factors are fixed at 1
# x1 is fixed at 0 for F2 and F3
# x5 is fixed at 0 for F1 and F3
# x20 is fixed at 0 for F1 and F2
# in total: 6 loadings and 3 variance are fixed

esem_model <-  'F1 =~ NA*x1+p2_1*x2+p3_1*x3+p4_1*x4+ 0*x5+ p6_1*x6+p7_1*x7+p8_1*x8+p9_1*x9+p10_1*x10+p11_1*x11+p12_1*x12+p13_1*x13+p14_1*x14+p15_1*x15+p16_1*x16+p17_1*x17+p18_1*x18+p19_1*x19+ 0*x20 +p21_1*x21+p22_1*x22+p23_1*x23+p24_1*x24+p25_1*x25+p26_1*x26+p27_1*x27+p28_1*x28+p29_1*x29+p30_1*x30 
F2 =~ 0*x1+p2_2*x2+p3_2*x3+p4_2*x4+ p5_2*x5 +p6_2*x6+p7_2*x7+p8_2*x8+p9_2*x9+p10_2*x10+p11_2*x11+p12_2*x12+p13_2*x13+p14_2*x14+p15_2*x15+p16_2*x16+p17_2*x17+p18_2*x18+p19_2*x19+ 0*x20+p21_2*x21+p22_2*x22+p23_2*x23+p24_2*x24+p25_2*x25+p26_2*x26+p27_2*x27+p28_2*x28+p29_2*x29+p30_2*x30 
F3 =~ 0*x1+p2_3*x2+p3_3*x3+p4_3*x4+ 0*x5+p6_3*x6+p7_3*x7+p8_3*x8+p9_3*x9+p10_3*x10+p11_3*x11+p12_3*x12+p13_3*x13+p14_3*x14+p15_3*x15+p16_3*x16+p17_3*x17+p18_3*x18+p19_3*x19+ p20_3*x20 +p21_3*x21+p22_3*x22+p23_3*x23+p24_3*x24+p25_3*x25+p26_3*x26+p27_3*x27+p28_3*x28+p29_3*x29+p30_3*x30 
 
y ~ F1
y ~ F2
y ~ F3

F1 ~~ 1*F1
F2 ~~ 1*F2
F3 ~~ 1*F3
' 

# Lavaan: SEM #
lav.out <- sem(esem_model, dat_train) 
  
  
coef(lav.out)[1:84]

lav_T <- predict(object = lav.out, newdata = dat_train)

lav_T %*% coef(lav.out)[85:87] > 0

(lav_T %*% coef(lav.out)[85:87] > 0) == y_train

sum((lav_T %*% coef(lav.out)[85:87] > 0) == y_train) # in-sample prediction

# predict on a test set
lav_test <- predict(object = lav.out, newdata = dat_test)

lav_test %*% coef(lav.out)[85:87] > 0

(lav_test %*% coef(lav.out)[85:87] > 0) == y_test

sum((lav_test %*% coef(lav.out)[85:87] > 0) == y_test) # out-sample prediction

# lavaan performs very well for regression fit and prediction. However, of course the loadings are not sparse


# RegSEM #
time1 <- Sys.time()
  
# I try many different lambda values (lasso penalties), to determine which lambda value leads to the 
# loadings that are closest to the defined parameters

reg.out <-cv_regsem(lav.out, n.lambda= 50, lambda.start = 0.50, type="lasso",jump=0.01, pars_pen="loadings", alpha = 0, gamma = 0, ncore = 20)

time2 <- Sys.time()

timetaken <- as.numeric(time2 - time1)

# determining which lambda value leads to the best-reflecting solutions
regsemcorrect <- c()

for (iii in 1:nrow(reg.out$fits)){
  
  regsem_coef <- reg.out$parameters[iii,]
  
  names(regsem_coef)
  
  # making the loadings into a P matrix
  regsem_P <- matrix(NA, nrow = 30, ncol = 3)
  
  for (i in 1:length(regsem_coef[1:84])){
    
    rowpart <- sapply(strsplit(names(regsem_coef)[i], " "), "[[", 3)
    
    columnpart <- sapply(strsplit(names(regsem_coef)[i], " "), "[[", 1)
    
    rowpart <- substr(rowpart, start = 2, stop = 3)
    
    columnpart <- substr(columnpart, start = 2, stop = 2)
    
    rowpart <- as.numeric(rowpart)
    
    columnpart <- as.numeric(columnpart)
    
    regsem_P[rowpart, columnpart] <- regsem_coef[i]
    
  }
  
  (regsem_P[is.na(regsem_P)] <- 0)
  
  regsemcorrect[iii] <- corrects(estimate = regsem_P, defined = truePx)
  
  
}


regsem_chosen <- which.max(regsemcorrect)

lambda_chosen <-  reg.out$fits[regsem_chosen,][1]

# making the RegSEM loadings into a P matrix
regsem_coef <- reg.out$parameters[regsem_chosen,]

regsem_P <- matrix(NA, nrow = 30, ncol = 3)

for (i in 1:length(regsem_coef[1:84])){
  
  rowpart <- sapply(strsplit(names(regsem_coef)[i], " "), "[[", 3)
  
  columnpart <- sapply(strsplit(names(regsem_coef)[i], " "), "[[", 1)
  
  rowpart <- substr(rowpart, start = 2, stop = 3)
  
  columnpart <- substr(columnpart, start = 2, stop = 2)
  
  rowpart <- as.numeric(rowpart)
  
  columnpart <- as.numeric(columnpart)
  
  regsem_P[rowpart, columnpart] <- regsem_coef[i]
  
}

(regsem_P[is.na(regsem_P)] <- 0)

regsem_P
  
  
  
regsem_T <- X_train %*% regsem_P %*% solve(t(regsem_P) %*% regsem_P) 

reg.out$final_pars[85:87] # regression coefficients

regsem_fit <- regsem_T %*% reg.out$final_pars[85:87]

regsem_fit <- sum((regsem_fit > 0) == y_train)

# SEM pred # 
X_test %*% regsem_P %*% solve(t(regsem_P) %*% regsem_P) 

regsem_pred <- X_test %*% regsem_P %*% solve(t(regsem_P) %*% regsem_P)  %*% reg.out$final_pars[85:87] > 0 

regsem_pred <- sum(regsem_pred == y_test)

# performance of in-sample and out-sample prediction are very good,
# but the loadings do not exactly reflect the defined weights
  
  
# 4. SCD-Cov-logR #####

R <- 3

blockcols <- c(J/2, J/2)

blockcols2 <- cumsum(blockcols)

blockindex <- list()

blockindex[[1]] <- 1:blockcols2[1]

for (i in 2:length(blockcols)){
  blockindex[[i]] <- (blockcols2[i-1]+1):blockcols2[i]
}

svdd <- svd(X_train)

nScree(svdd$d)

# a random try #
scd1 <- scd_cov_logr(X = X_train, y = y_train, blockcols = c(J/2, J/2), 
                                 R = 3, alpha = 0.3, lasso = c(40/2,40, 40),
                                 glasso = rep(1,3), ridge_y = 10, inits = "rational",
                                 nrstart = 1, MAXITER = 5000, seed = 1, include_rational = TRUE, 
                                 stop_value = 1e-5)

scd1$result$W


# cross validation for alpha and ridge #
alpha_range <- c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)

ridge_range <-  c(50, 30, 10, 5, 3, 1, 0.5, 0.1)

ranges <- expand.grid(alpha_range, ridge_range)

colnames(ranges) <- c("alpha", "ridge")

scd_cv_results <- c()
scd_cv_conv <- c()

scd_cv_results <- matrix(NA, ncol = 3+2, nrow = nrow(ranges))

for (i in 1:nrow(ranges)){
  
  scdcv1 <- scd_cv(X = X_train, y = y_train, R = R, blockcols = c(J/2, J/2),
                   alpha = ranges[i,]$alpha, lasso = rep(0,3), glasso = rep(0,3), 
                   ridge_y = ranges[i,]$ridge,
                   MAXITER = 500, nrFolds = 5, seed_cv = i+1, type_measure = "mse",
                   inits = "rational", seed_method = i+1, include_rational = TRUE)
  
  
  scd_cv_results[i,] <- c(scdcv1$cve, scdcv1$se, sum(scdcv1$convergence),
                          ranges[i,]$alpha, ranges[i,]$ridge)
  
  print(i)
}

time2 <- Sys.time()

# choosing the alpha and ridge through 1 standard error rule #
alpha_ridge_cv_results <- scd_cv_results

scd_cv_results <- alpha_ridge_cv_results

cv_min_index <- which.min(scd_cv_results[,1])

scd_cv_results[which.min(scd_cv_results[,1]),]

serule <- scd_cv_results[cv_min_index,1] + scd_cv_results[cv_min_index,2]

scd_possible <- scd_cv_results[scd_cv_results[,1] < serule,]

scd_possible <- scd_possible[order(scd_possible[,4], decreasing = F),]

if (is.vector(scd_possible)){
  scd_chosen <- scd_possible
} else {
  scd_chosen <- scd_possible[1,]
}

scd_alpha <- scd_chosen[4]
scd_ridge <- scd_chosen[5]
  
  
# cross validation for lasso and glasso 
lasso_range <- c(100, 45, 30, 15, 10, 7, 5, 1, 0.5)

glasso_range <-  c(10, 5, 2, 1, 0.5, 0.1)

ranges <- expand.grid(lasso_range, glasso_range)

colnames(ranges) <- c("lasso", "glasso")

scd_cv_results <- c()
scd_cv_conv <- c()

scd_cv_results <- matrix(NA, ncol = 3+2, nrow = nrow(ranges))


for (i in 1:nrow(ranges)){
  
  scdcv1 <- scd_cv(X = X_train, y = y_train, R = R, blockcols = c(J/2, J/2),
                   alpha = scd_alpha, lasso = c(ranges[i,]$lasso[1]/2, ranges[i,]$lasso[1], ranges[i,]$lasso[1]), 
                   glasso = rep(ranges[i,]$glasso,3), 
                   ridge_y = scd_ridge,
                   MAXITER = 500, nrFolds = 5, seed_cv  = i+1, seed_method = 1+1, 
                   inits = "rational", include_rational = TRUE, 
                   type_measure = "mse")
  
  
  scd_cv_results[i,] <- c(scdcv1$cve, scdcv1$se, sum(scdcv1$convergence),
                          ranges[i,]$lasso, ranges[i,]$glasso)
  
  print(i)
}

time2 <- Sys.time()

lasso_glasso_cv_results <- scd_cv_results

# scd_cv_results <- lasso_glasso_cv_results

cv_min_index <- which.min(scd_cv_results[,1])

serule <- scd_cv_results[cv_min_index,1] + scd_cv_results[cv_min_index,2]

scd_possible <- scd_cv_results[scd_cv_results[,1] < serule,]

scd_possible[order(scd_possible[,4], decreasing = T), ]

scd_cv_results[cv_min_index,]


# actual estimation #
scd1 <- scd_cov_logr(X = X_train, y = y_train, blockcols = c(J/2, J/2), 
                     R = R, alpha = scd_alpha, lasso = c(45/2, 45, 45),
                     glasso = rep(2,3), ridge_y = scd_ridge, inits = "rational",
                     nrstart = 1, MAXITER = 5000, seed = 50,
                     stop_value = 1e-5)

scd1$result$W

scd1$result$Py

# in-sample and out-sample prediction (both 92 observations correctly classified)
scd_fit <- X_train %*% scd1$result$W %*% scd1$result$Py + scd1$result$py0

scd_fit <- sum(((1 / (1 + exp(-scd_fit))) > 0.5)  == y_train)

scd_pred <- X_test %*% scd1$result$W %*% scd1$result$Py + scd1$result$py0

scd_pred <- sum(((1 / (1 + exp(-scd_pred))) > 0.5)  == y_test)
  
  
  
  
  
# 5. lasso logistic regression ####

# we don't know the true parameters for logistic regression, because the method does not assume
# the covariate structure.
# Therefore, we will do cross validation

set.seed(123)
cv.lasso <- cv.glmnet(X_train, y_train, alpha = 1, family = "binomial")
plot(cv.lasso)

cv.lasso$lambda.min

cv.lasso$lambda.1se

coef(cv.lasso, cv.lasso$lambda.1se)

logistic <- glmnet(X_train, y_train, family = "binomial", alpha = 1, lambda =  cv.lasso$lambda.1se)

coef(logistic)

logistic_prob_test <- predict(object = logistic, newx = X_test)

sum((logistic_prob_test > 0.5) == y_test)

logistic_prob <- predict(object = logistic, newx = X_train)

sum((logistic_prob > 0.5) == y_train)

# 5. PCR (SCaDS - logR) #####

# model selection for lasso and group lasso; 
# Just need to look for penalties that return the correct weights structure

groups <- rep(15, 2)

ncomp <- 3

# starting values for the weights matrix
Wstart <- svd(X_train)$v[,1:3]

pcr <- newAlgoCpp(X_train, 
                  ridge = rep(0,3),
                  lasso = rep(45,3),
                  constraints = matrix(1, ncol(X_train), 3),
                  grouplasso = rep(10,3),
                  elitistlasso = rep(0,3),
                  groups = groups,
                  Q = ncomp,
                  itr = 100000,
                  nStarts = 1,
                  printLoss = FALSE,
                  coorDec = FALSE,
                  Wstart = Wstart)



pcr$W

pcr_T <- X_train %*% pcr$W

# logistic regression step # 
pcr_log <- glm(y_train ~ pcr_T, family=binomial(link='logit'))

coef(pcr_log)

pcr_fit <- predict(object = pcr_log, newx = X_train)

sum(abs(fitted(pcr_log) - (1 / (1 + exp(-pcr_fit)))))

y_train - fitted(pcr_log)

sum(y_train == (fitted(pcr_log) > 0.5))


pcr_log$coefficients

pcr_fit <- pcr_T %*% pcr_log$coefficients[-1] + pcr_log$coefficient[1]

# pcr_ber <- ber(pred = ((1 / (1 + exp(-pcr_fit))) > 0.5), true = y_train)

pcr_pred <- X_test %*%  pcr$W %*% pcr_log$coefficients[-1] + pcr_log$coefficients[1]

sum(((1 / (1 + exp(-pcr_pred))) > 0.5) == y_test)
  
  
# 6. DIABLO ####
  
# No need for cross validation because oracle information is known already

R <- 3
# first arranging the objects for the diablo function #
block1_3 <- data.frame(X_train[,1:(J/2)]) # first datablock
block2_3 <- data.frame(X_train[,(J/2 + 1):J]) # second datablock
block3_3 <- data.frame(X_train) # superblock


plsX_3 <- list(block1_3 = block1_3, 
               block2_3 = block2_3,
               block3_3 = block3_3)

# design matrix specification through PLS and observing pairwise correlation
# (this is as recommended by the DIABLO researchers)

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

# resulting design matrix is the null design (all values are at 0)

if (cond$dimension == "low"){
  list.keepX_3 = list(block1_3 = c(4), 
                      block2_3 = c(4),
                      block3_3 = c(8))
} else {
  list.keepX_3 = list(block1_3 = c(25), 
                      block2_3 = c(25),
                      block3_3 = c(40))
  
}

# actual model fitting - multiomics approach #
# always 3 components #
diathree1 <- block.splsda(X = plsX_3, Y = y_train,
                          ncomp = c(1), keepX = list.keepX_3, 
                          design = design_cor)

# fit #
diathree_fit <- predict(object = diathree1, newdata = plsX_3, dist = "max.dist")

diathree_fit$WeightedVote$max.dist[,1]

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

diathree_pred <- sum(y_test == diathree_pred$WeightedVote$max.dist[,1])


dia_d1 <- c(diathree1$loadings$block1_3, rep(0,15))

dia_d2 <- c(rep(0,15), diathree1$loadings$block2_3)

diathree_loadings <- cbind(dia_d1, dia_d2, diathree1$loadings$block3_3)

diathree1$weights



# save(reg.out, regsem_P, regsem_chosen, X_train, X_test, y_train, y_test, 
#      pcr, diathree1, 
#      diathree_loadings, scd_possible, alpha_ridge_cv_results, lasso_glasso_cv_results, scd1, 
#      plsX_3, plsXtest_3, diathree_fit, diathree_pred,
#      logistic, cv.lasso,
#      file = "../toy example/toy_results_signal050_seed11.Rdata")
  
  
  
  