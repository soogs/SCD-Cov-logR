# scdcovlogr_cv  #
# initiated: 16 march 2021 #
# last modified: 20 may 2021 #

# cross validation for scd-cov-logR#

# update 30 april 2021: #
# I allow for the multistart approach by adding the input "inits", "seed_cv" and "seed_method"

# update 20 may 2021: #
# I allow for a different CVE criterion: number of cases correctly classified in the test fold


scd_cv <- function(X, y, R, alpha, lasso, 
                       glasso, ridge_y, nrFolds, seed_cv, seed_method, inits, include_rational,
                       blockcols, MAXITER, stop_value, type_measure){
  
  # Input #
  # X : predictor matrix
  # y : binary criterion variable (only values 0 or 1)
  # blockcols : vector specifying the number of variables that each data block has
  # R : number of covariates
  # alpha : weighting parameter between X and Y
  # lasso : lasso
  # glasso : group lasso
  # ridge_y : ridge penalty on Y
  # inits : the initialization strategy used for model-fitting. "rational", "multistart"
  # seed_cv : seed used for the CV. the manner of folds being divided is dependent on this
  # seed_method : seed used for the model-fitting. This matters if you want to do multistart; the random starting values are dependent on this
  # include_rational : if the rational start is included as part of the multistart procedure
  # MAXITER: maximum number of iterations
  # type_measure: the type of measure used for cross validation error. "mse", "deviance" or "hit"
  
  set.seed(seed_cv)
  
  # randomly shuffle the data
  sampled <- sample(nrow(X))
  
  X <- X[sampled,]
  
  y <- y[sampled]
  
  #Create equally sized folds
  folds <- cut(seq(1,nrow(X)),breaks=nrFolds,labels=FALSE)
  
  cve_k <- data.frame(error_y = NA)
  
  corrects_k <- data.frame(corrects = NA)
  
  convergence <- c()
  
  for(k in 1:nrFolds){
    
    # defining test and train data  
    test_index <- which(folds==k, arr.ind=TRUE)
    
    X_train <- X[-test_index,]
    X_test <- X[test_index, ]
    
    y_train <- y[-test_index]
    y_test <- y[test_index]
    
    # model fitting #
    scd1 <- scd_cov_logr(X = X_train, y = y_train, blockcols = blockcols, 
                                     R = R, alpha = alpha, lasso = lasso,
                                     glasso = glasso, ridge_y = ridge_y, 
                                     inits = inits,
                                     include_rational = include_rational,
                                     nrstart = 1, MAXITER = MAXITER, seed = seed_method,
                                     stop_value = 1e-5)
    
    if (length(scd1$result$loss_hist) < MAXITER){
      converged <- T
    } else {
      converged <- F
    }
    
    # out of sample prediction #
    pred <- X_test %*% scd1$result$W %*% scd1$result$Py + scd1$result$py0
    
    if(type_measure == "mse"){
      
      ssq <- sum(((1 / (1 + exp(-pred))) - y_test)^2)
      
      mse <- mean(((1 / (1 + exp(-pred))) - y_test)^2)
      
      prob <- 1 / (1 + exp(-pred))
      
      mse2 <- mean((y_test - (prob))^2 + ((y_test != 1) - (1-prob))^2)
      
      cve_k[k,1] <- mse2
      
    } else if (type_measure == "deviance"){
      prob <- 1 / (1 + exp(-pred))
      
      prob <- pmin(pmax(prob, 1e-5), 1 - 1e-5)
      
      lp <- (y_test == 0) * log(1-prob) + y_test * log(prob)
      
      devv <- -2 * lp
      
      devv <- mean(devv)
      
      cve_k[k,1] <- devv
      
    } else if (type_measure == "hits"){
      prob <- 1 / (1 + exp(-pred))
      
      pred_case <- prob >= 0.5
      
      hits <- pred_case == y_test
      
      hit <- sum(hits)
      
      cve_k[k,1] <- hit
      
    } else {stop("no / wrong type measure specified")}
    
    
    pred_class <- (1 / (1 + exp(-pred))) > 0.5
    
    corrects <- sum(as.numeric(pred_class) == y_test)
    
    corrects_k[k,1] <- corrects
    
    convergence[k] <- converged
    
    
  }
  
  cve <- colMeans(cve_k)
  # cve = mean cv error
  
  corrects_mean <- colMeans(corrects_k)
  
  se <- apply(cve_k,2,sd) / sqrt(nrFolds)
  # se = standard deviation of the mean cv errors per k / sqrt(nrFold)
  
  corrects_se <- apply(corrects_k,2,sd) / sqrt(nrFolds)
  
  return(list(cve = cve, se = se, cve_k = cve_k, 
              corrects_mean = corrects_mean, 
              corrects_se = corrects_se,
              corrects_k = corrects_k,
              convergence = convergence))
}

