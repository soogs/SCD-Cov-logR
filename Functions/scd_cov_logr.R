# SCD-Cov-logR #
# initiated: 19 May 2020
# last modified: 30 April 2021

# update 30 april 2021:
# making the function work with multistart approach 

scd_cov_logr <- function(X, y, blockcols, R, alpha, 
                             lasso, glasso, ridge_y, 
                             inits = c("rational", "oracle", "multistart"), 
                             nrstart = 10, 
                             include_rational = TRUE,
                             MAXITER = 10000,
                             seed,
                             stop_value = 1e-10){
  
  # Input #
  # X : predictor matrix
  # y : binary criterion variable (only values 0 or 1)
  # blockcols : vector specifying the number of variables that each data block has
  # R : number of covariates
  # alpha : weighting parameter between X and Y
  # lasso : lasso
  # glasso : group lasso
  # ridge_y : ridge penalty on Y
  # inits : starting value specification
  # nrstart : number of multiple starts
  # include_rational : whether rational start is included as a set of the multistart
  # MAXITER : Maximum number of iteration
  # stop_value : tolerance
  
  
  # blockindex definition #
  blockcols2 <- cumsum(blockcols)
  
  blockindex <- list()
  
  blockindex[[1]] <- 1:blockcols2[1]
  
  for (i in 2:length(blockcols)){
    blockindex[[i]] <- (blockcols2[i-1]+1):blockcols2[i]
  }
  
  # 1. define sub-functions ####
  
  updatePx <- function(wX, W, X){
    
    if (sum(colSums(W != 0) == 0) > 0){
      # print ("ERROR: W matrix has zero-columns. This disallows the P computation.")
      return (NA)
    }
    
    Tmat <- X %*% W
    K1 <- t(wX) %*% Tmat
    K <- t(Tmat) %*% wX %*% t(wX) %*% Tmat
    eigs <- eigen(K)
    V <- eigs$vectors
    Ssq <- eigs$values
    
    S <- diag(Ssq^(-0.5))
    
    if(ncol(W) == 1){
      S <- matrix(Ssq^(-0.5), ncol = 1, nrow = 1)
    }
    
    Px <- K1 %*% V %*% S %*% t(V)
    
    return (Px)
  }
  
  # updatePy #
  updatePy <- function(X, W, qq, z, ridge_y, py0, alpha){
    
    n <- nrow(X)
    
    R <- ncol(W)
    
    left <- MASS::ginv(t(X %*% W) %*% diag(c(qq)) %*% (X %*% W) + diag((2 / alpha) * ridge_y, R))
    
    right <- t(X %*% W) %*% diag(c(qq)) %*% z - (py0 * t(X %*% W) %*% qq)
    
    result <- left %*% right
    
    return (result)
  }
  
  # updatePy0 #
  updatePy0 <- function(X, W, qq, z, Py, alpha){
    
    result <- (1/ sum(qq)) * t(qq) %*% (z - X %*% W %*% Py)
    
    
    return (result)
    
    # i confirm that the above vector algebra is the same as the element-wise calculation:
    # hey <- c()
    # for(i in 1:n){
    #   hey[i] <- qq[i] * (z[i] - t(X[i,] %*% W %*% Py))
    # }
    # 
    # sum(hey) * (1 / (n * mean(qq)))
  }
  
  # loss calculation function #
  
  # the loss is made up of:
  # quadratic approximation of the negative log-likelihood of logistic regression
  # PCA loss and the penalties
  
  # first calculate q and z needed for the quadratic approximation
  quad_ingredients <- function(X, y, py0, Py, W){
    
    p <- 1 / (1 + exp(-(py0 + t(Py) %x% X %*% c(W))))
    
    qq <- p * (1 - p)
    
    # to avoid diverging problem (section 3 in Friedman paper) # 
    near1 <- abs(p - 1) < 1e-5
    
    p[near1] <- 1
    qq[near1] <- 1e-5
    
    near0 <- abs(p - 0) < 1e-5
    
    p[near0] <- 1e-5
    qq[near0] <- 1e-5
    
    z <- (py0 + t(Py) %x% X %*% c(W)) + (y - p) / qq
    
    result <- list(qq = qq, z = z)
    
    return(result)
  }
  
  
  losscal <- function(y, X, W, Py, py0, Px, qq, z, alpha, lasso, glasso, ridge_y, blockindex){
    
    # neg_lQ = negative log likelihood quadratic approximation
    neg_lQ <- (1 / 2) * t(qq) %*% (z - py0 - X %*% W %*% Py)^2
    
    # i confirm that neg_lQ is the same as element-wise calculation:
    # yo <- c()
    # 
    # for(i in 1:5){
    #   yo[i] <- qq[i] * (z[i] - py0 - t(Py) %x% t(X[i,]) %*% c(Wold))^2
    # }
    
    pca_loss <- sum((X - X %*% W %*% t(Px))^2) 
    
    # figuring out the penalties
    
    glasso_norm <- function(x, blockindex, R){
      l2norm <- 0
      
      for (r in 1:R){
        for (i in 1:length(blockindex)){
          l2norm <- l2norm + glasso[r] * sqrt(sum(x[blockindex[[i]],r]^2)) * sqrt(length(blockindex[[i]]))
        }
      }
      return (l2norm)
    }
    
    l2norm <- glasso_norm(x = W, blockindex = blockindex, R = ncol(W))
    
    result <- alpha * as.numeric(neg_lQ) + (1 - alpha) * pca_loss + 
      (l2norm) + as.numeric(t(lasso) %*% colSums(abs(W))) + (ridge_y * sum(Py^2))
    
    return (result)
  }
  
  # pcovr function from the bmc bioinformatics paper (Van Deun, Crompvoets and Ceulemans, 2018)
  pcovr <- function(X, Y, R, alpha){
    
    # if Y is provided as a vector.. #
    if (is.vector(Y)){
      Y <- matrix(data = Y, ncol = 1)
    }
    
    I <- nrow (X)
    Jx <- ncol (X)
    Jy <- ncol (Y)
    J <- Jx + Jy
    eps <- 1e-12
    
    # [Iy,Jy]=size(Y);
    # if Iy~=I, disp(' size Y and X do not match ');return;end;
    
    # weighting X and Y according to the alpha parameter
    w1 <- (I*J*alpha) / (sum(Y^2))
    w2 <- (I*J*(1-alpha)) / (sum(X^2))
    
    wY <- sqrt(w1) * Y
    wX <- sqrt(w2) * X
    
    # Z = [wY wX]
    Xstar <- cbind(wY, wX)
    
    # SVD data #
    if (J > I){
      XX <- X %*% t(X)
      eigs <- eigen(XX)
      Ux <- eigs$vectors
      Ssq <- eigs$values
      Sxsq <-Ssq
      Sxsq[Sxsq < 1e-6] <- 1e-6
      Sx <- sqrt(Sxsq)
      
      invSx <- Sx^(-1)
      
      Vx <- t(X) %*% Ux %*% diag(invSx)
      
      Uxstar <- t(Ux) %*% Xstar
      
      svd1 <- svd(Uxstar)
      U <- svd1$u[,1:R]
      S <- diag(svd1$d[1:R])
      V <- svd1$v[,1:R]
      
      
      if (R == 1){
        U <- matrix(svd1$u[,1:R], ncol = 1)
        S <- matrix(svd1$d[1:R], nrow = 1, ncol = 1)
        V <- matrix(svd1$v[,1:R], ncol = 1)
      }
      
    } else if (I >= J){
      XX <- t(X) %*% X
      eigs <- eigen(XX)
      
      Vx <- eigs$vectors
      Ssq <- eigs$values
      Sxsq <- Ssq
      Sxsq[Sxsq < 1e-6] <- 1e-6
      
      Sx <- sqrt(Sxsq)
      
      invSx <- Sx^(-1)
      
      Ux <- X %*% Vx %*% diag(invSx)
      Uxstar <- t(Ux) %*% Xstar
      
      svd1 <- svd(Uxstar)
      U <- svd1$u[,1:R]
      S <- diag(svd1$d[1:R])
      V <- svd1$v[,1:R]
      
      if (R == 1){
        U <- matrix(svd1$u[,1:R], ncol = 1)
        S <- matrix(svd1$d[1:R], nrow = 1, ncol = 1)
        V <- matrix(svd1$v[,1:R], ncol = 1)
      }
    }
    
    P <- V
    
    W <- Vx %*% diag(invSx) %*% U %*% S
    
    # reweighting step #
    Py <- P[1:Jy,] / (sqrt(w1))
    Px <- P[-c(1:Jy),] / (sqrt(w2))
    
    # if Py turns out to be a vector,
    # we need to make it into a matrix of one row 
    if (is.vector(Py)){
      Py <- matrix(data = Py, ncol = 1)
    }
    
    # fit measure #
    RsqX <- 1-sum(sum((X - X %*% W %*% t(Px))^2))/(sum(sum(X^2)))
    Rsqy <- 1-sum(sum((Y - X %*% W %*% t(Py))^2))/(sum(sum(Y^2)))
    
    return_list <- list(W = W, Px = Px, Py = Py, RsqX = RsqX, Rsqy = Rsqy)
    
    return (return_list)
  }
  
  # 2. define a few objects ####
  
  # blockindex definition #
  blockcols2 <- cumsum(blockcols)
  
  blockindex <- list()
  
  blockindex[[1]] <- 1:blockcols2[1]
  
  for (i in 2:length(blockcols)){
    blockindex[[i]] <- (blockcols2[i-1]+1):blockcols2[i]
  }
  
  # 3. initial values generated ####
  
  # initial value #
  # for both W and P, we have to provide initial values # 
  # borrowing the initial value system from the spca_adj function # 
  
  if (inits == "rational"){
    # W and P from pcovr
    
    y2 <- y == 0
    
    y_wide <- cbind(y, y2)
    
    y_wide <- scale(y_wide, T , T)
    
    pcovr_results <- pcovr(X = X, Y = y_wide, R = R, alpha = alpha)
    
    W <- pcovr_results$W
    # Px <- pcovr_results$Px
    # Py <- pcovr_results$Py
    # P <- rbind(Py, Px)
    
    nrstart <- 1
  } 
  
  if (inits == "oracle"){
    # user-specified matrices for W and P
    nrstart <- 1
  } 
  
  # vector to save results later
  multi_results <- list()
  multi_loss <- c()
  
  for (nr in 1:nrstart){
    
    if (inits == "multistart"){
      # initial values W, Px (orthogonal), Py
      
      if (nr == nrstart & include_rational){
        # the last set of starting values is always rational start
        
        y2 <- y == 0
        
        y_wide <- cbind(y, y2)
        
        y_wide <- scale(y_wide, T , T)
        
        pcovr_results <- pcovr(X = X, Y = y_wide, R = R, alpha = alpha)
        
        W <- pcovr_results$W
        
        svdX <- svd(X)
        Px <- svdX$v[,1:R]
        
        if (R == 1){
          Px <- matrix(svdX$v[,1:R], ncol = 1)
        }
        
        set.seed(seed)
        
        logregs <- glm(formula = y ~ X %*% W, family = "binomial")$coefficients
        
        Py <- matrix(logregs[-1], ncol = 1)
        
        py0 <- as.numeric(logregs[1])
        
      } else {
        # for the other sets of starting values, initial values are random
        
        set.seed(seed+nr-1)
        
        W <- matrix(stats::runif(n = ncol(X)*R, min = -1, max = 1), nrow = ncol(X), ncol = R)
        
        Px <- matrix(stats::runif(n = ncol(X)*R, min = -1, max = 1), nrow = ncol(X), ncol = R)
        Px <- qr.Q(qr(Px))
        
        Py <- matrix(stats::runif(n = R, min = -1, max = 1), nrow = R, ncol = 1)
        
        py0 <- runif(1, min = -1, max = 1)
      }
      
    }
    
    
    if (inits == "rational"){
      # initial values for P from SVD #
      svdX <- svd(X)
      Px <- svdX$v[,1:R]
      
      if (R == 1){
        Px <- matrix(svdX$v[,1:R], ncol = 1)
      }
      
      set.seed(seed)
      
      logregs <- glm(formula = y ~ X %*% W, family = "binomial")$coefficients
      
      Py <- matrix(logregs[-1], ncol = 1)
      
      py0 <- as.numeric(logregs[1])
      
    }
    
    # initial loss defined #
    loss0 <- 10000000
    
    Wold <- W + 100
    
    # convergence starting 
    conv <- 0
    iter <- 1
    
    loss_hist <- c(loss0)
    
    # 5. estimation ####
    while (conv == 0){
      
      quads <- quad_ingredients(X = X, y = y, py0 = py0, Py = Py, W = W)
      
      qq <- quads$qq
      
      z <- quads$z
      
      # after update of the quadratic approximation, new loss is calculated #
      loss0 <- losscal(y = y, X = X, W = W, Py = Py, py0 = py0, Px = Px,
                       qq = qq, z = z, alpha = alpha, lasso = lasso,
                       glasso = glasso, ridge_y = ridge_y, blockindex = blockindex)
      
      
      # W given P #
      
      Wnew <- checkit_reference(X = X, W = W, Px = Px, Py = Py, 
                                py0 = py0, lasso = lasso, glasso = glasso, ridge = ridge_y,
                                blockindex = blockindex, R = R, alpha = alpha, z = z, qq = qq)$Wnew
      
      if(anyNA(Wnew)){
        print("Wnew computation failed")
      }
      
      # print ("weights estimation complete")
      
      loss1 <- losscal(y = y, X = X, W = Wnew, Py = Py, py0 = py0, Px = Px,
                       qq = qq, z = z, alpha = alpha, lasso = lasso,
                       glasso = glasso, ridge_y = ridge_y, blockindex = blockindex)
      
      # loss should be non-increasing #
      if (loss0 < loss1){
        print ("ERROR: current loss > previous loss, after weights estimation")
      }
      
      # check convergence via loss #
      if (abs(loss0 - loss1) < stop_value){
        conv <- 1
      }
      
      # monitoring the convergence of the weights matrix #
      if(sum(abs(Wnew - W)) < 1e-5){
        conv <- 1
      }
      
      W <- Wnew
      
      iter <- iter + 1
      
      loss_hist[iter] <- loss1
      
      loss0 <- loss1
      
      # new calculation quads #
      quads <- quad_ingredients(X = X, y = y, py0 = py0, Py = Py, W = W)
      
      qq <- quads$qq
      
      z <- quads$z
      
      # after update of the quadratic approximation, new loss is calculated #
      loss0 <- losscal(y = y, X = X, W = W, Py = Py, py0 = py0, Px = Px,
                       qq = qq, z = z, alpha = alpha, lasso = lasso,
                       glasso = glasso, ridge_y = ridge_y, blockindex = blockindex)
      
      
      # Px given W #
      Px <- updatePx(wX = X, X = X, W = W)
      # in the SCD-CovDA version or ordinary SCD-CovR, 
      # wX = weighted X is given
      # but in the current function, i do not weight the X
      
      # if Px could not be calculated #
      if (anyNA(Px)){
        
        svdP <- svd(t(W) %*% t(X) %*% X)
        
        Px <- svdP$v %*% t(svdP$u)
        
      }
      
      # loss checking #
      loss1 <- losscal(y = y, X = X, W = W, Py = Py, py0 = py0, Px = Px,
                       qq = qq, z = z, alpha = alpha, lasso = lasso,
                       glasso = glasso, ridge_y = ridge_y, blockindex = blockindex)
      
      
      # loss should be non-incresing #
      if (loss0 < loss1){
        print ("ERROR: current loss > previous loss, after Px estimation")
      }
      
      # convergence checking via loss #
      if (abs(loss0 - loss1) < stop_value){
        conv <- 1
      }
      
      loss0 <- loss1
      
      # new calculation quads #
      quads <- quad_ingredients(X = X, y = y, py0 = py0, Py = Py, W = W)
      
      qq <- quads$qq
      
      z <- quads$z
      
      # after update of the quadratic approximation, new loss is calculated #
      loss0 <- losscal(y = y, X = X, W = W, Py = Py, py0 = py0, Px = Px,
                       qq = qq, z = z, alpha = alpha, lasso = lasso,
                       glasso = glasso, ridge_y = ridge_y, blockindex = blockindex)
      
      
      if (alpha != 0){
        # Py update #
        Py <- updatePy(X = X, W = W, qq = qq, z = z, ridge_y = ridge_y, py0 = py0, alpha = alpha)
        
        
        # loss check #
        loss1 <- losscal(y = y, X = X, W = W, Py = Py, py0 = py0, Px = Px,
                         qq = qq, z = z, alpha = alpha, lasso = lasso,
                         glasso = glasso, ridge_y = ridge_y, blockindex = blockindex)
        
        
        # loss should be non-incresing #
        if (loss0 < loss1){
          print ("ERROR: current loss > previous loss, after Py estimation")
        }
        
        # convergence checking via loss #
        if (abs(loss0 - loss1) < stop_value){
          conv <- 1
        }
        
        
        # new calculation quads #
        quads <- quad_ingredients(X = X, y = y, py0 = py0, Py = Py, W = W)
        
        qq <- quads$qq
        
        z <- quads$z
        
        # after update of the quadratic approximation, new loss is calculated #
        loss0 <- losscal(y = y, X = X, W = W, Py = Py, py0 = py0, Px = Px,
                         qq = qq, z = z, alpha = alpha, lasso = lasso,
                         glasso = glasso, ridge_y = ridge_y, blockindex = blockindex)
        
        # py0 update #
        py0 <- as.numeric(updatePy0(X = X, W = W, qq = qq, z = z, Py = Py, alpha = alpha))
        
        loss1 <- losscal(y = y, X = X, W = W, Py = Py, py0 = py0, Px = Px,
                         qq = qq, z = z, alpha = alpha, lasso = lasso,
                         glasso = glasso, ridge_y = ridge_y, blockindex = blockindex)
        
      }
      
      # loss should be non-incresing #
      if (loss0 < loss1){
        print ("ERROR: current loss > previous loss, after Py0 estimation")
      }
      
      iter <- iter + 1
      
      loss_hist[iter] <- loss1
      
      loss0 <- loss1
      
      # if MAXITER is reached #
      if (iter > MAXITER){
        print ("max iteration reached")
        conv <- 1
      }
      
    }
    
    # reweighting step #
    # Py <- matrix(P[1:Jy,] / (sqrt(w1)), ncol = R)
    # Px <- P[-c(1:Jy),] / (sqrt(w2))
    # 
    # if (R == 1){
    #   Px <- matrix(P[-c(1:Jy),] / (sqrt(w2)), ncol = 1)
    # }
    # 
    # P <- rbind(Py, Px)
    
    result_list <- list(W = W, Py = Py, Px = Px, py0 = py0, 
                        loss = loss1, loss_hist = loss_hist, 
                        iter = iter)
    
    multi_loss[nr] <- loss1
    multi_results[[nr]] <- result_list
    
  }
  
  lossmin <- which.min(multi_loss)
  multi_results_min <- multi_results[[lossmin]]
  
  # browser()
  
  return_results <- list(loss = multi_loss, result = multi_results_min, lossmin = lossmin)
  return (return_results)
}

print("please source updateW.cpp function")
