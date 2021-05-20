
# 500 Family data analysis #
# initiated: 3-may-2021 #
# last modified: 20-may-2021

# 20-may-2021 edit:
# I incorporate the leave-one-out CV for the family dataset, and use the number of correctly classified cases as CV error

# 0. Introduction ####
# The SCD-Cov-logR method is administered on the 500 family study dataset
# The dataset is publicly available at: https://www.icpsr.umich.edu/web/ICPSR/studies/4549?searchSource=revise&q=500+family+study#

# processing and preparing the dataset was done following the analysis done by Zhengguo Gu. 
# (Z Gu, K Van Deun - Behavior research methods, 2019)

# 1. data preparation #####
library(foreign)
library(psych)
library(RegularizedSCA)
library(nFactors)

Rcpp::sourceCpp("./functions/updateW_diffpen.cpp")

source("./functions/scd_cov_logR.R")

source("./functions/scd_cv.R")

DS3 = read.spss("E:\\Users\\park\\Desktop\\family data\\04549-0003-Data.sav", to.data.frame = TRUE)  # to be changed by user
DS4 = read.spss("E:\\Users\\park\\Desktop\\family data\\04549-0004-Data.sav", to.data.frame = TRUE) # to be changed by user

sim_impu <- function(DATA, M_values, repl){
  # M_values: missing values in the raw data. Some missing values are coded as 9, some are coded as 999 etc. 
  # repl = T, missing values are replaced with mean. 
  
  if(missing(repl)){
    repl = F
  }
  
  nrow <- dim(DATA)[1]
  ncol <- dim(DATA)[2]
  DATA <- data.matrix(DATA)
  
  for(i in 1:length(M_values)){
    DATA[which(DATA == M_values[i])] <- NA
  }
  
  if(repl == T){
    
    for(n in 1:nrow){
      for(m in 1:ncol){
        if(is.na(DATA[n, m])){
          DATA[n, m] <- mean(DATA[n, ], na.rm = T)
        }
      }
    }
  }
  
  return(DATA)
}

DS3[, c(352:367)] <- sim_impu(DS3[, c(352:367)], 9)
DS3[, c(353, 356, 359, 360, 363, 365)] <- 6 - DS3[, c(353, 356, 359, 360, 363, 365)] #reverse coding
DS3[, c(352:367)] <- sim_impu(DS3[, c(352:367)], 9, repl = T) #NA replaced with mean
summary(DS3[, c(352:367)])

DS3[, c(368:371)] <- sim_impu(DS3[, c(368:371)], 9)
DS3[, c(370, 371)] <- 6 - DS3[, c(370, 371)] #so the higher the value the calmer.
DS3[, c(368:371)] <- sim_impu(DS3[, c(368:371)], 9, repl = T)
summary(DS3[, c(368:371)])


DS3[, c(379:386)] <- sim_impu(DS3[, c(379:386)], c(9, -8), repl = T)
summary(DS3[, c(379:386)])

DS3[, c(387:404)] <- sim_impu(DS3[, c(387:404)], c(9, -8), repl = T)
summary(DS3[, c(387:404)])


DS3[, c(441:460)] <- sim_impu(DS3[, c(441:460)], c(9, -8), repl = T)
summary(DS3[, c(441:460)])

DS3[, c(540:543)] <- sim_impu(DS3[, c(540:543)], 9, repl = T)
DS3[, c(540:543)] <- 6 - DS3[, c(540:543)] #reverse coding
summary(DS3[, c(540:543)])

DS3[, c(544:549)] <- sim_impu(DS3[, c(544:549)], 9, repl = T)
summary(DS3[, c(544:549)])

DS3[, c(550:553)] <- sim_impu(DS3[, c(550:553)], 9)
DS3[, c(550, 551)] <- 5 - DS3[, c(550, 551)] #note, category 0 to 4
DS3[, c(550:553)] <- sim_impu(DS3[, c(550:553)], 9, repl = T)
DS3[, c(550:553)] <- 5 - DS3[, c(550:553)]  # reverse coding again 
summary(DS3[, c(550:553)])

DS3[, c(555:568)] <- sim_impu(DS3[, c(555:568)], 9)
DS3[, c(556:562, 565, 566)] <- 4 -  DS3[, c(556:562, 565, 566)]
DS3[, c(555:568)] <- sim_impu(DS3[, c(555:568)], 9, repl = T)
summary(DS3[, c(555:568)])


DS3_NEW <- rowMeans(DS3[, c(352:367)])
DS3_NEW <- cbind(DS3_NEW, rowMeans(DS3[, c(368:371)]))
DS3_NEW <- cbind(DS3_NEW, rowMeans(DS3[, c(379:386)]))
DS3_NEW <- cbind(DS3_NEW, rowMeans(DS3[, c(387:404)]))
DS3_NEW <- cbind(DS3_NEW, rowMeans(DS3[, c(441:460)]))
DS3_NEW <- cbind(DS3_NEW, rowMeans(DS3[, c(540:543)]))
DS3_NEW <- cbind(DS3_NEW, rowMeans(DS3[, c(544:549)]))
DS3_NEW <- cbind(DS3_NEW, rowMeans(DS3[, c(550:553)]))
DS3_NEW <- cbind(DS3_NEW, rowMeans(DS3[, c(555:568)]))
DS3_NEW <- cbind(DS3[, 1:3],DS3_NEW)
colnames(DS3_NEW)[4:12] <- c("Relationship with partners", "Argue with partners", "Childs bright future",
                             "Activities with children", "Family hassles", "Feeling about parenting", 
                             "Communation with children", "Argue with children", "Confidence about oneself")
summary(DS3_NEW)


DS4[, c(176:183)] <- sim_impu(DS4[, c(176:183)], 9)
DS4[, c(176, 181, 182)] <- 3 - DS4[, c(176, 181, 182)]  #recoding, answer category 0~3
DS4[, c(176:183)] <- sim_impu(DS4[, c(176:183)], 9, repl = T)
summary(DS4[, c(176:183)])



DS4[which(DS4[, 304] == 99), 304]<- NA #99:missing value
DS4[which(DS4[, 304] == 9), 304] <- NA #9: no grade

DS4[,304] <- sim_impu(as.data.frame(DS4[,304]), 9)

head(DS4[,304])

DS4[, 304] <- 9 - DS4[, 304]
summary(DS4[, 304])



DS4[, c(326:342)] <- sim_impu(DS4[, c(326:342)], 9, repl = T)
summary(DS4[, c(326:342)])


DS4[, c(359:367)] <- sim_impu(DS4[, c(359:367)], 9)
DS4[, c(361, 362)] <- 6 - DS4[, c(361, 362)]
DS4[, c(359:367)] <- sim_impu(DS4[, c(359:367)], 9, repl = T)
summary(DS4[, c(359:367)])


DS4[, c(381:397)] <- sim_impu(DS4[, c(381:397)], 9)
DS4[, c(383,384, 386:388, 390:394, 396,397)] <- 4 - DS4[, c(383,384, 386:388, 390:394, 396,397)] #reverse coding
DS4[, c(381:397)] <- sim_impu(DS4[, c(381:397)], 9, repl = T)
summary(DS4[, c(381:397)])

DS4[, c(399:418)] <- sim_impu(DS4[, c(399:418)], 9)
DS4[, c(399, 400, 402, 403, 405:411, 413, 414, 416:418)] <- 3 - DS4[, c(399, 400, 402, 403, 405:411, 413, 414, 416:418)]
DS4[, c(399:418)] <- sim_impu(DS4[, c(399:418)], 9, repl = T)
summary(DS4[, c(399:418)])


DS4[, c(468:480)] <- sim_impu(DS4[, c(468:480)], 9, repl = T)
summary(DS4[, c(468:480)])


DS4_NEW <- rowMeans(DS4[, 176:183]) #self confidence/esteem
DS4_NEW <- cbind(DS4_NEW, DS4[, 304]) #academic performance
DS4_NEW <- cbind(DS4_NEW, rowMeans(DS4[, 326:342])) # social/extracurricular
DS4_NEW <- cbind(DS4_NEW, rowMeans(DS4[, 359:367])) # importance of friendship
DS4_NEW <- cbind(DS4_NEW, rowMeans(DS4[, 381:397])) # self image
DS4_NEW <- cbind(DS4_NEW, rowMeans(DS4[, 399:418])) # happiness
DS4_NEW <- cbind(DS4_NEW, rowMeans(DS4[, 468:480])) # confidence about future
DS4_NEW <- cbind(DS4[, 1:3], DS4_NEW)
colnames(DS4_NEW)[4:10] <- c("Self confidence/esteem", "Academic performance", "Social life and extracurricular", "Importance of friendship",
                             "Self Image", "Happiness", "Confidence about the future")
summary(DS4_NEW)

DS3_NEW <- DS3_NEW[, -8] # removing family hassles

DS3_NEW_Final <- DS3_NEW[-sort(unique(which(is.na(DS3_NEW), arr.ind = T)[, 1])), ]
DS4_NEW_Final <- DS4_NEW[-sort(unique(which(is.na(DS4_NEW), arr.ind = T)[, 1])), ]
summary(DS3_NEW_Final)
summary(DS4_NEW_Final)

DS3_Mom <- DS3_NEW_Final[DS3_NEW_Final[, 2]=="mom", ]
DS3_Dad <- DS3_NEW_Final[DS3_NEW_Final[, 2]=="dad", ]
DS4_Kid <- DS4_NEW_Final
family_index <- intersect(intersect(DS3_Mom[, 3], DS3_Dad[, 3]), DS4_Kid[, 3])
DS3_Mom_Final <- DS3_Mom[DS3_Mom[, 3]==family_index[1], ]
for(i in 2:length(family_index)){
  DS3_Mom_Final <- rbind(DS3_Mom_Final, DS3_Mom[DS3_Mom[, 3]==family_index[i], ])
}
DS3_Dad_Final <- DS3_Dad[DS3_Dad[, 3]==family_index[1], ]
for(i in 2:length(family_index)){
  DS3_Dad_Final <- rbind(DS3_Dad_Final, DS3_Dad[DS3_Dad[, 3]==family_index[i], ])
}
DS4_Kid_Final <- DS4_Kid[DS4_Kid[, 3]==family_index[1], ]
for(i in 2:length(family_index)){
  DS4_Kid_Final <- rbind(DS4_Kid_Final, DS4_Kid[DS4_Kid[, 3]==family_index[i], ])
}
DS4_Kid_Final <- DS4_Kid_Final[!duplicated(DS4_Kid_Final[, 3]), ]  # a family can have more than 1 child. 

family_data <- list("mom" = DS3_Mom_Final[, -c(1:3)], "dad" = DS3_Dad_Final[, -c(1:3)], "child" = DS4_Kid_Final[, -c(1:3)])
colnames(family_data[[1]]) <- c("M: Relationship with partners", "M: Argue with partners", "M: Childs bright future",
                                "M: Activities with children", "M: Feeling about parenting", 
                                "M: Communation with children", "M: Argue with children", 
                                "M: Confidence about oneself")
colnames(family_data[[2]]) <- c("D: Relationship with partners", "D: Argue with partners", "D: Childs bright future",
                                "D: Activities with children", "D: Feeling about parenting", 
                                "D: Communation with children", "D: Argue with children", 
                                "D: Confidence about oneself")

dim(family_data)

dim(family_data$mom)


describe(family_data[[1]])  #mother
describe(family_data[[2]])  #father
describe(family_data[[3]])  #child

data <- cbind(family_data[[1]], family_data[[2]], family_data[[3]])

colnames(family_data[[3]])[2]

family_data[[3]][,2]

# academic achievement (most recent grade) as the outcome variable of interest
table(family_data[[3]][,2]) 
# 8 = most A
# 7 = half A half B
# 6 = most B
# 5 = half B half C

y <- as.numeric(family_data[[3]][,2] > 5)

# we will take below 6 as underachievement, and above 6 as overachievement

data<- cbind(pre_process(family_data[[1]]), pre_process(family_data[[2]]), pre_process(family_data[[3]][,-2]))
num_var <- cbind(dim(family_data[[1]])[2], dim(family_data[[2]])[2], dim(family_data[[3]])[2])

which(colnames(data) == "Academic performance")

y <- as.numeric(family_data[[3]][,2] > 5)

X <- cbind(pre_process(family_data[[1]]), pre_process(family_data[[2]]), pre_process(family_data[[3]][,-2]))

y1 <- which(y == 1)

# selecting a subset of observations to make the two groups even
set.seed(203)

selected <- sample(x = y1, size = sum(y == 0))

y0 <- which(y == 0)

selected <- c(selected, y0)

y <- y[selected]

X <- cbind(pre_process(family_data[[1]][selected,]),
           pre_process(family_data[[2]][selected,]),
           pre_process(family_data[[3]][selected,-2]))



# 2. Model selection for Rational starting values ####

# number of components
svdd <- svd(X)

nScree(svdd$d)

nScree(svdd$d^2)

plot(svdd$d^2)

plot(svdd$d)



colnames(X)

blockcols <- c(8, 8, 6)

blockcols2 <- cumsum(blockcols)

blockindex <- list()

blockindex[[1]] <- 1:blockcols2[1]

for (i in 2:length(blockcols)){
  blockindex[[i]] <- (blockcols2[i-1]+1):blockcols2[i]
}



# 2.1. cross validation for alpha and ridge (rational starting values)

# alpha and ridge CV #
alpha_range <- c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)
ridge_range <-  c(0.01, 0.05, 0.1, 0.5, 1, 2, 5, 10, 15, 20)
ridge_range <- sort(ridge_range, decreasing = TRUE) # bigger ridge is preferred in 1SE rule

alpha_ridge_ranges <- expand.grid(alpha_range, ridge_range)

colnames(alpha_ridge_ranges) <- c("alpha", "ridge")

pcs_cv_results <- matrix(NA, ncol = 7, nrow = nrow(alpha_ridge_ranges))

pcs_cv_results <- as.data.frame(pcs_cv_results)

colnames(pcs_cv_results) <- c("dev", "dev_se", "conv", "corrects", "corrects_se", 
                              "alpha", "ridge")

# generating random seed values for each combination of lasso and glasso
set.seed(7272)

pcs_seeds <- sample(1:100000, nrow(alpha_ridge_ranges))


for (i in 1:nrow(alpha_ridge_ranges)){
  
  scdcv1 <- scd_cv(X = X, y = y, R = 2, 
                       blockcols = blockcols,
                       alpha = alpha_ridge_ranges[i,]$alpha, 
                       lasso = rep(0,2), 
                       glasso = rep(0,2), 
                       ridge_y = alpha_ridge_ranges[i,]$ridge, seed_cv = pcs_seeds[i],
                       seed_method = 111, inits = "rational", include_rational = TRUE, 
                       MAXITER = 1000, nrFolds = 20, type_measure = "deviance")
  
  pcs_cv_results[i,] <- c(scdcv1$cve, scdcv1$se, 
                          sum(scdcv1$convergence),
                          scdcv1$corrects_mean,
                          scdcv1$corrects_se,
                          alpha_ridge_ranges[i,]$alpha, alpha_ridge_ranges[i,]$ridge)
  
  print(i)
}


pcs_cv_results

plot(pcs_cv_results[,1])

alpha_ridge_cv <- pcs_cv_results


# minimum error 
cv_min_index <- which.min(pcs_cv_results[,1])

# minimum error + 1 standard error
serule <- pcs_cv_results[cv_min_index,1] + pcs_cv_results[cv_min_index,2]

scd_possible <- pcs_cv_results[pcs_cv_results[,1] <= serule,]


scd_possible[order(scd_possible[,6]),]

alpha_chosen <- scd_possible[order(scd_possible[,6]),][1,]$alpha

ridge_chosen <-  scd_possible[order(scd_possible[,6]),][1,]$ridge

alpha_ridge_cve <-  scd_possible[order(scd_possible[,6]),][1,]$dev

alpha_ridge_se <-  scd_possible[order(scd_possible[,6]),][1,]$dev_se

alpha_ridge_conv <-  scd_possible[order(scd_possible[,6]),][1,]$conv



# 2.2. cross validation for lasso and glasso (rational starting values)

lasso_range <- c(10, 7, 5, 3, 1, 0.5, 0.3, 0.1, 0.05, 0)
glasso_range <- c(10, 7, 5, 3, 1, 0.5, 0.3, 0.1, 0.05, 0)

lasso_glasso_ranges <- expand.grid(lasso_range, glasso_range)

colnames(lasso_glasso_ranges) <- c("lasso", "glasso")

pcs_cv_results <- matrix(NA, ncol = 7, nrow = nrow(lasso_glasso_ranges))

pcs_cv_results <- as.data.frame(pcs_cv_results)

colnames(pcs_cv_results) <- c("dev", "dev_se", "conv", "corrects", "corrects_se", 
                              "lasso", "glasso")

set.seed(7272)

pcs_seeds <- sample(1:100000, nrow(lasso_glasso_ranges))

for (i in 1:nrow(lasso_glasso_ranges)){
  
  scdcv1 <- scd_cv(X = X, y = y, R = 2, 
                   blockcols = blockcols,
                   alpha = alpha_chosen, 
                   lasso = rep(lasso_glasso_ranges[i,]$lasso,2), 
                   glasso = rep(lasso_glasso_ranges[i,]$glasso,2), 
                   ridge_y = ridge_chosen,
                   MAXITER = 1000, nrFolds = 20, seed_cv = pcs_seeds[i], 
                   seed_method = 111, inits = "rational",
                   include_rational = TRUE, 
                   type_measure = "deviance")
  
  pcs_cv_results[i,] <- c(scdcv1$cve, scdcv1$se, 
                          sum(scdcv1$convergence),
                          scdcv1$corrects_mean,
                          scdcv1$corrects_se,
                          lasso_glasso_ranges[i,]$lasso, lasso_glasso_ranges[i,]$glasso)
  
  print(i)
  
}


lasso_glasso_cv <- pcs_cv_results

# minimum error 
cv_min_index <- which.min(pcs_cv_results[,1])

# minimum error + 1 standard error
serule <- pcs_cv_results[cv_min_index,1] + pcs_cv_results[cv_min_index,2]

scd_possible <- pcs_cv_results[pcs_cv_results[,1] < serule,]

scd_possible[order(scd_possible[,6], decreasing = T),]

lasso_chosen <- scd_possible[order(scd_possible[,6], decreasing = T),][1,]$lasso

glasso_chosen <- scd_possible[order(scd_possible[,6], decreasing = T),][1,]$glasso

lasso_glasso_cve <- scd_possible[order(scd_possible[,6], decreasing = T),][1,]$dev

lasso_glasso_se <- scd_possible[order(scd_possible[,6], decreasing = T),][1,]$dev_se

lasso_glasso_conv <-  scd_possible[order(scd_possible[,6], decreasing = T),][1,]$conv


# 2.3. final analysis with rational starting values 
scd_rational <- scd_cov_logr(X = X, y = y, blockcols = blockcols,
                                 R = 2, 
                                 alpha = alpha_chosen, 
                                 lasso = rep(lasso_chosen,2),
                                 glasso = rep(glasso_chosen,2),
                                 ridge_y = ridge_chosen, 
                                 inits = "rational", nrstart = 1, MAXITER = 10000,
                                 seed = 222, stop_value = 1e-5)


scd_rational$result$W
scd_rational$result$Py
scd_rational$result$py0

scd_rational$result$W

scd_rational$loss

pred <- 1 / (1 + exp(-(X %*% scd_rational$result$W %*% scd_rational$result$Py + scd_rational$result$py0)))

pred_class <- pred > 0.5 

insample_correct <- sum(pred_class ==  y)

insample_mse <- mean((y - (1 - pred))^2 + ((y != 1) - pred)^2)


# 3. Model selection for the multistart approach ####

# 3.1. alpha and ridge CV #
alpha_range <- c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)
ridge_range <-  c(0.01, 0.05, 0.1, 0.5, 1, 2, 5, 10, 15, 20)
ridge_range <- sort(ridge_range, decreasing = TRUE) # bigger ridge is preferred in 1SE rule

alpha_ridge_ranges <- expand.grid(alpha_range, ridge_range)

colnames(alpha_ridge_ranges) <- c("alpha", "ridge")

cv_multistart_alpha_ridge <- list()


for (nr in 1:50){
  
  family_cv_results <- matrix(NA, ncol = 7, nrow = nrow(alpha_ridge_ranges))
  
  
  for (i in 1:nrow(alpha_ridge_ranges)){
    
    
    scdcv1 <- scd_cv(X = X, y = y, R = 2, 
                     blockcols = blockcols,
                     alpha = alpha_ridge_ranges[i,]$alpha, 
                     lasso = rep(0,2), 
                     glasso = rep(0,2), 
                     ridge_y = alpha_ridge_ranges[i,]$ridge,
                     MAXITER = 1000, nrFolds = 20, seed_cv = i+1, seed_method = nr, 
                     inits = "multistart", include_rational = FALSE,
                     type_measure = "deviance")
    
    family_cv_results[i,] <- c(scdcv1$cve, scdcv1$se, 
                               sum(scdcv1$convergence),
                               scdcv1$corrects_mean,
                               scdcv1$corrects_se,
                               alpha_ridge_ranges[i,]$alpha, alpha_ridge_ranges[i,]$ridge)
    
  }
  
  colnames(family_cv_results) <- c("cve", "se", "conv", "corrects_mean", "corrects_se", "alpha", "ridge")
  
  cv_multistart_alpha_ridge[[nr]] <- family_cv_results
  
  print(nr)
}


# function to carry out the 1se rule <specific for family data> # 

se_select <- function(x){
  
  x <- as.data.frame(x)
  
  cv_min_index <- which.min(x$cve)
  
  serule <- x[cv_min_index,1] + x[cv_min_index,2]
  
  scd_possible <- x[x$cve < serule,]
  
  scd_possible <- scd_possible[order(scd_possible$alpha, decreasing = F),]
  
  if (is.vector(scd_possible)){
    scd_chosen <- scd_possible
  } else {
    scd_chosen <- scd_possible[1,]
  }
  
  return(scd_chosen)
}

alpha_ridge_cv_chosen <- lapply(cv_multistart_alpha_ridge, se_select)

time2 <- Sys.time()

# save(cv_multistart_alpha_ridge, file = "../family data/cv_multistart_alpha_ridge_again.Rdata")




# 3.2. cross validation for lasso and glasso

lasso_range <- c(10, 7, 5, 3, 1, 0.5, 0.3, 0.1, 0.05, 0)
glasso_range <- c(10, 7, 5, 3, 1, 0.5, 0.3, 0.1, 0.05, 0)

lasso_glasso_ranges <- expand.grid(lasso_range, glasso_range)

colnames(lasso_glasso_ranges) <- c("lasso", "glasso")


cv_multistart_lasso_glasso <- list()

for (nr in 1:50){
  
  family_cv_results <- matrix(NA, ncol = 7, nrow = nrow(lasso_glasso_ranges))
  
  
  for (i in 1:nrow(lasso_glasso_ranges)){
    
    
    scdcv1 <- scd_cv(X = X, y = y, R = 2,
                     blockcols = blockcols,
                     alpha = alpha_ridge_cv_chosen[[nr]]$alpha, 
                     lasso = rep(lasso_glasso_ranges[i,]$lasso,2), 
                     glasso = rep(lasso_glasso_ranges[i,]$glasso,2), 
                     ridge_y = alpha_ridge_cv_chosen[[nr]]$ridge,
                     MAXITER = 1000, nrFolds = 20, seed_cv = i+5, seed_method = nr, 
                     inits = "multistart", include_rational = FALSE,
                     type_measure = "deviance")

    
    family_cv_results[i,] <- c(scdcv1$cve, scdcv1$se, 
                               sum(scdcv1$convergence),
                               scdcv1$corrects_mean,
                               scdcv1$corrects_se,
                               lasso_glasso_ranges[i,]$lasso, lasso_glasso_ranges[i,]$glasso)
    
  }
  
  
  colnames(family_cv_results) <- c("cve", "se", "conv", "corrects_mean", "corrects_se", "lasso", "glasso")
  
  cv_multistart_lasso_glasso[[nr]] <- family_cv_results
  
  print(nr)
}


# save(cv_multistart_lasso_glasso, file = "../family data/cv_multistart_lasso_glasso_again.Rdata")

# function to carry out the 1se rule: lasso and glasso # 
se_select_lasso_glasso <- function(x){
  
  x <- as.data.frame(x)
  
  cv_min_index <- which.min(x$cve)
  
  serule <- x$cve[cv_min_index] + x$se[cv_min_index]
  
  scd_possible <- x[x$cve < serule,]
  
  scd_possible <- scd_possible[order(scd_possible$lasso, decreasing = T),]
  
  if (is.vector(scd_possible)){
    scd_chosen <- scd_possible
  } else {
    scd_chosen <- scd_possible[1,]
  }
  
  return(scd_chosen)
}



lasso_glasso_cv_chosen <- lapply(cv_multistart_lasso_glasso, se_select_lasso_glasso)




# 3.3. actual estimation: with the 50 multiple starts #

scd_multi <- list()

for (nr in 1:50){
  
  alpha_chosen <- alpha_ridge_cv_chosen[[nr]]$alpha
  
  ridge_chosen <- alpha_ridge_cv_chosen[[nr]]$ridge
  
  lasso_chosen <- lasso_glasso_cv_chosen[[nr]]$lasso
  
  glasso_chosen <- lasso_glasso_cv_chosen[[nr]]$glasso
  
  
  scd1 <- scd_cov_logr(X = X, y = y, blockcols = blockcols, 
                       R = 2, alpha = alpha_chosen, 
                       lasso = rep(lasso_chosen, 2),
                       glasso = rep(glasso_chosen,2), 
                       ridge_y = ridge_chosen, 
                       inits = "multistart", include_rational = FALSE,
                       nrstart = 1, MAXITER = 5000, seed = nr,
                       stop_value = 1e-5)
  scd_multi[[nr]] <- scd1
  
  print(nr)
}

# observing the objective attained by the methods
losses <- unlist(lapply(scd_multi, function(x){(x$loss)}))

which.min(losses)

lapply(scd_multi, function(x){(x$loss)})

scd_multi[[48]]$result$W

scd_multi[[48]]$result$Py

scd_multi[[48]]$result$py0
# smallest loss comes from rational start. so that is accepted, and reported in the paper.





# 4. interpreting the results from the rational starting values ####

rownames(scd_rational$result$W) <- colnames(X)

tmat <- as.data.frame(X %*% scd_rational$result$W)

tmat <- cbind(tmat, y)

tmat <- cbind(tmat, pred_class)

tmat <- as.data.frame(tmat)

tmat$y <- as.factor(tmat$y)

tmat$pred_class <- as.factor(tmat$pred_class)


data_long <- tidyr::gather(tmat, y, score, V1:V2, factor_key = TRUE)

library(ggplot2)

c_d1_plot <- ggplot(tmat, aes(x=V1, y=V2, color = y)) + 
  geom_point(size = 5) +
  theme_bw() +
  # xlim(c(-3, 3)) +
  # ylim(c(-4.5, 3)) +
  theme(legend.title = element_blank(),
        legend.position = "bottom",
        axis.title = element_text(size = 7),
        axis.text = element_text(size = 6),
        strip.text = element_text(size = 6),
        legend.text = element_text(size = 6),
        legend.box.margin=margin(-10,0,0,0))
# this plot is also reported in the paper



# 5. cv using this final model #####
scdcv1 <- scd_cv(X = X, y = y, R = 2, 
                     blockcols = blockcols,
                     alpha = alpha_chosen, 
                     lasso = rep(lasso_chosen,2), 
                     glasso = rep(glasso_chosen,2), 
                     ridge_y = ridge_chosen,
                     MAXITER = 10000, 
                     nrFolds = 20, seed_cv = 7272, seed_method = 111, inits = "rational", include_rational = TRUE, 
                     type_measure = "deviance")

scdcv1$cve # 1,29


# cv using hits as the criterion #
scdcv1 <- scd_cv(X = X, y = y, R = 2, 
                 blockcols = blockcols,
                 alpha = 0.2, 
                 lasso = rep(10,2), 
                 glasso = rep(10,2), 
                 ridge_y = 2,
                 MAXITER = 10000, 
                 nrFolds = 58, seed_cv = 7272, seed_method = 111, inits = "rational", include_rational = TRUE, 
                 type_measure = "hits")

scdcv1$cve # 0.7068
# 41 out-of-sample families are correctly classified 


# logistic regression cv (not reported in the paper) #
cv.logistic <- glmnet::cv.glmnet(x = X, y = y, family = "binomial", 
                                 type.measure = "deviance", nfolds = 20,
                                 alpha = 1)

plot(cv.logistic$cvm)

plot(cv.logistic)

cv.logistic$cvm[cv.logistic$index[2,]] # 1.35


cv.logistic$cvm[cv.logistic$index[1,]]

scdcv1$cve


library(glmnet)
# lambda.1se 
logistic <- glmnet(x = X, y = y, family = "binomial", alpha = 1, lambda = cv.logistic$lambda.1se)

coef(logistic)

sum((predict(object = logistic, newx = X, type = "response") > 0.5) == y)

# performance is therefore comparable to regularized logistic regression
