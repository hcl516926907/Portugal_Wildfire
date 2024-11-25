dir.data <- "/home/pgrad2/2448355h/My_PhD_Project/00_Dataset/Porgual_Wildfire"
dir.out <- "/home/pgrad2/2448355h/My_PhD_Project/01_Output/Portugal_Wildfire"

library(INLA)

inla.setOption(scale.model.default = TRUE)

library(mgcv)
library(ggplot2)
library(rgdal)
library(sf)
library(fmesher)
library(dplyr)
library(RColorBrewer)
library(terra)
library(rgeos)
library(spdep)
library(raster)
library(scoringRules)
library(fastDummies)
# bru_safe_sp(force = TRUE)

load(file.path(dir.data, "Data_For_Fitting_District.RData"))
data.fit <- data.fit.district
data.fit$NAME_1 <- as.character(data.fit$NAME_1)

data.fit$z <- as.vector((data.fit$y>0)+0)
data.fit$HVegTyp <- as.factor(data.fit$HVegTyp)
# levels(data.fit$HVegTyp) <-  c('NoCover','Evergreen_broadleaf_trees','Mixed_forest','Interrupted_forest')
data.fit$LVegTyp <- as.factor(data.fit$LVegTyp)
# levels(data.fit$LVegTyp) <- c('NoCover','Crops','Tall_grass','Semidesert','Evergreen_shrubs')
data.fit$month <- as.factor(data.fit$month)
data.fit <- dummy_cols(data.fit, 
                       select_columns = c('HVegTyp','LVegTyp','month','NAME_1'),remove_first_dummy=TRUE)

data.fit$year.af.2017.cnt = as.integer(data.fit$year>2017)
# covar.names <- colnames(data.fit)[c(-1,-2,-5,-6,-7,-8,-13,-16,-24)]
covar.names <- c('lon.grid','lat.grid','year',
                 'FWI','HVegCov','HVegLAI','LVegCov','LVegLAI',
                 'Pricp','RHumi','Temp','UComp','VComp','DewPoint','WindSpeed',
                 'LVegTyp_7','LVegTyp_16',
                 'month_2','month_3','month_4','month_5','month_6','month_7',
                 'month_8','month_9','month_10','month_11','month_12')
                 # "NAME_1_Beja",
                 # "NAME_1_Braga", "NAME_1_Bragança", "NAME_1_Castelo Branco",
                 # "NAME_1_Coimbra", "NAME_1_Évora", "NAME_1_Faro",
                 # "NAME_1_Guarda", "NAME_1_Leiria", "NAME_1_Lisboa",
                 # "NAME_1_Portalegre","NAME_1_Porto",
                 # "NAME_1_Santarém", "NAME_1_Setúbal", "NAME_1_Viana do Castelo",
                 # "NAME_1_Vila Real", "NAME_1_Viseu" )
                

data.rf <- data.fit[data.fit$year <= 2022, c(covar.names,'grid.idx','time.idx','log_ba','z','y')]
data.rf[is.na(data.rf$log_ba),'log_ba'] <- 0


library(xgboost)
library(rBayesianOptimization)
library(pROC)

set.seed(1234)

pos.loss <- mean(data.rf[data.rf$y>0,'log_ba'])
loss_eval <- function(pred,y,log_ba){
  pos.idx <- which(y>0)
  l1 <- sum((pred[pos.idx]*log_ba[pos.idx] - log_ba[pos.idx])^2)
  l2 <- sum((pred[-pos.idx]*pos.loss)^2)
  return((l1+l2)/length(pred))
}

xgb_z_cv   <- function(eta, max_depth, min_child_weight, subsample, colsample_bytree) {
  k <- 4
  loss_train_list <- numeric(k)
  loss_test_list <- numeric(k)

    data <- data.rf
  
  params <- list(
    objective = "binary:logistic",
    eval_metric = "auc",
    eta = eta,
    max_depth = as.integer(max_depth),
    min_child_weight = min_child_weight,
    subsample = subsample,
    colsample_bytree = colsample_bytree,
    scale_pos_weight = 1
  )
  
  for (i in 1:k) {
    select.year <- 2011:2013 + 3*(i-1)
    index <- which(data$year %in% select.year )
    train_data <- data[-index, covar.names ]
    train_target <- data[-index, 'z']
    test_data <- data[index, covar.names]
    test_target <- data[index, 'z' ]
    dtrain <- xgb.DMatrix(data = as.matrix(train_data), label = train_target)
    dtest <- xgb.DMatrix(data = as.matrix(test_data), label = test_target)
    
    set.seed(1234)
    model <- xgb.train(params = params, data = dtrain, nrounds = 100,
                       nthread = 30)
    
    train_pred <- predict(model, dtrain)
    test_pred <- predict(model, dtest)
  
    
    loss_train_list[i] <- loss_eval(train_pred,train_target,data[-index, 'log_ba'])
    loss_test_list[i] <-  loss_eval(test_pred,test_target,data[index, 'log_ba'])

  }
  
  mean_loss_train <- mean(loss_train_list)
  mean_loss_test <- mean(loss_test_list)
  
  
  return(list(Score=-mean_loss_test,Pred=NULL))
}

bounds_z <- list(
  eta = c(0.03, 0.1),
  max_depth = c(4L,10L),
  min_child_weight = c(1L,10L),
  subsample = c(0.5, 1.0),
  colsample_bytree = c(0.5, 1.0)
)

set.seed(1234)
t1 <- Sys.time()
xgb_z_params <- BayesianOptimization(
  FUN = xgb_z_cv,
  bounds = bounds_z,
  init_points = 50, # Number of initial random points
  n_iter = 30,      # Number of iterations for optimization
  acq = "ucb",      # Acquisition function, e.g., "ucb", "ei", "poi"
  kappa = 2.576,    # Parameter for the acquisition function
  eps = 0.0,        # Exploration-exploitation trade-off
  verbose = TRUE
)
t2 <- Sys.time()
print(t2-t1)

z_best_params <- xgb_z_params$Best_Par
print(z_best_params)
print( xgb_z_params$Best_Value)


tune_grid <- expand.grid(
  nrounds = seq(10,100,5),
  eta = z_best_params['eta'],
  max_depth =  z_best_params['max_depth'],
  gamma = 0,
  scale_pos_weight = 1,
  colsample_bytree = z_best_params['colsample_bytree'],
  min_child_weight = z_best_params['min_child_weight'],
  subsample = z_best_params['subsample']
)

results_df5 <- xgboost_tune_z(tune_grid)
print(results_df5[order(results_df5$mean_loss_test,decreasing=F),])

param_z <- list(
  nrounds = 90,
  eta = z_best_params['eta'],
  max_depth =  z_best_params['max_depth'],
  gamma = 0,
  scale_pos_weight = 1,
  colsample_bytree = z_best_params['colsample_bytree'],
  min_child_weight = z_best_params['min_child_weight'],
  subsample = z_best_params['subsample']
)

print(param_z)

#######################################
dextGP <- function(y,xi, sigma, kappa, log=FALSE){
  y1 <- 1+xi*y/sigma
  out.of.sup <- y1<=0 | y<=0
  y1[out.of.sup] <- NaN
  if (log){
    if (xi != 0){
      num <- log(kappa) + (kappa-1)*log(1-y1^(-1/xi))
      den <- log(sigma) + (1+1/xi)*log(y1)
    }else{
      num <- log(kappa) + (kappa-1)*log(1-exp(-y/sigma))
      den <- log(sigma) + y/sigma
    }
    res <- num-den
    res[which(out.of.sup)] <- -Inf
  }else{
    if(xi!=0){
      num <- kappa*(1-y1^(-1/xi))^(kappa-1)
      den <- sigma*y1^(1/xi +1)
    }else{
      num <- kappa*(1-exp(-y/sigma))^(kappa-1)
      den <- sigma*exp(y/sigma)
    }
    res <- num/den
    res[which(out.of.sup)] <- 0
  }
  return(res)
}

create_custom_objective <- function(xi, kappa, alpha) {
  function(preds, dtrain) {
    y <- getinfo(dtrain, "label")
    q_alpha <- exp(preds)
    
    # Compute sigma
    if (abs(xi) > 10^-6) {
      sigma <- xi * q_alpha / ( (1 - alpha^(1 / kappa))^(-xi) - 1 )
    } else {
      sigma <- -q_alpha / log(1 - alpha^(1 / kappa))
    }
    
    # Compute y1
    y1 <- 1 + xi * y / sigma
    
    # Handle invalid values
    invalid <- y1 <= 0 | y <= 0
    y1[invalid] <- NaN  # Avoid division by zero or log of non-positive numbers
    
    # Negative log-likelihood
    nll <- -dextGP(y, xi, sigma, kappa, log=TRUE)
    
    # Compute derivative of l w.r.t sigma
    if (abs(xi) > 10^-6) {
      # First derivative
      dl_dsigma <- (1 / sigma) - ((xi + 1) * y) / (sigma^2 * y1) + 
        ((kappa - 1) * y * y1^(-1 / xi - 1)) / (sigma^2 * (1 - y1^(-1 / xi)))
      
      # Second derivative (complex expression)
      # Compute components step by step
      term1 <- -1 / sigma^2
      
      term2 <- 2 * (xi + 1) * y / (sigma^3 * y1)
      term3 <- (xi + 1) * xi * y^2 / (sigma^4 * y1^2)
      
      term4_num <- 2 * y * y1^(-1 / xi - 1)
      term4_den <- sigma^3 * (1 - y1^(-1 / xi))
      term4 <- -(kappa - 1) * term4_num / term4_den
      
      term5_num <- (1 +  xi) * y^2 * y1^(-1 / xi - 2)
      term5_den <- sigma^4 * (1 - y1^(-1 / xi))
      term5 <- (kappa - 1) * term5_num / term5_den
      
      term6_num <-  y^2 * y1^(-2 / xi - 2)
      term6_den <- sigma^4 * (1 - y1^(-1 / xi))^2
      term6 <- (kappa - 1) * term6_num / term6_den
      
      d2l_dsigma2 <- term1 + term2 - term3 + term4 + term5 + term6
      
      # Compute gradient and hessian w.r.t pred

    } else {
      # For xi = 0, special case handling
      y2 <- exp(-y/sigma)
      dl_dsigma <- (1 / sigma) - y / sigma^2 + 
        ((kappa - 1) * y * y2) / (sigma^2 * (1 - y2))
      
      # Second derivative
      term1 <- -1 / sigma^2
      term2 <- 2 * y / sigma^3
      
      term3 <- (kappa - 1) * y2*y^2 / (sigma^4 * (1 - y2))
      
      term4_num <- 2*sigma*y2*(1-y2) - y2^2*y
      term4_den <- sigma^4*(1-y2)^2
      
      term4 <- (kappa-1)*y*term4_num/term4_den
      d2l_dsigma2 <- term1 + term2 + term3 - term4
      

    }
    
    grad <- sigma * dl_dsigma
    hess <- sigma^2 * d2l_dsigma2 + sigma * dl_dsigma
    
    # Handle NaNs or Infs
    grad[is.na(grad) | is.infinite(grad)] <- 0
    hess[is.na(hess) | is.infinite(hess)] <- 1  # To avoid division by zero in optimization
    
    return(list(grad = grad, hess = hess))
  }
}

loss_eval.egp <- function(preds, y, xi, kappa, alpha){
  q_alpha <- exp(preds)
  
  # Compute sigma based on xi
  sigma <- if (xi != 0) {
    xi * q_alpha / ( (1 - alpha^(1 / kappa))^(-xi) - 1 )
  } else {
    - q_alpha / log(1 - alpha^(1 / kappa))
  }
  
  return(-dextGP(y, xi, sigma, kappa, log=TRUE))
}




xgb_ba_cv   <- function(alpha, xi, kappa, eta, max_depth, min_child_weight, subsample, colsample_bytree) {
  k <- 4
  loss_train_list <- numeric(k)
  loss_test_list <- numeric(k)
  
  data <- data.rf
  custom_obj <- create_custom_objective(xi, kappa, alpha)
  
  
  params <- list(
    objective = custom_obj,
    eta = eta,
    max_depth = as.integer(max_depth),
    min_child_weight = min_child_weight,
    subsample = subsample,
    colsample_bytree = colsample_bytree
  )
  
  for (i in 1:k) {
    select.year <- 2011:2013 + 3*(i-1)
    index <- which(data$year %in% select.year )
    train_data <- data[-index, covar.names ]
    train_target <- data[-index, 'log_ba']
    pos.ind.train <- train_target>0
    train_data <- train_data[pos.ind.train,]
    train_target <- train_target[pos.ind.train]
    # print(table(data[-index,'grid.idx']))
    
    test_data <- data[index, covar.names]
    test_target <- data[index, 'log_ba' ]
    pos.ind.test <- test_target>0
    test_data <- test_data[pos.ind.test,]
    test_target <- test_target[pos.ind.test]
    
    dtrain <- xgb.DMatrix(data = as.matrix(train_data), label = train_target)
    dtest <- xgb.DMatrix(data = as.matrix(test_data), label = test_target)
    
    set.seed(1234)
    model <- xgb.train(params = params, data = dtrain, nrounds = 100,
                       nthread = 30)
    
    train_pred <- predict(model, dtrain)
    test_pred <- predict(model, dtest)
    
    loss_train_list[i] <- mean(loss_eval.egp(train_pred,train_target,xi,kappa, alpha))
    loss_test_list[i] <-  mean(loss_eval.egp(test_pred,test_target,xi,kappa, alpha))
    
  }
  
  mean_loss_train <- mean(loss_train_list)
  mean_loss_test <- mean(loss_test_list)
  
  if (is.infinite(mean_loss_test)){
    mean_loss_test <- 10
  }
  return(list(Score=-mean_loss_test,Pred=NULL))
}


bounds_ba <- list(
  alpha = c(0.2,0.8),
  xi = c(-0.05,0.3),
  kappa = c(0.5,5),
  eta = c(0.02, 0.1),
  max_depth = c(2L,6L),
  min_child_weight = c(1L,15L),
  subsample = c(0.5, 1.0),
  colsample_bytree = c(0.5, 1.0)
)

set.seed(1234)
t1 <- Sys.time()
xgb_ba_params <- BayesianOptimization(
  FUN = xgb_ba_cv,
  bounds = bounds_ba,
  init_points = 100, # Number of initial random points
  n_iter = 30,      # Number of iterations for optimization
  acq = "ucb",      # Acquisition function, e.g., "ucb", "ei", "poi"
  kappa = 2.576,    # Parameter for the acquisition function
  eps = 0.0,        # Exploration-exploitation trade-off
  verbose = TRUE
)
t2 <- Sys.time()
print(t2-t1)

ba_best_params <- xgb_ba_params$Best_Par
print(ba_best_params)
print( xgb_ba_params$Best_Value)



xi <- -0.0191887
kappa <- 2.72282415
alpha <- 0.5

custom_obj <- create_custom_objective(xi, kappa, alpha)




params.ba <- list(
  nrounds = 90,
  eta = ba_best_params['eta'],
  max_depth =  ba_best_params['max_depth'],
  gamma = 0,
  colsample_bytree = ba_best_params['colsample_bytree'],
  min_child_weight = ba_best_params['min_child_weight'],
  subsample = ba_best_params['subsample'],
  objective=custom_obj
)

print(params.ba)



i <- 1
data <- data.rf
select.year <- 2011:2013 + 3*(i-1)
index <- which(data$year %in% select.year )
train_data <- data[-index, covar.names ]
train_target <- data[-index, 'log_ba']
pos.ind.train <- train_target>0
train_data <- train_data[pos.ind.train,]
train_target <- train_target[pos.ind.train]
# print(table(data[-index,'grid.idx']))

test_data <- data[index, covar.names]
test_target <- data[index, 'log_ba' ]
pos.ind.test <- test_target>0
test_data <- test_data[pos.ind.test,]
test_target <- test_target[pos.ind.test]


dtrain <- xgb.DMatrix(data = as.matrix(train_data), label = train_target)
dtest <- xgb.DMatrix(data = as.matrix(test_data), label = test_target)

set.seed(1234)
model <- xgb.train(params = params, data = dtrain, nrounds = params$nrounds,nthread=30)

train_pred <- predict(model, dtrain)
test_pred <- predict(model, dtest)

summary(exp(train_pred))

loss_train <- mean(loss_eval.egp(train_pred,train_target,xi,kappa, alpha))
loss_test <-  mean(loss_eval.egp(test_pred,test_target,xi,kappa, alpha))
print(loss_test)


