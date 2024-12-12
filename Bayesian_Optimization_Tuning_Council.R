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

load(file.path(dir.data, "Data_For_Fitting_Council.RData"))
data.fit <- data.fit.council
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
                 'HVegTyp_6','HVegTyp_18','HVegTyp_19',
                 'LVegTyp_1','LVegTyp_7','LVegTyp_11','LVegTyp_16',
                 'month_2','month_3','month_4','month_5','month_6','month_7',
                 'month_8','month_9','month_10','month_11','month_12',
                 "NAME_1_Beja",
                 "NAME_1_Braga", "NAME_1_Bragança", "NAME_1_Castelo Branco",
                 "NAME_1_Coimbra", "NAME_1_Évora", "NAME_1_Faro",
                 "NAME_1_Guarda", "NAME_1_Leiria", "NAME_1_Lisboa",
                 "NAME_1_Portalegre","NAME_1_Porto",
                 "NAME_1_Santarém", "NAME_1_Setúbal", "NAME_1_Viana do Castelo",
                 "NAME_1_Vila Real", "NAME_1_Viseu"
)



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
    scale_pos_weight = 2
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

xgboost_tune_z <- function(tune_grid){
  
  k <- 4
  
  results_df <- data.frame(
    eta = numeric(0),
    max_depth = numeric(0),
    gamma = numeric(0),
    colsample_bytree = numeric(0),
    min_child_weight = numeric(0),
    subsample = numeric(0),
    nrounds = numeric(0),
    mean_auc_train = numeric(0),
    mean_auc_test = numeric(0)
  )
  for (i in 1:nrow(tune_grid)) {
    print(i)
    params <- list(
      eta = tune_grid$eta[i],
      max_depth = tune_grid$max_depth[i],
      gamma = tune_grid$gamma[i],
      colsample_bytree = tune_grid$colsample_bytree[i],
      min_child_weight = tune_grid$min_child_weight[i],
      subsample = tune_grid$subsample[i],
      nrounds = tune_grid$nrounds[i],
      objective = "binary:logistic",
      eval_metric = "auc",
      scale_pos_weight = tune_grid$scale_pos_weight[i],
      verbosity = 0
    )
    
    result <- kfold_cv(data.rf, 'z', covar.names, params, k)
    
    # Store results in the dataframe
    results_df <- rbind(results_df, data.frame(
      eta = tune_grid$eta[i],
      max_depth = tune_grid$max_depth[i],
      gamma = tune_grid$gamma[i],
      colsample_bytree = tune_grid$colsample_bytree[i],
      min_child_weight = tune_grid$min_child_weight[i],
      subsample = tune_grid$subsample[i],
      nrounds = tune_grid$nrounds[i],
      scale_pos_weight = tune_grid$scale_pos_weight[i],
      mean_auc_train = result$mean_auc_train,
      mean_auc_test = result$mean_auc_test,
      mean_loss_train = result$mean_loss_train,
      mean_loss_test = result$mean_loss_test
    ))
  }
  
  return(results_df)
}
# 

kfold_cv  <- function(data, target.name ,covar.names, params, k) {
  auc_train_list <- numeric(k)
  auc_test_list <- numeric(k)
  
  loss_train_list <- numeric(k)
  loss_test_list <- numeric(k)
  
  
  
  for (i in 1:k) {
    select.year <- 2011:2013 + 3*(i-1)
    index <- which(data$year %in% select.year )
    train_data <- data[-index, covar.names ]
    train_target <- data[-index, target.name]
    test_data <- data[index, covar.names]
    test_target <- data[index, target.name ]
    dtrain <- xgb.DMatrix(data = as.matrix(train_data), label = train_target)
    dtest <- xgb.DMatrix(data = as.matrix(test_data), label = test_target)
    
    set.seed(1234)
    model <- xgb.train(params = params, data = dtrain, nrounds = params$nrounds,nthread=30)
    
    train_pred <- predict(model, dtrain)
    test_pred <- predict(model, dtest)
    
    auc_train <- auc(train_target, train_pred, quiet = TRUE)
    auc_test <- auc(test_target, test_pred, quiet = TRUE)
    
    auc_train_list[i] <- auc_train
    auc_test_list[i] <- auc_test
    
    loss_train_list[i] <- loss_eval(train_pred,train_target,data[-index, 'log_ba'])
    loss_test_list[i] <-  loss_eval(test_pred,test_target,data[index, 'log_ba'])
  }
  
  mean_auc_train <- mean(auc_train_list)
  mean_auc_test <- mean(auc_test_list)
  
  mean_loss_train <- mean(loss_train_list)
  mean_loss_test <- mean(loss_test_list)
  
  return(list(mean_auc_test = mean_auc_test, mean_auc_train = mean_auc_train,
              mean_loss_train = mean_loss_train, mean_loss_test=mean_loss_test))
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
  scale_pos_weight = 2,
  colsample_bytree = z_best_params['colsample_bytree'],
  min_child_weight = z_best_params['min_child_weight'],
  subsample = z_best_params['subsample']
)

results_df5 <- xgboost_tune_z(tune_grid)
print(results_df5[order(results_df5$mean_loss_test,decreasing=F),])

param_z <- list(
  nrounds = 100,
  eta = z_best_params['eta'],
  max_depth =  z_best_params['max_depth'],
  gamma = 0,
  scale_pos_weight = 2,
  colsample_bytree = z_best_params['colsample_bytree'],
  min_child_weight = z_best_params['min_child_weight'],
  subsample = z_best_params['subsample']
)

save(param_z,z_best_params, file=file.path(dir.out,'XGBoost_z_council_hyperpara.RData'))

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




xgb_ba_cv   <- function(alpha, xi, kappa, eta, max_depth, min_child_weight=1, subsample=1, colsample_bytree=1) {
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

kfold_cv_pos <- function(data, target.name ,covar.names, params, k) {
  loss_train_list <- numeric(k)
  loss_test_list <- numeric(k)
  
  for (i in 1:k) {
    select.year <- 2011:2013 + 3*(i-1)
    index <- which(data$year %in% select.year )
    train_data <- data[-index, covar.names ]
    train_target <- data[-index, target.name]
    pos.ind.train <- train_target>0
    train_data <- train_data[pos.ind.train,]
    train_target <- train_target[pos.ind.train]
    # print(table(data[-index,'grid.idx']))
    
    test_data <- data[index, covar.names]
    test_target <- data[index, target.name ]
    pos.ind.test <- test_target>0
    test_data <- test_data[pos.ind.test,]
    test_target <- test_target[pos.ind.test]
    
    
    dtrain <- xgb.DMatrix(data = as.matrix(train_data), label = train_target)
    dtest <- xgb.DMatrix(data = as.matrix(test_data), label = test_target)
    
    set.seed(1234)
    model <- xgb.train(params = params, data = dtrain, nrounds = params$nrounds,nthread=30,)
    
    train_pred <- predict(model, dtrain)
    test_pred <- predict(model, dtest)
    
    
    loss_train_list[i] <- mean(loss_eval.egp(train_pred,train_target,xi,kappa, alpha))
    loss_test_list[i] <-  mean(loss_eval.egp(test_pred,test_target,xi,kappa, alpha))
    
  }
  
  mean_loss_train <- mean(loss_train_list)
  mean_loss_test <- mean(loss_test_list)
  
  
  return(list(mean_loss_test = mean_loss_test, mean_loss_train = mean_loss_train))
}

xgboost_tune_ba <- function(tune_grid){
  
  k <- 4
  
  results_df <- data.frame(
    eta = numeric(0),
    max_depth = numeric(0),
    gamma = numeric(0),
    colsample_bytree = numeric(0),
    min_child_weight = numeric(0),
    subsample = numeric(0),
    nrounds = numeric(0),
    mean_loss_train = numeric(0),
    mean_loss_test = numeric(0),
    objective=character(0)
    
  )
  for (i in 1:nrow(tune_grid)) {
    print(i)
    params <- list(
      eta = tune_grid$eta[i],
      max_depth = tune_grid$max_depth[i],
      gamma = tune_grid$gamma[i],
      colsample_bytree = tune_grid$colsample_bytree[i],
      min_child_weight = tune_grid$min_child_weight[i],
      subsample = tune_grid$subsample[i],
      nrounds = tune_grid$nrounds[i],
      objective = custom_obj,
      verbosity = 0
    )
    
    result <- kfold_cv_pos(data.rf, 'log_ba', covar.names,  params, k)
    
    # Store results in the dataframe
    results_df <- rbind(results_df, data.frame(
      eta = tune_grid$eta[i],
      max_depth = tune_grid$max_depth[i],
      gamma = tune_grid$gamma[i],
      colsample_bytree = tune_grid$colsample_bytree[i],
      min_child_weight = tune_grid$min_child_weight[i],
      subsample = tune_grid$subsample[i],
      nrounds = tune_grid$nrounds[i],
      mean_loss_train = result$mean_loss_train,
      mean_loss_test = result$mean_loss_test,
      objective = 'egp_loss'
    ))
  }
  
  return(results_df)
}



bounds_ba <- list(
  alpha = c(0.1,0.9),
  xi = c(-0.3,0.3),
  kappa = c(0.5,5),
  eta = c(0.02, 0.1),
  max_depth = c(2L,6L)
  # min_child_weight = c(1L,15L),
  # subsample = c(0.5, 1.0),
  # colsample_bytree = c(0.5, 1.0)
)

set.seed(1234)
t1 <- Sys.time()
xgb_ba_params <- BayesianOptimization(
  FUN = xgb_ba_cv,
  bounds = bounds_ba,
  init_points = 100, # Number of initial random points
  n_iter = 50,      # Number of iterations for optimization
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



xi <- ba_best_params['xi']
kappa <- ba_best_params['kappa']
alpha <- ba_best_params['alpha']

custom_obj <- create_custom_objective(xi, kappa, alpha)
tune_grid <- expand.grid(
  nrounds = 100,
  eta = ba_best_params['eta'],
  max_depth = ba_best_params['max_depth'],
  gamma = 0,
  colsample_bytree = seq(0.5,1,0.1),
  min_child_weight = 1:10,
  subsample = seq(0.5,1,0.1)
)

results_df_reg <- xgboost_tune_ba(tune_grid)

results_df_reg[order(results_df_reg$mean_loss_test,decreasing=F),]


tune_grid <- expand.grid(
  nrounds = seq(10,100,5),
  eta = ba_best_params['eta'],
  max_depth = ba_best_params['max_depth'],
  gamma = 0,
  colsample_bytree = 0.7,
  min_child_weight = 10,
  subsample = 0.9
)

results_df_reg2 <- xgboost_tune_ba(tune_grid)

results_df_reg2[order(results_df_reg2$mean_loss_test,decreasing=F),]


param_ba <- list(
  nrounds = 100,
  eta = ba_best_params['eta'],
  max_depth =  3,
  gamma = 0,
  colsample_bytree = 0.7,
  min_child_weight = 10,
  subsample = 0.9,
  objective=custom_obj
)

save(param_ba, ba_best_params,custom_obj, file=file.path(dir.out,'XGBoost_ba_council_hyperpara.RData'))

##########################################
trunc_poisson_loss <- function(y,preds){
  lambda <- exp(preds)
  return(mean(lambda - y*preds + log(1-exp(-lambda))))
}


trunc_poisson_obj <- function(preds, dtrain) {
  # preds: Vector of raw predictions from XGBoost (f(x))
  # dtrain: xgb.DMatrix with data
  # label: observed counts (must be > 0 for truncated Poisson)
  y <- getinfo(dtrain, "label")
  
  # Compute lambda = exp(pred)
  lambda <- exp(preds)
  
  # Compute gradient
  # g = lambda - y + (lambda * exp(-lambda)) / (1 - exp(-lambda))
  g <- lambda - y + (lambda * exp(-lambda)) / (1 - exp(-lambda))
  
  # Compute hessian
  # h = lambda + [lambda * exp(-lambda) * (1 - lambda - exp(-lambda))] / (1 - exp(-lambda))^2
  numerator <- lambda * exp(-lambda) * (1 - lambda - exp(-lambda))
  denominator <- (1 - exp(-lambda))^2
  h <- lambda + numerator / denominator
  return(list(grad = g, hess = h))
}


xgb_cnt_cv   <- function(eta, max_depth, min_child_weight, subsample, colsample_bytree) {
  k <- 4
  loss_train_list <- numeric(k)
  loss_test_list <- numeric(k)
  
  data <- data.rf
  target.name <- 'y'
  
  params <- list(
    objective = trunc_poisson_obj,
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
    train_target <- data[-index, target.name]
    pos.ind.train <- train_target>0
    train_data <- train_data[pos.ind.train,]
    train_target <- train_target[pos.ind.train]
    # print(table(data[-index,'grid.idx']))
    
    test_data <- data[index, covar.names]
    test_target <- data[index, target.name ]
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
    
    
    loss_train <- trunc_poisson_loss(train_target,train_pred)
    loss_test <- trunc_poisson_loss(test_target,test_pred)
    
    loss_train_list[i] <- loss_train
    loss_test_list[i] <- loss_test
    
  }
  
  mean_loss_train <- mean(loss_train_list)
  mean_loss_test <- mean(loss_test_list)
  
  
  return(list(Score=-mean_loss_test,Pred=NULL))
}

bounds_cnt <- list(
  eta = c(0.03, 0.1),
  max_depth = c(2L,6L),
  min_child_weight = c(1L,10L),
  subsample = c(0.5, 1.0),
  colsample_bytree = c(0.5, 1.0)
)

set.seed(1234)
t1 <- Sys.time()
xgb_cnt_params <- BayesianOptimization(
  FUN = xgb_cnt_cv,
  bounds = bounds_cnt,
  init_points = 50, # Number of initial random points
  n_iter = 30,      # Number of iterations for optimization
  acq = "ucb",      # Acquisition function, e.g., "ucb", "ei", "poi"
  kappa = 2.576,    # Parameter for the acquisition function
  eps = 0.0,        # Exploration-exploitation trade-off
  verbose = TRUE
)
t2 <- Sys.time()
print(t2-t1)

cnt_best_params <- xgb_cnt_params$Best_Par
print(cnt_best_params)
print( xgb_cnt_params$Best_Value)



kfold_cv_poisson <- function(data, target.name ,covar.names, params, k) {
  loss_train_list <- numeric(k)
  loss_test_list <- numeric(k)
  
  for (i in 1:k) {
    select.year <- 2011:2013 + 3*(i-1)
    index <- which(data$year %in% select.year )
    train_data <- data[-index, covar.names ]
    train_target <- data[-index, target.name]
    pos.ind.train <- train_target>0
    train_data <- train_data[pos.ind.train,]
    train_target <- train_target[pos.ind.train]
    # print(table(data[-index,'grid.idx']))
    
    test_data <- data[index, covar.names]
    test_target <- data[index, target.name ]
    pos.ind.test <- test_target>0
    test_data <- test_data[pos.ind.test,]
    test_target <- test_target[pos.ind.test]
    
    
    dtrain <- xgb.DMatrix(data = as.matrix(train_data), label = train_target)
    dtest <- xgb.DMatrix(data = as.matrix(test_data), label = test_target)
    
    set.seed(1234)
    model <- xgb.train(params = params, data = dtrain, nrounds = params$nrounds,nthread=30)
    
    train_pred <- predict(model, dtrain)
    test_pred <- predict(model, dtest)
    
    loss_train <- trunc_poisson_loss(train_target,train_pred)
    loss_test <- trunc_poisson_loss(test_target,test_pred)
    
    loss_train_list[i] <- loss_train
    loss_test_list[i] <- loss_test
  }
  
  mean_loss_train <- mean(loss_train_list)
  mean_loss_test <- mean(loss_test_list)
  
  
  return(list(mean_loss_test = mean_loss_test, mean_loss_train = mean_loss_train))
}

xgboost_tune_cnt <- function(tune_grid){
  
  k <- 4
  
  results_df <- data.frame(
    eta = numeric(0),
    max_depth = numeric(0),
    gamma = numeric(0),
    colsample_bytree = numeric(0),
    min_child_weight = numeric(0),
    subsample = numeric(0),
    nrounds = numeric(0),
    mean_loss_train = numeric(0),
    mean_loss_test = numeric(0),
    objective=character(0)
    
  )
  for (i in 1:nrow(tune_grid)) {
    print(i)
    params <- list(
      eta = tune_grid$eta[i],
      max_depth = tune_grid$max_depth[i],
      gamma = tune_grid$gamma[i],
      colsample_bytree = tune_grid$colsample_bytree[i],
      min_child_weight = tune_grid$min_child_weight[i],
      subsample = tune_grid$subsample[i],
      nrounds = tune_grid$nrounds[i],
      objective = trunc_poisson_obj,
      verbosity = 0
    )
    
    result <- kfold_cv_poisson(data.rf, 'y', covar.names,  params, k)
    
    # Store results in the dataframe
    results_df <- rbind(results_df, data.frame(
      eta = tune_grid$eta[i],
      max_depth = tune_grid$max_depth[i],
      gamma = tune_grid$gamma[i],
      colsample_bytree = tune_grid$colsample_bytree[i],
      min_child_weight = tune_grid$min_child_weight[i],
      subsample = tune_grid$subsample[i],
      nrounds = tune_grid$nrounds[i],
      mean_loss_train = result$mean_loss_train,
      mean_loss_test = result$mean_loss_test,
      objective = tune_grid$objective[i]
    ))
  }
  
  return(results_df)
}



tune_grid <- expand.grid(
  nrounds = seq(10,100,5),
  eta = cnt_best_params['eta'],
  max_depth = cnt_best_params['max_depth'],
  gamma = 0,
  colsample_bytree = cnt_best_params['colsample_bytree'],
  min_child_weight = cnt_best_params['min_child_weight'],
  subsample = cnt_best_params['subsample'],
  objective=c('count:poisson')
)


results_df_pois_1 <- xgboost_tune_cnt(tune_grid)

results_df_pois_1[order(results_df_pois_1$mean_loss_test,decreasing=F),]


param_cnt <- list(
  nrounds = 100,
  eta = cnt_best_params['eta'],
  max_depth = cnt_best_params['max_depth'],
  gamma = 0,
  colsample_bytree = cnt_best_params['colsample_bytree'],
  min_child_weight = cnt_best_params['min_child_weight'],
  subsample = cnt_best_params['subsample'],
  objective=trunc_poisson_obj
)

save(param_cnt, cnt_best_params, 
     file=file.path(dir.out,'XGBoost_cnt_council_hyperpara.RData'))