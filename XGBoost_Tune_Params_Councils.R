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


data.fit$z <- as.vector((data.fit$y>0)+0)
data.fit$HVegTyp <- as.factor(data.fit$HVegTyp)
# levels(data.fit$HVegTyp) <-  c('NoCover','Evergreen_broadleaf_trees','Mixed_forest','Interrupted_forest')
data.fit$LVegTyp <- as.factor(data.fit$LVegTyp)
# levels(data.fit$LVegTyp) <- c('NoCover','Crops','Tall_grass','Semidesert','Evergreen_shrubs')
data.fit$month <- as.factor(data.fit$month)
data.fit <- dummy_cols(data.fit, 
           select_columns = c('HVegTyp','LVegTyp','month'),remove_first_dummy=TRUE)

# covar.names <- colnames(data.fit)[c(-1,-2,-5,-6,-7,-8,-13,-16,-24)]
covar.names <- c('lon.grid','lat.grid','year',
                 'FWI','HVegCov','HVegLAI','LVegCov','LVegLAI',
                 'Pricp','RHumi','Temp','UComp','VComp','DewPoint','WindSpeed',
                 'HVegTyp_6','HVegTyp_18','HVegTyp_19',
                 'LVegTyp_1','LVegTyp_7','LVegTyp_11','LVegTyp_16',
                 'month_2','month_3','month_4','month_5','month_6','month_7','month_8','month_9','month_10','month_11','month_12'
                 )

data.rf <- data.fit[data.fit$year <= 2022, c(covar.names,'grid.idx','time.idx','log_ba','z','y')]
data.rf[is.na(data.rf$log_ba),'log_ba'] <- 0

library(xgboost)
library(pROC)

pos.loss <- mean(data.rf[data.rf$y>0,'log_ba'])
loss_eval <- function(pred,y,log_ba){
  pos.idx <- which(y>0)
  l1 <- sum((pred[pos.idx]*log_ba[pos.idx] - log_ba[pos.idx])^2)
  l2 <- sum((pred[-pos.idx]*pos.loss)^2)
  return((l1+l2)/length(pred))
}

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


##################Tune max_depth, min_child_weight, scale_pos_weight #########
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
tune_grid <- expand.grid(
  nrounds = 100,
  eta = seq(0.04,0.07,0.01),
  max_depth = c(4,5,6,7),
  gamma = 0,
  scale_pos_weight = c(7),
  colsample_bytree = 0.8,
  min_child_weight = 1,
  subsample = 0.7
)

tune_grid <- expand.grid(
  nrounds = 100,
  eta = seq(0.06,0.09,0.01),
  max_depth = c(6,7,8,9),
  gamma = 0,
  scale_pos_weight = c(6,7),
  colsample_bytree = 0.8,
  min_child_weight = 1,
  subsample = 0.7
)


tune_grid <- expand.grid(
  nrounds = 100,
  eta = seq(0.07,0.1,0.01),
  max_depth = c(5,6,7,8),
  gamma = 0,
  scale_pos_weight = c(2),
  colsample_bytree = 0.8,
  min_child_weight = 1,
  subsample = 0.7
)


results_df1 <- xgboost_tune_z(tune_grid)
print(results_df1[order(results_df1$mean_auc_test,decreasing=T),])
print(results_df1[order(results_df1$mean_loss_test,decreasing=F),])
#

tune_grid <- expand.grid(
  nrounds = 100,
  eta = 0.09,
  max_depth = 7,
  gamma = 0,
  scale_pos_weight = 2,
  colsample_bytree = 0.8,
  min_child_weight = 1:10,
  subsample = 0.7
)

results_df2 <- xgboost_tune_z(tune_grid)

print(results_df2[order(results_df2$mean_auc_test,decreasing=T),])
print(results_df2[order(results_df2$mean_loss_test,decreasing=F),])
#
# save(results_df3, file=file.path(dir.out, 'XGBoost_Param_z_3.RData'))


tune_grid <- expand.grid(
  nrounds = 100,
  eta = 0.09,
  max_depth = 7,
  gamma = 0,
  scale_pos_weight = 2,
  colsample_bytree = c(0.4,0.5,0.6,0.7,0.8,0.9,1),
  min_child_weight = 1,
  subsample = c(0.4,0.5,0.6,0.7,0.8,0.9,1)
)

results_df3 <- xgboost_tune_z(tune_grid)

print(results_df3[order(results_df3$mean_auc_test,decreasing=T),])
print(results_df3[order(results_df3$mean_loss_test,decreasing=F),])
# save(results_df4, file=file.path(dir.out, 'XGBoost_Param_z_4.RData'))


tune_grid <- expand.grid(
  nrounds = 100,
  eta = 0.09,
  max_depth = 7,
  gamma = seq(0,1,0.1),
  scale_pos_weight = 2,
  colsample_bytree = 0.4,
  min_child_weight = 1,
  subsample = 1
)

results_df4 <- xgboost_tune_z(tune_grid)

print(results_df4[order(results_df4$mean_auc_test,decreasing=T),])
print(results_df4[order(results_df4$mean_loss_test,decreasing=F),])
# save(results_df5, file=file.path(dir.out, 'XGBoost_Param_z_5.RData'))

tune_grid <- expand.grid(
  nrounds = seq(10,100,5),
  eta = 0.09,
  max_depth = 7,
  gamma = 0.7,
  scale_pos_weight = 2,
  colsample_bytree = 0.4,
  min_child_weight = 1,
  subsample = 1
)
#
results_df5 <- xgboost_tune_z(tune_grid)
print(results_df5[order(results_df5$mean_auc_test,decreasing=T),])
print(results_df5[order(results_df5$mean_loss_test,decreasing=F),])
# save(results_df6, file=file.path(dir.out, 'XGBoost_Param_z_6.RData'))


# load(file.path(dir.out, 'XGBoost_Param_z_6.RData'))

plot(results_df5$nrounds,results_df5$mean_auc_train,type='l',col='blue',ylim=c(0.85,0.95))
lines(results_df5$nrounds,results_df5$mean_auc_test,col='red',type='l')

tune_grid_z <- expand.grid(
  nrounds = 100,
  eta = 0.09,
  max_depth = 7,
  gamma = 0.7,
  scale_pos_weight = 2,
  colsample_bytree = 0.4,
  min_child_weight = 1,
  subsample = 1
)



################################## Log BA ###########################
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
    model <- xgb.train(params = params, data = dtrain, nrounds = params$nrounds,nthread=30)

    train_pred <- predict(model, dtrain)
    test_pred <- predict(model, dtest)

    loss_train <- mean((train_pred-train_target)^2)
    loss_test <- mean((test_pred-test_target)^2)

    loss_train_list[i] <- loss_train
    loss_test_list[i] <- loss_test
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
      objective = tune_grid$objective[i],
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
      objective = tune_grid$objective[i]
    ))
  }

  return(results_df)
}
# 
# custom_loss <- function(preds, dtrain, alpha) {
#   # Get the true labels
#   labels <- getinfo(dtrain, "label")
#   
#   # Calculate the gradient
#   grad <- 2 * alpha * (preds - labels)
#   
#   # Calculate the hessian
#   hess <- rep(2 * alpha, length(labels))
#   
#   # Return as a list
#   return(list(grad = grad, hess = hess))
# }
# 
# 
# train_with_custom_loss <- function(alpha, dtrain, params) {
#   # Define the custom objective function with alpha
#   custom_objective <- function(preds, dtrain) {
#     custom_loss(preds, dtrain, alpha)
#   }
#   # Train the XGBoost model
#   model <- xgboost(
#     data = dtrain,
#     params = params,
#     verbose = 0  # Turn off printing for cleaner output
#   )
#   return(model)
# }

tune_grid <- expand.grid(
  nrounds = 100,
  eta = seq(0.02,0.1,0.01),
  max_depth = c(2,3,4),
  gamma = 0,
  colsample_bytree = 0.8,
  min_child_weight = 1,
  subsample = 0.7,
  objective=c('reg:squarederror','reg:gamma')
)


results_df_reg <- xgboost_tune_ba(tune_grid)

results_df_reg[order(results_df_reg$mean_loss_test,decreasing=F),]




tune_grid <- expand.grid(
  nrounds = 100,
  eta = 0.05,
  max_depth = 3,
  gamma = 0,
  colsample_bytree = 0.8,
  min_child_weight = 1:10,
  subsample = 0.7,
  objective=c('reg:squarederror')
)

results_df_reg_1 <- xgboost_tune_ba(tune_grid)

results_df_reg_1[order(results_df_reg_1$mean_loss_test,decreasing=F),]



tune_grid <- expand.grid(
  nrounds = 100,
  eta = 0.05,
  max_depth = 3,
  gamma = 0,
  colsample_bytree = seq(0.4,1,0.1),
  min_child_weight = 5,
  subsample = seq(0.4,1,0.1),
  objective=c('reg:squarederror')
)


results_df_reg_2 <- xgboost_tune_ba(tune_grid)

results_df_reg_2[order(results_df_reg_2$mean_loss_test,decreasing=F),]




tune_grid <- expand.grid(
  nrounds = 100,
  eta = 0.05,
  max_depth = 3,
  gamma = seq(0,1,0.1),
  colsample_bytree = 0.4,
  min_child_weight = 5,
  subsample = 0.5,
  objective=c('reg:squarederror')
)


results_df_reg_4 <- xgboost_tune_ba(tune_grid)

results_df_reg_4[order(results_df_reg_4$mean_loss_test,decreasing=F),]


tune_grid <- expand.grid(
  nrounds = seq(10,100,5),
  eta = 0.05,
  max_depth = 3,
  gamma = 0.2,
  colsample_bytree = 0.4,
  min_child_weight = 5,
  subsample = 0.5,
  objective=c('reg:squarederror')
)


results_df_reg_5 <- xgboost_tune_ba(tune_grid)
results_df_reg_5
results_df_reg_5[order(results_df_reg_5$mean_loss_test,decreasing=F),]



plot(results_df_reg_5$nrounds,results_df_reg_5$mean_loss_train,type='l',col='blue')
lines(results_df_reg_5$nrounds,results_df_reg_5$mean_loss_test,col='red')

tune_grid <- expand.grid(
  nrounds = 100,
  eta = 0.05,
  max_depth = 3,
  gamma = 0.2,
  colsample_bytree = 0.4,
  min_child_weight = 5,
  subsample = 0.5,
  objective=c('reg:squarederror')
)


################################### Count ###############################

poisson_loss <- function(y,lambda){
  return(mean(-y*log(lambda)+lambda))
}

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
    
    loss_train <- poisson_loss(train_target,train_pred)
    loss_test <- poisson_loss(test_target,test_pred)
    
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
      objective = tune_grid$objective[i],
      verbosity = 0,
      eval_metric = "logloss"
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
  nrounds = 100,
  eta = seq(0.04,0.1,0.01),
  max_depth = 5:9,
  gamma = 0,
  colsample_bytree =0.8,
  min_child_weight = 1,
  subsample = 0.7,
  objective=c('count:poisson')
)


results_df_pois_1 <- xgboost_tune_cnt(tune_grid)

results_df_pois_1[order(results_df_pois_1$mean_loss_test,decreasing=F),]



tune_grid <- expand.grid(
  nrounds = 100,
  eta = 0.09,
  max_depth = 6,
  gamma = 0,
  colsample_bytree = 0.8,
  min_child_weight = 1:15,
  subsample = 0.7,
  objective=c('count:poisson')
)
results_df_pois_3 <- xgboost_tune_cnt(tune_grid)

results_df_pois_3[order(results_df_pois_3$mean_loss_test,decreasing=F),]


tune_grid <- expand.grid(
  nrounds = 100,
  eta = 0.09,
  max_depth = 6,
  gamma = 0,
  colsample_bytree = seq(0.4,1,0.1),
  min_child_weight = 1,
  subsample = seq(0.4,1,0.1),
  objective=c('count:poisson')
)
results_df_pois_4 <- xgboost_tune_cnt(tune_grid)

results_df_pois_4[order(results_df_pois_4$mean_loss_test,decreasing=F),]


tune_grid <- expand.grid(
  nrounds = 100,
  eta = 0.09,
  max_depth = 6,
  gamma = seq(0,1,0.1),
  colsample_bytree = 0.7,
  min_child_weight = 1,
  subsample = 0.5,
  objective=c('count:poisson')
)
results_df_pois_5 <- xgboost_tune_cnt(tune_grid)

results_df_pois_5[order(results_df_pois_5$mean_loss_test,decreasing=F),]


tune_grid <- expand.grid(
  nrounds = seq(10,100,5),
  eta = 0.09,
  max_depth = 6,
  gamma = 0,
  colsample_bytree = 0.7,
  min_child_weight = 1,
  subsample = 0.5,
  objective=c('count:poisson')
)
results_df_pois_6 <- xgboost_tune_cnt(tune_grid)

results_df_pois_6[order(results_df_pois_6$mean_loss_test,decreasing=F),]


plot(results_df_pois_6$nrounds,results_df_pois_6$mean_loss_train,type='l',col='blue')
lines(results_df_pois_6$nrounds,results_df_pois_6$mean_loss_test,col='red')

tune_grid <- expand.grid(
  nrounds = 95,
  eta = 0.09,
  max_depth = 6,
  gamma = 0,
  colsample_bytree = 0.7,
  min_child_weight = 1,
  subsample = 0.5,
  objective=c('count:poisson')
)
