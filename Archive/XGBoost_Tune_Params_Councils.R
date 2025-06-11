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
# data.fit$month <- as.factor(data.fit$month)
data.fit <- dummy_cols(data.fit, 
                       select_columns = c('HVegTyp','LVegTyp','NAME_1'),remove_first_dummy=TRUE)

data.fit$year.af.2017.cnt = as.integer(data.fit$year>2017)
# covar.names <- colnames(data.fit)[c(-1,-2,-5,-6,-7,-8,-13,-16,-24)]
covar.names <- c('lon.grid','lat.grid','year',
                 'FWI','HVegCov','HVegLAI','LVegCov','LVegLAI',
                 'Pricp','RHumi','Temp','UComp','VComp','DewPoint','WindSpeed',
                 'HVegTyp_6','HVegTyp_18','HVegTyp_19',
                 'LVegTyp_1','LVegTyp_7','LVegTyp_11','LVegTyp_16',
                 'month',
                 "NAME_1_Beja",
                 "NAME_1_Braga", "NAME_1_Bragança", "NAME_1_Castelo Branco",
                 "NAME_1_Coimbra", "NAME_1_Évora", "NAME_1_Faro",
                 "NAME_1_Guarda", "NAME_1_Leiria", "NAME_1_Lisboa",
                 "NAME_1_Portalegre","NAME_1_Porto",
                 "NAME_1_Santarém", "NAME_1_Setúbal", "NAME_1_Viana do Castelo",
                 "NAME_1_Vila Real", "NAME_1_Viseu"
)



data.rf <- data.fit[data.fit$year <= 2022, c(covar.names,'grid.idx','time.idx','z','y','area_ha','log_ba')]
data.rf[is.na(data.rf$area_ha),'area_ha'] <- 0
data.rf[is.na(data.rf$log_ba),'log_ba'] <- 0


library(xgboost)
library(rBayesianOptimization)
library(pROC)

set.seed(1234)


weighted_loss <- function(y, preds, type){
  if (type=='CNT'){
    thres.cnt <- c(0:15,17,19,21,23,25,30)
    loss <- rep(NA,length(thres.cnt))
    w <- 1 - (1+(thres.cnt+1)^2/1000)^{-1/4}
    w <- w/w[length(thres.cnt)]
    for (i in 1:length(thres.cnt)){
      thres <- thres.cnt[i]
      I.true <-  y <= thres
      I.pred <-  preds<= thres
      loss[i] <- sum(w[i]*(I.true - I.pred)^2) 
    }
  } else if (type=='BA'){
    thres.ba <- c( seq(0,100,10),150, 200, 300, 400, 500, 1000, 1500, 2000, 5000, 10000,20000)
    loss <- rep(NA,length(thres.ba))
    w <- 1 - (1+(thres.ba+1)/1000)^{-1/4}
    w <- w/w[length(thres.ba)]
    for (i in 1:length(thres.ba)){
      thres <- thres.ba[i]
      I.true <-  y <= thres
      I.pred <-  preds<= thres
      loss[i] <- sum(w[i]*(I.true - I.pred)^2) 
    }
    
  }else{
    loss <- 0
  }
  return(sum(loss)) 
}

median.ba <- as.numeric(quantile(data.rf[data.rf$y>0,'area_ha'],0.5))
median.cnt <- as.numeric(quantile(data.rf[data.rf$y>0,'y'],0.5))

weighted_loss_z <- function(y, preds, cnt, area_ha){
  pos.idx <- which(y>0)
  loss.cnt.1 <- weighted_loss(cnt[pos.idx], preds[pos.idx]*cnt[pos.idx], type='CNT')
  loss.ba.1 <- weighted_loss(area_ha[pos.idx], preds[pos.idx]*area_ha[pos.idx], type='BA')

  loss.cnt.2 <- weighted_loss(0, preds[-pos.idx]*median.cnt, type='CNT')
  loss.ba.2 <- weighted_loss(0, preds[-pos.idx]*median.ba, type='BA')

  return(loss.cnt.1 + loss.ba.1 + loss.cnt.2 + loss.ba.2)
}

kfold_cv  <- function(data, target.name ,covar.names, params, k) {
  auc_train_list <- numeric(k)
  auc_test_list <- numeric(k)
  
  loss_train_list <- numeric(k)
  loss_test_list <- numeric(k)
  
  for (i in 1:k) {
    select.year <- 2011 + (i-1)
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
    
    loss_train_list[i] <- weighted_loss_z(train_target, train_pred,data[-index, 'y'], data[-index, 'area_ha'])
    loss_test_list[i] <-  weighted_loss_z(test_target, test_pred,data[index, 'y'], data[index, 'area_ha'])
  }
  
  mean_auc_train <- mean(auc_train_list)
  mean_auc_test <- mean(auc_test_list)
  
  mean_loss_train <- sum(loss_train_list)/k
  mean_loss_test <- sum(loss_test_list)
  
  return(list(mean_auc_test = mean_auc_test, mean_auc_train = mean_auc_train,
              mean_loss_train = mean_loss_train, mean_loss_test=mean_loss_test))
}


##################Tune max_depth, min_child_weight, scale_pos_weight #########
xgboost_tune_z <- function(tune_grid){
  
  k <- 12
  
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
  eta = seq(0.03,0.1,0.01),
  max_depth = 3:10,
  gamma = 0,
  scale_pos_weight = 6.87,
  colsample_bytree = 0.6,
  min_child_weight = 1,
  subsample = 0.9
)



results_df1 <- xgboost_tune_z(tune_grid)
print(results_df1[order(results_df1$mean_auc_test,decreasing=T),])
print(results_df1[order(results_df1$mean_loss_test,decreasing=F),])
#

tune_grid <- expand.grid(
  nrounds = 100,
  eta = 0.07,
  max_depth = 7,
  gamma = 0,
  scale_pos_weight = 1:15,
  colsample_bytree = 0.7,
  min_child_weight = 1,
  subsample = 0.8
)

results_df1_1 <- xgboost_tune_z(tune_grid)
print(results_df1_1[order(results_df1_1$mean_auc_test,decreasing=T),])
print(results_df1_1[order(results_df1_1$mean_loss_test,decreasing=F),])
#

tune_grid <- expand.grid(
  nrounds = 100,
  eta = 0.1,
  max_depth = 6,
  gamma = 0,
  scale_pos_weight = 6.87,
  colsample_bytree = 0.6,
  min_child_weight = 1:25,
  subsample = 0.9
)

results_df2 <- xgboost_tune_z(tune_grid)

print(results_df2[order(results_df2$mean_auc_test,decreasing=T),])
print(results_df2[order(results_df2$mean_loss_test,decreasing=F),])
#
# save(results_df3, file=file.path(dir.out, 'XGBoost_Param_z_3.RData'))



tune_grid <- expand.grid(
  nrounds = 100,
  eta = 0.1,
  max_depth = 6,
  gamma = 0,
  scale_pos_weight = 6.87,
  colsample_bytree = c(0.4,0.5,0.6,0.7,0.8,0.9,1),
  min_child_weight = 18,
  subsample = c(0.4,0.5,0.6,0.7,0.8,0.9,1)
)

results_df3 <- xgboost_tune_z(tune_grid)

print(results_df3[order(results_df3$mean_auc_test,decreasing=T),])
print(results_df3[order(results_df3$mean_loss_test,decreasing=F),])
# save(results_df4, file=file.path(dir.out, 'XGBoost_Param_z_4.RData'))


tune_grid <- expand.grid(
  nrounds = 100,
  eta = 0.1,
  max_depth = 6,
  gamma = seq(0,2,0.1),
  scale_pos_weight = 6.87,
  colsample_bytree = 0.4,
  min_child_weight = 18,
  subsample = 0.9
)

results_df4 <- xgboost_tune_z(tune_grid)

print(results_df4[order(results_df4$mean_auc_test,decreasing=T),])
print(results_df4[order(results_df4$mean_loss_test,decreasing=F),])
# save(results_df5, file=file.path(dir.out, 'XGBoost_Param_z_5.RData'))

tune_grid <- expand.grid(
  nrounds = seq(10,100,5),
  eta = 0.1,
  max_depth = 6,
  gamma = 1.9,
  scale_pos_weight = 6.87,
  colsample_bytree = 0.4,
  min_child_weight = 18,
  subsample = 0.9
)
#
results_df5 <- xgboost_tune_z(tune_grid)
print(results_df5[order(results_df5$mean_auc_test,decreasing=T),])
print(results_df5[order(results_df5$mean_loss_test,decreasing=F),])
save(results_df5, file=file.path(dir.out, 'XGBoost_Param_z_5.RData'))


# load(file.path(dir.out, 'XGBoost_Param_z_6.RData'))

plot(results_df5$nrounds,results_df5$mean_auc_train,type='l',col='blue',ylim=c(0.85,0.95))
lines(results_df5$nrounds,results_df5$mean_auc_test,col='red',type='l')

params_z <- list(
  nrounds = 100,
  eta = 0.1,
  max_depth = 6,
  gamma = 1.9,
  scale_pos_weight = 6.87,
  colsample_bytree = 0.4,
  min_child_weight = 18,
  subsample = 0.9,
  objective = "binary:logistic",
  verbosity = 0
)


################################## Log BA ###########################
kfold_cv_pos <- function(data, target.name ,covar.names, params, k) {
  loss_train_list <- numeric(k)
  loss_test_list <- numeric(k)
  
  wloss_train_list <- numeric(k)
  wloss_test_list <- numeric(k)

  for (i in 1:k) {
    select.year <- 2011 + (i-1)
    index <- which(data$year %in% select.year )
    train_data <- data[-index, covar.names ]
    train_target <- data[-index, target.name]
    test_data <- data[index, covar.names]
    test_target <- data[index, target.name ]
    
    
    
    pos.ind.train <- train_target>0
    train_data <- train_data[pos.ind.train,]
    train_target <- train_target[pos.ind.train]
    # print(table(data[-index,'grid.idx']))

    pos.ind.test <- test_target>0
    test_data <- test_data[pos.ind.test,]
    test_target <- test_target[pos.ind.test]


    dtrain <- xgb.DMatrix(data = as.matrix(train_data), label = train_target)
    dtest <- xgb.DMatrix(data = as.matrix(test_data), label = test_target)

    set.seed(1234)
    model <- xgb.train(params = params, data = dtrain, nrounds = params$nrounds,nthread=30)

    train_pred <- predict(model, dtrain)
    test_pred <- predict(model, dtest)

    # loss_train <- weighted_loss(exp(train_target), exp(train_pred),'BA')
    # loss_test <- weighted_loss(exp(test_target), exp(test_pred), 'BA')
    
    loss_train <- sum((train_target - train_pred)^2)
    loss_test <- sum((test_target - test_pred)^2)
    
    wloss_train <- weighted_loss(exp(train_target),exp(train_pred),'BA')
    wloss_test <- weighted_loss(exp(test_target), exp(test_pred), 'BA')

    loss_train_list[i] <- loss_train
    loss_test_list[i] <- loss_test
    
    wloss_train_list[i] <- wloss_train
    wloss_test_list[i] <- wloss_test
  }

  mean_loss_train <- sum(loss_train_list)/k
  mean_loss_test <- sum(loss_test_list)
  
  mean_wloss_train <- sum(wloss_train_list)/k
  mean_wloss_test <- sum(wloss_test_list)


  return(list(mean_loss_test = mean_loss_test, mean_loss_train = mean_loss_train,
              weighted_loss_train = mean_wloss_train, weighted_loss_test =mean_wloss_test ))
}

xgboost_tune_ba <- function(tune_grid){

  k <- 12

  results_df <- data.frame(
    eta = numeric(0),
    max_depth = numeric(0),
    gamma = numeric(0),
    colsample_bytree = numeric(0),
    min_child_weight = numeric(0),
    subsample = numeric(0),
    nrounds = numeric(0),
    loss_train = numeric(0),
    loss_test = numeric(0),
    weighted_loss_train = numeric(0),
    weighed_loss_test = numeric(0),
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
      loss_train = result$mean_loss_train,
      loss_test = result$mean_loss_test,
      objective = tune_grid$objective[i],
      weighted_loss_train = result$weighted_loss_train,
      weighted_loss_test = result$weighted_loss_test
      
    ))
  }

  return(results_df)
}

tune_grid <- expand.grid(
  nrounds = 100,
  eta = seq(0.02,0.1,0.01),
  max_depth = c(2,3,4,5),
  gamma = 0,
  colsample_bytree = 0.9,
  min_child_weight = 1,
  subsample = 0.6,
  objective=c('reg:squarederror')
)


results_df_reg <- xgboost_tune_ba(tune_grid)

results_df_reg[order(results_df_reg$loss_test,decreasing=F),]




tune_grid <- expand.grid(
  nrounds = 100,
  eta = 0.07,
  max_depth = 3,
  gamma = 0,
  colsample_bytree = 0.9,
  min_child_weight = 1:10,
  subsample = 0.6,
  objective=c('reg:squarederror')
)

results_df_reg_1 <- xgboost_tune_ba(tune_grid)

results_df_reg_1[order(results_df_reg_1$loss_test,decreasing=F),]



tune_grid <- expand.grid(
  nrounds = 100,
  eta = 0.07,
  max_depth = 3,
  gamma = 0,
  colsample_bytree = seq(0.4,1,0.1),
  min_child_weight = 3,
  subsample = seq(0.4,1,0.1),
  objective=c('reg:squarederror')
)


results_df_reg_2 <- xgboost_tune_ba(tune_grid)

results_df_reg_2[order(results_df_reg_2$loss_test,decreasing=F),]




tune_grid <- expand.grid(
  nrounds = 100,
  eta = 0.07,
  max_depth = 3,
  gamma = 0:10,
  colsample_bytree = 0.9,
  min_child_weight = 3,
  subsample = 0.6,
  objective=c('reg:squarederror')
)


results_df_reg_4 <- xgboost_tune_ba(tune_grid)

results_df_reg_4[order(results_df_reg_4$loss_test,decreasing=F),]


tune_grid <- expand.grid(
  nrounds = seq(10,100,5),
  eta = 0.07,
  max_depth = 3,
  gamma = 2,
  colsample_bytree = 0.9,
  min_child_weight = 3,
  subsample = 0.6,
  objective=c('reg:squarederror')
)


results_df_reg_5 <- xgboost_tune_ba(tune_grid)
results_df_reg_5[order(results_df_reg_5$loss_test,decreasing=F),]



plot(results_df_reg_5$nrounds,results_df_reg_5$mean_loss_train,type='l',col='blue')
plot(results_df_reg_5$nrounds,results_df_reg_5$mean_loss_test,col='red')
save(results_df_reg_5, file=file.path(dir.out, 'XGBoost_Param_ba_5.RData'))


params_ba <- list(
  nrounds = 95,
  eta = 0.07,
  max_depth =  3,
  gamma = 2,
  colsample_bytree = 0.9,
  min_child_weight = 3,
  subsample = 0.6,
  objective='reg:squarederror'
)

################################### Count ###############################

trunc_poisson_loss <- function(y,preds){
  lambda <- exp(preds)
  return(sum(lambda - y*preds + log(1-exp(-lambda))))
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




kfold_cv_poisson <- function(data, target.name ,covar.names, params, k) {
  loss_train_list <- numeric(k)
  loss_test_list <- numeric(k)
  wloss_train_list <- numeric(k)
  wloss_test_list <- numeric(k)
  
  for (i in 1:k) {
    select.year <- 2011 + (i-1)
    index <- which(data$year %in% select.year )
    train_data <- data[-index, covar.names ]
    train_target <- data[-index, target.name]
    test_data <- data[index, covar.names]
    test_target <- data[index, target.name ]
    
    
    pos.ind.train <- train_target>0
    train_data <- train_data[pos.ind.train,]
    train_target <- train_target[pos.ind.train]
    # print(table(data[-index,'grid.idx']))
    
    pos.ind.test <- test_target>0
    test_data <- test_data[pos.ind.test,]
    test_target <- test_target[pos.ind.test]
    
    
    dtrain <- xgb.DMatrix(data = as.matrix(train_data), label = train_target)
    dtest <- xgb.DMatrix(data = as.matrix(test_data), label = test_target)
    
    set.seed(1234)
    model <- xgb.train(params = params, data = dtrain, nrounds = params$nrounds,nthread=30)
    
    train_pred <- predict(model, dtrain)
    test_pred <- predict(model, dtest)
    
    # loss_train <- weighted_loss(train_target,train_pred,type='CNT')
    # loss_test <- weighted_loss(test_target,test_pred,type='CNT')
    
    loss_train <- trunc_poisson_loss(train_target,train_pred)
    loss_test <- trunc_poisson_loss(test_target,test_pred)
    
    wloss_train <- weighted_loss((train_target),exp(train_pred),'CNT')
    wloss_test <- weighted_loss((test_target), exp(test_pred), 'CNT')
    
    loss_train_list[i] <- loss_train
    loss_test_list[i] <- loss_test
    
    wloss_train_list[i] <- wloss_train
    wloss_test_list[i] <- wloss_test
  }
  
  mean_loss_train <- sum(loss_train_list)/k
  mean_loss_test <- sum(loss_test_list)
  
  mean_wloss_train <- sum(wloss_train_list)/k
  mean_wloss_test <- sum(wloss_test_list)
  
  
  return(list(mean_loss_test = mean_loss_test, mean_loss_train = mean_loss_train,
              mean_wloss_train=mean_wloss_train, mean_wloss_test = mean_wloss_test ))
}

xgboost_tune_cnt <- function(tune_grid){

  k <- 12

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
    weighted_loss_train = numeric(0),
    weighted_loss_test = numeric(0),
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
      weighted_loss_train = result$mean_wloss_train,
      weighted_loss_test = result$mean_wloss_test,
      objective = 'truncated_poisson'
    ))
  }

  return(results_df)
}



tune_grid <- expand.grid(
  nrounds = 100,
  eta = seq(0.01,0.1,0.01),
  max_depth = 2:7,
  gamma = 0,
  colsample_bytree =0.7,
  min_child_weight = 1,
  subsample = 1
)


results_df_pois_1 <- xgboost_tune_cnt(tune_grid)

results_df_pois_1[order(results_df_pois_1$mean_loss_test,decreasing=F),]



tune_grid <- expand.grid(
  nrounds = 100,
  eta = 0.08,
  max_depth = 4,
  gamma = 0,
  colsample_bytree = 0.7,
  min_child_weight = 1:15,
  subsample = 1
)
results_df_pois_3 <- xgboost_tune_cnt(tune_grid)

results_df_pois_3[order(results_df_pois_3$mean_loss_test,decreasing=F),]


tune_grid <- expand.grid(
  nrounds = 100,
  eta = 0.08,
  max_depth = 4,
  gamma = 0,
  colsample_bytree = seq(0.4,1,0.1),
  min_child_weight = 1,
  subsample = seq(0.4,1,0.1)
)
results_df_pois_4 <- xgboost_tune_cnt(tune_grid)

results_df_pois_4[order(results_df_pois_4$mean_loss_test,decreasing=F),]


tune_grid <- expand.grid(
  nrounds = 100,
  eta = 0.08,
  max_depth = 4,
  gamma = seq(0,1,0.1),
  colsample_bytree = 0.7,
  min_child_weight = 1,
  subsample = 1
)
results_df_pois_5 <- xgboost_tune_cnt(tune_grid)

results_df_pois_5[order(results_df_pois_5$mean_loss_test,decreasing=F),]


tune_grid <- expand.grid(
  nrounds = seq(10,100,5),
  eta = 0.08,
  max_depth = 4,
  gamma = 0.4,
  colsample_bytree = 0.7,
  min_child_weight = 1,
  subsample = 1
)
results_df_pois_6 <- xgboost_tune_cnt(tune_grid)

results_df_pois_6[order(results_df_pois_6$mean_loss_test,decreasing=F),]
save(results_df_pois_6, file=file.path(dir.out, 'XGBoost_Param_cnt_6.RData'))

plot(results_df_pois_6$nrounds,results_df_pois_6$mean_loss_train,type='l',col='blue')
lines(results_df_pois_6$nrounds,results_df_pois_6$mean_loss_test,col='red')

params_cnt <- list(
  nrounds = 100,
  eta = 0.08,
  max_depth = 4,
  gamma = 0.4,
  colsample_bytree = 0.7,
  min_child_weight = 1,
  subsample = 1,
  objective=trunc_poisson_obj
)
