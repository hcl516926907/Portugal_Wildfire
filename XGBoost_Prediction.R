dir.data <- "/home/pgrad2/2448355h/My_PhD_Project/00_Dataset/Porgual_Wildfire"
dir.out <- '/home/pgrad2/2448355h/My_PhD_Project/01_Output/Portugal_Wildfire'

library(fastDummies)
library(dplyr)

load(file.path(dir.data, "Data_For_Fitting_Council.RData"))
data.fit <- data.fit.council


####################################################################
# INLA
####################################################################
district.factor <- data.fit$NAME_1
data.fit$NAME_1 <- droplevels(district.factor[district.factor != "Açores" & district.factor != "Madeira"])

data.fit$z <- as.vector((data.fit$y>0)+0)
data.fit$HVegTyp <- as.factor(data.fit$HVegTyp)
# levels(data.fit$HVegTyp) <-  c('NoCover','Evergreen_broadleaf_trees','Mixed_forest','Interrupted_forest')
data.fit$LVegTyp <- as.factor(data.fit$LVegTyp)
# levels(data.fit$LVegTyp) <- c('NoCover','Crops','Tall_grass','Semidesert','Evergreen_shrubs')
# data.fit$month <- as.factor(data.fit$month)
data.fit <- dummy_cols(data.fit, 
                       select_columns = c('HVegTyp','LVegTyp','NAME_1'),remove_first_dummy=TRUE)

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


data.rf <- data.fit[data.fit$year <= 2022, c(covar.names,'grid.idx','time.idx','area_ha','z','y')]
data.rf$log_ba <- log(data.rf$area_ha)
data.rf[is.na(data.rf$log_ba),'log_ba'] <- 0
data.rf[is.na(data.rf$area_ha),'area_ha'] <- 0

library(xgboost)
# library(caret)
library(pROC)


kfold_cv_pred  <- function(data, target.name ,covar.names, params, k) {
  auc_train_list <- numeric(k)
  auc_test_list <- numeric(k)
  
  for (i in 1:k) {
    # select.year <- 2011:2013 + 3*(i-1)
    select.year <- 2011 + (i-1)
    index <- which(data$year %in% select.year )
    train_data <- data[-index, covar.names ]
    train_target <- data[-index, target.name]
    test_data <- data[index, covar.names]
    test_target <- data[index, target.name ]
    
    dtrain <- xgb.DMatrix(data = as.matrix(train_data), label = train_target)
    dtest <- xgb.DMatrix(data = as.matrix(test_data), label = test_target)
    
    set.seed(1234)
    model <- xgb.train(params = params, data = dtrain, nrounds = params$nrounds)
    
    train_pred <- predict(model, dtrain)
    test_pred <- predict(model, dtest)
    
    data[index,'score'] <- test_pred
  }
  
  
  return(data$score)
}


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

# params_z <- list(
#   nrounds = 100,
#   eta = 0.1,
#   max_depth = 6,
#   gamma = 1.3,
#   scale_pos_weight = 6.87,
#   colsample_bytree = 0.6,
#   min_child_weight = 16,
#   subsample = 0.9,
#   objective = "binary:logistic",
#   verbosity = 0
# )


data.rf$score_z <- kfold_cv_pred(data.rf, 'z', covar.names ,params_z, 12)


TrainSet1 <- xgb.DMatrix(data=as.matrix(data.rf[, c(covar.names)]),label=data.rf[, c('z')])
set.seed(1234)
model.z <- xgb.train(params = params_z, data = TrainSet1, nrounds = params_z$nrounds)

xgb.importance(covar.names,model.z)

data.fit[data.fit$year <= 2022, 'score_z'] <- data.rf$score_z
data.fit[data.fit$year > 2022, 'score_z'] <- predict(model.z, as.matrix(data.fit[data.fit$year > 2022,c(covar.names)]))
data.fit %>% group_by(year)%>% summarize(AUC=as.numeric(auc(z,score_z,quiet=TRUE)))


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

# params_ba <- list(
#   nrounds = 95,
#   eta = 0.1,
#   max_depth =  2,
#   gamma = 9,
#   colsample_bytree = 0.9,
#   min_child_weight = 7,
#   subsample = 0.6,
#   objective='reg:squarederror'
# )

kfold_cv_pred_pos <- function(data, target.name ,covar.names, params, k) {
  loss_train_list <- numeric(k)
  loss_test_list <- numeric(k)
  for (i in 1:k) {
    select.year <- 2011 + (i-1)
    index <- which(data$year %in% select.year )
    train_data <- data[-index, covar.names ]
    train_target <- data[-index, target.name]
    pos.ind.train <- train_target>0
    train_data <- train_data[pos.ind.train,]
    train_target <- train_target[pos.ind.train]
    # print(table(data[-index,'grid.idx']))
    
    test_data <- data[index, covar.names]
    test_target <- data[index, target.name ]
    
    
    dtrain <- xgb.DMatrix(data = as.matrix(train_data), label = train_target)
    dtest <- xgb.DMatrix(data = as.matrix(test_data), label = test_target)
    
    set.seed(1234)
    model <- xgb.train(params = params, data = dtrain, nrounds = params$nrounds)
    
    train_pred <- predict(model, dtrain)
    test_pred <- predict(model, dtest)
    
    data[index, 'score'] <- test_pred
  }
  
  
  return(data$score)
}


data.rf$score_ba <- kfold_cv_pred_pos(data.rf,'log_ba',covar.names, params_ba, 12)

TrainSet2 <- xgb.DMatrix(data=as.matrix(data.rf[data.rf$y>0, c(covar.names)]),label=data.rf[data.rf$y>0, c('log_ba')])
#
set.seed(1234)
model.ba <- xgb.train(params = params_ba, data = TrainSet2, nrounds = params_ba$nrounds)


data.fit[data.fit$year <= 2022, 'score_ba'] <- data.rf$score_ba
data.fit[data.fit$year > 2022, 'score_ba'] <- predict(model.ba, as.matrix(data.fit[data.fit$year > 2022,c(covar.names)]))



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

# params_cnt <- list(
#   nrounds = 60,
#   eta = 0.08,
#   max_depth = 5,
#   gamma = 0,
#   colsample_bytree = 0.7,
#   min_child_weight = 1,
#   subsample = 1,
#   objective=trunc_poisson_obj
# )
data.rf$score_cnt <- kfold_cv_pred_pos(data.rf,'y',covar.names, params_cnt, 12)


TrainSet3 <- xgb.DMatrix(data=as.matrix(data.rf[data.rf$y>0, c(covar.names)]),label=data.rf[data.rf$y>0, c('y')])

set.seed(1234)
model.cnt <- xgb.train(params = params_cnt, data = TrainSet3, nrounds = params_cnt$nrounds)


data.fit[data.fit$year <= 2022, 'score_cnt'] <- exp(data.rf$score_cnt)
data.fit[data.fit$year > 2022 , 'score_cnt'] <- exp(predict(model.cnt, as.matrix(data.fit[data.fit$year > 2022 ,c(covar.names)])))

# 
# 
data.fit[is.na(data.fit$log_ba),'log_ba'] <- 0




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

data.fit %>% group_by(year)%>% summarize(AUC=as.numeric(auc(z,score_z,quiet=TRUE)))

data.fit[data.fit$y>0,]%>%
  group_by(year)%>% summarize(loss=weighted_loss(y,score_cnt,'CNT'))


data.fit[data.fit$y>0,]%>%
  group_by(year)%>% summarize(loss=weighted_loss(area_ha,exp(score_ba),'BA'))


save(data.fit, file=file.path(dir.out, 'XGBoost_Score_Council_2.RData'))


save(data.rf, covar.names,model.z, model.ba, model.cnt, 
     file=file.path(dir.out, 'XGBoost_Models.RData'))
