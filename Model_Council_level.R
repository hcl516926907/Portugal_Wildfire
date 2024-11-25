dir.data <- "/home/pgrad2/2448355h/My_PhD_Project/00_Dataset/Porgual_Wildfire"
dir.out <- '/home/pgrad2/2448355h/My_PhD_Project/01_Output/Portugal_Wildfire'


library(INLA)

inla.setOption(scale.model.default = TRUE)

library(mgcv)
library(ggplot2)
library(rgdal)
library(fmesher)
library(dplyr)
library(RColorBrewer)
library(terra)
# library(rgeos)
library(spdep)
library(raster)
# library(scoringRules)
library(fastDummies)
# bru_safe_sp(force = TRUE)
library(sf)

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
data.fit$month <- as.factor(data.fit$month)
data.fit <- dummy_cols(data.fit, 
                       select_columns = c('HVegTyp','LVegTyp','month','NAME_1'),remove_first_dummy=TRUE)

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


data.rf <- data.fit[data.fit$year <= 2022, c(covar.names,'grid.idx','time.idx','area_ha','z','y')]

#
library(xgboost)
# library(caret)
library(pROC)
data.rf$log_ba <- log(data.rf$area_ha)
data.rf[is.na(data.rf$log_ba),'log_ba'] <- 0

kfold_cv_pred  <- function(data, target.name ,covar.names, params, k) {
  auc_train_list <- numeric(k)
  auc_test_list <- numeric(k)
  
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
    model <- xgb.train(params = params, data = dtrain, nrounds = params$nrounds)
    
    train_pred <- predict(model, dtrain)
    test_pred <- predict(model, dtest)
    
    data[index,'score'] <- test_pred
  }
  
  
  return(data$score)
}


params_z <- list(
  nrounds = 100,
  eta = 0.09,
  max_depth = 7,
  gamma = 0.7,
  scale_pos_weight = 2,
  colsample_bytree = 0.4,
  min_child_weight = 1,
  subsample = 1,
  objective = "binary:logistic",
  eval_metric = "auc",
  verbosity = 0
)

#
data.rf$score_z <- kfold_cv_pred(data.rf, 'z', covar.names ,params_z, 4)


TrainSet1 <- xgb.DMatrix(data=as.matrix(data.rf[, c(covar.names)]),label=data.rf[, c('z')])
set.seed(1234)
model.z <- xgb.train(params = params_z, data = TrainSet1, nrounds = params_z$nrounds)


data.fit[data.fit$year <= 2022, 'score_z'] <- data.rf$score_z
data.fit[data.fit$year > 2022, 'score_z'] <- predict(model.z, as.matrix(data.fit[data.fit$year > 2022,c(covar.names)]))

#
#
# # with coordiantes
params_ba <- list(
  nrounds = 100,
  eta = 0.05,
  max_depth = 3,
  gamma = 0.2,
  colsample_bytree = 0.4,
  min_child_weight = 5,
  subsample = 0.5,
  objective=c('reg:squarederror')
)
# 
# 
kfold_cv_pred_pos <- function(data, target.name ,covar.names, params, k) {
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


data.rf$score_ba <- kfold_cv_pred_pos(data.rf,'log_ba',covar.names, params_ba, 4)

TrainSet2 <- xgb.DMatrix(data=as.matrix(data.rf[data.rf$y>0, c(covar.names)]),label=data.rf[data.rf$y>0, c('log_ba')])
#
set.seed(1234)
model.ba <- xgb.train(params = params_ba, data = TrainSet2, nrounds = params_ba$nrounds)


data.fit[data.fit$year <= 2022, 'score_ba'] <- data.rf$score_ba
data.fit[data.fit$year > 2022, 'score_ba'] <- predict(model.ba, as.matrix(data.fit[data.fit$year > 2022,c(covar.names)]))

# 
# 
params_cnt <- list(
  nrounds = 95,
  eta = 0.09,
  max_depth = 6,
  gamma = 0,
  colsample_bytree = 0.7,
  min_child_weight = 1,
  subsample = 0.5,
  objective=c('count:poisson')
)

data.rf$score_cnt <- kfold_cv_pred_pos(data.rf,'y',covar.names, params_ba, 4)


TrainSet3 <- xgb.DMatrix(data=as.matrix(data.rf[data.rf$y>0, c(covar.names)]),label=data.rf[data.rf$y>0, c('y')])

set.seed(1234)
model.cnt <- xgb.train(params = params_cnt, data = TrainSet3, nrounds = params_cnt$nrounds)


data.fit[data.fit$year <= 2022, 'score_cnt'] <- data.rf$score_cnt
data.fit[data.fit$year > 2022 , 'score_cnt'] <- predict(model.cnt, as.matrix(data.fit[data.fit$year > 2022 ,c(covar.names)]))

# 
# 
data.fit[is.na(data.fit$log_ba),'log_ba'] <- 0

data.fit %>% group_by(year)%>% summarize(AUC=as.numeric(auc(z,score_z,quiet=TRUE)))

data.fit[data.fit$y>0,]%>%mutate(square_loss=(log_ba-score_ba)^2,sum_loss=sum(square_loss))%>%
  group_by(year)%>% summarize(avg_loss=mean(square_loss))


poisson_loss <- function(y,lambda){
  return(mean(-y*log(lambda)+lambda))
}


data.fit[data.fit$y>0,]%>%group_by(year)%>% summarize(avg_loss=poisson_loss(y,score_cnt))



dist<-shapefile(file.path(dir.data,'shapefile', "distritos.shp"))
dist$ID_0 <- as.factor(iconv(as.character(dist$ID_0), "UTF-8"))
dist$ISO <- as.factor(iconv(as.character(dist$ISO), "UTF-8"))
dist$NAME_0 <- as.factor(iconv(as.character(dist$NAME_0), "UTF-8"))
dist$ID_1 <- as.factor(iconv(as.character(dist$ID_1), "UTF-8"))
dist$NAME_1 <- as.factor(iconv(as.character(dist$NAME_1), "UTF-8"))
dist$HASC_1 <- as.factor(iconv(as.character(dist$HASC_1),  "UTF-8"))
dist$CCN_1<- as.factor(iconv(as.character(dist$CCN_1),  "UTF-8"))
dist$CCA_1 <- as.factor(iconv(as.character(dist$CCA_1), "UTF-8"))
dist$TYPE_1 <- as.factor(iconv(as.character(dist$TYPE_1), "UTF-8"))
dist$ENGTYPE_1 <- as.factor(iconv(as.character(dist$ENGTYPE_1), "UTF-8"))
dist$NL_NAME_1 <- as.factor(iconv(as.character(dist$NL_NAME_1), "UTF-8"))
dist$VARNAME_1 <- as.factor(iconv(as.character(dist$VARNAME_1), "UTF-8"))
dist=dist[dist$NAME_1!="Açores",]
dist=dist[dist$NAME_1!="Madeira",]
sf_districts <- st_as_sf(dist)

# map <- conc[,'NAME_2']
conc<-shapefile(file.path(dir.data, "shapefile","concelhos.shp"))
conc$ID_0 <- as.factor(iconv(as.character(conc$ID_0),  "UTF-8"))
conc$ISO <- as.factor(iconv(as.character(conc$ISO), "UTF-8"))
conc$NAME_0 <- as.factor(iconv(as.character(conc$NAME_0), "UTF-8"))
conc$NAME_1 <- as.factor(iconv(as.character(conc$NAME_1), "UTF-8"))
conc$ID_2 <- as.factor(iconv(as.character(conc$ID_2), "UTF-8"))
conc$NAME_2 <- as.factor(iconv(as.character(conc$NAME_2), "UTF-8"))
conc$HASC_2 <- as.factor(iconv(as.character(conc$HASC_2),  "UTF-8"))
conc$CCN_2 <- as.factor(iconv(as.character(conc$CCN_2),  "UTF-8"))
conc$CCA_2 <- as.factor(iconv(as.character(conc$CCA_2),  "UTF-8"))
conc$TYPE_2 <- as.factor(iconv(as.character(conc$TYPE_2),"UTF-8"))
conc$ENGTYPE_2 <- as.factor(iconv(as.character(conc$ENGTYPE_2), "UTF-8"))
conc$NL_NAME_2 <- as.factor(iconv(as.character(conc$NL_NAME_2), "UTF-8"))
conc$VARNAME_2 <- as.factor(iconv(as.character(conc$VARNAME_2), "UTF-8"))
conc=conc[conc$NAME_1!="Azores",]
conc=conc[conc$NAME_1!="Madeira",]
conc$NAME_1<-as.factor(droplevels(conc$NAME_1))
conc$NAME_2<-as.factor(droplevels(conc$NAME_2))
sf_conc <- st_as_sf(conc)

map <- sf_conc[,'NAME_2']
rownames(map) <- sf_conc$NAME_2

nb <- poly2nb(map)
map$grid.idx <- 1:nrow(map)
centroids <- st_centroid(map)
centroid_coords <- st_coordinates(centroids)
map$lon.grid <- centroid_coords[,1]
map$lat.grid <- centroid_coords[,2]

ggplot(map) + geom_sf()+
  geom_text(aes(lon.grid, lat.grid, label = grid.idx), size=2,color='red') 


nb2INLA("map.adj", nb)
g.council <- inla.read.graph(filename = "map.adj")

# 
# 
map.district <- sf_districts[,'NAME_1']
rownames(map.district) <- sf_districts$NAME_1

nb.district <- poly2nb(map.district)
map.district$grid.idx.district <-  1:nrow(map.district)
centroids.district <- st_centroid(map.district)
centroid_coords.district <- st_coordinates(centroids.district)
map.district$lon.grid <- centroid_coords.district[,1]
map.district$lat.grid <- centroid_coords.district[,2]

ggplot(map.district) + geom_sf()+
  geom_text(aes(lon.grid, lat.grid, label = grid.idx.district), size=2,color='red')

nb2INLA("map.district.adj", nb.district)
g.district <- inla.read.graph(filename = "map.district.adj")



data.fit.reshape <- reshape(data.fit[,c('grid.idx','time.idx','y')],
                            timevar = "time.idx",
                            idvar = "grid.idx",
                            direction = "wide"
)

data.fit$month <- (data.fit$time.idx-1)%%12 + 1
data.fit$year <- (data.fit$time.idx-1)%/%12 + 2011


# B2.merge <- merge(B1, data.fit.reshape, by = "grid.idx")
# save(B2.merge, file=file.path(dir.out, 'grid_cell_map_0.1.RData'))
# 


###################merge district predictions###################
# load(file=file.path(dir.out, 'districts_prediction.RData'))
# 
# data.fit <- merge(data.fit, data.fit.dist.pred, by=c('NAME_1','time.idx'))
# 
# 
# 
data.fit <- merge(data.fit,st_drop_geometry(map.district[,c('NAME_1','grid.idx.district')]),by='NAME_1')


data.fit <- data.fit[order(data.fit$time.idx,data.fit$grid.idx),]
 
z <- as.vector((data.fit$y>0)+0)
log.ba <- as.vector(ifelse(data.fit$y>0, data.fit$log_ba, NA))
cnt <- as.vector(ifelse(data.fit$y>0, data.fit$y, NA))


#prepare for prediction
z[which(data.fit$year>2022)] <- NA
log.ba[which(data.fit$year>2022)] <- NA
cnt[which(data.fit$year>2022)] <- NA

# # 
# knots_cnt <- quantile(data.fit[data.fit$y>0,'score_cnt'],c(0.05,0.2,0.35,0.5,0.65,0.8,0.95))
# mesh_score_1 <- inla.mesh.1d(knots_cnt, boundary=c('free','free'))
# A1 <- inla.spde.make.A(mesh_score_1, loc=data.fit$score_cnt)
# spde_score_1 <-  inla.spde2.pcmatern(mesh_score_1,
#                                      prior.range = c(1, 0.1),
#                                      prior.sigma = c(1, 0.05))
# 
# 
# spde_score_1.idx <- inla.spde.make.index("score_1", n.spde = spde_score_1$n.spde)
# 
# knots_z <- quantile(data.fit[,'score_z'],c(0.05,seq(0.1,0.9,0.1),0.95))
# mesh_score_2 <- inla.mesh.1d(knots_z,boundary=c('free','free'))
# 
# A2 <- inla.spde.make.A(mesh_score_2, loc=data.fit$score_z)
# spde_score_2 <-  inla.spde2.pcmatern(mesh_score_2,
#                                      prior.range = c(0.1, 0.05),
#                                      prior.sigma = c(1, 0.05))
# spde_score_2.idx <- inla.spde.make.index("score_2", n.spde = spde_score_2$n.spde)
# 
# knots_ba <- quantile(data.fit[data.fit$y>0,'score_ba'],c(0.05,0.2,0.35,0.5,0.65,0.8,0.95))
# mesh_score_3 <- inla.mesh.1d(knots_ba, boundary=c('free','free'))
# 
# A3 <- inla.spde.make.A(mesh_score_3, loc=data.fit$score_ba)
# spde_score_3 <-  inla.spde2.pcmatern(mesh_score_3,
#                                      prior.range = c(10, 0.2),
#                                      prior.sigma = c(1, 0.1))
# spde_score_3.idx <- inla.spde.make.index("score_3", n.spde = spde_score_3$n.spde)




#####district score#########
# knots_cnt.dist <- quantile(data.fit[data.fit$y>0,'score_cnt.dist'],c(0.05,0.2,0.35,0.5,0.65,0.8,0.95))
# mesh_score_1.dist <- inla.mesh.1d(knots_cnt.dist, boundary=c('free','free'))
# A1.dist <- inla.spde.make.A(mesh_score_1.dist, loc=data.fit$score_cnt.dist)
# spde_score_1.dist <-  inla.spde2.pcmatern(mesh_score_1.dist,
#                                      prior.range = c(1, 0.1),
#                                      prior.sigma = c(1, 0.05))
# spde_score_1.dist.idx <- inla.spde.make.index("score_1.dist", n.spde = spde_score_1.dist$n.spde)
# 
# 
# knots_z.dist <- quantile(data.fit[,'score_z.dist'],c(0.05,seq(0.1,0.9,0.1),0.95))
# mesh_score_2.dist <- inla.mesh.1d(knots_z.dist,boundary=c('free','free'))
# A2.dist <- inla.spde.make.A(mesh_score_2.dist, loc=data.fit$score_z.dist)
# spde_score_2.dist <-  inla.spde2.pcmatern(mesh_score_2.dist,
#                                      prior.range = c(0.1, 0.05),
#                                      prior.sigma = c(1, 0.05))
# spde_score_2.dist.idx <- inla.spde.make.index("score_2.dist", n.spde = spde_score_2.dist$n.spde)
# 
# 
# knots_ba.dist <- quantile(data.fit[data.fit$y>0,'score_ba.dist'],c(0.05,0.2,0.35,0.5,0.65,0.8,0.95))
# mesh_score_3.dist <- inla.mesh.1d(knots_ba.dist, boundary=c('free','free'))
# A3.dist <- inla.spde.make.A(mesh_score_3.dist, loc=data.fit$score_ba.dist)
# spde_score_3.dist <-  inla.spde2.pcmatern(mesh_score_3.dist,
#                                      prior.range = c(10, 0.2),
#                                      prior.sigma = c(1, 0.1))
# spde_score_3.dist.idx <- inla.spde.make.index("score_3.dist", n.spde = spde_score_3.dist$n.spde)
# 

# 
# cnt.stack <- inla.stack(
#   data= list(Y.log=cbind(cnt,NA,NA)),
#   A <- list(1,1,A1, 1 , 1,1,1,1,1,1,1,1,1,1,1,1 ,1,1,1,1,1),
#   effect = list(Intercept1=rep(1,nrow(data.fit)), month.idx.cnt=data.fit$month, score_1=spde_score_1.idx,
#                 grid.idx.cnt=data.fit$grid.idx, 
#                 grid.idx.dist.cnt = data.fit$grid.idx.district,
#                 year.af.2017.cnt = as.integer(data.fit$year>2017),
#                 month_2.cnt = data.fit$month_2, month_3.cnt = data.fit$month_3,
#                 month_4.cnt = data.fit$month_4, month_5.cnt = data.fit$month_5,
#                 month_6.cnt = data.fit$month_6, month_7.cnt = data.fit$month_7,
#                 month_8.cnt = data.fit$month_8, month_9.cnt = data.fit$month_9,
#                 month_10.cnt = data.fit$month_10, month_11.cnt = data.fit$month_11,
#                 month_12.cnt = data.fit$month_12,
#                 year.idx.cnt=data.fit$year-2010,
#                 year.af.2017.cnt = as.numeric(data.fit$year>2017),
#                 eps.time.cnt = data.fit$time.idx,
#                 eps.dist.cnt = unclass(data.fit$NAME_1)),
#   tag='cnt'
# )
# 
# z.stack <- inla.stack(
#   data= list(Y.log=cbind(NA,z,NA)),
#   A <- list(1,1,A2, 1 , 1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1),
#   effect = list(Intercept2=rep(1,nrow(data.fit)), month.idx.z=data.fit$month, score_2=spde_score_2.idx,
#                 grid.idx.z=data.fit$grid.idx, 
#                 grid.idx.dist.z = data.fit$grid.idx.district,
#                 year.af.2017.z = as.integer(data.fit$year>2017),
#                 month_2.z = data.fit$month_2, month_3.z = data.fit$month_3,
#                 month_4.z = data.fit$month_4, month_5.z = data.fit$month_5,
#                 month_6.z = data.fit$month_6, month_7.z = data.fit$month_7,
#                 month_8.z = data.fit$month_8, month_9.z = data.fit$month_9,
#                 month_10.z = data.fit$month_10, month_11.z = data.fit$month_11,
#                 month_12.z = data.fit$month_12,
#                 year.idx.z=data.fit$year-2010,
#                 year.af.2017.z = as.numeric(data.fit$year>2017),
#                 eps.time.z = data.fit$time.idx,
#                 eps.dist.z = unclass(data.fit$NAME_1)),
#   tag='z'
# )
# 
# 
# ba.stack <- inla.stack(
#   data= list(Y.log=cbind(NA,NA,log.ba)),
#   A <- list(1,1, A3, 1 , 1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1),
#   effect = list(Intercept3=rep(1,nrow(data.fit)), month.idx.ba=data.fit$month, score_3=spde_score_3.idx,
#                 grid.idx.ba = data.fit$grid.idx,
#                 grid.idx.dist.ba = data.fit$grid.idx.district,
#                 year.af.2017.ba = as.integer(data.fit$year>2017),
#                 month_2.ba = data.fit$month_2, month_3.ba = data.fit$month_3,
#                 month_4.ba = data.fit$month_4, month_5.ba = data.fit$month_5,
#                 month_6.ba = data.fit$month_6, month_7.ba = data.fit$month_7,
#                 month_8.ba = data.fit$month_8, month_9.ba = data.fit$month_9,
#                 month_10.ba = data.fit$month_10, month_11.ba = data.fit$month_11,
#                 month_12.ba = data.fit$month_12,
#                 year.idx.ba=data.fit$year-2010,
#                 year.af.2017.z = as.numeric(data.fit$year>2017),
#                 eps.time.ba = data.fit$time.idx,
#                 eps.dist.ba = unclass(data.fit$NAME_1)),
#   tag='ba'
# )
# 
# 
# all.stack <- inla.stack(cnt.stack, z.stack, ba.stack )


n1 <- dim(data.fit)[1]
nothing <- rep(NA, n1)

cntNA <- as.vector(c(cnt,nothing,nothing))
zNA <- as.vector(c(nothing, z, nothing))
baNA.log = as.vector(c(nothing, nothing, log.ba))

Y.log <- matrix(c(cntNA,zNA, baNA.log), ncol=3)

grid.idx.cnt = c(data.fit$grid.idx, nothing, nothing)# fire ignition
grid.idx.z = c(nothing,data.fit$grid.idx, nothing)# BA ind
grid.idx.ba = c(nothing, nothing, data.fit$grid.idx)# BA

grid.idx.dist.cnt = c(data.fit$grid.idx.district, nothing, nothing)# fire ignition
grid.idx.dist.z = c(nothing,data.fit$grid.idx.district, nothing)# BA ind
grid.idx.dist.ba = c(nothing, nothing, data.fit$grid.idx.district)# BA


time.idx.cnt <- c(data.fit$time.idx, nothing, nothing)# fire ignition
time.idx.z <- c(nothing, data.fit$time.idx, nothing)# BA ind
time.idx.ba  <-  c(nothing, nothing, data.fit$time.idx)# BA

intercept.cnt <- c(rep(1,n1),nothing, nothing)
intercept.z <- c(nothing, rep(1,n1), nothing)
intercept.ba <- c(nothing, nothing, rep(1,n1))

score_cnt_grp <- inla.group(data.fit$score_cnt, n = 20, method = "quantile")
score.cnt <- c(score_cnt_grp,nothing, nothing)

score_z_grp <- inla.group(data.fit$score_z, n = 20, method = "quantile")
score.z <- c(nothing, score_z_grp, nothing)

score_ba_grp <- inla.group(data.fit$score_ba, n = 20, method = "quantile")
score.ba <- c(nothing, nothing, score_ba_grp)


month.idx.cnt <- c(data.fit$month, nothing, nothing)# fire ignition
month.idx.z <- c(nothing, data.fit$month, nothing)# BA ind
month.idx.ba  <-  c(nothing, nothing, data.fit$month)# BA

year.idx <- data.fit$year - 2010
year.idx.cnt <- c(year.idx, nothing, nothing)# fire ignition
year.idx.z <- c(nothing, year.idx, nothing)# BA ind
year.idx.ba  <-  c(nothing, nothing,year.idx)# BA

after.2017 <- as.integer(data.fit$year>2017)
year.af.2017.cnt <- c(after.2017, nothing, nothing)
year.af.2017.z <- c(nothing, after.2017, nothing)
year.af.2017.ba <- c(nothing, nothing,  after.2017)



data.list=list(Y.log=Y.log, 

          grid.idx.cnt=grid.idx.cnt, 
          grid.idx.z=grid.idx.z, 
          grid.idx.ba=grid.idx.ba,
          
          grid.idx.dist.cnt = grid.idx.dist.cnt,
          grid.idx.dist.z = grid.idx.dist.z,
          grid.idx.dist.ba = grid.idx.dist.ba,

          time.idx.cnt = time.idx.cnt,
          time.idx.z = time.idx.z,
          time.idx.ba = time.idx.ba,

          intercept.cnt = intercept.cnt,
          intercept.z = intercept.z,
          intercept.ba = intercept.ba,
          
          month.idx.cnt = month.idx.cnt,
          month.idx.z = month.idx.z,
          month.idx.ba = month.idx.ba,
          
          year.idx.cnt = year.idx.cnt,
          year.idx.z = year.idx.z,
          year.idx.ba = year.idx.ba,
          
          year.af.2017.cnt = year.af.2017.cnt,
          year.af.2017.z = year.af.2017.z,
          year.af.2017.ba = year.af.2017.ba,
          
          score.cnt = score.cnt,
          score.z = score.z,
          score.ba = score.ba
          
)





# formula <- Y.log ~  0  +
#   intercept.cnt + year.af.2017.cnt +
#   f(grid.idx.cnt, copy='grid.idx.z',fixed=F)  +
#   f(month.idx.cnt, model='seasonal', season.length=12 )+
#   f(year.idx.cnt, model='rw1')+ f(score.cnt, model='rw1' )+
#   f(grid.idx.dist.cnt, copy="grid.idx.dist.z",fixed=F, group=time.idx.cnt, control.group = list(model='iid'))  +
# 
#   intercept.z + year.af.2017.z +
#   f(grid.idx.z, model='bym2',graph=g.council)  +
#   f(month.idx.z, model='seasonal',season.length=12 )+
#   f(year.idx.z, model='rw1')+  f(score.z, model='rw1' )+
#   f(grid.idx.dist.z, model='bym2', graph=g.district,group=time.idx.z, control.group = list(model='iid')) +
# 
#   intercept.ba + year.af.2017.ba +
#   f(grid.idx.ba, copy='grid.idx.z',fixed=F)  +
#   f(month.idx.ba, model='seasonal',season.length=12 ) +
#   f(year.idx.ba, model='rw1') +  f(score.ba, model='rw1' ) +
#   f(grid.idx.dist.ba, copy="grid.idx.dist.z",fixed=F, group=time.idx.ba, control.group = list(model='iid'))  


formula <- Y.log ~  0  +
  intercept.cnt + year.af.2017.cnt +
  f(grid.idx.cnt, copy='grid.idx.z',fixed=F,group=month.idx.cnt,control.group = list(model='iid'))  +
  f(year.idx.cnt, model='iid')+ f(score.cnt, model='rw1' )+
  f(grid.idx.dist.cnt, copy="grid.idx.dist.z",fixed=F, group=time.idx.cnt, control.group = list(model='iid'))  +
  
  intercept.z + year.af.2017.z +
  f(grid.idx.z, model='bym2',graph=g.council,group=month.idx.z,control.group = list(model='iid'))  +
  f(year.idx.z, model='iid')+  f(score.z, model='rw1' )+
  f(grid.idx.dist.z, model='bym2', graph=g.district,group=time.idx.z, control.group = list(model='iid')) +
  
  intercept.ba + year.af.2017.ba +
  f(grid.idx.ba, copy='grid.idx.z',fixed=F,group=month.idx.ba,control.group = list(model='iid'))  +
  f(year.idx.ba, model='iid') +  f(score.ba, model='rw1' ) +
  f(grid.idx.dist.ba, copy="grid.idx.dist.z",fixed=F, group=time.idx.ba, control.group = list(model='iid'))  





# save.image(file=file.path(dir.out, 'tmp.RData'))
# load(file.path(dir.out,'tmp.RData'))
t1 <- Sys.time()
res <- inla(formula,
               family = c('poisson','binomial', 'weibull'), data = data.list,  Ntrials=1,
               control.predictor = list(compute = TRUE, link=c(rep(1,n1),rep(2,n1),rep(3,n1))),
               verbose=TRUE,
               control.compute=list(config = TRUE),
               control.family = list( list(), list(), list(control.link=list(quantile=0.5))),
               control.fixed = list(expand.factor.strategy = 'inla')
)
t2 <- Sys.time()
print(t2-t1)

summary(res)

save(res, file=file.path(dir.out, 'Model_weibull_spde_0.1.RData'))
