dir.data <- "/home/pgrad2/2448355h/My_PhD_Project/00_Dataset/Urban_Fires"
dir.out <- '/home/pgrad2/2448355h/My_PhD_Project/01_Output/Urban_Fire_Model'


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
# bru_safe_sp(force = TRUE)
load(file.path(dir.data, "burn_area","finaldata_wildfire.RData"))
load(file.path(dir.data, "burn_area","weather_covariates.RData"))


date.start <- as.Date('2012/01/01', format="%Y/%m/%d")
date.end <- as.Date('2020/12/31', format="%Y/%m/%d")

max.day.idx <- as.numeric(as.Date(date.end)-date.start) + 1

elapsed_months <- function(end_date, start_date) {
  ed <- as.POSIXlt(end_date)
  sd <- as.POSIXlt(start_date)
  12 * (ed$year - sd$year) + (ed$mon - sd$mon) + 1
}
elapsed_months(data.new[1,]$open,data.new[1,]$open)

data.new$day.idx <- as.numeric(as.Date(data.new$open)-date.start) + 1
data.new$year.idx <- data.new$year -2011
data.new$month.idx <- as.integer(format(data.new$open,format='%m'))
data.new$time.idx <- elapsed_months(data.new$open, date.start)



loc.data.utm <- st_as_sf(data.new, coords=c('x_utm_new','y_utm_new'), crs='+proj=utm +zone=29 +datum=WGS84 +units=km +no_defs' )

cond <- data.new$month.idx==9 & data.new$length > 24*60
coords <- SpatialPointsDataFrame(data.new[cond,],coords=data.new[cond,c('x_utm_new','y_utm_new')], 
                                 proj4string=CRS('+proj=utm +zone=29 +datum=WGS84 +units=km +no_defs' ))
coords$time.idx <- coords$year.idx
library(rnaturalearth)
map <- ne_countries(type = "countries", country = "Portugal",
                    scale = "medium", returnclass = "sf")
projutm <- "+proj=utm +zone=29 +datum=WGS84 +units=km +no_defs"
map <- st_transform(map, crs = projutm)
mainland_bbox <- st_as_sfc(st_bbox(c(xmin = 0, xmax = 1000, ymin = 4000 , ymax = 4800), crs = st_crs(map)))
map_mainland <- st_intersection(map, mainland_bbox)

ggplot() + geom_sf(data = map_mainland) +
  geom_sf(data = st_as_sf(loc.data.utm)) + coord_sf(datum = projutm)

loc.d <- cbind(st_coordinates(map_mainland)[, 1], st_coordinates(map_mainland)[, 2])



domain.polys <- Polygons(list(Polygon(loc.d)), '0')
domainSP <- SpatialPolygons(list(domain.polys),proj4string=CRS(projutm))



f.get.cov <- function(dataset, cov.name){
  time <- dataset$month.idx
  
  x <- dataset$x_utm_new
  y <- dataset$y_utm_new
  spp <- SpatialPoints(data.frame(x=x,y=y),
                       proj4string=CRS('+proj=utm +zone=29 +datum=WGS84 +units=km +no_defs' ))
  spp <- spTransform(spp, CRS("EPSG:4326"))
  v <- rep(NA, nrow(dataset))
  for (t in unique(time)){
    spdf <- SpatialPixelsDataFrame(points = grid_pixels,
                                   data = data.frame(var=as.vector(get(paste(cov.name,'.month',sep=''))[,,t]),
                                                     time=t))
    proj4string(spdf) <- CRS("EPSG:4326")
    idx <- which(time==t)
    v[idx] <- over(spp[idx,],spdf[,'var'])$var
  }
  return(v)
}


B <- SpatialPolygonsDataFrame(domainSP, data.frame('weight'=1), match.ID = F) 


grd <- SpatialGrid(GridTopology(cellcentre.offset = c(-10, 36),
                                cellsize = c(0.0625,0.0625),
                                cells.dim = c(68, 116)))

# grd <- SpatialGrid(GridTopology(cellcentre.offset = c(-10, 36),
#                                 cellsize = c(0.125,0.125),
#                                 cells.dim = c(34, 58)))

grd_sp <- SpatialPixelsDataFrame(points = grd, data = data.frame(grid.idx = 1:length(grd)),proj4string = CRS("+proj=longlat +datum=WGS84"))

grd_poly <- as(grd_sp, 'SpatialPolygonsDataFrame')
grd_poly <- spTransform(grd_poly, CRS(projutm))
grd_poly$lon.grid <- grd_sp@coords[,1]
grd_poly$lat.grid <- grd_sp@coords[,2]

B2 <- as_Spatial(st_intersection(st_as_sf(grd_poly),st_as_sf(B)))
coord.cent <- st_coordinates(st_centroid(st_as_sf(B2)))
B2$grid.cent.x.utm <- coord.cent[,1]
B2$grid.cent.y.utm <- coord.cent[,2]
B2$E <- area(B2)

B2$grid.idx <- 1:nrow(B2)

ggplot(st_as_sf(B2)) + geom_sf(aes(fill=grid.idx))

data.merge <- B2 |>  
  st_as_sf() |> # cast to sf
  mutate(grid_id = row_number()) |> # create unique ID
  st_join(loc.data.utm) |> # join the species dataset
  group_by(grid_id)


data.fit.ba <- data.merge %>% st_drop_geometry() %>% filter( time.idx<=108, length >=24*60 )%>% 
  group_by(grid.idx,grid_id, year.idx, month.idx, time.idx) %>%
  summarise(area_ha = sum(area_ha), 
            log_ba = log(area_ha),
            y = n(),
            lon.grid = mean(lon.grid),
            lat.grid = mean(lat.grid),
            x.utm = mean(grid.cent.x.utm),
            y.utm = mean(grid.cent.y.utm),
            E = mean(E))


data.fit2 <- do.call(rbind, lapply(1:108, function(x) {B2@data$time.idx = x 
return(B2@data)}))
print(dim(data.fit2))
data.fit2 <- merge(data.fit2, data.fit.ba[,c('grid.idx','time.idx','year.idx', 'month.idx','y','area_ha','log_ba')],
                   by=c('grid.idx','time.idx'),all.x=T)

# data.fit2 <- merge(data.fit2, data.fit.ba1[,c('grid.idx','time.idx.grp','y','area_ha','log_ba')],
#                    by.x=c('grid.idx','time.idx'),by.y = c('grid.idx','time.idx.grp'),all.x=T)

# data.fit2 <- merge(data.fit2, data.fit.ba[,c('grid.idx','year.idx','y','area_ha','log_ba')],
#                    by.x=c('grid.idx','time.idx'),
#                    by.y=c('grid.idx','year.idx'),all.x=T)


print(dim(data.fit2))
data.fit2[is.na(data.fit2$y),'y'] <- 0
summary(data.fit2$y)

data.fit2$month.idx <- (data.fit2$time.idx-1)%%12 + 1
data.fit2$year.idx <- (data.fit2$time.idx-1)%/%12 + 1
data.fit2$E1 <- 1
# 
f.get.cov <- function(dataset, cov.name){
  time <- dataset$time.idx
  
  x <- dataset$grid.cent.x.utm
  y <- dataset$grid.cent.y.utm
  
  
  spp <- SpatialPoints(data.frame(x=x,y=y),
                       proj4string=CRS('+proj=utm +zone=29 +datum=WGS84 +units=km +no_defs' ))
  spp <- spTransform(spp, CRS("EPSG:4326"))
  v <- rep(NA, nrow(dataset))
  for (t in unique(time)){
    spdf <- SpatialPixelsDataFrame(points = grid_pixels,
                                   data = data.frame(var=as.vector(get(paste(cov.name,'.month',sep=''))[,,t]),
                                                     time=t))
    proj4string(spdf) <- CRS("EPSG:4326")
    idx <- which(time==t)
    v[idx] <- over(spp[idx,],spdf[,'var'])$var
  }
  return(v)
}
# 
data.fit2$FWI <- f.get.cov(data.fit2,'fwi')
data.fit2$HVegCov <- f.get.cov(data.fit2,'HVegCov')
data.fit2$HVegLAI <- f.get.cov(data.fit2, 'HVegLAI')
data.fit2$HVegTyp <- as.factor(f.get.cov(data.fit2, 'HVegTyp'))
data.fit2$LVegCov <- f.get.cov(data.fit2,'LVegCov')
data.fit2$LVegLAI <- f.get.cov(data.fit2, 'LVegLAI')
data.fit2$LVegTyp <- as.factor(f.get.cov(data.fit2, 'LVegTyp'))
data.fit2$Pricp <- f.get.cov(data.fit2, 'pricp')
data.fit2$RHumi <- f.get.cov(data.fit2, 'rhumi')
data.fit2$Temp <- f.get.cov(data.fit2, 'temp')
data.fit2$UComp <- f.get.cov(data.fit2, 'u')
data.fit2$VComp <- f.get.cov(data.fit2, 'v')
data.fit2$sqrt_ba <- sqrt(data.fit2$area_ha)
data.fit2[is.na(data.fit2$sqrt_ba), 'sqrt_ba'] <- 0
data.fit2$z <- as.vector((data.fit2$y>0)+0)

covar.names <- c('FWI','HVegCov','HVegLAI','LVegCov','LVegLAI','Pricp','RHumi','Temp',
                 'UComp','VComp')
for (var in covar.names ){
  data.fit2[,var] <- (data.fit2[,var]-mean(data.fit2[,var]))/sd(data.fit2[,var])
}

covar.names <- c(covar.names ,#'LVegTyp','HVegTyp',
                 'lon.grid','lat.grid',
                 'month.idx')

data.rf <- data.fit2[data.fit2$time.idx <= 96, c(covar.names,'grid.idx','time.idx','year.idx','area_ha','z','y')]
# for (var in c('LVegTyp','HVegTyp','month.idx')){
#   data.rf[,var] <- as.factor(data.rf[,var])
# }

library(xgboost)
library(pROC)
data.rf$log_ba <- log(data.rf$area_ha)
data.rf[is.na(data.rf$log_ba),'log_ba'] <- 0




kfold_cv  <- function(data, target.name ,covar.names, params, k) {
  auc_train_list <- numeric(k)
  auc_test_list <- numeric(k)
  
  for (i in 1:k) {
    index <- which(data$year.idx==i)
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
  }
  
  mean_auc_train <- mean(auc_train_list)
  mean_auc_test <- mean(auc_test_list)
  
  return(list(mean_auc_test = mean_auc_test, mean_auc_train = mean_auc_train))
}


##################Tune max_depth, min_child_weight, scale_pos_weight #########
xgboost_tune_z <- function(tune_grid){
  
  k <- 8
  
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
      mean_auc_test = result$mean_auc_test
    ))
  }
  
  return(results_df)
}
# 
# tune_grid <- expand.grid(
#   nrounds = 100,
#   eta = c(0.05,0.06,0.07,0.08),
#   max_depth = 4,
#   gamma = 0,
#   scale_pos_weight = 250,
#   colsample_bytree = 1,
#   min_child_weight = 1,
#   subsample = 1
# )
# 
# 
# results_df1 <- xgboost_tune_z(tune_grid)
# print(results_df1[order(results_df1$mean_auc_test,decreasing=T),])
# save(results_df1, file=file.path(dir.out, 'XGBoost_Param_z_1.RData'))
# 250      0.9571726     0.8766037

# 
# tune_grid <- expand.grid(
#   nrounds = 100,
#   eta = 0.03,
#   max_depth = c(2,3,4,5,6,7),
#   gamma = 0,
#   scale_pos_weight = 250,
#   colsample_bytree = 1,
#   min_child_weight = 1,
#   subsample = 1
# )
# 
# results_df2 <- xgboost_tune_z(tune_grid)
# 
# print(results_df2[order(results_df2$mean_auc_test,decreasing=T),])
# 
# save(results_df2, file=file.path(dir.out, 'XGBoost_Param_z_2.RData'))
# load(file.path(dir.out, 'XGBoost_Param_z_2.RData'))


# 
# tune_grid <- expand.grid(
#   nrounds = 100,
#   eta = 0.06,
#   max_depth = 4,
#   gamma = 0,
#   scale_pos_weight = 250,
#   colsample_bytree = 1,
#   min_child_weight = c(1,3,5,7),
#   subsample = 1
# )
# 
# results_df3 <- xgboost_tune_z(tune_grid)
# 
# print(results_df3[order(results_df3$mean_auc_test,decreasing=T),])
# 
# save(results_df3, file=file.path(dir.out, 'XGBoost_Param_z_3.RData'))


# tune_grid <- expand.grid(
#   nrounds = 100,
#   eta = 0.06,
#   max_depth = 4,
#   gamma = 0,
#   scale_pos_weight = 250,
#   colsample_bytree = c(0.6,0.7,0.8,0.9,1),
#   min_child_weight = 1,
#   subsample = c(0.6,0.7,0.8,0.9,1)
# )
# 
# results_df4 <- xgboost_tune_z(tune_grid)
# 
# print(results_df4[order(results_df4$mean_auc_test,decreasing=T),])
# 
# save(results_df4, file=file.path(dir.out, 'XGBoost_Param_z_4.RData'))


# tune_grid <- expand.grid(
#   nrounds = 100,
#   eta = 0.06,
#   max_depth = 4,
#   gamma = c(0,0.1,0.2,0.3,0.4,0.5),
#   scale_pos_weight = 250,
#   colsample_bytree = 0.7,
#   min_child_weight = 1,
#   subsample = 0.7
# )
# 
# results_df5 <- xgboost_tune_z(tune_grid)
# 
# print(results_df5[order(results_df5$mean_auc_test,decreasing=T),])
# 
# save(results_df5, file=file.path(dir.out, 'XGBoost_Param_z_5.RData'))

# tune_grid <- expand.grid(
#   nrounds = seq(10,100,5),
#   eta = 0.06,
#   max_depth = 4,
#   gamma = 0,
#   scale_pos_weight = 250,
#   colsample_bytree = 0.7,
#   min_child_weight = 1,
#   subsample = 0.7
# )
# 
# results_df6 <- xgboost_tune_z(tune_grid)
# print(results_df6[order(results_df6$mean_auc_test,decreasing=T),])
# save(results_df6, file=file.path(dir.out, 'XGBoost_Param_z_6.RData'))


load(file.path(dir.out, 'XGBoost_Param_z_6.RData'))

plot(results_df6$nrounds,results_df6$mean_auc_train,type='l',col='blue',ylim=c(0.85,0.95))
plot(results_df6$nrounds,results_df6$mean_auc_test,col='red',type='l')

tune_grid_z <- expand.grid(
  nrounds = 70,
  eta = 0.06,
  max_depth = 4,
  gamma = 0,
  scale_pos_weight = 250,
  colsample_bytree = 0.7,
  min_child_weight = 1,
  subsample = 0.7
)



################################## Log BA ###########################
kfold_cv_pos <- function(data, target.name ,covar.names, params, k) {
  loss_train_list <- numeric(k)
  loss_test_list <- numeric(k)
  
  for (i in 1:k) {
    index <- which(data$year.idx==i)
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
    model <- xgb.train(params = params, data = dtrain, nrounds = params$nrounds,nthread=4)
    
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
  
  k <- 7
  
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



tune_grid <- expand.grid(
  nrounds = 100,
  eta = seq(0.01,0.1,0.01),
  max_depth = 2,
  gamma = 0,
  colsample_bytree = 1,
  min_child_weight = 1,
  subsample = 1,
  objective=c('reg:squarederror')
)


results_df_reg <- xgboost_tune_ba(tune_grid)

results_df_reg[order(results_df_reg$mean_loss_test,decreasing=F),]




tune_grid <- expand.grid(
  nrounds = 100,
  eta = 0.09,
  max_depth = c(2,3,4,5,6,7,8),
  gamma = 0,
  colsample_bytree = 1,
  min_child_weight = 1,
  subsample = 1,
  objective=c('reg:squarederror')
)

results_df_reg_1 <- xgboost_tune_ba(tune_grid)

results_df_reg_1[order(results_df_reg_1$mean_loss_test,decreasing=F),]



tune_grid <- expand.grid(
  nrounds = 100,
  eta = 0.09,
  max_depth = 2,
  gamma = 0,
  colsample_bytree = 1,
  min_child_weight = c(1,2,3,4,5,6,7,8,9,10,11,12,13),
  subsample = 1,
  objective=c('reg:squarederror')
)


results_df_reg_2 <- xgboost_tune_ba(tune_grid)

results_df_reg_2[order(results_df_reg_2$mean_loss_test,decreasing=F),]




tune_grid <- expand.grid(
  nrounds = 100,
  eta = 0.09,
  max_depth = 2,
  gamma = 0,
  colsample_bytree = seq(0.4,1,0.1),
  min_child_weight = 1,
  subsample = seq(0.4,1,0.1),
  objective=c('reg:squarederror')
)


results_df_reg_3 <- xgboost_tune_ba(tune_grid)

results_df_reg_3[order(results_df_reg_3$mean_loss_test,decreasing=F),]




tune_grid <- expand.grid(
  nrounds = 100,
  eta = 0.09,
  max_depth = 2,
  gamma = seq(0,1,0.1),
  colsample_bytree = 0.8,
  min_child_weight = 1,
  subsample = 0.9,
  objective=c('reg:squarederror')
)


results_df_reg_4 <- xgboost_tune_ba(tune_grid)

results_df_reg_4[order(results_df_reg_4$mean_loss_test,decreasing=F),]


tune_grid <- expand.grid(
  nrounds = seq(10,100,5),
  eta = 0.09,
  max_depth = 2,
  gamma = 0.4,
  colsample_bytree = 0.8,
  min_child_weight = 1,
  subsample = 0.9,
  objective=c('reg:squarederror')
)


results_df_reg_5 <- xgboost_tune_ba(tune_grid)
results_df_reg_5

plot(results_df_reg_5$nrounds,results_df_reg_5$mean_loss_train,type='l',col='blue')
lines(results_df_reg_5$nrounds,results_df_reg_5$mean_loss_test,col='red')

tune_grid <- expand.grid(
  nrounds = 60,
  eta = 0.09,
  max_depth = 2,
  gamma = 0.4,
  colsample_bytree = 0.8,
  min_child_weight = 1,
  subsample = 0.9,
  objective=c('reg:squarederror')
)


################################### Count ###############################



xgboost_tune_cnt <- function(tune_grid){
  
  k <- 7
  
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
    
    result <- kfold_cv_pos(data.rf, 'y', covar.names,  params, k)
    
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
  eta = seq(0.01,0.1,0.01),
  max_depth = 3,
  gamma = 0,
  colsample_bytree =1,
  min_child_weight = 1,
  subsample = 1,
  objective=c('count:poisson')
)


results_df_pois_1 <- xgboost_tune_cnt(tune_grid)

results_df_pois_1[order(results_df_pois_1$mean_loss_test,decreasing=F),]


tune_grid <- expand.grid(
  nrounds = 100,
  eta = 0.07,
  max_depth = c(2,3,4,5,6),
  gamma = 0,
  colsample_bytree =1,
  min_child_weight = 1,
  subsample = 1,
  objective=c('count:poisson')
)


results_df_pois_2 <- xgboost_tune_cnt(tune_grid)

results_df_pois_2[order(results_df_pois_2$mean_loss_test,decreasing=F),]


tune_grid <- expand.grid(
  nrounds = 100,
  eta = 0.07,
  max_depth = 3,
  gamma = 0,
  colsample_bytree = 1,
  min_child_weight = c(1,2,3,4,5,6,7,8,9),
  subsample = 1,
  objective=c('count:poisson')
)
results_df_pois_3 <- xgboost_tune_cnt(tune_grid)

results_df_pois_3[order(results_df_pois_3$mean_loss_test,decreasing=F),]


tune_grid <- expand.grid(
  nrounds = 100,
  eta = 0.07,
  max_depth = 3,
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
  eta = 0.07,
  max_depth = 3,
  gamma = seq(0,1,0.1),
  colsample_bytree = 0.7,
  min_child_weight = 1,
  subsample = 0.9,
  objective=c('count:poisson')
)
results_df_pois_5 <- xgboost_tune_cnt(tune_grid)

results_df_pois_5[order(results_df_pois_5$mean_loss_test,decreasing=F),]


tune_grid <- expand.grid(
  nrounds = seq(10,100,5),
  eta = 0.07,
  max_depth = 3,
  gamma = 0.7,
  colsample_bytree = 0.7,
  min_child_weight = 1,
  subsample = 0.9,
  objective=c('count:poisson')
)
results_df_pois_6 <- xgboost_tune_cnt(tune_grid)

results_df_pois_6[order(results_df_pois_6$mean_loss_test,decreasing=F),]


plot(results_df_pois_6$nrounds,results_df_pois_6$mean_loss_train,type='l',col='blue')
lines(results_df_pois_6$nrounds,results_df_pois_6$mean_loss_test,col='red')

tune_grid <- expand.grid(
  nrounds = 100,
  eta = 0.07,
  max_depth = 3,
  gamma = 0.7,
  colsample_bytree = 0.7,
  min_child_weight = 1,
  subsample = 0.9,
  objective=c('count:poisson')
)
