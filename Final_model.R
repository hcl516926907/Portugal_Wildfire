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
                                cellsize = c(0.25,0.25),
                                cells.dim = c(17, 29)))

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

data.rf <- data.fit2[data.fit2$time.idx <= 84, c(covar.names,'y','area_ha','z')]
# for (var in c('LVegTyp','HVegTyp','month.idx')){
#   data.rf[,var] <- as.factor(data.rf[,var])
# }

library(xgboost)
library(caret)
library(pROC)
data.rf$log_ba <- log(data.rf$area_ha)
data.rf[is.na(data.rf$log_ba),'log_ba'] <- 0

kfold_cv_pred  <- function(data, target, params, k) {
  auc_train_list <- numeric(k)
  auc_test_list <- numeric(k)
  
  for (i in 1:k) {
    index <- (1:(192*12))+(i-1)*(192*12)
    train_data <- data[-index, ]
    train_target <- target[-index]
    test_data <- data[index, ]
    test_target <- target[index]
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

##without coordiates
# params_z <- list(
#   eta = 0.02,
#   max_depth = 2,
#   gamma = 0.2,
#   colsample_bytree = 0.8,
#   min_child_weight = 3,
#   subsample = 0.8,
#   nrounds = 100,
#   objective = "binary:logistic",
#   eval_metric = "auc",
#   scale_pos_weight = 27,
#   verbosity = 0
# )

##with coordiates
params_z <- list(
  eta = 0.07,
  max_depth = 4,
  gamma = 0,
  colsample_bytree = 0.9,
  min_child_weight = 3,
  subsample = 0.8,
  nrounds = 35,
  objective = "binary:logistic",
  eval_metric = "auc",
  scale_pos_weight = 27,
  verbosity = 0
)

data.rf$score_z <- kfold_cv_pred(data.rf[, c(covar.names)], data.rf[, 'z'], params_z, 7)


TrainSet1 <- xgb.DMatrix(data=as.matrix(data.rf[, c(covar.names)]),label=data.rf[, c('z')])
set.seed(1234)
model.z <- xgb.train(params = params_z, data = TrainSet1, nrounds = params_z$nrounds)


data.fit2[data.fit2$time.idx <= 84, 'score_z'] <- data.rf$score_z
data.fit2[data.fit2$time.idx > 84, 'score_z'] <- predict(model.z, as.matrix(data.fit2[data.fit2$time.idx > 84,c(covar.names)]))


# params_ba <- list(
#   eta = 0.05,
#   max_depth = 3,
#   gamma = 1,
#   colsample_bytree = 0.7,
#   min_child_weight = 5,
#   subsample = 0.7,
#   nrounds = 150,
#   objective = 'reg:tweedie',
#   verbosity = 0
# )
# data.rf$score_ba <- kfold_cv_pred(data.rf[, c(covar.names,'score_z')], data.rf[, 'log_ba'], params_ba, 7)
# 
# 
# TrainSet2 <- xgb.DMatrix(data=as.matrix(data.rf[, c(covar.names,'score_z')]),label=data.rf[, c('log_ba')])
# 
# model.ba <- xgb.train(params = params_ba, data = TrainSet2, nrounds = params_ba$nrounds)

#without coordiantes
# params_ba <- list(
#   eta = 0.1,
#   max_depth = 2,
#   gamma = 1,
#   colsample_bytree = 0.7,
#   min_child_weight = 5,
#   subsample = 0.8,
#   nrounds = 50,
#   objective = 'reg:squarederror',
#   verbosity = 0
# )

# with coordiantes
params_ba <- list(
  nrounds = 60,
  eta = 0.05,
  max_depth = 2,
  gamma = 0.7,
  colsample_bytree = 0.9,
  min_child_weight = 8,
  subsample = 0.5,
  objective = 'reg:squarederror',
  verbosity = 0
)


kfold_cv_pred_1  <- function(data, target, params, k) {
  auc_train_list <- numeric(k)
  auc_test_list <- numeric(k)
  data$score <- NA
  for (i in 1:k) {
    index <- (1:(192*12))+(i-1)*(192*12)
    
    train_data <- data[-index, ]
    train_target <- target[-index]
    pos.ind.train <- train_target>0
    train_data <- train_data[pos.ind.train,]
    train_target <- train_target[pos.ind.train]
    
    test_data <- data[index, ]
    test_target <- target[index]
    pos.ind.test <- test_target>0
    test_data <- test_data[pos.ind.test,]
    test_target <- test_target[pos.ind.test]
    
    
    dtrain <- xgb.DMatrix(data = as.matrix(train_data), label = train_target)
    dtest <- xgb.DMatrix(data = as.matrix(test_data), label = test_target)
    
    set.seed(1234)
    model <- xgb.train(params = params, data = dtrain, nrounds = params$nrounds)
    
    train_pred <- predict(model, dtrain)
    test_pred <- predict(model, dtest)
    
    data[index,][pos.ind.test,'score'] <- test_pred
  }
  
  
  return(data$score)
}


data.rf$score_ba <- kfold_cv_pred_1(data.rf[, c(covar.names)], data.rf[, 'log_ba'], params_ba, 7)

TrainSet2 <- xgb.DMatrix(data=as.matrix(data.rf[data.rf$y>0, c(covar.names)]),label=data.rf[data.rf$y>0, c('log_ba')])
# 
set.seed(1234)
model.ba <- xgb.train(params = params_ba, data = TrainSet2, nrounds = params_ba$nrounds)


# data.fit2[data.fit2$time.idx <= 84, 'score_ba'] <- data.rf$score_ba
# data.fit2[data.fit2$time.idx > 84, 'score_ba'] <- predict(model.ba, as.matrix(data.fit2[data.fit2$time.idx > 84,c(covar.names,'score_z')]))

data.fit2[data.fit2$time.idx <= 84, 'score_ba'] <- data.rf$score_ba
data.fit2[data.fit2$time.idx > 84 & data.fit2$y>0, 'score_ba'] <- predict(model.ba, as.matrix(data.fit2[data.fit2$time.idx > 84 & data.fit2$y>0,c(covar.names)]))



#without coordinates
# params_cnt <- list(
#   eta = 0.1,
#   max_depth = 3,
#   gamma = 1,
#   colsample_bytree = 0.7,
#   min_child_weight = 3,
#   subsample = 0.7,
#   nrounds = 50,
#   objective = 'count:poisson',
#   verbosity = 0
# )

#with coordiantes
params_cnt <- list(
  eta = 0.09,
  max_depth = 5,
  gamma = 0.8,
  colsample_bytree = 0.7,
  min_child_weight = 7,
  subsample = 1,
  nrounds = 70,
  objective = 'count:poisson',
  verbosity = 0
)
TrainSet3 <- xgb.DMatrix(data=as.matrix(data.rf[, c(covar.names,'score_z')]),label=data.rf[, c('y')])

set.seed(1234)
model.cnt <- xgb.train(params = params_cnt, data = TrainSet3, nrounds = params_cnt$nrounds)
data.rf$score_cnt <- kfold_cv_pred(data.rf[, c(covar.names,'score_z')], data.rf[, 'y'], params_cnt, 7)


data.fit2[data.fit2$time.idx <= 84, 'score_cnt'] <- data.rf$score_cnt
data.fit2[data.fit2$time.idx > 84, 'score_cnt'] <- predict(model.cnt, as.matrix(data.fit2[data.fit2$time.idx > 84,c(covar.names,'score_z')]))





data.fit2[is.na(data.fit2$log_ba),'log_ba'] <- 0

auc(data.fit2[data.fit2$year.idx<=7,'z'], data.fit2[data.fit2$year.idx<=7,'score_z'], quiet = TRUE)
auc(data.fit2[data.fit2$year.idx>7,'z'], data.fit2[data.fit2$year.idx>7,'score_z'], quiet = TRUE)

# sum((data.fit2[data.fit2$year.idx<=7,'log_ba']-data.fit2[data.fit2$year.idx<=7,'score_ba'])^2)/sum(data.fit2$year.idx<=7)
# sum((data.fit2[data.fit2$year.idx>7,'log_ba']-data.fit2[data.fit2$year.idx>7,'score_ba'])^2)/sum(data.fit2$year.idx>7)

sum((data.fit2[data.fit2$year.idx<=7 & data.fit2$y>0,'log_ba']-data.fit2[data.fit2$year.idx<=7 & data.fit2$y>0,'score_ba'])^2)/sum(data.fit2$year.idx<=7)
sum((data.fit2[data.fit2$year.idx>7 & data.fit2$y>0,'log_ba']-data.fit2[data.fit2$year.idx>7 & data.fit2$y>0,'score_ba'])^2)/sum(data.fit2$year.idx>7)


sum((data.fit2[data.fit2$year.idx<=7,'y']-data.fit2[data.fit2$year.idx<=7,'score_cnt'])^2)/sum(data.fit2$year.idx<=7)
sum((data.fit2[data.fit2$year.idx>7,'y']-data.fit2[data.fit2$year.idx>7,'score_cnt'])^2)/sum(data.fit2$year.idx>7)



# for (var in c('score_cnt','score_ba','score_z') ){
#   data.fit2[,var] <- (data.fit2[,var]-mean(data.fit2[,var],na.rm=T))/sd(data.fit2[,var],na.rm=T)
# }

data.fit3 <- reshape(data.fit2[,c('grid.idx','time.idx','y')],
                     timevar = "time.idx",
                     idvar = "grid.idx",
                     direction = "wide"
)

B2.merge <- merge(B2, data.fit3, by = "grid.idx")


B2.adj <- poly2nb(B2)
nb2INLA("map.adj", B2.adj)
g <- inla.read.graph(filename = "map.adj")



n1 <- dim(data.fit2)[1]
#
nothing1 <- rep(NA, n1)
nothing2 <- rep(NA, n1)

data.fit2 <- data.fit2[order(data.fit2$time.idx),]

z <- as.vector((data.fit2$y>0)+0)
log.ba <- as.vector(ifelse(data.fit2$y>0, data.fit2$log_ba, NA))
cnt = as.vector(data.fit2$y)
# 

#prepare for prediction
z[which(data.fit2$time.idx>=85)] <- NA
log.ba[which(data.fit2$time.idx>=85)] <- NA
cnt[which(data.fit2$time.idx>=85)] <- NA



cntNA <- as.vector(c(cnt,nothing1,nothing2))
zNA <- as.vector(c(nothing1, z, nothing2))
baNA.log = as.vector(c(nothing1, nothing2, log.ba))


outcome.matrix.log <- matrix(c(cntNA,zNA, baNA.log), ncol=3)


Intercept1 <- c(rep(1,n1),nothing1, nothing2)
Intercept2 <- c(nothing1, rep(1,n1), nothing2)
Intercept3 <- c(nothing1, nothing2, rep(1,n1))


i.spat1 = c(data.fit2$grid.idx, nothing1, nothing2)# fire ignition
# i.spat.iid = c(data.fit2$grid.idx, nothing1, nothing2)# fire ignition
i.spat2 = c(nothing1,data.fit2$grid.idx, nothing2)# BA ind
i.spat3 = c(nothing1, nothing2, data.fit2$grid.idx)# BA
# i.spat3.iid = c(nothing1, nothing2, data.fit2$grid.idx)# BA


time.idx1 <- c(data.fit2$month.idx, nothing1, nothing2)# fire ignition
time.idx2 <- c(nothing1, data.fit2$month.idx, nothing2)# BA ind
time.idx3  <-  c(nothing1, nothing2, data.fit2$month.idx)# BA

score_cnt_grp <- inla.group(data.fit2$score_cnt, n = 10, method = "cut")
score_1 <- c(data.fit2$score_cnt, nothing1, nothing2)
score_1_grp <- c(score_cnt_grp, nothing1, nothing2)

score_z_grp <- inla.group(data.fit2$score_z, n = 20, method = "quantile")
score_2 <- c(nothing1,  data.fit2$score_z, nothing2)
score_2_grp <- c(nothing1, score_z_grp, nothing2)

score_ba_grp <- rep(NA, nrow(data.fit2))
score_ba_grp[which(data.fit2$y>0)] <- inla.group(data.fit2[data.fit2$y>0,]$score_ba, n = 20, method = "cut")
  
score_3 <- c(nothing1, nothing2, data.fit2$score_ba)
score_3_grp <- c(nothing1, nothing2, score_ba_grp)


# mesh_score_1 <- inla.mesh.1d(seq(0,2,by=0.25),boundary=c('dirichlet','free')) 
mesh_score_1 <- inla.mesh.1d(seq(0.02,1.6,by=0.5),boundary=c('free','free')) 
A1 <- inla.spde.make.A(mesh_score_1, loc=data.fit2$score_cnt)
spde_score_1 <-  inla.spde2.pcmatern(mesh_score_1, 
                                     prior.range = c(1, 0.05),
                                     prior.sigma = c(1, 0.05))
# spde_score_1 <-  inla.spde2.matern(mesh_score_1, constr = TRUE)

spde_score_1.idx <- inla.spde.make.index("score_1", n.spde = spde_score_1$n.spde)

mesh_score_2 <- inla.mesh.1d(seq(0.03, 0.9, by = 0.1),boundary=c('free','free'))
# mesh_score_2 <- inla.mesh.1d(seq(0, 1, by = 0.2),boundary=c('dirichlet','dirichlet'))
A2 <- inla.spde.make.A(mesh_score_2, loc=data.fit2$score_z)
spde_score_2 <-  inla.spde2.pcmatern(mesh_score_2, 
                                     prior.range = c(0.2, 0.05),
                                     prior.sigma = c(1, 0.05))
# spde_score_2 <-  inla.spde2.matern(mesh_score_2, constr = TRUE)
spde_score_2.idx <- inla.spde.make.index("score_2", n.spde = spde_score_2$n.spde)


data.fit2[is.na(data.fit2$score_ba),'score_ba'] <- 0
mesh_score_3 <- inla.mesh.1d(seq(2.6, 7, by = 1),boundary=c('free','free'))
# mesh_score_3 <- inla.mesh.1d(seq(0, 7, by = 1),boundary=c('dirichlet','free'))
A3 <- inla.spde.make.A(mesh_score_3, loc=data.fit2$score_ba)
spde_score_3 <-  inla.spde2.pcmatern(mesh_score_3, 
                                     prior.range = c(2, 0.05),
                                     prior.sigma = c(1, 0.05))
# spde_score_3 <-  inla.spde2.matern(mesh_score_3, constr = TRUE)
spde_score_3.idx <- inla.spde.make.index("score_3", n.spde = spde_score_3$n.spde)



cnt.stack <- inla.stack(
  data= list(Y.log=cbind(cnt,NA,NA)),
  A <- list(1,1,1,A1,1),
  effect = list(Intercept1=rep(1,nrow(data.fit2)), idarea1=data.fit2$grid.idx, time.idx1=data.fit2$month.idx, score_1=spde_score_1.idx,
                score_1_grp = score_cnt_grp),
  tag='cnt'
)

z.stack <- inla.stack(
  data= list(Y.log=cbind(NA,z,NA)),
  A <- list(1,1,1,A2,1),
  effect = list(Intercept2=rep(1,nrow(data.fit2)), idarea2=data.fit2$grid.idx, time.idx2=data.fit2$month.idx, score_2=spde_score_2.idx,
                score_2_grp=score_z_grp),
  tag='z'
)

ba.stack <- inla.stack(
  data= list(Y.log=cbind(NA,NA,log.ba)),
  A <- list(1,1,1, A3,1),
  effect = list(Intercept3=rep(1,nrow(data.fit2)), idarea3=data.fit2$grid.idx, time.idx3=data.fit2$month.idx, score_3=spde_score_3.idx,
                score_3_grp=score_ba_grp),
  tag='ba'
)

all.stack <- inla.stack(cnt.stack, z.stack, ba.stack )


data=list(Y.log=outcome.matrix.log, 

          idarea1=i.spat1, 
          idarea2=i.spat2, 
          idarea3=i.spat3,
          
          idarea1.iid = i.spat1, 
          idarea3.iid = i.spat3,
          
          time.idx1 = time.idx1,
          time.idx2 = time.idx2,
          time.idx3 = time.idx3,
          
          Intercept1 = Intercept1,
          Intercept2 = Intercept2,
          Intercept3 = Intercept3,

          score_1_grp = score_1_grp,
          score_2_grp = score_2_grp,
          score_3_grp = score_3_grp,
          
          score_1 = score_1,
          score_2 = score_2,
          score_3 = score_3
)

library("splines")
hyper.rw <-  list(prec = list(prior="loggamma",param=c(1,1)))
# formula1 <- Y.log ~  -1  +
#   Intercept1 +  f(idarea1, copy='idarea2',fixed=F) + f(time.idx1, model='rw1') + ns(score_1, df=10)+
#   Intercept2 +  f(idarea2, model='bym2', graph=g) + f(time.idx2, model='rw1') + ns(score_2, df=10)+
#   Intercept3 + f(idarea3,  copy='idarea2',fixed=F) + f(time.idx3, model='rw1') + ns(score_3, df=10)

formula1 <- Y.log ~  -1  +
  Intercept1 +  f(idarea1, copy='idarea2',fixed=F) + f(time.idx1, model='rw1') + f(score_1_grp, model='rw1',constr = TRUE )+
  Intercept2 +  f(idarea2, model='bym2', graph=g) + f(time.idx2, model='rw1') + f(score_2_grp, model='rw1',constr = TRUE )+
  Intercept3 + f(idarea3,  copy='idarea2',fixed=F) + f(time.idx3, model='rw1',hyper=hyper.rw) + f(score_3_grp, model='rw1',constr = TRUE , hyper=hyper.rw)

formula2 <- Y.log ~  -1  +
  Intercept1 +  f(idarea1, copy='idarea2',fixed=F) +  f(time.idx1, model='rw1',hyper=hyper.rw) + f(score_1, model=spde_score_1 )+
  Intercept2 +  f(idarea2, model='bym2', graph=g) +  f(time.idx2, model='rw1',hyper=hyper.rw) + f(score_2, model=spde_score_2 )+
  Intercept3 + f(idarea3,  copy='idarea2',fixed=F) +  f(time.idx3, model='rw1',hyper=hyper.rw) +f(score_3, model=spde_score_3 )


t1 <- Sys.time()
res1 <- inla(formula1,
               family = c('poisson','binomial', 'gamma'), data = data,  Ntrials=1,
               control.predictor = list(compute = TRUE, link=c(rep(1,n1),rep(2,n1),rep(3,n1))),
               verbose=TRUE,
               control.compute=list(config = TRUE),
               control.family = list( list(), list(), list()),
               control.fixed = list(expand.factor.strategy = 'inla')
)
t2 <- Sys.time()
print(t2-t1)
summary(res1)

save(res1, file=file.path(dir.out,'Final_Model_1.RData'))



t1 <- Sys.time()
res1.1 <- inla(formula2,
             family = c('poisson','binomial', 'gamma'), data = inla.stack.data(all.stack),  Ntrials=1,
             control.predictor = list(A = inla.stack.A(all.stack), compute = TRUE, link=c(rep(1,n1),rep(2,n1),rep(3,n1))),
             verbose=TRUE,
             control.compute=list(config = TRUE),
             control.family = list( list(), list(), list()),
             control.fixed = list(expand.factor.strategy = 'inla')
)
t2 <- Sys.time()
print(t2-t1)
summary(res1.1)



t1 <- Sys.time()
res2 <- inla(formula1,
             family = c('poisson','binomial', 'gp'), data = data,  Ntrials=1,
             control.predictor = list(compute = TRUE, link=c(rep(1,n1),rep(2,n1),rep(3,n1))),
             verbose=TRUE,
             control.compute=list(config = TRUE),
             control.family = list( list(), list(), list(control.link=list(quantile=0.5))),
             control.fixed = list(expand.factor.strategy = 'inla')
)
t2 <- Sys.time()
print(t2-t1)
summary(res2)

save(res2, file=file.path(dir.out,'Final_Model_2.RData'))

t1 <- Sys.time()
res2.1 <- inla(formula2,
               family = c('poisson','binomial', 'gp'), data = inla.stack.data(all.stack),  Ntrials=1,
               control.predictor = list(A = inla.stack.A(all.stack), compute = TRUE, link=c(rep(1,n1),rep(2,n1),rep(3,n1))),
               verbose=TRUE,
               control.compute=list(config = TRUE),
               control.family = list( list(), list(), list(control.link=list(quantile=0.5))),
               control.fixed = list(expand.factor.strategy = 'inla')
)
t2 <- Sys.time()
print(t2-t1)
summary(res2.1)

t1 <- Sys.time()
res3 <- inla(formula1,
             family = c('poisson','binomial', 'weibull'), data = data,  Ntrials=1,
             control.predictor = list(compute = TRUE, link=c(rep(1,n1),rep(2,n1),rep(3,n1))),
             verbose=TRUE,
             control.compute=list(config = TRUE),
             control.family = list( list(), list(), list()),
             control.fixed = list(expand.factor.strategy = 'inla'))
t2 <- Sys.time()
print(t2-t1)
summary(res3)

save(res3, file=file.path(dir.out,'Final_Model_3.RData'))



t1 <- Sys.time()
res3.1 <- inla(formula2,
               family = c('poisson','binomial', 'weibull'), data = inla.stack.data(all.stack),  Ntrials=1,
               control.predictor = list(A = inla.stack.A(all.stack), compute = TRUE, link=c(rep(1,n1),rep(2,n1),rep(3,n1))),
               verbose=TRUE,
               control.compute=list(config = TRUE),
               control.family = list( list(), list(), list()),
               control.fixed = list(expand.factor.strategy = 'inla')
)
t2 <- Sys.time()
print(t2-t1)
summary(res3.1)

save(res3.1, file=file.path(dir.out,'Final_Model_3.1.RData'))


post.pred.gpd.parallel <- function(samples, idx.pred.pois, idx.pred.z, idx.pred.ba, alpha,n.samples=200 ){
  rgp = function(n, sigma, eta, alpha, xi = 0.001)
  {
    if (missing(sigma)) {
      stopifnot(!missing(eta) && !missing(alpha))
      sigma = exp(eta) * xi / ((1.0 - alpha)^(-xi) -1.0)
    }
    return (sigma / xi * (runif(n)^(-xi) -1.0))
  }
  t1 <- Sys.time()
  res <- foreach(j = 1:nrow(data.fit2)) %dopar%{
    set.seed(j)
    pred.cnt <- rep(NA, n.samples)
    pred.z <- rep(NA, n.samples)
    pred.ba <- rep(NA, n.samples)
    for (i in 1:n.samples){
      eta.pois <- samples[[i]]$latent[idx.pred.pois,1]
      
      eta.z <- samples[[i]]$latent[idx.pred.z,1]
      p <- exp(eta.z)/(1 + exp(eta.z))
      
      
      eta.ba <- samples[[i]]$latent[idx.pred.ba,1]
      xi <- samples[[i]]$hyperpar[1]
      
      lambda <- exp(eta.pois)
      
      pred.cnt[i] <- rpois(1, lambda[j] )
      
      z <- rbinom(1, size=1, prob=p[j])
      pred.z[i] <- z
      if (z==1){
        pred.ba[i] <- rgp(1,eta=eta.ba[j], alpha=alpha, xi=xi)
      }else{
        pred.ba[i] <- 0
      }
    }
    res.list <- list('pred.cnt'=pred.cnt, 'pred.z'=pred.z,'pred.ba'=pred.ba)
    
  }
  t2 <- Sys.time()
  print(t2-t1)
  pred.cnt <- list()
  pred.z <- list()
  pred.ba <- list()
  for (i in 1:nrow(data.fit2)){
    pred.cnt[[i]] <- res[[i]][['pred.cnt']]
    pred.z[[i]] <- res[[i]][['pred.z']]
    pred.ba[[i]] <- res[[i]][['pred.ba']]
  }
  
  
  return(list('pred.ba'=pred.ba,'pred.cnt'=pred.cnt, 'pred.z'=pred.z))
}



post.pred.gamma.parallel <- function(samples, idx.pred.pois, idx.pred.z, idx.pred.ba,n.samples=200 ){

  t1 <- Sys.time()
  res <- foreach(j = 1:nrow(data.fit2)) %dopar%{
    set.seed(j)
    pred.cnt <- rep(NA, n.samples)
    pred.z <- rep(NA, n.samples)
    pred.ba <- rep(NA, n.samples)
    for (i in 1:n.samples){
      eta.pois <- samples[[i]]$latent[idx.pred.pois,1]
      
      eta.z <- samples[[i]]$latent[idx.pred.z,1]
      p <- exp(eta.z)/(1 + exp(eta.z))
      
      
      eta.ba <- samples[[i]]$latent[idx.pred.ba,1]
      prec.par <- samples[[i]]$hyperpar[1]
      
      a = prec.par
      b = eta.ba / a
      
      
      lambda <- exp(eta.pois)
      
      pred.cnt[i] <- rpois(1, lambda[j] )
      
      z <- rbinom(1, size=1, prob=p[j])
      pred.z[i] <- z
      if (z==1){
        pred.ba[i] <- rgamma(1, shape = a, scale = b[j])
      }else{
        pred.ba[i] <- 0
      }
    }
    res.list <- list('pred.cnt'=pred.cnt, 'pred.z'=pred.z,'pred.ba'=pred.ba)
    
  }
  t2 <- Sys.time()
  print(t2-t1)
  pred.cnt <- list()
  pred.z <- list()
  pred.ba <- list()
  for (i in 1:nrow(data.fit2)){
    pred.cnt[[i]] <- res[[i]][['pred.cnt']]
    pred.z[[i]] <- res[[i]][['pred.z']]
    pred.ba[[i]] <- res[[i]][['pred.ba']]
  }
  
  
  return(list('pred.ba'=pred.ba,'pred.cnt'=pred.cnt, 'pred.z'=pred.z))
}


post.pred.weibull.parallel <- function(samples, idx.pred.pois, idx.pred.z, idx.pred.ba,n.samples=200 ){
  t1 <- Sys.time()
  res <- foreach(j = 1:nrow(data.fit2)) %dopar%{
    set.seed(j)
    pred.cnt <- rep(NA, n.samples)
    pred.z <- rep(NA, n.samples)
    pred.ba <- rep(NA, n.samples)
    for (i in 1:n.samples){
      eta.pois <- samples[[i]]$latent[idx.pred.pois,1]
      
      eta.z <- samples[[i]]$latent[idx.pred.z,1]
      p <- exp(eta.z)/(1 + exp(eta.z))
      
      
      eta.ba <- samples[[i]]$latent[idx.pred.ba,1]
      
      alpha.ba <- samples[[i]]$hyperpar[1]
      
      lambda.ba = exp(eta.ba)
      
      
      lambda <- exp(eta.pois)
      
      pred.cnt[i] <- rpois(1, lambda[j] )
      
      z <- rbinom(1, size=1, prob=p[j])
      pred.z[i] <- z
      if (z==1){
        # pred.ba[i] <- rgamma(1, shape = a, scale = b[j])
        pred.ba[i] <- rweibull(1, shape = alpha.ba, scale = lambda.ba[j]^(-1/alpha.ba))
      }else{
        pred.ba[i] <- 0
      }
    }
    res.list <- list('pred.cnt'=pred.cnt, 'pred.z'=pred.z,'pred.ba'=pred.ba)
    
  }
  t2 <- Sys.time()
  print(t2-t1)
  pred.cnt <- list()
  pred.z <- list()
  pred.ba <- list()
  for (i in 1:nrow(data.fit2)){
    pred.cnt[[i]] <- res[[i]][['pred.cnt']]
    pred.z[[i]] <- res[[i]][['pred.z']]
    pred.ba[[i]] <- res[[i]][['pred.ba']]
  }
  
  
  return(list('pred.ba'=pred.ba,'pred.cnt'=pred.cnt, 'pred.z'=pred.z))
}


library(doParallel)
cl <- makeCluster(3)
registerDoParallel(cl)

n.samples = 200

n1 <- 192*108


# idx.pred.pois <- 1:n1
# idx.pred.z <- (n1+1):(2*n1)
# idx.pred.ba <- (2*n1+1):(3*n1)
# result <- res0
# 
# 
# # idx.pred.pois <- 1:n1
# # idx.pred.ba <- (n1+1):(2*n1)
# result <- res0
# 
samples = inla.posterior.sample(1, result = res3.1, seed=1234)

idx.pred.pois <- 1:n1
idx.pred.z <- (n1+1):(2*n1)
idx.pred.ba <- (2*n1+1):(3*n1)
idx.time.idx1 <- (3*n1+1):(3*n1+12)

idx.idarea2.u <- (3*n1+13):(3*n1+204)
idx.idarea2.v <- (3*n1+205):(3*n1+396)

idx.time.idx2 <- (3*n1+ 397):(3*n1+408)
idx.time.idx3 <- (3*n1+ 409):(3*n1+420)

idx.idarea1.u <- (3*n1+421):(3*n1+612)
idx.idarea1.v <- (3*n1+613):(3*n1+804)

idx.idarea3.u <- (3*n1+805):(3*n1+996)
idx.idarea3.v <- (3*n1+997):(3*n1+1188)

idx.Intercept1 <- 3*n1+1189
idx.ns_score1 <- (3*n1+1190):(3*n1+1199)

idx.Intercept2 <- 3*n1+1200
idx.ns_score2 <- (3*n1+1201):(3*n1+1210)

idx.Intercept3 <- 3*n1+1211
idx.ns_score3 <- (3*n1+1212):(3*n1+1221)


effect <- c()
for (i in 1:n.samples){
  effect <- rbind(effect, samples[[i]]$latent[idx.ns_score1,])
}

attr(ns(score_1, df = 10), "knots")
attr(ns(score_1, df = 10), "Boundary.knots")

plot(attr(ns(score_2, df = 10), "knots"), colMeans(effect),type='l',ylim=c(-10,10))
lines(1:10, apply(effect,2,quantile,0.975),col='red')
lines(1:10, apply(effect,2,quantile,0.025),col='red')



summary((res3.1$summary.fitted.values[idx.pred.pois,'mean']-data.fit2$y)[which(data.fit2$y>0)])
summary((res3$summary.fitted.values[idx.pred.pois,'mean']-data.fit2$y)[which(data.fit2$y>0)])



summary(res3$summary.fitted.values[idx.pred.pois,'mean'])

set.seed(123)
test <- rnorm(100)  # Generate some normally distributed data

# Fit the natural spline with df = 10
spline_fit <- ns(test, df = 10)

# Extract the internal knots
internal_knots <- attr(spline_fit, "knots")

# Display the internal knots
print(internal_knots)

# Plot the histogram of the data
hist(score_1, breaks = 20, main = "Data Distribution and Knot Placement", xlab = "score_1", col = "lightblue", border = "white")

# Add vertical lines for the internal knots
abline(v = internal_knots, col = "red", lwd = 2, lty = 2)
legend("topright", legend = "Internal Knots", col = "red", lty = 2, lwd = 2)




n <- 10
d <- data.frame(x = 1:n, y =  colMeans(effect))
ggplot(d,aes(x,y)) + geom_point() + 
  geom_line(data=data.frame(spline(d, n=100*n)))

# load(file=file.path(dir.out,'Final_Model_2.RData'))


load(file=file.path(dir.out,'Final_Model_3.1_pred.sp.RData'))
pred.cnt <- pred.sp$pred.cnt
pred.ba <- pred.sp$pred.ba
pred.z <- pred.sp$pred.z
# result <- res1

n1 <- 192*108
idx.pred.pois <- 1:n1
idx.pred.z <- (n1+1):(2*n1)
idx.pred.ba <- (2*n1+1):(3*n1)

crps.ba <- rep(NA, nrow(data.fit2))
for (i in 1:nrow(data.fit2) ){
  y <- log(data.fit2$area_ha[i])
  if (is.na(y)) y <- 0
  crps.ba[i] <- crps_sample(
    y,
    pred.ba[[i]],
    method = "edf")
}

ins.idx <- 192*84
round(sum(crps.ba[1:ins.idx])/ins.idx,4)
round(sum(crps.ba[(ins.idx+1):nrow(data.fit2)])/(nrow(data.fit2)-ins.idx),4)
round(sum(crps.ba)/length(crps.ba),4)

crps.cnt <- rep(NA, nrow(data.fit2))
for (i in 1:nrow(data.fit2) ){
  y <- data.fit2$y[i]
  if (is.na(y)) y <- 0
  crps.cnt[i] <- crps_sample(
    y,
    pred.ba[[i]],
    method = "edf")
}

ins.idx <- 192*84
round(sum(crps.cnt[1:ins.idx])/ins.idx,4)
round(sum(crps.cnt[(ins.idx+1):nrow(data.fit2)])/(nrow(data.fit2)-ins.idx),4)
round(sum(crps.cnt)/length(crps.cnt),4)




data.fit2[,'Estimated_Lambda'] <- sapply(pred.cnt,mean)
data.fit2$Lower_CNT <- sapply(pred.cnt,quantile,0.025)
data.fit2$Upper_CNT <- sapply(pred.cnt,quantile,0.975)
data.fit2$Scaling_Residuals <- data.fit2$y - data.fit2$Estimated_Lambda
summary(data.fit2$Scaling_Residuals)

data.fit2[,'Estimated_Z'] <- sapply(pred.z,mean)
data.fit2$Lower_Z <- sapply(pred.z,quantile,0.025)
data.fit2$Upper_Z <- sapply(pred.z,quantile,0.975)

data.fit2[,'Estimated_BA'] <- sapply(pred.ba,mean)
data.fit2[,'Lower_BA'] <- sapply(pred.ba,quantile,0.025)
data.fit2[,'Upper_BA'] <- sapply(pred.ba,quantile,0.975)




########################   Visualization #################################
dist<-shapefile(file.path("/home/pgrad2/2448355h/My_PhD_Project/00_Dataset/Urban_Fires","distritos.shp"))
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


grid.cell.coord <- st_as_sf(data.fit2, coords = c("lon.grid", "lat.grid"), crs = 4326)

# merged_sf <- st_join(grid.cell.coord, sf_districts[,'NAME_1'], join = st_within)
merged_sf1 <- st_join(grid.cell.coord, sf_districts[,'NAME_1'], join = st_nearest_feature)

merged_sf1 %>% group_by(NAME_1) %>% summarize(n_cnt = sum(y))




#----------------------boxplot of CNT/BA in single district --------------
ggplot() +geom_sf(data=merged_sf1, aes(fill=NAME_1,col=NAME_1)) + 
  geom_sf(data = sf_districts,color = "black", fill = NA)


dist.name <- 'Castelo Branco'
joint.post.sp <- function(x) Reduce("+", x)
df.dist <- data.frame(month=rep(1:108, each=200))
df.boxplot.true <- data.frame(month=1:108, cnt=rep(NA, 108), log.ba= rep(NA,108))
for (t in 1:108){
  idx <- which(merged_sf1$NAME_1==dist.name & merged_sf1$time.idx==t)
  cnt.true <- sum(data.fit2[idx,'y'])
  ba.true <- sum(data.fit2[idx,'log_ba'])
  
  df.boxplot.true[df.boxplot.true$month==t,c('cnt','log.ba')] <- c(cnt.true, ba.true)
  
  pred.cnt.dist <- joint.post.sp(pred.cnt[idx])
  pred.ba.dist <- joint.post.sp(pred.ba[idx])
  df.dist[df.dist$month==t,'sample_cnt'] <- pred.cnt.dist
  df.dist[df.dist$month==t,'sample_ba'] <- pred.ba.dist
}


ggplot(df.dist[df.dist$month>=85,], aes(x = factor(month), y = sample_ba)) +
  geom_boxplot() + 
  geom_line(data=df.boxplot.true[df.boxplot.true$month>=85,], aes(x=factor(month),y=log.ba,group = 1), col='red',linewidth=1)+
  labs(x = "Month", y = "Posterior Predictive Sample", title = "Posterior Predictive Samples over Time") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))  # Rotate x-axis labels for clarity


ggplot(df.dist[df.dist$month>=85,], aes(x = factor(month), y = sample_cnt)) +
  geom_boxplot() + 
  geom_line(data=df.boxplot.true[df.boxplot.true$month>=85,], aes(x=factor(month),y=cnt,group = 1), col='red',linewidth=1)+
  labs(x = "Month", y = "Posterior Predictive Sample", title = "Posterior Predictive Samples over Time") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))  # Rotate x-axis labels for clarity


ggplot(df.dist[df.dist$month<=84  & df.dist$month>=61 ,], aes(x = factor(month), y = sample_ba)) +
  geom_boxplot() + 
  geom_line(data=df.boxplot.true[df.boxplot.true$month<=84 & df.boxplot.true$month>=61,], aes(x=factor(month),y=log.ba,group = 1), col='red',linewidth=1)+
  labs(x = "Month", y = "Posterior Predictive Sample", title = "Posterior Predictive Samples over Time") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))  # Rotate x-axis labels for clarity


ggplot(df.dist[df.dist$month<=84  & df.dist$month>=61,], aes(x = factor(month), y = sample_cnt)) +
  geom_boxplot() + 
  geom_line(data=df.boxplot.true[df.boxplot.true$month<=84 & df.boxplot.true$month>=61,], aes(x=factor(month),y=cnt,group = 1), col='red',linewidth=1)+
  labs(x = "Month", y = "Posterior Predictive Sample", title = "Posterior Predictive Samples over Time") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))  # Rotate x-axis labels for clarity

########################################################################

head(merged_sf1)




B2_sf <- st_as_sf(B2.merge)
library(tidyr)
B2_sf <- gather(B2_sf, y.time.idx, y, paste0("y.", 1:108))
B2_sf$time.idx <- as.integer(substring(B2_sf$y.time.idx, 3, 5))


B2_sf <- merge(
  B2_sf[,c('grid.idx','time.idx')], data.fit2,
  
  by.x=c('grid.idx','time.idx'),
  by.y=c('grid.idx','time.idx'),
)

Year <- 6
# data.fit2$year.idx <- 1
# B2_sf$year.idx <- 1

lprange.scale <- c(max(abs(data.fit2[data.fit2$year.idx==Year, ]$y)))*c(-1,1)

csc.scale <- scale_fill_gradientn(colours = brewer.pal(3, "RdBu"), limits = lprange.scale)
# csc.scale <- scale_fill_gradientn(colours = c("#F7F7F7","#67A9CF"), limits = lprange.scale)
ggplot() + geom_sf(data=B2_sf[B2_sf$year.idx==Year,],aes(fill = y),lwd = 0.1) +
  # geom_sf(data=st_as_sf(district),alpha = 0)+
  facet_wrap(~time.idx, dir = "h", ncol =3) +
  ggtitle(paste("Actual y in 2012")) +
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)
        # axis.text.y = element_blank(),
        # axis.ticks = element_blank()
  )+csc.scale

lprange.scale <- c(max(abs(data.fit2[data.fit2$year.idx==Year, ]$sqrt_ba)))*c(-1,1)
csc.scale <- scale_fill_gradientn(colours = brewer.pal(3, "RdBu"), limits = lprange.scale)
# csc.scale <- scale_fill_gradientn(colours = c("#F7F7F7","#67A9CF"), limits = lprange.scale)
ggplot() + geom_sf(data=B2_sf[B2_sf$year.idx==Year,],aes(fill = sqrt_ba),lwd = 0.1) +
  # geom_sf(data=st_as_sf(district),alpha = 0)+
  facet_wrap(~time.idx, dir = "h", ncol =3) +
  ggtitle(paste("Estimated Burn area in Sep from 2012-2020")) +
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)
        # axis.text.y = element_blank(),
        # axis.ticks = element_blank()
  )+csc.scale

lprange.scale <- max(
  abs(data.fit2[data.fit2$year.idx==Year, ]$Latent_Effect_Z),na.rm=T
)*c(-1,1)

csc.scale <- scale_fill_gradientn(colours = brewer.pal(3, "RdBu"), limits = lprange.scale)

ggplot() + geom_sf(data=B2_sf[B2_sf$year.idx==Year,],aes(fill = Latent_Effect_Z),lwd = 0.1) +
  # geom_sf(data=st_as_sf(district),alpha = 0)+
  facet_wrap(~time.idx, dir = "h", ncol =3) +
  ggtitle(paste("Estimated Latent Z in 2012")) +
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)
        # axis.text.y = element_blank(),
        # axis.ticks = element_blank()
  )+csc.scale


lprange.scale <- max(
  abs(data.fit2[data.fit2$year.idx==Year, ]$Latent_Effect_BA),na.rm=T
)*c(-1,1)

csc.scale <- scale_fill_gradientn(colours = brewer.pal(3, "RdBu"), limits = lprange.scale)

ggplot() + geom_sf(data=B2_sf[B2_sf$year.idx==Year,],aes(fill = Latent_Effect_BA),lwd = 0.1) +
  # geom_sf(data=st_as_sf(district),alpha = 0)+
  facet_wrap(~time.idx, dir = "h", ncol =3) +
  ggtitle(paste("Estimated Latent BA in 2012")) +
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)
        # axis.text.y = element_blank(),
        # axis.ticks = element_blank()
  )+csc.scale

lprange.scale <- max(
  abs(data.fit2[data.fit2$year.idx==Year, ]$Estimated_BA),na.rm=T
)*c(-1,1)

csc.scale <- scale_fill_gradientn(colours = brewer.pal(3, "RdBu"), limits = lprange.scale)

ggplot() + geom_sf(data=B2_sf[B2_sf$year.idx==Year,],aes(fill = Estimated_BA),lwd = 0.1) +
  # geom_sf(data=st_as_sf(district),alpha = 0)+
  facet_wrap(~time.idx, dir = "h", ncol =3) +
  ggtitle(paste("Estimated BA^{1/3} Means in 2012 with covariates")) +
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)
        # axis.text.y = element_blank(),
        # axis.ticks = element_blank()
  )+csc.scale


lprange.scale <- max(
  abs((df.plot[df.plot$year.idx==Year, ]$area_ha)^{1/3}),na.rm=T
)*c(-1,1)

csc.scale <- scale_fill_gradientn(colours = brewer.pal(3, "RdBu"), limits = lprange.scale)

B2_sf[is.na(B2_sf$area_ha),'area_ha'] <- 0
B2_sf$area_ha2 <- (B2_sf$area_ha)^{1/3}
ggplot() + geom_sf(data=B2_sf[B2_sf$year.idx==Year,],aes(fill = area_ha2),lwd = 0.1) +
  # geom_sf(data=st_as_sf(district),alpha = 0)+
  facet_wrap(~time.idx, dir = "h", ncol =3) +
  ggtitle(paste("BA^{1/3} in 2012")) +
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)
        # axis.text.y = element_blank(),
        # axis.ticks = element_blank()
  )+csc.scale

df.plot <- data.fit2
df.plot[is.na(df.plot$area_ha),'area_ha'] <- 1
df.plot[,'cover_ind'] <- df.plot[,'Upper_BA']>=log(df.plot[,'area_ha']) & df.plot[,'Lower_BA']<=log(df.plot[,'area_ha'])

View(df.plot %>% group_by(year.idx) %>% summarize(coverate=mean(cover_ind)))
View(df.plot[df.plot$area_ha>1,] %>% group_by(year.idx) %>% summarize(coverate=mean(cover_ind)))

View(df.plot %>% group_by(year.idx) %>% summarize(coverate=mean(crps)))

ggplot(df.plot[!is.na(df.plot$area_ha)&df.plot$year.idx==8,], aes(x = grid.idx, y = log(area_ha))) +
  geom_point() +
  geom_line(aes(y = Estimated_BA), linetype = "dashed", color = "blue") +
  geom_line(aes(y = Upper_BA), linetype = "dashed", color = "red") +
  geom_line(aes(y = Lower_BA), linetype = "dashed", color = "red") +
  labs(title = "Estimated Burn area with 95% credible band in Year 2020 (with predicator) ",
       x = "Grid Index",
       y = "Transformed BA") +
  facet_wrap(~time.idx, dir = "h", ncol =3) +
  theme_minimal()

df_long <- melt(df.plot, id.vars = "y")

# Create scatter plots with facet_wrap

var <- 'score_ba'
lprange.scale <- max(
  abs(df.plot[df.plot$year.idx==1, var ]),na.rm=T
)*c(-1,1)

csc.scale <- scale_fill_gradientn(colours = brewer.pal(3, "RdBu"), limits = lprange.scale)
ggplot() + geom_sf(data=B2_sf[B2_sf$year.idx==1,],aes(fill = .data[[var]]),lwd = 0.1) +
  # geom_sf(data=st_as_sf(district),alpha = 0)+
  facet_wrap(~time.idx, dir = "h", ncol =3) +
  ggtitle(paste(var,"in 2012")) +
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)
        # axis.text.y = element_blank(),
        # axis.ticks = element_blank()
  )+csc.scale

ggplot(df.plot, aes(x = .data[[var]], y =(area_ha)^{1/3})) +
  geom_point() + 
  labs(x = var, y = "BA^{1/3}") +
  geom_smooth(method = 'lm', color = 'red', se = TRUE) +
  theme_minimal()

ggplot(df.plot, aes(x = .data[[var]], y =z)) +
  geom_point() + 
  labs(x = var, y = "BA^{1/3}") +
  theme_minimal()


ggplot(df.plot[df.plot$area_ha>0 & df.plot$score<100 ,], aes(x = .data[[var]], y =(area_ha))) +
  geom_point() + 
  labs(x = var, y = "BA^{1/3}") +
  geom_smooth(method = 'lm', color = 'red', se = TRUE) +
  theme_minimal()



ggplot(df.plot[df.plot$area_ha>0,], aes(x = pca_comp1, y =pca_comp3, color=(area_ha)^{1/3})) +
  geom_point() +
  labs(x = "pca_comp1", y = "pca_comp3") +
  theme_minimal() + scale_color_gradientn(colours=c("#F7F7F7","#67A9CF"))

plot(df.plot[,'VComp'],df.plot[,'log_ba'])
res.pca$scores[,2]
plot(res.pca$scores[,2],df.plot[,'log_ba'])
