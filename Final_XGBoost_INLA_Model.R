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


# grd <- SpatialGrid(GridTopology(cellcentre.offset = c(-10, 36),
#                                 cellsize = c(0.0625,0.0625),
#                                 cells.dim = c(68, 116)))

grd <- SpatialGrid(GridTopology(cellcentre.offset = c(-10, 36),
                                cellsize = c(0.125,0.125),
                                cells.dim = c(34, 58)))

# grd <- SpatialGrid(GridTopology(cellcentre.offset = c(-10, 36),
#                                 cellsize = c(0.25,0.25),
#                                 cells.dim = c(17, 29)))

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

#####################################################################
# XGBoost
#####################################################################
# # 
# data.fit2 <- do.call(rbind, lapply(1:108, function(x) {B2@data$time.idx = x
# return(B2@data)}))
# print(dim(data.fit2))
# data.fit2 <- merge(data.fit2, data.fit.ba[,c('grid.idx','time.idx','year.idx', 'month.idx','y','area_ha','log_ba')],
#                    by=c('grid.idx','time.idx'),all.x=T)
# 
# 
# print(dim(data.fit2))
# data.fit2[is.na(data.fit2$y),'y'] <- 0
# summary(data.fit2$y)
# 
# 
# 
# data.fit2$month.idx <- (data.fit2$time.idx-1)%%12 + 1
# data.fit2$year.idx <- (data.fit2$time.idx-1)%/%12 + 1
# data.fit2$E1 <- 1
# #
# 
# f.get.cov <- function(dataset, cov.name){
#   time <- dataset$time.idx
# 
#   x <- dataset$grid.cent.x.utm
#   y <- dataset$grid.cent.y.utm
# 
# 
#   spp <- SpatialPoints(data.frame(x=x,y=y),
#                        proj4string=CRS('+proj=utm +zone=29 +datum=WGS84 +units=km +no_defs' ))
#   spp <- spTransform(spp, CRS("EPSG:4326"))
#   v <- rep(NA, nrow(dataset))
#   for (t in unique(time)){
#     spdf <- SpatialPixelsDataFrame(points = grid_pixels,
#                                    data = data.frame(var=as.vector(get(paste(cov.name,'.month',sep=''))[,,t]),
#                                                      time=t))
#     proj4string(spdf) <- CRS("EPSG:4326")
#     idx <- which(time==t)
#     v[idx] <- over(spp[idx,],spdf[,'var'])$var
#   }
#   return(v)
# }
# #
# data.fit2$FWI <- f.get.cov(data.fit2,'fwi')
# data.fit2$HVegCov <- f.get.cov(data.fit2,'HVegCov')
# data.fit2$HVegLAI <- f.get.cov(data.fit2, 'HVegLAI')
# data.fit2$HVegTyp <- as.factor(f.get.cov(data.fit2, 'HVegTyp'))
# data.fit2$LVegCov <- f.get.cov(data.fit2,'LVegCov')
# data.fit2$LVegLAI <- f.get.cov(data.fit2, 'LVegLAI')
# data.fit2$LVegTyp <- as.factor(f.get.cov(data.fit2, 'LVegTyp'))
# data.fit2$Pricp <- f.get.cov(data.fit2, 'pricp')
# data.fit2$RHumi <- f.get.cov(data.fit2, 'rhumi')
# data.fit2$Temp <- f.get.cov(data.fit2, 'temp')
# data.fit2$UComp <- f.get.cov(data.fit2, 'u')
# data.fit2$VComp <- f.get.cov(data.fit2, 'v')
# data.fit2$sqrt_ba <- sqrt(data.fit2$area_ha)
# data.fit2[is.na(data.fit2$sqrt_ba), 'sqrt_ba'] <- 0
# data.fit2$z <- as.vector((data.fit2$y>0)+0)
# #
# covar.names <- c('FWI','HVegCov','HVegLAI','LVegCov','LVegLAI','Pricp','RHumi','Temp',
#                  'UComp','VComp')
# for (var in covar.names ){
#   data.fit2[,var] <- (data.fit2[,var]-mean(data.fit2[,var]))/sd(data.fit2[,var])
# }
# 
# covar.names <- c(covar.names ,#'LVegTyp','HVegTyp',
#                  'lon.grid','lat.grid',
#                  'month.idx')
# 
# data.rf <- data.fit2[data.fit2$time.idx <= 96, c(covar.names,'grid.idx','time.idx','year.idx','area_ha','z','y')]
# 
# 
# 
# 
# #
# library(xgboost)
# # library(caret)
# library(pROC)
# data.rf$log_ba <- log(data.rf$area_ha)
# data.rf[is.na(data.rf$log_ba),'log_ba'] <- 0
# 
# kfold_cv_pred  <- function(data, target.name ,covar.names, params, k) {
#   auc_train_list <- numeric(k)
#   auc_test_list <- numeric(k)
# 
#   for (i in 1:k) {
#     index <- which(data$year.idx==i)
#     train_data <- data[-index, covar.names ]
#     train_target <- data[-index, target.name]
#     test_data <- data[index, covar.names]
#     test_target <- data[index, target.name ]
# 
#     dtrain <- xgb.DMatrix(data = as.matrix(train_data), label = train_target)
#     dtest <- xgb.DMatrix(data = as.matrix(test_data), label = test_target)
# 
#     set.seed(1234)
#     model <- xgb.train(params = params, data = dtrain, nrounds = params$nrounds)
# 
#     train_pred <- predict(model, dtrain)
#     test_pred <- predict(model, dtest)
# 
#     data[index,'score'] <- test_pred
#   }
# 
# 
#   return(data$score)
# }
# 
# params_z <- list(
#   nrounds = 100,
#   eta = 0.06,
#   max_depth = 3,
#   gamma = 0.2,
#   scale_pos_weight = 75,
#   colsample_bytree = 0.7,
#   min_child_weight = 7,
#   subsample = 0.9,
#   objective = "binary:logistic",
#   eval_metric = "auc",
#   verbosity = 0
# )
# #
# data.rf$score_z <- kfold_cv_pred(data.rf, 'z', covar.names ,params_z, 8)
# 
# 
# TrainSet1 <- xgb.DMatrix(data=as.matrix(data.rf[, c(covar.names)]),label=data.rf[, c('z')])
# set.seed(1234)
# model.z <- xgb.train(params = params_z, data = TrainSet1, nrounds = params_z$nrounds)
# 
# 
# data.fit2[data.fit2$time.idx <= 96, 'score_z'] <- data.rf$score_z
# data.fit2[data.fit2$time.idx > 96, 'score_z'] <- predict(model.z, as.matrix(data.fit2[data.fit2$time.idx > 96,c(covar.names)]))
# 
# #
# #
# # # with coordiantes
# params_ba <- list(
#   nrounds = 60,
#   eta = 0.06,
#   max_depth = 2,
#   gamma = 0.6,
#   colsample_bytree = 0.9,
#   min_child_weight = 2,
#   subsample = 0.6,
#   objective=c('reg:squarederror')
# )
# 
# 
# kfold_cv_pred_pos <- function(data, target.name ,covar.names, params, k) {
#   loss_train_list <- numeric(k)
#   loss_test_list <- numeric(k)
#   for (i in 1:k) {
#     index <- which(data$year.idx==i)
#     train_data <- data[-index, covar.names ]
#     train_target <- data[-index, target.name]
#     pos.ind.train <- train_target>0
#     train_data <- train_data[pos.ind.train,]
#     train_target <- train_target[pos.ind.train]
#     # print(table(data[-index,'grid.idx']))
# 
#     test_data <- data[index, covar.names]
#     test_target <- data[index, target.name ]
# 
# 
#     dtrain <- xgb.DMatrix(data = as.matrix(train_data), label = train_target)
#     dtest <- xgb.DMatrix(data = as.matrix(test_data), label = test_target)
# 
#     set.seed(1234)
#     model <- xgb.train(params = params, data = dtrain, nrounds = params$nrounds)
# 
#     train_pred <- predict(model, dtrain)
#     test_pred <- predict(model, dtest)
# 
#     data[index, 'score'] <- test_pred
#   }
# 
# 
#   return(data$score)
# }
# 
# 
# data.rf$score_ba <- kfold_cv_pred_pos(data.rf,'log_ba',covar.names, params_ba, 8)
# 
# TrainSet2 <- xgb.DMatrix(data=as.matrix(data.rf[data.rf$y>0, c(covar.names)]),label=data.rf[data.rf$y>0, c('log_ba')])
# #
# set.seed(1234)
# model.ba <- xgb.train(params = params_ba, data = TrainSet2, nrounds = params_ba$nrounds)
# 
# 
# data.fit2[data.fit2$time.idx <= 96, 'score_ba'] <- data.rf$score_ba
# data.fit2[data.fit2$time.idx > 96, 'score_ba'] <- predict(model.ba, as.matrix(data.fit2[data.fit2$time.idx > 96,c(covar.names)]))
# 
# 
# 
# params_cnt <- list(
#   nrounds = 100,
#   eta = 0.05,
#   max_depth = 3,
#   gamma = 0,
#   colsample_bytree = 0.5,
#   min_child_weight = 12,
#   subsample = 0.9,
#   objective=c('count:poisson')
# )
# 
# data.rf$score_cnt <- kfold_cv_pred_pos(data.rf,'y',covar.names, params_ba, 8)
# 
# 
# TrainSet3 <- xgb.DMatrix(data=as.matrix(data.rf[data.rf$y>0, c(covar.names)]),label=data.rf[data.rf$y>0, c('y')])
# 
# set.seed(1234)
# model.cnt <- xgb.train(params = params_cnt, data = TrainSet3, nrounds = params_cnt$nrounds)
# 
# 
# data.fit2[data.fit2$time.idx <= 96, 'score_cnt'] <- data.rf$score_cnt
# data.fit2[data.fit2$time.idx > 96 , 'score_cnt'] <- predict(model.cnt, as.matrix(data.fit2[data.fit2$time.idx > 96 ,c(covar.names)]))
# 
# 
# 
# data.fit2[is.na(data.fit2$log_ba),'log_ba'] <- 0
# 
# auc(data.fit2[data.fit2$year.idx<=8,'z'], data.fit2[data.fit2$year.idx<=8,'score_z'], quiet = TRUE)
# auc(data.fit2[data.fit2$year.idx>8,'z'], data.fit2[data.fit2$year.idx>8,'score_z'], quiet = TRUE)
# 
# sum((data.fit2[data.fit2$year.idx<=8 & data.fit2$y>0,'log_ba']-data.fit2[data.fit2$year.idx<=8 & data.fit2$y>0,'score_ba'])^2)/sum(data.fit2$year.idx<=8)
# sum((data.fit2[data.fit2$year.idx>8 & data.fit2$y>0,'log_ba']-data.fit2[data.fit2$year.idx>8 & data.fit2$y>0,'score_ba'])^2)/sum(data.fit2$year.idx>8)
# 
# 
# sum((data.fit2[data.fit2$year.idx<=8 & data.fit2$y>0,'y']-data.fit2[data.fit2$year.idx<=8 & data.fit2$y>0,'score_cnt'])^2)/sum(data.fit2$year.idx<=8)
# sum((data.fit2[data.fit2$year.idx>8 & data.fit2$y>0,'y']-data.fit2[data.fit2$year.idx>8 & data.fit2$y>0,'score_cnt'])^2)/sum(data.fit2$year.idx>8)
# 
# 
# data.fit2 <- data.fit2[order(data.fit2$time.idx),]


# data.fit2$score_cnt_bin <- cut(data.fit2$score_cnt, 4,lables=1:4)


# save(data.fit2, file=file.path(dir.out, 'data.fit2.score_0.0625.RData'))
# load(file.path(dir.out, 'data.fit2.score_0.0625.RData'))


# save(data.fit2, file=file.path(dir.out, 'data.fit2.score_0.125.RData'))
load(file.path(dir.out, 'data.fit2.score_0.125.RData'))

# load(file.path(dir.out, 'data.fit2.score_0.25.RData'))



table(cut(data.fit2[data.fit2$y>0,'score_ba'], breaks=8))
data.fit2$score_ba_bin <- cut(data.fit2$score_ba, breaks=c(0,3.19,3.64,4.09,4.54,4.99,100),labels=1:6)
table((data.fit2[data.fit2$y>0,'score_ba_bin']),useNA='always')


table(cut(data.fit2[data.fit2$y>0,'score_cnt'], breaks=8))
data.fit2$score_cnt_bin <- cut(data.fit2$score_cnt, breaks=c(0, 1.03,1.12,1.21,100),labels=1:4)
table((data.fit2[data.fit2$y>0,'score_cnt_bin']),useNA='always')


data.fit2$score_z_bin <- cut(data.fit2$score_z, breaks=c(seq(0,0.8,0.1),1),labels=1:9)
table((data.fit2[,'score_z_bin']),useNA='always')
####################################################################
# INLA
####################################################################

data.fit3 <- reshape(data.fit2[,c('grid.idx','time.idx','y')],
                     timevar = "time.idx",
                     idvar = "grid.idx",
                     direction = "wide"
)

B2.merge <- merge(B2, data.fit3, by = "grid.idx")

save(B2.merge, file=file.path(dir.out, 'grid_cell_map_0.125.RData'))

B2.adj <- poly2nb(B2)
nb2INLA("map.adj", B2.adj)
g <- inla.read.graph(filename = "map.adj")


z <- as.vector((data.fit2$y>0)+0)
log.ba <- as.vector(ifelse(data.fit2$y>0, data.fit2$log_ba, NA))
# cnt = as.vector(data.fit2$y)
cnt <- as.vector(ifelse(data.fit2$y>0, data.fit2$y, NA))
#

#prepare for prediction
z[which(data.fit2$time.idx>=97)] <- NA
log.ba[which(data.fit2$time.idx>=97)] <- NA
cnt[which(data.fit2$time.idx>=97)] <- NA





coo <- as.matrix(data.fit2[,c('grid.cent.x.utm','grid.cent.y.utm')])

# mesh.spat <- inla.mesh.2d(
#   # loc = coo,
#   boundary = domainSP,
#   offset = c(20, 30),
#   cutoff = 1, max.edge = c(40, 60)
# )

mesh.spat <- inla.mesh.2d(
  # loc = coo,
  boundary = domainSP,
  offset = c(20, 30),
  cutoff = 1, max.edge = c(20, 30)
)

plot(mesh.spat)
points(coo, col = "red",pch=19,cex=0.1)

mesh.spat$n

spde.spat <- inla.spde2.pcmatern(mesh = mesh.spat,
                                 # PC-prior on range: P(practic.range < 0.05) = 0.01
                                 prior.range = c(60, 0.05),
                                 # PC-prior on sigma: P(sigma > 1) = 0.01
                                 prior.sigma = c(1, 0.01))


A.spat <- inla.spde.make.A(mesh = mesh.spat, loc = coo)




resolution <- "0.125"


if (resolution=='0.125'){
  #SPDE model
  # mesh_score_1 <- inla.mesh.1d(c(0.93, 0.99, 1.03,  1.06, 1.5 , 1.8),boundary=c('free','free'))
  #BYM2 model
  mesh_score_1 <- inla.mesh.1d(c(1,1.05, 1.1, 1.15,1.2, 1.25 ),boundary=c('free','free')) 
  A1 <- inla.spde.make.A(mesh_score_1, loc=data.fit2$score_cnt)
  spde_score_1 <-  inla.spde2.pcmatern(mesh_score_1, 
                                       prior.range = c(0.1, 0.05),
                                       prior.sigma = c(1, 0.05))
  # spde_score_1 <-  inla.spde2.matern(mesh_score_1, constr = TRUE)
  
  spde_score_1.idx <- inla.spde.make.index("score_1", n.spde = spde_score_1$n.spde)
  
  mesh_score_2 <- inla.mesh.1d(c(0.01, 0.07,0.2,0.3,0.5,0.8),boundary=c('free','free'))
  # mesh_score_2 <- inla.mesh.1d(seq(0, 1, by = 0.2),boundary=c('dirichlet','dirichlet'))
  A2 <- inla.spde.make.A(mesh_score_2, loc=data.fit2$score_z)
  spde_score_2 <-  inla.spde2.pcmatern(mesh_score_2, 
                                       prior.range = c(0.2, 0.05),
                                       prior.sigma = c(1, 0.05))
  # spde_score_2 <-  inla.spde2.matern(mesh_score_2, constr = TRUE)
  spde_score_2.idx <- inla.spde.make.index("score_2", n.spde = spde_score_2$n.spde)
  
  
  mesh_score_3 <- inla.mesh.1d(c(3.14, 3.6, 4, 4.4, 4.8),boundary=c('free','free'))
  # mesh_score_3 <- inla.mesh.1d(seq(0, 7, by = 1),boundary=c('dirichlet','free'))
  A3 <- inla.spde.make.A(mesh_score_3, loc=data.fit2$score_ba)
  spde_score_3 <-  inla.spde2.pcmatern(mesh_score_3, 
                                       prior.range = c(0.5, 0.05),
                                       prior.sigma = c(1, 0.1))
  # spde_score_3 <-  inla.spde2.matern(mesh_score_3, constr = TRUE)
  spde_score_3.idx <- inla.spde.make.index("score_3", n.spde = spde_score_3$n.spde)
}

if (resolution =='0.25'){
  mesh_score_1 <- inla.mesh.1d(c(0.9,1.1,1.2, 1.3,1.6,1.9),boundary=c('free','free')) 
  A1 <- inla.spde.make.A(mesh_score_1, loc=data.fit2$score_cnt)
  spde_score_1 <-  inla.spde2.pcmatern(mesh_score_1, 
                                       prior.range = c(0.3, 0.3),
                                       prior.sigma = c(1, 0.05))
  # spde_score_1 <-  inla.spde2.matern(mesh_score_1, constr = TRUE)
  
  spde_score_1.idx <- inla.spde.make.index("score_1", n.spde = spde_score_1$n.spde)
  
  mesh_score_2 <- inla.mesh.1d(c(0.03,0.05, 0.1,0.2,0.4,0.6,0.8),boundary=c('free','free'))
  # mesh_score_2 <- inla.mesh.1d(seq(0, 1, by = 0.2),boundary=c('dirichlet','dirichlet'))
  A2 <- inla.spde.make.A(mesh_score_2, loc=data.fit2$score_z)
  spde_score_2 <-  inla.spde2.pcmatern(mesh_score_2, 
                                       prior.range = c(0.2, 0.3),
                                       prior.sigma = c(1, 0.05))
  # spde_score_2 <-  inla.spde2.matern(mesh_score_2, constr = TRUE)
  spde_score_2.idx <- inla.spde.make.index("score_2", n.spde = spde_score_2$n.spde)
  
  
  mesh_score_3 <- inla.mesh.1d(c(2.2, 3, 3.4, 4, 5),boundary=c('free','free'))
  # mesh_score_3 <- inla.mesh.1d(seq(0, 7, by = 1),boundary=c('dirichlet','free'))
  A3 <- inla.spde.make.A(mesh_score_3, loc=data.fit2$score_ba)
  spde_score_3 <-  inla.spde2.pcmatern(mesh_score_3, 
                                       prior.range = c(1, 0.3),
                                       prior.sigma = c(1, 0.05))
  # spde_score_3 <-  inla.spde2.matern(mesh_score_3, constr = TRUE)
  spde_score_3.idx <- inla.spde.make.index("score_3", n.spde = spde_score_3$n.spde)
  
}




cnt.stack <- inla.stack(
  data= list(Y.log=cbind(cnt,NA,NA)),
  A <- list(1,1,1,1,A1, A.spat, A.spat,1),
  effect = list(Intercept1=rep(1,nrow(data.fit2)), idarea1=data.fit2$grid.idx, 
                idarea.cnt=data.fit2$grid.idx, time.idx1=data.fit2$month.idx, score_1=spde_score_1.idx,
                spat.share.cnt=1:spde.spat$n.spde, spat.cnt=1:spde.spat$n.spde, score_cnt_bin=data.fit2$score_cnt_bin  ),
  tag='cnt'
)

z.stack <- inla.stack(
  data= list(Y.log=cbind(NA,z,NA)),
  A <- list(1,1,1,A2, A.spat,1),
  effect = list(Intercept2=rep(1,nrow(data.fit2)), idarea2=data.fit2$grid.idx, time.idx2=data.fit2$month.idx, score_2=spde_score_2.idx,
                spat.z=1:spde.spat$n.spde, score_z_bin = data.fit2$score_z_bin  ),
  tag='z'
)


ba.stack <- inla.stack(
  data= list(Y.log=cbind(NA,NA,log.ba)),
  A <- list(1,1,1,1, A3, A.spat, A.spat,1),
  effect = list(Intercept3=rep(1,nrow(data.fit2)), idarea3=data.fit2$grid.idx,
                idarea.ba=data.fit2$grid.idx, time.idx3=data.fit2$month.idx, score_3=spde_score_3.idx,
                spat.share.ba=1:spde.spat$n.spde, spat.ba = 1:spde.spat$n.spde, score_ba_bin=data.fit2$score_ba_bin),
  tag='ba'
)


all.stack <- inla.stack(cnt.stack, z.stack, ba.stack )



########-----------------------Prepare prediction data------------------------
grd <- SpatialGrid(GridTopology(cellcentre.offset = c(-10, 36),
                                cellsize = c(0.125,0.125),
                                cells.dim = c(34, 58)))


grd_sp <- SpatialPixelsDataFrame(points = grd, data = data.frame(grid.idx = 1:length(grd)),proj4string = CRS("+proj=longlat +datum=WGS84"))

grd_poly <- as(grd_sp, 'SpatialPolygonsDataFrame')
grd_poly <- spTransform(grd_poly, CRS(projutm))
grd_poly$lon.grid <- grd_sp@coords[,1]
grd_poly$lat.grid <- grd_sp@coords[,2]

B2.pred <- as_Spatial(st_intersection(st_as_sf(grd_poly),st_as_sf(B)))
coord.cent <- st_coordinates(st_centroid(st_as_sf(B2.pred)))
B2.pred$grid.cent.x.utm <- coord.cent[,1]
B2.pred$grid.cent.y.utm <- coord.cent[,2]
B2.pred$E <- area(B2.pred)

B2.pred$grid.idx <- 1:nrow(B2.pred)


ggplot(st_as_sf(B2.pred)) + geom_sf(aes(fill=grid.idx))




data.merge.pred <- B2.pred |>  
  st_as_sf() |> # cast to sf
  mutate(grid_id = row_number()) |> # create unique ID
  st_join(loc.data.utm) |> # join the species dataset
  group_by(grid_id)


data.fit.ba.pred <- data.merge.pred %>% st_drop_geometry() %>% filter( time.idx<=108, length >=24*60 )%>% 
  group_by(grid.idx,grid_id, year.idx, month.idx, time.idx) %>%
  summarise(area_ha = sum(area_ha), 
            log_ba = log(area_ha),
            y = n(),
            lon.grid = mean(lon.grid),
            lat.grid = mean(lat.grid),
            x.utm = mean(grid.cent.x.utm),
            y.utm = mean(grid.cent.y.utm),
            E = mean(E))



data.fit2.pred <- do.call(rbind, lapply(97:108, function(x) {B2.pred@data$time.idx = x
return(B2.pred@data)}))
print(dim(data.fit2.pred))
data.fit2.pred <- merge(data.fit2.pred, data.fit.ba.pred[,c('grid.idx','time.idx','year.idx', 'month.idx','y','area_ha','log_ba')],
                        by=c('grid.idx','time.idx'),all.x=T)


print(dim(data.fit2.pred))
data.fit2.pred[is.na(data.fit2.pred$y),'y'] <- 0
summary(data.fit2.pred$y)

data.fit2.pred$month.idx <- (data.fit2.pred$time.idx-1)%%12 + 1
data.fit2.pred$year.idx <- (data.fit2.pred$time.idx-1)%/%12 + 1
data.fit2.pred$E1 <- 1

data.fit2.pred.sf <- st_as_sf(data.fit2.pred, coords=c('grid.cent.x.utm','grid.cent.y.utm'),crs='+proj=utm +zone=29 +datum=WGS84 +units=km +no_defs')
data.fit2.sf <- st_as_sf(data.fit2, coords=c('grid.cent.x.utm','grid.cent.y.utm'),crs='+proj=utm +zone=29 +datum=WGS84 +units=km +no_defs')


for (time.idx in 97:108){
  data.fit2.pred[data.fit2.pred$time.idx==time.idx,
                 c('score_z','score_ba','score_cnt',
                   'score_z_bin','score_ba_bin','score_cnt_bin')] <- st_join(data.fit2.pred.sf[data.fit2.pred.sf$time.idx==time.idx,], 
                                                                 data.fit2.sf[data.fit2.sf$time.idx==time.idx,c('score_z','score_ba','score_cnt',
                                                                                                                'score_z_bin','score_ba_bin','score_cnt_bin')], 
                                                                 join = st_nearest_feature) |> st_drop_geometry()|> dplyr::select(score_z,score_ba,score_cnt,
                                                                                                                                  score_z_bin,score_ba_bin,score_cnt_bin)
}

data.fit2.pred <- data.fit2.pred[order(data.fit2.pred$time.idx),]

data.fit2.pred[is.na(data.fit2.pred$log_ba),'log_ba'] <- 0

# ggplot(data.fit2.sf[data.fit2.sf$time.idx==104,]) + geom_sf(aes(fill=score_z,color=score_z))
# ggplot(data.fit2.pred.sf[data.fit2.pred.sf$time.idx==104,]) + geom_sf(aes(fill=score_z,color=score_z))
# 

coo.pred <- as.matrix(data.fit2.pred[,c('grid.cent.x.utm','grid.cent.y.utm')])


A.spat.1.pred <- inla.spde.make.A(mesh = mesh.spat, loc = coo.pred)

A.spat.2.pred <- inla.spde.make.A(mesh = mesh.spat, loc = coo.pred)

A.spat.3.pred <- inla.spde.make.A(mesh = mesh.spat, loc = coo.pred)


A1.pred <- inla.spde.make.A(mesh_score_1, loc=data.fit2.pred$score_cnt)

A2.pred <- inla.spde.make.A(mesh_score_2, loc=data.fit2.pred$score_z)

A3.pred <- inla.spde.make.A(mesh_score_3, loc=data.fit2.pred$score_ba)


data.fit3.pred <- reshape(data.fit2.pred[,c('grid.idx','time.idx','y')],
                          timevar = "time.idx",
                          idvar = "grid.idx",
                          direction = "wide"
)

B2.merge.pred <- merge(B2.pred, data.fit3.pred, by = "grid.idx")

# save(B2.merge.pred, file=file.path(dir.out, 'grid_cell_map_0.0625.RData'))
# save(data.fit2.pred, A.spat.1.pred, A.spat.2.pred, 
#      A.spat.3.pred, A1.pred, A2.pred, A3.pred,
#      file=file.path(dir.out, 'data.fit2.pred_0.0625.RData'))

# save(B2.merge.pred, file=file.path(dir.out, 'grid_cell_map_0.125.RData'))
# save(data.fit2.pred, A.spat.1.pred, A.spat.2.pred,
#      A.spat.3.pred, A1.pred, A2.pred, A3.pred,
#      file=file.path(dir.out, 'data.fit2.pred_0.125.RData'))
#--------------------------------------------------------------


n1 <- dim(data.fit2)[1]


hyper.rw <-  list(prec = list(prior="loggamma",param=c(1,1)))
hyper.time.rw <- list(prec = list(prior="loggamma",param=c(80,20)))


hyper.bym2 = list (
  phi = list(
  prior = "pc",
  param = c(0.5 , 2/3) ,
  initial = 3) ,
  # prec = list(
  #   prior = "pc. prec",
  #   param = c (3 , 0.1) ,
  #   initial = 5)
  prec = list(prior="loggamma",param=c(1,0.5))
  )
# 
# formula <- Y.log ~  -1  +
#   Intercept1 + f(spat.cnt, copy='spat.z',fixed=F) +  f(time.idx1, model='rw1',hyper=hyper.rw) + f(score_1, model=spde_score_1 )+
#   Intercept2 + f(spat.z, model=spde.spat) +  f(time.idx2, model='rw1',hyper=hyper.rw) + f(score_2, model=spde_score_2 )+
#   Intercept3 + f(spat.ba, copy='spat.z',fixed=F) +  f(time.idx3, model='rw1',hyper=hyper.rw) +f(score_3, model=spde_score_3 )

formula <- Y.log ~  -1  +
  Intercept1 + f(spat.cnt, copy='spat.z',fixed=F)  + f(score_1, model=spde_score_1 )+
  Intercept2 + f(spat.z, model=spde.spat)  + f(score_2, model=spde_score_2 )+
  Intercept3 + f(spat.ba, copy='spat.z',fixed=F)  +f(score_3, model=spde_score_3 )


# formula <- Y.log ~  -1  +
#   Intercept1 + f(spat.cnt, copy='spat.z',fixed=F) +  f(time.idx1, model='rw1',hyper=hyper.time.rw) + f(score_cnt_bin, model='rw1',hyper=hyper.rw )+
#   Intercept2 + f(spat.z, model=spde.spat) +  f(time.idx2, model='rw1',hyper=hyper.time.rw) + f(score_z_bin, model='rw1',hyper=hyper.rw )+
#   Intercept3 + f(spat.ba, copy='spat.z',fixed=F) +  f(time.idx3, model='rw1',hyper=hyper.time.rw) + f(score_ba_bin, model='rw1',hyper=hyper.rw )

# formula <- Y.log ~  -1  +
#   Intercept1 + f(spat.cnt, copy='spat.z',fixed=F) +   f(score_cnt_bin, model='rw1',hyper=hyper.rw )+
#   Intercept2 + f(spat.z, model=spde.spat) + f(score_z_bin, model='rw1',hyper=hyper.rw )+
#   Intercept3 + f(spat.ba, copy='spat.z',fixed=F) + f(score_ba_bin, model='rw1',hyper=hyper.rw )


# formula <- Y.log ~  -1  +
#   Intercept1 +  f(idarea1, copy='idarea2',fixed=F) +  f(time.idx1, model='rw1',hyper=hyper.rw) + f(score_1, model=spde_score_1 )+
#   Intercept2 +  f(idarea2, model='bym2', graph=g) +  f(time.idx2, model='rw1',hyper=hyper.rw) + f(score_2, model=spde_score_2 )+
#   Intercept3 + f(idarea3,  copy='idarea2',fixed=F) +  f(time.idx3, model='rw1',hyper=hyper.rw) +f(score_3, model=spde_score_3 )

# formula <- Y.log ~  -1  +
#   Intercept1 + f(spat.share.cnt, copy='spat.z',fixed=F) + f(spat.cnt, model=spde.spat) +
#   f(time.idx1, model='rw1',hyper=hyper.rw) + f(score_1, model=spde_score_1 )+
# 
#   Intercept2 + f(spat.z, model=spde.spat) +  f(time.idx2, model='rw1',hyper=hyper.rw) + f(score_2, model=spde_score_2 )+
# 
#   Intercept3 + f(spat.share.ba, copy='spat.z',fixed=F) +  f(spat.ba, model=spde.spat)  +
#     f(time.idx3, model='rw1',hyper=hyper.rw) +f(score_3, model=spde_score_3 )

# formula <- Y.log ~  -1  +
#   Intercept1 +  f(idarea1, copy='idarea2',fixed=F) +  f(time.idx1, model='rw1',hyper=hyper.rw) + f(score_1, model=spde_score_1 )+
#   Intercept2 +  f(idarea2, model='bym2', graph=g) +  f(time.idx2, model='rw1', hyper=hyper.rw) + f(score_2, model=spde_score_2 )+
#   Intercept3 + f(idarea3,  copy='idarea2',fixed=F) + f(idarea.ba, model='bym2', graph=g, hyper=hyper.bym2)+ f(time.idx3, model='rw1',hyper=hyper.rw) +f(score_3, model=spde_score_3 )

# formula <- Y.log ~  -1  +
#   Intercept1 +  f(idarea1, copy='idarea2',fixed=F) +  f(time.idx1, model='rw1',hyper=hyper.rw) +
#   Intercept2 +  f(idarea2, model='bym2', graph=g) +  f(time.idx2, model='rw1', hyper=hyper.rw) + f(score_2, model=spde_score_2 )+
#   Intercept3 + f(idarea3,  copy='idarea2',fixed=F) + f(idarea.ba, model='bym2', graph=g, hyper=hyper.bym2)+ f(time.idx3, model='rw1',hyper=hyper.rw)

# formula <- Y.log ~  -1  +
#   Intercept1 +  f(idarea1, copy='idarea2',fixed=F) + f(time.idx1, model='rw1',hyper=hyper.rw) + f(score_cnt_bin, model='rw1',hyper=hyper.rw )+
#   Intercept2 +  f(idarea2, model='bym2', graph=g) +  f(time.idx2, model='rw1', hyper=hyper.rw) + f(score_z_bin, model='rw1' )+
#   Intercept3 + f(idarea3,  copy='idarea2',fixed=F) + f(time.idx3, model='rw1',hyper=hyper.rw) +f(score_ba_bin, model='rw1' ,hyper=hyper.rw)


t1 <- Sys.time()
res <- inla(formula,
               family = c('poisson','binomial', 'gamma'), data = inla.stack.data(all.stack),  Ntrials=1,
               control.predictor = list(A = inla.stack.A(all.stack), compute = TRUE, link=c(rep(1,n1),rep(2,n1),rep(3,n1))),
               verbose=TRUE,
               control.compute=list(config = TRUE),
               control.family = list( list(), list(), list()),
               control.fixed = list(expand.factor.strategy = 'inla')
)
# res <- inla(formula,
#                family = c('poisson','binomial', 'weibull'), data = inla.stack.data(all.stack),  Ntrials=1,
#                control.predictor = list(A = inla.stack.A(all.stack), compute = TRUE, link=c(rep(1,n1),rep(2,n1),rep(3,n1))),
#                verbose=TRUE,
#                control.compute=list(config = TRUE),
#                control.family = list( list(), list(), list()),
#                control.fixed = list(expand.factor.strategy = 'inla')
# )
# res <- inla(formula,
#             family = c('poisson','binomial', 'egp'), data = inla.stack.data(all.stack),  Ntrials=1,
#             control.predictor = list(A = inla.stack.A(all.stack), compute = TRUE, link=c(rep(1,n1),rep(2,n1),rep(3,n1))),
#             verbose=TRUE,
#             control.compute=list(config = TRUE),
#             control.family = list( list(), list(), list(control.link = list(quantile = 0.5),
#                                                         hyper = list(
#                                                           tail = list(
#                                                           ##initial = xi.intern,
#                                                           fixed = !TRUE,
#                                                           prior = "pc.egptail",
#                                                           param = c(5, -1, 1)),
#                                                           shape = list(
#                                                             ##initial = kappa.intern,
#                                                             fixed = !TRUE,
#                                                             prior = "loggamma",
#                                                             param = c(100, 100)
#                                                             )
#                                                           )
#                                                         )),
#             control.fixed = list(expand.factor.strategy = 'inla')
# )
t2 <- Sys.time()
print(t2-t1)
summary(res)

# save(res, file=file.path(dir.out, 'Model_weibull_0.0625.RData'))
# save(res, file=file.path(dir.out, 'Model_egp_spde_0.125.RData'))
# save(res, file=file.path(dir.out, 'Model_weibull_bym2_0.125.RData'))
# save(res, file=file.path(dir.out, 'Model_weibull_spde_0.125.RData'))
save(res, file=file.path(dir.out, 'Model_gamma_spde_0.125.RData'))



