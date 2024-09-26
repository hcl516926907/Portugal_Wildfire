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

data.rf <- data.fit2[data.fit2$time.idx <= 84, c(covar.names,'grid.idx','time.idx','year.idx','area_ha','z','y')]
# for (var in c('LVegTyp','HVegTyp','month.idx')){
#   data.rf[,var] <- as.factor(data.rf[,var])
# }

library(xgboost)
# library(caret)
library(pROC)
data.rf$log_ba <- log(data.rf$area_ha)
data.rf[is.na(data.rf$log_ba),'log_ba'] <- 0

kfold_cv_pred  <- function(data, target.name ,covar.names, params, k) {
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
    model <- xgb.train(params = params, data = dtrain, nrounds = params$nrounds)
    
    train_pred <- predict(model, dtrain)
    test_pred <- predict(model, dtest)
    
    data[index,'score'] <- test_pred
  }
  
  
  return(data$score)
}

params_z <- list(
  nrounds = 55,
  eta = 0.05,
  max_depth = 4,
  gamma = 0.1,
  scale_pos_weight = 25,
  colsample_bytree = 0.8,
  min_child_weight = 1,
  subsample = 0.8,
  objective = "binary:logistic",
  eval_metric = "auc",
  scale_pos_weight = 25,
  verbosity = 0
)

data.rf$score_z <- kfold_cv_pred(data.rf, 'z', covar.names ,params_z, 7)


TrainSet1 <- xgb.DMatrix(data=as.matrix(data.rf[, c(covar.names)]),label=data.rf[, c('z')])
set.seed(1234)
model.z <- xgb.train(params = params_z, data = TrainSet1, nrounds = params_z$nrounds)


data.fit2[data.fit2$time.idx <= 84, 'score_z'] <- data.rf$score_z
data.fit2[data.fit2$time.idx > 84, 'score_z'] <- predict(model.z, as.matrix(data.fit2[data.fit2$time.idx > 84,c(covar.names)]))



# with coordiantes
params_ba <- list(
  nrounds = 80,
  eta = 0.04,
  max_depth = 2,
  gamma = 0,
  colsample_bytree = 1,
  min_child_weight = 8,
  subsample = 0.8,
  objective=c('reg:squarederror')
)


kfold_cv_pred_pos <- function(data, target.name ,covar.names, params, k) {
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


data.rf$score_ba <- kfold_cv_pred_pos(data.rf,'log_ba',covar.names, params_ba, 7)

TrainSet2 <- xgb.DMatrix(data=as.matrix(data.rf[data.rf$y>0, c(covar.names)]),label=data.rf[data.rf$y>0, c('log_ba')])
# 
set.seed(1234)
model.ba <- xgb.train(params = params_ba, data = TrainSet2, nrounds = params_ba$nrounds)


# data.fit2[data.fit2$time.idx <= 84, 'score_ba'] <- data.rf$score_ba
# data.fit2[data.fit2$time.idx > 84, 'score_ba'] <- predict(model.ba, as.matrix(data.fit2[data.fit2$time.idx > 84,c(covar.names,'score_z')]))

data.fit2[data.fit2$time.idx <= 84, 'score_ba'] <- data.rf$score_ba
data.fit2[data.fit2$time.idx > 84, 'score_ba'] <- predict(model.ba, as.matrix(data.fit2[data.fit2$time.idx > 84,c(covar.names)]))



params_cnt <- list(
  nrounds = 80,
  eta = 0.04,
  max_depth = 3,
  gamma = 0.4,
  colsample_bytree = 0.4,
  min_child_weight = 8,
  subsample = 1,
  objective=c('count:poisson')
)

data.rf$score_cnt <- kfold_cv_pred_pos(data.rf,'y',covar.names, params_ba, 7)


TrainSet3 <- xgb.DMatrix(data=as.matrix(data.rf[data.rf$y>0, c(covar.names)]),label=data.rf[data.rf$y>0, c('y')])

set.seed(1234)
model.cnt <- xgb.train(params = params_cnt, data = TrainSet3, nrounds = params_cnt$nrounds)


data.fit2[data.fit2$time.idx <= 84, 'score_cnt'] <- data.rf$score_cnt
data.fit2[data.fit2$time.idx > 84 , 'score_cnt'] <- predict(model.cnt, as.matrix(data.fit2[data.fit2$time.idx > 84 ,c(covar.names)]))




data.fit2[is.na(data.fit2$log_ba),'log_ba'] <- 0

auc(data.fit2[data.fit2$year.idx<=7,'z'], data.fit2[data.fit2$year.idx<=7,'score_z'], quiet = TRUE)
auc(data.fit2[data.fit2$year.idx>7,'z'], data.fit2[data.fit2$year.idx>7,'score_z'], quiet = TRUE)

# sum((data.fit2[data.fit2$year.idx<=7,'log_ba']-data.fit2[data.fit2$year.idx<=7,'score_ba'])^2)/sum(data.fit2$year.idx<=7)
# sum((data.fit2[data.fit2$year.idx>7,'log_ba']-data.fit2[data.fit2$year.idx>7,'score_ba'])^2)/sum(data.fit2$year.idx>7)

sum((data.fit2[data.fit2$year.idx<=7 & data.fit2$y>0,'log_ba']-data.fit2[data.fit2$year.idx<=7 & data.fit2$y>0,'score_ba'])^2)/sum(data.fit2$year.idx<=7)
sum((data.fit2[data.fit2$year.idx>7 & data.fit2$y>0,'log_ba']-data.fit2[data.fit2$year.idx>7 & data.fit2$y>0,'score_ba'])^2)/sum(data.fit2$year.idx>7)


sum((data.fit2[data.fit2$year.idx<=7 & data.fit2$y>0,'y']-data.fit2[data.fit2$year.idx<=7 & data.fit2$y>0,'score_cnt'])^2)/sum(data.fit2$year.idx<=7)
sum((data.fit2[data.fit2$year.idx>7 & data.fit2$y>0,'y']-data.fit2[data.fit2$year.idx>7 & data.fit2$y>0,'score_cnt'])^2)/sum(data.fit2$year.idx>7)




data.fit3 <- reshape(data.fit2[,c('grid.idx','time.idx','y')],
                     timevar = "time.idx",
                     idvar = "grid.idx",
                     direction = "wide"
)

B2.merge <- merge(B2, data.fit3, by = "grid.idx")

save(B2.merge, file=file.path(dir.out, 'grid_cell_map.RData'))

B2.adj <- poly2nb(B2)
nb2INLA("map.adj", B2.adj)
g <- inla.read.graph(filename = "map.adj")



n1 <- dim(data.fit2)[1]
#
nothing1 <- rep(NA, n1)
nothing2 <- rep(NA, n1)

data.fit2 <- data.fit2[order(data.fit2$time.idx),]


coo <- as.matrix(data.fit2[,c('grid.cent.x.utm','grid.cent.y.utm')])
# mesh.spat <- inla.mesh.2d(
#   loc = coo, offset = c(20, 30),
#   cutoff = 1, max.edge = c(60, 80)
# )

mesh.spat <- inla.mesh.2d(
  # loc = coo,
  boundary = domainSP,
  offset = c(20, 30),
  cutoff = 1, max.edge = c(40, 60)
)

plot(mesh.spat)
points(coo, col = "red")

mesh.spat$n

spde.spat <- inla.spde2.pcmatern(mesh = mesh.spat,
                            # PC-prior on range: P(practic.range < 0.05) = 0.01
                            prior.range = c(60, 0.05),
                            # PC-prior on sigma: P(sigma > 1) = 0.01
                            prior.sigma = c(1, 0.01)) 

indexs.1 <- inla.spde.make.index("spat1", spde.spat$n.spde)
A.spat.1 <- inla.spde.make.A(mesh = mesh.spat, loc = coo)

indexs.2 <- inla.spde.make.index("spat2", spde.spat$n.spde)
A.spat.2 <- inla.spde.make.A(mesh = mesh.spat, loc = coo)

indexs.3 <- inla.spde.make.index("spat3", spde.spat$n.spde)
A.spat.3 <- inla.spde.make.A(mesh = mesh.spat, loc = coo)

z <- as.vector((data.fit2$y>0)+0)
log.ba <- as.vector(ifelse(data.fit2$y>0, data.fit2$log_ba, NA))
# cnt = as.vector(data.fit2$y)
cnt <- as.vector(ifelse(data.fit2$y>0, data.fit2$y, NA))
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
mesh_score_1 <- inla.mesh.1d(c(0.9,1.1,1.3,1.6,1.9),boundary=c('free','free')) 
A1 <- inla.spde.make.A(mesh_score_1, loc=data.fit2$score_cnt)
spde_score_1 <-  inla.spde2.pcmatern(mesh_score_1, 
                                     prior.range = c(0.3, 0.05),
                                     prior.sigma = c(1, 0.05))
# spde_score_1 <-  inla.spde2.matern(mesh_score_1, constr = TRUE)

spde_score_1.idx <- inla.spde.make.index("score_1", n.spde = spde_score_1$n.spde)

mesh_score_2 <- inla.mesh.1d(c(0.1,0.2,0.4,0.6,0.8),boundary=c('free','free'))
# mesh_score_2 <- inla.mesh.1d(seq(0, 1, by = 0.2),boundary=c('dirichlet','dirichlet'))
A2 <- inla.spde.make.A(mesh_score_2, loc=data.fit2$score_z)
spde_score_2 <-  inla.spde2.pcmatern(mesh_score_2, 
                                     prior.range = c(0.2, 0.05),
                                     prior.sigma = c(1, 0.05))
# spde_score_2 <-  inla.spde2.matern(mesh_score_2, constr = TRUE)
spde_score_2.idx <- inla.spde.make.index("score_2", n.spde = spde_score_2$n.spde)


mesh_score_3 <- inla.mesh.1d(c(3,4,5,6),boundary=c('free','free'))
# mesh_score_3 <- inla.mesh.1d(seq(0, 7, by = 1),boundary=c('dirichlet','free'))
A3 <- inla.spde.make.A(mesh_score_3, loc=data.fit2$score_ba)
spde_score_3 <-  inla.spde2.pcmatern(mesh_score_3, 
                                     prior.range = c(0.5, 0.05),
                                     prior.sigma = c(1, 0.05))
# spde_score_3 <-  inla.spde2.matern(mesh_score_3, constr = TRUE)
spde_score_3.idx <- inla.spde.make.index("score_3", n.spde = spde_score_3$n.spde)



cnt.stack <- inla.stack(
  data= list(Y.log=cbind(cnt,NA,NA)),
  A <- list(1,1,1,A1,1,A.spat.1),
  effect = list(Intercept1=rep(1,nrow(data.fit2)), idarea1=data.fit2$grid.idx, time.idx1=data.fit2$month.idx, score_1=spde_score_1.idx,
                score_1_grp = score_cnt_grp, spat1=indexs.1),
  tag='cnt'
)

dim(cnt.stack$A)
1 + 192 + 12 + ncol(A1) + length(table(score_cnt_grp))

z.stack <- inla.stack(
  data= list(Y.log=cbind(NA,z,NA)),
  A <- list(1,1,1,A2,1,A.spat.2),
  effect = list(Intercept2=rep(1,nrow(data.fit2)), idarea2=data.fit2$grid.idx, time.idx2=data.fit2$month.idx, score_2=spde_score_2.idx,
                score_2_grp=score_z_grp, spat2=indexs.2),
  tag='z'
)

dim(z.stack$A)
1 + 192 + 12 + ncol(A2) + length(table(score_z_grp))


ba.stack <- inla.stack(
  data= list(Y.log=cbind(NA,NA,log.ba)),
  A <- list(1,1,1, A3,1,A.spat.3),
  effect = list(Intercept3=rep(1,nrow(data.fit2)), idarea3=data.fit2$grid.idx, time.idx3=data.fit2$month.idx, score_3=spde_score_3.idx,
                score_3_grp=score_ba_grp, spat3=indexs.3),
  tag='ba'
)

dim(ba.stack$A)
1 + 192 + 12 + ncol(A3) + length(table(score_ba_grp))


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

formula3 <- Y.log ~  -1  +
  Intercept1 + f(spat1, copy='spat2',fixed=F) +  f(time.idx1, model='rw1',hyper=hyper.rw) + f(score_1, model=spde_score_1 )+
  Intercept2 + f(spat2, model=spde.spat) +  f(time.idx2, model='rw1',hyper=hyper.rw) + f(score_2, model=spde_score_2 )+
  Intercept3 + f(spat3, copy='spat2',fixed=F) +  f(time.idx3, model='rw1',hyper=hyper.rw) +f(score_3, model=spde_score_3 )


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
save(res1.1, file=file.path(dir.out,'Final_Model_1.1.RData'))


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
save(res2.1, file=file.path(dir.out,'Final_Model_2.1.RData'))

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

# hurdle poisson and weibull
t1 <- Sys.time()
res3.2 <- inla(formula2,
               family = c('poisson','binomial', 'weibull'), data = inla.stack.data(all.stack),  Ntrials=1,
               control.predictor = list(A = inla.stack.A(all.stack), compute = TRUE, link=c(rep(1,n1),rep(2,n1),rep(3,n1))),
               verbose=TRUE,
               control.compute=list(config = TRUE),
               control.family = list( list(), list(), list()),
               control.fixed = list(expand.factor.strategy = 'inla')
)
t2 <- Sys.time()
print(t2-t1)
summary(res3.2)

save(res3.2, file=file.path(dir.out,'Final_Model_3.2.RData'))

# load(file.path(dir.out,'Final_Model_3.2.RData'))

t1 <- Sys.time()
res3.3 <- inla(formula3,
               family = c('poisson','binomial', 'weibull'), data = inla.stack.data(all.stack),  Ntrials=1,
               control.predictor = list(A = inla.stack.A(all.stack), compute = TRUE, link=c(rep(1,n1),rep(2,n1),rep(3,n1))),
               verbose=TRUE,
               control.compute=list(config = TRUE),
               control.family = list( list(), list(), list()),
               control.fixed = list(expand.factor.strategy = 'inla')
)
t2 <- Sys.time()
print(t2-t1)
summary(res3.3)

save(res3.3, file=file.path(dir.out,'Final_Model_3.3.RData'))


n1 <- 192*108



load(file=file.path(dir.out,'Final_Model_3.1_pred.sp_200.RData'))
# load(file=file.path(dir.out,'Final_Model_3.2_pred.sp_1000.RData'))
# load(file=file.path(dir.out,'Final_Model_3.3_pred.sp_200.RData'))
pred.cnt <- pred.sp$pred.cnt
pred.ba <- pred.sp$pred.ba
pred.z <- pred.sp$pred.z

pred.cnt.weibull <- pred.sp$pred.cnt
pred.ba.weibull <- pred.sp$pred.ba
pred.z.weibull <- pred.sp$pred.z

load(file=file.path(dir.out,'Final_Model_2.1_pred.sp.RData'))
pred.cnt.gpd <- pred.sp$pred.cnt
pred.ba.gpd <- pred.sp$pred.ba
pred.z.gpd <- pred.sp$pred.z

load(file=file.path(dir.out,'Final_Model_1.1_pred.sp.RData'))
pred.cnt.gamma <- pred.sp$pred.cnt
pred.ba.gamma <- pred.sp$pred.ba
pred.z.gamma <- pred.sp$pred.z



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


eta.z <- res3.2$summary.linear.predictor[idx.pred.z,1]

pred.z.prob <- exp(eta.z)/(1 + exp(eta.z))


A<- verify(data.fit2$z, pred.z.prob, frcst.type = "prob", obs.type = "binary")
reliability.plot(A, titl = "Alternative plot")

leq.u <- function(x, u){
  return(sum(x<=u))
}
eq.u <- function(x,u){
  return(sum(x==u))
}
bet.u <- function(x,u){
  return(sum(x>=u[1] & x<u[2]))
}

u.cnt <- c(0,1,2,3,4,5,9)
pred.cnt.weibull.prob <- rep(NA,length(u.cnt)-1)
pred.cnt.gamma.prob <- rep(NA,length(u.cnt)-1)
pred.cnt.gpd.prob <- rep(NA,length(u.cnt)-1)
obs.cnt.prob <- rep(NA,length(u.cnt)-1)
u.cnt <- c(0,1,2,3,4,5,9)
u.cnt.label <- c('[0,1)','[1,2)', '[2,3)', '[3,4)','[4,5)','[5,+infty)')
for (i in 1:(length(u.cnt)-1)){
  pred.cnt.weibull.prob[i] <- sum(sapply(pred.cnt.weibull, bet.u, u.cnt[i:(i+1)]))/200/length(pred.cnt.weibull)
  pred.cnt.gamma.prob[i] <- sum(sapply(pred.cnt.gamma, bet.u, u.cnt[i:(i+1)]))/200/length(pred.cnt.gamma)
  pred.cnt.gpd.prob[i] <- sum(sapply(pred.cnt.gpd, bet.u, u.cnt[i:(i+1)]))/200/length(pred.cnt.gpd)
  obs.cnt.prob[i] <- sum(bet.u(data.fit2$y, u.cnt[i:(i+1)]))/nrow(data.fit2)
}
plot(log(pred.cnt.prob),log(obs.cnt.prob))
abline(a=0,b=1)

df.realibity.cnt <- data.frame('pred.cnt.prob'=c(pred.cnt.weibull.prob, pred.cnt.gamma.prob,pred.cnt.gpd.prob),
                               'obs.cnt.prob'=rep(obs.cnt.prob,3),
                               'model' = rep(c('weibull','gamma','gpd'),each=length(obs.cnt.prob)),
                               'cnt'=rep(u.cnt.label,3))

cc.cnt <- scales::seq_gradient_pal("#FABC3F", "#E85C0D", "Lab")(seq(0,1,length.out=length(u.cnt.label)))

ggplot(df.realibity.cnt[df.realibity.cnt$model=='weibull',], aes(x=log(pred.cnt.prob),y=log(obs.cnt.prob), color=cnt)) + 
  geom_point(size=2.5, shape=17)+
  # geom_line(group=1,color='black')+
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black") +
  scale_color_manual(values = cc.cnt)  

ggplot(df.realibity.cnt, aes(x=cnt,y=(pred.cnt.prob)-(obs.cnt.prob),color=model,group = model)) + 
  geom_point(size=2.5, shape=17)+
  geom_line(alpha=0.5)


cnt.weibull.unlist <- unlist(pred.cnt.weibull)
cnt.gamma.unlist <- unlist(pred.cnt.gamma)
cnt.gpd.unlist <- unlist(pred.cnt.gpd)


df.sharp.cnt.weibull <- as.data.frame(table(cnt.weibull.unlist[cnt.weibull.unlist>0]))
df.sharp.cnt.weibull$model <- 'weibull'

df.sharp.cnt.gamma <- as.data.frame(table(cnt.gamma.unlist[cnt.gamma.unlist>0]))
df.sharp.cnt.gamma$model <- 'gamma'

df.sharp.cnt.gpd <- as.data.frame(table(cnt.gpd.unlist[cnt.gpd.unlist>0]))
df.sharp.cnt.gpd$model <- 'gpd'

df.sharp.cnt <- rbind(df.sharp.cnt.weibull,df.sharp.cnt.gamma,df.sharp.cnt.gpd)
colnames(df.sharp.cnt) <- c('cnt','freq','model')

ggplot(data=df.sharp.cnt, aes(x=cnt, y=freq, fill=model)) +
  geom_bar(stat="identity",position=position_dodge())

ggplot(data=df.sharp.cnt[df.sharp.cnt$freq<200,], aes(x=cnt, y=freq, fill=model)) +
  geom_bar(stat="identity",position=position_dodge())

u.ba <- c(1,seq(10,100,10),150,200,250,300,400,500,1000,1500,2000,5000,10000,20000,30000,40000,100000)
u.log.ba <- c(0,1,log(u.ba)[-1]) 
u.ba.label <- as.factor(c('[0,1)','[1,10)', '[10,20)', '[20,30)', '[30,40)','[40,50)',
                '[50,60)','[60,70)','[70,80)','[80,90)','[90,100)',
                '[100,150)','[150,200)','[200,250)','[250,300)',
                '[300,400)','[400,500)','[500,1000)','[1000,1500)',
                '[1500,2000)','[2000,5000)','[5000,10000)','[10000,20000)',
                '[20000,30000)','[30000,40000)','[40000,+infty)'))
levels(u.ba.label) <- c('[0,1)','[1,10)', '[10,20)', '[20,30)', '[30,40)','[40,50)',
                        '[50,60)','[60,70)','[70,80)','[80,90)','[90,100)',
                        '[100,150)','[150,200)','[200,250)','[250,300)',
                        '[300,400)','[400,500)','[500,1000)','[1000,1500)',
                        '[1500,2000)','[2000,5000)','[5000,10000)','[10000,20000)',
                        '[20000,30000)','[30000,40000)','[40000,+infty)')
pred.ba.weibull.prob <- rep(NA, length(u.log.ba)-1 )
pred.ba.gamma.prob <- rep(NA, length(u.log.ba)-1 )
pred.ba.gpd.prob <- rep(NA, length(u.log.ba)-1 )
obs.ba.prob <- rep(NA, length(u.log.ba)-1 )
for (i in 1:(length(u.log.ba)-1)){
  pred.ba.weibull.prob[i] <- sum(sapply(pred.ba.weibull, bet.u, u.log.ba[i:(i+1)]))/200/length(pred.ba.weibull)
  pred.ba.gamma.prob[i] <- sum(sapply(pred.ba.gamma, bet.u, u.log.ba[i:(i+1)]))/200/length(pred.ba.gamma)
  pred.ba.gpd.prob[i] <- sum(sapply(pred.ba.gpd, bet.u, u.log.ba[i:(i+1)]))/200/length(pred.ba.gpd)
  obs.ba.prob[i] <- sum(bet.u(data.fit2$log_ba, u.log.ba[i:(i+1)]))/nrow(data.fit2)
}
  
plot(log(pred.ba.prob),log(obs.ba.prob))
abline(a=0,b=1)


cc.ba <- scales::seq_gradient_pal("#FABC3F", "#E85C0D", "Lab")(seq(0,1,length.out=length(u.ba.label)))


df.realibity.ba <- data.frame('pred.ba.prob'=c(pred.ba.weibull.prob, pred.ba.gamma.prob,pred.ba.gpd.prob),
                              'model' = rep(c('weibull','gamma','gpd'),each=length(obs.ba.prob)),
                              'obs.ba.prob'=rep(obs.ba.prob,3),
                              'ba'=u.ba.label)

ggplot(df.realibity.ba[df.realibity.ba$model=='weibull',], aes(x=log(pred.ba.prob),y=log(obs.ba.prob), color=ba)) + 
  geom_point(size=2.5, shape=17)+
  # geom_line(group=1,color='black')+
  xlim(c(-11,0))+
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black") +
  scale_color_manual(values = cc.ba)  


ggplot(df.realibity.ba, aes(x=ba,y=(pred.ba.prob)-(obs.ba.prob))) + 
  geom_point(size=2.5, shape=17)+
  geom_line(group=1)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
  # geom_line(group=1,color='black')+
  # geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black") 

ggplot(df.realibity.ba, aes(x=ba,y=(pred.ba.prob)-(obs.ba.prob),color=model,group = model)) + 
  geom_point(size=2.5, shape=17)+
  geom_line(alpha=0.5)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))




ggplot(df.realibity.ba, aes(x=cnt,y=pred.ba.prob)) + 
  geom_bar(stat="identity")

ggplot(df.realibity.ba[df.realibity.ba$cnt!='[0,1)',], aes(x=cnt,y=pred.ba.prob)) + 
  geom_bar(stat="identity")




ba.weibull.unlist <- unlist(pred.ba.weibull)
ba.gamma.unlist <- unlist(pred.ba.gamma)
ba.gpd.unlist <- unlist(pred.ba.gpd)


df.sharp.ba.weibull <- data.frame('log.ba'=ba.weibull.unlist[ba.weibull.unlist>0])
df.sharp.ba.weibull$model <- 'weibull'

df.sharp.ba.gamma <- data.frame('log.ba'=ba.gamma.unlist[ba.gamma.unlist>0])
df.sharp.ba.gamma$model <- 'gamma'

df.sharp.ba.gpd <- data.frame('log.ba'=ba.gpd.unlist[ba.gpd.unlist>0])
df.sharp.ba.gpd$model <- 'gpd'

df.sharp.ba <- rbind(df.sharp.ba.weibull,df.sharp.ba.gamma,df.sharp.ba.gpd)

ggplot(df.sharp.ba, aes(x=model, y=log.ba, fill=model)) + 
  geom_boxplot() 


u.z <- seq(0,1,0.1)
pred.z.prob <- rep(NA,length(u.z)-1)
obs.z.prob <- rep(NA,length(u.z)-1)

for (i in 1:length(u.z)-1){
  pred.z.prob[i] <- sum(sapply(pred.z, bet.u, u.z[i:(i+1)]))/200/length(pred.z)
  obs.z.prob[i] <- sum(bet.u(data.fit2$z, u.z[i:(i+1)]))/nrow(data.fit2)
}
plot(pred.z.prob,(obs.z.prob))
abline(a=0,b=1)

plot(0:9, (obs.cnt.prob)-(pred.cnt.prob),type='l')
plot(0:9, log(obs.cnt.prob)/log(pred.cnt.prob),type='l')


  
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
sf_districts <- st_as_sf(dist)

grid.cell.coord <- st_as_sf(data.fit2, coords = c("lon.grid", "lat.grid"), crs = 4326)

# merged_sf <- st_join(grid.cell.coord, sf_districts[,'NAME_1'], join = st_within)
merged_sf1 <- st_join(grid.cell.coord, sf_districts[,'NAME_1'], join = st_nearest_feature)

merged_sf1 %>% group_by(NAME_1) %>% summarize(n_cnt = sum(y)) %>% arrange(desc(n_cnt), ascending=F)

merged_sf1 %>% group_by(NAME_1) %>% summarize(n_ba = sum(log_ba)) %>% arrange(desc(n_ba), ascending=F)


merged_sf1 %>% group_by(month.idx) %>% summarize(n_cnt = sum(y))


#----------------------boxplot of CNT/BA in single district --------------
ggplot() +geom_sf(data=merged_sf1, aes(fill=NAME_1,col=NAME_1)) + 
  geom_sf(data = sf_districts,color = "black", fill = NA)


# dist.name <- 'Castelo Branco'
# dist.name <- 'Guarda'
dist.name <- 'Vila Real'
joint.post.sp <- function(x) Reduce("+", x)
df.dist <- data.frame(month=rep(1:108, each=length(pred.cnt[[1]])))
df.dist$year_label <- 2012 + (df.dist$month-1)%/%12
df.dist$month_label <- sprintf("%02d", (df.dist$month-1)%%12+1)
df.dist$time_label <- paste(df.dist$year_label ,df.dist$month_label  ,sep='_')

df.boxplot.true <- data.frame(month=1:108, cnt=rep(NA, 108), log.ba= rep(NA,108))
df.boxplot.true$year_label <- 2012 + (df.boxplot.true$month-1)%/%12
df.boxplot.true$month_label <- sprintf("%02d", (df.boxplot.true$month-1)%%12+1)
df.boxplot.true$time_label <- paste(df.boxplot.true$year_label ,df.boxplot.true$month_label  ,sep='_')

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

ggplot(df.dist[df.dist$month>=85,], aes(x = factor(time_label), y = sample_cnt)) +
  geom_boxplot() + 
  geom_line(data=df.boxplot.true[df.boxplot.true$month>=85,], aes(x=factor(time_label),y=cnt,group = 1), col='red',linewidth=1)+
  labs(x = "Time", y = "fire count", title = paste("Predictive Total Fire Count in", dist.name,"During 2019-2020", sep=' ')) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        axis.title.x = element_blank())  # Rotate x-axis labels for clarity+


ggplot(df.dist[df.dist$month>=85,], aes(x = factor(time_label), y = sample_ba)) +
  geom_boxplot() + 
  geom_line(data=df.boxplot.true[df.boxplot.true$month>=85,], aes(x=factor(time_label),y=log.ba,group = 1), col='red',linewidth=1)+
  labs(x = "Time", y = "log burn area", title = paste("Predictive Total Burn Area in", dist.name,"During 2019-2020", sep=' ')) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        axis.title.x = element_blank())  # Rotate x-axis labels for clarity+




# 
# ggplot(df.dist[df.dist$month<=84  & df.dist$month>=61 ,], aes(x = factor(month), y = sample_ba)) +
#   geom_boxplot() + 
#   geom_line(data=df.boxplot.true[df.boxplot.true$month<=84 & df.boxplot.true$month>=61,], aes(x=factor(month),y=log.ba,group = 1), col='red',linewidth=1)+
#   labs(x = "Month", y = "Posterior Predictive Sample", title = "Posterior Predictive Samples over Time") +
#   theme_minimal() +
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))  # Rotate x-axis labels for clarity
# 
# 
# ggplot(df.dist[df.dist$month<=84  & df.dist$month>=61,], aes(x = factor(month), y = sample_cnt)) +
#   geom_boxplot() + 
#   geom_line(data=df.boxplot.true[df.boxplot.true$month<=84 & df.boxplot.true$month>=61,], aes(x=factor(month),y=cnt,group = 1), col='red',linewidth=1)+
#   labs(x = "Month", y = "Posterior Predictive Sample", title = "Posterior Predictive Samples over Time") +
#   theme_minimal() +
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))  # Rotate x-axis labels for clarity

#----------------------boxplot of CNT/BA in single month --------------


fire.in.month <- st_drop_geometry(merged_sf1) %>% group_by(NAME_1, time.idx) %>% 
             summarize(cnt.true = sum(y),
                       ba.true = sum(log_ba),
                       )
pred.in.month <- c()
for (i in 1:nrow(fire.in.month)){
  
  idx <- which(merged_sf1$NAME_1==as.character(fire.in.month[i,]$NAME_1) & merged_sf1$time.idx==fire.in.month[i,]$time.idx)
  pred.cnt.dist <- joint.post.sp(pred.cnt[idx])
  pred.ba.dist <- joint.post.sp(pred.ba[idx])
  
  pred.in.month <- rbind(pred.in.month,fire.in.month[i,] %>% 
    mutate(cnt.pred = quantile(pred.cnt.dist,0.5),
           ba.pred = quantile(pred.ba.dist,0.5),
           
           cnt.upp = quantile(pred.cnt.dist,0.975),
           ba.upp = quantile(pred.ba.dist,0.975),
           
           cnt.low = quantile(pred.cnt.dist,0.025),
           ba.low = quantile(pred.ba.dist,0.025)
  
           ) 
  )
}


pred.in.month.1 <- merge(sf_districts[,'NAME_1'],pred.in.month,  by='NAME_1')


st_drop_geometry(pred.in.month)

year <- 2019
month <- 8
idx <- 12*(year-2012) + month


lprange.scale.cnt <- c(0, max(pred.in.month[pred.in.month$time.idx>=85,c('cnt.true','cnt.pred')]))

csc.scale.cnt <- scale_fill_gradient( low = "#F7F7F7", high = "#E4003A", limits = lprange.scale.cnt)

ggplot() + 
  geom_sf(data=pred.in.month.1[pred.in.month.1$time.idx==idx,],aes(fill = cnt.pred),lwd = 0.1)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        axis.title.x = element_blank())+csc.scale.cnt+
  labs(title=paste('PM of Fire Count in', 'Aug', year,sep=' '))
  

ggplot() + geom_sf(data=pred.in.month.1[pred.in.month.1$time.idx==idx,],aes(fill = cnt.upp),lwd = 0.1)
ggplot() + geom_sf(data=pred.in.month.1[pred.in.month.1$time.idx==idx,],aes(fill = cnt.low),lwd = 0.1)


ggplot() + 
  geom_sf(data=pred.in.month.1[pred.in.month.1$time.idx==idx,],aes(fill = cnt.true),lwd = 0.1)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        axis.title.x = element_blank())+csc.scale.cnt+
  labs(title='Actual Fire Count in Aug 2019')



lprange.scale.ba <- c(0, max(pred.in.month[pred.in.month$time.idx>=85,c('ba.true','ba.pred')]))

csc.scale.ba <- scale_fill_gradient( low = "#F7F7F7", high = "#E4003A", limits = lprange.scale.ba)

ggplot() + 
  geom_sf(data=pred.in.month.1[pred.in.month.1$time.idx==idx,],aes(fill = ba.pred),lwd = 0.1)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        axis.title.x = element_blank())+csc.scale.ba+
  labs(title='PM of Burn Area in Aug 2019')

ggplot() +
  geom_sf(data=pred.in.month.1[pred.in.month.1$time.idx==idx,],aes(fill = ba.true),lwd = 0.1)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        axis.title.x = element_blank())+csc.scale.ba+
  labs(title='Actual Burn Area in Aug 2019')



#-------------------------------effect of covariates-------------------------

samples = inla.posterior.sample(200, result = res3.2, seed=1234)

n.grid <- 192
idx.pred.pois <- 1:n1
idx.pred.z <- (n1+1):(2*n1)
idx.pred.ba <- (2*n1+1):(3*n1)

idx.predicator <- (3*n1+1):(3*n1+679) # 679 columns in the A matrix (all.stack$A). 


# only for bym models
n2 <- 3*n1+679

idx.time.idx1 <- (n2+1):(n2+12)
idx.score_1 <- (idx.time.idx1[length(idx.time.idx1)]+1):(idx.time.idx1[length(idx.time.idx1)]+ spde_score_1$n.spde)

idx.idarea2.u <- (idx.score_1[length(idx.score_1)]+1):(idx.score_1[length(idx.score_1)] + n.grid)
idx.idarea2.v <- (idx.idarea2.u[length(idx.idarea2.u)]+1):(idx.idarea2.u[length(idx.idarea2.u)] + n.grid)

idx.time.idx2 <- (idx.idarea2.v[length(idx.idarea2.v)]+1):(idx.idarea2.v[length(idx.idarea2.v)]+12)
idx.score_2 <- (idx.time.idx2[length(idx.time.idx2)]+1):(idx.time.idx2[length(idx.time.idx2)]+ spde_score_2$n.spde)

idx.time.idx3 <- (idx.score_2[length(idx.score_2)]+1):(idx.score_2[length(idx.score_2)]+12)
idx.score_3 <- (idx.time.idx3[length(idx.time.idx3)]+1):(idx.time.idx3[length(idx.time.idx3)]+ spde_score_3$n.spde)

idx.idarea1.u <- (idx.score_3[length(idx.score_3)]+1):(idx.score_3[length(idx.score_3)]+n.grid)
idx.idarea1.v <- (idx.idarea1.u[length(idx.idarea1.u)]+1):(idx.idarea1.u[length(idx.idarea1.u)]+n.grid)

idx.idarea3.u <- (idx.idarea1.v[length(idx.idarea1.v)]+1):(idx.idarea1.v[length(idx.idarea1.v)]+n.grid)
idx.idarea3.v <- (idx.idarea3.u[length(idx.idarea3.u)]+1):(idx.idarea3.u[length(idx.idarea3.u)]+n.grid)

idx.Intercept1 <- idx.idarea3.v[length(idx.idarea3.v)]+1
idx.Intercept2 <- idx.Intercept1 + 1
idx.Intercept3 <- idx.Intercept1 + 2

samples[[1]]$latent[idx.pred.pois,][193:200]

# (samples[[1]]$latent[idx.Intercept1,] + samples[[1]]$latent[idx.idarea1.u,] + 
# samples[[1]]$latent[idx.time.idx1[2],] + as.vector(A1%*%matrix(samples[[1]]$latent[idx.score_1,])))[193:200]


pred.manual <- samples[[1]]$latent[idx.Intercept1,] + rowSums(expand.grid(samples[[1]]$latent[idx.idarea1.u,],samples[[1]]$latent[idx.time.idx1,]))    +
 +  as.vector(A1%*%matrix(samples[[1]]$latent[idx.score_1,]))




df.score.cnt <- data.frame('score' = data.fit2$score_cnt)
mat.score.cnt <- c()
mat.score.z <- c()
mat.score.ba <- c()

mat.time.cnt <- c()
mat.time.z <- c()
mat.time.ba <- c()

mat.spat.cnt <- c()
mat.spat.z <- c()
mat.spat.ba <- c()

for (i in 1:length(samples)){
  lp.score.cnt <- as.vector(A1%*%matrix(samples[[i]]$latent[idx.score_1,]))
  lp.score.z <- as.vector(A2%*%matrix(samples[[i]]$latent[idx.score_2,]))
  lp.score.ba <- as.vector(A3%*%matrix(samples[[i]]$latent[idx.score_3,]))
  
  mat.score.cnt <- rbind(mat.score.cnt, lp.score.cnt)
  mat.score.z <- rbind(mat.score.z, lp.score.z)
  mat.score.ba <- rbind(mat.score.ba, lp.score.ba)
  
  mat.time.cnt <- rbind(mat.time.cnt,samples[[i]]$latent[idx.time.idx1,])
  mat.time.z <- rbind(mat.time.z,samples[[i]]$latent[idx.time.idx2,])
  mat.time.ba <- rbind(mat.time.ba,samples[[i]]$latent[idx.time.idx3,])
  
  mat.spat.cnt <- rbind(mat.spat.cnt,samples[[i]]$latent[idx.Intercept1,] + 
                        rowSums(expand.grid(samples[[i]]$latent[idx.idarea1.u,],samples[[i]]$latent[idx.time.idx1,])))
  mat.spat.z <- rbind(mat.spat.z,samples[[i]]$latent[idx.Intercept2,] + 
                          rowSums(expand.grid(samples[[i]]$latent[idx.idarea2.u,],samples[[i]]$latent[idx.time.idx2,])))
  mat.spat.ba <- rbind(mat.spat.ba,samples[[i]]$latent[idx.Intercept3,] + 
                          rowSums(expand.grid(samples[[i]]$latent[idx.idarea3.u,],samples[[i]]$latent[idx.time.idx3,])))
}

data.fit2$effect_score_cnt_mean <- colMeans(mat.score.cnt)
data.fit2$effect_score_cnt_upp <- apply(mat.score.cnt,2 ,quantile,0.975)
data.fit2$effect_score_cnt_low <- apply(mat.score.cnt,2 ,quantile,0.025)

data.fit2$effect_score_z_mean <- colMeans(mat.score.z)
data.fit2$effect_score_z_upp <- apply(mat.score.z,2 ,quantile,0.975)
data.fit2$effect_score_z_low <- apply(mat.score.z,2 ,quantile,0.025)

data.fit2$effect_score_ba_mean <- colMeans(mat.score.ba)
data.fit2$effect_score_ba_upp <- apply(mat.score.ba,2 ,quantile,0.975)
data.fit2$effect_score_ba_low <- apply(mat.score.ba,2 ,quantile,0.025)

data.fit2$effect_spat_cnt_mean <- colMeans(mat.spat.cnt)
data.fit2$effect_spat_cnt_upp <- apply(mat.spat.cnt,2 ,quantile,0.975)
data.fit2$effect_spat_cnt_low <- apply(mat.spat.cnt,2 ,quantile,0.025)

data.fit2$effect_spat_z_mean <- colMeans(mat.spat.z)
data.fit2$effect_spat_z_upp <- apply(mat.spat.z,2 ,quantile,0.975)
data.fit2$effect_spat_z_low <- apply(mat.spat.z,2 ,quantile,0.025)

data.fit2$effect_spat_ba_mean <- colMeans(mat.spat.ba)
data.fit2$effect_spat_ba_upp <- apply(mat.spat.ba,2 ,quantile,0.975)
data.fit2$effect_spat_ba_low <- apply(mat.spat.ba,2 ,quantile,0.025)

ggplot(data.fit2) +
  stat_smooth(aes(x = score_cnt, y = effect_score_cnt_mean), col='#EB5B00') + 
  stat_smooth(aes(x = score_cnt, y = effect_score_cnt_upp),linetype='dashed',col='#EB5B00') + 
  stat_smooth(aes(x = score_cnt, y = effect_score_cnt_low),linetype='dashed',col='#EB5B00')+
  geom_point(aes(x = score_cnt, y = -5),size=0.8) + 
  labs(y='linear predicator')
 

ggplot(data.fit2) +
  stat_smooth(aes(x = score_z, y = effect_score_z_mean), col='#EB5B00') + 
  stat_smooth(aes(x = score_z, y = effect_score_z_upp),linetype='dashed',col='#EB5B00') + 
  stat_smooth(aes(x = score_z, y = effect_score_z_low),linetype='dashed',col='#EB5B00') +
  geom_point(aes(x = score_z, y = -6),size=0.8)+
  labs(y='linear predicator')


ggplot(data.fit2) +
  stat_smooth(aes(x = score_ba, y = effect_score_ba_mean),col='#EB5B00') + 
  stat_smooth(aes(x = score_ba, y = effect_score_ba_upp),linetype='dashed',col='#EB5B00') + 
  stat_smooth(aes(x = score_ba, y = effect_score_ba_low),linetype='dashed',col='#EB5B00') + 
  geom_point(aes(x = score_ba, y = -1.5),size=0.8)+
  labs(y='linear predicator')

  



B2_sf <- st_as_sf(B2.merge)
library(tidyr)
B2_sf <- gather(B2_sf, y.time.idx, y, paste0("y.", 1:108))
B2_sf$time.idx <- as.integer(substring(B2_sf$y.time.idx, 3, 5))

B2_sf <- merge(
  B2_sf[,c('grid.idx','time.idx')], data.fit2,
  
  by.x=c('grid.idx','time.idx'),
  by.y=c('grid.idx','time.idx'),
)

year <- 2019
month <- 1
idx <- 12*(year-2012) + month


lprange.scale.effect <- c(0, 0.9)

csc.scale.effect <- scale_fill_gradient( low = "#F7F7F7", high = "#E4003A", limits = lprange.scale.effect)


ggplot() + geom_sf(data=B2_sf[B2_sf$time.idx==idx,],aes(fill = exp(effect_spat_z_mean)/(1+exp(effect_spat_z_mean))),lwd = 0.1) +
  # geom_sf(data=st_as_sf(district),alpha = 0)+
  # facet_wrap(~time.idx, dir = "h", ncol =3) +
  # ggtitle(paste("Actual y in 2012")) +
  theme(legend.position = "right",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)
        # axis.text.y = element_blank(),
        # axis.ticks = element_blank()
  )+ csc.scale.effect+
  labs(fill  = "BYM2 + Month Effect",title='Mean of the linear predicator of BYM2 + Month' )

ggplot() + geom_sf(data=B2_sf[B2_sf$time.idx==idx,],aes(fill = score_z),lwd = 0.1) +
  # geom_sf(data=st_as_sf(district),alpha = 0)+
  # facet_wrap(~time.idx, dir = "h", ncol =3) +
  # ggtitle(paste("Actual y in 2012")) +
  theme(legend.position = "right",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)
        # axis.text.y = element_blank(),
        # axis.ticks = element_blank()
  )+csc.scale.effect+
  labs(title='Score z')


ggplot() + geom_sf(data=B2_sf[B2_sf$time.idx==103,],aes(fill = exp(effect_score_cnt_mean)),lwd = 0.1) +
  # geom_sf(data=st_as_sf(district),alpha = 0)+
  # facet_wrap(~time.idx, dir = "h", ncol =3) +
  # ggtitle(paste("Actual y in 2012")) +
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)
        # axis.text.y = element_blank(),
        # axis.ticks = element_blank()
  )+csc.scale.effect

ggplot() + geom_sf(data=B2_sf[B2_sf$time.idx==103,],aes(fill = score_cnt),lwd = 0.1) +
  # geom_sf(data=st_as_sf(district),alpha = 0)+
  # facet_wrap(~time.idx, dir = "h", ncol =3) +
  # ggtitle(paste("Actual y in 2012")) +
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)
        # axis.text.y = element_blank(),
        # axis.ticks = element_blank()
  )+csc.scale.effect
#------------------------- feature importance for xgboost model-------------
library(SHAPforxgboost)

shap_values <- shap.values(xgb_model = model.z, X_train = TrainSet1)
# The ranked features by mean |SHAP|
shap_values$mean_shap_score


# To prepare the long-format data:
shap_long <- shap.prep(xgb_model = model.z, X_train = as.matrix(data.rf[, c(covar.names)]))
# is the same as: using given shap_contrib
shap_long <- shap.prep(shap_contrib = shap_values$shap_score, X_train = as.matrix(data.rf[, c(covar.names)]))
# (Notice that there will be a data.table warning from `melt.data.table` due to `dayint` coerced from
# integer to double)

# **SHAP summary plot**
shap.plot.summary(shap_long)

# sometimes for a preview, you want to plot less data to make it faster using `dilute`
shap.plot.summary(shap_long, x_bound  = 1.2, dilute = 10)

# Alternatives options to make the same plot:
# option 1: start with the xgboost model
shap.plot.summary.wrap1(mod, X = dataX)



fig_list = lapply(names(shap_values$mean_shap_score)[1:6], shap.plot.dependence, 
                  data_long = shap_long, dilute = 5)
gridExtra::grid.arrange(grobs = fig_list, ncol = 2)


# choose to show top 4 features by setting `top_n = 4`, set 6 clustering groups.  
plot_data <- shap.prep.stack.data(shap_contrib = shap_values$shap_score, top_n = 4, n_groups = 6)

# choose to zoom in at location 500, set y-axis limit using `y_parent_limit`  
# it is also possible to set y-axis limit for zoom-in part alone using `y_zoomin_limit`  
shap.plot.force_plot(plot_data, zoom_in_location = 500, y_parent_limit = c(-1,1))

# plot by each cluster
shap.plot.force_plot_bygroup(plot_data)

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
