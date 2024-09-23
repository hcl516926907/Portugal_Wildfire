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
# library(rgeos)
# library(spdep)
library(raster)
# library(scoringRules)
# bru_safe_sp(force = TRUE)

library(doParallel)
cl <- detectCores()
print(cl)
cl <- 10

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
                 # 'lon.grid','lat.grid',
                 'month.idx')

data.rf <- data.fit2[data.fit2$time.idx <= 84, c(covar.names,'lon.grid','lat.grid','area_ha','z','y')]
# for (var in c('LVegTyp','HVegTyp','month.idx')){
#   data.rf[,var] <- as.factor(data.rf[,var])
# }

library(xgboost)
library(rBayesianOptimization)
library(pROC)
data.rf$log_ba <- log(data.rf$area_ha)
data.rf[is.na(data.rf$log_ba),'log_ba'] <- 0

set.seed(100)

xgb_z_cv   <- function(eta, max_depth, min_child_weight, gamma, subsample, colsample_bytree,
                       scale_pos_weight) {
  k <- 7
  auc_train_list <- numeric(k)
  auc_test_list <- numeric(k)
  data <- data.rf[, c(covar.names)]
  target <- data.rf[, 'z']
  
  params <- list(
    objective = "binary:logistic",
    eval_metric = "auc",
    eta = eta,
    max_depth = as.integer(max_depth),
    min_child_weight = min_child_weight,
    gamma = gamma,
    subsample = subsample,
    colsample_bytree = colsample_bytree,
    scale_pos_weight = scale_pos_weight
  )
  
  for (i in 1:k) {
    index <- (1:(192*12))+(i-1)*(192*12)
    train_data <- data[-index, ]
    train_target <- target[-index]
    test_data <- data[index, ]
    test_target <- target[index]
    dtrain <- xgb.DMatrix(data = as.matrix(train_data), label = train_target)
    dtest <- xgb.DMatrix(data = as.matrix(test_data), label = test_target)
    
    set.seed(1234)
    model <- xgb.train(params = params, data = dtrain, nrounds = 100,
                       nthread = cl)
    
    train_pred <- predict(model, dtrain)
    test_pred <- predict(model, dtest)
    
    auc_train <- auc(train_target, train_pred, quiet = TRUE)
    auc_test <- auc(test_target, test_pred, quiet = TRUE)
    
    auc_train_list[i] <- auc_train
    auc_test_list[i] <- auc_test
  }
  
  mean_auc_train <- mean(auc_train_list)
  mean_auc_test <- mean(auc_test_list)
  
  return(list(Score=mean_auc_test,Pred=NULL))
}

bounds_z <- list(
  eta = c(0.01, 0.1),
  max_depth = c(2L, 3L, 4L, 5L, 6L,7L,8L),
  min_child_weight = c(1,10),
  gamma = c(0,1),
  subsample = c(0.4, 1.0),
  colsample_bytree = c(0.4, 1.0),
  scale_pos_weight = c(20,40)
)

set.seed(1234)
t1 <- Sys.time()
xgb_z_params <- BayesianOptimization(
  FUN = xgb_z_cv,
  bounds = bounds_z,
  init_points = 10, # Number of initial random points
  n_iter = 50,      # Number of iterations for optimization
  acq = "ucb",      # Acquisition function, e.g., "ucb", "ei", "poi"
  kappa = 2.576,    # Parameter for the acquisition function
  eps = 0.0,        # Exploration-exploitation trade-off
  verbose = TRUE,
  kernel = list(type = "matern", nu = 5/2)
)
t2 <- Sys.time()
print(t2-t1)

z_best_params <- xgb_z_params$Best_Par
print(z_best_params)



xgb_z_lonlat_cv   <- function(eta, max_depth, min_child_weight, gamma, subsample, colsample_bytree,
                       scale_pos_weight,nrounds) {
  k <- 7
  auc_train_list <- numeric(k)
  auc_test_list <- numeric(k)
  data <- data.rf[, c(covar.names,'lon.grid','lat.grid')]
  target <- data.rf[, 'z']
  
  params <- list(
    objective = "binary:logistic",
    eval_metric = "auc",
    eta = eta,
    max_depth = as.integer(max_depth),
    min_child_weight = min_child_weight,
    gamma = gamma,
    subsample = subsample,
    colsample_bytree = colsample_bytree,
    scale_pos_weight = scale_pos_weight,
    nrounds = nrounds
  )
  
  for (i in 1:k) {
    index <- (1:(192*12))+(i-1)*(192*12)
    train_data <- data[-index, ]
    train_target <- target[-index]
    test_data <- data[index, ]
    test_target <- target[index]
    dtrain <- xgb.DMatrix(data = as.matrix(train_data), label = train_target)
    dtest <- xgb.DMatrix(data = as.matrix(test_data), label = test_target)
    
    set.seed(1234)
    model <- xgb.train(params = params, data = dtrain, nrounds = 100,
                       nthread = cl)
    
    train_pred <- predict(model, dtrain)
    test_pred <- predict(model, dtest)
    
    auc_train <- auc(train_target, train_pred, quiet = TRUE)
    auc_test <- auc(test_target, test_pred, quiet = TRUE)
    
    auc_train_list[i] <- auc_train
    auc_test_list[i] <- auc_test
  }
  
  mean_auc_train <- mean(auc_train_list)
  mean_auc_test <- mean(auc_test_list)
  
  return(list(Score=mean_auc_test,Pred=NULL))
}
set.seed(1234)
t1 <- Sys.time()
xgb_z_params_1 <- BayesianOptimization(
  FUN = xgb_z_lonlat_cv,
  bounds = bounds_z,
  init_points = 10, # Number of initial random points
  n_iter = 50,      # Number of iterations for optimization
  acq = "ucb",      # Acquisition function, e.g., "ucb", "ei", "poi"
  kappa = 2.576,    # Parameter for the acquisition function
  eps = 0.0,        # Exploration-exploitation trade-off
  verbose = TRUE,
  kernel = list(type = "matern", nu = 5/2)
)
t2 <- Sys.time()
print(t2-t1)

z_best_params_1 <- xgb_z_params_1$Best_Par
print(z_best_params_1)

