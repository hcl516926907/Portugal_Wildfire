dir.data <- "/home/pgrad2/2448355h/My_PhD_Project/00_Dataset/Porgual_Wildfire"
dir.out <- '/home/pgrad2/2448355h/My_PhD_Project/01_Output/Portugal_Wildfire'


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
load(file.path(dir.data, "Data_For_Fitting.RData"))



#####################################################################
# XGBoost
#####################################################################


data.fit$z <- as.vector((data.fit$y>0)+0)
data.fit$HVegTyp <- as.factor(data.fit$HVegTyp)
levels(data.fit$HVegTyp) <-  c('NoCover','Evergreen_broadleaf_trees','Mixed_forest','Interrupted_forest')
data.fit$LVegTyp <- as.factor(data.fit$LVegTyp)
levels(data.fit$LVegTyp) <- c('NoCover','Crops','Tall_grass','Semidesert','Evergreen_shrubs')
data.fit$month <- as.factor(data.fit$month)
data.fit <- dummy_cols(data.fit, 
                       select_columns = c('HVegTyp','LVegTyp','month'),remove_first_dummy=TRUE)

covar.names <- colnames(data.fit)[c(-1,-2,-5,-6,-7,-8,-13,-16,-24)]

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
  eta = 0.07,
  max_depth = 5,
  gamma = 1,
  scale_pos_weight = 18,
  colsample_bytree = 0.6,
  min_child_weight = 17,
  subsample = 0.8,
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
  nrounds = 85,
  eta = 0.08,
  max_depth = 2,
  gamma = 3,
  colsample_bytree = 0.7,
  min_child_weight = 6,
  subsample = 0.8,
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
  nrounds = 90,
  eta = 0.08,
  max_depth = 3,
  gamma = 0.2,
  colsample_bytree = 0.9,
  min_child_weight = 5,
  subsample = 0.6,
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
 
save(model.z,model.ba,model.cnt,file=file.path(dir.out,'XGBoost_Models.RData'))
####################################################################
# INLA
####################################################################


loc.data.utm <- st_transform(st_as_sf(data.fit, coords=c('lon.grid','lat.grid'), crs=4326 ),
                             crs='+proj=utm +zone=29 +datum=WGS84 +units=km +no_defs')

# cond <- data.new$month.idx==9 & data.new$length > 24*60
# coords <- SpatialPointsDataFrame(data.new[cond,],coords=data.new[cond,c('x_utm_new','y_utm_new')], 
#                                  proj4string=CRS('+proj=utm +zone=29 +datum=WGS84 +units=km +no_defs' ))
# coords$time.idx <- coords$year.idx
library(rnaturalearth)
map <- ne_countries(type = "countries", country = "Portugal",
                    scale = "medium", returnclass = "sf")
projutm <- "+proj=utm +zone=29 +datum=WGS84 +units=km +no_defs"
map <- st_transform(map, crs = projutm)
mainland_bbox <- st_as_sfc(st_bbox(c(xmin = 0, xmax = 1000, ymin = 4000 , ymax = 4800), crs = st_crs(map)))
map_mainland <- st_intersection(map, mainland_bbox)
mainland <- st_transform(map_mainland,crs=4326)
# ggplot() + geom_sf(data = map_mainland) +
#   geom_sf(data = st_as_sf(loc.data.utm)) + coord_sf(datum = projutm)

loc.d <- cbind(st_coordinates(map_mainland)[, 1], st_coordinates(map_mainland)[, 2])



domain.polys <- Polygons(list(Polygon(loc.d)), '0')
domainSP <- SpatialPolygons(list(domain.polys),proj4string=CRS(projutm))


B <- SpatialPolygonsDataFrame(domainSP, data.frame('weight'=1), match.ID = F) 

grid_era5 <- expand.grid(lon=seq(-10,-6,0.1), lat=seq(43,36,-0.1))

grd <- st_intersection(st_as_sf(grid_era5,coords=c('lon','lat'),crs=4326),mainland) 

grd <- SpatialPixelsDataFrame(points = st_coordinates(grd), data = data.frame(grid.idx = 1:nrow(grd)),proj4string = CRS("+proj=longlat +datum=WGS84"))

B1 <- as(grd, 'SpatialPolygonsDataFrame')
B1$lon.grid <- grd@coords[,1]
B1$lat.grid <- grd@coords[,2]

data.fit.reshape <- reshape(data.fit[,c('grid.idx','time.idx','y')],
                     timevar = "time.idx",
                     idvar = "grid.idx",
                     direction = "wide"
)

B2.merge <- merge(B1, data.fit.reshape, by = "grid.idx")

# save(B2.merge, file=file.path(dir.out, 'grid_cell_map_0.1.RData'))
# 


z <- as.vector((data.fit$y>0)+0)
log.ba <- as.vector(ifelse(data.fit$y>0, data.fit$log_ba, NA))
# cnt = as.vector(data.fit$y)
cnt <- as.vector(ifelse(data.fit$y>0, data.fit$y, NA))
#

#prepare for prediction
z[which(data.fit$year>2022)] <- NA
log.ba[which(data.fit$year>2022)] <- NA
cnt[which(data.fit$year>2022)] <- NA





coo <- as.matrix(st_coordinates(loc.data.utm))

# mesh.spat <- inla.mesh.2d(
#   # loc = coo,
#   boundary = domainSP,
#   offset = c(20, 30),
#   cutoff = 1, max.edge = c(40, 60)
# )

mesh.spat <- inla.mesh.2d(
  # loc = unique(coo),
  boundary = domainSP,
  offset = c(20, 30),
  cutoff = 1, max.edge = c(20, 30)
)

plot(mesh.spat)
points(unique(coo), col = "red",pch=19,cex=0.1)

mesh.spat$n


data.fit$month <- as.integer(data.fit$month)


spde.spat <- inla.spde2.pcmatern(mesh = mesh.spat,
                                 # PC-prior on range: P(practic.range < 0.05) = 0.01
                                 prior.range = c(60, 0.05),
                                 # PC-prior on sigma: P(sigma > 1) = 0.01
                                 prior.sigma = c(1, 0.01))

# spat.z <- inla.spde.make.index('spat.z', n.spde = spde.spat$n.spde,
#                                n.group = 12)
spat.z <- inla.spde.make.index('spat.z', n.spde = spde.spat$n.spde)

# A.spat <- inla.spde.make.A(mesh = mesh.spat, loc = coo,  group=data.fit$month)
A.spat <- inla.spde.make.A(mesh = mesh.spat, loc = coo)



# spat.cnt <-  inla.spde.make.index('spat.cnt', n.spde = spde.spat$n.spde,
#                                   n.group = 12)
# spat.ba <-  inla.spde.make.index('spat.ba', n.spde = spde.spat$n.spde,
#                                  n.group = 12)
spat.cnt <-  inla.spde.make.index('spat.cnt', n.spde = spde.spat$n.spde)
spat.ba <-  inla.spde.make.index('spat.ba', n.spde = spde.spat$n.spde)


knots_cnt <- quantile(data.fit[data.fit$y>0,'score_cnt'],c(0.05,0.2,0.35,0.5,0.65,0.8,0.95))
mesh_score_1 <- inla.mesh.1d(knots_cnt, boundary=c('free','free'))
A1 <- inla.spde.make.A(mesh_score_1, loc=data.fit$score_cnt)
spde_score_1 <-  inla.spde2.pcmatern(mesh_score_1,
                                     prior.range = c(0.5, 0.3),
                                     prior.sigma = c(1, 0.05))


spde_score_1.idx <- inla.spde.make.index("score_1", n.spde = spde_score_1$n.spde)

knots_z <- quantile(data.fit[,'score_z'],c(0.05,seq(0.1,0.9,0.1),0.95))
mesh_score_2 <- inla.mesh.1d(knots_z,boundary=c('free','free'))

A2 <- inla.spde.make.A(mesh_score_2, loc=data.fit$score_z)
spde_score_2 <-  inla.spde2.pcmatern(mesh_score_2,
                                     prior.range = c(0.1, 0.05),
                                     prior.sigma = c(1, 0.05))
spde_score_2.idx <- inla.spde.make.index("score_2", n.spde = spde_score_2$n.spde)

knots_ba <- quantile(data.fit[data.fit$y>0,'score_ba'],c(0.05,0.2,0.35,0.5,0.65,0.8,0.95))
mesh_score_3 <- inla.mesh.1d(knots_ba, boundary=c('free','free'))

A3 <- inla.spde.make.A(mesh_score_3, loc=data.fit$score_ba)
spde_score_3 <-  inla.spde2.pcmatern(mesh_score_3,
                                     prior.range = c(0.3, 0.2),
                                     prior.sigma = c(1, 0.1))
spde_score_3.idx <- inla.spde.make.index("score_3", n.spde = spde_score_3$n.spde)




cnt.stack <- inla.stack(
  data= list(Y.log=cbind(cnt,NA,NA)),
  A <- list(1,1,A1, A.spat ),
  effect = list(Intercept1=rep(1,nrow(data.fit)), month.idx.cnt=data.fit$month, score_1=spde_score_1.idx,
                spat.cnt=spat.cnt ),
  tag='cnt'
)

z.stack <- inla.stack(
  data= list(Y.log=cbind(NA,z,NA)),
  A <- list(1,1,A2, A.spat),
  effect = list(Intercept2=rep(1,nrow(data.fit)), month.idx.z=data.fit$month, score_2=spde_score_2.idx,
                spat.z=spat.z ),
  tag='z'
)


ba.stack <- inla.stack(
  data= list(Y.log=cbind(NA,NA,log.ba)),
  A <- list(1,1, A3, A.spat),
  effect = list(Intercept3=rep(1,nrow(data.fit)), month.idx.ba=data.fit$month, score_3=spde_score_3.idx,
                spat.ba = spat.ba),
  tag='ba'
)


all.stack <- inla.stack(cnt.stack, z.stack, ba.stack )



# ########-----------------------Prepare prediction data------------------------
# grd <- SpatialGrid(GridTopology(cellcentre.offset = c(-10, 36),
#                                 cellsize = c(0.125,0.125),
#                                 cells.dim = c(34, 58)))
# 
# 
# grd_sp <- SpatialPixelsDataFrame(points = grd, data = data.frame(grid.idx = 1:length(grd)),proj4string = CRS("+proj=longlat +datum=WGS84"))
# 
# grd_poly <- as(grd_sp, 'SpatialPolygonsDataFrame')
# grd_poly <- spTransform(grd_poly, CRS(projutm))
# grd_poly$lon.grid <- grd_sp@coords[,1]
# grd_poly$lat.grid <- grd_sp@coords[,2]
# 
# B2.pred <- as_Spatial(st_intersection(st_as_sf(grd_poly),st_as_sf(B)))
# coord.cent <- st_coordinates(st_centroid(st_as_sf(B2.pred)))
# B2.pred$grid.cent.x.utm <- coord.cent[,1]
# B2.pred$grid.cent.y.utm <- coord.cent[,2]
# B2.pred$E <- area(B2.pred)
# 
# B2.pred$grid.idx <- 1:nrow(B2.pred)
# 
# 
# ggplot(st_as_sf(B2.pred)) + geom_sf(aes(fill=grid.idx))
# 
# 
# 
# 
# data.merge.pred <- B2.pred |>  
#   st_as_sf() |> # cast to sf
#   mutate(grid_id = row_number()) |> # create unique ID
#   st_join(loc.data.utm) |> # join the species dataset
#   group_by(grid_id)
# 
# 
# data.fit.ba.pred <- data.merge.pred %>% st_drop_geometry() %>% filter( time.idx<=108, length >=24*60 )%>% 
#   group_by(grid.idx,grid_id, year.idx, month.idx, time.idx) %>%
#   summarise(area_ha = sum(area_ha), 
#             log_ba = log(area_ha),
#             y = n(),
#             lon.grid = mean(lon.grid),
#             lat.grid = mean(lat.grid),
#             x.utm = mean(grid.cent.x.utm),
#             y.utm = mean(grid.cent.y.utm),
#             E = mean(E))
# 
# 
# 
# data.fit.pred <- do.call(rbind, lapply(97:108, function(x) {B2.pred@data$time.idx = x
# return(B2.pred@data)}))
# print(dim(data.fit.pred))
# data.fit.pred <- merge(data.fit.pred, data.fit.ba.pred[,c('grid.idx','time.idx','year.idx', 'month.idx','y','area_ha','log_ba')],
#                         by=c('grid.idx','time.idx'),all.x=T)
# 
# 
# print(dim(data.fit.pred))
# data.fit.pred[is.na(data.fit.pred$y),'y'] <- 0
# summary(data.fit.pred$y)
# 
# data.fit.pred$month.idx <- (data.fit.pred$time.idx-1)%%12 + 1
# data.fit.pred$year.idx <- (data.fit.pred$time.idx-1)%/%12 + 1
# data.fit.pred$E1 <- 1
# 
# data.fit.pred.sf <- st_as_sf(data.fit.pred, coords=c('grid.cent.x.utm','grid.cent.y.utm'),crs='+proj=utm +zone=29 +datum=WGS84 +units=km +no_defs')
# data.fit.sf <- st_as_sf(data.fit, coords=c('grid.cent.x.utm','grid.cent.y.utm'),crs='+proj=utm +zone=29 +datum=WGS84 +units=km +no_defs')
# 
# 
# for (time.idx in 97:108){
#   data.fit.pred[data.fit.pred$time.idx==time.idx,
#                  c('score_z','score_ba','score_cnt',
#                    'score_z_bin','score_ba_bin','score_cnt_bin')] <- st_join(data.fit.pred.sf[data.fit.pred.sf$time.idx==time.idx,], 
#                                                                  data.fit.sf[data.fit.sf$time.idx==time.idx,c('score_z','score_ba','score_cnt',
#                                                                                                                 'score_z_bin','score_ba_bin','score_cnt_bin')], 
#                                                                  join = st_nearest_feature) |> st_drop_geometry()|> dplyr::select(score_z,score_ba,score_cnt,
#                                                                                                                                   score_z_bin,score_ba_bin,score_cnt_bin)
# }
# 
# data.fit.pred <- data.fit.pred[order(data.fit.pred$time.idx),]
# 
# data.fit.pred[is.na(data.fit.pred$log_ba),'log_ba'] <- 0
# 
# # ggplot(data.fit.sf[data.fit.sf$time.idx==104,]) + geom_sf(aes(fill=score_z,color=score_z))
# # ggplot(data.fit.pred.sf[data.fit.pred.sf$time.idx==104,]) + geom_sf(aes(fill=score_z,color=score_z))
# # 
# 
# coo.pred <- as.matrix(data.fit.pred[,c('grid.cent.x.utm','grid.cent.y.utm')])
# 
# 
# A.spat.1.pred <- inla.spde.make.A(mesh = mesh.spat, loc = coo.pred)
# 
# A.spat.2.pred <- inla.spde.make.A(mesh = mesh.spat, loc = coo.pred)
# 
# A.spat.3.pred <- inla.spde.make.A(mesh = mesh.spat, loc = coo.pred)
# 
# 
# A1.pred <- inla.spde.make.A(mesh_score_1, loc=data.fit.pred$score_cnt)
# 
# A2.pred <- inla.spde.make.A(mesh_score_2, loc=data.fit.pred$score_z)
# 
# A3.pred <- inla.spde.make.A(mesh_score_3, loc=data.fit.pred$score_ba)
# 
# 
# data.fit3.pred <- reshape(data.fit.pred[,c('grid.idx','time.idx','y')],
#                           timevar = "time.idx",
#                           idvar = "grid.idx",
#                           direction = "wide"
# )
# 
# B2.merge.pred <- merge(B2.pred, data.fit3.pred, by = "grid.idx")
# 
# # save(B2.merge.pred, file=file.path(dir.out, 'grid_cell_map_0.0625.RData'))
# # save(data.fit.pred, A.spat.1.pred, A.spat.2.pred, 
# #      A.spat.3.pred, A1.pred, A2.pred, A3.pred,
# #      file=file.path(dir.out, 'data.fit.pred_0.0625.RData'))
# 
# # save(B2.merge.pred, file=file.path(dir.out, 'grid_cell_map_0.125.RData'))
# # save(data.fit.pred, A.spat.1.pred, A.spat.2.pred,
# #      A.spat.3.pred, A1.pred, A2.pred, A3.pred,
# #      file=file.path(dir.out, 'data.fit.pred_0.125.RData'))
# #--------------------------------------------------------------


n1 <- dim(data.fit)[1]


hyper.rw <-  list(prec = list(prior="loggamma",param=c(1,1)))
hyper.time.rw <- list(prec = list(prior="loggamma",param=c(80,20)))

t1 <- Sys.time()

#
formula <- Y.log ~  -1  +
  Intercept1 + f(spat.cnt, copy='spat.z',fixed=F)  + f(score_1, model=spde_score_1 )+
  Intercept2 + f(spat.z, model=spde.spat)  + f(score_2, model=spde_score_2 )+
  Intercept3 + f(spat.ba, copy='spat.z',fixed=F)  + f(score_3, model=spde_score_3 )

# formula <- Y.log ~  -1  +
#   Intercept1 + f(spat.cnt, model=spde.spat)  + f(score_1, model=spde_score_1 )+
#   Intercept2 + f(spat.z, model=spde.spat)  + f(score_2, model=spde_score_2 )+
#   Intercept3 + f(spat.ba, copy='spat.z',fixed=F)  + f(score_3, model=spde_score_3 )
# 

# formula <- Y.log ~  -1  +
#   Intercept1 + f(spat.cnt, model=spde.spat)  + f(score_1, model=spde_score_1 )+
#   Intercept2 + f(spat.z, model=spde.spat)  + f(score_2, model=spde_score_2 )+
#   Intercept3 + f(spat.ba, model=spde.spat)  + f(score_3, model=spde_score_3 )
# 

# formula <- Y.log ~  -1  +
#   Intercept1 + f(spat.cnt, copy='spat.z',fixed=F, group=spat.cnt.group)  + f(score_1, model=spde_score_1 )+
#   Intercept2 + f(spat.z, model=spde.spat, group = spat.z.group, control.group = list(model = "ar", order=1))  + f(score_2, model=spde_score_2 )+
#   Intercept3 + f(spat.ba, copy='spat.z',fixed=F, group=spat.ba.group)  + f(score_3, model=spde_score_3 )
# 

res <- inla(formula,
            family = c('poisson','binomial', 'weibull'), data = inla.stack.data(all.stack),  Ntrials=1,
            control.predictor = list(A = inla.stack.A(all.stack), compute = TRUE, link=c(rep(1,n1),rep(2,n1),rep(3,n1))),
            verbose=TRUE,
            control.compute=list(config = TRUE),
            control.family = list( list(), list(), list()),
            control.fixed = list(expand.factor.strategy = 'inla')
)

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
# t2 <- Sys.time()
# print(t2-t1)
summary(res)


# save(res, file=file.path(dir.out, 'Model_gamma_spde_0.1.RData'))
save(res, file=file.path(dir.out, 'Model_weibull_spde_0.1.RData'))
# save(res, file=file.path(dir.out, 'Model_egp_spde_0.1.RData'))

