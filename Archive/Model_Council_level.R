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
library(units)

load(file.path(dir.out, 'XGBoost_Score_Council_2.RData'))

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

save(map, map.district,  file=file.path(dir.out, 'map_council_district.RData'))

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
load(file=file.path(dir.out, 'districts_prediction.RData'))

data.fit <- merge(data.fit, data.fit.dist.pred, by=c('NAME_1','time.idx'))
# 
# 
# 
data.fit <- merge(data.fit,st_drop_geometry(map.district[,c('NAME_1','grid.idx.district')]),by='NAME_1')

###########################################################################
n.group <- 40


data.fit$z <- as.vector((data.fit$y>0)+0)

log.ba.agg <- data.fit[data.fit$year<=2022,] %>% group_by(NAME_2) %>%
  summarise(cum_log_ba=sum(log_ba,na.rm=TRUE)) %>%
  mutate(log.ba.idx=inla.group(cum_log_ba,n=25,method='quantile',idx.only=T))

cnt.agg <- data.fit[data.fit$year<=2022,] %>% group_by(NAME_2) %>%
  summarise(cum_cnt=sum(y,na.rm=TRUE)) %>%
  mutate(cnt.idx=inla.group(cum_cnt,n=25,method='quantile',idx.only=T))

z.agg <- data.fit[data.fit$year<=2022,] %>% group_by(NAME_2) %>%
  summarise(cum_z=sum(z,na.rm=TRUE)) %>%
  mutate(z.idx=inla.group(cum_z,n=25,method='quantile',idx.only=T))


data.fit <- merge(data.fit,log.ba.agg[,c("NAME_2",'log.ba.idx')],by='NAME_2')
data.fit <- merge(data.fit, cnt.agg[,c('NAME_2','cnt.idx')], by='NAME_2')
data.fit <- merge(data.fit, z.agg[,c('NAME_2','z.idx')], by='NAME_2')


# high.risk.log.ba <- data.fit[data.fit$year<=2022,] %>% group_by(NAME_2) %>%
#   summarise(cum_log_ba=sum(log_ba,na.rm=TRUE)) %>%
#   mutate(bins=cut(cum_log_ba,breaks=c(-1,quantile(cum_log_ba,0.7),max(cum_log_ba))))
# map.log.ba <- st_make_valid(map)
# map.log.ba <- merge(map.log.ba,high.risk.log.ba,by='NAME_2')
# map.log.ba.dissolve <- map.log.ba %>% group_by(bins) %>% summarise(do_union = TRUE) %>%
#   st_cast("POLYGON") %>%  # Separate multipart polygons into individual polygons
#   mutate(polygon_id = row_number())
# 
# area_threshold <- set_units(20000000, "m^2")
# 
# # Calculate areas and filter
# map.log.ba.dissolve.1 <- map.log.ba.dissolve %>%
#   filter(st_area(.) > area_threshold)
# 
# 
# nb.diss.ba <- poly2nb(map.log.ba.dissolve.1)
# map.log.ba.dissolve.1$diss.idx.ba <- 1:nrow(map.log.ba.dissolve.1)
# centroids <- st_centroid(map.log.ba.dissolve.1)
# centroid_coords <- st_coordinates(centroids)
# map.log.ba.dissolve.1$lon.grid <- centroid_coords[,1]
# map.log.ba.dissolve.1$lat.grid <- centroid_coords[,2]
# 
# ggplot(map.log.ba.dissolve.1) + geom_sf()+
#   geom_text(aes(lon.grid, lat.grid, label = diss.idx.ba), size=2,color='red')
# 
# nb2INLA("map.diss.ba.adj", nb.diss.ba)
# g.diss.ba <- inla.read.graph(filename = "map.diss.ba.adj")
# 
# 
# 
# 
# 
# high.risk.z <- data.fit[data.fit$year<=2022,] %>% group_by(NAME_2) %>%
#   summarise(cum_z=sum(z,na.rm=TRUE)) %>%
#   mutate(bins=cut(cum_z,breaks=c(-1,quantile(cum_z,0.7),max(cum_z))))
# map.z <- st_make_valid(map)
# map.z <- merge(map.z,high.risk.z,by='NAME_2')
# map.z.dissolve <- map.z %>% group_by(bins) %>% summarise(do_union = TRUE) %>%
#   st_cast("POLYGON") %>%  # Separate multipart polygons into individual polygons
#   mutate(polygon_id = row_number())
# 
# area_threshold <- set_units(20000000, "m^2")
# 
# # Calculate areas and filter
# map.z.dissolve.1 <- map.z.dissolve %>%
#   filter(st_area(.) > area_threshold)
# # 
# 
# 
# 
# 
# 
# nb.diss.z <- poly2nb(map.z.dissolve.1)
# map.z.dissolve.1$diss.idx.z <- 1:nrow(map.z.dissolve.1)
# 
# centroids <- st_centroid(map.z.dissolve.1)
# centroid_coords <- st_coordinates(centroids)
# map.z.dissolve.1$lon.grid <- centroid_coords[,1]
# map.z.dissolve.1$lat.grid <- centroid_coords[,2]
# 
# ggplot(map.z.dissolve.1) + geom_sf()+
#   geom_text(aes(lon.grid, lat.grid, label = diss.idx.z), size=2,color='red')
# 
# nb2INLA("map.diss.z.adj", nb.diss.z)
# g.diss.z <- inla.read.graph(filename = "map.diss.z.adj")
# 
# 
# 
# 
# high.risk.cnt <- data.fit[data.fit$year<=2022,] %>% group_by(NAME_2) %>%
#   summarise(cum_cnt=sum(y,na.rm=TRUE)) %>%
#   mutate(bins=cut(cum_cnt,breaks=c(-1,quantile(cum_cnt,0.7),max(cum_cnt))))
# map.cnt <- st_make_valid(map)
# map.cnt <- merge(map.cnt,high.risk.cnt,by='NAME_2')
# map.cnt.dissolve <- map.cnt %>% group_by(bins) %>% summarise(do_union = TRUE) %>%
#   st_cast("POLYGON") %>%  # Separate multipart polygons into individual polygons
#   mutate(polygon_id = row_number())
# 
# area_threshold <- set_units(20000000, "m^2")
# 
# # Calculate areas and filter
# map.cnt.dissolve.1 <- map.cnt.dissolve %>%
#   filter(st_area(.) > area_threshold)
# # 
# 
# nb.diss.cnt <- poly2nb(map.cnt.dissolve.1)
# map.cnt.dissolve.1$diss.idx.cnt <- 1:nrow(map.cnt.dissolve.1)
# 
# centroids <- st_centroid(map.cnt.dissolve.1)
# centroid_coords <- st_coordinates(centroids)
# map.cnt.dissolve.1$lon.grid <- centroid_coords[,1]
# map.cnt.dissolve.1$lat.grid <- centroid_coords[,2]
# 
# ggplot(map.cnt.dissolve.1) + geom_sf()+
#   geom_text(aes(lon.grid, lat.grid, label = diss.idx.cnt), size=2,color='red')
# 
# nb2INLA("map.diss.cnt.adj", nb.diss.cnt)
# g.diss.cnt <- inla.read.graph(filename = "map.diss.cnt.adj")
# 
# 
# 
# data.fit <-  data.fit |> 
#   st_as_sf(coords = c("lon.grid", "lat.grid"), crs = 4326, remove=FALSE) |>
#   st_join(map.cnt.dissolve.1[,'diss.idx.cnt'], join = st_within) |>
#   st_drop_geometry()
# 
# data.fit <-  data.fit |> 
#   st_as_sf(coords = c("lon.grid", "lat.grid"), crs = 4326, remove=FALSE) |>
#   st_join(map.z.dissolve.1[,'diss.idx.z'], join = st_within) |>
#   st_drop_geometry()
# 
# data.fit <-  data.fit |> 
#   st_as_sf(coords = c("lon.grid", "lat.grid"), crs = 4326, remove=FALSE) |>
#   st_join(map.log.ba.dissolve.1[,'diss.idx.ba'], join = st_within) |>
#   st_drop_geometry()

data.fit <- data.fit[order(data.fit$time.idx,data.fit$grid.idx),]
 
z <- as.vector((data.fit$y>0)+0)
log.ba <- as.vector(ifelse(data.fit$y>0, exp(data.fit$log_ba)^{1/2}, NA))
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


rw.group <- function(x1,x2,n){
  q.bins <- c(0, ppoints(n-1), 1)
  q.bins.value <- lead((q.bins+ lag(q.bins))/2)[1:n]
  bins <- unique(quantile(x1, probs = q.bins))
  bins.value <-  unique(quantile(x1, probs = q.bins.value))
  bins[1] <- min(bins[1], min(x2))
  bins[length(bins)] <- max(bins[length(bins)], max(x2))
  x2.bins <- cut(x2, breaks = bins, include.lowest = TRUE)
  levels(x2.bins) <- bins.value
  x2.bins <- as.numeric(as.character(x2.bins))
  return(x2.bins)
}


score_cnt_grp <- rw.group(data.fit[data.fit$y>0,'score_cnt'],data.fit$score_cnt,n.group)
score.cnt <- c(score_cnt_grp,nothing, nothing)

score_z_grp <- rw.group(data.fit$score_z,data.fit$score_z,n.group)
score.z <- c(nothing, score_z_grp, nothing)

score_ba_grp <- rw.group(data.fit[data.fit$y>0,'score_ba'],data.fit$score_ba,n.group)
score.ba <- c(nothing, nothing, score_ba_grp)


score_cnt_dist_group <-  rw.group(data.fit[data.fit$y>0,'score_cnt.dist'],data.fit$score_cnt.dist,n.group)
score.cnt.dist <- c(score_cnt_dist_group, nothing, nothing)

score_z_dist_group <-  rw.group(data.fit$score_z.dist,data.fit$score_z.dist,n.group)
score.z.dist <- c(nothing, score_z_dist_group, nothing)

score_ba_dist_group <- rw.group(data.fit[data.fit$y>0,'score_ba.dist'],data.fit$score_ba.dist,n.group)
score.ba.dist <- c(nothing, nothing, score_ba_dist_group)

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


log.ba.idx <- c(nothing, nothing,  data.fit$log.ba.idx)
cnt.idx <- c(data.fit$cnt.idx, nothing, nothing)
z.idx <- c(nothing, data.fit$z.idx, nothing)
# 
# diss.idx.cnt <- c(data.fit$diss.idx.cnt, nothing, nothing)
# diss.idx.z <- c(nothing, data.fit$diss.idx.z, nothing)
# diss.idx.ba <- c(nothing, nothing, data.fit$diss.idx.ba)
print(summary(baNA.log))
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
          
          score.cnt.dist = score.cnt.dist,
          score.z.dist = score.z.dist,
          score.ba.dist = score.z.dist,
          
          score.cnt = score.cnt,
          score.z = score.z,
          score.ba = score.ba,
          
          log.ba.idx = log.ba.idx,
          cnt.idx = cnt.idx,
          z.idx = z.idx
          # 
          # diss.idx.cnt = diss.idx.cnt,
          # diss.idx.z = diss.idx.z,
          # diss.idx.ba = diss.idx.ba
          
)




# hyper.theta <-  list(prec = list(prior="loggamma",param=c(1,0.1)))


formula <- Y.log ~  0  +
  intercept.cnt + 
  # year.af.2017.cnt +
  f(grid.idx.cnt, copy='grid.idx.z',fixed=F)  +
  f(year.idx.cnt, model='iid')+ f(score.cnt, model='rw1' )+
  f(grid.idx.dist.cnt, copy="grid.idx.dist.z",fixed=F, group=time.idx.cnt, control.group = list(model='iid'))  +


  intercept.z + 
  # year.af.2017.z +
  f(grid.idx.z, model='bym2',graph=g.council)  +
  f(year.idx.z, model='iid')+  f(score.z, model='rw1' )+
  f(grid.idx.dist.z, model='bym2', graph=g.district,group=time.idx.z, control.group = list(model='iid')) +

  intercept.ba + 
  # year.af.2017.ba +
  f(grid.idx.ba, copy='grid.idx.z',fixed=F)  +
  f(year.idx.ba, model='iid') +  f(score.ba, model='rw1' ) +
  f(grid.idx.dist.ba, copy="grid.idx.dist.z",fixed=F, group=time.idx.ba, control.group = list(model='iid'))


# formula <- Y.log ~  0  +
#   intercept.cnt + year.af.2017.cnt +
#   f(grid.idx.cnt, copy='grid.idx.z',fixed=F)  +
#   f(year.idx.cnt, model='iid')+ f(score.cnt, model='rw1' )+
#   f(grid.idx.dist.cnt, copy="grid.idx.dist.z",fixed=F, group=time.idx.cnt, control.group = list(model='iid'))  +
#   
#   
#   intercept.z + year.af.2017.z +
#   f(grid.idx.z, model='bym2',graph=g.council)  +
#   f(year.idx.z, model='iid')+  f(score.z, model='rw1' )+
#   f(grid.idx.dist.z, model='bym2', graph=g.district,group=time.idx.z, control.group = list(model='iid')) +
#   
#   intercept.ba + year.af.2017.ba +
#   f(grid.idx.ba, model='bym2',graph=g.council)  +
#   f(year.idx.ba, model='iid') +  f(score.ba, model='rw1' ) +
#   f(grid.idx.dist.ba, model='bym2', graph=g.district, group=time.idx.ba, control.group = list(model='iid'))
# 
# 



# save.image(file=file.path(dir.out, 'tmp.RData'))
# load(file.path(dir.out,'tmp.RData'))
t1 <- Sys.time()
res <- inla(formula,
               family = c('nzpoisson','binomial', 'gamma'), data = data.list,  Ntrials=1,
               control.predictor = list(compute = TRUE, link=c(rep(1,n1),rep(2,n1),rep(3,n1))),
               verbose=TRUE,
               control.compute=list(config = TRUE),
               control.family = list( list(), list(), list(control.link=list(quantile=0.5))),
               control.fixed = list(expand.factor.strategy = 'inla')
)

# res <- inla(formula,
#             family = c('nzpoisson','binomial', 'egp'), data = data.list,  Ntrials=1,
#             control.predictor = list(compute = TRUE, link=c(rep(1,n1),rep(2,n1),rep(3,n1))),
#             
#             verbose=TRUE,
#             control.compute=list(config = TRUE),
#             control.family = list( list(), list(), list(control.link = list(quantile = 0.5),
#                                                         hyper = list(
#                                                           tail = list(
#                                                           ##initial = xi.intern,
#                                                           fixed = !TRUE,
#                                                           prior = "pc.egptail",
#                                                           param = c(10, -0.5, 0.5)),
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

# save(res, file=file.path(dir.out, 'Model_egp_bym2.RData'))
# save(res, file=file.path(dir.out, 'Model_weibull_bym2.RData'))
save(res, file=file.path(dir.out, 'Model_gamma_bym2.RData'))

