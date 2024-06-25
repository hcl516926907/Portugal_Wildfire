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

covar.names <- c(covar.names ,'LVegTyp','HVegTyp','month.idx','lon.grid','lat.grid')
library(randomForest)

data.rf <- data.fit2[data.fit2$time.idx <= 84, c(covar.names,'area_ha','z')]
# for (var in c('LVegTyp','HVegTyp','month.idx')){
#   data.rf[,var] <- as.factor(data.rf[,var])
# }


data.rf$area_ha <- log(data.rf$area_ha)
data.rf[is.na(data.rf$area_ha),'area_ha'] <- 0
data.rf$z <- as.factor(data.rf$z)
set.seed(100)

train <- sample(nrow(data.rf), nrow(data.rf), replace = FALSE)
TrainSet1 <- data.rf[train, c(covar.names,'z')]

model.z <- randomForest(z ~ ., data = TrainSet1, importance = TRUE)
model.z

weight <- ifelse( data.rf$z==1, 3, 1 )
model.z1 <- randomForest(z ~ ., data = TrainSet1, weights=weight, importance = TRUE)
model.z1
data.rf$score_z <- as.vector(predict(model.z, data.rf,type = "prob")[,2])

# 1. Area under curve


TrainSet2 <- data.rf[train, c(covar.names, 'score_z','area_ha')]
model.ba <- randomForest(area_ha ~ ., data = TrainSet2, importance = TRUE)
model.ba
# # 
# save(model1,model2, file=file.path(dir.out,'RF_model.RData'))

# load(file=file.path(dir.out,'RF_model.RData'))
data.fit2$score_z <-  as.vector(predict(model.z, data.fit2,type = "prob")[,2])
data.fit2$score_ba <-  predict(model.ba, data.fit2)


for (var in c('score_ba','score_z') ){
  data.fit2[,var] <- (data.fit2[,var]-mean(data.fit2[,var]))/sd(data.fit2[,var])
}

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
baNA1 = as.vector(c(nothing1, nothing2, data.fit2$area_ha))
baNA2 = as.vector(c(nothing1, nothing2, sqrt(data.fit2$area_ha)))
baNA3 = as.vector(c(nothing1, nothing2, (data.fit2$area_ha)^{1/3}))
baNA4 = as.vector(c(nothing1, nothing2, (data.fit2$area_ha)^{1/4}))

outcome.matrix.log <- matrix(c(cntNA,zNA, baNA.log), ncol=3)
outcome.matrix1 <- matrix(c(cntNA,zNA, baNA1), ncol=3)
outcome.matrix2 <- matrix(c(cntNA,zNA, baNA2), ncol=3)
outcome.matrix3 <- matrix(c(cntNA,zNA, baNA3), ncol=3)
outcome.matrix4 <- matrix(c(cntNA,zNA, baNA4), ncol=3)



Intercept1 <- c(rep(1,n1),nothing1, nothing2)
Intercept2 <- c(nothing1, rep(1,n1), nothing2)
Intercept3 <- c(nothing1, nothing2, rep(1,n1))


i.spat1 = c(data.fit2$grid.idx, nothing1, nothing2)# fire ignition
i.spat.iid = c(data.fit2$grid.idx, nothing1, nothing2)# fire ignition
i.spat2 = c(nothing1,data.fit2$grid.idx, nothing2)# BA ind
i.spat3 = c(nothing1, nothing2, data.fit2$grid.idx)# BA
i.spat3.iid = c(nothing1, nothing2, data.fit2$grid.idx)# BA


time.idx1 <- c(data.fit2$month.idx, nothing1, nothing2)# fire ignition
time.idx2 <- c(nothing1, data.fit2$month.idx, nothing2)# BA ind
time.idx3  <-  c(nothing1, nothing2, data.fit2$month.idx)# BA


RHumi_1 <- c(data.fit2$RHumi,nothing1, nothing2)
RHumi_2 <- c(nothing1,  data.fit2$RHumi, nothing2)
RHumi_3 <- c(nothing1, nothing2, data.fit2$RHumi)

RHumi_grp <- inla.group(data.fit2$RHumi, n = 20, method = "quantile")
RHumi_grp_1 <- c(RHumi_grp, nothing1, nothing2)
RHumi_grp_2 <- c(nothing1, RHumi_grp, nothing2)
# RHumi_grp_3 <- c(nothing1, nothing2, RHumi_grp)
RHumi_grp_3 <- c(nothing1, nothing2, RHumi_grp)
# 
Temp_1 <- c(data.fit2$Temp,nothing1, nothing2)
Temp_2 <- c(nothing1,  data.fit2$Temp, nothing2)
Temp_3 <- c(nothing1, nothing2, data.fit2$Temp)

Temp_grp <- inla.group(data.fit2$Temp, n = 20, method = "quantile")
Temp_grp_1 <- c(Temp_grp,nothing1, nothing2)
Temp_grp_2 <- c(nothing1,  Temp_grp, nothing2)
Temp_grp_3 <- c(nothing1, nothing2, Temp_grp)

# 
UComp_1 <- c(data.fit2$UComp,nothing1, nothing2)
UComp_2 <- c(nothing1,  data.fit2$UComp, nothing2)
UComp_3 <- c(nothing1, nothing2, data.fit2$UComp)

UComp_grp <- inla.group(data.fit2$UComp, n = 20, method = "quantile")
UComp_grp_1 <- c(UComp_grp,nothing1, nothing2)
UComp_grp_2 <- c(nothing1,  UComp_grp, nothing2)
UComp_grp_3 <- c(nothing1, nothing2, UComp_grp)

VComp_1 <- c(data.fit2$VComp,nothing1, nothing2)
VComp_2 <- c(nothing1,  data.fit2$VComp, nothing2)
VComp_3 <- c(nothing1, nothing2, data.fit2$VComp)

VComp_grp <- inla.group(data.fit2$VComp, n = 20, method = "quantile")
VComp_grp_1 <- c(VComp_grp,nothing1, nothing2)
VComp_grp_2 <- c(nothing1,  VComp_grp, nothing2)
VComp_grp_3 <- c(nothing1, nothing2, VComp_grp)

score_z_grp <- inla.group(data.fit2$score_z, n = 20, method = "quantile")
score_2 <- c(nothing1,  data.fit2$score_z, nothing2)
score_2_grp <- c(nothing1, score_z_grp, nothing2)

score_ba_grp <- inla.group(data.fit2$score_ba, n = 40, method = "quantile")
score_ba_grp[which(score_ba_grp<0.1)] <- 0
score_3 <- c(nothing1, nothing2, data.fit2$score_ba)
score_3_grp <- c(nothing1, nothing2, score_ba_grp)

inner.length <- 0.7
outer.length <- 1
Cov.mesh <- inla.mesh.2d(loc.domain = data.fit2[,c('UComp','VComp')], max.edge=c(inner.length,outer.length))
plot(Cov.mesh)
points(data.fit2[,c('UComp','VComp')],pch=19, col='blue', cex=0.1)

plot(Cov.mesh)
points(data.fit2[data.fit2$y>0,c('UComp','VComp')],pch=19, col='blue', cex=0.1)

A.Cov <- inla.spde.make.A(mesh=Cov.mesh, loc=as.matrix(data.fit2[,c('UComp','VComp')]))
SPDE.Cov <- inla.spde2.pcmatern(mesh = Cov.mesh, 
                                prior.range = c(0.9, 0.3), # P(range < 50) = 0.1
                                prior.sigma = c(1, 0.1)) # P(sigma > 1) = 0.01
s.index <-  inla.spde.make.index(name = "Wind",
                                 n.spde = SPDE.Cov$n.spde)

cnt.stack <- inla.stack(
  data= list(Y=cbind(cnt,NA,NA)),
  A <- list(1,1,1),
  effect = list(Intercept1=rep(1,nrow(data.fit2)), idarea1=data.fit2$grid.idx, time.idx1=data.fit2$time.idx ),
  tag='pois'
)

z.stack <- inla.stack(
  data= list(Y=cbind(NA,z,NA)),
  A <- list(1,1,1),
  effect = list(Intercept2=rep(1,nrow(data.fit2)), idarea2=data.fit2$grid.idx, time.idx2=data.fit2$time.idx ),
  tag='bern'
)

ba.stack <- inla.stack(
  data= list(Y=cbind(NA,NA,(data.fit2$area_ha)^{1/3})),
  A <- list(1,1,1, A.Cov),
  effect = list(Intercept3=rep(1,nrow(data.fit2)), idarea3=data.fit2$grid.idx, time.idx3=data.fit2$time.idx, i.cov = Cov.mesh$n ),
  tag='BYM'
)

all.stack <- inla.stack(cnt.stack, z.stack, ba.stack )

data=list(Y.log=outcome.matrix.log, 
          Y1 =outcome.matrix1,
          Y2 =outcome.matrix2,
          Y3 =outcome.matrix3,
          Y4 =outcome.matrix4,
          idarea1=i.spat1, idarea2=i.spat2, idarea3=i.spat3,
                            idarea1.iid = i.spat1, idarea3.iid = i.spat3,
                    time.idx1 = time.idx1,
                    time.idx2 = time.idx2,
                    # time.idx2.group = time.idx2.group,
                    time.idx3 = time.idx3,
                    # time.idx3.group = time.idx3.group,
                    Intercept1 = Intercept1,
                    Intercept2 = Intercept2,
                    Intercept3 = Intercept3,
                    RHumi_1 = RHumi_1,
                    RHumi_2 = RHumi_2,
                    RHumi_3 = RHumi_3,
                    # 
                    RHumi_grp_1 = RHumi_grp_1,
                    RHumi_grp_2 = RHumi_grp_2,
                    RHumi_grp_3 = RHumi_grp_3,
          
                    Temp_1 = Temp_1,
                    Temp_2 = Temp_2,
                    Temp_3 = Temp_3,
          
                    Temp_grp_1 = Temp_grp_1,
                    Temp_grp_2 = Temp_grp_2,
                    Temp_grp_3 = Temp_grp_3,

                    UComp_1 = UComp_1,
                    UComp_2 = UComp_2,
                    UComp_3 = UComp_3,
          
                    UComp_grp_1 = UComp_grp_1,
                    UComp_grp_2 = UComp_grp_2,
                    UComp_grp_3 = UComp_grp_3,
          
                    VComp_1 = VComp_1,
                    VComp_2 = VComp_2,
                    VComp_3 = VComp_3,
          
                    VComp_grp_1 = VComp_grp_1,
                    VComp_grp_2 = VComp_grp_2,
                    VComp_grp_3 = VComp_grp_3,
          
                    score_2_grp = score_2_grp,
                    score_3_grp = score_3_grp,
          
                    score_2 = score_2,
                    score_3 = score_3
          
)


formula0 <- Y3 ~  -1  + Intercept1 +  f(idarea1, model='bym2',graph = g, group=time.idx1) +
  Intercept2 +  f(idarea2, model='bym2',graph = g, group=time.idx2) +
  Intercept3 + f(idarea3,  copy='idarea2',fixed=F,group=time.idx3) 

t1 <- Sys.time()
res0 <- inla(formula0,
             family = c('poisson','binomial', 'gp'), data = data,  Ntrials=1,
             control.predictor = list(compute = TRUE, link=c(rep(1,n1),rep(2,n1),rep(3,n1))),
             verbose=TRUE,
             control.compute=list(config = TRUE),
             control.family = list( list(), list(), list(control.link=list(quantile=0.5))),
             control.fixed = list(expand.factor.strategy = 'inla')
)
t2 <- Sys.time()
print(t2-t1)
summary(res0)

save(res0, file=file.path(dir.out,'Temp0.RData'))

formula1.1 <- Y3 ~  -1 + Intercept1 +  f(idarea1, model='bym2',graph = g, group=time.idx1) +
  
  Intercept2 + f(idarea2, model='bym2',graph = g,group=time.idx2) + 
  f(score_2, model='linear', mean.linear=1, prec.linear=1) +
  
  Intercept3 + f(idarea3,  copy='idarea2',fixed=F, group=time.idx3) +
  f(score_3, model='linear',mean.linear=1, prec.linear=1)


t1 <- Sys.time()
res1.1 <- inla(formula1.1,
             family = c('poisson','binomial', 'gp'), data = data,  Ntrials=1,
             control.predictor = list(compute = TRUE, link=c(rep(1,n1),rep(2,n1),rep(3,n1))),
             verbose=TRUE,
             control.compute=list(config = TRUE),
             control.family = list( list(), list(), list(control.link=list(quantile=0.5))),
             control.fixed = list(expand.factor.strategy = 'inla')
)
t2 <- Sys.time()
print(t2-t1)
summary(res1.1)

hyper.bym2 = list (phi = list (
  prior = "pc",
  param = c(0.5 , 0.01) ,
  initial = 0.8) ,
  prec = list (
    prior = "pc.prec",
    param = c (6.5 , 0.9) ,
    initial = 5))

formula1.2 <- Y3 ~  -1 + Intercept1 +  f(idarea1, model='bym2',graph = g, group=time.idx1) +
  Intercept2 +  f(idarea2, model='bym2',graph = g,group=time.idx2, hyper=hyper.bym2) + 
  Intercept3 +  f(idarea3,  copy='idarea2',fixed=F, group=time.idx3) +
  f(score_3_grp, model='rw1', hyper = list(theta = list(prior="pc.prec", param=c(1,0.01))))


t1 <- Sys.time()
res1.2 <- inla(formula1.2,
               family = c('poisson','binomial', 'gp'), data = data,  Ntrials=1,
               control.predictor = list(compute = TRUE, link=c(rep(1,n1),rep(2,n1),rep(3,n1))),
               verbose=TRUE,
               control.compute=list(config = TRUE),
               control.family = list( list(), list(), list(control.link=list(quantile=0.5))),
               control.fixed = list(expand.factor.strategy = 'inla')
)
t2 <- Sys.time()
print(t2-t1)
summary(res1.2)
save(res1.2, file=file.path(dir.out,'Temp1.2.RData'))



formula1.3 <- Y.log ~  -1 + Intercept1 +  f(idarea1, model='bym2',graph = g, group=time.idx1) +
  Intercept2 +  f(idarea2, model='bym2',graph = g,group=time.idx2, hyper=hyper.bym2) + 
  f(score_2_grp, model='rw1',hyper = list(theta = list(prior="pc.prec", param=c(0.1,0.01))))+
  Intercept3 +  f(idarea3,  copy='idarea2',fixed=F, group=time.idx3) +
  f(score_3_grp, model='rw1', hyper = list(theta = list(prior="pc.prec", param=c(1,0.01))))


t1 <- Sys.time()
res1.3 <- inla(formula1.3,
               family = c('poisson','binomial', 'gp'), data = data,  Ntrials=1,
               control.predictor = list(compute = TRUE, link=c(rep(1,n1),rep(2,n1),rep(3,n1))),
               verbose=TRUE,
               control.compute=list(config = TRUE),
               control.family = list( list(), list(), list(control.link=list(quantile=0.5))),
               control.fixed = list(expand.factor.strategy = 'inla')
)
t2 <- Sys.time()
print(t2-t1)
summary(res1.3)

save(res1.3, file=file.path(dir.out,'Temp1.3.RData'))


formula2 <- Y ~  -1 + Intercept1 +  f(idarea1, model='bym2',graph = g)+ f(time.idx1, model='ar1') +
  Intercept2 +  f(idarea2, model='bym2',graph = g) + f(time.idx2, model='ar1') +
  Intercept3 +  f(idarea3,  copy='idarea2',fixed=F) + f(time.idx3, model='ar1') + f(i.cov, model=SPDE.Cov)


t1 <- Sys.time()
res2 <- inla(formula2,
             family = c('poisson','binomial', 'gp'), data = inla.stack.data(all.stack),  Ntrials=1,
             control.predictor = list(A=inla.stack.A(all.stack),compute=FALSE),
             verbose=TRUE,
             control.compute=list(config = TRUE),
             control.family = list( list(), list(), list(control.link=list(quantile=0.5))),
             control.fixed = list(expand.factor.strategy = 'inla')
)
t2 <- Sys.time()
print(t2-t1)
summary(res2)




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


library(doParallel)
cl <- makeCluster(3)
registerDoParallel(cl)

n.samples = 200


# idx.pred.pois <- 1:n1
# idx.pred.z <- (n1+1):(2*n1)
# idx.pred.ba <- (2*n1+1):(3*n1)
# result <- res0
# 
# 
# # idx.pred.pois <- 1:n1
# # idx.pred.ba <- (n1+1):(2*n1)
# # result <- res0
# 
# samples = inla.posterior.sample(n.samples, result = result, seed=1234)
# 
# 
# # pred.sp <- post.pred.gamma(samples, idx.pred.pois, idx.pred.z, idx.pred.ba)
# t1 <- Sys.time()
# 
# pred.sp <- post.pred.gpd.parallel(samples, idx.pred.pois, idx.pred.z, idx.pred.ba, alpha=0.5, n.samples=n.samples)
# t2 <- Sys.time()
# print(t2-t1)

load(file=file.path(dir.out,'Temp1.3.RData'))
load(file=file.path(dir.out,'Temp1.3_pred.sp.RData'))
pred.cnt <- pred.sp$pred.cnt
pred.ba <- pred.sp$pred.ba
pred.z <- pred.sp$pred.z
result <- res1.3

idx.pred.pois <- 1:n1
idx.pred.z <- (n1+1):(2*n1)
idx.pred.ba <- (2*n1+1):(3*n1)

crps.val <- rep(NA, nrow(data.fit2))
for (i in 1:nrow(data.fit2) ){
  y <- log(data.fit2$area_ha[i])
  if (is.na(y)) y <- 0
  crps.val[i] <- crps_sample(
    y,
    pred.ba[[i]],
    method = "edf")
}

ins.idx <- 192*84
round(sum(crps.val[1:ins.idx])/ins.idx,4)
round(sum(crps.val[(ins.idx+1):nrow(data.fit2)])/(nrow(data.fit2)-ins.idx),4)
round(sum(crps.val)/length(crps.val),4)

data.fit2$crps <- crps.val
############## vector name in the samples$latent:#################################
# Predicator:1, ... Predicator:172 means linear predicator of grid_1_1:9, 
# grid_2_1:9,...grid_192_1:9

# grid.idx:1,...grid.idx5184:  u+v of grid_1:192 in group1, u of grid_1:192 in group1
                             # u+v of grid_1:192 in group2, u of grid_1:192 in group2 ...

# idx.pred.pois <- 1:n1
# idx.pred.ba <- (n1+1):(2*n1)
# idx.latent.id1 <- 3457:6912
# idx.latent.id2 <- 6913:10368
# idx.intercept1 <- 10369
# idx.intercept2 <- 10370

# group1 <- 1729: (1728 + 192*2)
# group2 <- 2113: (2112 + 192*2)
# group1 <- 2497: (2496 + 192*2)
# 
# samples[[1]]$latent[idx.pred.pois,1][1:10]
# samples[[1]]$latent[5185,1] # intercept
# 
# samples[[1]]$latent[group1,1][1] + samples[[1]]$latent[5185,1]
# samples[[1]]$latent[group2,1][1] + samples[[1]]$latent[5185,1]
# samples[[1]]$latent[group3,1][1] + samples[[1]]$latent[5185,1]
# 
# 






data.fit2[,'z'] <- as.vector((data.fit2$y>0)+0)
data.fit2[,'Latent_Effect_Cnt'] <- result$summary.fitted.values$mean[idx.pred.pois]
data.fit2[,'Latent_Effect_BA'] <- result$summary.fitted.values$mean[idx.pred.ba]
data.fit2[,'Latent_Effect_Z'] <- result$summary.fitted.values$mean[idx.pred.z]

data.fit2[,'Estimated_Lambda'] <- sapply(pred.cnt,mean)
data.fit2$Lower_Bound <- sapply(pred.cnt,quantile,0.025)
data.fit2$Upper_Bound <- sapply(pred.cnt,quantile,0.975)
data.fit2$Scaling_Residuals <- data.fit2$y - data.fit2$Estimated_Lambda
summary(data.fit2$Scaling_Residuals)

data.fit2[,'Estimated_Z'] <- sapply(pred.z,mean)

data.fit2[,'Estimated_BA'] <- sapply(pred.ba,mean)
data.fit2[,'Lower_BA'] <- sapply(pred.ba,quantile,0.025)
data.fit2[,'Upper_BA'] <- sapply(pred.ba,quantile,0.975)




# data.fit2$BA_Residuals <- log(data.fit2$area_ha) - data.fit2$Estimated_BA
# summary(data.fit2$BA_Residuals)



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
  labs(title = "Estimated Burn area with 95% credible band in Year 2017 (with predicator) ",
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
