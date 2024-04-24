dir.data <- "/home/pgrad2/2448355h/My_PhD_Project/00_Dataset/Urban_Fires"
dir.out <- '/home/pgrad2/2448355h/My_PhD_Project/01_Output/Urban_Fire_Model'


library(inlabru)
library(INLA)
library(mgcv)
library(ggplot2)
library(rgdal)
library(sf)
library(fmesher)
library(dplyr)
library(RColorBrewer)
library(terra)
library(rgeos)
# bru_safe_sp(force = TRUE)
load(file.path(dir.data, "burn_area","finaldata_ruralfire_fwi.RData"))
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
data.new$month.idx <- elapsed_months(data.new$open, date.start)

time_cf$date <- paste(time_cf$year,time_cf$month,time_cf$day,sep='/')
time_cf$month.idx <- elapsed_months(time_cf$date, date.start)

loc.data.utm <- st_as_sf(data.new, coords=c('x_utm_new','y_utm_new'), crs='+proj=utm +zone=29 +datum=WGS84 +units=km +no_defs' )


coords <- SpatialPointsDataFrame(data.new, coords=data.new[,c('x_utm_new','y_utm_new')], 
                                 proj4string=CRS('+proj=utm +zone=29 +datum=WGS84 +units=km +no_defs' ))

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

# inner.length <- 20
# outer.length <- 50
# cutoff <- 6

inner.length <- 40
outer.length <- 70
cutoff <- 15

mesh <- inla.mesh.2d(loc.domain = loc.d, max.edge=c(inner.length,outer.length), cutoff=cutoff,
                     crs=projutm)
mesh$n
plot(mesh)
# points(coords, col = 4, pch = 19, cex=0.2)

domain.polys <- Polygons(list(Polygon(loc.d)), '0')
domainSP <- SpatialPolygons(list(domain.polys),proj4string=CRS(projutm))

#no of mesh points...
nv <- mesh$n 

ggplot() +
  gg(mesh) +
  gg(domainSP) +
  # gg(mrsea$samplers) +
  gg(loc.data.utm, size = 0.1) +
  facet_wrap(~year.idx) +
  ggtitle("Fire Occurance by year")



r <- 50
p1 <- 0.1
sigma <- 8
p2 <- 0.1

spde <- inla.spde2.pcmatern(mesh = mesh,
                            # PC-prior on range: P(practic.range < r) = p1
                            prior.range = c(r, p1),
                            # PC-prior on sigma: P(sigma > sigma) = p2
                            prior.sigma = c(sigma, p2))


mesh1D.fwi <- fm_mesh_1d(seq(0,100,by=5), boundary = "free")
spde1D.fwi <- inla.spde2.pcmatern(mesh1D.fwi,
                              prior.range = c(1, 0.1),
                              prior.sigma = c(2, 0.1)
)
ggplot() +
  gg(mesh1D)

mesh1D.wspd <- fm_mesh_1d(seq(0,15,by=0.5), boundary = "free")
spde1D.wspd <- inla.spde2.pcmatern(mesh1D.wspd,
                              prior.range = c(1, 0.5),
                              prior.sigma = c(2, 0.5)
)

mesh1D.pricp <- fm_mesh_1d(seq(0,0.5,by=0.01), boundary = "free")
spde1D.pricp <- inla.spde2.pcmatern(mesh1D.pricp,
                                   prior.range = c(1, 0.5),
                                   prior.sigma = c(2, 0.5)
)

mesh1D.temp <- fm_mesh_1d(seq(0,35,by=1), boundary = "free")
spde1D.temp <- inla.spde2.pcmatern(mesh1D.temp,
                                    prior.range = c(1, 0.5),
                                    prior.sigma = c(2, 0.5)
)

mesh1D.rhumi <- fm_mesh_1d(seq(0,100,by=5), boundary = "free")
spde1D.rhumi <- inla.spde2.pcmatern(mesh1D.rhumi,
                                   prior.range = c(1, 0.5),
                                   prior.sigma = c(2, 0.5)
)

temproal_spdf_fwi <- function(time){
  time <- as.integer(time)
  m <- which(time_cf$month.idx==time)
  # average over the month.idx
  spdf <- SpatialPixelsDataFrame(points = grid_pixels, 
                                 data = data.frame(var=as.vector(rowMeans(fwi_array[,,m],dims=2)),
                                                   time=time))
  proj4string(spdf) <- CRS("EPSG:4326")
  return(spdf)
}


#.data. conatins coordinates in the data and coordiantes of the dual mesh within the boundary.
f.fwi.utm <- function(where) {
  x <- where@coords[,1]
  y <- where@coords[,2]
  time <- where$month.idx
  spp <- SpatialPoints(data.frame(x=x,y=y),
                       proj4string=CRS('+proj=utm +zone=29 +datum=WGS84 +units=km +no_defs' ))
  spp <- spTransform(spp, CRS("EPSG:4326"))
  v <- rep(NA, nrow(where))
  for (t in unique(time)){
    temp_spdf <- temproal_spdf_fwi(t)
    idx <- which(time==t)
    v[idx] <- over(spp[idx,],temp_spdf[,'var'])$var
    if (any(is.na(v[idx]))) {
      v[idx] <- bru_fill_missing(temp_spdf, spp[idx,], v[idx])
    }
  }
  return(v)
}

f.time <- function(where) {
  v <- where$month.idx
  return(v)
}

temproal_spdf_rh <- function(time){
  m <- as.integer(time)
  # average over the month.idx
  spdf <- SpatialPixelsDataFrame(points = grid_pixels, 
                                 data = data.frame(var=as.vector(rhumi.month[,,m]),
                                                   time=time))
  proj4string(spdf) <- CRS("EPSG:4326")
  return(spdf)
}

f.rh.utm <- function(where) {
  x <- where@coords[,1]
  y <- where@coords[,2]
  time <- where$month.idx
  spp <- SpatialPoints(data.frame(x=x,y=y),
                       proj4string=CRS('+proj=utm +zone=29 +datum=WGS84 +units=km +no_defs' ))
  spp <- spTransform(spp, CRS("EPSG:4326"))
  v <- rep(NA, nrow(where))
  for (t in unique(time)){
    temp_spdf <- temproal_spdf_rh(t)
    idx <- which(time==t)
    v[idx] <- over(spp[idx,],temp_spdf[,'var'])$var
    if (any(is.na(v[idx]))) {
      v[idx] <- bru_fill_missing(temp_spdf, spp[idx,], v[idx])
    }
  }
  return(v)
}


temproal_spdf_wspd <- function(time){
  m <- as.integer(time)
  # average over the month.idx
  spdf <- SpatialPixelsDataFrame(points = grid_pixels, 
                                 data = data.frame(var=as.vector(wind.spd.month[,,m]),
                                                   time=time))
  proj4string(spdf) <- CRS("EPSG:4326")
  return(spdf)
}

f.wspd.utm <- function(where) {
  x <- where@coords[,1]
  y <- where@coords[,2]
  time <- where$month.idx
  spp <- SpatialPoints(data.frame(x=x,y=y),
                       proj4string=CRS('+proj=utm +zone=29 +datum=WGS84 +units=km +no_defs' ))
  spp <- spTransform(spp, CRS("EPSG:4326"))
  v <- rep(NA, nrow(where))
  for (t in unique(time)){
    temp_spdf <- temproal_spdf_wspd(t)
    idx <- which(time==t)
    v[idx] <- over(spp[idx,],temp_spdf[,'var'])$var
    if (any(is.na(v[idx]))) {
      v[idx] <- bru_fill_missing(temp_spdf, spp[idx,], v[idx])
    }
  }
  return(v)
}
f.wspd.utm(coords[1:10,])

temproal_spdf_wdir <- function(time){
  m <- as.integer(time)
  # average over the month.idx
  spdf <- SpatialPixelsDataFrame(points = grid_pixels, 
                                 data = data.frame(var=as.vector(wind.dirc.month[,,m]),
                                                   time=time))
  proj4string(spdf) <- CRS("EPSG:4326")
  return(spdf)
}

f.wdir.utm <- function(where) {
  x <- where@coords[,1]
  y <- where@coords[,2]
  time <- where$month.idx
  spp <- SpatialPoints(data.frame(x=x,y=y),
                       proj4string=CRS('+proj=utm +zone=29 +datum=WGS84 +units=km +no_defs' ))
  spp <- spTransform(spp, CRS("EPSG:4326"))
  v <- rep(NA, nrow(where))
  for (t in unique(time)){
    temp_spdf <- temproal_spdf_wdir(t)
    idx <- which(time==t)
    v[idx] <- over(spp[idx,],temp_spdf[,'var'])$var
    if (any(is.na(v[idx]))) {
      v[idx] <- bru_fill_missing(temp_spdf, spp[idx,], v[idx])
    }
  }
  return(v)
}
f.wdir.utm(coords[1:10,])

temproal_spdf_temp <- function(time){
  m <- as.integer(time)
  # average over the month.idx
  spdf <- SpatialPixelsDataFrame(points = grid_pixels, 
                                 data = data.frame(var=as.vector(temp.month[,,m]),
                                                   time=time))
  proj4string(spdf) <- CRS("EPSG:4326")
  return(spdf)
}

f.temp.utm <- function(where) {
  x <- where@coords[,1]
  y <- where@coords[,2]
  time <- where$month.idx
  spp <- SpatialPoints(data.frame(x=x,y=y),
                       proj4string=CRS('+proj=utm +zone=29 +datum=WGS84 +units=km +no_defs' ))
  spp <- spTransform(spp, CRS("EPSG:4326"))
  v <- rep(NA, nrow(where))
  for (t in unique(time)){
    temp_spdf <- temproal_spdf_temp(t)
    idx <- which(time==t)
    v[idx] <- over(spp[idx,],temp_spdf[,'var'])$var
    if (any(is.na(v[idx]))) {
      v[idx] <- bru_fill_missing(temp_spdf, spp[idx,], v[idx])
    }
  }
  return(v)
}
f.temp.utm(coords[1:10,])

temproal_spdf_prcip <- function(time){
  m <- as.integer(time)
  # average over the month.idx
  spdf <- SpatialPixelsDataFrame(points = grid_pixels, 
                                 data = data.frame(var=as.vector(pricp.month[,,m]),
                                                   time=time))
  proj4string(spdf) <- CRS("EPSG:4326")
  return(spdf)
}

f.pricp.utm <- function(where) {
  x <- where@coords[,1]
  y <- where@coords[,2]
  time <- where$month.idx
  spp <- SpatialPoints(data.frame(x=x,y=y),
                       proj4string=CRS('+proj=utm +zone=29 +datum=WGS84 +units=km +no_defs' ))
  spp <- spTransform(spp, CRS("EPSG:4326"))
  v <- rep(NA, nrow(where))
  for (t in unique(time)){
    temp_spdf <- temproal_spdf_prcip(t)
    idx <- which(time==t)
    v[idx] <- over(spp[idx,],temp_spdf[,'var'])$var
    if (any(is.na(v[idx]))) {
      v[idx] <- bru_fill_missing(temp_spdf, spp[idx,], v[idx])
    }
  }
  return(v)
}
f.pricp.utm(coords[1:100,])

system.time(test1 <- f.fwi.utm(coords))


ppxl <- fm_pixels(mesh, mask = domainSP, format = "sp",dims=c(30,90))
ppxl_all <- fm_cprod(ppxl, data.frame(month.idx = seq_len(108)))


# model 1
cmp1 <- coordinates + month.idx ~  Intercept(1) +
  mySmooth(coordinates , model = spde) +
  fwi(f.fwi.utm(.data.), model = "linear")


t1 <-  Sys.time()
fit1 <- lgcp(cmp1, data=coords, samplers = domainSP,
             domain = list(coordinates  = mesh, month.idx = seq_len(108)),
             options=list(verbose=TRUE))
t2 <-  Sys.time()
print(t2-t1)
summary(fit1)

lambda1 <- predict(
  fit1,
  ppxl_all,
  ~ data.frame(month.idx = month.idx, lambda = exp( Intercept+ fwi + mySmooth ))
)

cmp1.1 <- coordinates + month.idx ~   Intercept(1) +
        mySmooth(coordinates , model = spde) + 
        fwi(f.fwi.utm(.data.), model = "linear") + 
        TimeEffect(f.time(.data.), model='rw1')


t1 <-  Sys.time()
fit1.1 <- lgcp(cmp1.1, data=coords, samplers = domainSP,
            domain = list(coordinates  = mesh, month.idx = seq_len(108)),
            options=list(verbose=TRUE))
t2 <-  Sys.time()
print(t2-t1)
summary(fit1.1)

lambda1.1 <- predict(
  fit1.1,
  ppxl_all,
  ~ data.frame(month.idx = month.idx, lambda = exp( Intercept+ fwi + mySmooth + TimeEffect ))
)


cmp1.2 <- coordinates + month.idx ~   Intercept(1) +
  mySmooth(coordinates , model = spde) + 
  fwi(f.fwi.utm(.data.), model = 'linear') + 
  TimeEffect(f.time(.data.), model='rw1')


t1 <-  Sys.time()
fit1.2 <- lgcp(cmp1.2, data=coords, samplers = domainSP,
               domain = list(coordinates  = mesh, month.idx = seq_len(108)),
               options=list(verbose=TRUE))
t2 <-  Sys.time()
print(t2-t1)
summary(fit1.2)

lambda1.2 <- predict(
  fit1.2,
  ppxl_all,
  ~ data.frame(month.idx = month.idx, lambda = exp( fwi + Intercept))
)



cmp1.3 <- coordinates + month.idx ~   Intercept(1) +
  mySmooth(coordinates , model = spde) + 
  fwi(f.fwi.utm(.data.), model = spde1D ) + 
  TimeEffect(f.time(.data.), model='rw1')


t1 <-  Sys.time()
fit1.3 <- lgcp(cmp1.3, data=coords, samplers = domainSP,
               domain = list(coordinates  = mesh, month.idx = seq_len(108)),
               options=list(verbose=TRUE))
t2 <-  Sys.time()
print(t2-t1)
summary(fit1.3)

lambda1.3 <- predict(
  fit1.3,
  ppxl_all,
  ~ data.frame(month.idx = month.idx, lambda = exp( fwi + Intercept ))
)


cmp1.4 <- coordinates + month.idx ~   Intercept(1) +
  mySmooth(coordinates , model = spde) + 
  fwi(f.fwi.utm(.data.), model = spde1D ) + 
  TimeEffect(f.time(.data.), model='rw1')


t1 <-  Sys.time()
fit1.4 <- lgcp(cmp1.4, data=coords, samplers = domainSP,
               domain = list(coordinates  = mesh, month.idx = seq_len(108)),
               options=list(verbose=TRUE))
t2 <-  Sys.time()
print(t2-t1)
summary(fit1.4)

cmp2 <- coordinates + month.idx ~  Intercept(1) +
  mySmooth(coordinates , model = spde) +
  fwi(f.fwi.utm(.data.), model = "linear") + 
  rh(f.rh.utm(.data.), model='linear') + 
  temp(f.temp.utm(.data.), model='linear') +
  wspd(f.wspd.utm(.data.), model='linear') + 
  wdir(f.wdir.utm(.data.), model='linear') + 
  prcip(f.pricp.utm(.data.), model='linear')



t1 <-  Sys.time()
fit2 <- lgcp(cmp2, data=coords, samplers = domainSP,
             domain = list(coordinates  = mesh, month.idx = seq_len(108)),
             options=list(verbose=TRUE))
t2 <-  Sys.time()
print(t2-t1)
summary(fit2)



lambda2 <- predict(
  fit2,
  ppxl_all,
  ~ data.frame(month.idx = month.idx, lambda = exp(Intercept+ fwi + mySmooth +  rh + 
                                                     temp + wspd + wdir + prcip))
)


cmp2.1 <- coordinates + month.idx ~  Intercept(1) +
  mySmooth(coordinates , model = spde) +
  fwi(f.fwi.utm(.data.), model = "linear") + 
  rh(f.rh.utm(.data.), model='linear') + 
  temp(f.temp.utm(.data.), model='linear') +
  wspd(f.wspd.utm(.data.), model='linear') + 
  wdir(f.wdir.utm(.data.), model='linear') + 
  prcip(f.pricp.utm(.data.), model='linear') + 
  TimeEffect(f.time(.data.), model='rw1')



t1 <-  Sys.time()
fit2.1 <- lgcp(cmp2.1, data=coords, samplers = domainSP,
             domain = list(coordinates  = mesh, month.idx = seq_len(108)),
             options=list(verbose=TRUE))
t2 <-  Sys.time()
print(t2-t1)
summary(fit2.1)

lambda2.1 <- predict(
  fit2.1,
  ppxl_all,
  ~ data.frame(month.idx = month.idx, lambda = exp(Intercept+ fwi + mySmooth +  rh + 
                                                     temp + wspd + wdir + prcip+
                                                     TimeEffect))
)


cmp2.3 <- coordinates + month.idx ~  Intercept(1) +
  mySmooth(coordinates , model = spde) +
  fwi(f.fwi.utm(.data.), model = spde1D) + 
  # rh(f.rh.utm(.data.), model='linear') + 
  # temp(f.temp.utm(.data.), model='linear') +
  wspd(f.wspd.utm(.data.), model=spde1D.wspd) + 
  # wdir(f.wdir.utm(.data.), model='linear') + 
  # prcip(f.pricp.utm(.data.), model='linear') + 
  TimeEffect(f.time(.data.), model='rw1')

t1 <-  Sys.time()
fit2.3 <- lgcp(cmp2.3, data=coords, samplers = domainSP,
               domain = list(coordinates  = mesh, month.idx = seq_len(108)),
               options=list(verbose=TRUE))
t2 <-  Sys.time()
print(t2-t1)
summary(fit2.3)

cmp2.4 <- coordinates + month.idx ~  Intercept(1) +
  mySmooth(coordinates , model = spde) +
  fwi(f.fwi.utm(.data.), model = spde1D.fwi) + 
  rh(f.rh.utm(.data.), model=spde1D.rhumi) +
  temp(f.temp.utm(.data.), model=spde1D.temp) +
  wspd(f.wspd.utm(.data.), model=spde1D.wspd) + 
  # wdir(f.wdir.utm(.data.), model='linear') + 
  prcip(f.pricp.utm(.data.), model=spde1D.pricp) +
  TimeEffect(f.time(.data.), model='rw1')

t1 <-  Sys.time()
fit2.4 <- lgcp(cmp2.4, data=coords, samplers = domainSP,
               domain = list(coordinates  = mesh, month.idx = seq_len(108)),
               options=list(verbose=TRUE))
t2 <-  Sys.time()
print(t2-t1)
summary(fit2.4)


lambda2.1 <- predict(
  fit2.1,
  ppxl_all,
  ~ data.frame(month.idx = month.idx, lambda = exp(Intercept+ fwi + mySmooth +  rh + 
                                                     temp + wspd + wdir + prcip+
                                                     TimeEffect))
)



cmp3 <- coordinates + month.idx ~  Intercept(1) +
  mySmooth(coordinates , model = spde, group=(month.idx-1)%/%12+1,ngroup=9,
           control.group=list(model="ar1")) +
  fwi(f.fwi.utm(.data.), model = spde1D )

t1 <-  Sys.time()
fit3 <- lgcp(cmp3, data=coords, samplers = domainSP,
               domain = list(coordinates  = mesh, month.idx = seq_len(108)),
               options=list(verbose=TRUE))
t2 <-  Sys.time()
print(t2-t1)
summary(fit3)

lambda3 <- predict(
  fit3,
  ppxl_all,
  ~ data.frame(month.idx = month.idx, lambda = exp(Intercept+ mySmooth+fwi))
)


lprange <- range(lambda1.2$mean,lambda1.3$mean)
csc <- scale_fill_gradientn(colours = brewer.pal(9, "YlOrRd"), limits = lprange)


elapsed_months('2016/01/01','2012/01/01')
elapsed_months('2016/12/01','2012/01/01')

pl <- ggplot() +
  # gg(as(lambda1[lambda1$month.idx %in% 49:60,], "SpatialPixelsDataFrame"), aes(fill = mean)) +
  gg(loc.data.utm[loc.data.utm$month.idx  %in% 49:60, ], size = 0.1) +
  csc+
  facet_wrap(~month.idx) +
  coord_sf()
pl




pl1 <- ggplot() +
  gg(as(lambda1[lambda1$month.idx %in% 49:60,], "SpatialPixelsDataFrame"), aes(fill = mean)) +
  # gg(loc.data.utm[loc.data.utm$month.idx  %in% 49:60, ], size = 0.1) +
  csc+
  facet_wrap(~month.idx) +
  coord_sf()
pl1

pl1.1 <- ggplot() +
  gg(as(lambda1.1[lambda1.1$month.idx %in% 49:60,], "SpatialPixelsDataFrame"), aes(fill = mean)) +
  # gg(loc.data.utm[loc.data.utm$month.idx  %in% 49:60, ], size = 0.1) +
  csc+
  facet_wrap(~month.idx) +
  coord_sf()
pl1.1
multiplot(pl1,pl1.1,cols=2)

pl1.2 <- ggplot() +
  gg(as(lambda1.2[lambda1.2$month.idx %in% 49:60,], "SpatialPixelsDataFrame"), aes(fill = mean)) +
  # gg(loc.data.utm[loc.data.utm$month.idx  %in% 49:60, ], size = 0.1) +
  csc+
  facet_wrap(~month.idx) +
  coord_sf()
pl1.2

pl1.3 <- ggplot() +
  gg(as(lambda1.3[lambda1.3$month.idx %in% 49:60,], "SpatialPixelsDataFrame"), aes(fill = mean)) +
  # gg(loc.data.utm[loc.data.utm$month.idx  %in% 49:60, ], size = 0.1) +
  csc+
  facet_wrap(~month.idx) +
  coord_sf()
pl1.3
multiplot(pl1.2,pl1.3,cols=2)


pl2 <- ggplot() +
  gg(as(lambda2[lambda2$month.idx %in% 49:60,], "SpatialPixelsDataFrame"), aes(fill = mean)) +
  # gg(loc.data.utm, size = 0.1) +
  # csc+
  facet_wrap(~month.idx) +
  coord_sf()
pl2

pl2.1 <- ggplot() +
  gg(as(lambda2.1[lambda2.1$month.idx %in% 49:60,], "SpatialPixelsDataFrame"), aes(fill = mean)) +
  # gg(loc.data.utm, size = 0.1) +
  csc+
  facet_wrap(~month.idx) +
  coord_sf()
pl2
multiplot(pl2,pl2.1, cols=2)


pl3 <- ggplot() +
  gg(as(lambda3[lambda3$month.idx %in% 49:60,], "SpatialPixelsDataFrame"), aes(fill = mean)) +
  # gg(loc.data.utm, size = 0.1) +
  # csc+
  facet_wrap(~month.idx) +
  coord_sf()
pl3

prepare_residual_calculations <- function(samplers, domain, observations, time.idx) {
  observations <- observations[observations$month.idx==time.idx,]
  
  # Calculate the integration weights for A_integrate
  ips <- fm_int(domain = domain, samplers = samplers)
  
  # Set-up the A_integrate matrix
  # A_integrate has as many rows as polygons in the samplers,
  # as many columns as mesh points
  A_integrate <- inla.spde.make.A(
    mesh = domain, ips, weights = ips$weight,
    block = ips$.block, block.rescale = "none"
  )
  
  
  # Set-up the A_sum matrix
  # A_sum has as many rows as polygons in the samplers,
  # as many columns as observed points
  # each row has 1s for the points in the corresponding polygon
  idx <- sf::st_within(sf::st_as_sf(observations), sf::st_as_sf(samplers), sparse = TRUE)
  A_sum <- sparseMatrix(
    i = unlist(idx),
    j = rep(
      seq_len(nrow(observations)),
      vapply(idx, length, 1L)
    ),
    x = rep(1, length(unlist(idx))),
    dims = c(nrow(samplers), nrow(observations))
  )
  
  
  
  # Setting up the data frame for calculating residuals
  observations$obs <- TRUE
  df <- SpatialPointsDataFrame(
    coords = rbind(domain$loc[, 1:2], coordinates(observations)),
    data = bind_rows(data.frame(obs = rep(FALSE, domain$n), month.idx=time.idx), observations@data),
    proj4string = fm_CRS(domain)
  )
  
  # Return A-sum, A_integrate and the data frame for predicting the residuals
  list(A_sum = A_sum, A_integrate = A_integrate, df = df)
}



partition <- function(samplers, resolution = NULL, nrows = NULL, ncols = NULL) {
  # Create a grid for the given boundary
  if (is.null(resolution)) {
    grid <- rast(terra::ext(samplers),
                 crs = proj4string(samplers),
                 nrows = nrows, ncols = ncols
    )
  }
  
  if (is.null(c(nrows, ncols))) {
    grid <- rast(terra::ext(samplers),
                 crs = proj4string(samplers),
                 resolution = resolution
    )
  }
  
  gridPolygon <- terra::as.polygons(grid)
  
  # Extract the boundary with subpolygons only
  sf::as_Spatial(sf::st_as_sf(terra::intersect(gridPolygon, terra::vect(samplers))))
}



# 
residual_df <- function(model, df, expr, A_sum, A_integrate) {
  # Compute residuals
  res <- predict(
    object = model,
    newdata = df,
    ~ {
      lambda <- eval(expr)
      h1 <- lambda * 0 + 1
      h2 <- 1 / lambda
      h3 <- 1 / sqrt(lambda)
      data.frame(
        Estimate_Lambda = as.vector(A_integrate %*% (h1 * lambda)[!obs]),
        Scaling_Residuals =
          as.vector(A_sum %*% h1[obs]) -
          as.vector(A_integrate %*% (h1 * lambda)[!obs]),
        Inverse_Residuals =
          as.vector(A_sum %*% h2[obs]) -
          as.vector(A_integrate %*% (h2 * lambda)[!obs]),
        Pearson_Residuals =
          as.vector(A_sum %*% h3[obs]) -
          as.vector(A_integrate %*% (h3 * lambda)[!obs])
      )
    },
    used = bru_used(expr)
  )
  # Label the three types of residuals
  res$Estimate_Lambda$Type <- "Estimated Lambda"
  res$Scaling_Residuals$Type <- "Scaling Residuals"
  res$Inverse_Residuals$Type <- "Inverse Residuals"
  res$Pearson_Residuals$Type <- "Pearson Residuals"
  do.call(rbind, res)
}


B <- SpatialPolygonsDataFrame(domainSP, data.frame('weight'=1), match.ID = F) 

B1 <- partition(samplers = B, nrows = 20, ncols = 10)
plot(B1, main = "Grid partitions of B")




prepare_residual_calculations(samplers=B1,domain=mesh,observations=coords, time.idx=1)

residual_df_temproal <- function(time.idx, samplers, domain, observations, model, expr){
  As <- prepare_residual_calculations(
    samplers = samplers, domain = domain,
    observations = observations, time.idx=time.idx
  )
  
  res <- residual_df(
    model, As$df, expr,
    As$A_sum, As$A_integrate
  )
  res$month.idx <- time.idx
  return(res)
}


res1.3.1 <- do.call(rbind, lapply(49:60, residual_df_temproal, samplers=B1, domain=mesh, observations=coords,
                                  model=fit1.3,
                                  expr=expression(exp(fwi))))

res1 <- do.call(rbind, lapply(49:60, residual_df_temproal, samplers=B1, domain=mesh, observations=coords,
                                                        model=fit1,
                                                        expr=expression(exp(mySmooth + Intercept + fwi))))
res1.1 <- do.call(rbind, lapply(49:60, residual_df_temproal, samplers=B1, domain=mesh, observations=coords,
                              model=fit1.1,
                              expr=expression(exp(mySmooth + Intercept+ fwi + TimeEffect))))

res2 <- do.call(rbind, lapply(49:60, residual_df_temproal, samplers=B1, domain=mesh, observations=coords,
                              model=fit2,
                              expr=expression(exp(Intercept+ fwi + mySmooth +  rh + 
                                                    temp + wspd + wdir + prcip))))
res2.1 <- do.call(rbind, lapply(49:60, residual_df_temproal, samplers=B1, domain=mesh, observations=coords,
                              model=fit2.1,
                              expr=expression(exp(Intercept+ fwi + mySmooth +  rh + 
                                                    temp + wspd + wdir + prcip + 
                                                    TimeEffect))))


set_csc <- function(residuals, col_theme) {
  # Store data for the colour scale of the plots for each type of residual
  cscrange <- data.frame(
    residuals %>%
      group_by(Type) %>%
      summarise(maxabs = max(abs(mean)))
  )
  
  # Set the colour scale for all three types of residuals
  scaling_csc <-
    scale_fill_gradientn(
      colours = brewer.pal(9, col_theme[1]),
      name = "Scaling Residual",
      limits =
        cscrange[cscrange$Type == "Scaling Residuals", 2] *
        c(-1, 1)
    )
  
  inverse_csc <-
    scale_fill_gradientn(
      colours = brewer.pal(9, col_theme[2]),
      name = "Inverse Residual",
      limits =
        cscrange[cscrange$Type == "Inverse Residuals", 2] *
        c(-1, 1)
    )
  
  pearson_csc <-
    scale_fill_gradientn(
      colours = brewer.pal(9, col_theme[3]),
      name = "Pearson Residual",
      limits =
        cscrange[cscrange$Type == "Pearson Residuals", 2] *
        c(-1, 1)
    )
  
  list("Scaling" = scaling_csc, "Inverse" = inverse_csc, "Pearson" = pearson_csc)
}




grd <- SpatialGrid(GridTopology(cellcentre.offset = c(-10, 36),
                                cellsize = c(0.25,0.25), 
                                cells.dim = c(17, 29)))
# sgrd <- SpatialGridDataFrame(grd, data = data.frame(val = runif(240)), proj4string = CRS("+proj=longlat +datum=WGS84"))
grd_sp <- SpatialPixelsDataFrame(points = grd, data = data.frame(id = 1:length(grd)),proj4string = CRS("+proj=longlat +datum=WGS84"))

grd_poly <- as(grd_sp, 'SpatialPolygons')
grd_poly <- spTransform(grd_poly, CRS(projutm))


B2 <- as_Spatial(st_intersection(st_as_sf(grd_poly),st_as_sf(B)))
B2.Cent <- SpatialPointsDataFrame( gCentroid(B2,byid=TRUE), data=data.frame(weight=rep(1,length(B2))))

B2.Cent$month.idx <- 1



samplers.all <-  do.call(rbind, lapply(49:60, function(x) {
                                              B2@data$month.idx=x 
                                              return(B2)
                                              }))
# test1 <- residual_df_temproal(1:108, samplers=B1, domain=mesh, observations=coords,
#                               model=fit1.3,
#                               expr=expression(exp(Intercept+ fwi)))

res1.2.1 <- do.call(rbind, lapply(49:60, residual_df_temproal, samplers=B2, domain=mesh, observations=coords,
                                  model=fit1.2,
                                  expr=expression(exp(fwi))))
res1.3.1 <- do.call(rbind, lapply(49:60, residual_df_temproal, samplers=B2, domain=mesh, observations=coords,
                                  model=fit1.3,
                                  expr=expression(exp(fwi))))

res1.3.2 <- do.call(rbind, lapply(49:60, residual_df_temproal, samplers=B2, domain=mesh, observations=coords,
                                  model=fit1.3,
                                  expr=expression(exp(fwi+Intercept))))

res1.2 <- do.call(rbind, lapply(49:60, residual_df_temproal, samplers=B2, domain=mesh, observations=coords,
                                model=fit1.2,
                                expr=expression(exp(fwi+Intercept + mySmooth + TimeEffect))))
res1.3 <- do.call(rbind, lapply(49:60, residual_df_temproal, samplers=B2, domain=mesh, observations=coords,
                                  model=fit1.3,
                                  expr=expression(exp(fwi+Intercept + mySmooth + TimeEffect))))

res2.3 <- do.call(rbind, lapply(49:60, residual_df_temproal, samplers=B2, domain=mesh, observations=coords,
                                  model=fit2.3,
                                  expr=expression(exp(fwi+wspd+Intercept+mySmooth + TimeEffect))))

res2.4 <- do.call(rbind, lapply(49:60, residual_df_temproal, samplers=B2, domain=mesh, observations=coords,
                                model=fit2.4,
                                expr=expression(exp(fwi+rh+temp+wspd+prcip+Intercept+mySmooth + TimeEffect))))


samplers.all$Lambda <- res1.2.1 %>%
  filter(Type == "Estimated Lambda" ) %>%
  pull(mean)

samplers.all$Scaling_Residual <- res1.2 %>%
  filter(Type == "Scaling Residuals" ) %>%
  pull(mean)

# samplers.all$Lambda_l0.025 <- res1.3.2 %>%
#   filter(Type == "Estimated Lambda" ) %>%
#   pull(q0.025)
# samplers.all$Lambda_u0.975 <- res1.3.2 %>%
#   filter(Type == "Estimated Lambda" ) %>%
#   pull(q0.975)
samplers.all$FWI <- f.fwi.utm(do.call(rbind, lapply(49:60, function(x) {
                              B2.Cent$month.idx=x 
                                return(B2.Cent)
                              })))
samplers.all$Wind_Spd <- f.wspd.utm(do.call(rbind, lapply(49:60, function(x) {
  B2.Cent$month.idx=x 
  return(B2.Cent)
})))
samplers.all$RHumi <- f.rh.utm(do.call(rbind, lapply(49:60, function(x) {
  B2.Cent$month.idx=x 
  return(B2.Cent)
})))
samplers.all$Temp <- f.temp.utm(do.call(rbind, lapply(49:60, function(x) {
  B2.Cent$month.idx=x 
  return(B2.Cent)
})))
samplers.all$Pricp <- f.pricp.utm(do.call(rbind, lapply(49:60, function(x) {
  B2.Cent$month.idx=x 
  return(B2.Cent)
})))
samplers.all$Actual_Cnt <- 0

for (i in 49:60){
  points_sf <- st_as_sf(coords[coords$month.idx==i,])
  polygons_sf <- st_as_sf(samplers.all[samplers.all$month.idx==i,])

  
  samplers.all[samplers.all$month.idx==i,'Actual_Cnt'] <-  lengths(st_intersects(polygons_sf, points_sf))
}

ggplot() +
  geom_point() +
  # ggtitle(paste(Months[j], ", COX-FWI Relationship", sep='')) +
  # geom_line(aes(x = x, y = y), linetype = 1, data = tab.spde1) +
  # geom_ribbon(aes(x = x, y = y, ymin = ll95, ymax = ul95), alpha = 0.2, data = tab.spde1) +
  # geom_point(data=data.plot, aes(y=-2, x=temp), size=0.1 , col='black') +
  geom_point(data=samplers.all@data , aes(y=log(Lambda), x=FWI), size=0.1, col='blue') +
  # geom_point(data=samplers.all@data , aes(y=Actual_Cnt, x=FWI), size=0.1, col='red') +
  facet_wrap(. ~ month.idx, ncol=3) +
  ggtitle('Model fit1.2, log(lambda)- FWI ')


ggplot() +
  geom_point() +
  # ggtitle(paste(Months[j], ", COX-FWI Relationship", sep='')) +
  # geom_line(aes(x = x, y = y), linetype = 1, data = tab.spde1) +
  # geom_ribbon(aes(x = x, y = y, ymin = ll95, ymax = ul95), alpha = 0.2, data = tab.spde1) +
  # geom_point(data=data.plot, aes(y=-2, x=temp), size=0.1 , col='black') +
  geom_point(data=samplers.all@data , aes(y=FWI, x=Pricp), size=0.1, col='blue') +
  facet_wrap(. ~ month.idx, ncol=3)  +
  ggtitle("FWI-Precipitation Relationship")


my_residual_plot_temporal <- function(samplers.all, residuals, model_name) {
  samplers.all$Lambda <- residuals %>%
    filter(Type == "Estimated Lambda" ) %>%
    pull(mean)
  
  samplers.all$Scaling_Residual <- residuals %>%
    filter(Type == "Scaling Residuals" ) %>%
    pull(mean)
  
  scaling <- ggplot() +
    gg(samplers.all, aes(fill = Scaling_Residual), alpha = 1, colour = NA) +
    # csc["Scaling"] +
    theme(legend.position = "bottom") +
    facet_wrap(~month.idx) +
    labs(subtitle = paste(model_name, "Scaling"))

  origin_cnt <- ggplot() +
    gg(samplers.all, aes(fill = Actual_Cnt), alpha = 1, colour = NA) +
    # csc["Actual_Cnt"] +
    theme(legend.position = "bottom") +
    facet_wrap(~month.idx) +
    labs(subtitle = paste(model_name, "Actual Cnt"))
  list(
    Scaling = scaling, Origin_cnt=origin_cnt
  )
}


plot1.2 <- my_residual_plot_temporal(samplers.all, res1.2, model_name='fit1.2')
plot1.3 <- my_residual_plot_temporal(samplers.all, res1.3, model_name='fit1.3')
plot2.3 <- my_residual_plot_temporal(samplers.all, res2.3, model_name='fit2.3')
plot2.4 <- my_residual_plot_temporal(samplers.all, res2.4, model_name='fit2.4')
lprange.scale <- range(res1.3 %>%
                         filter(Type == "Scaling Residuals" ) %>%
                         pull(mean),
                       res2.3 %>%
                         filter(Type == "Scaling Residuals" ) %>%
                         pull(mean),
                       res2.4 %>%
                         filter(Type == "Scaling Residuals" ) %>%
                         pull(mean))
csc.scale <- scale_fill_gradientn(colours = brewer.pal(9, "YlOrRd"), limits = lprange.scale)

p1 <- plot1.3$Scaling+csc.scale
p2 <- plot2.3$Scaling+csc.scale
p3 <- plot2.4$Scaling+csc.scale

multiplot(p1,p2,p3, cols=3)

pc <- prcomp(samplers.all@data[,c('FWI','Wind_Spd','RHumi','Temp','Pricp')],
             center = TRUE,
             scale. = TRUE)
attributes(pc)




csc <- set_csc(res, rep("RdBu", 3))
model_name <- 'spatial_temporal'

 
residual_plot_temporal <- function(samplers.all, residuals, csc, model_name) {
  # Initialise the scaling residuals plot
  samplers.all$Residual <- residuals %>%
    filter(Type == "Scaling Residuals") %>%
    pull(mean)
  scaling <- ggplot() +
    gg(samplers.all, aes(fill = Residual), alpha = 1, colour = NA) +
    csc["Scaling"] +
    theme(legend.position = "bottom") +
    facet_wrap(~month.idx) +
    labs(subtitle = paste(model_name, "Scaling"))
  
  # Initialise the inverse residuals plot
  samplers.all$Residual <- residuals %>%
    filter(Type == "Inverse Residuals") %>%
    pull(mean)
  inverse <- ggplot() +
    gg(samplers.all, aes(fill = Residual), alpha = 1, colour = NA) +
    csc["Inverse"] +
    theme(legend.position = "bottom") +
    facet_wrap(~month.idx) +
    labs(subtitle = paste(model_name, "Inverse"))
  
  # Initialise the Pearson residuals plot
  samplers.all$Residual <- residuals %>%
    filter(Type == "Pearson Residuals") %>%
    pull(mean)
  pearson <- ggplot() +
    gg(samplers.all, aes(fill = Residual), alpha = 1, colour = NA) +
    csc["Pearson"] +
    theme(legend.position = "bottom") +
    facet_wrap(~month.idx) +
    labs(subtitle = paste(model_name, "Pearson"))
  
  # Return the three plots in a list
  list(
    Scaling = scaling, Inverse = inverse,
    Pearson = pearson
  )
  
}




fit_csc <- set_csc(res, rep("RdBu", 3))
# Store plots
plotB1_1 <- residual_plot_temporal(samplers.all, res1, fit_csc, "SPDE + FWI")
plotB1_1.1 <- residual_plot_temporal(samplers.all, res1.1, fit_csc, "SPDE + FWI + TimeIdx")

plotB1_2 <- residual_plot_temporal(samplers.all, res2, fit_csc, "SPDE + FWI + Weather")

plotB1_2.1 <- residual_plot_temporal(samplers.all, res2.1, fit_csc, "SPDE + FWI + Weather + TimeIdx")


((plotB1_1$Scaling | plotB1_1.1$Scaling )/
   ( plotB1_2$Scaling | plotB1_2.1$Scaling )) +
  plot_annotation(title = "Scaling Residuals") +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

((plotB1_1$Inverse | plotB1_1.1$Inverse )/
    ( plotB1_2$Inverse | plotB1_2.1$Inverse )) +
  plot_annotation(title = "Inverse Residuals") +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

((plotB1_1$Pearson | plotB1_1.1$Pearson )/
    ( plotB1_2$Pearson | plotB1_2.1$Pearson )) +
  plot_annotation(title = "Pearson Residuals") +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

res.summary <- cbind(
res1 %>% group_by(Type) %>% summarise(sum_res2= sum(mean^2)),
(res1.1 %>% group_by(Type) %>% summarise(sum_res2= sum(mean^2)))[,2],
(res2 %>% group_by(Type) %>% summarise(sum_res2= sum(mean^2)))[,2],
(res2.1 %>% group_by(Type) %>% summarise(sum_res2= sum(mean^2)))[,2])
colnames(res.summary) <- c('Type','Model1','Model1.1','Model2','Model2.1')
round(res.summary[,2:5],2)
# pred<- predict(
#   fit,
#   fm_pixels(mesh, mask = domainSP, format = "sp"),
#   ~ data.frame(
#     lambda = exp(mySmooth + Intercept),
#     loglambda = mySmooth + Intercept
#   )
# )
# 
# pl1 <- ggplot() +
#   gg(pred$lambda) +
#   gg(domainSP) +
#   ggtitle("LGCP fit to Points", subtitle = "(Response Scale)")
# 
# pl2 <- ggplot() +
#   gg(pred$loglambda) +
#   gg(domainSP, alpha = 0) +
#   ggtitle("LGCP fit to Points", subtitle = "(Linear Predictor Scale)")
# 
# multiplot(pl1, pl2, cols = 2)
# 
# 
# ggplot() +
#   gg(cbind(pred$lambda, data.frame(property = "q0.500")), aes(fill = median)) +
#   gg(cbind(pred$lambda, data.frame(property = "q0.025")), aes(fill = q0.025)) +
#   gg(cbind(pred$lambda, data.frame(property = "q0.975")), aes(fill = q0.975)) +
#   coord_equal() +
#   facet_wrap(~property)
# 
# 
int.plot <- plot(fit, "Intercept")
spde.range <- spde.posterior(fit, "mySmooth", what = "range")
spde.logvar <- spde.posterior(fit, "mySmooth", what = "log.variance")
range.plot <- plot(spde.range)
var.plot <- plot(spde.logvar)
# 
multiplot(range.plot, var.plot, int.plot)
# 
# 
# corplot <- plot(spde.posterior(fit, "mySmooth", what = "matern.correlation"))
# covplot <- plot(spde.posterior(fit, "mySmooth", what = "matern.covariance"))
# multiplot(covplot, corplot)
# 
# 
# Lambda <- predict(
#   fit,
#   fm_int(mesh, domainSP),
#   ~ sum(weight * exp(mySmooth + Intercept))
# )
# Lambda
# 
# Nest <- predict(
#   fit, fm_int(mesh, domainSP),
#   ~ data.frame(
#     N = 9000:10000,
#     dpois(9000:10000,
#           lambda = sum(weight * exp(mySmooth + Intercept))
#     )
#   )
# )
# 
# inla.qmarginal(c(0.025, 0.5, 0.975), marginal = list(x = Nest$N, y = Nest$mean))
# inla.emarginal(identity, marginal = list(x = Nest$N, y = Nest$mean))
# Nest$plugin_estimate <- dpois(Nest$N, lambda = Lambda$mean)
# ggplot(data = Nest) +
#   geom_line(aes(x = N, y = mean, colour = "Posterior")) +
#   geom_line(aes(x = N, y = plugin_estimate, colour = "Plugin"))
