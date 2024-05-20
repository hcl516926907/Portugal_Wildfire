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
load(file.path(dir.data, 'burn_area', 'Vegetation.RData'))


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
data.new$month <- as.integer(format(data.new$open,format='%m'))
data.new$month.idx <- elapsed_months(data.new$open, date.start)

# time_cf$date <- paste(time_cf$year,time_cf$month,time_cf$day,sep='/')
# time_cf$month.idx <- elapsed_months(time_cf$date, date.start)

loc.data.utm <- st_as_sf(data.new, coords=c('x_utm_new','y_utm_new'), crs='+proj=utm +zone=29 +datum=WGS84 +units=km +no_defs' )


coords <- SpatialPointsDataFrame(data.new[data.new$month==9,], coords=data.new[data.new$month==9,c('x_utm_new','y_utm_new')], 
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

inner.length <- 20
outer.length <- 50
cutoff <- 6

# inner.length <- 40
# outer.length <- 70
# cutoff <- 15

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


f.time <- function(where) {
  v <- where$time.idx
  return(v)
}

f.time.factor <- function(where) {
  v <- as.factor(where$time.idx)
  return(v)
}


coords$time.idx <- coords$year.idx
f.NDVI.utm <- function(where) {
  x <- where@coords[,1]
  y <- where@coords[,2]
  time <- where$time.idx
  spp <- SpatialPoints(data.frame(x=x,y=y),
                       proj4string=CRS('+proj=utm +zone=29 +datum=WGS84 +units=km +no_defs' ))
  v <- rep(NA, nrow(where))
  for (t in unique(time)){
    temp_spdf <- NDVI.list[[2011+t]]
    inside <- !is.na(over(temp_spdf, domainSP))
    
    # Subset SpatialPixelsDataFrame to only include points inside the polygon
    temp_spdf <- temp_spdf[inside, ]
    temp_spdf$CMG.0.05.Deg.Monthly.NDVI <- (temp_spdf$CMG.0.05.Deg.Monthly.NDVI- mean(temp_spdf$CMG.0.05.Deg.Monthly.NDVI))
    temp_spdf$CMG.0.05.Deg.Monthly.NDVI <- temp_spdf$CMG.0.05.Deg.Monthly.NDVI/max(abs(temp_spdf$CMG.0.05.Deg.Monthly.NDVI))
    idx <- which(time==t)
    v[idx] <- over(spp[idx,],temp_spdf[,'CMG.0.05.Deg.Monthly.NDVI'])$CMG.0.05.Deg.Monthly.NDVI
    if (any(is.na(v[idx]))) {
      v[idx] <- bru_fill_missing(temp_spdf, spp[idx,], v[idx])
    }
  }
  return(v)
}

summary(f.NDVI.utm(coords))

f.NDVI.utm.discrete <- function(where) {
  x <- where@coords[,1]
  y <- where@coords[,2]
  time <- where$time.idx
  spp <- SpatialPoints(data.frame(x=x,y=y),
                       proj4string=CRS('+proj=utm +zone=29 +datum=WGS84 +units=km +no_defs' ))
  v <- rep(NA, nrow(where))
  for (t in unique(time)){
    temp_spdf <- NDVI.list[[2011+t]]
    inside <- !is.na(over(temp_spdf, domainSP))
    # Subset SpatialPixelsDataFrame to only include points inside the polygon
    temp_spdf <- temp_spdf[inside, ]
    temp_spdf$CMG.0.05.Deg.Monthly.NDVI <- (temp_spdf$CMG.0.05.Deg.Monthly.NDVI- mean(temp_spdf$CMG.0.05.Deg.Monthly.NDVI))
    temp_spdf$CMG.0.05.Deg.Monthly.NDVI <- temp_spdf$CMG.0.05.Deg.Monthly.NDVI/max(abs(temp_spdf$CMG.0.05.Deg.Monthly.NDVI))
    
    idx <- which(time==t)
    v[idx] <- over(spp[idx,],temp_spdf[,'CMG.0.05.Deg.Monthly.NDVI'])$CMG.0.05.Deg.Monthly.NDVI
    if (any(is.na(v[idx]))){
      v[idx] <- bru_fill_missing(temp_spdf, spp[idx,], v[idx])
    }
    breaks <- seq(-1, 1, length.out = 21)  # Creates 20 intervals
    # Use cut to assign each value to a bin
    bins <- cut(v[idx], breaks, include.lowest = TRUE)
    v[idx] <- as.numeric(bins)
  }
  return(v)
}
table(f.NDVI.utm.discrete(coords))
summary(f.NDVI.utm(coords))


f.EVI.utm <- function(where) {
  x <- where@coords[,1]
  y <- where@coords[,2]
  time <- where$time.idx
  spp <- SpatialPoints(data.frame(x=x,y=y),
                       proj4string=CRS('+proj=utm +zone=29 +datum=WGS84 +units=km +no_defs' ))
  v <- rep(NA, nrow(where))
  for (t in unique(time)){
    temp_spdf <- EVI.list[[2011+t]]
    inside <- !is.na(over(temp_spdf, domainSP))
    # Subset SpatialPixelsDataFrame to only include points inside the polygon
    temp_spdf <- temp_spdf[inside, ]
    temp_spdf$CMG.0.05.Deg.Monthly.EVI <- (temp_spdf$CMG.0.05.Deg.Monthly.EVI- mean(temp_spdf$CMG.0.05.Deg.Monthly.EVI))
    temp_spdf$CMG.0.05.Deg.Monthly.EVI <- temp_spdf$CMG.0.05.Deg.Monthly.EVI/max(abs(temp_spdf$CMG.0.05.Deg.Monthly.EVI))
    idx <- which(time==t)
    v[idx] <- over(spp[idx,],temp_spdf[,'CMG.0.05.Deg.Monthly.EVI'])$CMG.0.05.Deg.Monthly.EVI
    if (any(is.na(v[idx]))) {
      v[idx] <- bru_fill_missing(temp_spdf, spp[idx,], v[idx])
    }
  }
  return(v)
}

summary(f.EVI.utm(coords))


mesh1D.NDVI <- fm_mesh_1d(seq(-1.2,1.2,by=0.05), boundary = "free")
spde1D.NDVI <- inla.spde2.pcmatern(mesh1D.NDVI,
                                   prior.range = c(0.005, 0.1),
                                   prior.sigma = c(0.1, 0.1)
)


cmp2.3 <- coordinates + time.idx ~  Intercept(1) +
  CSmooth(coordinates , model = spde) + 
  TSmooth(coordinates , model = spde, group = time.idx, ngroup = 9) +
  NDVI(f.NDVI.utm(.data.), model = spde1D.NDVI)+
  EVI(f.EVI.utm(.data.), model = spde1D.NDVI)



t1 <-  Sys.time()
fit2.3 <- lgcp(cmp2.3, data=coords, samplers = domainSP,
               domain = list(coordinates  = mesh, time.idx = seq_len(9)),
               options=list(verbose=TRUE))
t2 <-  Sys.time()
print(t2-t1)
summary(fit2.3)

save(fit2.3, file=file.path(dir.out, 'LGCP_2SPDE_wi_cova.RData'))
