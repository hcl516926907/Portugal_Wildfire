library(ncdf4)
library(CFtime)
library(lattice)
library(sp)
library(dplyr)
library(abind)
library(sf)
library(ggplot2)
library(cffdrs)
library(raster)
library(DescTools)
dir.data <- "/home/pgrad2/2448355h/My_PhD_Project/00_Dataset/Porgual_Wildfire"

####################################################################
#Load ERA5 Land covariates (2011–2023)
####################################################################

ncin <- nc_open(file.path(dir.data,'covariates', 'ERA5_Land_Daily_2011.nc'))

lon <- round(ncvar_get(ncin,"longitude"),1)
nlon <- dim(lon)
lat <- round(ncvar_get(ncin,"latitude"),1)
nlat <- dim(lat)
print(c(nlon,nlat))


time <- ncvar_get(ncin,"valid_time")
tunits <- ncatt_get(ncin,"valid_time","units")
tunits

cf <- CFtime(tunits$value, calendar = "proleptic_gregorian", time) # convert time to CFtime class
cf

timestamps <- CFtimestamp(cf) # get character-string times
time_cf <- CFparse(cf, timestamps) # parse the string into date components
time_cf

u10 <- ncvar_get(ncin,'u10')
v10 <- ncvar_get(ncin,'v10')
d2m <- ncvar_get(ncin,'d2m')
t2m <- ncvar_get(ncin,'t2m')
lai_hv <- ncvar_get(ncin,'lai_hv')
lai_lv <- ncvar_get(ncin,'lai_lv')

timestamp.all <- timestamps
time_cf.all <- time_cf

u10.all <- u10
v10.all <- v10
d2m.all <- d2m
t2m.all <- t2m 
lai_hv.all <- lai_hv
lai_lv.all <- lai_lv


for (year in 2012:2023){
  print(year)
  filename <- paste('ERA5_Land_Daily_', year, '.nc',sep='')
  ncin <- nc_open(file.path(dir.data,'covariates', filename))
  time <- ncvar_get(ncin,"valid_time")
  tunits <- ncatt_get(ncin,"valid_time","units")
  
  cf <- CFtime(tunits$value, calendar = "proleptic_gregorian", time) # convert time to CFtime class
  timestamps <- CFtimestamp(cf) # get character-string times
  time_cf <- CFparse(cf, timestamps) # parse the string into date components
  
  
  timestamp.all <- c(timestamp.all, timestamp)
  time_cf.all <- rbind(time_cf.all, time_cf)
  
  u10 <- ncvar_get(ncin,'u10')
  v10 <- ncvar_get(ncin,'v10')
  d2m <- ncvar_get(ncin,'d2m')
  t2m <- ncvar_get(ncin,'t2m')
  lai_hv <- ncvar_get(ncin,'lai_hv')
  lai_lv <- ncvar_get(ncin,'lai_lv')
  
  u10.all <- abind(u10.all, u10, along=3)
  v10.all <- abind(v10.all, v10, along=3)
  d2m.all <- abind(d2m.all, d2m, along=3)
  t2m.all <- abind(t2m.all, t2m, along=3)
  lai_hv.all <- abind(lai_hv.all, lai_hv, along=3)
  lai_lv.all <- abind(lai_lv.all, lai_lv, along=3)
  
}

u10.all[is.na(u10.all)] <- 0
v10.all[is.na(v10.all)] <- 0
d2m.all[is.na(d2m.all)] <- 0
t2m.all[is.na(t2m.all)] <- 0
lai_hv.all[is.na(lai_hv.all)] <- 0
lai_lv.all[is.na(lai_lv.all)] <- 0

dim(time_cf.all)


# Load vegetation data
ncin <- nc_open(file.path(dir.data,'covariates', 'low_veg_type.nc'))
lon.veg <- round(ncvar_get(ncin,"longitude"),1)
lat.veg <- round(ncvar_get(ncin,"latitude"),1)
summary(lon.veg)
summary(lat.veg)
type_lv.all <- ncvar_get(ncin,'tvl')
dim(type_lv.all)

type_lv.all <- round(type_lv.all[match(lon+360,lon.veg),match(lat,lat.veg)])


ncin <- nc_open(file.path(dir.data,'covariates', 'high_veg_type.nc'))
lon.veg <- round(ncvar_get(ncin,"longitude"),1)
lat.veg <- round(ncvar_get(ncin,"latitude"),1)
summary(lon.veg)
summary(lat.veg)
type_hv.all <- ncvar_get(ncin,'tvh')
dim(type_hv.all)

type_hv.all <- round(type_hv.all[match(lon+360,lon.veg),match(lat,lat.veg)])



ncin <- nc_open(file.path(dir.data,'covariates', 'low_veg_cover.nc'))
lon.veg <- round(ncvar_get(ncin,"longitude"),1)
lat.veg <- round(ncvar_get(ncin,"latitude"),1)
summary(lon.veg)
summary(lat.veg)
cover_lv.all <- ncvar_get(ncin,'cvl')
dim(cover_lv.all)

cover_lv.all <- cover_lv.all[match(lon+360,lon.veg),match(lat,lat.veg)]



ncin <- nc_open(file.path(dir.data,'covariates', 'high_veg_cover.nc'))
lon.veg <- round(ncvar_get(ncin,"longitude"),1)
lat.veg <- round(ncvar_get(ncin,"latitude"),1)
summary(lon.veg)
summary(lat.veg)
cover_hv.all <- ncvar_get(ncin,'cvh')
dim(cover_hv.all)

cover_hv.all <- cover_hv.all[match(lon+360,lon.veg),match(lat,lat.veg)]
dim(cover_hv.all)

# load daily total pricipation data

ncin <- nc_open(file.path(dir.data,'covariates', 'ERA5_Land_Daily_Total_Pricp_2011.nc'))
lon.pricp <- round(ncvar_get(ncin,"longitude"),1)
lat.pricp <- round(ncvar_get(ncin,"latitude"),1)


time <- ncvar_get(ncin,"valid_time")
tunits <- ncatt_get(ncin,"valid_time","units")

cf <- CFtime(tunits$value, calendar = "proleptic_gregorian", time) # convert time to CFtime class
timestamps <- CFtimestamp(cf) # get character-string times
time_cf_pricp <- CFparse(cf, timestamps) # parse the string into date components

time_cf_pricp$ind <- 0
time_cf_pricp[time_cf_pricp$hour==13,'ind'] <- 1
time_cf_pricp[1,'ind'] <- 1
time_cf_pricp$daily.index <- cumsum(time_cf_pricp$ind)

tp <- ncvar_get(ncin,'tp')
n.days <- nrow(unique(time_cf_pricp[,c('month','day')]))
tp.daily <- array(NA, dim=c(41,71,n.days))
for (i in 1:n.days){
  tp.daily[,,i] <- rowSums(tp[,,which(time_cf_pricp$daily.index==i)], dims=2)*10^3
  tp.daily[is.na(tp.daily)] <- 0
}
rm(tp)


tp.daily.all <- tp.daily
for (year in 2012:2023){
  print(year)
  filename <- paste('ERA5_Land_Daily_Total_Pricp_', year, '.nc',sep='')
  ncin <- nc_open(file.path(dir.data,'covariates', filename))
  time <- ncvar_get(ncin,"valid_time")
  tunits <- ncatt_get(ncin,"valid_time","units")
  
  cf <- CFtime(tunits$value, calendar = "proleptic_gregorian", time) # convert time to CFtime class
  timestamps <- CFtimestamp(cf) # get character-string times
  time_cf_pricp <- CFparse(cf, timestamps) # parse the string into date components
  
  
  time_cf_pricp$ind <- 0
  time_cf_pricp[time_cf_pricp$hour==13,'ind'] <- 1
  time_cf_pricp[1,'ind'] <- 1
  time_cf_pricp$daily.index <- cumsum(time_cf_pricp$ind)
  
  tp <- ncvar_get(ncin,'tp')
  n.days <- nrow(unique(time_cf_pricp[,c('month','day')]))
  tp.daily <- array(NA, dim=c(41,71,n.days))
  for (i in 1:n.days){
    tp.daily[,,i] <- rowSums(tp[,,which(time_cf_pricp$daily.index==i)], dims=2)*10^3
    tp.daily[is.na(tp.daily)] <- 0
  }
  rm(tp)
  
  
  tp.daily.all <- abind(tp.daily.all, tp.daily, along=3)
  
}

dim(tp.daily.all)


####################################################################
#Create relative humidity and fire weather index
####################################################################

relative_humi <- function(t,td){
  E <- exp(17.67*td/(td + 243.5))
  ES <- exp(17.67*t/(t + 243.5))
  return(100*E/ES)
}

fwi.daily.all <- array(NA, dim=dim(t2m.all))
ws.daily.all <- array(NA, dim=dim(t2m.all))
rh.daily.all <- array(NA, dim=dim(t2m.all))
for (i in 1:length(lon)){
  for (j in 1:length(lat)){
    print(c(i,j))
    temp <- t2m.all[i,j,] - 273.15
    dewpoint <- d2m.all[i,j,] - 273.15
    rh <- relative_humi(temp, dewpoint)
    ws <- sqrt(u10.all[i,j,]^2 + v10.all[i,j,]^2)
    
    rh.daily.all[i,j,] <- rh
    ws.daily.all[i,j,] <- ws
    
    prec <- tp.daily.all[i,j,]
    fwi.df <- data.frame(long=lon[i],lat=lat[j],
                         yr=time_cf.all$year,mon=time_cf.all$month,
                         day=time_cf.all$day, temp=temp, ws=ws,
                         rh = rh, prec=prec)
    fwi.daily.all[i,j,] <- round(fwi(
      fwi.df,
      init = data.frame(ffmc = 80, dmc = 10, dc = 16, lat = lat[j])
    )$FWI,2)
  }
}

#####create monthly covariates array###############
n.months <- nrow(unique(time_cf.all[,c('year','month')]))
time_cf.all$month.index <- time_cf.all$month + 12*(time_cf.all$year-2011) 
fwi.monthly.all <- array(NA, dim=c(41,71,n.months))
rh.monthly.all <- array(NA, dim=c(41,71,n.months))
ws.monthly.all <- array(NA, dim=c(41,71,n.months))
u10.monthly.all <- array(NA, dim=c(41,71,n.months))
v10.monthly.all <- array(NA, dim=c(41,71,n.months))
d2m.monthly.all <- array(NA, dim=c(41,71,n.months))
t2m.monthly.all <- array(NA, dim=c(41,71,n.months))
lai_hv.monthly.all <- array(NA, dim=c(41,71,n.months))
lai_lv.monthly.all <- array(NA, dim=c(41,71,n.months))
tp.monthly.all <- array(NA, dim=c(41,71,n.months))
fill.na <- function(x, i=5) {
  if( is.na(x)[i] ) {
    return( round(mean(x, na.rm=TRUE),2) )
  } else {
    return( round(x[i],2) )
  }
}

fill.na.mode <- function(x, i=5) {
  if( is.na(x)[i] ) {
    return( Mode(x, na.rm=TRUE)[1] )
  } else {
    return(x[i] )
  }
}

padding_boundary <- function(mat,method='mean'){

  r <- raster(mat)
  if (method=='mean'){
    r2 <- focal(r, w = matrix(1,3,3), fun = fill.na, 
                pad = TRUE, na.rm = FALSE )
  }else{
    r2 <- focal(r, w = matrix(1,3,3), fun = fill.na.mode, 
                pad = TRUE, na.rm = FALSE )
  }

  return(as.matrix(r2))
}

for (i in 1:n.months){
  fwi.monthly.all[,,i] <- rowMeans(fwi.daily.all[,,which(time_cf.all$month.index==i)], dims=2)
  rh.monthly.all[,,i] <- rowMeans(rh.daily.all[,,which(time_cf.all$month.index==i)], dims=2)
  ws.monthly.all[,,i] <- rowMeans(ws.daily.all[,,which(time_cf.all$month.index==i)], dims=2)
  u10.monthly.all[,,i] <- rowMeans(u10.all[,,which(time_cf.all$month.index==i)], dims=2)
  v10.monthly.all[,,i] <- rowMeans(v10.all[,,which(time_cf.all$month.index==i)], dims=2)
  lai_hv.monthly.all[,,i] <- rowMeans(lai_hv.all[,,which(time_cf.all$month.index==i)], dims=2)
  lai_lv.monthly.all[,,i] <- rowMeans(lai_lv.all[,,which(time_cf.all$month.index==i)], dims=2)
  t2m.monthly.all[,,i] <- rowMeans(t2m.all[,,which(time_cf.all$month.index==i)], dims=2) - 273.15
  d2m.monthly.all[,,i] <- rowMeans(d2m.all[,,which(time_cf.all$month.index==i)], dims=2) - 273.15
  tp.monthly.all[,,i] <- rowMeans(tp.daily.all[,,which(time_cf.all$month.index==i)], dims=2)

  # Below is to padding the missing value around the boundary of Portugal.
  na.index <- abs(t2m.monthly.all[,,i]+273.15)<10^-6
  
  fwi.monthly.all[,,i][na.index] <- NA
  fwi.monthly.all[,,i] <- padding_boundary(fwi.monthly.all[,,i])
  
  rh.monthly.all[,,i][na.index] <- NA
  rh.monthly.all[,,i] <- padding_boundary(rh.monthly.all[,,i])
  
  ws.monthly.all[,,i][na.index] <- NA
  ws.monthly.all[,,i] <- padding_boundary(ws.monthly.all[,,i])
  
  u10.monthly.all[,,i][na.index] <- NA
  u10.monthly.all[,,i] <- padding_boundary(u10.monthly.all[,,i])
  
  v10.monthly.all[,,i][na.index] <- NA
  v10.monthly.all[,,i] <- padding_boundary(v10.monthly.all[,,i])
  
  lai_hv.monthly.all[,,i][na.index] <- NA
  lai_hv.monthly.all[,,i] <- padding_boundary(lai_hv.monthly.all[,,i])
  
  lai_lv.monthly.all[,,i][na.index] <- NA
  lai_lv.monthly.all[,,i] <- padding_boundary(lai_lv.monthly.all[,,i])
  
  t2m.monthly.all[,,i][na.index] <- NA
  t2m.monthly.all[,,i] <- padding_boundary(t2m.monthly.all[,,i])
  
  d2m.monthly.all[,,i][na.index] <- NA
  d2m.monthly.all[,,i] <- padding_boundary(d2m.monthly.all[,,i])
  
  tp.monthly.all[,,i][na.index] <- NA
  tp.monthly.all[,,i] <- padding_boundary(d2m.monthly.all[,,i])
  
}



type_lv.all[na.index] <- NA
type_lv.all <- padding_boundary(type_lv.all,method='mode')
type_hv.all[na.index] <- NA
type_hv.all <- padding_boundary(type_hv.all,method='mode')

cover_lv.all[na.index] <- NA
cover_lv.all <- padding_boundary(cover_lv.all)
cover_hv.all[na.index] <- NA
cover_hv.all <- padding_boundary(cover_hv.all)




#################################################################
#Value Check
################################################################



grid_era5 <- expand.grid(lon=lon, lat=lat)
library(rnaturalearth)
world <- ne_countries(scale = "medium", returnclass = "sf")
portugal <- world[world$admin == "Portugal", ]

grid_era5$data <- as.vector(rep(1,nrow(grid_era5)))
grid_sf <- st_as_sf(grid_era5, coords = c("lon", "lat"), crs = 4326)

mainland_bbox <- st_as_sfc(st_bbox(c(xmin = -10, xmax = -6, ymin = 35 , ymax = 43), crs = st_crs(portugal)))
mainland <- st_intersection(portugal, mainland_bbox)
# Clip grid to the mainland of Portugal (conceptually - ensure mainland shapefile for real use)
grid_clipped <- st_intersection(grid_sf, mainland)
grid_within_portugal <- st_coordinates(grid_clipped)
grid_within_portugal <- round(grid_within_portugal,2)


ggplot() + geom_sf(data = mainland) +
  geom_sf(data = grid_clipped)
 

#data check
cov.name <- 'tp'
t <- 1
spdf <- SpatialPixelsDataFrame(points = grid_era5[,1:2],
                               data = data.frame(var=as.vector(get(paste(cov.name,'.monthly.all',sep=''))[,,t]),
                                                 time=t))
proj4string(spdf) <- CRS("EPSG:4326")
ggplot() +
  geom_sf(data = mainland) +
  geom_sf(data = st_intersection(st_as_sf(spdf),mainland),aes(color=var))




cov.name <- 'type_hv'
spdf <- SpatialPixelsDataFrame(points = grid_era5[,1:2],
                               data = data.frame(var=as.vector(get(paste(cov.name,'.all',sep='')))))
proj4string(spdf) <- CRS("EPSG:4326")
ggplot() +
  geom_sf(data = mainland) +
  geom_sf(data = st_intersection(st_as_sf(spdf),mainland),aes(color=var))


save(type_lv.all,
     type_hv.all,
     cover_lv.all,
     cover_hv.all,
     fwi.monthly.all,
     rh.monthly.all,
     ws.monthly.all,
     u10.monthly.all,
     v10.monthly.all,
     lai_hv.monthly.all,
     lai_lv.monthly.all,
     t2m.monthly.all,
     d2m.monthly.all,
     tp.monthly.all,
     file=file.path(dir.data, 'Covariates.RData')
     )


####################################################################
# Merge Covariates with fire data
####################################################################

load(file.path(dir.data, 'Covariates.RData'))
load(file.path(dir.data,"Wildfire.RData"))


wildfire <- as.data.frame(data[data$AreaTotal_ha>1&data$length>3 ,])
# wildfire <- as.data.frame(data[data$AreaTotal_ha>1 & data$length>6  ,])
# wildfire <- as.data.frame(data[data$AreaPov_ha>1  ,])
dim(wildfire)
# Area>1, length>6 5010 rows
# Area>1, length>12 1787 rows
# Area>1, length>24 788 rows

wildfire$year <- as.integer(format(wildfire$open,format='%Y'))
wildfire$month <- as.integer(format(wildfire$open,format='%m'))
wildfire$time.idx <- (wildfire$year-min(wildfire$year))*12 + wildfire$month

wildfire%>%group_by(year)%>%summarise(burn_area=sum(AreaTotal_ha),
                                      n=n())


grd <- st_intersection(st_as_sf(grid_era5,coords=c('lon','lat'),crs=4326),mainland) 

grd <- SpatialPixelsDataFrame(points = st_coordinates(grd), data = data.frame(grid.idx = 1:nrow(grd)),proj4string = CRS("+proj=longlat +datum=WGS84"))

B1 <- as(grd, 'SpatialPolygonsDataFrame')
B1$lon.grid <- grd@coords[,1]
B1$lat.grid <- grd@coords[,2]


data.merge <- B1 |>  
  st_as_sf() |> # cast to sf
  st_join(st_as_sf(wildfire,coords=c('lon','lat'),crs=4326)) 

data.agg <- data.merge %>% st_drop_geometry() %>% 
  group_by(grid.idx, time.idx) %>%
  summarise(area_ha = sum(AreaTotal_ha), 
            log_ba = log(area_ha),
            y = n())

data.fit <- do.call(rbind, lapply(1:n.months, function(x) {B1@data$time.idx = x
return(B1@data)}))
print(dim(data.fit))
data.fit <- merge(data.fit, data.agg[,c('grid.idx','time.idx','y','area_ha','log_ba')],
                   by=c('grid.idx','time.idx'),all.x=T)

print(dim(data.fit))


data.fit[is.na(data.fit$y),'y'] <- 0
table(data.fit$y)


data.fit$month <- (data.fit$time.idx-1)%%12 + 1
data.fit$year <- (data.fit$time.idx-1)%/%12 + min(wildfire$year)

data.fit <- data.fit[order(data.fit$time.idx,data.fit$grid.idx),]

ggplot(data=st_as_sf(data.fit[data.fit$time.idx==152,], coords = c('lon.grid','lat.grid'),crs=4326))+
  geom_sf(aes(color=log_ba))


f.get.cov <- function(dataset, cov.name){
  time <- dataset$time.idx

  x <- dataset$lon.grid
  y <- dataset$lat.grid
  spp <- SpatialPoints(data.frame(x=x,y=y),
                       proj4string=CRS("EPSG:4326" ))
  v <- rep(NA, nrow(dataset))
  for (t in unique(time)){
    if(cov.name %in% c('type_hv.all','type_lv.all',
                       'cover_hv.all','cover_lv.all')){
      spdf <- SpatialPixelsDataFrame(points = grid_era5[,1:2],
                                     data = data.frame(var=as.vector(get(cov.name)),
                                                       time=t))
    }else{
      spdf <- SpatialPixelsDataFrame(points = grid_era5[,1:2],
                                     data = data.frame(var=as.vector(get(cov.name)[,,t]),
                                                       time=t))
    }

    proj4string(spdf) <- CRS("EPSG:4326")
    idx <- which(time==t)
    v[idx] <- over(spp[idx,],spdf[,'var'])$var
  }
  return(v)
}

data.fit$FWI <- f.get.cov(data.fit,'fwi.monthly.all')
data.fit$HVegCov <- f.get.cov(data.fit,'cover_hv.all')
data.fit$HVegLAI <- f.get.cov(data.fit, 'lai_hv.monthly.all')
data.fit$HVegTyp <- as.factor(f.get.cov(data.fit, 'type_hv.all'))
data.fit$LVegCov <- f.get.cov(data.fit,'cover_lv.all')
data.fit$LVegLAI <- f.get.cov(data.fit, 'lai_lv.monthly.all')
data.fit$LVegTyp <- as.factor(f.get.cov(data.fit, 'type_lv.all'))
data.fit$Pricp <- f.get.cov(data.fit, 'tp.monthly.all')
data.fit$RHumi <- f.get.cov(data.fit, 'rh.monthly.all')
data.fit$Temp <- f.get.cov(data.fit, 't2m.monthly.all')
data.fit$UComp <- f.get.cov(data.fit, 'u10.monthly.all')
data.fit$VComp <- f.get.cov(data.fit, 'v10.monthly.all')
data.fit$DewPoint <- f.get.cov(data.fit, 'd2m.monthly.all')
data.fit$WindSpeed <- f.get.cov(data.fit, 'ws.monthly.all')



ggplot(data=st_as_sf(data.fit[data.fit$time.idx==31,], coords = c('lon.grid','lat.grid'),crs=4326))+
  geom_sf(aes(color=Temp))

save(data.fit,file=file.path(dir.data,'Data_For_Fitting.RData'))




####################################################################
# Aggregate the covariates to council level
####################################################################

load(file.path(dir.data,"Wildfire.RData"))

dist<-shapefile(file.path(dir.data,'shapefile', "distritos.shp"))
# Municipalities
conc<-shapefile(file.path(dir.data, "shapefile","concelhos.shp"))

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


wildfire <- as.data.frame(data[data$AreaTotal_ha>1&data$length>3 ,])
# wildfire <- as.data.frame(data[data$AreaTotal_ha>1 & data$length>6  ,])
# wildfire <- as.data.frame(data[data$AreaPov_ha>1  ,])
dim(wildfire)
# Area>1, length>6 5010 rows
# Area>1, length>12 1787 rows
# Area>1, length>24 788 rows

wildfire$year <- as.integer(format(wildfire$open,format='%Y'))
wildfire$month <- as.integer(format(wildfire$open,format='%m'))
wildfire$time.idx <- (wildfire$year-min(wildfire$year))*12 + wildfire$month

wildfire%>%group_by(year)%>%summarise(burn_area=sum(AreaTotal_ha),
                                      n=n())

map <- sf_conc[,'NAME_2']
rownames(map) <- sf_conc$NAME_2

library(spdep)
nb <- poly2nb(map)
map$grid.idx <- 1:nrow(map)
centroids <- st_centroid(map)
centroid_coords <- st_coordinates(centroids)
map$lon.grid <- centroid_coords[,1]
map$lat.grid <- centroid_coords[,2]
ggplot(map) + geom_sf()+
  geom_text(aes(lon.grid, lat.grid, label = grid.idx), size=2,color='red') 


data.merge <- map |>
  st_join(st_as_sf(wildfire,coords=c('lon','lat'),crs=4326)) 

data.agg <- data.merge %>% st_drop_geometry() %>% 
  group_by(grid.idx, time.idx) %>%
  summarise(area_ha = sum(AreaTotal_ha), 
            log_ba = log(area_ha),
            y = n())

B1 <- st_drop_geometry(map)

n.months <- 156

data.fit.council <- do.call(rbind, lapply(1:n.months, function(x) {B1$time.idx = x
return(B1)}))
print(dim(data.fit.council))
data.fit.council <- merge(data.fit.council, data.agg[,c('grid.idx','time.idx','y','area_ha','log_ba')],
                  by=c('grid.idx','time.idx'),all.x=T)

print(dim(data.fit.council))


data.fit.council[is.na(data.fit.council$y),'y'] <- 0
table(data.fit.council$y)


data.fit.council$month <- (data.fit.council$time.idx-1)%%12 + 1
data.fit.council$year <- (data.fit.council$time.idx-1)%/%12 + min(wildfire$year)

data.fit.council <- data.fit.council[order(data.fit.council$time.idx,data.fit.council$grid.idx),]

ggplot(data=st_as_sf(data.fit.council[data.fit.council$time.idx==152,], coords = c('lon.grid','lat.grid'),crs=4326))+
  geom_sf(aes(color=log_ba))


load(file=file.path(dir.data,'Data_For_Fitting.RData'))
data.fit.sf <- data.fit |> 
  st_as_sf(coords = c("lon.grid", "lat.grid"), crs = 4326) |>
  st_join(sf_conc[,'NAME_2'], join = st_nearest_feature) |>
  st_drop_geometry()
           
my.mode <- function(x){
  return(as.integer(tail(names(sort(table(x))), 1)))
}

cov.council <- data.fit.sf %>%
  group_by(NAME_2,time.idx) %>%
  summarise(
    FWI = mean(FWI),
    HVegCov = mean(HVegCov),
    HVegLAI = mean(HVegLAI),
    HVegTyp = my.mode(HVegTyp),
    LVegCov = mean(LVegCov),
    LVegLAI = mean(LVegLAI),
    LVegTyp = my.mode(LVegTyp),
    Pricp = mean(Pricp),
    RHumi = mean(RHumi),
    Temp = mean(Temp),
    UComp = mean(UComp),
    VComp = mean(VComp),
    DewPoint = mean(DewPoint),
    WindSpeed = mean(WindSpeed)
  )


data.fit.council <- merge(data.fit.council,cov.council,by=c("NAME_2",'time.idx'),all.x=TRUE)

cov.names <- c('FWI','HVegCov','HVegLAI','HVegTyp','LVegCov','LVegLAI','LVegTyp',
               'Pricp','RHumi','Temp','UComp','VComp','DewPoint','WindSpeed')

missing.council <- unique(data.fit.council[is.na(data.fit.council$FWI),'NAME_2'] )

if(length(missing.council)>0){
  for (t in 1:n.months){
    print(t)
    for (council in missing.council){
      cond1 <- data.fit.council$time.idx==t & data.fit.council$NAME_2==council
      cond2 <- data.fit$time.idx==t 
      tmp <- st_as_sf(data.fit[cond2, c("lon.grid", "lat.grid",cov.names)],coords = c("lon.grid", "lat.grid"), crs = 4326)
      tmp$HVegTyp <- as.integer(as.character(tmp$HVegTyp))
      tmp$LVegTyp <- as.integer(as.character(tmp$LVegTyp))
      
      data.fit.council[cond1,cov.names] <- 
        data.fit.council[cond1,c('lon.grid','lat.grid')] |>
        st_as_sf(coords = c("lon.grid", "lat.grid"), crs = 4326) |>
        st_join(tmp, join = st_nearest_feature) |>
        st_drop_geometry()
    }
  }
}

  
ggplot(data=st_as_sf(data.fit.council[data.fit.council$time.idx==31,], coords = c('lon.grid','lat.grid'),crs=4326))+
  geom_sf(aes(color=Temp))



data.fit.council <- data.fit.council |> st_as_sf(coords=c('lon.grid','lat.grid'),crs=4326,remove=FALSE) |>
  st_join(sf_districts[,'NAME_1']) |> st_drop_geometry()



data.fit.council <- data.fit.council[order(data.fit.council$time.idx,data.fit.council$grid.idx),]


save(data.fit.council,file=file.path(dir.data,'Data_For_Fitting_Council.RData'))

