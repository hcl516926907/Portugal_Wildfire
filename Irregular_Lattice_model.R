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
library(spdep)
library(raster)
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
data.new$month.idx <- as.integer(format(data.new$open,format='%m'))
data.new$time.idx <- elapsed_months(data.new$open, date.start)


# time_cf$date <- paste(time_cf$year,time_cf$month,time_cf$day,sep='/')
# time_cf$month.idx <- elapsed_months(time_cf$date, date.start)

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

# inner.length <- 20
# outer.length <- 50
# cutoff <- 6
inner.length <- 40
outer.length <- 70
cutoff <- 15
# 
mesh <- inla.mesh.2d(loc.domain = loc.d, max.edge=c(inner.length,outer.length), cutoff=cutoff,
                     crs=projutm)
mesh$n
plot(mesh)
# points(coords, col = 4, pch = 19, cex=0.2)



domain.polys <- Polygons(list(Polygon(loc.d)), '0')
domainSP <- SpatialPolygons(list(domain.polys),proj4string=CRS(projutm))

# #no of mesh points...
# nv <- mesh$n 
# 
# ggplot() +
#   gg(mesh) +
#   gg(domainSP) +
#   # gg(mrsea$samplers) +
#   gg(loc.data.utm[cond,], size = 0.1) +
#   facet_wrap(~year.idx) +
#   ggtitle("Fire Occurance by year")


# 
r <- 50
p1 <- 0.1
sigma <- 8
p2 <- 0.1
# 
spde <- inla.spde2.pcmatern(mesh = mesh,
                            # PC-prior on range: P(practic.range < r) = p1
                            prior.range = c(r, p1),
                            # PC-prior on sigma: P(sigma > sigma) = p2
                            prior.sigma = c(sigma, p2))


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

# 
grd <- SpatialGrid(GridTopology(cellcentre.offset = c(-10, 36),
                                cellsize = c(0.25,0.25),
                                cells.dim = c(17, 29)))

# grd <- SpatialGrid(GridTopology(cellcentre.offset = c(-10, 36),
#                                 cellsize = c(0.125,0.125),
#                                 cells.dim = c(34, 58)))

# sgrd <- SpatialGridDataFrame(grd, data = data.frame(val = runif(240)), proj4string = CRS("+proj=longlat +datum=WGS84"))

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

ggplot(st_as_sf(B2)) + geom_sf(aes(fill=E))

data.merge <- B2 |>  
  st_as_sf() |> # cast to sf
  mutate(grid_id = row_number()) |> # create unique ID
  st_join(loc.data.utm) |> # join the species dataset
  group_by(grid_id)
#######################################################################
#Fewer points after join. Boundary issues.
####################################################################

data.fit <- data.merge %>% st_drop_geometry() %>% filter( time.idx<=108, length >=24*60 )%>% 
  group_by(grid.idx,grid_id, year.idx, month.idx, time.idx) %>% 
  summarise(y = n(), 
            lon.grid = mean(lon.grid),
            lat.grid = mean(lat.grid),
            x.utm = mean(grid.cent.x.utm),
            y.utm = mean(grid.cent.y.utm),
            E = mean(E))

data.fit2 <- do.call(rbind, lapply(1:108, function(x) {B2@data$time.idx = x 
return(B2@data)}))
print(dim(data.fit2))
data.fit2 <- merge(data.fit2, data.fit[,c('grid.idx','time.idx','y')],
                   by=c('grid.idx','time.idx'),all.x=T)
# data.fit2 <- merge(data.fit2, data.fit[,c('grid.idx','month.idx','y')],
#                    by.x=c('grid.idx','time.idx'),
#                    by.y=c('grid.idx','month.idx'),all.x=T)

print(dim(data.fit2))
data.fit2[is.na(data.fit2$y),'y'] <- 0
summary(data.fit2$y)




data.fit3 <- reshape(data.fit2[,c('grid.idx','time.idx','y')],
                     timevar = "time.idx",
                     idvar = "grid.idx",
                     direction = "wide"
)

B2.merge <- merge(B2, data.fit3, by = "grid.idx")



B2.adj <- poly2nb(B2)
nb2INLA("map.adj", B2.adj)
g <- inla.read.graph(filename = "map.adj")

data.fit2$idarea <- as.numeric(as.factor(data.fit2$grid.idx))
data.fit2$idarea1 <- data.fit2$idarea

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

data.fit2$fwi <- f.get.cov(data.fit2,'fwi')
data.fit2$HVegCov <- f.get.cov(data.fit2,'HVegCov')
data.fit2$HVegHAI <- f.get.cov(data.fit2, 'HVegLAI')
# data.fit2$HVegTyp <- as.factor(f.get.cov(data.fit2, 'HVegTyp'))
data.fit2$LVegCov <- f.get.cov(data.fit2,'LVegCov')
data.fit2$LVegHAI <- f.get.cov(data.fit2, 'LVegLAI')
data.fit2$LVegTyp <- as.factor(f.get.cov(data.fit2, 'LVegTyp'))
data.fit2$month.idx <- (data.fit2$time.idx-1)%%12 + 1
data.fit2$year.idx <- (data.fit2$time.idx-1)%/%12 + 1
data.fit2$E1 <- 1
data.fit2$LVegCov.grp <- inla.group(data.fit2$LVegCov, n = 20, method = "quantile")
data.fit2$HVegCov.grp <- inla.group(data.fit2$HVegCov, n = 20, method = "quantile")

set.seed(1234)
zero.cnt.month <- data.fit2%>%group_by(time.idx)%>% summarize(y_cnt=sum(y)) %>% filter(y_cnt==0)
data.fit2[data.fit2$time.idx %in% zero.cnt.month$time.idx,'y'] <- rbinom(8832,1,0.2)*10^-5

points_sf <- st_as_sf(data.fit2[data.fit2$time.idx==1,], coords = c("lon.grid", "lat.grid"), crs = "+proj=longlat +datum=WGS84 +no_defs")
points_sf <- st_transform(points_sf, projutm)
points_sf
points_sf <- points_sf[,c('grid.idx')] %>%
  mutate(
    geometry = st_buffer(geometry, dist = 13.5, endCapStyle = "SQUARE")
  )
plot(points_sf)
points_sf[points_sf$grid.idx==74,]

points_sf <- st_transform(points_sf, crs='+proj=longlat +datum=WGS84 +no_defs')
intersections <- st_intersection(points_sf, st_as_sf(district))

# Calculate the area of each intersection
intersections$area <- st_area(intersections)

intersections <- intersections %>% group_by(grid.idx)%>%
   arrange(desc(area)) %>% 
   distinct(grid.idx, .keep_all=T) %>%
   ungroup()

intersections[intersections$grid.idx==93,]
intersections$district <- intersections$NAME_1

# points_in_districts <- st_join(points_sf, st_as_sf(district))


# data.fit2$district <- as.character(points_in_districts$NAME_1)
# data.fit2[is.na(data.fit2$district),'district'] <- 'Out of Boundary'
data.fit2 <- merge(data.fit2,intersections[,c('grid.idx','district')], by='grid.idx')
data.fit2$district <- as.factor(data.fit2$district)

ggplot(intersections) + geom_sf(aes(fill=NAME_1))

# iset <- inla.spde.make.index('i', n.spde = spde$n.spde,
#                              n.group = 9)
# 
# A <- inla.spde.make.A(mesh = mesh,
#                       loc = cbind(data.fit2$grid.cent.x.utm, data.fit2$grid.cent.y.utm), 
#                       group = data.fit2$year.idx) 
# 
# 
# sdat <- inla.stack(
#   data = list(y = data.fit2$y), 
#   A = list(A), 
#   effects = list(iset),
#   tag = 'stdata') 
# 
# formulae <- y ~ 0 +  f(i, model = spde, group = i.group, 
#                           control.group = list(model = 'ar1')) 
# 
# # Model fitting
# res2 <- inla(formulae, family='poisson',  data = inla.stack.data(sdat), 
#             control.predictor = list(compute = TRUE,
#                                      A = inla.stack.A(sdat)), 
#             control.fixed = list(expand.factor.strategy = 'inla'),
#             verbose=TRUE)
# idat <- inla.stack.index(sdat, 'stdata')$data
# 
# data.fit2$Estimated_Lambda <- res$summary.fitted.values[idat, "mean"]
# data.fit2$Lower_Bound <- res$summary.fitted.values[idat, "0.025quant"]
# data.fit2$Upper_Bound <- res$summary.fitted.values[idat, "0.975quant"]
# data.fit2$Scaling_Residuals <- data.fit2$y - data.fit2$Estimated_Lambda
# summary(data.fit2$Scaling_Residuals)
# 
prior <- list(
  prec = list(
    prior = "pc.prec",
    # P(sqrt(sd) > 0.3)=0.01
    param = c(2, 0.01)),
  phi = list(
    prior = "pc",
    # P(phi < 0.5) = 0.1
    param = c(0.9, 0.01))
)




formula <- y ~ 1 +   f(idarea, model = "bym", graph = g,
                       group = time.idx, control.group = list(model = "ar", order=2))


                 
t1 <- Sys.time()
res <- inla(formula,
            family = "poisson", data = data.fit2, 
            control.predictor = list(compute = TRUE),
            verbose=TRUE,
            control.fixed = list(expand.factor.strategy = 'inla')
)
t2 <- Sys.time()
print(t2-t1)
summary(res)


formula1.1 <- y ~  1 +  f(idarea, model = "bym2", graph = g, group = time.idx, control.group = list(model = "ar", order=2),
                       hyper = prior, scale.model=TRUE)
t1 <- Sys.time()
res1.1 <- inla(formula1.1,
            family = "poisson", data = data.fit2, 
            control.predictor = list(compute = TRUE),
            verbose=TRUE,
            control.fixed = list(expand.factor.strategy = 'inla'),
            control.inla = list(strategy='adaptive', int.strategy = 'eb'),
)
t2 <- Sys.time()
print(t2-t1)
summary(res1.1)

formula1.2 <- y ~  1 +  f(idarea, model = "besag", graph = g, group = time.idx, control.group = list(model = "ar", order=2))
t1 <- Sys.time()
res1.2 <- inla(formula1.2,
               family = "poisson", data = data.fit2, 
               control.predictor = list(compute = TRUE),
               verbose=TRUE,
               control.fixed = list(expand.factor.strategy = 'inla')
)
t2 <- Sys.time()
print(t2-t1)
summary(res1.2)


result <- res1.1

data.fit2$Estimated_Lambda <- result$summary.fitted.values[, "mean"]
data.fit2$Lower_Bound <- result$summary.fitted.values[, "0.025quant"]
data.fit2$Upper_Bound <- result$summary.fitted.values[, "0.975quant"]
data.fit2$Scaling_Residuals <- data.fit2$y - data.fit2$Estimated_Lambda
summary(data.fit2$Scaling_Residuals)

B2_sf <- st_as_sf(B2.merge)
library(tidyr)
B2_sf <- gather(B2_sf, y.time.idx, y, paste0("y.", 1:9))
B2_sf$time.idx <- as.integer(substring(B2_sf$y.time.idx, 3, 5))

B2_sf <- merge(
  B2_sf, data.fit2[,c('grid.idx','time.idx','month.idx','year.idx','Estimated_Lambda','Lower_Bound','Upper_Bound','Scaling_Residuals')],
  by.x=c('grid.idx','time.idx'),
  by.y=c('grid.idx','time.idx'),
)


lprange.scale <- max(
  abs(data.fit2$Scaling_Residuals)
)*c(-1,1)

csc.scale <- scale_fill_gradientn(colours = brewer.pal(3, "RdBu"), limits = lprange.scale)

ggplot() + geom_sf(data=B2_sf,aes(fill = Scaling_Residuals),lwd = 0.1) +
  # geom_sf(data=st_as_sf(district),alpha = 0)+
  facet_wrap(~time.idx, dir = "h", ncol =3) +
  ggtitle(paste("Pois Residuals in")) +
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)
        # axis.text.y = element_blank(),
        # axis.ticks = element_blank()
  )+csc.scale



formula1 <- y ~  1 +  f(idarea, model = "bym", graph = g, 
                        group = time.idx, control.group = list(model = "ar", order=2)) 


res1 <- inla(formula1,
            family = "poisson", data = data.fit2, 
            control.predictor = list(compute = TRUE),
            verbose=TRUE,
            control.fixed = list(expand.factor.strategy = 'inla')
)


formula2 <- y ~  1 +  f(idarea, model = "bym", graph = g, 
                        group = time.idx, control.group = list(model = "ar", order=2)) + 
                       f(LVegCov, model='linear') + f(HVegCov, model='linear')
                      
                          



res2 <- inla(formula2,
             family = "poisson", data = data.fit2, 
             control.predictor = list(compute = TRUE),
             verbose=TRUE,
             control.fixed = list(expand.factor.strategy = 'inla')
)


data.fit2$Estimated_Lambda <- res2$summary.fitted.values[, "mean"]
data.fit2$Lower_Bound <- res2$summary.fitted.values[, "0.025quant"]
data.fit2$Upper_Bound <- res2$summary.fitted.values[, "0.975quant"]
data.fit2$Scaling_Residuals <- data.fit2$y - data.fit2$Estimated_Lambda
summary(data.fit2$Scaling_Residuals)


formula2.1 <- y ~  1 +  f(idarea, model = "bym", graph = g, 
                        group = time.idx, control.group = list(model = "ar", order=2)) + 
  f(LVegHAI, model='linear') + f(HVegHAI, model='linear') 





res2.1 <- inla(formula2.1,
             family = "poisson", data = data.fit2, 
             control.predictor = list(compute = TRUE),
             verbose=TRUE,
             control.fixed = list(expand.factor.strategy = 'inla')
)

data.fit2$Estimated_Lambda <- res2.1$summary.fitted.values[, "mean"]
data.fit2$Lower_Bound <- res2.1$summary.fitted.values[, "0.025quant"]
data.fit2$Upper_Bound <- res2.1$summary.fitted.values[, "0.975quant"]
data.fit2$Scaling_Residuals <- data.fit2$y - data.fit2$Estimated_Lambda
summary(data.fit2$Scaling_Residuals)




formula2.2 <- y ~  1 +  f(idarea, model = "bym", graph = g, 
                          group = time.idx, control.group = list(model = "ar", order=2)) + 
  f(LVegHAI, model='linear') + f(HVegHAI, model='linear') + 
  f(LVegCov, model='linear') + f(HVegCov, model='linear')





res2.2 <- inla(formula2.2,
               family = "poisson", data = data.fit2, 
               control.predictor = list(compute = TRUE),
               verbose=TRUE,
               control.fixed = list(expand.factor.strategy = 'inla')
)

data.fit2$Estimated_Lambda <- res2.2$summary.fitted.values[, "mean"]
data.fit2$Lower_Bound <- res2.2$summary.fitted.values[, "0.025quant"]
data.fit2$Upper_Bound <- res2.2$summary.fitted.values[, "0.975quant"]
data.fit2$Scaling_Residuals <- data.fit2$y - data.fit2$Estimated_Lambda
summary(data.fit2$Scaling_Residuals)


formula2.3 <- y ~  0 +  f(idarea, model = "bym", graph = g, 
                          group = time.idx, control.group = list(model = "ar", order=2), hyper=list()) + 
                         district
                          




res2.3 <- inla(formula2.3,
               family = "poisson", data = data.fit2, 
               control.predictor = list(compute = TRUE),
               verbose=TRUE,
               control.fixed = list(expand.factor.strategy = 'inla')
)



data.fit2$Estimated_Lambda <- res2.3$summary.fitted.values[, "mean"]
data.fit2$Lower_Bound <- res2.3$summary.fitted.values[, "0.025quant"]
data.fit2$Upper_Bound <- res2.3$summary.fitted.values[, "0.975quant"]
data.fit2$Scaling_Residuals <- data.fit2$y - data.fit2$Estimated_Lambda
summary(data.fit2$Scaling_Residuals)


mesh1D.Cov <- fm_mesh_1d(seq(-0.2,1.2,by=0.05), boundary = "free")
spde1D.Cov <- inla.spde2.pcmatern(mesh1D.Cov,
                                   prior.range = c(0.005, 0.1),
                                   prior.sigma = c(0.1, 0.1)
)


formula3 <- y ~  1 +  f(idarea, model = "bym", graph = g, 
                        group = time.idx, control.group = list(model = "ar", order=2)) + 
  f(LVegCov.grp,  model='rw1') + f(HVegCov.grp, model='rw1')




res3 <- inla(formula3,
            family = "poisson", data = data.fit2, E=E1,
            control.predictor = list(compute = TRUE),
            verbose=TRUE,
            control.fixed = list(expand.factor.strategy = 'inla')
)

data.fit2$Estimated_Lambda <- res3$summary.fitted.values[, "mean"]
data.fit2$Lower_Bound <- res3$summary.fitted.values[, "0.025quant"]
data.fit2$Upper_Bound <- res3$summary.fitted.values[, "0.975quant"]
data.fit2$Scaling_Residuals <- data.fit2$y - data.fit2$Estimated_Lambda
summary(data.fit2$Scaling_Residuals)



formula4<- y ~  1 +  f(idarea, model = "bym", graph = g, 
                        group = time.idx, control.group = list(model = "ar", order=2)) + 
  f(LVegCov, model='linear') + f(HVegCov, model='linear')



res4 <- inla(formula4,
             family = "zeroinflatedpoisson0", data = data.fit2, E=E1,
             control.predictor = list(compute = TRUE),
             verbose=TRUE,
             control.fixed = list(expand.factor.strategy = 'inla')
)

data.fit2$Estimated_Lambda <- res4$summary.fitted.values[, "mean"]
data.fit2$Lower_Bound <- res4$summary.fitted.values[, "0.025quant"]
data.fit2$Upper_Bound <- res4$summary.fitted.values[, "0.975quant"]
data.fit2$Scaling_Residuals <- data.fit2$y - data.fit2$Estimated_Lambda
summary(data.fit2$Scaling_Residuals)

formula5<- y ~  1 +  f(idarea, model = "bym", graph = g, 
                       group = time.idx, control.group = list(model = "ar", order=2)) + 
  f(LVegCov, model='linear') + f(HVegCov, model='linear')



res5 <- inla(formula5,
             family = "zeroinflatedpoisson1", data = data.fit2, E=E1,
             control.predictor = list(compute = TRUE),
             verbose=TRUE,
             control.fixed = list(expand.factor.strategy = 'inla')
)

formula6<- y ~  1 +  f(idarea, model = "bym", graph = g, 
                       group = time.idx, control.group = list(model = "ar", order=2)) 



res6 <- inla(formula6,
             family = "nbinomial", data = data.fit2, E=E1,
             control.predictor = list(compute = TRUE),
             verbose=TRUE,
             control.fixed = list(expand.factor.strategy = 'inla'),
             control.family = list(prior="gaussian", param = c(0,0.01))
)

####################################################################

for (t in  c(1,seq(3,30,by=3))){
  data.fit <- data.merge %>% st_drop_geometry() %>% filter(month.idx==9 & length >=t*60 )%>% 
              group_by(grid.idx,grid_id, year.idx) %>% 
              summarise(y = n(), 
                      lon.grid = mean(lon.grid),
                      lat.grid = mean(lat.grid),
                      x.utm = mean(grid.cent.x.utm),
                      y.utm = mean(grid.cent.y.utm),
                      E = mean(E))
  
  data.fit2 <- do.call(rbind, lapply(1:9, function(x) {B2@data$year.idx = x 
                                                        return(B2@data)}))
  print(dim(data.fit2))
  data.fit2 <- merge(data.fit2, data.fit[,c('grid.idx','year.idx','y')],
                     by=c('grid.idx','year.idx'),all.x=T)
  print(dim(data.fit2))
  data.fit2[is.na(data.fit2$y),'y'] <- 0
  summary(data.fit2$y)
  
  
  
  data.fit3 <- reshape(data.fit2[,c('grid.idx','year.idx','y')],
                timevar = "year.idx",
                idvar = "grid.idx",
                direction = "wide"
  )
  
  B2.merge <- merge(B2, data.fit3, by = "grid.idx")
  # ggplot(st_as_sf(B2)) + geom_sf(aes(fill=y.7))
  
  B2_sf <- st_as_sf(B2.merge)
  library(tidyr)
  B2_sf <- gather(B2_sf, year.idx, y, paste0("y.", 1:9))
  B2_sf$time.idx <- as.integer(substring(B2_sf$year.idx, 3, 3))
  
  # ggplot(B2_sf) + geom_sf(aes(fill = y)) +
  #   facet_wrap(~year.idx, dir = "h", ncol = 3) +
  #   ggtitle("Fire ignition") + theme_bw() +
  #   theme(
  #     axis.text.x = element_blank(),
  #     axis.text.y = element_blank(),
  #     axis.ticks = element_blank()
  #   ) +
  #   scale_fill_gradient2(
  #     midpoint = 1, low = "blue", mid = "white", high = "red"
  #   )
  
  
  B2.adj <- poly2nb(B2)
  nb2INLA("map.adj", B2.adj)
  g <- inla.read.graph(filename = "map.adj")
  
  data.fit2$idarea <- as.numeric(as.factor(data.fit2$grid.idx))
  data.fit2$idarea1 <- data.fit2$idarea
  
  formula <- y ~ 1 + f(idarea, model = "bym", graph = g, group = year.idx, control.group = list(model = "ar1"))
  
  res <- inla(formula,
              family = "poisson", data = data.fit2, 
              control.predictor = list(compute = TRUE),
              verbose=TRUE
  )
  
  
  data.fit2$Estimated_Lambda <- res2$summary.fitted.values[, "mean"]
  data.fit2$Lower_Bound <- res2$summary.fitted.values[, "0.025quant"]
  data.fit2$Upper_Bound <- res2$summary.fitted.values[, "0.975quant"]
  data.fit2$Scaling_Residuals <- data.fit2$y - data.fit2$Estimated_Lambda
  
  result <- data.fit2[,c('grid.idx','year.idx','y','Estimated_Lambda','Lower_Bound','Upper_Bound','Scaling_Residuals')]
  result$cutoff <- t
  
  if (t==1){
    result.all <- result
  }else{
    result.all <- rbind(result.all, result)
  }
}


# 
data.fit2$Estimated_Lambda <- res2$summary.fitted.values[, "mean"]
data.fit2$Lower_Bound <- res2$summary.fitted.values[, "0.025quant"]
data.fit2$Upper_Bound <- res2$summary.fitted.values[, "0.975quant"]

# data.fit2$Estimated_Lambda <- res2$summary.fitted.values[idat, "mean"]
# data.fit2$Lower_Bound <- res2$summary.fitted.values[idat, "0.025quant"]
# data.fit2$Upper_Bound <- res2$summary.fitted.values[idat, "0.975quant"]

data.fit2$Scaling_Residuals <- data.fit2$y - data.fit2$Estimated_Lambda
summary(data.fit2$Scaling_Residuals)


B2_sf <- st_as_sf(B2.merge)
library(tidyr)
B2_sf <- gather(B2_sf, y.time.idx, y, paste0("y.", 1:108))
B2_sf$time.idx <- as.integer(substring(B2_sf$y.time.idx, 3, 5))

B2_sf <- merge(
  B2_sf, data.fit2[,c('grid.idx','time.idx','month.idx','year.idx','Estimated_Lambda','Lower_Bound','Upper_Bound','Scaling_Residuals')],
  by.x=c('grid.idx','time.idx'),
  by.y=c('grid.idx','time.idx'),
)


# B2_sf <- st_as_sf(B2.merge)
# library(tidyr)
# B2_sf <- gather(B2_sf, y.time.idx, y, paste0("y.", 1:9))
# B2_sf$time.idx <- as.integer(substring(B2_sf$y.time.idx, 3, 5))
# 
# 
# B2_sf <- merge(
#   B2_sf, data.fit2[,c('grid.idx','year.idx','Estimated_Lambda','Lower_Bound','Upper_Bound','Scaling_Residuals')],
#   by.x=c('grid.idx','time.idx'),
#   by.y=c('grid.idx','year.idx'),
# )
# 




year <- 2
lprange.scale <- max(
  abs(data.fit2$Scaling_Residuals)
  # abs(data.fit2[data.fit2$year.idx==year, ]$y - res2$summary.fitted.values[data.fit2$year.idx==year, "mean"] )
  # abs(data.fit2[data.fit2$year.idx==year, ]$y )
)*c(-1,1)

# require(cowplot)
# lprange.scale <- max(
#   abs(res_LGCP_1[res_LGCP_1$Type=='Scaling Residuals',]$mean),
#   abs(data.fit2$Scaling_Residuals )
# )*c(-1,1)

csc.scale <- scale_fill_gradientn(colours = brewer.pal(3, "RdBu"), limits = lprange.scale)

ggplot() + geom_sf(data=B2_sf[B2_sf$year.idx==year,],aes(fill = Scaling_Residuals),lwd = 0.1) +
  geom_sf(data=st_as_sf(district),alpha = 0)+
  facet_wrap(~month.idx, dir = "h", ncol =3) +
  ggtitle(paste("Pois Residuals in", year + 2011)) +
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)
    # axis.text.y = element_blank(),
    # axis.ticks = element_blank()
  )+
  csc.scale

lprange.scale <- max(
  # abs(data.fit2$y - res$summary.fitted.values[, "mean"])
  # abs(data.fit2[data.fit2$year.idx==year, ]$y - res2$summary.fitted.values[data.fit2$year.idx==year, "mean"] )
  abs(data.fit2[data.fit2$year.idx==year, ]$y )
)*c(-1,1)
csc.scale <- scale_fill_gradientn(colours = brewer.pal(3, "RdBu"), limits = lprange.scale)
ggplot() + geom_sf(data=B2_sf[B2_sf$year.idx==year,],aes(fill = y)) +
  geom_sf(data=st_as_sf(district),alpha = 0)+
  facet_wrap(~month.idx, dir = "h", ncol =3) +
  ggtitle(paste("Observations in", year + 2011)) +
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)
        # axis.text.y = element_blank(),
        # axis.ticks = element_blank()
  )+
  csc.scale



# B2_sf$Residual <- B2_sf$Scaling_Residuals
# ggplot(B2_sf) + geom_sf(aes(fill = Residual)) +
#   facet_wrap(~time.idx, dir = "h", ncol =3) +
#   ggtitle("Residuals of BYM model") +
#   theme(legend.position = "bottom",
#         axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)
#     # axis.text.y = element_blank(),
#     # axis.ticks = element_blank()
#   )+
#   csc.scale

samplers.all$Residual <- samplers.all$Scaling_Residual
p2 <- ggplot(st_as_sf(samplers.all)) +
  geom_sf(aes(fill = Residual)) +
  # csc["Scaling"] +
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)
        # axis.text.y = element_blank(),
        # axis.ticks = element_blank()
  )+
  facet_wrap(~time.idx) + 
  ggtitle("Residuals of SPDE model")  + csc.scale 

summary(B2_sf[B2_sf$year.idx==year,'y'])


ggplot(B2_sf[B2_sf$year.idx==year,]) + geom_sf(aes(fill = y)) +
  facet_wrap(~month.idx, dir = "h", ncol =3) +
  ggtitle("Scaling Residuals of SPDE model on the lattices") + theme_bw() +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_blank()
  )+
  csc.scale
summary(B2_sf[B2_sf$year.idx==year,'y'])

ggplot(B2_sf[B2_sf$year.idx==year,]) + geom_sf(aes(fill = Estimated_Lambda)) +
  facet_wrap(~month.idx, dir = "h", ncol =3) +
  ggtitle("Scaling Residuals of SPDE model on the lattices") + theme_bw() +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_blank()
  )+  csc.scale


result.all$cutoff <- as.factor(result.all$cutoff)
result.all$Pearson_Residuals <- result.all$y/sqrt(result.all$Estimated_Lambda) - sqrt(result.all$Estimated_Lambda)
ggplot(result.all, aes(x=cutoff,y=Pearson_Residuals,group=cutoff)) + 
  geom_boxplot()+
  scale_x_discrete(limits=factor(c(1,seq(3,30,by=3))) , labels=factor(c(1,seq(3,30,by=3))))

result.all.sorted <- result.all %>% group_by(cutoff) %>%
  arrange(year.idx, y, .by_group = TRUE)


result.all.sorted$plot.idx <- rep(1:1728, 11) 
result.all.sorted$Estimated_Lambda <- result.all.sorted$y - result.all.sorted$Scaling_Residuals

cf <- 18
ggplot(data=result.all.sorted[result.all.sorted$cutoff==cf,] )+
  geom_line(aes(x=plot.idx, y=Estimated_Lambda)) + 
  geom_line(aes(x=plot.idx, y=y), col='red') + 
  geom_vline(xintercept = seq(192,length.out=9, by=192), linetype="dotted", 
             color = "blue", size=0.5)+
  geom_text(aes(x= 96, label="2012", y=-5), colour="red")+
  geom_text(aes(x= 288, label="2013", y=-5), colour="red")+
  geom_text(aes(x= 480, label="2014", y=-5), colour="red")+
  geom_text(aes(x= 672, label="2015", y=-5), colour="red")+
  geom_text(aes(x= 864, label="2016", y=-5), colour="red")+
  geom_text(aes(x= 1056, label="2017", y=-5), colour="red")+
  geom_text(aes(x= 1248, label="2018", y=-5), colour="red")+
  geom_text(aes(x= 1440, label="2019", y=-5), colour="red")+
  geom_text(aes(x= 1632, label="2020", y=-5), colour="red") + 
  ylim(c(-10,40))+
  ggtitle(paste('Cutoff = ', cf))




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

district <- dist[dist$NAME_0=="Portugal",]