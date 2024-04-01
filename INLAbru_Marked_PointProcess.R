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
library(terra)
library(RColorBrewer)
library(patchwork)
library(raster)
library(rgdal)

# bru_safe_sp(force = TRUE)
load(file.path(dir.data, "burn_area","finaldata_urbanfire.RData"))

data.new$year.idx <- data.new$year -2011
loc.data.utm <- st_as_sf(data.new, coords=c('x_utm_new','y_utm_new'), crs='+proj=utm +zone=29 +datum=WGS84 +units=km +no_defs' )

# coords <-  st_coordinates(st_as_sf(loc.data.utm))
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



r <- 50
p1 <- 0.1
sigma <- 8
p2 <- 0.1

point_matern <- inla.spde2.pcmatern(mesh = mesh,
                            # PC-prior on range: P(practic.range < r) = p1
                            prior.range = c(r, p1),
                            # PC-prior on sigma: P(sigma > sigma) = p2
                            prior.sigma = c(sigma, p2))

mark_matern <- inla.spde2.pcmatern(mesh = mesh,
                                  # PC-prior on range: P(practic.range < r) = p1
                                  prior.range = c(r, p1),
                                  # PC-prior on sigma: P(sigma > sigma) = p2
                                  prior.sigma = c(sigma, p2))



cmp <-   ~  -1 + 
  point_field(coordinates, model=point_matern) + 
  mark_field(coordinates, model=mark_matern) + 
  Inter_point(1) + Inter_mark(1) 

point_lik <- like("cp",
                     formula=coordinates ~ point_field + Inter_point,
                     include=c("point_field","Inter_point"),
                     data=coords,
                     domain=list(coordinates=mesh),
                     samplers=domainSP)
mark_lik <- like("lognormal",
                    formula=area_ha ~ Inter_mark + mark_field + point_field,
                    include=c("Inter_mark","point_field","mark_field"),
                    data=coords)

t1 <-  Sys.time()
fit <- bru(cmp, 
           point_lik,
           mark_lik,
          options=list(verbose=TRUE))
t2 <-  Sys.time()
print(t2-t1)

pred_point<- predict(
  fit,
  fm_pixels(mesh, mask = domainSP, format = "sp"),
  ~ exp(point_field + Inter_point)
)
pred_mark<- predict(
  fit,
  fm_pixels(mesh, mask = domainSP, format = "sp"),
  ~ exp(Inter_mark + mark_field + point_field)
  )

pl1 <- ggplot() +
  gg(pred_mark, aes(fill = q0.975)) +
  gg(domainSP) +
  ggtitle("Burn area fitted in Mared PP", subtitle = "(Response Scale)")


mark_lik1 <- like("lognormal",
                 formula=area_ha ~ Inter_mark + mark_field,
                 include=c("Inter_mark","mark_field"),
                 data=coords)
t1 <-  Sys.time()
fit1 <- bru(cmp, 
           # point_lik,
           mark_lik1,
           options=list(verbose=TRUE))
t2 <-  Sys.time()
print(t2-t1)
pred_mark1<- predict(
  fit1,
  fm_pixels(mesh, mask = domainSP, format = "sp"),
  ~ exp(Inter_mark + mark_field )
)

pl2 <- ggplot() +
  gg(pred_mark1, aes(fill = q0.975)) +
  gg(domainSP) +
  ggtitle("Burn area fitted only", subtitle = "(Response Scale)")
pl2

multiplot(pl1, pl2, cols = 2)
