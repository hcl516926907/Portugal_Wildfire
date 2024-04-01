dir.data <- "/home/pgrad2/2448355h/My_PhD_Project/00_Dataset/Urban_Fires"
dir.out <- '/home/pgrad2/2448355h/My_PhD_Project/01_Output/Urban_Fire_Model'


library(INLA)
library(mgcv)
library(ggplot2)
library(rgdal)
library(sf)
library(tibble)
library(dplyr)

book.mesh.dual <- function(mesh) {
  if (mesh$manifold=='R2') {
    ce <- t(sapply(1:nrow(mesh$graph$tv), function(i)
      colMeans(mesh$loc[mesh$graph$tv[i, ], 1:2])))
    library(parallel)
    pls <- mclapply(1:mesh$n, function(i) {
      p <- unique(Reduce('rbind', lapply(1:3, function(k) {
        j <- which(mesh$graph$tv[,k]==i)
        if (length(j)>0) 
          return(rbind(ce[j, , drop=FALSE],
                       cbind(mesh$loc[mesh$graph$tv[j, k], 1] +
                               mesh$loc[mesh$graph$tv[j, c(2:3,1)[k]], 1], 
                             mesh$loc[mesh$graph$tv[j, k], 2] +
                               mesh$loc[mesh$graph$tv[j, c(2:3,1)[k]], 2])/2))
        else return(ce[j, , drop=FALSE])
      })))
      j1 <- which(mesh$segm$bnd$idx[,1]==i)
      j2 <- which(mesh$segm$bnd$idx[,2]==i)
      if ((length(j1)>0) | (length(j2)>0)) {
        p <- unique(rbind(mesh$loc[i, 1:2], p,
                          mesh$loc[mesh$segm$bnd$idx[j1, 1], 1:2]/2 +
                            mesh$loc[mesh$segm$bnd$idx[j1, 2], 1:2]/2, 
                          mesh$loc[mesh$segm$bnd$idx[j2, 1], 1:2]/2 +
                            mesh$loc[mesh$segm$bnd$idx[j2, 2], 1:2]/2))
        yy <- p[,2]-mean(p[,2])/2-mesh$loc[i, 2]/2
        xx <- p[,1]-mean(p[,1])/2-mesh$loc[i, 1]/2
      }
      else {
        yy <- p[,2]-mesh$loc[i, 2]
        xx <- p[,1]-mesh$loc[i, 1]
      }
      Polygon(p[order(atan2(yy,xx)), ])
    })
    return(SpatialPolygons(lapply(1:mesh$n, function(i)
      Polygons(list(pls[[i]]), i))))
  }
  else stop("It only works for R2!")
}

load(file.path(dir.data, "burn_area","finaldata_urbanfire.RData"))

# time.sp <- as.character(deframe(data[,'open']))>'2017-01-01'
# data.flt <- data[time.sp,c("lon", "lat",'cent.lon','cent.lat')]
# dim(data.flt %>% count(lat,lon))
# data.cox <- data.flt %>% count(lon,lat)

# n1 <- nrow(data %>% count(lon,lat))
# 
# set.seed(1234)
# data$year <- format(data$open, format="%Y")
# data$month <- format(data$open, format="%m")
# data %>% count(year, month,district)
# 
# data.flt <- data %>% group_by(year, district) %>% sample_frac(0.4)
# data.cox <- data.flt %>% count(lon,lat)
# n2 <- nrow(data.cox);n2
# 
# loc.data <- st_as_sf(data.cox, coords = c("lon", "lat"), crs = 4326)

# loc.data.utm <- st_transform(st_as_sf(loc.data), "+proj=utm +zone=29 +datum=WGS84 +units=km +no_defs")
loc.data.utm <- st_as_sf(data.new, coords=c('x_utm_new','y_utm_new'), crs='+proj=utm +zone=29 +datum=WGS84 +units=km +no_defs' )

coo <- st_coordinates(st_as_sf(loc.data.utm))

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
library(terra)

# raster grid covering map
grid.pred <- terra::rast(map_mainland, nrows = 100, ncols = 60)
# coordinates of all cells
xy <- terra::xyFromCell(grid.pred, 1:ncell(grid.pred))

dp <- st_as_sf(as.data.frame(xy), coords = c("x", "y"),
               crs = st_crs(map))
indicespointswithin <- which(st_intersects(dp, map_mainland,
                                           sparse = FALSE))

# points within the map
dp <- st_filter(dp, map_mainland)

ggplot() + geom_sf(data = map_mainland) +
  geom_sf(data = dp) + coord_sf(datum = projutm)

coop <- st_coordinates(dp)

# coords <- coordinates(loc.data)
# hull_indices <- chull(coords)
# hull_coords <- coords[hull_indices, ]
# 
# # Create a Polygon and then a SpatialPolygons object
# poly <- Polygon(hull_coords)
# sp_poly <- SpatialPolygons(list(Polygons(list(poly), "convex_hull")))



# inner.length <- 8
# outer.length <- 20
# cutoff <- 4

# inner.length <- 16
# outer.length <- 40
# cutoff <- 4

# inner.length <- 20
# outer.length <- 50
# cutoff <- 6

inner.length <- 40
outer.length <- 70
cutoff <- 15

# plot points
# mesh <- inla.mesh.2d(loc.domain = loc.data, boundary=sp_poly, max.edge=c(inner.length,outer.length), cutoff=0.1)
mesh <- inla.mesh.2d(loc.domain = loc.d, max.edge=c(inner.length,outer.length), cutoff=cutoff)
mesh$n
plot(mesh, asp = 1, main = '')
points(loc.data.utm, col = 4, pch = 19, cex=0.2)

# Number of vertices in the mesh
nv <- mesh$n
# Number of points in the data
n <- nrow(coo)

r <- 50
p1 <- 0.1
sigma <- 8
p2 <- 0.1

spde <- inla.spde2.pcmatern(mesh = mesh,
                            prior.range = c(r, p1),

                            prior.sigma = c(sigma, p2)) 



# library("deldir")
# library("SDraw")

dmesh <- book.mesh.dual(mesh)
plot(dmesh)
axis(1)
axis(2)

domain.polys <- Polygons(list(Polygon(loc.d)), '0')
domainSP <- SpatialPolygons(list(domain.polys))

# Voronoi polygons (as SpatialPolygons)
# mytiles <- voronoi.polygons(SpatialPoints(mesh$loc[, 1:2]))


require(rgeos)

w <- sapply(1:length(dmesh), function(i) {
  if (gIntersects(dmesh[i, ], domainSP))
    return(gArea(gIntersection(dmesh[i, ], domainSP)))
  else return(0)
})

sum(w)
st_area(map_mainland)
gArea(domainSP)

plot(mesh)
points(mesh$loc[which(w > 0), 1:2], col = "black", pch = 20, cex=0.2)
points(mesh$loc[which(w == 0), 1:2], col = "red", pch = 20, cex=0.2)

k <- 9
mesh.t <- inla.mesh.1d(2012:2020)
st.vol <- rep(w, k) * rep(diag(inla.mesh.fem(mesh.t)$c0), nv)

y.pp = rep(0:1, c(k*nv, n))
# y.pp = c(rep(0,nv), data.cox$n)
e.pp = c(st.vol, rep(0, n))

# Projection matrix for the integration points (mesh vertices)
# A.int <- Diagonal(nv, rep(1, nv))
A.int <- Diagonal(k*nv, rep(1, nv))
# Projection matrix for observed points (event locations)
# A.y <- inla.spde.make.A(mesh = mesh, loc = coo)
A.y <- inla.spde.make.A(mesh = mesh, loc = coo,
                        n.group = length(mesh.t$n),
                        group= data.new$year,group.mesh = mesh.t)
print(dim(A.y))

idx <- inla.spde.make.index('s', spde$n.spde, n.group = mesh.t$n)

# Projection matrix for mesh vertices and event locations
A.pp <- rbind(A.int, A.y)
Ap.pp <- inla.spde.make.A(mesh = mesh, loc = do.call("rbind", rep(list(coop), k)),
                          n.group = length(mesh.t$n),
                          group= rep(2012:2020, each=nrow(coop)),
                          group.mesh = mesh.t)

# stack for estimation
stk.e.pp <- inla.stack(tag = "est.pp",
                       data = list(y = y.pp, e = e.pp), 
                       A = list( A.pp,1),
                       effects = list(idx, list(b0 = rep(1, k*nv + n))))


# stack for prediction stk.p
stk.p.pp <- inla.stack(tag = "pred.pp",
                       data = list(y = rep(NA, nrow(coop)*k), e = rep(0, nrow(coop)*k)),
                       A = list( Ap.pp,1),
                       effects = list(idx, data.frame(b0 = rep(1, nrow(coop)*k)))
                       )
                                      

stk.full.pp <- inla.stack(stk.e.pp, stk.p.pp)

#P(cor > 0.7) =0.7
pcrho <- list(prior = 'pc.cor1', param = c(0.7, 0.7))
formula <- y ~ 0 + b0 + f(s, model = spde, group = s.group,
                          control.group = list(model = 'ar1',
                                               hyper = list(rho = pcrho)))

t1 <- Sys.time()
res.inla <- inla(formula,  family = 'poisson',
            data = inla.stack.data(stk.full.pp),
            control.predictor = list(compute = TRUE, link = 1,
                                     A = inla.stack.A(stk.full.pp)),
            E = inla.stack.data(stk.full.pp)$e,
            control.compute = list(waic=TRUE, config=TRUE),
            verbose=TRUE)
t2 <- Sys.time()
print(t2-t1)

# range is changed when using different mesh?


res <- res.inla
r0 <- diff(range(loc.d[, 1])) / diff(range(loc.d[, 2]))
prj <- inla.mesh.projector(mesh, xlim = range(loc.d[, 1]),
                           ylim = range(loc.d[, 2]), dims = c(100, 100 / r0)) 
ov <- over(SpatialPoints(prj$lattice$loc), domainSP)
m.prj <- lapply(1:k, function(j) {
  r <- inla.mesh.project(prj,
                         res$summary.ran$s$mean[1:nv + (j - 1) * nv])
  r[is.na(ov)] <- NA
  return(r) 
})

index <- inla.stack.index(stk.full.pp, tag = "pred.pp")$data
pred_mean <- res$summary.fitted.values[index, "mean"]
pred_ll <- res$summary.fitted.values[index, "0.025quant"]
pred_ul <- res$summary.fitted.values[index, "0.975quant"]


# 
# grid$mean <- NA
# grid$ll <- NA
# grid$ul <- NA
# 
# grid$mean[indicespointswithin] <- pred_mean
# grid$ll[indicespointswithin] <- pred_ll
# grid$ul[indicespointswithin] <- pred_ul

grid.t <- list()
for (i in 1:k){
  grid.t[[i]] <- grid.pred
  grid.t[[i]]$mean <- NA
  grid.t[[i]]$ll <- NA
  grid.t[[i]]$ul <- NA
  
  pred.idx <- which(rep(2012:2020, each=nrow(coop)) == (i+2011))
  grid.t[[i]]$mean[indicespointswithin] <- pred_mean[pred.idx]
  grid.t[[i]]$ll[indicespointswithin] <- pred_ll[pred.idx]
  grid.t[[i]]$ul[indicespointswithin] <- pred_ul[pred.idx]
  
}


library(rasterVis)
levelplot(raster::brick(grid.t[[1]]), layout = c(3, 1),
          names.attr = c("Mean", "2.5 percentile", "97.5 percentile"),par.settings= BuRdTheme())

r_df <- as.data.frame(grid.t[[1]], xy=TRUE)

# Use ggplot
ggplot(data = r_df, aes(x = x, y = y, fill = mean)) +
  geom_tile() +
  scale_fill_viridis_c() +  # Using viridis color scale
  coord_fixed() +  # Keep aspect ratio
  labs(fill = "Value", x = "Longitude", y = "Latitude") +
  theme_minimal()

for (i in 1:k){
  r_df <- as.data.frame(grid.t[[i]], xy=TRUE)
  
  # Use ggplot
  assign(paste('p',i, sep=''),  ggplot(data = r_df, aes(x = x, y = y, fill = mean)) +
    geom_tile() +
    scale_fill_viridis_c() +  # Using viridis color scale
    coord_fixed() +  # Keep aspect ratio
    labs(fill = "Value", x = "Longitude", y = "Latitude") +
    ggtitle(paste('Year', i+2011))+
  theme_minimal() )

}
library(gridExtra)
grid.arrange(p1,p2,p3,
             p4,p5,p6,
             p7,p8,p9, nrow = 3)
# library(fields)
# local.plot.field = function(field, mesh, xlim=c(-15,-4), ylim=c(35,43), ...){
#   stopifnot(length(field) == mesh$n)
#   proj = inla.mesh.projector(mesh, xlim = xlim, 
#                              ylim = ylim, dims=c(300, 300))
#   field.proj = inla.mesh.project(proj, field)
#   n.col = 20
#   image.plot(list(x = proj$x, y=proj$y, z = field.proj), 
#              xlim = xlim, ylim = ylim, col = plasma(n.col), nlevel=n.col+1, ...)
# }
# local.plot.field(res$summary.random[['s']][['mean']], mesh)
# len = res$summary.hyperpar[2, '0.5quant']
# # - the posterior median range
# arrows(5-0.5*len, 5, 5+0.5*len, 5, length=0.05, angle=90, code=3, lwd=3)
# 
# 
# 
# # Extract the estimated intensity values
# estimated_intensity <- pp.res0$summary.random$spatial.field[, "mean"]
# 
# # Create a data frame for plotting
# plot_data <- data.frame(coords = mesh$loc[, 1:2], intensity = estimated_intensity)
# 
# # Convert to a spatial object (if necessary) and then to an sf object for ggplot2 plotting
# spatial_points <- SpatialPointsDataFrame(coords = plot_data[,1:2], data = plot_data)
# plot_data_sf <- st_as_sf(spatial_points, coords = c("coords.X1", "coords.X2"), crs = 4326)
# 
# # Plot using ggplot2
# ggplot(data = plot_data_sf) +
#   geom_sf(aes(color = intensity, size = intensity), alpha = 0.6) +
#   scale_color_viridis_c() + # Use a color scale that represents intensity well
#   scale_size(range = c(1, 5)) + # Adjust size range based on your data
#   theme_minimal() +
#   labs(title = "Estimated Intensity of Log Gaussian Cox Model",
#        color = "Intensity",
#        size = "Intensity")
# 
# 
# 
# spde.est <- inla.spde2.result(inla = res, name = "spatial.field",
#                               spde = spde, do.transf = TRUE)
# inla.zmarginal(spde.est$marginals.variance.nominal[[1]])
# inla.zmarginal(spde.est$marginals.range.nominal[[1]])
