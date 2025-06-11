dir.data <- "/home/pgrad2/2448355h/My_PhD_Project/00_Dataset/Urban_Fires"
dir.out <- '/home/pgrad2/2448355h/My_PhD_Project/01_Output/Urban_Fire_Model'


library(INLA)
library(mgcv)
library(ggplot2)
library(rgdal)
library(sf)
library(tibble)

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

load(file.path(dir.data,"finaldata.RData"))

time.sp <- as.character(deframe(data[,'open']))>'2018-01-01'
data.flt <- data[time.sp,c("lon", "lat")]
# loc.data <- SpatialPointsDataFrame(coords = data.flt[, c("lon", "lat")], data = data.flt, proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0 +units=km"))
loc.data <- st_as_sf(data.flt, coords = c("lon", "lat"), crs = 4326)

loc.data.utm <- st_transform(st_as_sf(loc.data), "+proj=utm +zone=29 +datum=WGS84 +units=m +no_defs")

coo <- st_coordinates(st_as_sf(loc.data.utm))

library(rnaturalearth)
map <- ne_countries(type = "countries", country = "Portugal",
                           scale = "medium", returnclass = "sf")
projutm <- "+proj=utm +zone=29 +datum=WGS84 +units=m +no_defs"
map <- st_transform(map, crs = projutm)
mainland_bbox <- st_as_sfc(st_bbox(c(xmin = 0, xmax = 1000000, ymin = 4000000 , ymax = 4800000), crs = st_crs(map)))
map_mainland <- st_intersection(map, mainland_bbox)

ggplot() + geom_sf(data = map_mainland) +
  geom_sf(data = st_as_sf(loc.data.utm)) + coord_sf(datum = projutm)

loc.d <- cbind(st_coordinates(map_mainland)[, 1], st_coordinates(map_mainland)[, 2])
library(terra)

# raster grid covering map
grid <- terra::rast(map_mainland, nrows = 50, ncols = 30)
# coordinates of all cells
xy <- terra::xyFromCell(grid, 1:ncell(grid))

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



inner.length <- 2800
outer.length <- 10000


# plot points
# mesh <- inla.mesh.2d(loc.domain = loc.data, boundary=sp_poly, max.edge=c(inner.length,outer.length), cutoff=0.1)
mesh <- inla.mesh.2d(loc.domain = loc.d, max.edge=c(inner.length,outer.length))
mesh$n
plot(mesh, asp = 1, main = '')
points(loc.data.utm, col = 4, pch = 19, cex=0.2)

# Number of vertices in the mesh
nv <- mesh$n
# Number of points in the data
n <- nrow(coo)

spde <- inla.spde2.pcmatern(mesh = mesh,
                            # PC-prior on range: P(practic.range < 0.05) = 0.01
                            prior.range = c(5000, 0.05),
                            # PC-prior on sigma: P(sigma > 1) = 0.01
                            # n=725, 2022-12-01, sigma= 5-10, 0.05
                            # n=6774 2022-01-01, 
                            prior.sigma = c(5, 0.05)) 

# A.global <- inla.spde.make.A(mesh = mesh, loc = coords)
# s.index <- inla.spde.make.index(name = "spatial.field",
#                                 n.spde = spde$n.spde)



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

y.pp = rep(0:1, c(nv, n))
e.pp = c(w, rep(0, n))

# Projection matrix for the integration points (mesh vertices)
A.int <- Diagonal(nv, rep(1, nv))
# Projection matrix for observed points (event locations)
A.y <- inla.spde.make.A(mesh = mesh, loc = coo)
# Projection matrix for mesh vertices and event locations
A.pp <- rbind(A.int, A.y)
Ap.pp <- inla.spde.make.A(mesh = mesh, loc = coop)

# stack for estimation
stk.e.pp <- inla.stack(tag = "est.pp",
                       data = list(y = y.pp, e = e.pp), 
                       A = list(1, A.pp),
                       effects = list(list(b0 = rep(1, nv + n)), list(s = 1:nv)))


# stack for prediction stk.p
stk.p.pp <- inla.stack(tag = "pred.pp",
                       data = list(y = rep(NA, nrow(coop)), e = rep(0, nrow(coop))),
                       A = list(1, Ap.pp),
                       effects = list(data.frame(b0 = rep(1, nrow(coop))),
                                      list(s = 1:nv)))

stk.full.pp <- inla.stack(stk.e.pp, stk.p.pp)

formula <- y ~ 0 + b0 + f(s, model = spde)
res <- inla(formula,  family = 'poisson',
            data = inla.stack.data(stk.full.pp),
            control.predictor = list(compute = TRUE, link = 1,
                                     A = inla.stack.A(stk.full.pp)),
                                      E = inla.stack.data(stk.full.pp)$e,
            # control.inla=list(diagonal=100),
            verbose=TRUE)

index <- inla.stack.index(stk.full.pp, tag = "pred.pp")$data
pred_mean <- res$summary.fitted.values[index, "mean"]
pred_ll <- res$summary.fitted.values[index, "0.025quant"]
pred_ul <- res$summary.fitted.values[index, "0.975quant"]


inla.spde2.result(inla = res, name = "spatial.field",
                  spde = spde, do.transf = TRUE)

grid$mean <- NA
grid$ll <- NA
grid$ul <- NA

grid$mean[indicespointswithin] <- pred_mean
grid$ll[indicespointswithin] <- pred_ll
grid$ul[indicespointswithin] <- pred_ul


library(rasterVis)
levelplot(raster::brick(grid), layout = c(3, 1),
          names.attr = c("Mean", "2.5 percentile", "97.5 percentile"))

library(fields)
local.plot.field = function(field, mesh, xlim=c(-15,-4), ylim=c(35,43), ...){
  stopifnot(length(field) == mesh$n)
  proj = inla.mesh.projector(mesh, xlim = xlim, 
                             ylim = ylim, dims=c(300, 300))
  field.proj = inla.mesh.project(proj, field)
  n.col = 20
  image.plot(list(x = proj$x, y=proj$y, z = field.proj), 
             xlim = xlim, ylim = ylim, col = plasma(n.col), nlevel=n.col+1, ...)
}
local.plot.field(res$summary.random[['s']][['mean']], mesh)
len = res$summary.hyperpar[2, '0.5quant']
# - the posterior median range
arrows(5-0.5*len, 5, 5+0.5*len, 5, length=0.05, angle=90, code=3, lwd=3)



# Extract the estimated intensity values
estimated_intensity <- pp.res0$summary.random$spatial.field[, "mean"]

# Create a data frame for plotting
plot_data <- data.frame(coords = mesh$loc[, 1:2], intensity = estimated_intensity)

# Convert to a spatial object (if necessary) and then to an sf object for ggplot2 plotting
spatial_points <- SpatialPointsDataFrame(coords = plot_data[,1:2], data = plot_data)
plot_data_sf <- st_as_sf(spatial_points, coords = c("coords.X1", "coords.X2"), crs = 4326)

# Plot using ggplot2
ggplot(data = plot_data_sf) +
  geom_sf(aes(color = intensity, size = intensity), alpha = 0.6) +
  scale_color_viridis_c() + # Use a color scale that represents intensity well
  scale_size(range = c(1, 5)) + # Adjust size range based on your data
  theme_minimal() +
  labs(title = "Estimated Intensity of Log Gaussian Cox Model",
       color = "Intensity",
       size = "Intensity")



spde.est <- inla.spde2.result(inla = pp.res0, name = "spatial.field",
                                  spde = spde, do.transf = TRUE)
inla.zmarginal(spde.est$marginals.variance.nominal[[1]])
inla.zmarginal(spde.est$marginals.range.nominal[[1]])
