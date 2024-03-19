library("sf")
library("spocc")

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

df <- occ(query = "solanum", from = "gbif",
          date = c("2015-01-01", "2022-12-31"),
          gbifopts = list(country = "BO"),
          has_coords = TRUE, limit = 1000)
d <- occ2df(df)

d <- st_as_sf(d[, 2:3], coords = c("longitude", "latitude"))
st_crs(d) <- "EPSG:4326"

st_crs("EPSG:5356")$proj4string
projUTM <- "+proj=utm +zone=19 +south +ellps=GRS80
+towgs84=0,0,0,0,0,0,0 +units=km +no_defs"
d <- st_transform(d, crs = projUTM)

library(rnaturalearth)
map <- ne_countries(type = "countries", country = "Bolivia",
                    scale = "medium", returnclass = "sf")
map <- st_transform(map, crs = projUTM)

library("ggplot2")
ggplot() + geom_sf(data = map) +
  geom_sf(data = d) + coord_sf(datum = projUTM)

coo <- st_coordinates(d)

library(sf)
library(terra)

# raster grid covering map
grid <- terra::rast(map, nrows = 100, ncols = 100)
# coordinates of all cells
xy <- terra::xyFromCell(grid, 1:ncell(grid))

# transform points to a sf object
dp <- st_as_sf(as.data.frame(xy), coords = c("x", "y"),
               crs = st_crs(map))

# indices points within the map
indicespointswithin <- which(st_intersects(dp, map,
                                           sparse = FALSE))

# points within the map
dp <- st_filter(dp, map)

ggplot() + geom_sf(data = map) +
  geom_sf(data = dp) + coord_sf(datum = projUTM)

coop <- st_coordinates(dp)


library(INLA)

summary(dist(coo)) # summary of distances between event locations


loc.d <- cbind(st_coordinates(map)[, 1], st_coordinates(map)[, 2])

mesh <- inla.mesh.2d(loc.domain = loc.d, max.edge = c(50, 100),
                     offset = c(50, 100), cutoff = 1)

plot(mesh)
points(coo, col = "red")
axis(1)
axis(2)

nv <- mesh$n
n <- nrow(coo)
spde <- inla.spde2.matern(mesh = mesh, alpha = 2, constr = TRUE)

dmesh <- book.mesh.dual(mesh)
plot(dmesh)
axis(1)
axis(2)

# Domain polygon is converted into a SpatialPolygons
domain.polys <- Polygons(list(Polygon(loc.d)), '0')
domainSP <- SpatialPolygons(list(domain.polys))

# Because the mesh is larger than the study area, we need to
# compute the intersection between each polygon
# in the dual mesh and the study area
library(rgeos)

w <- sapply(1:length(dmesh), function(i) {
  if (gIntersects(dmesh[i, ], domainSP))
    return(gArea(gIntersection(dmesh[i, ], domainSP)))
  else return(0)
})

sum(w) # sum weights
st_area(map) # area of the study region

plot(mesh)
points(mesh$loc[which(w > 0), 1:2], col = "black", pch = 20)
points(mesh$loc[which(w == 0), 1:2], col = "red", pch = 20)

y.pp <- rep(0:1, c(nv, n))
e.pp <- c(w, rep(0, n))

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

# stk.full has stk.e and stk.p
stk.full.pp <- inla.stack(stk.e.pp, stk.p.pp)

formula <- y ~ 0 + b0 + f(s, model = spde)

res <- inla(formula,  family = 'poisson',
            data = inla.stack.data(stk.full.pp),
            control.predictor = list(compute = TRUE, link = 1,
                                     A = inla.stack.A(stk.full.pp)),
            E = inla.stack.data(stk.full.pp)$e)

index <- inla.stack.index(stk.full.pp, tag = "pred.pp")$data
pred_mean <- res$summary.fitted.values[index, "mean"]
pred_ll <- res$summary.fitted.values[index, "0.025quant"]
pred_ul <- res$summary.fitted.values[index, "0.975quant"]

grid$mean <- NA
grid$ll <- NA
grid$ul <- NA

grid$mean[indicespointswithin] <- pred_mean
grid$ll[indicespointswithin] <- pred_ll
grid$ul[indicespointswithin] <- pred_ul

library(rasterVis)
levelplot(raster::brick(grid), layout = c(3, 1),
          names.attr = c("Mean", "2.5 percentile", "97.5 percentile"))