library("spatstat")
data(clmfires)

clmfires0407 <- clmfires[clmfires$marks$date >= "2004-01-01"]

# Set `urban` instead of `artifgreen`
idx <- which(clmfires.extra$clmcov100$landuse$v %in% c("artifgreen",
                                                       "farm"))
clmfires.extra$clmcov100$landuse$v[idx] <- "urban"

# Convert to factor
clmfires.extra$clmcov100$landuse$v <- 
  factor(as.character(clmfires.extra$clmcov100$landuse$v))
# Set right dimension of raster object

dim(clmfires.extra$clmcov100$landuse$v) <- c(100, 100)

clmfires.extra$clmcov100$landuse

#In addition, we will rescale `elevation` (to express it in kilometers) and
#`orientation` (to be in radians) so that fixed effects are better estimated:
clmfires.extra$clmcov100$elevation <- 
  clmfires.extra$clmcov100$elevation / 1000
clmfires.extra$clmcov100$orientation <- 
  clmfires.extra$clmcov100$orientation * pi / 180 


clmfires0407 <- clmfires0407[clmfires0407$marks$cause == "lightning", ]


clm.bdy <- do.call(cbind, clmfires0407$window$bdry[[1]])
#Define mesh
clm.mesh <- inla.mesh.2d(loc.domain = clm.bdy, max.edge = c(15, 50),
                         offset = c(10, 10))

#Points
clm.pts <- as.matrix(coords(clmfires0407))
clm.mesh.pts <- as.matrix(clm.mesh$loc[, 1:2])
allpts <- rbind(clm.mesh.pts, clm.pts)

# Number of vertices in the mesh
nv <- clm.mesh$n
# Number of points in the data
n <- nrow(clm.pts)

#Create SPDE
clm.spde <- inla.spde2.pcmatern(mesh = clm.mesh, alpha = 2,
                                prior.range = c(50, 0.9), # P(range < 50) = 0.9
                                prior.sigma = c(1, 0.01) # P(sigma > 10) = 0.01
)

library("deldir")
library("SDraw")

# Voronoi polygons (as SpatialPolygons)
mytiles <- voronoi.polygons(SpatialPoints(clm.mesh$loc[, 1:2]))

# C-LM bounday as SpatialPolygons
clmbdy.sp <- SpatialPolygons(list(Polygons(list(Polygon (clm.bdy)),
                                           ID = "1"))) 

#Compute weights
require(rgeos)

w <- sapply(1:length(mytiles), function(p) {
  aux <- mytiles[p, ]  
  
  if(gIntersects(aux, clmbdy.sp) ) {
    return(gArea(gIntersection(aux, clmbdy.sp)))
  } else {
    return(0)
  }
})

# Sum of weights
sum(w)


# Area of study region
gArea(clmbdy.sp)

plot(clmbdy.sp)


#Prepare data
y.pp = rep(0:1, c(nv, n))
e.pp = c(w, rep(0, n))

lmat <- inla.spde.make.A(clm.mesh, clm.pts)
imat <- Diagonal(nv, rep(1, nv))

A.pp <-rbind(imat, lmat)

clm.spde.index <- inla.spde.make.index(name = "spatial.field",
                                       n.spde = clm.spde$n.spde)



#Covariates
allpts.ppp <- ppp(allpts[, 1], allpts[, 2], owin(xrange = c(-15.87, 411.38), 
                                                 yrange = c(-1.44, 405.19)))

# Assign values of covariates to points using value of nearest pixel
covs100 <- lapply(clmfires.extra$clmcov100, function(X){
  pixels <- nearest.pixel(allpts.ppp$x, allpts.ppp$y, X)
  sapply(1:npoints(allpts.ppp), function(i) {
    X[pixels$row[i], pixels$col[i]]
  })
})

# Intercept for spatial model
covs100$b0 <- rep(1, nv + n)



#Create data structure
clm.stack <- inla.stack(data = list(y = y.pp, e = e.pp),
                        A = list(A.pp, 1), 
                        effects = list(clm.spde.index, covs100),
                        tag = "pp")


#Data structure for prediction
library("maptools")
sgdf <- as(clmfires.extra$clmcov100$elevation, "SpatialGridDataFrame")
sp.bdy <- as(clmfires$window, "SpatialPolygons")
idx <- over(sgdf, sp.bdy)
spdf <- as(sgdf[!is.na(idx), ], "SpatialPixelsDataFrame")

pts.pred <- coordinates(spdf)
n.pred <- nrow(pts.pred)

#Get covariates (using subsetting operator in spatstat)
ppp.pred <- ppp(pts.pred[, 1], pts.pred[, 2], window = clmfires0407$window)
covs100.pred <- lapply(clmfires.extra$clmcov100, function(X) {
  X[ppp.pred]
})

# Intercept for spatial model
covs100.pred$b0 <- rep(1, n.pred)

#Prediction points
A.pred <- inla.spde.make.A (mesh = clm.mesh, loc = pts.pred)
clm.stack.pred <- inla.stack(data = list(y = NA),
                             A = list(A.pred, 1),
                             effects = list(clm.spde.index, covs100.pred), 
                             tag = "pred")


#Join data
join.stack <- inla.stack(clm.stack, clm.stack.pred)

pp.res0 <- inla(y ~ 1 + 
                  f(spatial.field, model = clm.spde), 
                family = "poisson", data = inla.stack.data(join.stack),
                control.predictor = list(A = inla.stack.A(join.stack), compute = TRUE,
                                         link = 1),
                control.inla = list(int.strategy = "eb"),
                E = inla.stack.data(join.stack)$e)


pp.res <- inla(y ~ 1 + landuse + elevation + orientation + slope +
                 f(spatial.field, model = clm.spde), 
               family = "poisson", data = inla.stack.data(join.stack),
               control.predictor = list(A = inla.stack.A(join.stack), compute = TRUE,
                                        link = 1), verbose = TRUE,
               control.inla = list(int.strategy = "eb"),
               control.fixed = list(expand.factor.strategy="inla"),
               E = inla.stack.data(join.stack)$e)

summary(pp.res)

#Prediction
idx <- inla.stack.index(join.stack, 'pred')$data
# MOdel with no covariates
spdf$SPDE0 <- pp.res0$summary.fitted.values[idx, "mean"]
#Model with covariates
spdf$SPDE <- pp.res$summary.fitted.values[idx, "mean"]


#Compute statistics in terms or range and variance
spde.clm.est <- inla.spde2.result(inla = pp.res, name = "spatial.field",
                                  spde = clm.spde, do.transf = TRUE)

#Variance
inla.zmarginal(spde.clm.est$marginals.variance.nominal[[1]])

#Range
inla.zmarginal(spde.clm.est$marginals.range.nominal[[1]])
