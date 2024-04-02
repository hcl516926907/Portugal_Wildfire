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

# bru_safe_sp(force = TRUE)
load(file.path(dir.data, "burn_area","finaldata_urbanfire.RData"))

data.new$year.idx <- data.new$year -2011
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

ips <- fm_int(
  domain = list(coordinates = mesh, year = 1:9),
  samplers =domainSP
)



r <- 50
p1 <- 0.1
sigma <- 8
p2 <- 0.1

spde <- inla.spde2.pcmatern(mesh = mesh,
                            # PC-prior on range: P(practic.range < r) = p1
                            prior.range = c(r, p1),
                            # PC-prior on sigma: P(sigma > sigma) = p2
                            prior.sigma = c(sigma, p2))




# cmp <- coordinates ~  Intercept(1) +
#   mySmooth(coordinates , model = spde)
# 
# t1 <-  Sys.time()
# fit <- lgcp(cmp, data=coords, samplers = domainSP,
#             domain = list(coordinates  = mesh))
# t2 <-  Sys.time()
# print(t2-t1)

cmp <- coordinates + year.idx ~  Intercept(1) +
mySmooth(coordinates , model = spde, group = year.idx, ngroup = 9)

t1 <-  Sys.time()
fit <- lgcp(cmp, data=coords, samplers = domainSP,
             domain = list(coordinates  = mesh, year.idx = seq_len(9)),
             options=list(verbose=TRUE))
t2 <-  Sys.time()
print(t2-t1)



ppxl <- fm_pixels(mesh, mask = domainSP, format = "sp")
ppxl_all <- fm_cprod(ppxl, data.frame(year.idx = seq_len(9)))

lambda1 <- predict(
  fit,
  ppxl_all,
  ~ data.frame(year.idx = year.idx, lambda = exp(mySmooth + Intercept))
)

pl1 <- ggplot() +
  gg(as(lambda1, "SpatialPixelsDataFrame"), aes(fill = q0.5)) +
  # gg(loc.data.utm, size = 0.1) +
  facet_wrap(~year.idx) +
  coord_sf()
pl1


prepare_residual_calculations <- function(samplers, domain, observations, year) {
  observations <- observations[observations$year.idx==year,]
  
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
    data = bind_rows(data.frame(obs = rep(FALSE, domain$n), year.idx=year), observations@data),
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




residual_df <- function(model, df, expr, A_sum, A_integrate, year) {
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
  res$Scaling_Residuals$Type <- "Scaling Residuals"
  res$Inverse_Residuals$Type <- "Inverse Residuals"
  res$Pearson_Residuals$Type <- "Pearson Residuals"
  do.call(rbind, res)
}


B <- SpatialPolygonsDataFrame(domainSP, data.frame('weight'=1), match.ID = F) 

B1 <- partition(samplers = B, nrows = 20, ncols = 10)
plot(B1, main = "Grid partitions of B")

residual_df_temproal <- function(year, samplers, domain, observations, expr){
  As <- prepare_residual_calculations(
    samplers = samplers, domain = domain,
    observations = observations, year=year
  )
  
  res <- residual_df(
    fit, As$df, expr,
    As$A_sum, As$A_integrate, year=year
  )
  res$year.idx <- year
  return(res)
}



res <- do.call(rbind, lapply(1:9, residual_df_temproal, samplers=B1, domain=mesh, observations=coords,
                                                        expr=expression(exp(mySmooth + Intercept))))

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

csc <- set_csc(res, rep("RdBu", 3))
model_name <- 'spatial_temporal'
samplers.all <-  fm_cprod(B1, data.frame(year.idx = seq_len(9)))

residual_plot_temporal <- function(samplers.all, residuals, csc, model_name) {
  # Initialise the scaling residuals plot
  samplers.all$Residual <- res %>%
    filter(Type == "Scaling Residuals") %>%
    pull(mean)
  scaling <- ggplot() +
    gg(samplers.all, aes(fill = Residual), alpha = 1, colour = NA) +
    csc["Scaling"] +
    theme(legend.position = "bottom") +
    facet_wrap(~year.idx) +
    labs(subtitle = paste(model_name, "Scaling"))
  
  # Initialise the inverse residuals plot
  samplers.all$Residual <- residuals %>%
    filter(Type == "Inverse Residuals") %>%
    pull(mean)
  inverse <- ggplot() +
    gg(samplers.all, aes(fill = Residual), alpha = 1, colour = NA) +
    csc["Inverse"] +
    theme(legend.position = "bottom") +
    facet_wrap(~year.idx) +
    labs(subtitle = paste(model_name, "Inverse"))
  
  # Initialise the Pearson residuals plot
  samplers.all$Residual <- residuals %>%
    filter(Type == "Pearson Residuals") %>%
    pull(mean)
  pearson <- ggplot() +
    gg(samplers.all, aes(fill = Residual), alpha = 1, colour = NA) +
    csc["Pearson"] +
    theme(legend.position = "bottom") +
    facet_wrap(~year.idx) +
    labs(subtitle = paste(model_name, "Pearson"))
  
  # Return the three plots in a list
  list(
    Scaling = scaling, Inverse = inverse,
    Pearson = pearson
  )
  
}




fit_csc <- set_csc(res, rep("RdBu", 3))
# Store plots
plotB1 <- residual_plot_temporal(samplers.all, res, fit_csc, "SPDE Model")



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
