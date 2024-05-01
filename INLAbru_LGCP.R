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
load(file.path(dir.data, "burn_area","finaldata_ruralfire_fwi.RData"))
load(file.path(dir.data, "burn_area","weather_covariates.RData"))
data.new$month <- as.integer(format(data.new$open,format='%m'))
data.new <- data.new[data.new$year==2020 & data.new$month==9,]

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


inner.length <- 40
outer.length <- 70
cutoff <- 15


mesh <- inla.mesh.2d(loc.domain = loc.d, max.edge=c(inner.length,outer.length), cutoff=cutoff,
                     crs=projutm)
mesh$n
plot(mesh)
# points(coords, col = 4, pch = 19, cex=0.2)

inner.length <- 20
outer.length <- 50
cutoff <- 6

mesh2 <- fm_mesh_2d_inla(loc.domain = loc.d,
                      max.edge=c(inner.length,outer.length),
                      cutoff=cutoff,
                     crs=projutm)
mesh2$n
plot(mesh2)

domain.polys <- Polygons(list(Polygon(loc.d)), '0')
domainSP <- SpatialPolygons(list(domain.polys),proj4string=CRS(projutm))


#no of mesh points...
nv <- mesh$n 



# r <- 50
# p1 <- 0.1
# sigma <- 8
# p2 <- 0.1

r <- 100
p1 <- 0.01
sigma <- 2
p2 <- 0.01

spde <- inla.spde2.pcmatern(mesh = mesh,
                            # PC-prior on range: P(practic.range < r) = p1
                            prior.range = c(r, p1),
                            # PC-prior on sigma: P(sigma > sigma) = p2
                            prior.sigma = c(sigma, p2))

spde2 <- inla.spde2.pcmatern(mesh = mesh2,
                            # PC-prior on range: P(practic.range < r) = p1
                            prior.range = c(r, p1),
                            # PC-prior on sigma: P(sigma > sigma) = p2
                            prior.sigma = c(sigma, p2))


mesh1D.LVegCov <- fm_mesh_1d(seq(0,1,by=0.05), boundary = "free")
spde1D.LVegCov <- inla.spde2.pcmatern(mesh1D.LVegCov,constr = TRUE,
                                   prior.range = c(0.001, 0.01),
                                   prior.sigma = c(0.02, 0.01)
)




# cmp <- coordinates ~  Intercept(1) +
#   mySmooth(coordinates , model = spde)
# 
# t1 <-  Sys.time()
# fit <- lgcp(cmp, data=coords, samplers = domainSP,
#             domain = list(coordinates  = mesh))
# t2 <-  Sys.time()
# print(t2-t1)

cmp1 <- coordinates  ~  Intercept(1) +
  mySmooth(coordinates , model = spde)
t1 <-  Sys.time()
fit1 <- lgcp(cmp1, data=coords, samplers = domainSP,
            domain = list(coordinates  = mesh),
            options=list(verbose=TRUE))
t2 <-  Sys.time()
print(t2-t1)


temproal_spdf_LVegTyp <- function(time){
  m <- as.integer(time)
  # average over the time.idx
  spdf <- SpatialPixelsDataFrame(points = grid_pixels,
                                 data = data.frame(var=as.factor(LVegTyp.month[,,m]),
                                                   time=time))
  proj4string(spdf) <- CRS("EPSG:4326")
  return(spdf)
}

f.LVegTyp.utm <- function(where) {
  x <- where@coords[,1]
  y <- where@coords[,2]
  spp <- SpatialPoints(data.frame(x=x,y=y),
                       proj4string=CRS('+proj=utm +zone=29 +datum=WGS84 +units=km +no_defs' ))
  spp <- spTransform(spp, CRS("EPSG:4326"))
  v <- rep(NA, nrow(where))
    temp_spdf <- temproal_spdf_LVegTyp(105)
    v<- over(spp,temp_spdf[,'var'])$var
    if (any(is.na(v))) {
      v <- bru_fill_missing(temp_spdf, spp, v)
    }
  return(v)
}

f.LVegTyp.utm(coords)

temproal_spdf_LVegCov <- function(time){
  m <- as.integer(time)
  # average over the time.idx
  spdf <- SpatialPixelsDataFrame(points = grid_pixels,
                                 data = data.frame(var=as.vector(LVegCov.month[,,m]),
                                                   time=time))
  proj4string(spdf) <- CRS("EPSG:4326")
  return(spdf)
}

f.LVegCov.utm <- function(where) {
  x <- where@coords[,1]
  y <- where@coords[,2]
  spp <- SpatialPoints(data.frame(x=x,y=y),
                       proj4string=CRS('+proj=utm +zone=29 +datum=WGS84 +units=km +no_defs' ))
  spp <- spTransform(spp, CRS("EPSG:4326"))
  v <- rep(NA, nrow(where))
  temp_spdf <- temproal_spdf_LVegCov(105)
  v<- over(spp,temp_spdf[,'var'])$var
  if (any(is.na(v))) {
    v <- bru_fill_missing(temp_spdf, spp, v)
  }
  return(v)
}


f.LVegTyp.utm_Combined <- function(where) {
  x <- where@coords[,1]
  y <- where@coords[,2]
  spp <- SpatialPoints(data.frame(x=x,y=y),
                       proj4string=CRS('+proj=utm +zone=29 +datum=WGS84 +units=km +no_defs' ))
  spp <- spTransform(spp, CRS("EPSG:4326"))
  v <- rep(NA, nrow(where))
  temp_spdf <- temproal_spdf_LVegTyp(105)
  v<- over(spp,temp_spdf[,'var'])$var
  v[which(v!='7')] <- '0'
  if (any(is.na(v))) {
    v <- bru_fill_missing(temp_spdf, spp, v)
    v[which(v!='7')] <- '0'
  }
  return(v)
}


f.LVegTyp.utm_Combined(coords)
popu_2015$Population_Rate <- (popu_2015$Population_Rate - mean(popu_2015$Population_Rate, na.rm = TRUE))/sd(popu_2015$Population_Rate)

f.Popu.utm <- function(where) {
  x <- where@coords[,1]
  y <- where@coords[,2]
  spp <- SpatialPoints(data.frame(x=x,y=y),
                       proj4string=CRS('+proj=utm +zone=29 +datum=WGS84 +units=km +no_defs' ))
  spp <- spTransform(spp, CRS("+proj=longlat +datum=WGS84 +no_defs"))
  v <- over(spp,popu_2015[,'Population_Rate'])$Population_Rate
  if (any(is.na(v))) {
    v <- bru_fill_missing(popu_2015, spp, v)
  }
  return(v)
}
summary(f.Popu.utm(coords))

cmp2.1 <- coordinates  ~  
  mySmooth(coordinates , model = spde)+
  LVegCov(f.LVegCov.utm(.data.),model='linear')+
  # Popu(f.Popu.utm(.data.),model='linear') +
  Intercept(1)
  # LVegTyp(f.LVegTyp.utm_Combined(.data.),model='factor_full') -1 

t1 <-  Sys.time()

fit2.1 <- lgcp(cmp2.1, data=coords, samplers = domainSP,
             domain = list(coordinates  = mesh),
             options=list(verbose=TRUE))
t2 <-  Sys.time()
print(t2-t1)
print(summary(fit2.1))


cmp2.2 <- coordinates  ~  
  mySmooth(coordinates , model = spde)+
  LVegCov(f.LVegCov.utm(.data.),model=spde1D.LVegCov)+
  # Popu(f.Popu.utm(.data.),model='linear') +
  Intercept(1)
# LVegTyp(f.LVegTyp.utm_Combined(.data.),model='factor_full') -1 

t1 <-  Sys.time()

fit2.2 <- lgcp(cmp2.2, data=coords, samplers = domainSP,
               domain = list(coordinates  = mesh),
               options=list(verbose=TRUE))
t2 <-  Sys.time()
print(t2-t1)
print(summary(fit2.2))

cmp2.3 <- coordinates  ~  
  mySmooth(coordinates , model = spde)+
  # LVegCov(f.LVegCov.utm(.data.),model=spde1D.LVegCov)+
  Popu(f.Popu.utm(.data.),model='linear') +
  Intercept(1)
# LVegTyp(f.LVegTyp.utm_Combined(.data.),model='factor_full') -1 

t1 <-  Sys.time()

fit2.3 <- lgcp(cmp2.3, data=coords, samplers = domainSP,
               domain = list(coordinates  = mesh),
               options=list(verbose=TRUE))
t2 <-  Sys.time()
print(t2-t1)
print(summary(fit2.3))

cmp2.4 <- coordinates  ~  
  mySmooth(coordinates , model = spde)+
  # LVegCov(f.LVegCov.utm(.data.),model=spde1D.LVegCov)+
  # Popu(f.Popu.utm(.data.),model='linear') +
  # Intercept(1)
LVegTyp(f.LVegTyp.utm_Combined(.data.),model='factor_full') -1

t1 <-  Sys.time()

fit2.4 <- lgcp(cmp2.4, data=coords, samplers = domainSP,
               domain = list(coordinates  = mesh),
               options=list(verbose=TRUE))
t2 <-  Sys.time()
print(t2-t1)
print(summary(fit2.4))

cmp3 <- coordinates  ~  Intercept(1) +
  mySmooth(coordinates, model = spde2) 

t1 <-  Sys.time()
fit3 <- lgcp(cmp3, data=coords, samplers = domainSP,
             domain = list(coordinates  = mesh2),
             options=list(verbose=TRUE))
t2 <-  Sys.time()
print(t2-t1)

cmp4 <- coordinates  ~
  mySmooth(coordinates, model = spde2) + 
  # LVegTyp(f.LVegTyp.utm(.data.),model='factor_full') -1
  LVegCov(f.LVegCov.utm(.data.),model='linear') + Intercept(1)

t1 <-  Sys.time()
fit4 <- lgcp(cmp4, data=coords, samplers = domainSP,
             domain = list(coordinates  = mesh2),
             options=list(verbose=TRUE))
t2 <-  Sys.time()
print(t2-t1)



new_data <- fm_pixels(mesh, mask = domainSP, format = "sp",dims=c(20,50))
new_data$open <- rep(data.new[1,]$open,nrow(new_data))
new_data$LVegTyp <- f.LVegTyp.utm(new_data)
new_data$LVegCov <-  f.LVegCov.utm(new_data)
new_data$LVegTyp_combined <-  f.LVegTyp.utm_Combined(new_data)
new_data$Populaiton <- f.Popu.utm(new_data)

ggplot() + gg(new_data,aes(fill=Populaiton)) 

lprange <- range(new_data$LVegCov)
csc <- scale_fill_gradientn(colours = brewer.pal(9, "YlOrRd"), limits = lprange)
ggplot() +
  # gg(mesh) +
  gg(new_data,aes(fill=LVegCov)) +
  # geom_sf(data = boundary, alpha = 0.1, fill = "blue") +
  geom_sf(data = st_as_sf(coords)) + 
  ggtitle('Low Vegetation Coverate')+
  csc


pred1<- predict(
  fit1,
  new_data,
  ~ data.frame(
    lambda = exp(mySmooth + Intercept)
  )
)

pred2<- predict(
  fit2,
  new_data,
  ~ data.frame(
    lambda = exp(mySmooth  + LVegCov + Intercept)
  )
)

pred3<- predict(
  fit3,
  new_data,
  ~ data.frame(
    lambda = exp(mySmooth   + Intercept)
  )
)

pred4<- predict(
  fit4,
  new_data,
  ~ data.frame(
    lambda = exp(mySmooth  + LVegTyp)
  )
)

# new_data$temp <- f.temp.utm(new_data)
# pred3<- predict(
#   fit3,
#   new_data,
#   ~ data.frame(
#     lambda = exp(mySmooth + Intercept + temp)
#   )
# )

lprange <- range(pred1$mean,pred2$mean,pred3$mean)
csc <- scale_fill_gradientn(colours = brewer.pal(9, "YlOrRd"), limits = lprange)
pl1 <- ggplot() +
  gg(pred1,aes(fill=mean)) +
  csc+
  gg(domainSP) +
  ggtitle("LGCP fit to Points", subtitle = "Only SPDE")

pl2 <- ggplot() +
  gg(pred2, aes(fill=mean)) +
  csc+
  gg(domainSP, alpha = 0) +
  ggtitle("LGCP fit to Points", subtitle = "SPDE + Covariates")

multiplot(pl1,pl2,cols=2)


pl3 <- ggplot() +
  gg(pred3, aes(fill=mean))  +
  csc+
  gg(domainSP, alpha = 0) +
  ggtitle("LGCP fit to Points", subtitle = "Only Covariates")

multiplot(pl1,pl2, pl3,cols = 3)


ggplot() +
  gg(cbind(pred$lambda, data.frame(property = "q0.500")), aes(fill = median)) +
  gg(cbind(pred$lambda, data.frame(property = "q0.025")), aes(fill = q0.025)) +
  gg(cbind(pred$lambda, data.frame(property = "q0.975")), aes(fill = q0.975)) +
  coord_equal() +
  facet_wrap(~property)


int.plot <- plot(fit, "Intercept")
spde.range <- spde.posterior(fit, "mySmooth", what = "range")
spde.logvar <- spde.posterior(fit, "mySmooth", what = "log.variance")
range.plot <- plot(spde.range)
var.plot <- plot(spde.logvar)
# 
multiplot(range.plot, var.plot, int.plot)
# 

prepare_residual_calculations <- function(samplers, domain, observations) {
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
    data = bind_rows(data.frame(obs = rep(FALSE, domain$n)), observations@data),
    proj4string = fm_CRS(domain)
  )
  
  # Return A-sum, A_integrate and the data frame for predicting the residuals
  list(A_sum = A_sum, A_integrate = A_integrate, df = df)
}




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



#' ---------------
#' residual_plot
#' ---------------
#'
#' plots the three types of residuals for each polygon
#'
#' Input:
#' @param samplers A SpatialPolygonsDataFrame containing partitions for which
#' residuals are to be calculated
#' @param residuals frame containing residual information for each of the
#' partitions of the subset 'B'
#' @param csc list of three colour scales for the three types of residuals
#' @param model_name string containing the name of the model being assessed
#'
#' Output:
#' @return a list of three subplots scaling, inverse and Pearson residuals
#' for the different partitions of samplers


residual_plot <- function(samplers, residuals, csc, model_name) {
  # Initialise the scaling residuals plot
  samplers$Residual <- residuals %>%
    filter(Type == "Scaling Residuals") %>%
    pull(mean)
  scaling <- ggplot() +
    gg(samplers, aes(fill = Residual), alpha = 1, colour = NA) +
    csc["Scaling"] +
    theme(legend.position = "bottom") +
    labs(subtitle = paste(model_name, "Scaling"))
  
  # Initialise the inverse residuals plot
  samplers$Residual <- residuals %>%
    filter(Type == "Inverse Residuals") %>%
    pull(mean)
  inverse <- ggplot() +
    gg(samplers, aes(fill = Residual), alpha = 1, colour = NA) +
    csc["Inverse"] +
    theme(legend.position = "bottom") +
    labs(subtitle = paste(model_name, "Inverse"))
  
  # Initialise the Pearson residuals plot
  samplers$Residual <- residuals %>%
    filter(Type == "Pearson Residuals") %>%
    pull(mean)
  pearson <- ggplot() +
    gg(samplers, aes(fill = Residual), alpha = 1, colour = NA) +
    csc["Pearson"] +
    theme(legend.position = "bottom") +
    labs(subtitle = paste(model_name, "Pearson"))
  
  # Return the three plots in a list
  list(
    Scaling = scaling, Inverse = inverse,
    Pearson = pearson
  )
}




#' ------------------
#' partition
#' ------------------
#'
#' Partitions the region based on the given criteria for calculating residuals
#' in each partition. Parts of this function are taken from concepts in
#' https://rpubs.com/huanfaChen/grid_from_polygon
#'
#' Input:
#' @param samplers A SpatialPolygonsDataFrame containing region for which
#' partitions need to be created
#' @param resolution resolution of the grids that are required
#' @param nrows number of rows of grids that are required
#' @param ncols number of columns of grids that are required
#'
#' Output:
#' @return a partitioned SpatialPolygonsDataFrame as required
#'
#'
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


#' ------------------
#' edit_df
#' ------------------
#'
#' Edits the residual data frames to remove columns that need not be displayed
#'
#' Input:
#' @param df the data frame which needs to be edited
#' @param columns a vector of columns that need to be deleted from df
#'
#' Output:
#' @return the edited data frame with only the desired columns
#'
#'
edit_df <- function(df, columns) {
  # Remove the columns that are not required
  df[, !(colnames(df) %in% columns)]
}

B <- SpatialPolygonsDataFrame(domainSP, data.frame('weight'=1), match.ID = F) 

As <- prepare_residual_calculations(
  samplers = B, domain = mesh,
  observations = coords
)

res <- residual_df(
  fit, As$df, expression(exp(mySmooth + Intercept)),
  As$A_sum, As$A_integrate
)
knitr::kable(edit_df(res, c("Type", "mean.mc_std_err", "sd.mc_std_err", "median")))



B1 <- partition(samplers = B, nrows = 20, ncols = 10)
plot(B1, main = "Grid partitions of B")


As1 <- prepare_residual_calculations(
  samplers = B1, domain = mesh,
  observations = coords
)
# Residuals for the vegetation model
res1 <- residual_df(
  fit1, As1$df, expression(exp(mySmooth + Intercept)),
  As1$A_sum, As1$A_integrate
)
knitr::kable(edit_df(res1, c(
  "Type", "mean.mc_std_err",
  "sd.mc_std_err", "median"
)))

grd <- SpatialGrid(GridTopology(cellcentre.offset = c(-10, 36),
                                cellsize = c(0.25,0.25), 
                                cells.dim = c(17, 29)))
# sgrd <- SpatialGridDataFrame(grd, data = data.frame(val = runif(240)), proj4string = CRS("+proj=longlat +datum=WGS84"))
grd_sp <- SpatialPixelsDataFrame(points = grd, data = data.frame(id = 1:length(grd)),proj4string = CRS("+proj=longlat +datum=WGS84"))

grd_poly <- as(grd_sp, 'SpatialPolygons')
grd_poly <- spTransform(grd_poly, CRS(projutm))


B2 <- as_Spatial(st_intersection(st_as_sf(grd_poly),st_as_sf(B)))

As2 <- prepare_residual_calculations(
  samplers = B2, domain = mesh,
  observations = coords
)
res2.1 <- residual_df(
  fit1, As2$df, expression(exp(mySmooth + Intercept)),
  As2$A_sum, As2$A_integrate
)
res2.2.1 <- residual_df(
  fit2.1, As2$df, expression(exp(mySmooth + LVegCov+Intercept)),
  As2$A_sum, As2$A_integrate
)

res2.2.2 <- residual_df(
  fit2.2, As2$df, expression(exp(mySmooth + LVegCov+Intercept)),
  As2$A_sum, As2$A_integrate
)

res2.2.3 <- residual_df(
  fit2.3, As2$df, expression(exp(mySmooth + Popu+Intercept)),
  As2$A_sum, As2$A_integrate
)

res2.2.4 <- residual_df(
  fit2.4, As2$df, expression(exp(mySmooth + LVegTyp)),
  As2$A_sum, As2$A_integrate
)

As2.3 <- prepare_residual_calculations(
  samplers = B2, domain = mesh2,
  observations = coords
)
res3 <- residual_df(
  fit3, As2.3$df, expression(exp(mySmooth + Intercept)),
  As2.3$A_sum, As2.3$A_integrate
)



fit_csc1 <- set_csc(res2.2.2, rep("RdBu", 3))
# Store plots
plotB2.1 <- residual_plot(B2, res2.1, fit_csc1, "fit1")
plotB2.2.1 <- residual_plot(B2, res2.2.1, fit_csc1, "fit2.1")
plotB2.2.2 <- residual_plot(B2, res2.2.2, fit_csc1, "fit2.2")
plotB2.2.3 <- residual_plot(B2, res2.2.3, fit_csc1, "fit2.3")
plotB2.2.4 <- residual_plot(B2, res2.2.4, fit_csc1, "fit2.4")
plotB2.3 <- residual_plot(B2, res2.3, fit_csc1, "fit3")

((plotB2.1$Pearson | plotB2.2.1$Pearson | plotB2.2.2$Pearson) ) +
  plot_annotation(title = "Vegetation Model") +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

((plotB2.2.3$Scaling | plotB2.2.3$Pearson ) )+
  plot_annotation(title = "Vegetation Model") +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

((plotB2.2.3$Pearson | plotB2.2.4$Pearson | plotB2.3$Pearson) ) +
  plot_annotation(title = "Vegetation Model") +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")



((plotB1$Scaling | plotB1$Inverse | plotB1$Pearson) ) +
  plot_annotation(title = "Vegetation Model") +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

((plotB2$Scaling | plotB2$Inverse | plotB2$Pearson) ) +
  plot_annotation(title = "Vegetation Model") +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")


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
district.utm <- spTransform(district, CRS(projutm))

As2 <- prepare_residual_calculations(
  samplers = district.utm, domain = mesh,
  observations = coords
)
# Residuals for the vegetation model
res2 <- residual_df(
  fit, As2$df, expression(exp(mySmooth + Intercept)),
  As2$A_sum, As2$A_integrate
)
knitr::kable(edit_df(res2, c(
  "Type", "mean.mc_std_err",
  "sd.mc_std_err", "median"
)))

fit_csc <- set_csc(res2, rep("RdBu", 3))
# Store plots
plot.district <- residual_plot(district.utm, res2, fit_csc, "SPDE Model")

((plot.district$Scaling | plot.district$Inverse | plot.district$Pearson) ) +
  plot_annotation(title = "Vegetation Model") +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")
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
