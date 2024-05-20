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


domain.polys <- Polygons(list(Polygon(loc.d)), '0')
domainSP <- SpatialPolygons(list(domain.polys),proj4string=CRS(projutm))



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


cmp <- coordinates + time.idx ~   Intercept(1) +
  mySmooth(coordinates , model = spde, group = time.idx, ngroup = 9)


for (t in c(1,seq(3,30,by=3))){
  cond <- data.new$month==9 & data.new$length > t*60
  coords <- SpatialPointsDataFrame(data.new[cond,],coords=data.new[cond,c('x_utm_new','y_utm_new')], 
                                   proj4string=CRS('+proj=utm +zone=29 +datum=WGS84 +units=km +no_defs' ))
  coords$time.idx <- coords$year.idx
  
  t1 <-  Sys.time()
  fit <- lgcp(cmp, data=coords, samplers = domainSP,
                 domain = list(coordinates  = mesh, time.idx = seq_len(9)),
                 options=list(verbose=TRUE))
  t2 <-  Sys.time()
  print(t2-t1)
  print(t)
  summary(fit)
  
  
  residuals <- do.call(rbind, lapply(1:9, residual_df_temproal, samplers=B2, domain=mesh, observations=coords,
                                  model=fit,
                                  expr=expression(exp(Intercept + mySmooth))))
  
  result <- data.frame(year.idx = rep(1:9,each=192))
  result$Actual_Cnt <- residuals[residuals$Type=='Scaling Residuals','mean'] + residuals[residuals$Type=='Estimated Lambda','mean'] 
  result$Scaling_Residuals <- residuals[residuals$Type=='Scaling Residuals','mean']
  result$Pearson_Residuals <- residuals[residuals$Type=='Pearson Residuals','mean']
  result$cutoff <- t
  
  if (t==1){
    result.all <- result
  }else{
    result.all <- rbind(result.all, result)
  }

}

summary(result.all[result.all$cutoff==2,'Scaling_Residuals'])

summary(result.all[result.all$cutoff==1,'Actual_Cnt'])


result.all$cutoff <- as.factor(result.all$cutoff)
ggplot(result.all, aes(x=cutoff,y=Actual_Cnt,group=cutoff)) + 
  geom_boxplot()+
  scale_x_discrete(limits=factor(c(1,seq(3,30,by=3))) , labels=factor(c(1,seq(3,30,by=3))))

ggplot(result.all, aes(x=cutoff,y=Scaling_Residuals,group=cutoff)) + 
  geom_boxplot()+
  scale_x_discrete(limits=factor(c(1,seq(3,30,by=3))) , labels=factor(c(1,seq(3,30,by=3))))

ggplot(result.all, aes(x=cutoff,y=Pearson_Residuals,group=cutoff)) + 
  geom_boxplot()+
  scale_x_discrete(limits=factor(c(1,seq(3,30,by=3))) , labels=factor(c(1,seq(3,30,by=3))))

result.all.sorted <- result.all %>% group_by(cutoff) %>%
  arrange(year.idx, Actual_Cnt, .by_group = TRUE)

result.all.sorted$plot.idx <- rep(1:1728, 11) 
result.all.sorted$Estimated_Lambda <- result.all.sorted$Actual_Cnt - result.all.sorted$Scaling_Residuals

cf <- 1
ggplot(data=result.all.sorted[result.all.sorted$cutoff==cf,] )+
  geom_line(aes(x=plot.idx, y=Estimated_Lambda)) + 
  geom_line(aes(x=plot.idx, y=Actual_Cnt), col='red') + 
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
  

plot(1:1728,result.all.sorted[result.all.sorted$cutoff==3,]$Actual_Cnt,type='l',ylim=c(-10,40))
lines(1:1728,result.all.sorted[result.all.sorted$cutoff==3,]$Scaling_Residuals,col='red')





samplers.all <-  do.call(rbind, lapply(1:9, function(x) {
  B2@data$time.idx=x 
  return(B2)
}))

samplers.all$Lambda <- residuals %>%
  filter(Type == "Estimated Lambda" ) %>%
  pull(mean)

samplers.all$Scaling_Residual <- residuals %>%
  filter(Type == "Scaling Residuals" ) %>%
  pull(mean)

my_residual_plot_temporal <- function(samplers.all, residuals, model_name) {
  samplers.all$Lambda <- residuals %>%
    filter(Type == "Estimated Lambda" ) %>%
    pull(mean)
  
  samplers.all$Scaling_Residual <- residuals %>%
    filter(Type == "Scaling Residuals" ) %>%
    pull(mean)
  
  samplers.all$Pearson_Residuals <- residuals %>%
    filter(Type == "Pearson Residuals" ) %>%
    pull(mean)
  
  scaling <- ggplot() +
    gg(samplers.all, aes(fill = Scaling_Residual), alpha = 1, colour = NA) +
    # csc["Scaling"] +
    theme(legend.position = "bottom") +
    facet_wrap(~time.idx) +
    labs(subtitle = paste(model_name, "Scaling"))
  
  pearson <- ggplot() +
    gg(samplers.all, aes(fill = Pearson_Residuals), alpha = 1, colour = NA) +
    # csc["Scaling"] +
    theme(legend.position = "bottom") +
    facet_wrap(~time.idx) +
    labs(subtitle = paste(model_name, "Pearson"))
  
  origin_cnt <- ggplot() +
    gg(samplers.all, aes(fill = Actual_Cnt), alpha = 1, colour = NA) +
    # csc["Actual_Cnt"] +
    theme(legend.position = "bottom") +
    facet_wrap(~time.idx) +
    labs(subtitle = paste(model_name, "Actual Cnt"))
  list(
    Scaling = scaling, Origin_cnt=origin_cnt,
    Pearson = pearson
  )
}

res_LGCP_24 <- residuals
res_LGCP_1 <- residuals
plot1 <- my_residual_plot_temporal(samplers.all, res_LGCP_1, model_name='Scaling Residuals of SPDE model')

###################################################
#appended functions, need to run them first
###################################################

prepare_residual_calculations <- function(samplers, domain, observations, time.idx) {
  observations <- observations[observations$time.idx==time.idx,]
  
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
    data = bind_rows(data.frame(obs = rep(FALSE, domain$n), time.idx=time.idx), observations@data),
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



# 
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
        Estimate_Lambda = as.vector(A_integrate %*% (h1 * lambda)[!obs]),
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
  res$Estimate_Lambda$Type <- "Estimated Lambda"
  res$Scaling_Residuals$Type <- "Scaling Residuals"
  res$Inverse_Residuals$Type <- "Inverse Residuals"
  res$Pearson_Residuals$Type <- "Pearson Residuals"
  do.call(rbind, res)
}


B <- SpatialPolygonsDataFrame(domainSP, data.frame('weight'=1), match.ID = F) 




residual_df_temproal <- function(time.idx, samplers, domain, observations, model, expr){
  As <- prepare_residual_calculations(
    samplers = samplers, domain = domain,
    observations = observations, time.idx=time.idx
  )
  
  res <- residual_df(
    model, As$df, expr,
    As$A_sum, As$A_integrate
  )
  res$time.idx <- time.idx
  return(res)
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




grd <- SpatialGrid(GridTopology(cellcentre.offset = c(-10, 36),
                                cellsize = c(0.25,0.25), 
                                cells.dim = c(17, 29)))
# sgrd <- SpatialGridDataFrame(grd, data = data.frame(val = runif(240)), proj4string = CRS("+proj=longlat +datum=WGS84"))
grd_sp <- SpatialPixelsDataFrame(points = grd, data = data.frame(id = 1:length(grd)),proj4string = CRS("+proj=longlat +datum=WGS84"))

grd_poly <- as(grd_sp, 'SpatialPolygons')
grd_poly <- spTransform(grd_poly, CRS(projutm))


B2 <- as_Spatial(st_intersection(st_as_sf(grd_poly),st_as_sf(B)))
B2.Cent <- SpatialPointsDataFrame( gCentroid(B2,byid=TRUE), data=data.frame(weight=rep(1,length(B2))))

B2.Cent$time.idx <- 1