dir.data <- "/home/pgrad2/2448355h/My_PhD_Project/00_Dataset/Urban_Fires"
dir.out <- '/home/pgrad2/2448355h/My_PhD_Project/01_Output/Urban_Fire_Model'


library(INLA)
library(mgcv)
library(ggplot2)
library(rgdal)
library(sf)
library(tibble)
library(dplyr)
library(rasterVis)

library(rnaturalearth)
map <- ne_countries(type = "countries", country = "Portugal",
                    scale = "medium", returnclass = "sf")
projutm <- "+proj=utm +zone=29 +datum=WGS84 +units=m +no_defs"
map <- st_transform(map, crs = projutm)
mainland_bbox <- st_as_sfc(st_bbox(c(xmin = 0, xmax = 1000000, ymin = 4000000 , ymax = 4800000), crs = st_crs(map)))
map_mainland <- st_intersection(map, mainland_bbox)


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


res_extract <- function(filename){
  load(file.path(dir.out, filename))

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
  
  levelplot(raster::brick(grid), layout = c(3, 1),
            names.attr = c("Mean", "2.5 percentile", "97.5 percentile"),par.settings= BuRdTheme())
  
  b0 <- res$marginals.fixed[[1]]
  ggplot(data.frame(inla.smarginal(b0)), aes(x, y)) +
    geom_line() +
    theme_bw()
  
  range <- res$marginals.hyperpar$`Range for s`
  ggplot(data.frame(inla.smarginal(range)), aes(x, y)) +
    geom_line() +
    theme_bw()
  
  
  stdev <- res$marginals.hyperpar$`Stdev for s`
  ggplot(data.frame(inla.smarginal(stdev)), aes(x, y)) +
    geom_line() +
    theme_bw()
  return(list('grid'=grid,'b0'=b0, 'range'=range, 'stdev'=stdev))
}



# filename1 <- "res_frac0.4_2200_10000_4000_0.05_8_0.05_.RData"
# filename2 <- "res_frac0.4_3500_10000_4000_0.05_8_0.05_.RData"
# filename3 <- "res_frac0.4_4000_10000_4000_0.05_8_0.05_.RData"
# filename4 <- "res_frac0.4_2800_20000_4000_0.05_8_0.05_.RData"

filename <- "res_frac0.4_8000_20000_9000_0.05_8_0.05_.RData"

res <- res_extract(filename)

levelplot(raster::brick(res$grid), layout = c(3, 1),
          names.attr = c("Mean", "2.5 percentile", "97.5 percentile"),par.settings= BuRdTheme())
ggplot(data.frame(inla.smarginal(res$b0)), aes(x, y)) +
  geom_line() +
  ggtitle('Posterior of b0') + 
  theme(plot.title = element_text(hjust = 0.5))

ggplot(data.frame(inla.smarginal(res$range)), aes(x, y)) +
  geom_line() +
  ggtitle('Posterior of range') + 
  theme(plot.title = element_text(hjust = 0.5))

ggplot(data.frame(inla.smarginal(res$stdev)), aes(x, y)) +
  geom_line() +
  ggtitle('Posterior of stdev') + 
  theme(plot.title = element_text(hjust = 0.5))

