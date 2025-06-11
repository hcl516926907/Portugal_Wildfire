dir.data <- "/home/pgrad2/2448355h/My_PhD_Project/00_Dataset/Porgual_Wildfire"
dir.out <- '/home/pgrad2/2448355h/My_PhD_Project/01_Output/Portugal_Wildfire'


library(INLA)

inla.setOption(scale.model.default = TRUE)

library(mgcv)
library(ggplot2)
# library(rgdal)
# library(fmesher)
library(dplyr)
library(RColorBrewer)
# library(terra)
# library(rgeos)
library(spdep)
# library(raster)
# library(scoringRules)
# library(fastDummies)
# bru_safe_sp(force = TRUE)
library(sf)
library(units)

load(file.path(dir.out, 'AutoRegressive_XGBoost_Predictions.RData'))


dist<-read_sf(file.path(dir.data,'shapefile', "distritos.shp"))
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
sf_districts <- st_make_valid(dist)

# map <- conc[,'NAME_2']
conc<-read_sf(file.path(dir.data, "shapefile","concelhos.shp"))
conc$ID_0 <- as.factor(iconv(as.character(conc$ID_0),  "UTF-8"))
conc$ISO <- as.factor(iconv(as.character(conc$ISO), "UTF-8"))
conc$NAME_0 <- as.factor(iconv(as.character(conc$NAME_0), "UTF-8"))
conc$NAME_1 <- as.factor(iconv(as.character(conc$NAME_1), "UTF-8"))
conc$ID_2 <- as.factor(iconv(as.character(conc$ID_2), "UTF-8"))
conc$NAME_2 <- as.factor(iconv(as.character(conc$NAME_2), "UTF-8"))
conc$HASC_2 <- as.factor(iconv(as.character(conc$HASC_2),  "UTF-8"))
conc$CCN_2 <- as.factor(iconv(as.character(conc$CCN_2),  "UTF-8"))
conc$CCA_2 <- as.factor(iconv(as.character(conc$CCA_2),  "UTF-8"))
conc$TYPE_2 <- as.factor(iconv(as.character(conc$TYPE_2),"UTF-8"))
conc$ENGTYPE_2 <- as.factor(iconv(as.character(conc$ENGTYPE_2), "UTF-8"))
conc$NL_NAME_2 <- as.factor(iconv(as.character(conc$NL_NAME_2), "UTF-8"))
conc$VARNAME_2 <- as.factor(iconv(as.character(conc$VARNAME_2), "UTF-8"))
conc=conc[conc$NAME_1!="Azores",]
conc=conc[conc$NAME_1!="Madeira",]
conc$NAME_1<-as.factor(droplevels(conc$NAME_1))
conc$NAME_2<-as.factor(droplevels(conc$NAME_2))
sf_conc <- st_make_valid(conc)

map <- sf_conc[,'NAME_2']
rownames(map) <- sf_conc$NAME_2

nb <- poly2nb(map)
map$grid.idx <- 1:nrow(map)
centroids <- st_centroid(map)
centroid_coords <- st_coordinates(centroids)
map$lon.grid <- centroid_coords[,1]
map$lat.grid <- centroid_coords[,2]
# 
# ggplot(map) + geom_sf()+
#   geom_text(aes(lon.grid, lat.grid, label = grid.idx), size=2,color='red')
# 

nb2INLA("map.adj", nb)
g.council <- inla.read.graph(filename = "map.adj")

# 
# 
map.district <- sf_districts[,'NAME_1']
rownames(map.district) <- sf_districts$NAME_1

nb.district <- poly2nb(map.district)
map.district$grid.idx.district <-  1:nrow(map.district)
centroids.district <- st_centroid(map.district)
centroid_coords.district <- st_coordinates(centroids.district)
map.district$lon.grid <- centroid_coords.district[,1]
map.district$lat.grid <- centroid_coords.district[,2]

# ggplot(map.district) + geom_sf()+
#   geom_text(aes(lon.grid, lat.grid, label = grid.idx.district), size=2,color='red')

nb2INLA("map.district.adj", nb.district)
g.district <- inla.read.graph(filename = "map.district.adj")

save(map, map.district,  file=file.path(dir.out, 'map_council_district.RData'))


data.fit <- subset(data.fit, select=-c(grid.idx))

data.fit <- merge(data.fit, as.data.frame(st_drop_geometry(map)[,c('NAME_2','grid.idx')]), on ='NAME_2')
data.fit <- merge(data.fit, as.data.frame(st_drop_geometry(map.district)[,c('NAME_1','grid.idx.district')]), on = 'NAME_1')



data.fit <- data.fit[order(data.fit$time.idx,data.fit$grid.idx),]

save(data.fit, file=file.path(dir.out, 'dataset_perf_evaluate.RData'))

z <- as.vector((data.fit$y>0)+0)
ba <- as.vector(ifelse(data.fit$y>0, (data.fit$area_ha)^{1/2}, NA))
cnt <- as.vector(ifelse(data.fit$y>0, data.fit$y, NA))


#prepare for prediction
z[which(data.fit$year>2022)] <- NA
ba[which(data.fit$year>2022)] <- NA
cnt[which(data.fit$year>2022)] <- NA



n.group <- 20

n1 <- dim(data.fit)[1]
nothing <- rep(NA, n1)

cntNA <- as.vector(c(cnt,nothing,nothing))
zNA <- as.vector(c(nothing, z, nothing))
baNA = as.vector(c(nothing, nothing, ba))

Y <- matrix(c(cntNA,zNA, baNA), ncol=3)

grid.idx.cnt = c(data.fit$grid.idx, nothing, nothing)# fire ignition
grid.idx.z = c(nothing,data.fit$grid.idx, nothing)# BA ind
grid.idx.ba = c(nothing, nothing, data.fit$grid.idx)# BA

grid.idx.dist.cnt = c(data.fit$grid.idx.district, nothing, nothing)# fire ignition
grid.idx.dist.z = c(nothing,data.fit$grid.idx.district, nothing)# BA ind
grid.idx.dist.ba = c(nothing, nothing, data.fit$grid.idx.district)# BA

month.idx.cnt <-  c(data.fit$month, nothing, nothing)
month.idx.z <-  c(nothing, data.fit$month, nothing)
month.idx.ba <-  c(nothing, nothing, data.fit$month)

year.idx.cnt <-  c(data.fit$year, nothing, nothing)
year.idx.z <-  c(nothing, data.fit$year, nothing)
year.idx.ba <-  c(nothing, nothing, data.fit$year)

intercept.cnt <- c(rep(1,n1),nothing, nothing)
intercept.z <- c(nothing, rep(1,n1), nothing)
intercept.ba <- c(nothing, nothing, rep(1,n1))

time.idx <- data.fit$time.idx - min(data.fit$time.idx) + 1
time.idx.cnt <- c(time.idx, nothing, nothing)
time.idx.z <- c(nothing, time.idx, nothing)
time.idx.ba <- c(nothing, nothing, time.idx)

rw.group <- function(x1,x2,n){
  q.bins <- c(0, ppoints(n-1), 1)
  q.bins.value <- lead((q.bins+ lag(q.bins))/2)[1:n]
  bins <- unique(quantile(x1, probs = q.bins))
  bins.value <-  unique(quantile(x1, probs = q.bins.value))
  bins[1] <- min(bins[1], min(x2))
  bins[length(bins)] <- max(bins[length(bins)], max(x2))
  x2.bins <- cut(x2, breaks = bins, include.lowest = TRUE)
  levels(x2.bins) <- bins.value
  x2.bins <- as.numeric(as.character(x2.bins))
  return(x2.bins)
}


score_cnt_grp <- rw.group(data.fit[data.fit$y>0,'score_cnt_h1'], data.fit$score_cnt_h1,n.group)
score.cnt <- c(score_cnt_grp,nothing, nothing)

# score_z_grp <- rw.group(data.fit$score_z_h1, data.fit$score_z_h1, n.group)
# score.z <- c(nothing, score_z_grp, nothing)
 
score_z_grp_cnt <- rw.group(data.fit$score_cnt_h1, data.fit$score_cnt_h1, n.group)
score.z_cnt <- c(nothing, score_z_grp_cnt, nothing)

score_z_grp_ba <- rw.group(data.fit$score_ba_h1, data.fit$score_ba_h1, n.group)
score.z_ba <- c(nothing, score_z_grp_ba, nothing)


score_ba_grp <- rw.group(data.fit[data.fit$y>0,'score_ba_h1'],data.fit$score_ba_h1,n.group)
score.ba <- c(nothing, nothing, score_ba_grp)


month.idx.cnt <- c(data.fit$month, nothing, nothing)# fire ignition
month.idx.z <- c(nothing, data.fit$month, nothing)# BA ind
month.idx.ba  <-  c(nothing, nothing, data.fit$month)# BA

year.idx <- data.fit$year - min(data.fit$year) + 1
year.idx.cnt <- c(year.idx, nothing, nothing)# fire ignition
year.idx.z <- c(nothing, year.idx, nothing)# BA ind
year.idx.ba  <-  c(nothing, nothing,year.idx)# BA


data.list=list(Y=Y, 

          grid.idx.cnt=grid.idx.cnt, 
          grid.idx.z=grid.idx.z, 
          grid.idx.ba=grid.idx.ba,
          
          grid.idx.dist.cnt = grid.idx.dist.cnt,
          grid.idx.dist.z = grid.idx.dist.z,
          grid.idx.dist.ba = grid.idx.dist.ba,
 
          intercept.cnt = intercept.cnt,
          intercept.z = intercept.z,
          intercept.ba = intercept.ba,
          
          month.idx.cnt = month.idx.cnt,
          month.idx.z = month.idx.z,
          month.idx.ba = month.idx.ba,
          
          year.idx.cnt = year.idx.cnt,
          year.idx.z = year.idx.z,
          year.idx.ba = year.idx.ba,
          
          time.idx.cnt = time.idx.cnt,
          time.idx.z = time.idx.z,
          time.idx.ba = time.idx.ba,
          
          
          score.cnt = score.cnt,
          # score.z = score.z,
          score.z_cnt = score.z_cnt,
          score.z_ba = score.z_ba,
          score.ba = score.ba
          
)




prec.prior <-  list(prec = list(prior="loggamma",param=c(0.1,0.1)))


 
formula <- Y ~  0  +
  intercept.cnt +
  f(grid.idx.cnt, copy='grid.idx.z',fixed=F,group=month.idx.cnt, control.group = list(model='iid'))  +
  f(year.idx.cnt, model='iid', hyper = prec.prior)+
  f(score.cnt, model='rw1', hyper = prec.prior )+
  f(grid.idx.dist.cnt, copy="grid.idx.dist.z",fixed=F, group=time.idx.cnt, control.group = list(model='iid'))  +


  intercept.z +
  f(grid.idx.z, model='bym2',graph=g.council,group=month.idx.z, control.group = list(model='iid'))  +
  f(year.idx.z, model='iid',hyper = prec.prior)+
  # f(score.z, model='rw1')+
  f(score.z_cnt, model='rw1', hyper = prec.prior) +  f(score.z_ba, model='rw1', hyper = prec.prior)+
  f(grid.idx.dist.z, model='bym2', graph=g.district,group=time.idx.z, control.group = list(model='iid')) +


  intercept.ba +
  f(grid.idx.ba, copy='grid.idx.z',fixed=F,group=month.idx.ba, control.group = list(model='iid'))  +
  f(year.idx.ba, model='iid',hyper = prec.prior) +
  f(score.ba, model='rw1', hyper = prec.prior ) +
  f(grid.idx.dist.ba, copy="grid.idx.dist.z",fixed=F, group=time.idx.ba, control.group = list(model='iid'))



#### Pure spatial-temporal effect
# formula <- Y ~  0  +
#   intercept.cnt + 
#   f(grid.idx.cnt, copy='grid.idx.z',fixed=F,group=month.idx.cnt, control.group = list(model='iid'))  +
#   f(year.idx.cnt, model='iid', hyper = prec.prior)+ 
#   # f(score.cnt, model='rw1', hyper = prec.prior )+
#   f(grid.idx.dist.cnt, copy="grid.idx.dist.z",fixed=F, group=time.idx.cnt, control.group = list(model='iid'))  +
#   
#   
#   intercept.z + 
#   f(grid.idx.z, model='bym2',graph=g.council,group=month.idx.z, control.group = list(model='iid'))  +
#   f(year.idx.z, model='iid',hyper = prec.prior)+  
#   # f(score.z_cnt, model='rw1', hyper = prec.prior) +  f(score.z_ba, model='rw1', hyper = prec.prior)+
#   f(grid.idx.dist.z, model='bym2', graph=g.district,group=time.idx.z, control.group = list(model='iid')) +
#   
#   
#   intercept.ba + 
#   f(grid.idx.ba, copy='grid.idx.z',fixed=F,group=month.idx.ba, control.group = list(model='iid'))  +
#   f(year.idx.ba, model='iid',hyper = prec.prior) +  
#   # f(score.ba, model='rw1', hyper = prec.prior ) +
#   f(grid.idx.dist.ba, copy="grid.idx.dist.z",fixed=F, group=time.idx.ba, control.group = list(model='iid'))



# save.image(file=file.path(dir.out, 'tmp.RData'))
# load(file.path(dir.out,'tmp.RData'))
# t1 <- Sys.time()
res <- inla(formula,
               family = c('nzpoisson','binomial', 'gamma'), data = data.list,  Ntrials=1,
               control.predictor = list(compute = TRUE, link=c(rep(1,n1),rep(2,n1),rep(3,n1))),
               verbose=TRUE,
               control.compute=list(config = TRUE),
               control.fixed = list(expand.factor.strategy = 'inla')
)

# res <- inla(formula,
#             family = c('nzpoisson','binomial', 'egp'), data = data.list,  Ntrials=1,
#             control.predictor = list(compute = TRUE, link=c(rep(1,n1),rep(2,n1),rep(3,n1))),
# 
#             verbose=TRUE,
#             control.compute=list(config = TRUE),
#             control.family = list( list(), list(), list(control.link = list(quantile = 0.5),
#                                                         hyper = list(
#                                                           tail = list(
#                                                           ##initial = xi.intern,
#                                                           fixed = !TRUE,
#                                                           prior = "pc.egptail",
#                                                           param = c(10, -0.5, 0.5)),
#                                                           shape = list(
#                                                             ##initial = kappa.intern,
#                                                             fixed = !TRUE,
#                                                             prior = "loggamma",
#                                                             param = c(100, 100)
#                                                             )
#                                                           )
#                                                         )),
#             control.fixed = list(expand.factor.strategy = 'inla')
# )
# 
# t2 <- Sys.time()
# print(t2-t1)

summary(res)

# save(res, file=file.path(dir.out, 'Model_egp_bym2.RData'))
# save(res, file=file.path(dir.out, 'Model_weibull_bym2.RData'))
# save(res, file=file.path(dir.out, 'Model_gamma_bym2.RData'))
gc()
n1 = nrow(data.fit)
# save(n1, res, file=file.path(dir.out, 'Model_egp_bym2_h1.RData'))
save(n1, res, file=file.path(dir.out, 'Model_gamma_bym2_h1.RData'))
# save(n1, res, file=file.path(dir.out, 'Model_weibull_bym2_h1.RData'))
# save(n1, res, file=file.path(dir.out, 'Model_egp_bym2_h1_noxgb.RData'))
gc()