dir.data <- "/home/pgrad2/2448355h/My_PhD_Project/00_Dataset/Porgual_Wildfire"
dir.out <- '/home/pgrad2/2448355h/My_PhD_Project/01_Output/Portugal_Wildfire'


library(INLA)

inla.setOption(scale.model.default = TRUE)

library(mgcv)
library(ggplot2)
library(rgdal)
library(sf)
library(fmesher)
library(dplyr)
library(RColorBrewer)
library(terra)
library(rgeos)
library(spdep)
library(raster)
library(scoringRules)
library(fastDummies)
# bru_safe_sp(force = TRUE)


conc<-shapefile(file.path(dir.data, "shapefile","concelhos.shp"))
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

post.pred.weibull.hurdle.parallel <- function(samples, idx.pred.pois, idx.pred.z, idx.pred.ba,n.samples=200 ){
  t1 <- Sys.time()
  # res <- foreach(j = 1:length(idx.pred.pois) ) %dopar%{
  pred.cnt.mat <- c()
  pred.z.mat <- c()
  pred.ba.mat <- c()
  param.cnt <- list()
  param.z <- list()
  param.ba <- list()
  hyper.ba <- list()
  
  for (i in 1:n.samples){
    eta.pois <- samples[[i]]$latent[idx.pred.pois,1]
    lambda.pois <- exp(eta.pois)
    
    
    eta.z <- samples[[i]]$latent[idx.pred.z,1]
    p <- exp(eta.z)/(1 + exp(eta.z))
    
    
    eta.ba <- samples[[i]]$latent[idx.pred.ba,1]
    alpha.ba <- samples[[i]]$hyperpar[1]
    
    lambda.ba  <-  exp(eta.ba)
    scale.ba <- lambda.ba^(-1/alpha.ba)
    
    pred.z <- rbinom(length(idx.pred.pois), size=1, prob=p)
    pred.cnt <- rpois(length(idx.pred.z), lambda.pois )
    pred.ba <- rweibull(length(idx.pred.ba), shape = alpha.ba, scale = scale.ba)
    
    zero.idx <- which(pred.z==0)
    pred.cnt[zero.idx] <- 0
    pred.ba[zero.idx] <- 0
    
    
    pred.cnt.mat <- cbind(pred.cnt.mat, pred.cnt)
    pred.z.mat <- cbind(pred.z.mat, pred.z)
    pred.ba.mat <- cbind(pred.ba.mat, pred.ba)
    
    param.cnt[[i]] <- lambda.pois
    param.z[[i]] <- p
    param.ba[[i]] <- scale.ba
    hyper.ba[[i]] <- alpha.ba
  }
  
  
  t2 <- Sys.time()
  print(t2-t1)
  pred.cnt <- list()
  pred.z <- list()
  pred.ba <- list()
  for (i in 1:length(idx.pred.pois)){
    pred.cnt[[i]] <- as.numeric(pred.cnt.mat[i,])
    pred.z[[i]] <- as.numeric(pred.z.mat[i,])
    pred.ba[[i]] <- as.numeric(pred.ba.mat[i,])
  }
  
  return(list('pred.ba'=pred.ba, 'pred.cnt'=pred.cnt, 'pred.z'=pred.z ))
}


post.pred.egp.hurdle.parallel <- function(samples, idx.pred.pois, idx.pred.z, idx.pred.ba,n.samples=200 ){
  t1 <- Sys.time()
  # res <- foreach(j = 1:length(idx.pred.pois) ) %dopar%{
  pred.cnt.mat <- c()
  pred.z.mat <- c()
  pred.ba.mat <- c()
  param.cnt <- list()
  param.z <- list()
  param.ba <- list()
  hyper.xi <- list()
  hyper.kappa <- list()
  
  regp <- function(n, eta, kappa, xi,alpha=0.5){
    a <- ((1-alpha^(1/kappa))^(-xi) -1)
    # x <- rnorm(n, sd = 0.3)
    # eta <- 0.9 + 1.1 * x
    q <- exp(eta)
    sigma <- xi * q / a
    y <- - (sigma / xi) * (1- (1-runif(n)^(1/kappa))^(-xi))
    return(y)
  }
  
  
  for (i in 1:n.samples){
    eta.pois <- samples[[i]]$latent[idx.pred.pois,1]
    lambda.pois <- exp(eta.pois)
    
    
    eta.z <- samples[[i]]$latent[idx.pred.z,1]
    p <- exp(eta.z)/(1 + exp(eta.z))
    
    
    eta.ba <- samples[[i]]$latent[idx.pred.ba,1]
    # xi <- samples[[i]]$hyperpar[1]
    xi <- -0.465
    kappa <- samples[[i]]$hyperpar[2]
    
    pred.z <- rbinom(length(idx.pred.pois), size=1, prob=p)
    pred.cnt <- rpois(length(idx.pred.z), lambda.pois )
    pred.ba <- regp(length(idx.pred.ba), eta = eta.ba, kappa=kappa, xi=xi)
    
    zero.idx <- which(pred.z==0)
    pred.cnt[zero.idx] <- 0
    pred.ba[zero.idx] <- 0
    
    
    pred.cnt.mat <- cbind(pred.cnt.mat, pred.cnt)
    pred.z.mat <- cbind(pred.z.mat, pred.z)
    pred.ba.mat <- cbind(pred.ba.mat, pred.ba)
    
    param.cnt[[i]] <- lambda.pois
    param.z[[i]] <- p
    param.ba[[i]] <- eta.ba
    hyper.xi[[i]] <- xi
    hyper.kappa[[i]] <- kappa
  }
  
  
  t2 <- Sys.time()
  print(t2-t1)
  pred.cnt <- list()
  pred.z <- list()
  pred.ba <- list()
  for (i in 1:length(idx.pred.pois)){
    pred.cnt[[i]] <- as.numeric(pred.cnt.mat[i,])
    pred.z[[i]] <- as.numeric(pred.z.mat[i,])
    pred.ba[[i]] <- as.numeric(pred.ba.mat[i,])
  }
  
  return(list('pred.ba'=pred.ba, 'pred.cnt'=pred.cnt, 'pred.z'=pred.z,
              'param.cnt'=param.cnt, 'param.z'=param.z, 'param.ba'=param.ba, 'hyper.xi'=hyper.xi,
              'hyper.kappa'=hyper.kappa ))
}


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



# district <- 'Viana do Castelo'
# district <- 'Viseu'
load(file.path(dir.out, "XGBoost_Score_Council.RData"))

data.fit$grid.idx.copy <- data.fit$grid.idx
pred.sp.all <- list()
####################################################################
# INLA
####################################################################
summary.name <- file.path(dir.out, 'INLA_Summary_Weibull.txt')
if (file.exists(summary.name)) {
  file.remove(summary.name)
}
for (district in unique(data.fit$NAME_1)){
  data.fit.dist <- data.fit[data.fit$NAME_1==district,]
  sf_conc <- st_as_sf(conc[conc$NAME_1==district,])
  
  map <- sf_conc[,'NAME_2']
  rownames(map) <- sf_conc$NAME_2
  
  nb <- poly2nb(map)
  map$grid.idx <- 1:nrow(map)
  centroids <- st_centroid(map)
  centroid_coords <- st_coordinates(centroids)
  map$lon.grid <- centroid_coords[,1]
  map$lat.grid <- centroid_coords[,2]
  
  ggplot(map) + geom_sf()+
    geom_text(aes(lon.grid, lat.grid, label = grid.idx), size=2,color='red') 
  
  
  nb2INLA("map.adj", nb)
  g.council <- inla.read.graph(filename = "map.adj")
  
  
  data.fit.dist <- data.fit.dist |> subset(select=-grid.idx) |>
    st_as_sf(coords=c('lon.grid','lat.grid'),crs=4326,remove=FALSE)|>
    st_join(map[,'grid.idx'])|>st_drop_geometry()
  
  
  
  data.fit.dist$month <- (data.fit.dist$time.idx-1)%%12 + 1
  data.fit.dist$year <- (data.fit.dist$time.idx-1)%/%12 + 2011
  
  
  n.group <- 5
  
  
  data.fit.dist$z <- as.vector((data.fit.dist$y>0)+0)
  
  log.ba.agg <- data.fit.dist[data.fit.dist$year<=2022,] %>% group_by(NAME_2) %>%
    summarise(cum_log_ba=sum(log_ba,na.rm=TRUE)) %>%
    mutate(log.ba.idx=inla.group(cum_log_ba,n=n.group,method='quantile',idx.only=T))
  
  cnt.agg <- data.fit.dist[data.fit.dist$year<=2022,] %>% group_by(NAME_2) %>%
    summarise(cum_cnt=sum(y,na.rm=TRUE)) %>%
    mutate(cnt.idx=inla.group(cum_cnt,n=n.group,method='quantile',idx.only=T))
  
  z.agg <- data.fit.dist[data.fit.dist$year<=2022,] %>% group_by(NAME_2) %>%
    summarise(cum_z=sum(z,na.rm=TRUE)) %>%
    mutate(z.idx=inla.group(cum_z,n=n.group,method='quantile',idx.only=T))
  
  
  data.fit.dist <- merge(data.fit.dist,log.ba.agg[,c("NAME_2",'log.ba.idx')],by='NAME_2')
  data.fit.dist <- merge(data.fit.dist, cnt.agg[,c('NAME_2','cnt.idx')], by='NAME_2')
  data.fit.dist <- merge(data.fit.dist, z.agg[,c('NAME_2','z.idx')], by='NAME_2')
  
  
  
  data.fit.dist <- data.fit.dist[order(data.fit.dist$time.idx,data.fit.dist$grid.idx.copy),]
  
  z <- as.vector((data.fit.dist$y>0)+0)
  log.ba <- as.vector(ifelse(data.fit.dist$y>0, data.fit.dist$log_ba, NA))
  cnt <- as.vector(ifelse(data.fit.dist$y>0, data.fit.dist$y, NA))
  
  
  #prepare for prediction
  z[which(data.fit.dist$year>2022)] <- NA
  log.ba[which(data.fit.dist$year>2022)] <- NA
  cnt[which(data.fit.dist$year>2022)] <- NA
  
  n1 <- dim(data.fit.dist)[1]
  nothing <- rep(NA, n1)
  
  cntNA <- as.vector(c(cnt,nothing,nothing))
  zNA <- as.vector(c(nothing, z, nothing))
  baNA.log = as.vector(c(nothing, nothing, log.ba))
  
  Y.log <- matrix(c(cntNA,zNA, baNA.log), ncol=3)
  
  grid.idx.cnt = c(data.fit.dist$grid.idx, nothing, nothing)# fire ignition
  grid.idx.z = c(nothing,data.fit.dist$grid.idx, nothing)# BA ind
  grid.idx.ba = c(nothing, nothing, data.fit.dist$grid.idx)# BA
  
  
  time.idx.cnt <- c(data.fit.dist$time.idx, nothing, nothing)# fire ignition
  time.idx.z <- c(nothing, data.fit.dist$time.idx, nothing)# BA ind
  time.idx.ba  <-  c(nothing, nothing, data.fit.dist$time.idx)# BA
  
  intercept.cnt <- c(rep(1,n1),nothing, nothing)
  intercept.z <- c(nothing, rep(1,n1), nothing)
  intercept.ba <- c(nothing, nothing, rep(1,n1))
  

  
  score_cnt_grp <- rw.group(data.fit.dist[data.fit.dist$y>0,'score_cnt'],data.fit.dist$score_cnt,40)
  score.cnt <- c(score_cnt_grp,nothing, nothing)
  
  score_z_grp <- rw.group(data.fit.dist$score_z,data.fit.dist$score_z,40)
  score.z <- c(nothing, score_z_grp, nothing)
  
  score_ba_grp <- rw.group(data.fit.dist[data.fit.dist$y>0,'score_ba'],data.fit.dist$score_ba,40)
  score.ba <- c(nothing, nothing, score_ba_grp)
  
  month.idx.cnt <- c(data.fit.dist$month, nothing, nothing)# fire ignition
  month.idx.z <- c(nothing, data.fit.dist$month, nothing)# BA ind
  month.idx.ba  <-  c(nothing, nothing, data.fit.dist$month)# BA
  
  year.idx <- data.fit.dist$year - 2010
  year.idx.cnt <- c(year.idx, nothing, nothing)# fire ignition
  year.idx.z <- c(nothing, year.idx, nothing)# BA ind
  year.idx.ba  <-  c(nothing, nothing,year.idx)# BA
  
  after.2017 <- as.integer(data.fit.dist$year>2017)
  year.af.2017.cnt <- c(after.2017, nothing, nothing)
  year.af.2017.z <- c(nothing, after.2017, nothing)
  year.af.2017.ba <- c(nothing, nothing,  after.2017)
  
  log.ba.idx <- c(nothing, nothing,  data.fit.dist$log.ba.idx)
  cnt.idx <- c(data.fit.dist$cnt.idx, nothing, nothing)
  z.idx <- c(nothing, data.fit.dist$z.idx, nothing)
  
  
  data.list=list(Y.log=Y.log, 
                 
                 grid.idx.cnt=grid.idx.cnt, 
                 grid.idx.z=grid.idx.z, 
                 grid.idx.ba=grid.idx.ba,
                 
                 time.idx.cnt = time.idx.cnt,
                 time.idx.z = time.idx.z,
                 time.idx.ba = time.idx.ba,
                 
                 time.idx.cnt1 = time.idx.cnt,
                 time.idx.z1 = time.idx.z,
                 time.idx.ba1 = time.idx.ba,
                 
                 intercept.cnt = intercept.cnt,
                 intercept.z = intercept.z,
                 intercept.ba = intercept.ba,
                 
                 month.idx.cnt = month.idx.cnt,
                 month.idx.z = month.idx.z,
                 month.idx.ba = month.idx.ba,
                 
                 year.idx.cnt = year.idx.cnt,
                 year.idx.z = year.idx.z,
                 year.idx.ba = year.idx.ba,
                 
                 year.af.2017.cnt = year.af.2017.cnt,
                 year.af.2017.z = year.af.2017.z,
                 year.af.2017.ba = year.af.2017.ba,
                 
                 log.ba.idx = log.ba.idx,
                 cnt.idx = cnt.idx,
                 z.idx = z.idx,
                 
                 score.cnt = score.cnt,
                 score.z = score.z,
                 score.ba = score.ba
                 
  )
  
  
  formula <- Y.log ~  0  +
    intercept.cnt + year.af.2017.cnt +
    f(grid.idx.cnt, copy='grid.idx.z',fixed=F)  +
    f(score.cnt, model='rw1')+
    f(time.idx.cnt, model='iid')+ 
    f(year.idx.cnt, model='iid') + 
    # f(cnt.idx, model='rw1')+
    
    intercept.z + year.af.2017.z +
    f(grid.idx.z, model='bym2',graph=g.council)  +
    f(score.z, model='rw1' )+
    f(time.idx.z, model='iid')+ 
    f(year.idx.z, model='iid') + 
    # f(z.idx, model='rw1')+
    
    
    intercept.ba + year.af.2017.ba +
    f(grid.idx.ba, copy='grid.idx.z',fixed=F)  +
    f(score.ba, model='rw1' ) +
    f(time.idx.ba, model='iid') + 
    f(year.idx.ba, model='iid') 
    # f(log.ba.idx, model='rw1')
  
  hyper.theta <-  list(prec = list(prior="loggamma",param=c(10,10)))
  
  formula.egp <- Y.log ~  0  +
    intercept.cnt + year.af.2017.cnt +
    f(grid.idx.cnt, copy='grid.idx.z',fixed=F)  +
    f(score.cnt, model='rw1')+
    f(time.idx.cnt, model='iid')+ 
    
    intercept.z + year.af.2017.z +
    f(grid.idx.z, model='bym2',graph=g.council)  +
    f(score.z, model='rw1' )+
    f(time.idx.z, model='iid')+ 
    
    
    intercept.ba + year.af.2017.ba +
    f(grid.idx.ba, copy='grid.idx.z',fixed=F)  +
    f(score.ba, model='rw1' ) +
    f(time.idx.ba, model='iid')
  
  
  n1 <- dim(data.fit.dist)[1]
  
  res <- inla(formula,
              family = c('poisson','binomial', 'weibull'), data = data.list,  Ntrials=1,
              control.predictor = list( compute = TRUE, link=c(rep(1,n1),rep(2,n1),rep(3,n1))),
              verbose=TRUE,
              control.compute=list(config = TRUE),
              control.family = list( list(), list(), list()),
              control.fixed = list(expand.factor.strategy = 'inla')
  )
  
  # 
  # res <- inla(formula.egp,
  #             family = c('poisson','binomial', 'egp'), data = data.list,  Ntrials=1,
  #             control.predictor = list( compute = TRUE, link=c(rep(1,n1),rep(2,n1),rep(3,n1))),
  #             verbose=TRUE,
  #             control.compute=list(config = TRUE),
  #             control.family = list( list(), list(), list(control.link = list(quantile = 0.6),
  #                                                         hyper = list(
  #                                                           tail = list(
  #                                                             ##initial = xi.intern,
  #                                                             fixed = !TRUE,
  #                                                             prior = "pc.egptail",
  #                                                             param = c(5, -0.5, 0.5)),
  #                                                           shape = list(
  #                                                             ##initial = kappa.intern,
  #                                                             fixed = !TRUE,
  #                                                             prior = "loggamma",
  #                                                             param = c(40, 20)
  #                                                           )
  #                                                         )
  #             )),
  #             control.fixed = list(expand.factor.strategy = 'inla')
  # )
  
  summary(res)
  # summary.name <- paste("INLA Summary Weibull ",district,'.txt',sep='')

  sink(summary.name,append = TRUE)
  print(summary(res))
  sink() 
  
  
  n.samples <- 2000
  samples = inla.posterior.sample(n.samples, result = res, seed=1234)
  
  

  idx.pred.pois <- 1:n1
  idx.pred.z <- (n1+1):(2*n1)
  idx.pred.ba <- (2*n1+1):(3*n1)
  
  
  
  
  t1 <- Sys.time()
  pred.sp <- post.pred.weibull.hurdle.parallel(samples, idx.pred.pois, idx.pred.z, idx.pred.ba, n.samples=n.samples)
  # pred.sp <- post.pred.egp.hurdle.parallel(samples, idx.pred.pois, idx.pred.z, idx.pred.ba, n.samples=n.samples)
  # pred.sp <- post.pred.gamma.hurdle.parallel(samples, idx.pred.pois, idx.pred.z, idx.pred.ba, n.samples=n.samples)
  
  t2 <- Sys.time()
  print(t2-t1)
  
  
  pred.sp.all[[district]] <- pred.sp
  print(paste('Done',district))
}
save(pred.sp.all, file=file.path(dir.out,'single_district_weibull_pred_sp.RData'))

###############################################################################
load(file=file.path(dir.out,'single_district_weibull_pred_sp.RData'))
pred.cnt <- c()
pred.ba <- c()
pred.z <- c()
for (t in unique(data.fit$time.idx)){
  for (district in unique(data.fit$NAME_1)){
    n.council <- length(unique(data.fit[data.fit$NAME_1==district,'NAME_2']))
    time.cond <- ((t-1)*n.council+1) : (t*n.council)
    pred.cnt <- c(pred.cnt,pred.sp.all[[district]]$pred.cnt[time.cond] )
    pred.ba <- c(pred.ba,pred.sp.all[[district]]$pred.ba[time.cond] )
    pred.z <- c(pred.z,pred.sp.all[[district]]$pred.z[time.cond] )
  }

}

rm(pred.sp.all)



crps.ba <- rep(NA, nrow(data.fit))
for (i in 1:nrow(data.fit) ){
  y <- data.fit$log_ba[i]
  if (is.na(y)) y <- 0
  crps.ba[i] <- crps_sample(
    y,
    pred.ba[[i]],
    method = "edf")
}
round(sum(crps.ba)/length(crps.ba),4)

df.crps.ba <- data.frame(crps=crps.ba,year=data.fit$year,district=data.fit$NAME_1)
df.crps.ba$label <- 'Train'
df.crps.ba[which(df.crps.ba$year>2022),'label'] <- 'Test'
df.crps.ba$label <- factor(df.crps.ba$label, level=c('Train','Test'))

print(df.crps.ba %>% group_by(district,label) %>% summarize(avg_crps_ba=mean(crps)),n=50)
print(df.crps.ba %>% group_by(label) %>% summarize(avg_crps_ba=mean(crps)),n=50)


crps.t.ba <- rep(NA, nrow(data.fit))
for (i in 1:nrow(data.fit) ){
  y <- data.fit$y[i]
  if (is.na(y)) y <- 0
  if (y>log(50)){
    crps.t.ba[i] <- crps_sample(
      y,
      pred.ba[[i]],
      method = "edf")
  }
}

df.crps.t.ba <- data.frame(crps.t=crps.t.ba,year=data.fit$year,district=data.fit$NAME_1)
df.crps.t.ba$label <- 'Train'
df.crps.t.ba[which(df.crps.t.ba$year>2022),'label'] <- 'Test'
df.crps.t.ba$label <- factor(df.crps.t.ba$label, level=c('Train','Test'))

print(df.crps.t.ba %>% group_by(district,label) %>% summarize(avg_crps_t_ba=mean(crps.t,na.rm=T)),n=50)
print(df.crps.t.ba %>% group_by(label) %>% summarize(avg_crps_t_ba=mean(crps.t,na.rm=T)),n=50)




#
#
crps.cnt <- rep(NA, nrow(data.fit))
for (i in 1:nrow(data.fit) ){
  y <- data.fit$y[i]
  if (is.na(y)) y <- 0
  crps.cnt[i] <- crps_sample(
    y,
    pred.cnt[[i]],
    method = "edf")
}

df.crps.cnt <- data.frame(crps=crps.cnt,year=data.fit$year,district=data.fit$NAME_1)
df.crps.cnt$label <- 'Train'
df.crps.cnt[which(df.crps.cnt$year>2022),'label'] <- 'Test'
df.crps.cnt$label <- factor(df.crps.cnt$label, level=c('Train','Test'))

print(df.crps.cnt %>% group_by(district,label) %>% summarize(avg_crps_cnt=mean(crps)),n=50)
print(df.crps.cnt %>% group_by(label) %>% summarize(avg_crps_cnt=mean(crps)),n=50)


crps.t.cnt <- rep(NA, nrow(data.fit))
for (i in 1:nrow(data.fit) ){
  y <- data.fit$y[i]
  if (is.na(y)) y <- 0
  if (y>5){
    crps.t.cnt[i] <- crps_sample(
      y,
      pred.cnt[[i]],
      method = "edf")
  }
}

df.crps.t.cnt <- data.frame(crps.t=crps.t.cnt,year=data.fit$year,district=data.fit$NAME_1)
df.crps.t.cnt$label <- 'Train'
df.crps.t.cnt[which(df.crps.t.cnt$year>2022),'label'] <- 'Test'
df.crps.t.cnt$label <- factor(df.crps.t.cnt$label, level=c('Train','Test'))

print(df.crps.t.cnt %>% group_by(district,label) %>% summarize(avg_crps_t_cnt=mean(crps.t,na.rm=T)),n=50)
print(df.crps.t.cnt %>% group_by(label) %>% summarize(avg_crps_t_cnt=mean(crps.t,na.rm=T)),n=50)


data.fit[data.fit$y==0,'log_ba'] <- 0
data.fit$z <- as.integer(data.fit$y>0)
data.fit$label <- 'Train'
data.fit[which(data.fit$year>2022),'label'] <- 'Test'
data.fit$label <- factor(data.fit$label, level=c('Train','Test'))


data.fit[,'Estimated_Lambda'] <- sapply(pred.cnt,mean)
data.fit$Lower_CNT <- sapply(pred.cnt,quantile,0.25)
data.fit$Upper_CNT <- sapply(pred.cnt,quantile,0.75)
data.fit$residual_cnt <- data.fit$y - data.fit$Estimated_Lambda
summary(data.fit$residual_cnt)
data.fit$cover_cnt <- data.fit$Lower_CNT <= data.fit$y & data.fit$Upper_CNT>=data.fit$y

data.fit[,'Estimated_Z'] <- sapply(pred.z,mean)
data.fit$Lower_Z <- sapply(pred.z,quantile,0.25)
data.fit$Upper_Z <- sapply(pred.z,quantile,0.75)


data.fit[,'Estimated_BA'] <- sapply(pred.ba,mean)
data.fit[,'Lower_BA'] <- sapply(pred.ba,quantile,0.25)
data.fit[,'Upper_BA'] <- sapply(pred.ba,quantile,0.75)
data.fit$residual_ba <- data.fit$log_ba - data.fit$Estimated_BA
summary(data.fit$residual_ba)
data.fit$cover_ba <- data.fit$Lower_BA <= data.fit$log_ba & data.fit$Upper_BA>=data.fit$log_ba

library(pROC)
auc(data.fit[data.fit$label=='Train',  'z'], data.fit[data.fit$label=='Train',  'Estimated_Z'], quiet = TRUE)
auc(data.fit[data.fit$label=='Test',  'z'], data.fit[data.fit$label=='Test',  'Estimated_Z'], quiet = TRUE)


print(data.fit %>% group_by(label) %>% summarize(avg_res_cnt=mean((residual_cnt)^2)),n=50)
print(data.fit %>% group_by(NAME_1,label) %>% summarize(avg_res_cnt=mean((residual_cnt)^2)),n=50)

print(data.fit %>% group_by(label) %>% summarize(avg_res_ba=mean((residual_ba)^2)),n=50)
print(data.fit %>% group_by(NAME_1, label) %>% summarize(avg_res_ba=mean((residual_ba)^2)),n=50)



print(data.fit[data.fit$y>5,] %>% group_by(label) %>% summarize(avg_res_cnt=mean((residual_cnt)^2)),n=50)
print(data.fit[data.fit$y>5,] %>% group_by(NAME_1,label) %>% summarize(avg_res_cnt=mean((residual_cnt)^2)),n=50)

print(data.fit[data.fit$log_ba>log(50),] %>% group_by(label) %>% summarize(avg_res_ba=mean((residual_ba)^2)),n=50)
print(data.fit[data.fit$log_ba>log(50),] %>% group_by(NAME_1, label) %>% summarize(avg_res_ba=mean((residual_ba)^2)),n=50)





merged_sf1 <- st_as_sf(data.fit, coords = c("lon.grid", "lat.grid"), crs = 4326)

dist.name <-  district
joint.post.sp <- function(x) Reduce("+", x)
df.dist <- data.frame(time.idx=rep(1:156, each=length(pred.cnt[[1]])))
df.dist$year_label <- rep(2011:2023,each=length(pred.cnt[[1]])*12)
df.dist$month_label <- sprintf("%02d", (df.dist$time.idx-1)%%12+1)
df.dist$time_label <- paste(df.dist$year_label ,df.dist$month_label  ,sep='_')

df.boxplot.true <- data.frame(time.idx=1:156, cnt=rep(NA, 12*13), log.ba= rep(NA,12*13))
df.boxplot.true$year_label <- rep(2011:2023,each=12)
df.boxplot.true$month_label <- sprintf("%02d", (df.boxplot.true$time.idx-1)%%12+1)
df.boxplot.true$time_label <- paste(df.boxplot.true$year_label ,df.boxplot.true$month_label  ,sep='_')

for (t in 1:156){
  idx <- which(merged_sf1$NAME_1==dist.name & merged_sf1$time.idx==t)
  cnt.true <- sum(data.fit[idx,'y'])
  ba.true <- sum(data.fit[idx,'log_ba'],na.rm=T)

  df.boxplot.true[df.boxplot.true$time.idx==t,c('cnt','log.ba')] <- c(cnt.true, ba.true)

  pred.cnt.dist <- joint.post.sp(pred.cnt[idx])
  pred.ba.dist <- joint.post.sp(pred.ba[idx])
  df.dist[df.dist$time.idx==t,'sample_cnt'] <- pred.cnt.dist
  df.dist[df.dist$time.idx==t,'sample_ba'] <- pred.ba.dist
}

ggplot(df.dist, aes(x = factor(time_label), y = sample_cnt)) +
  geom_boxplot(outlier.shape=NA) +
  geom_line(data=df.boxplot.true, aes(x=factor(time_label),y=cnt,group = 1), col='red',linewidth=1)+
  labs(x = "Time", y = "fire count", title = paste("Predictive Total Fire Count in", dist.name, sep=' ')) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        axis.title.x = element_blank())  # Rotate x-axis labels for clarity+


ggplot(df.dist, aes(x = factor(time_label), y = sample_ba)) +
  geom_boxplot(outlier.shape=NA) +
  geom_line(data=df.boxplot.true, aes(x=factor(time_label),y=log.ba,group = 1), col='red',linewidth=1)+
  labs(x = "Time", y = "log burn area", title = paste("Predictive Total Burn Area in", dist.name, sep=' ')) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        axis.title.x = element_blank())  # Rotate x-axis labels for clarity+


#----------------------------------boxplot by district----------------------

# time.idx <- 82
time.idx <-  154

fire.dist.month <- st_drop_geometry(merged_sf1) %>% group_by(NAME_1, time.idx) %>% 
  summarize(cnt.true = sum(y),
            ba.true = sum(log_ba,na.rm=T),
  )

dist.vector <- fire.dist.month[fire.dist.month$time.idx>=97,] %>% group_by(NAME_1) %>% 
  summarize(cnt.true = sum(cnt.true),
            ba.true = sum(ba.true),
  )%>% arrange(-ba.true)

dist.vector <- as.character(dist.vector$NAME_1)

df.dist.month <- data.frame(dist=factor(rep(dist.vector, each=length(pred.cnt[[1]]))),
                            time.idx=time.idx, month=(time.idx-1)%%12 + 1)
levels(df.dist.month$dist) <- dist.vector


df.dist.month.true <- data.frame(dist=factor(dist.vector),time.idx=time.idx,
                                 month=(time.idx-1)%%12 + 1)
levels(df.dist.month.true$dist) <- dist.vector

##INLA output
for (dist in unique(df.dist.month$dist)){
  idx <- which(merged_sf1$NAME_1==dist & merged_sf1$time.idx==time.idx)
  
  pred.cnt.dist <- joint.post.sp(pred.cnt[idx])
  pred.ba.dist <- joint.post.sp(pred.ba[idx])
  df.dist.month[df.dist.month$dist==dist,'sample_cnt'] <- pred.cnt.dist
  df.dist.month[df.dist.month$dist==dist,'sample_ba'] <- pred.ba.dist
}

##XGBoost output
# cutoff <- 0.7
# xgb.pred.z <- (as.integer(data.fit$score_z>cutoff))
# xgb.pred.cnt <- xgb.pred.z*data.fit$score_cnt
# xgb.pred.ba <- xgb.pred.z*data.fit$score_ba
# for (dist in unique(df.dist.month$dist)){
#   idx <- which(merged_sf1$NAME_1==dist & merged_sf1$time.idx==time.idx)
#   
#   pred.cnt.dist <- sum(xgb.pred.cnt[idx])
#   pred.ba.dist <- sum(xgb.pred.ba[idx])
#   df.dist.month[df.dist.month$dist==dist,'sample_cnt'] <- pred.cnt.dist
#   df.dist.month[df.dist.month$dist==dist,'sample_ba'] <- pred.ba.dist
# }

df.dist.month.true <- merge(df.dist.month.true,fire.dist.month,
                            by.x=c('dist','time.idx'),
                            by.y=c('NAME_1','time.idx'))



ggplot(df.dist.month, aes(x = dist, y = sample_cnt)) +
  geom_boxplot(outlier.shape=NA) + 
  geom_line(data=df.dist.month.true, aes(x=dist,y=cnt.true,group = 1), col='red',linewidth=1)+
  # labs(x = "Time", y = "fire count", title = paste("Predictive Total Fire Count in", dist.name,"During 2019-2020", sep=' ')) +
  theme_minimal() +
  ggtitle(paste("time index",time.idx))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        axis.title.x = element_blank())  # Rotate x-axis labels for clarity+


ggplot(df.dist.month, aes(x = factor(dist), y = sample_ba)) +
  geom_boxplot(outlier.shape=NA) + 
  geom_line(data=df.dist.month.true, aes(x=factor(dist),y=ba.true,group = 1), col='red',linewidth=1)+
  # labs(x = "Time", y = "fire count", title = paste("Predictive Total Fire Count in", dist.name,"During 2019-2020", sep=' ')) +
  theme_minimal() +  
  ggtitle(paste("time index",time.idx))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        axis.title.x = element_blank())  # Rotate x-axis labels for clarity+


