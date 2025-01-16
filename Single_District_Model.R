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

# district <- 'Viana do Castelo'
district <- 'Évora'
load(file.path(dir.out, "XGBoost_Score_Council.RData"))
data.fit <- data.fit[data.fit$NAME_1==district,]


####################################################################
# INLA
####################################################################




# map <- conc[,'NAME_2']
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


data.fit <- data.fit |> subset(select=-grid.idx) |>
  st_as_sf(coords=c('lon.grid','lat.grid'),crs=4326,remove=FALSE)|>
  st_join(map[,'grid.idx'])|>st_drop_geometry()



data.fit$month <- (data.fit$time.idx-1)%%12 + 1
data.fit$year <- (data.fit$time.idx-1)%/%12 + 2011


n.group <- 5


data.fit$z <- as.vector((data.fit$y>0)+0)

log.ba.agg <- data.fit[data.fit$year<=2022,] %>% group_by(NAME_2) %>%
  summarise(cum_log_ba=sum(log_ba,na.rm=TRUE)) %>%
  mutate(log.ba.idx=inla.group(cum_log_ba,n=n.group,method='quantile',idx.only=T))

cnt.agg <- data.fit[data.fit$year<=2022,] %>% group_by(NAME_2) %>%
  summarise(cum_cnt=sum(y,na.rm=TRUE)) %>%
  mutate(cnt.idx=inla.group(cum_cnt,n=n.group,method='quantile',idx.only=T))

z.agg <- data.fit[data.fit$year<=2022,] %>% group_by(NAME_2) %>%
  summarise(cum_z=sum(z,na.rm=TRUE)) %>%
  mutate(z.idx=inla.group(cum_z,n=n.group,method='quantile',idx.only=T))


data.fit <- merge(data.fit,log.ba.agg[,c("NAME_2",'log.ba.idx')],by='NAME_2')
data.fit <- merge(data.fit, cnt.agg[,c('NAME_2','cnt.idx')], by='NAME_2')
data.fit <- merge(data.fit, z.agg[,c('NAME_2','z.idx')], by='NAME_2')



data.fit <- data.fit[order(data.fit$time.idx,data.fit$grid.idx),]

z <- as.vector((data.fit$y>0)+0)
log.ba <- as.vector(ifelse(data.fit$y>0, data.fit$log_ba, NA))
cnt <- as.vector(ifelse(data.fit$y>0, data.fit$y, NA))


#prepare for prediction
z[which(data.fit$year>2022)] <- NA
log.ba[which(data.fit$year>2022)] <- NA
cnt[which(data.fit$year>2022)] <- NA

n1 <- dim(data.fit)[1]
nothing <- rep(NA, n1)

cntNA <- as.vector(c(cnt,nothing,nothing))
zNA <- as.vector(c(nothing, z, nothing))
baNA.log = as.vector(c(nothing, nothing, log.ba))

Y.log <- matrix(c(cntNA,zNA, baNA.log), ncol=3)

grid.idx.cnt = c(data.fit$grid.idx, nothing, nothing)# fire ignition
grid.idx.z = c(nothing,data.fit$grid.idx, nothing)# BA ind
grid.idx.ba = c(nothing, nothing, data.fit$grid.idx)# BA


time.idx.cnt <- c(data.fit$time.idx, nothing, nothing)# fire ignition
time.idx.z <- c(nothing, data.fit$time.idx, nothing)# BA ind
time.idx.ba  <-  c(nothing, nothing, data.fit$time.idx)# BA

intercept.cnt <- c(rep(1,n1),nothing, nothing)
intercept.z <- c(nothing, rep(1,n1), nothing)
intercept.ba <- c(nothing, nothing, rep(1,n1))


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


score_cnt_grp <- rw.group(data.fit[data.fit$y>0,'score_cnt'],data.fit$score_cnt,40)
score.cnt <- c(score_cnt_grp,nothing, nothing)

score_z_grp <- rw.group(data.fit$score_z,data.fit$score_z,40)
score.z <- c(nothing, score_z_grp, nothing)

score_ba_grp <- rw.group(data.fit[data.fit$y>0,'score_ba'],data.fit$score_ba,40)
score.ba <- c(nothing, nothing, score_ba_grp)

month.idx.cnt <- c(data.fit$month, nothing, nothing)# fire ignition
month.idx.z <- c(nothing, data.fit$month, nothing)# BA ind
month.idx.ba  <-  c(nothing, nothing, data.fit$month)# BA

year.idx <- data.fit$year - 2010
year.idx.cnt <- c(year.idx, nothing, nothing)# fire ignition
year.idx.z <- c(nothing, year.idx, nothing)# BA ind
year.idx.ba  <-  c(nothing, nothing,year.idx)# BA

after.2017 <- as.integer(data.fit$year>2017)
year.af.2017.cnt <- c(after.2017, nothing, nothing)
year.af.2017.z <- c(nothing, after.2017, nothing)
year.af.2017.ba <- c(nothing, nothing,  after.2017)

log.ba.idx <- c(nothing, nothing,  data.fit$log.ba.idx)
cnt.idx <- c(data.fit$cnt.idx, nothing, nothing)
z.idx <- c(nothing, data.fit$z.idx, nothing)


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
  f(cnt.idx, model='rw1')+
  
  intercept.z + year.af.2017.z +
  f(grid.idx.z, model='bym2',graph=g.council)  +
  f(score.z, model='rw1' )+
  f(time.idx.z, model='iid')+ 
  f(year.idx.z, model='iid') + 
  f(z.idx, model='rw1')+
  
  
  intercept.ba + year.af.2017.ba +
  f(grid.idx.ba, copy='grid.idx.z',fixed=F)  +
  f(score.ba, model='rw1' ) +
  f(time.idx.ba, model='iid') + 
  f(year.idx.ba, model='iid') + 
  f(log.ba.idx, model='rw1')

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


n1 <- dim(data.fit)[1]

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



n.samples <- 2000
samples = inla.posterior.sample(n.samples, result = res, seed=1234)



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
  
  return(list('pred.ba'=pred.ba, 'pred.cnt'=pred.cnt, 'pred.z'=pred.z,
              'param.cnt'=param.cnt, 'param.z'=param.z, 'param.ba'=param.ba, 'hyper.ba'=hyper.ba))
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





# 
idx.pred.pois <- 1:n1
idx.pred.z <- (n1+1):(2*n1)
idx.pred.ba <- (2*n1+1):(3*n1)




t1 <- Sys.time()
pred.sp <- post.pred.weibull.hurdle.parallel(samples, idx.pred.pois, idx.pred.z, idx.pred.ba, n.samples=n.samples)
# pred.sp <- post.pred.egp.hurdle.parallel(samples, idx.pred.pois, idx.pred.z, idx.pred.ba, n.samples=n.samples)
# pred.sp <- post.pred.gamma.hurdle.parallel(samples, idx.pred.pois, idx.pred.z, idx.pred.ba, n.samples=n.samples)

t2 <- Sys.time()
print(t2-t1)


pred.cnt <- pred.sp$pred.cnt
pred.ba <- pred.sp$pred.ba
pred.z <- pred.sp$pred.z



data.fit[,'Estimated_Lambda'] <- sapply(pred.cnt,mean)
data.fit$Lower_CNT <- sapply(pred.cnt,quantile,0.025)
data.fit$Upper_CNT <- sapply(pred.cnt,quantile,0.975)
data.fit$Scaling_Residuals <- data.fit$y - data.fit$Estimated_Lambda
summary(data.fit$Scaling_Residuals)

data.fit[,'Estimated_Z'] <- sapply(pred.z,mean)
data.fit$Lower_Z <- sapply(pred.z,quantile,0.025)
data.fit$Upper_Z <- sapply(pred.z,quantile,0.975)

data.fit[,'Estimated_BA'] <- sapply(pred.ba,mean)
data.fit[,'Lower_BA'] <- sapply(pred.ba,quantile,0.025)
data.fit[,'Upper_BA'] <- sapply(pred.ba,quantile,0.975)



crps.ba <- rep(NA, n1)
for (i in 1:n1 ){
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


# 
# 
crps.cnt <- rep(NA, n1)
for (i in 1:n1 ){
  y <- data.fit$y[i]
  if (is.na(y)) y <- 0
  crps.cnt[i] <- crps_sample(
    y,
    pred.ba[[i]],
    method = "edf")
}

df.crps.cnt <- data.frame(crps=crps.cnt,year=data.fit$year,district=data.fit$NAME_1)
df.crps.cnt$label <- 'Train'
df.crps.cnt[which(df.crps.cnt$year<=2022),'label'] <- 'Test'
df.crps.cnt$label <- factor(df.crps.cnt$label, level=c('Train','Test'))

print(df.crps.cnt %>% group_by(district,label) %>% summarize(avg_crps_cnt=mean(crps)),n=50)
print(df.crps.cnt %>% group_by(label) %>% summarize(avg_crps_cnt=mean(crps)),n=50)




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

