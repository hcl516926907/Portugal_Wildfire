dir.data <- "/home/pgrad2/2448355h/My_PhD_Project/00_Dataset/Urban_Fires"
dir.out <- '/home/pgrad2/2448355h/My_PhD_Project/01_Output/Urban_Fire_Model'


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
# bru_safe_sp(force = TRUE)

# load(file.path(dir.out, 'data.fit2.score_0.0625.RData'))
# load(file=file.path(dir.out,'Model_weibull_0.0625_pred_sp_200.RData'))

# load(file.path(dir.out, 'data.fit2.score_0.25.RData'))
# load(file=file.path(dir.out,'Model_weibull_0.25_pred_sp_200.RData'))
 
# load(file.path(dir.out, 'data.fit2.score_0.125.RData'))
# load(file=file.path(dir.out,'Model_weibull_0.125_pred_sp_200.RData'))

load(file.path(dir.out, 'data.fit2.pred_0.0625.RData'))
# load(file.path(dir.out, 'data.fit2.pred_0.125.RData'))
load(file=file.path(dir.out,'Model_weibull_0.125_pred_0.0625_pred_sp_200.RData'))
# load(file=file.path(dir.out,'Model_weibull_0.25_pred_0.0625_pred_sp_200.RData'))
# load(file=file.path(dir.out,'Model_weibull_0.125_pred_0.125_pred_sp_200.RData'))


pred.cnt <- pred.sp$pred.cnt
pred.ba <- pred.sp$pred.ba
pred.z <- pred.sp$pred.z



# n.grid <- 2554
# n.grid <- 192
# n.grid <- 681
# n1 <- n.grid*108
# n1 <- n.grid*12
# pred.time <- c(97, 108)
# subset.idx <- (1+(pred.time[1]-1)*n.grid) : (n.grid + (pred.time[2]-1)*n.grid )

n1 <- nrow(data.fit2.pred)

# idx.pred.pois <- 1:n1
# idx.pred.z <- (n1+1):(2*n1)
# idx.pred.ba <- (2*n1+1):(3*n1)
# 
# idx.pred.cnt <- idx.pred.pois[subset.idx]
# idx.pred.z <- idx.pred.z[subset.idx]
# idx.pred.ba <- idx.pred.ba[subset.idx]



crps.ba <- rep(NA, n1)
for (i in 1:n1 ){
  y <- log(data.fit2.pred$area_ha[i])
  if (is.na(y)) y <- 0
  crps.ba[i] <- crps_sample(
    y,
    pred.ba[[i]],
    method = "edf")
}
round(sum(crps.ba)/length(crps.ba),4)



crps.cnt <- rep(NA, n1)
for (i in 1:n1 ){
  y <- data.fit2.pred$y[i]
  if (is.na(y)) y <- 0
  crps.cnt[i] <- crps_sample(
    y,
    pred.ba[[i]],
    method = "edf")
}

round(sum(crps.cnt)/length(crps.cnt),4)


data.fit2.2020 <- data.fit2.pred
data.fit2.2020[is.na(data.fit2.2020$log_ba),'log_ba'] <- 0

data.fit2.2020[,'Estimated_Lambda'] <- sapply(pred.cnt,mean)
data.fit2.2020$Lower_CNT <- sapply(pred.cnt,quantile,0.025)
data.fit2.2020$Upper_CNT <- sapply(pred.cnt,quantile,0.975)
data.fit2.2020$Scaling_Residuals <- data.fit2.2020$y - data.fit2.2020$Estimated_Lambda
summary(data.fit2.2020$Scaling_Residuals)

data.fit2.2020[,'Estimated_Z'] <- sapply(pred.z,mean)
data.fit2.2020$Lower_Z <- sapply(pred.z,quantile,0.025)
data.fit2.2020$Upper_Z <- sapply(pred.z,quantile,0.975)

data.fit2.2020[,'Estimated_BA'] <- sapply(pred.ba,mean)
data.fit2.2020[,'Lower_BA'] <- sapply(pred.ba,quantile,0.025)
data.fit2.2020[,'Upper_BA'] <- sapply(pred.ba,quantile,0.975)




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
sf_districts <- st_as_sf(dist)


conc<-shapefile("/home/pgrad2/2448355h/My_PhD_Project/00_Dataset/Urban_Fires/concelhos.shp")
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
sf_conc <- st_as_sf(conc)

grid.cell.coord <- st_as_sf(data.fit2.2020, coords = c("lon.grid", "lat.grid"), crs = 4326)

# merged_sf <- st_join(grid.cell.coord, sf_districts[,'NAME_1'], join = st_within)
merged_sf1 <- st_join(grid.cell.coord, sf_districts[,'NAME_1'], join = st_nearest_feature)

merged_sf1 %>% group_by(NAME_1) %>% summarize(n_cnt = sum(y)) %>% arrange(desc(n_cnt), ascending=F)

merged_sf1 %>% group_by(NAME_1) %>% summarize(n_ba = sum(log_ba)) %>% arrange(desc(n_ba), ascending=F)


merged_sf1 %>% group_by(month.idx) %>% summarize(n_cnt = sum(y))


#----------------------boxplot of CNT/BA in single district --------------
ggplot() +geom_sf(data=merged_sf1, aes(fill=NAME_1,col=NAME_1)) + 
  geom_sf(data = sf_districts,color = "black", fill = NA)


dist.name <- 'Castelo Branco'
# dist.name <- 'Guarda'
# dist.name <- 'Vila Real'
joint.post.sp <- function(x) Reduce("+", x)
df.dist <- data.frame(month=rep(97:108, each=length(pred.cnt[[1]])))
df.dist$year_label <- 2020
df.dist$month_label <- sprintf("%02d", (df.dist$month-1)%%12+1)
df.dist$time_label <- paste(df.dist$year_label ,df.dist$month_label  ,sep='_')

df.boxplot.true <- data.frame(month=97:108, cnt=rep(NA, 12), log.ba= rep(NA,12))
df.boxplot.true$year_label <- 2020
df.boxplot.true$month_label <- sprintf("%02d", (df.boxplot.true$month-1)%%12+1)
df.boxplot.true$time_label <- paste(df.boxplot.true$year_label ,df.boxplot.true$month_label  ,sep='_')

for (t in 97:108){
  idx <- which(merged_sf1$NAME_1==dist.name & merged_sf1$time.idx==t)
  cnt.true <- sum(data.fit2.2020[idx,'y'])
  ba.true <- sum(data.fit2.2020[idx,'log_ba'])
  
  df.boxplot.true[df.boxplot.true$month==t,c('cnt','log.ba')] <- c(cnt.true, ba.true)
  
  pred.cnt.dist <- joint.post.sp(pred.cnt[idx])
  pred.ba.dist <- joint.post.sp(pred.ba[idx])
  df.dist[df.dist$month==t,'sample_cnt'] <- pred.cnt.dist
  df.dist[df.dist$month==t,'sample_ba'] <- pred.ba.dist
}

ggplot(df.dist[df.dist$month>=85,], aes(x = factor(time_label), y = sample_cnt)) +
  geom_boxplot() + 
  geom_line(data=df.boxplot.true[df.boxplot.true$month>=85,], aes(x=factor(time_label),y=cnt,group = 1), col='red',linewidth=1)+
  labs(x = "Time", y = "fire count", title = paste("Predictive Total Fire Count in", dist.name,"in 2020", sep=' ')) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        axis.title.x = element_blank())  # Rotate x-axis labels for clarity+


ggplot(df.dist[df.dist$month>=85,], aes(x = factor(time_label), y = sample_ba)) +
  geom_boxplot() + 
  geom_line(data=df.boxplot.true[df.boxplot.true$month>=85,], aes(x=factor(time_label),y=log.ba,group = 1), col='red',linewidth=1)+
  labs(x = "Time", y = "log burn area", title = paste("Predictive Total Burn Area in", dist.name,"in 2020", sep=' ')) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        axis.title.x = element_blank())  # Rotate x-axis labels for clarity+



#----------------------------------boxplot by district----------------------

time.idx <- 104

merged_sf2 <- st_join(grid.cell.coord, sf_conc[,'NAME_2'], join = st_nearest_feature)


fire.dist.month <- st_drop_geometry(merged_sf1) %>% group_by(NAME_1, time.idx) %>% 
  summarize(cnt.true = sum(y),
            ba.true = sum(log_ba),
  )

dist.vector <- fire.dist.month %>% group_by(NAME_1) %>% 
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

for (dist in unique(df.dist.month$dist)){
  idx <- which(merged_sf1$NAME_1==dist & merged_sf1$time.idx==time.idx)
  
  pred.cnt.dist <- joint.post.sp(pred.cnt[idx])
  pred.ba.dist <- joint.post.sp(pred.ba[idx])
  df.dist.month[df.dist.month$dist==dist,'sample_cnt'] <- pred.cnt.dist
  df.dist.month[df.dist.month$dist==dist,'sample_ba'] <- pred.ba.dist
}


df.dist.month.true <- merge(df.dist.month.true,fire.dist.month,
                            by.x=c('dist','time.idx'),
                            by.y=c('NAME_1','time.idx'))



ggplot(df.dist.month, aes(x = dist, y = sample_cnt)) +
  geom_boxplot() + 
  geom_line(data=df.dist.month.true, aes(x=dist,y=cnt.true,group = 1), col='red',linewidth=1)+
  # labs(x = "Time", y = "fire count", title = paste("Predictive Total Fire Count in", dist.name,"During 2019-2020", sep=' ')) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        axis.title.x = element_blank())  # Rotate x-axis labels for clarity+


ggplot(df.dist.month, aes(x = factor(dist), y = sample_ba)) +
  geom_boxplot() + 
  geom_line(data=df.dist.month.true, aes(x=factor(dist),y=ba.true,group = 1), col='red',linewidth=1)+
  # labs(x = "Time", y = "fire count", title = paste("Predictive Total Fire Count in", dist.name,"During 2019-2020", sep=' ')) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        axis.title.x = element_blank())  # Rotate x-axis labels for clarity+



####################################################################
ggplot() +geom_sf(data=merged_sf2, aes(fill=NAME_2,col=NAME_2)) + 
  geom_sf(data = sf_conc,color = "black", fill = NA)+
  theme(legend.position="none")



time.idx <- 104



fire.conc.month <- st_drop_geometry(merged_sf2) %>% group_by(NAME_2, time.idx) %>% 
  summarize(cnt.true = sum(y),
            ba.true = sum(log_ba),
  )

conc.vector <- fire.conc.month[fire.conc.month$time.idx==time.idx,] %>% group_by(NAME_2) %>% 
  summarize(cnt.true = sum(cnt.true),
            ba.true = sum(ba.true),
  )%>% arrange(-ba.true) 
conc.vector <- as.character(conc.vector$NAME_2)

df.conc.month <- data.frame(conc=factor(rep(conc.vector, each=length(pred.cnt[[1]]))),
                            time.idx=time.idx, month=(time.idx-1)%%12 + 1)
levels(df.conc.month$conc) <- conc.vector


df.conc.month.true <- data.frame(conc=factor(conc.vector),time.idx=time.idx,
                                 month=(time.idx-1)%%12 + 1)
levels(df.conc.month.true$conc) <- conc.vector

for (conc in unique(df.conc.month$conc)){
  idx <- which(merged_sf2$NAME_2==conc & merged_sf2$time.idx==time.idx)
  
  pred.cnt.conc <- joint.post.sp(pred.cnt[idx])
  pred.ba.conc <- joint.post.sp(pred.ba[idx])
  df.conc.month[df.conc.month$conc==conc,'sample_cnt'] <- pred.cnt.conc
  df.conc.month[df.conc.month$conc==conc,'sample_ba'] <- pred.ba.conc
}


df.conc.month.true <- merge(df.conc.month.true,fire.conc.month,
                            by.x=c('conc','time.idx'),
                            by.y=c('NAME_2','time.idx'))



ggplot(df.conc.month, aes(x = conc, y = sample_cnt)) +
  geom_boxplot() + 
  geom_line(data=df.conc.month.true, aes(x=conc,y=cnt.true,group = 1), col='red',linewidth=1)+
  # labs(x = "Time", y = "fire count", title = paste("Predictive Total Fire Count in", conc.name,"During 2019-2020", sep=' ')) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        axis.title.x = element_blank())  # Rotate x-axis labels for clarity+


ggplot(df.conc.month, aes(x = factor(conc), y = sample_ba)) +
  geom_boxplot() + 
  geom_line(data=df.conc.month.true, aes(x=factor(conc),y=ba.true,group = 1), col='red',linewidth=1)+
  # labs(x = "Time", y = "fire count", title = paste("Predictive Total Fire Count in", conc.name,"During 2019-2020", sep=' ')) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        axis.title.x = element_blank())  # Rotate x-axis labels for clarity+

