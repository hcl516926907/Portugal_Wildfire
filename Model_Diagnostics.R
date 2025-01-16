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
library(ggridges)
# bru_safe_sp(force = TRUE)


# load(file.path(dir.out, 'data.fit.score_0.1.RData'))
load(file.path(dir.data, "Data_For_Fitting_Council.RData"))
data.fit <- data.fit.council
rm(data.fit.council)

data.fit[data.fit$y==0,'log_ba'] <- 0
data.fit[data.fit$y==0,'area_ha'] <- 0
data.fit$transformed_ba <- data.fit$area_ha^{1/2}


# load(file=file.path(dir.out,'Model_gamma_bym2_pred.RData'))
# load(file=file.path(dir.out,'Model_weibull_bym2_pred.RData'))
load(file=file.path(dir.out,'Model_egp_bym2_pred.RData'))

pred.cnt <- pred.sp$pred.cnt
pred.ba <- pred.sp$pred.ba
pred.z <- pred.sp$pred.z

rm(pred.sp)

n1 <- nrow(data.fit)
# n1 <- n.grid*12
# pred.time <- c(97, 108)
# subset.idx <- (1+(pred.time[1]-1)*n.grid) : (n.grid + (pred.time[2]-1)*n.grid )

# n1 <- nrow(data.fit.pred)

idx.pred.pois <- 1:n1
idx.pred.z <- (n1+1):(2*n1)
idx.pred.ba <- (2*n1+1):(3*n1)

# idx.pred.cnt <- idx.pred.pois[subset.idx]
# idx.pred.z <- idx.pred.z[subset.idx]
# idx.pred.ba <- idx.pred.ba[subset.idx]
district <- 'Viseu'
district.idx <- data.fit$NAME_1==district


crps.ba <- rep(NA, nrow(data.fit))
for (i in 1:nrow(data.fit) ){
  y <- data.fit$transformed_ba[i]
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


crps.ba.org <- rep(NA, nrow(data.fit))
for (i in 1:nrow(data.fit) ){
  y <- data.fit$area_ha[i]
  if (is.na(y)) y <- 0
  crps.ba.org[i] <- crps_sample(
    y,
    pred.ba[[i]]^2,
    method = "edf")
}
round(sum(crps.ba.org)/length(crps.ba.org),4)


df.crps.ba.org <- data.frame(crps=crps.ba.org,year=data.fit$year,district=data.fit$NAME_1)
df.crps.ba.org$label <- 'Train'
df.crps.ba.org[which(df.crps.ba$year>2022),'label'] <- 'Test'
df.crps.ba.org$label <- factor(df.crps.ba$label, level=c('Train','Test'))

print(df.crps.ba.org %>% group_by(district,label) %>% summarize(avg_crps_ba_org=mean(crps)),n=50)
print(df.crps.ba.org %>% group_by(label) %>% summarize(avg_crps_ba_org=mean(crps)),n=50)

# 
# 
# crps.t.ba <- rep(NA, nrow(data.fit))
# for (i in 1:nrow(data.fit) ){
#   y <- data.fit$y[i]
#   if (is.na(y)) y <- 0
#   if (y>log(50)){
#     crps.t.ba[i] <- crps_sample(
#       y,
#       pred.ba[[i]],
#       method = "edf")
#   }
# }
# 
# df.crps.t.ba <- data.frame(crps.t=crps.t.ba,year=data.fit$year,district=data.fit$NAME_1)
# df.crps.t.ba$label <- 'Train'
# df.crps.t.ba[which(df.crps.t.ba$year>2022),'label'] <- 'Test'
# df.crps.t.ba$label <- factor(df.crps.t.ba$label, level=c('Train','Test'))
# 
# print(df.crps.t.ba %>% group_by(district,label) %>% summarize(avg_crps_t_ba=mean(crps.t,na.rm=T)),n=50)
# print(df.crps.t.ba %>% group_by(label) %>% summarize(avg_crps_t_ba=mean(crps.t,na.rm=T)),n=50)




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

# 
# crps.t.cnt <- rep(NA, nrow(data.fit))
# for (i in 1:nrow(data.fit) ){
#   y <- data.fit$y[i]
#   if (is.na(y)) y <- 0
#   if (y>5){
#     crps.t.cnt[i] <- crps_sample(
#       y,
#       pred.cnt[[i]],
#       method = "edf")
#   }
# }
# 
# df.crps.t.cnt <- data.frame(crps.t=crps.t.cnt,year=data.fit$year,district=data.fit$NAME_1)
# df.crps.t.cnt$label <- 'Train'
# df.crps.t.cnt[which(df.crps.t.cnt$year>2022),'label'] <- 'Test'
# df.crps.t.cnt$label <- factor(df.crps.t.cnt$label, level=c('Train','Test'))
# 
# print(df.crps.t.cnt %>% group_by(district,label) %>% summarize(avg_crps_t_cnt=mean(crps.t,na.rm=T)),n=50)
# print(df.crps.t.cnt %>% group_by(label) %>% summarize(avg_crps_t_cnt=mean(crps.t,na.rm=T)),n=50)


 


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
data.fit$residual_ba <- data.fit$transformed_ba - data.fit$Estimated_BA
summary(data.fit$residual_ba)
data.fit$cover_ba <- data.fit$Lower_BA <= data.fit$transformed_ba & data.fit$Upper_BA>=data.fit$transformed_ba

data.fit$z <- as.integer(data.fit$y > 0)
data.fit$label <- 'Train'
data.fit[which(data.fit$year > 2022), 'label'] <- 'Test'
data.fit$label <- factor(data.fit$label, levels = c('Train', 'Test'))

library(pROC)
auc(data.fit[data.fit$label=='Train',  'z'], data.fit[data.fit$label=='Train',  'Estimated_Z'], quiet = TRUE)
auc(data.fit[data.fit$label=='Test',  'z'], data.fit[data.fit$label=='Test',  'Estimated_Z'], quiet = TRUE)


# print(df.crps.t.cnt %>% group_by(district,label) %>% summarize(avg_crps_t_cnt=mean(crps.t,na.rm=T)),n=50)
# print(df.crps.t.cnt %>% group_by(label) %>% summarize(avg_crps_t_cnt=mean(crps.t,na.rm=T)),n=50)


weighted_loss <- function(y, preds, type) {
  if (type == 'CNT') {
    thres.cnt <- c(0:15, 17, 19, 21, 23, 25, 30)
    loss <- rep(NA, length(thres.cnt))
    w <- 1 - (1 + (thres.cnt + 1)^2 / 1000)^(-1/4)
    w <- w / w[length(thres.cnt)]
    for (i in 1:length(thres.cnt)) {
      thres <- thres.cnt[i]
      I.true <- y <= thres
      I.pred <- mean(preds <= thres)
      loss[i] <- sum(w[i] * (I.true - I.pred)^2) 
    }
  } else if (type == 'BA') {
    thres.ba <- c(seq(0, 100, 10), 150, 200, 300, 400, 500, 1000, 1500, 2000, 5000, 10000, 20000)
    loss <- rep(NA, length(thres.ba))
    w <- 1 - (1 + (thres.ba + 1) / 1000)^(-1/4)
    w <- w / w[length(thres.ba)]
    for (i in 1:length(thres.ba)) {
      thres <- thres.ba[i]
      I.true <- y <= thres
      I.pred <- mean(preds <= thres)
      loss[i] <- sum(w[i] * (I.true - I.pred)^2) 
    }
  } else {
    loss <- 0
  }
  return(sum(loss)) 
}



loss_cnt <- rep(NA, nrow(data.fit))
loss_ba <- rep(NA, nrow(data.fit))

for (i in 1:nrow(data.fit)) {
  loss_cnt[i] = weighted_loss(data.fit[i, 'y'], pred.cnt[[i]], 'CNT')
  loss_ba[i] = weighted_loss(data.fit[i, 'area_ha'], (pred.ba[[i]])^2, 'BA')
}

data.fit$loss_cnt = loss_cnt
data.fit$loss_ba = loss_ba

data.fit %>% 
  group_by(label) %>% 
  summarize(loss_cnt = sum(loss_cnt))

data.fit %>% 
  group_by(label) %>% 
  summarize(loss_ba = sum(loss_ba))

# print(data.fit %>% group_by(label) %>% summarize(avg_res_cnt=mean((residual_cnt)^2)),n=50)
# print(data.fit %>% group_by(NAME_1,label) %>% summarize(avg_res_cnt=mean((residual_cnt)^2)),n=50)
# 
# print(data.fit %>% group_by(label) %>% summarize(avg_res_ba=mean((residual_ba)^2)),n=50)
# print(data.fit %>% group_by(NAME_1, label) %>% summarize(avg_res_ba=mean((residual_ba)^2)),n=50)
# 
# 
# print(data.fit[data.fit$y>5,] %>% group_by(label) %>% summarize(avg_res_cnt=mean((residual_cnt)^2)),n=50)
# print(data.fit[data.fit$y>5,] %>% group_by(NAME_1,label) %>% summarize(avg_res_cnt=mean((residual_cnt)^2)),n=50)
# 
# print(data.fit[data.fit$transformed_ba>log(50),] %>% group_by(label) %>% summarize(avg_res_ba=mean((residual_ba)^2)),n=50)
# print(data.fit[data.fit$transformed_ba>log(50),] %>% group_by(NAME_1, label) %>% summarize(avg_res_ba=mean((residual_ba)^2)),n=50)
# 
#-----------------------------------------------------------------
# Posterior Predictiv Check

thres.ba <- c(10, 50, 100, 150, 200, 300, 400, 500, 1000,2000)
# ba.break <- c(-1, 100,200,500,10^10)
df.bin.ba <- data.frame(
  threshold = thres.ba,
  prob = sapply(thres.ba, function(thresh) {
    mean(data.fit$area_ha > thresh, na.rm = TRUE)
  })
)
ba.rep.all <- data.frame()
for (i in 1:1000){
  ba.rep <- data.frame(ba_pred = sapply(pred.ba, function(vec) vec[i]^2))
  ba.rep$ba_pred_bin <- data.fit$ba_bins
  ba.rep$year <- data.fit$year
  ba.rep <- ba.rep[ba.rep$year==2023,]
  df.plot <- data.frame(
    threshold = thres.ba,
    prob = sapply(thres.ba, function(thresh) {
      mean(ba.rep$ba_pred > thresh, na.rm = TRUE)
    })
  )
  ba.rep.all <- rbind(ba.rep.all, df.plot)
  print(i)
}
p <- ggplot(ba.rep.all, aes(x = prob, y = as.factor(threshold))) +
  geom_density_ridges(stat = "binline",scale = 2, rel_min_height = 0.001,binwidth=0.001,
                      fill = "darkgoldenrod2") +
  geom_segment(data = df.bin.ba,
               aes(x = prob, y = as.factor(threshold), xend = prob,
                   yend = as.numeric(as.factor(threshold))+1 + c(rep(0,length(threshold)-1),1)),
               color = "red",
               linetype = "dashed",
               size=1,show.legend=F)+
  theme_minimal() +
  # xlim(-5, 1000)+
  labs( x = "Exceedance Probability", y = "Burn Area (ha)")+
  theme(
    axis.line = element_line(colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.ticks = element_line(color = "black"),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.text.x=element_text(size=15),
    axis.text.y=element_text(size=15),
    axis.title.x=element_text(size=15),
    axis.title.y=element_text(size=15)) 
print(p)

png(filename = file.path(dir.out, "Threshold_Exceedance_BA.png"), width = 3000, height = 1500, res=300)
print(p)
dev.off()




thres.cnt<- c(5,7, 9, 11,13,15,20)
# ba.break <- c(-1, 100,200,500,10^10)
df.bin.cnt <- data.frame(
  threshold = thres.cnt,
  prob = sapply(thres.cnt, function(thresh) {
    mean(data.fit$y > thresh, na.rm = TRUE)
  })
)
# geom_segment(data = df.bin.org, 
#              aes(x = ba, y = ba_bins, xend = ba, http://130.209.66.182:8787/graphics/38ce3b99-b76b-4fcf-9bf8-5a37dddda197.png
#                  yend = as.numeric(ba_bins)+1),
#              color = "red", 
#              linetype = "dashed",
#              size=1,show.legend=F)+

cnt.rep.all <- data.frame()
for (i in 1:1000){
  cnt.rep <- data.frame(cnt_pred = sapply(pred.cnt, function(vec) vec[i]))
  cnt.rep$year <- data.fit$year
  cnt.rep <- cnt.rep[cnt.rep$year==2023,]
  df.plot.cnt <- data.frame(
    threshold = thres.cnt,
    prob = sapply(thres.cnt, function(thresh) {
      mean(cnt.rep$cnt_pred > thresh, na.rm = TRUE)
    })
  )
  cnt.rep.all <- rbind(cnt.rep.all, df.plot.cnt)
  print(i)
}

p <- ggplot(cnt.rep.all, aes(x = prob, y = as.factor(threshold))) +
  geom_density_ridges(stat = "binline",scale = 2, rel_min_height = 0.001,binwidth=0.0003,
                      fill = "darkgoldenrod2") +
  geom_segment(data = df.bin.cnt,
               aes(x = prob, y = as.factor(threshold), xend = prob,
                   yend = as.numeric(as.factor(threshold))+1 + c(rep(0,length(threshold)-1),1)),
               color = "red",
               linetype = "dashed",
               size=1,show.legend=F)+
  scale_fill_viridis_c(name = "Temp. [F]", option = "C") +
  theme_minimal() +
  # scale_x_continuous(breaks = seq(0,0.2,0.005) )+
  labs( x = "Exceedance Probability", y = "Fire Count")+
  theme(
    axis.line = element_line(colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.ticks = element_line(color = "black"),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.text.x=element_text(size=15),
    axis.text.y=element_text(size=15),
    axis.title.x=element_text(size=15),
    axis.title.y=element_text(size=15)) 
print(p)

png(filename = file.path(dir.out, "Threshold_Exceedance_CNT.png"), width = 3000, height = 1500, res=300)
print(p)
dev.off()



#-----------------------------------------------------------------
dist<-shapefile(file.path(dir.data,'shapefile', "distritos.shp"))
# Municipalities
conc<-shapefile(file.path(dir.data, "shapefile","concelhos.shp"))

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


merged_sf1 <- st_as_sf(data.fit, coords = c("lon.grid", "lat.grid"), crs = 4326)




#----------------------boxplot of CNT/BA in single district --------------
ggplot() +geom_sf(data=merged_sf1, aes(fill=NAME_1,col=NAME_1)) + 
  geom_sf(data = sf_districts,color = "black", fill = NA)


# dist.name <- 'Castelo Branco'
# dist.name <- 'Guarda'
# dist.name <- 'Vila Real'
# dist.name <-  'Viseu'
# dist.name <- 'Porto'
joint.post.sp <- function(x) Reduce("+", x)
df.dist <- data.frame(time.idx=rep(1:156, each=length(pred.cnt[[1]])))
df.dist$year_label <- rep(2011:2023,each=length(pred.cnt[[1]])*12)
df.dist$month_label <- sprintf("%02d", (df.dist$time.idx-1)%%12+1)
df.dist$time_label <- paste(df.dist$year_label ,df.dist$month_label  ,sep='_')
df.dist$label <- ifelse(df.dist$ year_label=='2023','test','train')

df.boxplot.true <- data.frame(time.idx=1:156, cnt=rep(NA, 12*13), log.ba= rep(NA,12*13))
df.boxplot.true$year_label <- rep(2011:2023,each=12)
df.boxplot.true$month_label <- sprintf("%02d", (df.boxplot.true$time.idx-1)%%12+1)
df.boxplot.true$time_label <- paste(df.boxplot.true$year_label ,df.boxplot.true$month_label  ,sep='_')

for (t in 1:156){
  idx <- which(merged_sf1$time.idx==t)
  cnt.true <- sum(data.fit[idx,'y'])
  ba.true <- sum((data.fit[idx,'area_ha']),na.rm=T)
  
  df.boxplot.true[df.boxplot.true$time.idx==t,c('cnt','ba')] <- c(sqrt(cnt.true), log(1+ba.true))
  
  pred.cnt.dist <- joint.post.sp(pred.cnt[idx])
  pred.ba.dist <- joint.post.sp(lapply(pred.ba[idx], function(x) x^2))
  df.dist[df.dist$time.idx==t,'sample_cnt'] <- sqrt(pred.cnt.dist)
  df.dist[df.dist$time.idx==t,'sample_ba'] <- log(1 + pred.ba.dist )
  
}

subset.label <- c(paste("2017_",sprintf("%02d", 1:12),sep=''),
                  paste("2023_",sprintf("%02d", 1:12),sep=''))
p <- ggplot(df.dist[df.dist$time_label %in% subset.label,], aes(x = factor(time_label), y = sample_cnt)) +
  geom_boxplot( aes(fill=label)) + 
  scale_fill_manual(values=c('train'= 'cornflowerblue', 'test'='darkgoldenrod2'))+
  geom_line(data=df.boxplot.true[df.boxplot.true$time_label %in% subset.label ,], aes(x=factor(time_label),y=cnt,group = 1), col='red',linewidth=0.8)+
  geom_point(data=df.boxplot.true[df.boxplot.true$time_label %in% subset.label ,], aes(x=factor(time_label),y=cnt,group = 1), col='red', size=2.5)+
  # geom_vline(xintercept=12.5, linetype="dashed", color = "bisque4",linewidth=1.2)+
  annotate("rect", xmin = 12.5, xmax = Inf, ymin = 0, ymax = Inf,
           fill = "cornsilk3", alpha = 0.2) +  # Define color and transparency
  
  labs(x = "Time", y = expression(sqrt("Fire Count"))) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1,size=15),
        axis.title.x = element_blank(),
        legend.position="none",
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks = element_line(color = "black"),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text.y=element_text(size=15),
        axis.title.y=element_text(size=15))  # Rotate x-axis labels for clarity+

print(p)

png(filename = file.path(dir.out, "Posterior_Prediction_CNT.png"), width = 4000, height = 2000, res=300)
print(p)
dev.off()


p <- ggplot(df.dist[df.dist$time_label %in% subset.label,], aes(x = factor(time_label), y = sample_ba)) +
  geom_boxplot( aes(fill=label)) + 
  scale_fill_manual(values=c('train'= 'cornflowerblue', 'test'='darkgoldenrod2'))+
  geom_line(data=df.boxplot.true[df.boxplot.true$time_label %in% subset.label ,], aes(x=factor(time_label),y=ba,group = 1), col='red',linewidth=0.8)+
  geom_point(data=df.boxplot.true[df.boxplot.true$time_label %in% subset.label ,], aes(x=factor(time_label),y=ba,group = 1), col='red', size=2.5)+
  # geom_vline(xintercept=12.5, linetype="dashed", color = "bisque4",linewidth=1.2)+
  annotate("rect", xmin = 12.5, xmax = Inf, ymin = 0, ymax = Inf,
           fill = "cornsilk3", alpha = 0.2) +  # Define color and transparency
  labs(x = "Time", y = "log (1 + Burn Area)") +
  # scale_y_continuous(limits = c(0, 500))+
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1,size=15),
        axis.title.x = element_blank(),
        legend.position="none",
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks = element_line(color = "black"),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text.y=element_text(size=15),
        axis.title.y=element_text(size=15))  # Rotate x-axis labels for clarity+
print(p)

png(filename = file.path(dir.out, "Posterior_Prediction_BA.png"), width = 4000, height = 2000, res=300)
print(p)
dev.off()





###################################SHAP value##########################################
load(file.path(dir.out, 'XGBoost_Models.RData'))
library(shapr)



explainer <- shapr(as.matrix(data.fit[data.fit$year <= 2022, covar.names]), model.z, n_combinations = 10000)
#> The specified model provides feature classes that are NA. The classes of data are taken as the truth.

# Specifying the phi_0, i.e. the expected prediction without any features
p <- mean(as.matrix(data.fit[data.fit$year <= 2022, 'z']))

# Computing the actual Shapley values with kernelSHAP accounting for feature dependence using
# the empirical (conditional) distribution approach with bandwidth parameter sigma = 0.1 (default)
explanation <- explain(
  as.matrix(data.fit[1:10, covar.names]),
  approach = "empirical",
  explainer = explainer,
  prediction_zero = p
)

# Printing the Shapley values for the test data.
# For more information about the interpretation of the values in the table, see ?shapr::explain.
print(explanation$dt)


# Plot the resulting explanations for observations 1 and 6
plot(explanation, plot_phi0 = FALSE, index_x_test = c(1, 6))






















#----------------------------------boxplot by district----------------------

time.idx <- 82
# time.idx <-  152

fire.dist.month <- st_drop_geometry(merged_sf1) %>% group_by(NAME_1, time.idx) %>% 
  summarize(cnt.true = sum(y),
            ba.true = sum(area_ha,na.rm=T),
  )
fire.dist.month$ba.true = ifelse(fire.dist.month$ba.true>0,log(fire.dist.month$ba.true),0)


dist.vector <- data.fit %>% group_by(NAME_1) %>% 
  summarize(lat = mean(lat.grid)
  )%>% arrange(-lat)

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
  pred.ba.dist <- joint.post.sp(lapply(pred.ba[idx], function(x) x^2))
  df.dist.month[df.dist.month$dist==dist,'sample_cnt'] <- pred.cnt.dist
  df.dist.month[df.dist.month$dist==dist,'sample_ba'] <- ifelse(pred.ba.dist>1,log(pred.ba.dist),0)
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
  # scale_y_continuous(limits = c(0, 100))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        axis.title.x = element_blank())  # Rotate x-axis labels for clarity+



####################################################################
ggplot() +geom_sf(data=merged_sf2, aes(fill=NAME_2,col=NAME_2)) + 
  geom_sf(data = sf_conc,color = "black", fill = NA)+
  theme(legend.position="none")



time.idx <- 104



fire.conc.month <- st_drop_geometry(merged_sf2) %>% group_by(NAME_2, time.idx) %>% 
  summarize(cnt.true = sum(y),
            ba.true = sum(transformed_ba),
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

