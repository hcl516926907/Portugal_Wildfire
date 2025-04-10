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



# load(file.path(dir.data, "XGBoost_Score_Council_2.RData"))
load(file.path(dir.out, 'XGBoost_Score_Council_2.RData'))
# data.fit <- data.fit.council
# rm(data.fit.council)

data.fit[data.fit$y==0,'log_ba'] <- 0
data.fit[data.fit$y==0,'area_ha'] <- 0
data.fit$transformed_ba <- data.fit$area_ha^{1/2}


############################################explorary analysis ##########33333333
df_ba <- data.frame(x = data.fit[data.fit$year!=2023,'log_ba'])
ggplot(df_ba, aes(x = x)) +
  geom_histogram(binwidth = 0.3, fill = "lightblue", color = "black") +
  labs(title = "Histogram of x", x = "x", y = "Frequency")




df <- data.frame(x = log(1+data.fit[data.fit$year!=2023,'area_ha'])) %>%
  mutate(group = if_else(x == 0, "zero", "nonzero"))

scale_factor <- 20

color = "cornflowerblue"
p <- ggplot() +
  # Histogram for Group A (primary axis)
  geom_histogram(data = subset(df, group == "zero"), 
                 aes(x = x), 
                 binwidth = 0.2, 
                 fill = color, 
                 alpha = 1) +
  # Histogram for Group B (transformed counts for secondary axis)
  geom_histogram(data = subset(df, group == "nonzero"), 
                 aes(x = x, y = after_stat(count) * scale_factor), 
                 binwidth = 0.2, 
                 fill = color, 
                 alpha = 1) +
  scale_y_continuous(
    name = "Frequency for zero",
    sec.axis = sec_axis(~ . / scale_factor, name = "Frequency for positive values")
  ) +
  labs(x = "log (1 + Burn Area)") +
  theme_minimal() +
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

png(filename = file.path(dir.out, "Histgram_burn_area.pdf"), width =4000 , height = 2000, res=300)
print(p)
dev.off()




df <- data.frame(x = sqrt(data.fit[data.fit$year!=2023,'y'])) %>%
  mutate(group = if_else(x == 0, "zero", "nonzero"))

scale_factor <- 3

color = "cornflowerblue"
p <- ggplot() +
  # Histogram for Group A (primary axis)
  geom_histogram(data = subset(df, group == "zero"), 
                 aes(x = x), 
                 binwidth = 0.2, 
                 fill = color, 
                 alpha = 1) +
  # Histogram for Group B (transformed counts for secondary axis)
  geom_histogram(data = subset(df, group == "nonzero"), 
                 aes(x = x, y = after_stat(count) * scale_factor), 
                 binwidth = 0.2, 
                 fill = color, 
                 alpha = 1) +
  scale_y_continuous(
    name = "Frequency for zero",
    sec.axis = sec_axis(~ . / scale_factor, name = "Frequency for positive values")
  ) +
  labs(x = expression(sqrt("Fire Count"))) +
  theme_minimal() +
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

png(filename = file.path(dir.out, "Histgram_fire_count.pdf"), width =4000 , height = 2000, res=300)
print(p)
dev.off()










# load(file=file.path(dir.out,'Model_gamma_bym2_pred.RData'))
# load(file=file.path(dir.out,'Model_weibull_bym2_pred.RData'))
load(file=file.path(dir.out,'Model_egp_bym2_pred.RData'))

pred.cnt <- pred.sp$pred.cnt
pred.ba <- pred.sp$pred.ba
pred.z <- pred.sp$pred.z



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
library(xgboost)
library(SHAPforxgboost)
# Calculate SHAP values

rename_var <- function(x){
  mapping <- c("lat.grid" = 'Lat',
               'FWI'= 'FWI',
               "RHumi"="RHumi",
               "month" = "Month",
               "year" = "Year",
               "Temp" = 'Tmp',
               "HVegCov" = "HVegCov",
               "LVegCov"  = "LVegCov" ,
               "Pricp" = "Pricp",
               "lon.grid"  = 'Lon',
               "UComp" = "UComp",
               "HVegLAI" = "HVegLAI",
               "LVegLAI" = "LVegLAI",
               "WindSpeed" = "WindSpeed",
               "DewPoint"   = "DewPoint",
               "LVegTyp_7"  = "VegTyp_Tall_Grass",
               'VComp' = 'VComp',
               'NAME_1_Santarém' = 'Dst_Santarém',
               'LVegTyp_1' = 'VegTyp_Crops',
               'NAME_1_Coimbra' = 'Dst_Coimbra',
               'LVegTyp_16' = "VegTyp_Evergreen_Shrubs",
               "NAME_1_Portalegre" = "Dst_Portalegre",
               'NAME_1_Guarda' = 'Dst_Guarda',
               "NAME_1_Viana do Castelo" = "Dst_Viana_do_Castelo",
               "NAME_1_Évora"  = "Dst_Évora",
               "NAME_1_Porto" = "Dst_Porto",
               "NAME_1_Viseu" = "Dst_Viseu",
               "NAME_1_Lisboa"  = "Dst_Lisboa",
               "NAME_1_Faro" = "Dst_Faro",
               
               "NAME_1_Vila Real" = "Dst_Vila_Real",
               "NAME_1_Leiria" = "Dst_Leiria",
               "NAME_1_Braga" = "Dst_Braga",
               "NAME_1_Setúbal" = "Dst_Setúbal",
               
               "NAME_1_Bragança" = "Dst_Bragança",
               "NAME_1_Castelo Branco" = "Dst_Castelo_Branco",
               "HVegTyp_19"  = "VegTyp_Interrupted_Forest",
               "NAME_1_Beja" = "Dst_Beja",
               
               "HVegTyp_6"  = "VegTyp_Evergreen_Broadleaf_Trees",
               "HVegTyp_18"  = "VegTyp_Mixed_Forest",
               "LVegTyp_11" = "VegTyp_Semidesert"
               )
  return(unname(mapping[x]))
}


X_train_z <- data.rf[,covar.names]
shap_values_z <- shap.values(xgb_model = model.z, X_train = xgb.DMatrix(as.matrix(X_train_z)) )

shap_long_z <- shap.prep(shap_contrib = shap_values_z$shap_score, X_train = as.matrix(X_train_z))
levels(shap_long_z$variable)

levels(shap_long_z$variable) <- rename_var(levels(shap_long_z$variable))

top20_z <- shap_long_z[shap_long_z$variable %in% levels(shap_long_z$variable)[1:10],] %>% droplevels()


my_shap_plot_summary <- function(data_long, x_bound = NULL, dilute = FALSE, scientific = FALSE, 
                                 my_format = NULL,
                                  min_color_bound = '#3399ff', max_color_bound = "#ff5050",
                                 kind = c("sina", "bar")){
  {
    kind <- match.arg(kind)
    if (kind == "bar") {
      imp <- shap.importance(data_long)
      p <- ggplot(imp, aes(x = variable, y = mean_abs_shap)) + 
        geom_bar(stat = "identity", fill = max_color_bound) + 
        coord_flip() + scale_x_discrete(limits = rev(levels(imp[["variable"]]))) + 
        theme_bw() + theme(axis.title.x = element_text(size = 10)) + 
        labs(x = element_blank(), y = "Avg(|SHAP|)")
      return(p)
    }
    if (scientific) {
      label_format = "%.1e"
    }
    else {
      label_format = "%.3f"
    }
    if (!is.null(my_format)) 
      label_format <- my_format
    N_features <- setDT(data_long)[, uniqueN(variable)]
    if (is.null(dilute)) 
      dilute = FALSE
    nrow_X <- nrow(data_long)/N_features
    if (dilute != 0) {
      dilute <- ceiling(min(nrow_X/10, abs(as.numeric(dilute))))
      set.seed(1234)
      data_long <- data_long[sample(nrow(data_long), min(nrow(data_long)/dilute, 
                                                         nrow(data_long)/2))]
    }
    x_bound <- if (is.null(x_bound)) 
      max(abs(data_long$value)) * 1.1
    else as.numeric(abs(x_bound))
    plot1 <- ggplot(data = data_long) + coord_flip(ylim = c(-x_bound, 
                                                            x_bound)) + geom_hline(yintercept = 0) + ggforce::geom_sina(aes(x = variable, 
                                                                                                                            y = value, color = stdfvalue), method = "counts", maxwidth = 0.7, 
                                                                                                                        alpha = 0.7) + geom_text(data = unique(data_long[, c("variable", 
                                                                                                                                                                             "mean_value")]), aes(x = variable, y = -Inf, label = sprintf(label_format, 
                                                                                                                                                                                                                                          mean_value)), size = 5, alpha = 0.7, hjust = -0.2, fontface = "bold", 
                                                                                                                                                 check_overlap = TRUE) + scale_color_gradient(low = min_color_bound, 
                                                                                                                                                                                              high = max_color_bound, breaks = c(0, 1), labels = c(" Low", 
                                                                                                                                                                                                                                                   "High "), guide = guide_colorbar(barwidth = 12, barheight = 0.3)) + 
      theme_bw() + theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(), 
                         legend.position = "bottom",
                         axis.text.x=element_text(size=15),
                         axis.text.y=element_text(size=15),
                         axis.title.x=element_text(size=15),
                         axis.title.y=element_blank(),
                         legend.text = element_text(size=12),
                         legend.title = element_text(size=15,margin = margin(r = 20))) + 
      scale_x_discrete(limits = rev(levels(data_long$variable)), 
                       labels = (rev(levels(data_long$variable)))) +
      labs(color = 'Feature value',y)
    return(plot1)
  }  
  
}
p <- my_shap_plot_summary(top20_z, x_bound  = 2.5, dilute = 20,
                       min_color_bound = '#3399ff',
                       max_color_bound  = "#ff5050") +
  labs(y = "SHAP Value for Fire Indicator XGBoost")
p


png(filename = file.path(dir.out, "SHAP_z.png"), width = 4000, height = 2000, res=300)
print(p)
dev.off()


X_train_ba <- data.rf[data.rf$y>0,covar.names]
shap_values_ba <- shap.values(xgb_model = model.ba, X_train = xgb.DMatrix(as.matrix(X_train_ba)) )

shap_long_ba <- shap.prep(shap_contrib = shap_values_ba$shap_score, X_train = as.matrix(X_train_ba))
levels(shap_long_ba$variable)

levels(shap_long_ba$variable) <- rename_var(levels(shap_long_ba$variable))

top20_ba <- shap_long_ba[shap_long_ba$variable %in% levels(shap_long_ba$variable)[1:10],] %>% droplevels()

p <- my_shap_plot_summary(top20_ba, x_bound  = 1.3, dilute = 20,
                       min_color_bound = '#3399ff',
                       max_color_bound  = "#ff5050") + 
  labs(y = "SHAP Value for Burn Area XGBoost")
p
png(filename = file.path(dir.out, "SHAP_ba.png"), width = 4000, height = 2000, res=300)
print(p)
dev.off()

X_train_cnt <- data.rf[data.rf$y>0,covar.names]
shap_values_cnt <- shap.values(xgb_model = model.cnt, X_train = xgb.DMatrix(as.matrix(X_train_cnt)) )

shap_long_cnt <- shap.prep(shap_contrib = shap_values_cnt$shap_score, X_train = as.matrix(X_train_cnt))
levels(shap_long_cnt$variable)

levels(shap_long_cnt$variable) <- rename_var(levels(shap_long_cnt$variable))

top20_cnt <- shap_long_cnt[shap_long_cnt$variable %in% levels(shap_long_cnt$variable)[1:10],] %>% droplevels()

p <- my_shap_plot_summary(top20_cnt, x_bound  = 2.3, dilute = 20,
                       min_color_bound = '#3399ff',
                       max_color_bound  = "#ff5050")+
  labs(y = "SHAP Value for Fire Count XGBoost") 
p
png(filename = file.path(dir.out, "SHAP_cnt.png"), width = 4000, height = 2000, res=300)
print(p)
dev.off()



################################latent effect###############

samples[[1]]

extract_var_name <- function(x){
  return (as.vector(strsplit(x, split=':')[[1]][1]))
}

effect_name <- sapply (rownames(samples[[1]]$latent), extract_var_name)

idx.score.ba = which(effect_name == 'score.ba')
idx.score.cnt = which(effect_name == 'score.cnt')
idx.score.z = which(effect_name == 'score.z')

idx.year.idx.ba = which(effect_name == 'year.idx.ba')
idx.year.idx.cnt = which(effect_name == 'year.idx.cnt')
idx.year.idx.z = which(effect_name == 'year.idx.z')

idx.intercept.ba = which(effect_name == 'intercept.ba')
idx.intercept.cnt = which(effect_name == 'intercept.cnt')
idx.intercept.z = which(effect_name == 'intercept.z')

idx.grid.idx.ba = which(effect_name == 'grid.idx.ba')
idx.grid.idx.cnt = which(effect_name == 'grid.idx.cnt')
idx.grid.idx.z = which(effect_name == 'grid.idx.z')

idx.grid.idx.dist.ba = which(effect_name == 'grid.idx.dist.ba')
idx.grid.idx.dist.cnt = which(effect_name == 'grid.idx.dist.cnt')
idx.grid.idx.dist.z = which(effect_name == 'grid.idx.dist.z')

test1 = samples[[1]]$latent[idx.pred.pois,]

test2 = samples[[1]]$latent[idx.intercept.cnt,] + samples[[1]]$latent[idx.year.idx.cnt,][data.fit$year-2010] + 
  samples[[1]]$latent[idx.grid.idx.cnt,][data.fit$grid.idx] + samples[[1]]$latent[idx.grid.idx.dist.cnt,][3] + 
  samples[[1]]$latent[idx.score.cnt,][as.integer(as.factor(score_cnt_grp))]

which(abs(test1 - test2)>0.0001)


dist_effect_ba_array = array(dim = c(18,156,1000))
council_effect_ba_array = array(dim = c(278,1000))

dist_effect_cnt_array = array(dim = c(18,156,1000))
council_effect_cnt_array = array(dim = c(278,1000))

dist_effect_z_array = array(dim = c(18,156,1000))
council_effect_z_array = array(dim = c(278,1000))

score_ba_array = array(dim=c(length(idx.score.ba), 1000))
score_cnt_array = array(dim=c(length(idx.score.cnt), 1000))
score_z_array = array(dim=c(length(idx.score.z), 1000))

year_effect_ba_array = array(dim=c(length(idx.year.idx.ba), 1000))
year_effect_cnt_array = array(dim=c(length(idx.year.idx.cnt), 1000))
year_effect_z_array = array(dim=c(length(idx.year.idx.z), 1000))

for (i in 1:length(samples)){
  dist_effect_cnt_array[,,i] =  matrix(samples[[i]]$latent[idx.grid.idx.dist.cnt,][1:(length(idx.grid.idx.dist.cnt)/2)],nrow=18,ncol=156)
  council_effect_cnt_array[,i] = matrix(samples[[i]]$latent[idx.grid.idx.cnt,][1:(length(idx.grid.idx.cnt)/2)],nrow=278,ncol=1)

  dist_effect_ba_array[,,i] =  matrix(samples[[i]]$latent[idx.grid.idx.dist.ba,][1:(length(idx.grid.idx.dist.ba)/2)],nrow=18,ncol=156)
  council_effect_ba_array[,i] = matrix(samples[[i]]$latent[idx.grid.idx.ba,][1:(length(idx.grid.idx.ba)/2)],nrow=278,ncol=1)

  dist_effect_z_array[,,i] =  matrix(samples[[i]]$latent[idx.grid.idx.dist.z,][1:(length(idx.grid.idx.dist.z)/2)],nrow=18,ncol=156)
  council_effect_z_array[,i] = matrix(samples[[i]]$latent[idx.grid.idx.z,][1:(length(idx.grid.idx.z)/2)],nrow=278,ncol=1)
  
  score_ba_array[,i] = samples[[i]]$latent[idx.score.ba]
  score_cnt_array[,i] = samples[[i]]$latent[idx.score.cnt]
  score_z_array[,i] = samples[[i]]$latent[idx.score.z]
  
  year_effect_ba_array[,i] = samples[[i]]$latent[idx.year.idx.ba]
  year_effect_cnt_array[,i] = samples[[i]]$latent[idx.year.idx.cnt]
  year_effect_z_array[,i] = samples[[i]]$latent[idx.year.idx.z]
  }

dist_effect_ba_mean <- apply(dist_effect_ba_array, c(1, 2), mean)
dist_effect_ba_q975 <-  apply(dist_effect_ba_array, c(1, 2), quantile, probs = 0.975)
dist_effect_ba_q025 <-  apply(dist_effect_ba_array, c(1, 2), quantile, probs = 0.025)

council_effect_ba_mean <- apply(council_effect_ba_array, 1, mean)
council_effect_ba_q975 <-  apply(council_effect_ba_array, 1, quantile, probs = 0.975)
council_effect_ba_q025 <-  apply(council_effect_ba_array, 1, quantile, probs = 0.025)


dist_effect_cnt_mean <- apply(dist_effect_cnt_array, c(1, 2), mean)
dist_effect_cnt_q975 <-  apply(dist_effect_cnt_array, c(1, 2), quantile, probs = 0.975)
dist_effect_cnt_q025 <-  apply(dist_effect_cnt_array, c(1, 2), quantile, probs = 0.025)

council_effect_cnt_mean <- apply(council_effect_cnt_array, 1, mean)
council_effect_cnt_q975 <-  apply(council_effect_cnt_array, 1, quantile, probs = 0.975)
council_effect_cnt_q025 <-  apply(council_effect_cnt_array, 1, quantile, probs = 0.025)


dist_effect_z_mean <- apply(dist_effect_z_array, c(1, 2), mean)
dist_effect_z_q975 <-  apply(dist_effect_z_array, c(1, 2), quantile, probs = 0.975)
dist_effect_z_q025 <-  apply(dist_effect_z_array, c(1, 2), quantile, probs = 0.025)

council_effect_z_mean <- apply(council_effect_z_array, 1, mean)
council_effect_z_q975 <-  apply(council_effect_z_array, 1, quantile, probs = 0.975)
council_effect_z_q025 <-  apply(council_effect_z_array, 1, quantile, probs = 0.025)


score_ba_mean = apply(score_ba_array, 1, mean)
score_ba_q975 = apply(score_ba_array, 1, quantile, probs = 0.975)
score_ba_q025 = apply(score_ba_array, 1, quantile, probs = 0.025)

score_cnt_mean = apply(score_cnt_array, 1, mean)
score_cnt_q975 = apply(score_cnt_array, 1, quantile, probs = 0.975)
score_cnt_q025 = apply(score_cnt_array, 1, quantile, probs = 0.025)

score_z_mean = apply(score_z_array, 1, mean)
score_z_q975 = apply(score_z_array, 1, quantile, probs = 0.975)
score_z_q025 = apply(score_z_array, 1, quantile, probs = 0.025)

year_effect_ba_mean = apply(year_effect_ba_array, 1, mean)
year_effect_ba_q975 = apply(year_effect_ba_array, 1, quantile, probs = 0.975)
year_effect_ba_q025 = apply(year_effect_ba_array, 1, quantile, probs = 0.025)

year_effect_cnt_mean = apply(year_effect_cnt_array, 1, mean)
year_effect_cnt_q975 = apply(year_effect_cnt_array, 1, quantile, probs = 0.975)
year_effect_cnt_q025 = apply(year_effect_cnt_array, 1, quantile, probs = 0.025)

year_effect_z_mean = apply(year_effect_z_array, 1, mean)
year_effect_z_q975 = apply(year_effect_z_array, 1, quantile, probs = 0.975)
year_effect_z_q025 = apply(year_effect_z_array, 1, quantile, probs = 0.025)

# 
# colnames(dist_effect_ba_mean) <- 1:156
# colnames(dist_effect_q975) <- 1:156
# colnames(dist_effect_q025) <- 1:156


load(file=file.path(dir.out, 'map_council_district.RData'))

effect_min_value = min(dist_effect_ba_q025,
                       council_effect_ba_q025,
                       dist_effect_cnt_q025,
                       council_effect_cnt_q025,
                       dist_effect_z_q025,
                       council_effect_z_q025
                       )
effect_max_value = max(dist_effect_ba_q975,
                       council_effect_ba_q975,
                       dist_effect_cnt_q975,
                       council_effect_cnt_q975,
                       dist_effect_z_q975,
                       council_effect_z_q975
)
selected_month = (1:12)+ 6*12
csc_mean <- scale_fill_gradientn(
  colours = c("#1a9850", '#91cf60', "#d9ef8b","#fee08b" ,'#fc8d59','#d73027'),
  name = "Mean",
  limits = c(min(dist_effect_ba_mean,
                 council_effect_ba_mean,
                 dist_effect_cnt_mean,
                 council_effect_cnt_mean,
                 dist_effect_z_mean,
                 council_effect_z_mean), 
             
             max(dist_effect_ba_mean,
                 council_effect_ba_mean,
                 dist_effect_cnt_mean,
                 council_effect_cnt_mean,
                 dist_effect_z_mean,
                 council_effect_z_mean
             ))
)

plot_district_effect <- function(df, selected_month, csc){
  
  map_district = map.district
  map_district$value = df[,selected_month[1]]
  map_district$month = month.name[1]
  for (i in 2:12){
    tmp = map.district 
    tmp$value = df[,selected_month[i]]
    tmp$month = month.name[i]
    map_district = rbind(map_district,tmp )
  }
  
  
  map_district$month = factor(map_district$month,levels = month.name )
  
  df_name = deparse(substitute(df))
  if (grepl('ba', df_name)){
    category = 'Burn Area'
  }else if(grepl('cnt', df_name)){
    category = 'Fire Count'
  }else if (grepl('z', df_name)){
    category = 'Fire Indicator'
  }else{
    category = ' '
  }
  title =  paste( "Effect of Districts for ", category,' in 2017', sep='')
  
  p = ggplot(data = map_district) +
    geom_sf(aes(fill = value)) +
    facet_wrap(~ month, ncol = 4) +
    csc +
    labs(title = title,
         fill = "Value") +
    theme(
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 15),
      axis.title.x = element_blank(),
      legend.position = "right",
      legend.text = element_text(size = 10),
      legend.title = element_text(size = 20),        # Increased legend title size
      axis.line = element_line(colour = "black"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.ticks = element_line(color = "black"),
      panel.border = element_blank(),
      panel.background = element_blank(),
      axis.text.y = element_text(size = 15),
      axis.title.y = element_text(size = 15),
      plot.title = element_text(hjust = 0.5, size = 25),  # Center and increase title size
      strip.text = element_text(size = 20)                # Increase facet header size
    )
  return(p)
}


p <- plot_district_effect(dist_effect_cnt_mean, selected_month, csc_mean)
png(filename = file.path(dir.out, "district_effect_cnt_mean.pdf"), width =5000 , height = 6666, res=300)
print(p)
dev.off()



plot_district_effect <- function(df,csc){
  map_countil = map 
  map_countil$value = df
  p <- ggplot(data = map_countil) +
    geom_sf(aes(fill=value)) +
    csc +
    labs(title = "Latent Effect of Councils",
         fill = "Value") +
    theme(
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 15),
      axis.title.x = element_blank(),
      legend.position = "right",
      legend.text = element_text(size = 10),
      legend.title = element_text(size = 20),        # Increased legend title size
      axis.line = element_line(colour = "black"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.ticks = element_line(color = "black"),
      panel.border = element_blank(),
      panel.background = element_blank(),
      axis.text.y = element_text(size = 15),
      axis.title.y = element_text(size = 15),
      plot.title = element_text(hjust = 0.5, size = 25)
    )
  return(p)
}





library(reshape2)

rownames(year_effect_cnt_array) <- years
rownames(year_effect_ba_array)  <- years
rownames(year_effect_z_array)   <- years

# Convert each matrix into a long-format data frame and add a column for the effect type
df_cnt <- melt(year_effect_cnt_array, varnames = c("Year", "Sample"), value.name = "Effect")
df_cnt$Type <- "cnt"

df_ba <- melt(year_effect_ba_array, varnames = c("Year", "Sample"), value.name = "Effect")
df_ba$Type <- "ba"

df_z <- melt(year_effect_z_array, varnames = c("Year", "Sample"), value.name = "Effect")
df_z$Type <- "z"

# Combine the three data frames
df <- rbind(df_cnt, df_ba, df_z)

# Filter out the year 2023
df <- subset(df, Year != "2023")

# Ensure Year is a factor with levels from 2011 to 2022
df$Year <- factor(df$Year, levels = as.character(2011:2022))

ggplot(df, aes(x = Year, y = Effect, fill = Type)) +
  geom_boxplot(position = position_dodge(width = 0.6),width = 0.5) +
  # Assign different colors to the data types (feel free to adjust these colors)
  scale_fill_manual(values = c("cnt" = "#a8dadc",
                               "ba"  = "#ed6a5a",
                               "z"   = "#F4D35E")) +

  labs(y = "Year Effect")+
  theme(
    axis.text.x = element_text( vjust = 0.5, size = 15),
    axis.title.x = element_blank(),
    legend.position = c(0.9, 0.9),
    legend.text = element_text(size = 15),
    legend.title = element_text(size = 20),        # Increased legend title size
    axis.line = element_line(colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.ticks = element_line(color = "black"),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.text.y = element_text(size = 15),
    axis.title.y = element_text(size = 15),
    plot.title = element_text(hjust = 0.5, size = 25))






index = 1:40

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


score_cnt_grp <- rw.group(data.fit[data.fit$y>0,'score_cnt'],data.fit$score_cnt,n.group)

score_z_grp <- rw.group(data.fit$score_z,data.fit$score_z,n.group)

score_ba_grp <- rw.group(data.fit[data.fit$y>0,'score_ba'],data.fit$score_ba,n.group)





rownames(score_cnt_array) <- sort(unique(score_cnt_grp))
rownames(score_ba_array)  <-  sort(unique(score_ba_grp))
rownames(score_z_array)   <-  sort(unique(score_z_grp))

# Convert each matrix into a long-format data frame and add a column for the effect type
df_cnt <- melt(score_cnt_array, varnames = c("index", "Sample"), value.name = "Effect")
df_cnt$Type <- "cnt"

df_ba <- melt(score_ba_array, varnames = c("index", "Sample"), value.name = "Effect")
df_ba$Type <- "ba"

df_z <- melt(score_z_array, varnames = c("index", "Sample"), value.name = "Effect")
df_z$Type <- "z"

# Combine the three data frames


# Ensure Year is a factor with levels from 2011 to 2022
df_cnt$index <- factor(sprintf("%0.2f",df_cnt$index))
df_ba$index <- factor(sprintf("%0.2f",df_ba$index))
df_z$index <- factor(sprintf("%0.2f",df_z$index))


ggplot(df_cnt, aes(x = index, y = Effect)) +
  geom_boxplot(fill = "#a8dadc",position = position_dodge(width = 0.6),width = 0.5) +
  # Assign different colors to the data types (feel free to adjust these colors)
  # scale_fill_manual(values = c("cnt" = "#0D3B66",
  #                              "ba"  = "#ed6a5a",
  #                              "z"   = "#F4D35E")) +
  
  labs(y = "Effect of Score_cnt")+
  ylim(-2, 2)+
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 15),
    axis.title.x = element_blank(),
    legend.position = "None",
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 20),        # Increased legend title size
    axis.line = element_line(colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.ticks = element_line(color = "black"),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.text.y = element_text(size = 15),
    axis.title.y = element_text(size = 15),
    plot.title = element_text(hjust = 0.5, size = 25)
  )



merge(map_district, data.fit[data.fit$time.idx==selected_month[1],c('grid.idx.district','score_cnt')], by = 'grid.idx.district')



map_district = map.district
map_district <- merge(map_district, data.fit[data.fit$time.idx==selected_month[1],c('grid.idx.district','score_cnt')], by = 'grid.idx.district')
map_district$month = month.name[1]
for (i in 2:12){
  tmp = map.district 
  tmp <- merge(tmp, data.fit[data.fit$time.idx==selected_month[i],c('grid.idx.district','score_cnt')], by = 'grid.idx.district')
  tmp$month = month.name[i]
  map_district = rbind(map_district,tmp )
}


map_district$month = factor(map_district$month,levels = month.name )

df_name = deparse(substitute(df))
if (grepl('ba', df_name)){
  category = 'Burn Area'
}else if(grepl('cnt', df_name)){
  category = 'Fire Count'
}else if (grepl('z', df_name)){
  category = 'Fire Indicator'
}else{
  category = ' '
}
title =  paste( "Effect of Districts for ", category,' in 2017', sep='')

p = ggplot(data = map_district) +
  geom_sf(aes(fill = score_cnt)) +
  facet_wrap(~ month, ncol = 4) +
  csc +
  labs(title = title,
       fill = "Value") +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 15),
    axis.title.x = element_blank(),
    legend.position = "right",
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 20),        # Increased legend title size
    axis.line = element_line(colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.ticks = element_line(color = "black"),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.text.y = element_text(size = 15),
    axis.title.y = element_text(size = 15),
    plot.title = element_text(hjust = 0.5, size = 25),  # Center and increase title size
    strip.text = element_text(size = 20)                # Increase facet header size
  )








df.score.cnt <- data.frame('score' = data.fit2$score_cnt)
mat.score.cnt <- c()
mat.score.z <- c()
mat.score.ba <- c()

mat.time.cnt <- c()
mat.time.z <- c()
mat.time.ba <- c()

mat.spat.cnt <- c()
mat.spat.z <- c()
mat.spat.ba <- c()

for (i in 1:length(samples)){
  lp.score.cnt <- as.vector(A1%*%matrix(samples[[i]]$latent[idx.score_1,]))
  lp.score.z <- as.vector(A2%*%matrix(samples[[i]]$latent[idx.score_2,]))
  lp.score.ba <- as.vector(A3%*%matrix(samples[[i]]$latent[idx.score_3,]))
  
  mat.score.cnt <- rbind(mat.score.cnt, lp.score.cnt)
  mat.score.z <- rbind(mat.score.z, lp.score.z)
  mat.score.ba <- rbind(mat.score.ba, lp.score.ba)
  
  mat.time.cnt <- rbind(mat.time.cnt,samples[[i]]$latent[idx.time.idx1,])
  mat.time.z <- rbind(mat.time.z,samples[[i]]$latent[idx.time.idx2,])
  mat.time.ba <- rbind(mat.time.ba,samples[[i]]$latent[idx.time.idx3,])
  
  mat.spat.cnt <- rbind(mat.spat.cnt,samples[[i]]$latent[idx.Intercept1,] + 
                          rowSums(expand.grid(samples[[i]]$latent[idx.idarea1.u,],samples[[i]]$latent[idx.time.idx1,])))
  mat.spat.z <- rbind(mat.spat.z,samples[[i]]$latent[idx.Intercept2,] + 
                        rowSums(expand.grid(samples[[i]]$latent[idx.idarea2.u,],samples[[i]]$latent[idx.time.idx2,])))
  mat.spat.ba <- rbind(mat.spat.ba,samples[[i]]$latent[idx.Intercept3,] + 
                         rowSums(expand.grid(samples[[i]]$latent[idx.idarea3.u,],samples[[i]]$latent[idx.time.idx3,])))
}

data.fit2$effect_score_cnt_mean <- colMeans(mat.score.cnt)
data.fit2$effect_score_cnt_upp <- apply(mat.score.cnt,2 ,quantile,0.975)
data.fit2$effect_score_cnt_low <- apply(mat.score.cnt,2 ,quantile,0.025)

data.fit2$effect_score_z_mean <- colMeans(mat.score.z)
data.fit2$effect_score_z_upp <- apply(mat.score.z,2 ,quantile,0.975)
data.fit2$effect_score_z_low <- apply(mat.score.z,2 ,quantile,0.025)

data.fit2$effect_score_ba_mean <- colMeans(mat.score.ba)
data.fit2$effect_score_ba_upp <- apply(mat.score.ba,2 ,quantile,0.975)
data.fit2$effect_score_ba_low <- apply(mat.score.ba,2 ,quantile,0.025)

data.fit2$effect_spat_cnt_mean <- colMeans(mat.spat.cnt)
data.fit2$effect_spat_cnt_upp <- apply(mat.spat.cnt,2 ,quantile,0.975)
data.fit2$effect_spat_cnt_low <- apply(mat.spat.cnt,2 ,quantile,0.025)

data.fit2$effect_spat_z_mean <- colMeans(mat.spat.z)
data.fit2$effect_spat_z_upp <- apply(mat.spat.z,2 ,quantile,0.975)
data.fit2$effect_spat_z_low <- apply(mat.spat.z,2 ,quantile,0.025)

data.fit2$effect_spat_ba_mean <- colMeans(mat.spat.ba)
data.fit2$effect_spat_ba_upp <- apply(mat.spat.ba,2 ,quantile,0.975)
data.fit2$effect_spat_ba_low <- apply(mat.spat.ba,2 ,quantile,0.025)

ggplot(data.fit2) +
  stat_smooth(aes(x = score_cnt, y = effect_score_cnt_mean), col='#EB5B00') + 
  stat_smooth(aes(x = score_cnt, y = effect_score_cnt_upp),linetype='dashed',col='#EB5B00') + 
  stat_smooth(aes(x = score_cnt, y = effect_score_cnt_low),linetype='dashed',col='#EB5B00')+
  geom_point(aes(x = score_cnt, y = -5),size=0.8) + 
  labs(y='linear predicator')




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

