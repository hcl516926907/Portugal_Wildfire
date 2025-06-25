dir.data <- "/home/pgrad2/2448355h/My_PhD_Project/00_Dataset/Porgual_Wildfire"
dir.out <- '/home/pgrad2/2448355h/My_PhD_Project/01_Output/Portugal_Wildfire'


library(INLA)
library(tidyr)

inla.setOption(scale.model.default = TRUE)

library(mgcv)
library(ggplot2)
library(sf)
library(dplyr)
library(RColorBrewer)
library(scoringRules)
library(ggridges)
library(data.table)


load(file.path(dir.out, 'dataset_perf_evaluate.RData'))


data.fit$transformed_ba <- data.fit$area_ha^{1/2}


########################## explorary analysis #############################
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
    axis.text.x=element_text(size=30),
    axis.text.y=element_text(size=30),
    axis.title.x=element_text(size=30),
    axis.title.y=element_text(size=25)) 

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
    axis.text.x=element_text(size=30),
    axis.text.y=element_text(size=30),
    axis.title.x=element_text(size=30),
    axis.title.y=element_text(size=25)) 

png(filename = file.path(dir.out, "Histgram_fire_count.pdf"), width =4000 , height = 2000, res=300)
print(p)
dev.off()




####################################################################
# Load posterior predictive samples 
####################################################################

# load(file=file.path(dir.out,'Model_gamma_bym2_pred_h1.RData'))
load(file=file.path(dir.out,'Model_egp_bym2_pred_h1.RData'))
# load(file=file.path(dir.out,'Model_weibull_bym2_pred_h1.RData'))
# load(file=file.path(dir.out,'Model_egp_bym2_pred_h1_noxgb.RData'))

pred.cnt <- pred.sp$pred.cnt
pred.ba <- pred.sp$pred.ba
pred.z <- pred.sp$pred.z



n1 <- nrow(data.fit)


idx.pred.pois <- 1:n1
idx.pred.z <- (n1+1):(2*n1)
idx.pred.ba <- (2*n1+1):(3*n1)


district <- 'Viseu'
district.idx <- data.fit$NAME_1==district

################################# CRPS ######################################
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

df.crps.ba <- data.frame(crps=crps.ba,year=data.fit$year,district=data.fit$NAME_1,transformed_ba= data.fit$transformed_ba)
df.crps.ba$label <- 'Train'
df.crps.ba[which(df.crps.ba$year>2022),'label'] <- 'Test'
df.crps.ba$label <- factor(df.crps.ba$label, level=c('Train','Test'))


print(df.crps.ba %>% filter(transformed_ba>0)%>% group_by(label) %>% summarize(avg_crps_ba=mean(crps)),n=50)



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


################################## AUC ##############################
library(pROC)
auc(data.fit[data.fit$label=='Train',  'z'], data.fit[data.fit$label=='Train',  'Estimated_Z'], quiet = TRUE)
auc(data.fit[data.fit$label=='Test',  'z'], data.fit[data.fit$label=='Test',  'Estimated_Z'], quiet = TRUE)


############################### weighted score ########################
weighted_loss <- function(y, preds, type) {
  if (type == 'CNT') {
    # thres.cnt <- c(0:15, 17, 19, 21, 23, 25, 30)
    thres.cnt <- c(0:10,15,20,25,30)
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
    # thres.ba <- c(seq(0, 100, 10), 150, 200, 300, 400, 500, 1000, 1500, 2000, 5000, 10000, 20000)
    thres.ba <- c(seq(0, 100, 20), 200, 300, 400, 500, 1000, 2000,5000,10000,20000,50000)
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



equal_weighted_loss <- function(y, preds, type) {
  if (type == 'CNT') {
    # thres.cnt <- c(0:15, 17, 19, 21, 23, 25, 30)
    thres.cnt <- c(0:10,15,20,25,30)
    loss <- rep(NA, length(thres.cnt))
    for (i in 1:length(thres.cnt)) {
      thres <- thres.cnt[i]
      I.true <- y <= thres
      I.pred <- mean(preds <= thres)
      loss[i] <- sum((I.true - I.pred)^2) 
    }
  } else if (type == 'BA') {
    # thres.ba <- c(seq(0, 100, 10), 150, 200, 300, 400, 500, 1000, 1500, 2000, 5000, 10000, 20000)
    thres.ba <- c(seq(0, 100, 20), 200, 300, 400, 500, 1000, 2000,5000,10000,20000,50000)
    loss <- rep(NA, length(thres.ba))
    for (i in 1:length(thres.ba)) {
      thres <- thres.ba[i]
      I.true <- y <= thres
      I.pred <- mean(preds <= thres)
      loss[i] <- sum( (I.true - I.pred)^2) 
    }
  } else {
    loss <- 0
  }
  return(sum(loss)) 
}


loss_cnt_uw <- rep(NA, nrow(data.fit))
loss_ba_uw <- rep(NA, nrow(data.fit))

for (i in 1:nrow(data.fit)) {
  loss_cnt_uw[i] = equal_weighted_loss(data.fit[i, 'y'], pred.cnt[[i]], 'CNT')
  loss_ba_uw[i] = equal_weighted_loss(data.fit[i, 'area_ha'], (pred.ba[[i]])^2, 'BA')
}

data.fit$loss_cnt_uw = loss_cnt_uw
data.fit$loss_ba_uw = loss_ba_uw

data.fit %>% 
  group_by(label) %>% 
  summarize(loss_cnt = sum(loss_cnt_uw))

data.fit %>% 
  group_by(label) %>% 
  summarize(loss_ba = sum(loss_ba_uw))


########################### Posterior Predictiv Check ######################

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
    axis.text.x=element_text(size=20),
    axis.text.y=element_text(size=20),
    axis.title.x=element_text(size=20),
    axis.title.y=element_text(size=20)) 
print(p)

png(filename = file.path(dir.out, "Threshold_Exceedance_BA.pdf"), width = 3000, height = 1500, res=300)
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
    axis.text.x=element_text(size=20),
    axis.text.y=element_text(size=20),
    axis.title.x=element_text(size=20),
    axis.title.y=element_text(size=20)) 
print(p)

png(filename = file.path(dir.out, "Threshold_Exceedance_CNT.pdf"), width = 3000, height = 1500, res=300)
print(p)
dev.off()



################## Raw fire count and burnt area in spatial #####################

dist<-read_sf(file.path(dir.data,'shapefile', "distritos.shp"))
# Municipalities
conc<-read_sf(file.path(dir.data, "shapefile","concelhos.shp"))

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


merged_sf1 <- st_as_sf(data.fit, coords = c("lon.grid", "lat.grid"), crs = 4326)


df_wildfire_analysis = data.fit%>% group_by(NAME_2) %>% summarize(mean_cnt = sqrt(mean(y)),
                                           mean_ba = log(1+mean(area_ha))) |> as.data.frame()

df_wildfire = sf_conc %>% left_join(df_wildfire_analysis, by='NAME_2')



csc_cnt <- scale_fill_gradientn(
  colours = c("#1a9850", "#91cf60", "#d9ef8b", "#fee08b", "#fc8d59", "#d73027"),
  name   = expression(sqrt(Count_mean)),
  limits = range(df_wildfire$mean_cnt)
)


csc_ba <- scale_fill_gradientn(
  colours = c("#1a9850", '#91cf60', "#d9ef8b","#fee08b" ,'#fc8d59','#d73027'),
  name =  expression(log(1+BurnArea_mean)),
  limits = c(min(df_wildfire$mean_ba),max(df_wildfire$mean_ba))
)



p <- ggplot() + geom_sf(data=df_wildfire, aes(fill=mean_cnt,col=mean_cnt))+
  csc_cnt +  guides(colour = "none") + 
  labs(title= expression(sqrt(Count_mean)))+
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 50),
    axis.title.x = element_blank(),
    legend.position = "right",
    legend.key.size = unit(4, "lines"),
    legend.text = element_text(size = 50),
    legend.title = element_blank(),        # Increased legend title size
    axis.line = element_line(colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.ticks = element_line(color = "black"),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.text.y = element_text(size = 50),
    axis.title.y = element_text(size = 50),
    plot.title = element_text(hjust = 0.5, size = 60),  # Center and increase title size
    strip.text = element_text(size = 50)                # Increase facet header size
  )

png(filename = file.path(dir.out, "raw_fire_count_spatial.pdf"), width =5000 , height = 6666, res=300)
print(p)
dev.off()


p <- ggplot() + geom_sf(data=df_wildfire, aes(fill=mean_ba,col=mean_ba))+
csc_ba +  guides(colour = "none") +  
  labs(title= expression(log(1+BurnArea_mean)))+
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 50),
    axis.title.x = element_blank(),
    legend.position = "right",
    legend.key.size = unit(4, "lines"),
    legend.text = element_text(size = 50),
    legend.title = element_blank(),        # Increased legend title size
    axis.line = element_line(colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.ticks = element_line(color = "black"),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.text.y = element_text(size = 50),
    axis.title.y = element_text(size = 50),
    plot.title = element_text(hjust = 0.5, size = 60),  # Center and increase title size
    strip.text = element_text(size = 50)                # Increase facet header size
  )


png(filename = file.path(dir.out, "raw_burn_area_spatial.pdf"), width =5000 , height = 6666, res=300)
print(p)
dev.off()


################## Raw fire count and burnt area by month #####################

df_wildfire_analysis_2 = data.fit%>% group_by(month) %>% summarize(mean_cnt = 5*sqrt(mean(y)),
                                                                  mean_ba = log(1+mean(area_ha))) |> as.data.frame()


df_long <- df_wildfire_analysis_2 %>%
  pivot_longer(
    cols      = c(mean_cnt, mean_ba),
    names_to  = "metric",
    values_to = "value"
  ) %>%
  mutate(
    month = factor(month,
                   levels = 1:12,
                   labels = month.abb),
    # **here** we relabel the factor levels to the exact plotmath calls
    metric = factor(metric,
                    levels = c("mean_cnt","mean_ba"),
                    labels = c("sqrt(Count_mean)",   # → \bar(n)
                               "log(1+BurnArea_mean)"   # → \bar(A)
                    )
    )
  )

p <- ggplot(df_long, aes(x = month, y = value)) +
  geom_col(fill = "cornflowerblue") +
  facet_wrap(~ metric,
             scales   = "free_y",
             ncol     = 1,
             labeller = label_parsed   # <— parse those "bar(n)" strings!
  ) +
  labs(
    x     = "Month",
    y     = NULL,
  ) +
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
    axis.title.y=element_text(size=15),
    strip.text = element_text(size = 15)) 

png(filename = file.path(dir.out, "Temporal_analysis.pdf"), width = 3000, height = 1500, res=300)
print(p)
dev.off()


############# boxplot of forecasted fire count and burnt area over Portugal #####################
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

for (t in (min(data.fit$time.idx)): (max(data.fit$time.idx))){
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
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1,size=25),
        axis.title.x = element_blank(),
        legend.position="none",
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks = element_line(color = "black"),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text.y=element_text(size=25),
        axis.title.y=element_text(size=25))  # Rotate x-axis labels for clarity+

print(p)

png(filename = file.path(dir.out, "Posterior_Prediction_CNT.pdf"), width = 4000, height = 2000, res=300)
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
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1,size=25),
        axis.title.x = element_blank(),
        legend.position="none",
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks = element_line(color = "black"),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text.y=element_text(size=25),
        axis.title.y=element_text(size=25))  # Rotate x-axis labels for clarity+
print(p)

png(filename = file.path(dir.out, "Posterior_Prediction_BA.pdf"), width = 4000, height = 2000, res=300)
print(p)
dev.off()





########################### SHAP value ###################################

load(file.path(dir.out,'AutoRegressive_XGBoost_Hyperparameters_CNT.RData'))
load(file.path(dir.out,'AutoRegressive_XGBoost_Hyperparameters_BA.RData'))
library(xgboost)
library(SHAPforxgboost)
# Calculate SHAP values

rename_var <- function(x){
  mapping <- c("lat.grid" = 'Lat',
               'FWI'= 'FWI',
               "RHumi"="RHumi",
               "month" = "Month",
               "year" = "Year",
               "Temp" = 'Temp',
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
               "LVegTyp_11" = "VegTyp_Semidesert",
               
               "hist_cnt_h1_3m" = 'conc_fc_hist_3',
               "hist_cnt_h1_5m" = 'conc_fc_hist_5',
               "cnt_ma_24" = "conc_fc_ma_24",
               'dist_cnt_lag0' = 'dist_fc_lag_1',
               'dist_hist_cnt_h1_1m' = 'dist_fc_hist_1',
               
               'hist_ba_h1_3m' = 'conc_ba_hist_3',
               'dist_ba_lag0' = 'dist_ba_lag_1',
               'hist_ba_h1_1m' = 'conc_ba_hist_1',
               'ba_ma_3' = 'conc_ba_ma_3',
               'hist_ba_h1_5m' = 'conc_ba_hist_5',
               'dist_ba_lag1' = 'dist_ba_lag_2',
               'ba_lag0' = 'conc_ba_lag_1',
               'month_cos' = 'month_cos',
               'dist_hist_ba_h1_1m' = 'dist_ba_hist_1'
               
               )
  return(unname(mapping[x]))
}

final_model <- function(data, covar.names, target.name, tuning_res){
  
  label <- data[,target.name]
  features <- data[, covar.names]
  dtrain <- xgb.DMatrix(data = as.matrix(features), label = label)
  
  model <- xgb.train(
    params = tuning_res$best_params,
    data = dtrain,
    nrounds = tuning_res$best_nrounds ,
    verbose = 0  )
  return(model)
}
w = 9
ba.covar.names.h1 <- c( paste('ba_lag',0:(w-1),sep=''),
                        paste('dist_ba_lag',0:(w-1),sep=''),
                        'ba_ma_3', 'ba_ma_6','ba_ma_9','ba_ma_12','ba_ma_24','ba_ma_36',
                        'dist_ba_ma_3', 'dist_ba_ma_6','dist_ba_ma_9','dist_ba_ma_12','dist_ba_ma_24','dist_ba_ma_36',
                        'hist_ba_h1_1m','hist_ba_h1_3m','hist_ba_h1_5m',
                        'dist_hist_ba_h1_1m','dist_hist_ba_h1_3m','dist_hist_ba_h1_5m',
                        
                        'lon.grid','lat.grid','year',
                        'month_sin','month_cos',
                        'FWI','HVegCov','HVegLAI','LVegCov','LVegLAI',
                        'Pricp','RHumi','Temp','UComp','VComp'
)
# df[,ba.covar.names.h1]

cnt.covar.names.h1 <- c( paste('cnt_lag',0:(w-1),sep=''),
                         paste('dist_cnt_lag',0:(w-1),sep=''),
                         'cnt_ma_3', 'cnt_ma_6','cnt_ma_9','cnt_ma_12','cnt_ma_24','cnt_ma_36',
                         'dist_cnt_ma_3', 'dist_cnt_ma_6','dist_cnt_ma_9','dist_cnt_ma_12','dist_cnt_ma_24','dist_cnt_ma_36',
                         'hist_cnt_h1_1m','hist_cnt_h1_3m','hist_cnt_h1_5m',
                         'dist_hist_cnt_h1_1m','dist_hist_cnt_h1_3m','dist_hist_cnt_h1_5m',
                         
                         'lon.grid','lat.grid','year',
                         'month_sin','month_cos',
                         'FWI','HVegCov','HVegLAI','LVegCov','LVegLAI',
                         'Pricp','RHumi','Temp','UComp','VComp'
)

# Optionally, inspect your final model or its performance.

data.fit.train <- data.fit[data.fit$year<2023,]

set.seed(123)
model_cnt_h1 <- final_model(as.data.frame(data.fit.train), cnt.covar.names.h1, 'xgb_cnt_h1',res_cnt_h1)
model_ba_h1 <- final_model(as.data.frame(data.fit.train), ba.covar.names.h1, 'xgb_ba_h1', res_ba_h1)

X_train_cnt = data.fit.train[,cnt.covar.names.h1]
shap_values_cnt <- shap.values(xgb_model = model_cnt_h1, X_train = xgb.DMatrix(as.matrix(X_train_cnt)) )


shap_long_cnt <- shap.prep(shap_contrib = shap_values_cnt$shap_score, X_train = as.matrix(X_train_cnt))
levels(shap_long_cnt$variable)

# levels(shap_long_cnt$variable) <- rename_var(levels(shap_long_cnt$variable))

top20_cnt <- shap_long_cnt[shap_long_cnt$variable %in% levels(shap_long_cnt$variable)[1:10],] %>% droplevels()


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
                                                                                                                                                                                                                                          mean_value)), size = 7, alpha = 0.7, hjust = -0.2, fontface = "bold", 
                                                                                                                                                 check_overlap = TRUE) + scale_color_gradient(low = min_color_bound, 
                                                                                                                                                                                              high = max_color_bound, breaks = c(0, 1), labels = c(" Low", 
                                                                                                                                                                                                                                                   "High "), guide = guide_colorbar(barwidth = 12, barheight = 0.3)) + 
      theme_bw() + theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(), 
                         legend.position = "bottom",
                         axis.text.x=element_text(size=25),
                         axis.text.y=element_text(size=25),
                         axis.title.x=element_text(size=25),
                         axis.title.y=element_blank(),
                         legend.text = element_text(size=25),
                         legend.title = element_text(size=25,margin = margin(r = 20))) + 
      scale_x_discrete(limits = rev(levels(data_long$variable)), 
                       labels = (rev(levels(data_long$variable)))) +
      labs(color = 'Covariate Value ')
    return(plot1)
  }  
  
}

levels(top20_cnt$variable) = rename_var(levels(top20_cnt$variable))
p <- my_shap_plot_summary(top20_cnt, x_bound  = 2, dilute = 20,
                       min_color_bound = '#3399ff',
                       max_color_bound  = "#ff5050") +
  labs(y = "SHAP Value for Fire Count XGBoost")
p


png(filename = file.path(dir.out, "SHAP_Fire_Count.pdf"), width = 4000, height = 2500, res=300)
print(p)
dev.off()


X_train_ba <- data.fit.train[,ba.covar.names.h1]
shap_values_ba <- shap.values(xgb_model = model_ba_h1, X_train = xgb.DMatrix(as.matrix(X_train_ba)) )

shap_long_ba <- shap.prep(shap_contrib = shap_values_ba$shap_score, X_train = as.matrix(X_train_ba))
levels(shap_long_ba$variable)


top20_ba <- shap_long_ba[shap_long_ba$variable %in% levels(shap_long_ba$variable)[1:10],] %>% droplevels()
levels(top20_ba$variable) = rename_var(levels(top20_ba$variable))
p <- my_shap_plot_summary(top20_ba, x_bound  = 5, dilute = 20,
                       min_color_bound = '#3399ff',
                       max_color_bound  = "#ff5050") + 
  labs(y = "SHAP Value for Burn Area XGBoost")
p
png(filename = file.path(dir.out, "SHAP_Burn_Area.pdf"), width = 4000, height = 2500, res=300)
print(p)
dev.off()


######################## latent effect  #####################

#---------------------extract all latent effects----------------
samples[[1]]

extract_var_name <- function(x){
  return (as.vector(strsplit(x, split=':')[[1]][1]))
}

effect_name <- sapply (rownames(samples[[1]]$latent), extract_var_name)

idx.score.ba = which(effect_name == 'score.ba')
idx.score.cnt = which(effect_name == 'score.cnt')
idx.score.z_cnt = which(effect_name == 'score.z_cnt')
idx.score.z_ba = which(effect_name == 'score.z_ba')

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


n.group <- 20

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

# score_z_grp <- rw.group(data.fit$score_z_h1, data.fit$score_z_h1, n.group)
# score.z <- c(nothing, score_z_grp, nothing)

score_z_grp_cnt <- rw.group(data.fit$score_cnt_h1, data.fit$score_cnt_h1, n.group)

score_z_grp_ba <- rw.group(data.fit$score_ba_h1, data.fit$score_ba_h1, n.group)


score_ba_grp <- rw.group(data.fit[data.fit$y>0,'score_ba_h1'],data.fit$score_ba_h1,n.group)



test2 = samples[[1]]$latent[idx.intercept.cnt,] + samples[[1]]$latent[idx.year.idx.cnt,][data.fit$year-min(data.fit$year)+1] + 
  samples[[1]]$latent[idx.grid.idx.cnt,][data.fit$grid.idx] + samples[[1]]$latent[idx.grid.idx.dist.cnt,][data.fit$grid.idx.district] + 
  samples[[1]]$latent[idx.score.cnt,][as.integer(as.factor(score_cnt_grp))]

# check the order of the group effect. The order is idx in group1, idx in group 2, ...
which(abs(test1 - test2)>0.0001)

total_time_idx = length(unique(data.fit$time.idx))
dist_effect_ba_array = array(dim = c(18,total_time_idx,1000))
council_effect_ba_array = array(dim = c(278,12,1000))

dist_effect_cnt_array = array(dim = c(18,total_time_idx,1000))
council_effect_cnt_array = array(dim = c(278,12,1000))

dist_effect_z_array = array(dim = c(18,total_time_idx,1000))
council_effect_z_array = array(dim = c(278,12, 1000))

score_ba_array = array(dim=c(length(idx.score.ba), 1000))
score_cnt_array = array(dim=c(length(idx.score.cnt), 1000))
score_z_cnt_array = array(dim=c(length(idx.score.z_cnt), 1000))
score_z_ba_array = array(dim=c(length(idx.score.z_ba), 1000))

year_effect_ba_array = array(dim=c(length(idx.year.idx.ba), 1000))
year_effect_cnt_array = array(dim=c(length(idx.year.idx.cnt), 1000))
year_effect_z_array = array(dim=c(length(idx.year.idx.z), 1000))

xi_array <- rep(NA,1000)
kappa_array <- rep(NA,1000)

for (i in 1:length(samples)){
  dist_effect_cnt_array[,,i] =  matrix(samples[[i]]$latent[idx.grid.idx.dist.cnt,][1:(length(idx.grid.idx.dist.cnt)/2)],nrow=18,ncol=total_time_idx)
  council_effect_cnt_array[,,i] = matrix(samples[[i]]$latent[idx.grid.idx.cnt,][1:(length(idx.grid.idx.cnt)/2)],nrow=278,ncol=12)

  dist_effect_ba_array[,,i] =  matrix(samples[[i]]$latent[idx.grid.idx.dist.ba,][1:(length(idx.grid.idx.dist.ba)/2)],nrow=18,ncol=total_time_idx)
  council_effect_ba_array[,,i] = matrix(samples[[i]]$latent[idx.grid.idx.ba,][1:(length(idx.grid.idx.ba)/2)],nrow=278,ncol=12)

  dist_effect_z_array[,,i] =  matrix(samples[[i]]$latent[idx.grid.idx.dist.z,][1:(length(idx.grid.idx.dist.z)/2)],nrow=18,ncol=total_time_idx)
  council_effect_z_array[,,i] = matrix(samples[[i]]$latent[idx.grid.idx.z,][1:(length(idx.grid.idx.z)/2)],nrow=278,ncol=12)
  
  score_ba_array[,i] = samples[[i]]$latent[idx.score.ba]
  score_cnt_array[,i] = samples[[i]]$latent[idx.score.cnt]
  score_z_cnt_array[,i] = samples[[i]]$latent[idx.score.z_cnt]
  score_z_ba_array[,i] = samples[[i]]$latent[idx.score.z_ba]
  
  year_effect_ba_array[,i] = samples[[i]]$latent[idx.year.idx.ba]
  year_effect_cnt_array[,i] = samples[[i]]$latent[idx.year.idx.cnt]
  year_effect_z_array[,i] = samples[[i]]$latent[idx.year.idx.z]
  
  xi_array[i] <- samples[[i]]$hyperpar[1]
  kappa_array[i] <- samples[[i]]$hyperpar[2]
  }

dist_effect_ba_mean <- apply(dist_effect_ba_array, c(1, 2), mean)
dist_effect_ba_q975 <-  apply(dist_effect_ba_array, c(1, 2), quantile, probs = 0.975)
dist_effect_ba_q025 <-  apply(dist_effect_ba_array, c(1, 2), quantile, probs = 0.025)

council_effect_ba_mean <- apply(council_effect_ba_array, c(1,2), mean)
council_effect_ba_q975 <-  apply(council_effect_ba_array, c(1,2), quantile, probs = 0.975)
council_effect_ba_q025 <-  apply(council_effect_ba_array, c(1,2), quantile, probs = 0.025)


dist_effect_cnt_mean <- apply(dist_effect_cnt_array, c(1, 2), mean)
dist_effect_cnt_q975 <-  apply(dist_effect_cnt_array, c(1, 2), quantile, probs = 0.975)
dist_effect_cnt_q025 <-  apply(dist_effect_cnt_array, c(1, 2), quantile, probs = 0.025)

council_effect_cnt_mean <- apply(council_effect_cnt_array, c(1,2), mean)
council_effect_cnt_q975 <-  apply(council_effect_cnt_array, c(1,2), quantile, probs = 0.975)
council_effect_cnt_q025 <-  apply(council_effect_cnt_array, c(1,2), quantile, probs = 0.025)


dist_effect_z_mean <- apply(dist_effect_z_array, c(1, 2), mean)
dist_effect_z_q975 <-  apply(dist_effect_z_array, c(1, 2), quantile, probs = 0.975)
dist_effect_z_q025 <-  apply(dist_effect_z_array, c(1, 2), quantile, probs = 0.025)

council_effect_z_mean <- apply(council_effect_z_array, c(1,2), mean)
council_effect_z_q975 <-  apply(council_effect_z_array, c(1,2), quantile, probs = 0.975)
council_effect_z_q025 <-  apply(council_effect_z_array, c(1,2), quantile, probs = 0.025)


score_ba_mean = apply(score_ba_array, 1, mean)
score_ba_q975 = apply(score_ba_array, 1, quantile, probs = 0.975)
score_ba_q025 = apply(score_ba_array, 1, quantile, probs = 0.025)

score_cnt_mean = apply(score_cnt_array, 1, mean)
score_cnt_q975 = apply(score_cnt_array, 1, quantile, probs = 0.975)
score_cnt_q025 = apply(score_cnt_array, 1, quantile, probs = 0.025)

score_z_cnt_mean = apply(score_z_cnt_array, 1, mean)
score_z_cnt_q975 = apply(score_z_cnt_array, 1, quantile, probs = 0.975)
score_z_cnt_q025 = apply(score_z_cnt_array, 1, quantile, probs = 0.025)

score_z_ba_mean = apply(score_z_ba_array, 1, mean)
score_z_ba_q975 = apply(score_z_ba_array, 1, quantile, probs = 0.975)
score_z_ba_q025 = apply(score_z_ba_array, 1, quantile, probs = 0.025)



year_effect_ba_mean = apply(year_effect_ba_array, 1, mean)
year_effect_ba_q975 = apply(year_effect_ba_array, 1, quantile, probs = 0.975)
year_effect_ba_q025 = apply(year_effect_ba_array, 1, quantile, probs = 0.025)

year_effect_cnt_mean = apply(year_effect_cnt_array, 1, mean)
year_effect_cnt_q975 = apply(year_effect_cnt_array, 1, quantile, probs = 0.975)
year_effect_cnt_q025 = apply(year_effect_cnt_array, 1, quantile, probs = 0.025)

year_effect_z_mean = apply(year_effect_z_array, 1, mean)
year_effect_z_q975 = apply(year_effect_z_array, 1, quantile, probs = 0.975)
year_effect_z_q025 = apply(year_effect_z_array, 1, quantile, probs = 0.025)


#------------------------Council level spatio-temporal effect-------------------

load(file=file.path(dir.out, 'map_council_district.RData'))


csc_mean_council <- scale_fill_gradientn(
  colours = c("#1a9850", '#91cf60', "#d9ef8b","#fee08b" ,'#fc8d59','#d73027'),
  name = "Mean",
  limits = c(min(
                 #  dist_effect_ba_mean,
                 # council_effect_ba_mean,
                 # dist_effect_cnt_mean,
                 # council_effect_cnt_mean,
                 # dist_effect_z_mean,
                 
                 council_effect_z_mean), 
             
             max(
               # dist_effect_ba_mean,
               #   council_effect_ba_mean,
               #   dist_effect_cnt_mean,
               #   council_effect_cnt_mean,
               #   dist_effect_z_mean,
                 council_effect_z_mean
             ))
)


map_plot = map
map_plot$value = council_effect_z_mean[,1]
map_plot$month = month.name[1]
for (i in 2:12){
  tmp = map 
  tmp$value = council_effect_z_mean[,i]
  tmp$month = month.name[i]
  map_plot = rbind(map_plot,tmp )
}


map_plot$month = factor(map_plot$month,levels = month.name )



p <-  ggplot(data = map_plot) +
  geom_sf(aes(fill = value)) +
  facet_wrap(~ month, ncol = 4) +
  csc_mean_council +
  labs(title = 'Council-Level Effect',
       fill = "Value") +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 30),
    axis.title.x = element_blank(),
    legend.position = "right",
    legend.key.size = unit(3, "lines"),
    legend.text = element_text(size = 30),
    legend.title = element_text(size = 30),        # Increased legend title size
    axis.line = element_line(colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.ticks = element_line(color = "black"),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.text.y = element_text(size = 30),
    axis.title.y = element_text(size = 30),
    plot.title = element_text(hjust = 0.5, size = 40),  # Center and increase title size
    strip.text = element_text(size = 30)                # Increase facet header size
  )
p

png(filename = file.path(dir.out, "council_effect_z_mean.pdf"), width =5000 , height = 6666, res=300)
print(p)
dev.off()

 
#------------------------District level spatio-temporal effect-------------------

map_plot = map.district
map_plot$mean = apply(dist_effect_z_mean,1, mean)
map_plot$sd = apply(dist_effect_z_mean,1, sd)


csc_mean_district <- scale_fill_gradientn(
  colours = c("#1a9850", '#91cf60', "#d9ef8b","#fee08b" ,'#fc8d59','#d73027'),
  name = "Mean",
  limits = c(min(
    map_plot$mean
  ), 
  
  max(
    map_plot$mean
  ))
)



p <-  ggplot(data = map_plot) +
  geom_sf(aes(fill = mean)) +
  csc_mean_district +
  labs(title = 'District-Level Effect',
       fill = "Value") +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 50),
    axis.title.x = element_blank(),
    legend.position = "right",
    legend.key.size = unit(3, "lines"),
    legend.text = element_text(size = 30),
    legend.title = element_blank(),        # Increased legend title size
    axis.line = element_line(colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.ticks = element_line(color = "black"),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.text.y = element_text(size = 50),
    axis.title.y = element_text(size = 50),
    plot.title = element_text(hjust = 0.5, size = 40),  # Center and increase title size
    strip.text = element_text(size = 50)                # Increase facet header size
  )
p

png(filename = file.path(dir.out, "district_effect_z_mean.pdf"), width =5000 , height = 6666, res=300)
print(p)
dev.off()


csc_sd_district <- scale_fill_gradientn(
  colours = c("#1a9850", '#91cf60', "#d9ef8b","#fee08b" ,'#fc8d59','#d73027'),
  name = "Mean",
  limits = c(min(
    map_plot$sd
  ), 
  
  max(
    map_plot$sd
  ))
)

p <-  ggplot(data = map_plot) +
  geom_sf(aes(fill = sd)) +
  csc_sd_district +
  labs(title = 'Council-Level Effect',
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
p

#------------------------Year effect-------------------

library(reshape2)

years <- 2014:2023
rownames(year_effect_cnt_array) <- years
rownames(year_effect_ba_array)  <- years
rownames(year_effect_z_array)   <- years

# Convert each matrix into a long-format data frame and add a column for the effect type
df_cnt <- melt(year_effect_cnt_array, varnames = c("Year", "Sample"), value.name = "Effect")
df_cnt$Type <- "Fire Count"

df_ba <- melt(year_effect_ba_array, varnames = c("Year", "Sample"), value.name = "Effect")
df_ba$Type <- "Burn Area"

df_z <- melt(year_effect_z_array, varnames = c("Year", "Sample"), value.name = "Effect")
df_z$Type <- "Fire Presence"

# Combine the three data frames
df <- rbind(df_cnt, df_ba, df_z)

# Filter out the year 2023
df <- subset(df, Year != "2023")

# Ensure Year is a factor with levels from 2011 to 2022
df$Year <- factor(df$Year, levels = as.character(2011:2022))


#------------------------XGBoost prediction effect-------------------


df_errorbar <- df%>% group_by(Year, Type)%>% summarize(mean = mean(Effect),
                                        lb = quantile(Effect,0.025),
                                        ub = quantile(Effect,0.975))
p <- ggplot(df_errorbar, aes(x = Year, y = mean)) +
  # Error bars still mapped to Type → dodge + colored by Type
  geom_errorbar(
    aes(ymin = lb, ymax = ub, color = Type),
    position = position_dodge(width = 0.6),
    width    = 0.5,
    linewidth     = 1.5
  ) +
  # Points share the same dodge‐grouping, but are drawn in black
  geom_point(
    aes(group = Type),       # ensures dodge offset matches the error bars
    color    = "black",      # force point color
    position = position_dodge(width = 0.6),
    size     = 2
  ) +
  scale_color_manual(values = c(
    "Fire Count"    = "#a8dadc",
    "Burn Area"     = "#ed6a5a",
    "Fire Presence" = "#F4D35E"
  )) +
  labs(y = "Year Effect") +
  theme(
    axis.text.x      = element_text(vjust = 0.5, size = 20),
    axis.title.x     = element_blank(),
    legend.position  = c(0.75, 0.9),
    legend.text      = element_text(size = 20),
    legend.title     = element_text(size = 20),
    axis.line        = element_line(colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.ticks       = element_line(color = "black"),
    panel.border     = element_blank(),
    panel.background = element_blank(),
    axis.text.y      = element_text(size = 20),
    axis.title.y     = element_text(size = 20),
    plot.title       = element_text(hjust = 0.5, size = 25)
  )

png(filename = file.path(dir.out, "Effect_Year.pdf"), width = 4000, height = 2000, res=300)
print(p)
dev.off()


index = 1:20


rownames(score_cnt_array) <- sort(unique(score_cnt_grp))
rownames(score_ba_array)  <-  sort(unique(score_ba_grp))
rownames(score_z_cnt_array)   <-  sort(unique(score_z_grp_cnt))
rownames(score_z_ba_array)   <-  sort(unique(score_z_grp_ba))

# Convert each matrix into a long-format data frame and add a column for the effect type
df_cnt <- melt(score_cnt_array, varnames = c("index", "Sample"), value.name = "Effect")
df_cnt$Type <- "cnt"

df_ba <- melt(score_ba_array, varnames = c("index", "Sample"), value.name = "Effect")
df_ba$Type <- "ba"

df_z_cnt <- melt(score_z_cnt_array, varnames = c("index", "Sample"), value.name = "Effect")
df_z_cnt$Type <- "z_cnt"

df_z_ba <- melt(score_z_ba_array, varnames = c("index", "Sample"), value.name = "Effect")
df_z_ba$Type <- "z_ba"

# Combine the three data frames


# Ensure Year is a factor with levels from 2011 to 2022
df_cnt$index <- factor(sprintf("%0.2f",df_cnt$index))
df_ba$index <- factor(sprintf("%0.2f",df_ba$index),levels=sprintf("%0.2f",unique(df_ba$index)))
# levels(df_ba$index) <-  sort(as.numeric(levels(df_ba$index)) )

df_z_cnt$index <- factor(sprintf("%0.2f",df_z_cnt$index))
df_z_ba$index <- factor(sprintf("%0.2f",df_z_ba$index))

df_z <- rbind(df_z_cnt, df_z_ba)



df_cnt_errorbar <- df_cnt%>% group_by(index)%>% summarize(mean = mean(Effect),
                                                           lb = quantile(Effect,0.025),
                                                           ub = quantile(Effect,0.975))
p <- ggplot(df_cnt_errorbar, aes(x = index, y = mean)) +
  geom_errorbar(
    aes(ymin = lb, ymax = ub),
    color = '#a8dadc',
    position = position_dodge(width = 0.6),
    width = 0.5,
    linewidth     = 1.5
  ) +
  # Draw a point at the mean
  geom_point(
    position = position_dodge(width = 0.6),
    size = 2
  ) +

  labs(y = "Effect of XGB_FC_Pred")+
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 30),
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
    axis.text.y = element_text(size = 30),
    axis.title.y = element_text(size = 30),
    plot.title = element_text(hjust = 0.5, size = 25)
  )
png(filename = file.path(dir.out, "Effect_Fire_Count.pdf"), width = 4000, height = 2000, res=300)
print(p)
dev.off()




df_ba_errorbar <- df_ba%>% group_by(index)%>% summarize(mean = mean(Effect),
                                                          lb = quantile(Effect,0.025),
                                                          ub = quantile(Effect,0.975))
p <- ggplot(df_ba_errorbar, aes(x = index, y = mean)) +
  geom_errorbar(
    aes(ymin = lb, ymax = ub),
    color = '#ed6a5a',
    position = position_dodge(width = 0.6),
    width = 0.5,
    linewidth     = 1.5
  ) +
  # Draw a point at the mean
  geom_point(
    position = position_dodge(width = 0.6),
    size = 2
  ) +
  
  labs(y = "Effect of XGB_BA_Pred")+
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 30),
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
    axis.text.y = element_text(size = 30),
    axis.title.y = element_text(size = 30),
    plot.title = element_text(hjust = 0.5, size = 25)
  )
png(filename = file.path(dir.out, "Effect_Burn_Area.pdf"), width = 4000, height = 2000, res=300)
print(p)
dev.off()



ggplot(df_z_cnt, aes(x = index, y = Effect)) +
  geom_boxplot(fill = "#a8dadc",position = position_dodge(width = 0.6),width = 0.5) +
  # Assign different colors to the data types (feel free to adjust these colors)
  # scale_fill_manual(values = c("cnt" = "#0D3B66",
  #                              "ba"  = "#ed6a5a",
  #                              "z"   = "#F4D35E")) +
  
  labs(y = "Effect of XGB_FC_Pred")+
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


ggplot(df_z_ba, aes(x = index, y = Effect)) +
  geom_boxplot(fill = "#ed6a5a",position = position_dodge(width = 0.6),width = 0.5) +
  # Assign different colors to the data types (feel free to adjust these colors)
  # scale_fill_manual(values = c("cnt" = "#0D3B66",
  #                              "ba"  = "#ed6a5a",
  #                              "z"   = "#F4D35E")) +
  
  labs(y = "Effect of XGB_BA_Pred")+
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

####################### Posterior of xi and kappa ##################### 
df <- data.frame(xi = xi_array)

# prior parameters
lambda <- 10       
xi_L   <- -0.5
xi_U   <-  0.5

# analytic normalizing constant
Z <- 2 * (1 - exp(-lambda * xi_U))

# truncated Laplace density
dtrunc_laplace <- function(x) {
  dens <- lambda * exp(-lambda * abs(x)) / Z
  dens[x < xi_L | x > xi_U] <- 0
  dens
}


p <- ggplot(df, aes(x = xi)) +
  # posterior density
  stat_density(aes(color = "Posterior"),
               geom   = "line",
               bw     = 0.01,   # or adjust = 1.5
               size   = 2,
               linetype = 'dashed'
  ) +
  # prior density
  stat_function(
    fun = dtrunc_laplace, 
    aes(color = "Prior"),
    size = 2, 
    linetype = 'dotdash',
    n = 512
  ) +
  scale_color_manual(
    name = "", 
    values = c("Posterior" = "#5459AC", "Prior" = "#CA7842")
  ) +
  labs(
    x = expression(xi), 
    y = expression("Density of " * xi)
  ) +
  xlim(-0.5, 0.5)+
  theme(
    axis.text.x = element_text(size = 30),
    axis.title.x = element_blank(),
    legend.position  = c(0.2, 0.8),
    legend.text = element_text(size = 30),
    legend.key.size = unit(3, "lines"),
    legend.title = element_text(size = 30),        # Increased legend title size
    axis.line = element_line(colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.ticks = element_line(color = "black"),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.text.y = element_text(size = 30),
    axis.title.y = element_text(size = 30),
    plot.title = element_text(hjust = 0.5, size = 25)
  )
png(filename = file.path(dir.out, "Posterior_density_xi.pdf"), width = 4000, height = 2000, res=300)
print(p)
dev.off()


dgamma_prior <- function(x){
  dens <- dgamma(x,shape=100,rate=100)
}
df <- data.frame(kappa = kappa_array)

p <- ggplot(df, aes(x = kappa)) +
  # posterior density
  stat_density(aes(color = "Posterior"),
               geom   = "line",
               adjust = 15,
               size   = 2,
               linetype = 'dashed'
  ) +
  # prior density
  stat_function(
    fun = dgamma_prior, 
    aes(color = "Prior"),
    size = 2, 
    linetype = 'dotdash',
    n = 512
  ) +
  scale_color_manual(
    name = "", 
    values = c("Posterior" = "#5459AC", "Prior" = "#CA7842")
  ) +
  labs(
    x = expression(xi), 
    y = expression("Density of " * kappa)
  ) +
  xlim(0, 4)+
  theme(
    axis.text.x = element_text(size = 30),
    axis.title.x = element_blank(),
    legend.position  = 'none',
    legend.text =  element_blank(),
    legend.key.size = unit(3, "lines"),
    legend.title = element_text(size = 30),        # Increased legend title size
    axis.line = element_line(colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.ticks = element_line(color = "black"),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.text.y = element_text(size = 30),
    axis.title.y = element_text(size = 30),
    plot.title = element_text(hjust = 0.5, size = 25)
  )

png(filename = file.path(dir.out, "Posterior_density_kappa.pdf"), width = 4000, height = 2000, res=300)
print(p)
dev.off()

######################## Estimated Likelihood ##############################
x <- seq(0,5,0.01)

kappa <- 2.739
xi <- 0.305

weibull_shape <- 1.394
lambda.ba <- 0.95
weibull_scale <- lambda.ba^(-1/weibull_shape)

median(rweibull(100000,shape = weibull_shape, scale = weibull_scale))

y_weibull <- dweibull(x,shape = weibull_shape, scale = weibull_scale)
plot(x,y_weibull, type='l')



mu.ba <- 0.97
pricision_gamma <- 1.868
a = pricision_gamma
b = mu.ba / a
median(rgamma(100000, shape = a, scale = b))

y_gamma <- dgamma(x, shape = a, scale = b)
plot(x,y_gamma, type='l')

dextGP <- function(y,xi, sigma, kappa, log=FALSE){
  y1 <- 1+xi*y/sigma
  out.of.sup <- y1<=0 | y<=0
  if (log){
    if (xi != 0){
      num <- log(kappa) + (kappa-1)*log(1-y1^(-1/xi))
      den <- log(sigma) + (1+1/xi)*log(y1)
    }else{
      num <- log(kappa) + (kappa-1)*log(1-exp(-y/sigma))
      den <- log(sigma) + y/sigma
    }
    res <- num-den
    res[which(out.of.sup)] <- -Inf
  }else{
    if(xi!=0){
      num <- kappa*(1-y1^(-1/xi))^(kappa-1)
      den <- sigma*y1^(1/xi +1)
    }else{
      num <- kappa*(1-exp(-y/sigma))^(kappa-1)
      den <- sigma*exp(y/sigma)
    }
    res <- num/den
    res[which(out.of.sup)] <- 0
  }
  return(res)
}

rextGP <- function(n, xi, sigma, kappa){
  U <- runif(n)
  U.G <- U^{1/kappa}
  if (xi!=0){
    return(sigma/xi*((1-U.G)^(-xi)-1))
  }else{
    return(-sigma*log(1-U.G))
  }
  
}
sigma <- 0.421
median(rextGP(100000,  xi, sigma, kappa))
y_egp <- dextGP(x,xi,sigma,kappa)
plot(x, y_egp, type='l')
lines(x, y_weibull, col='red')
lines(x, y_gamma, col='blue' )


df <- data.frame(
  x           = x,
  ExtendedGP  = y_egp,
  Weibull     = y_weibull,
  Gamma       = y_gamma
)

# 2. Pivot to “long” format
df_long <- pivot_longer(
  df,
  cols = c(ExtendedGP, Weibull, Gamma),
  names_to  = "Distribution",
  values_to = "Density"
)

# 3. Plot with ggplot2
p <- ggplot(df_long, aes(x = x, y = Density, 
                    color = Distribution, 
                    linetype = Distribution)) +
  geom_line(size = 2) +
  labs(
    y     = "Density",
  ) +
  theme(
    axis.text.x = element_text(size = 30),
    axis.title.x = element_blank(),
    legend.position  = c(0.7, 0.7),
    legend.text =  element_text(size = 30),
    legend.key.size = unit(3, "lines"),
    legend.title = element_text(size = 30),      # Increased legend title size
    axis.line = element_line(colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.ticks = element_line(color = "black"),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.text.y = element_text(size = 30),
    axis.title.y = element_text(size = 30),
    plot.title = element_text(hjust = 0.5, size = 30)
  )
png(filename = file.path(dir.out, "Density_likelihood.pdf"), width = 4000, height = 2000, res=300)
print(p)
dev.off()
 
