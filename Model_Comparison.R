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
load(file.path(dir.out, "XGBoost_Score_Council.RData"))


get_plot_data <- function(file.name){

  if (file.name =='/home/pgrad2/2448355h/My_PhD_Project/01_Output/Portugal_Wildfire/single_district_weibull_pred_sp.RData'){
    load(file=file.name)
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
    
  }else{
    load(file=file.name)
  }

  
  
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
  df.crps.ba$type <- 'CRPS_Log_BA'

  
  crps.t.ba <- rep(NA, nrow(data.fit))
  for (i in 1:nrow(data.fit) ){
    y <- data.fit$log_ba[i]
    if (is.na(y)) y <- 0
    if (y>log(50)){
      crps.t.ba[i] <- crps_sample(
        y,
        pred.ba[[i]],
        method = "edf")
    }
  }
  
  df.crps.t.ba <- data.frame(crps=crps.t.ba,year=data.fit$year,district=data.fit$NAME_1)
  df.crps.t.ba$label <- 'Train'
  df.crps.t.ba[which(df.crps.t.ba$year>2022),'label'] <- 'Test'
  df.crps.t.ba$label <- factor(df.crps.t.ba$label, level=c('Train','Test'))
  df.crps.t.ba$type <- 'CRPS_Thres_Log_BA'


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
  df.crps.cnt$type <- 'CRPS_CNT'
  
  crps.t.cnt <- rep(NA, nrow(data.fit))
  for (i in 1:nrow(data.fit) ){
    y <- data.fit$y[i]
    if (is.na(y)) y <- 0
    if (y>3){
      crps.t.cnt[i] <- crps_sample(
        y,
        pred.cnt[[i]],
        method = "edf")
    }
  }
  
  df.crps.t.cnt <- data.frame(crps=crps.t.cnt,year=data.fit$year,district=data.fit$NAME_1)
  df.crps.t.cnt$label <- 'Train'
  df.crps.t.cnt[which(df.crps.t.cnt$year>2022),'label'] <- 'Test'
  df.crps.t.cnt$label <- factor(df.crps.t.cnt$label, level=c('Train','Test'))
  df.crps.t.cnt$type <- 'CRPS_Thres_CNT'

  crps.z <- rep(NA, nrow(data.fit))
  for (i in 1:nrow(data.fit) ){
    y <- data.fit$z[i]
    if (is.na(y)) y <- 0
    crps.z[i] <- crps_sample(
      y,
      pred.z[[i]],
      method = "edf")
  }
  round(sum(crps.z)/length(crps.z),4)
  
  df.crps.z <- data.frame(crps=crps.z,year=data.fit$year,district=data.fit$NAME_1)
  df.crps.z$label <- 'Train'
  df.crps.z[which(df.crps.z$year>2022),'label'] <- 'Test'
  df.crps.z$label <- factor(df.crps.z$label, level=c('Train','Test'))
  df.crps.z$type <- 'CRPS_Z'
  
  
  crps.t.z <- rep(NA, nrow(data.fit))
  for (i in 1:nrow(data.fit) ){
    z <- data.fit$z[i]
    y <- data.fit$log_ba[i]
    if (is.na(y)) y <- 0
    if (y>log(50)){
      crps.t.z[i] <- crps_sample(
        z,
        pred.z[[i]],
        method = "edf")
    }
  }
  
  df.crps.t.z <- data.frame(crps=crps.t.z,year=data.fit$year,district=data.fit$NAME_1)
  df.crps.t.z$label <- 'Train'
  df.crps.t.z[which(df.crps.t.z$year>2022),'label'] <- 'Test'
  df.crps.t.z$label <- factor(df.crps.t.z$label, level=c('Train','Test'))
  df.crps.t.z$type <- 'CRPS_Thres_Z'
  
  
  
  
  
  
  data.fit[data.fit$y==0,'log_ba'] <- 0
  data.fit$z <- as.integer(data.fit$y>0)
  data.fit$label <- 'Train'
  data.fit[which(data.fit$year>2022),'label'] <- 'Test'
  data.fit$label <- factor(data.fit$label, level=c('Train','Test'))
  
  
  data.fit[,'Estimated_CNT'] <- sapply(pred.cnt,mean)

  data.fit[,'Estimated_Z'] <- sapply(pred.z,mean)
  
  
  data.fit[,'Estimated_BA'] <- sapply(pred.ba,mean)

  library(pROC)
  auc_train <- auc(data.fit[data.fit$label=='Train',  'z'], data.fit[data.fit$label=='Train',  'Estimated_Z'], quiet = TRUE)
  auc_test <- auc(data.fit[data.fit$label=='Test',  'z'], data.fit[data.fit$label=='Test',  'Estimated_Z'], quiet = TRUE)
  
  
  return(list('CRPS'=rbind(df.crps.ba,df.crps.t.ba,df.crps.cnt,df.crps.t.cnt,
                           df.crps.z, df.crps.t.z),
              'AUC' = c(auc_train,auc_test),
              'MSE' = data.fit[,c('NAME_1','label','year','y','z','log_ba',
                                  'Estimated_CNT','Estimated_Z','Estimated_BA')]))
}

file.name.single <- file.path(dir.out,'single_district_weibull_pred_sp.RData')
signle_model_list <- get_plot_data(file.name.single)

file.name.joint <- file.path(dir.out,'Model_weibull_bym2_base_pred_sp_2000_subset.RData')
joint_model_list <- get_plot_data(file.name.joint)


df.crps.single <- signle_model_list$CRPS
df.crps.single$model <- "single"
df.crps.single%>% group_by(type,label) %>% summarise(crps=mean(crps,na.rm=T))

df.crps.joint <- joint_model_list$CRPS
df.crps.joint$model <- 'joint'
df.crps.joint%>% group_by(type,label) %>% summarise(crps=mean(crps,na.rm=T))

df.crps <- rbind(df.crps.single,df.crps.joint)

ggplot(df.crps[(df.crps$type=='CRPS_CNT')&(df.crps$label=='Train'),], 
       aes(x = district, y = crps,color=model)) +
  geom_boxplot()+
  ggtitle('CRPS of Fire Count in Training set')+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        axis.title.x = element_blank(),
        plot.title = element_text(hjust = 0.5)) 

ggplot(df.crps[(df.crps$type=='CRPS_CNT')&(df.crps$label=='Test'),], 
       aes(x = district, y = crps,color=model)) +
  geom_boxplot()+
  ggtitle('CRPS of Fire Count in Test set')+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        axis.title.x = element_blank(),
        plot.title = element_text(hjust = 0.5)) 


ggplot(df.crps[(df.crps$type=='CRPS_Thres_CNT')&(df.crps$label=='Train'),], 
       aes(x = district, y = crps,color=model)) +
  geom_boxplot()+
  ggtitle('CRPS of Count >3 in Training set')+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        axis.title.x = element_blank(),
        plot.title = element_text(hjust = 0.5)) 


ggplot(df.crps[(df.crps$type=='CRPS_Thres_CNT')&(df.crps$label=='Test'),], 
       aes(x = district, y = crps,color=model)) +
  geom_boxplot()+
  ggtitle('CRPS of Count >3 in Test set')+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        axis.title.x = element_blank(),
        plot.title = element_text(hjust = 0.5)) 


ggplot(df.crps[(df.crps$type=='CRPS_Log_BA')&(df.crps$label=='Train'),], 
       aes(x = district, y = crps,color=model)) +
  geom_boxplot()+
  ggtitle('CRPS of Log Burn Area in Training set')+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        axis.title.x = element_blank(),
        plot.title = element_text(hjust = 0.5)) 


ggplot(df.crps[(df.crps$type=='CRPS_Log_BA')&(df.crps$label=='Test'),], 
       aes(x = district, y = crps,color=model)) +
  geom_boxplot()+
  ggtitle('CRPS of Log Burn Area in Test set')+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        axis.title.x = element_blank(),
        plot.title = element_text(hjust = 0.5)) 


ggplot(df.crps[(df.crps$type=='CRPS_Thres_Log_BA')&(df.crps$label=='Train'),], 
       aes(x = district, y = crps,color=model)) +
  geom_boxplot()+
  ggtitle('CRPS of Log Burn Area > log(50) in Training set')+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        axis.title.x = element_blank(),
        plot.title = element_text(hjust = 0.5)) 

ggplot(df.crps[(df.crps$type=='CRPS_Thres_Log_BA')&(df.crps$label=='Test'),], 
       aes(x = district, y = crps,color=model)) +
  geom_boxplot()+
  ggtitle('CRPS of Log Burn Area > log(50) in Test set')+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        axis.title.x = element_blank(),
        plot.title = element_text(hjust = 0.5)) 






df.mse.single <- signle_model_list$MSE
df.mse.single$model <- "single"

df.mse.joint <- joint_model_list$MSE
df.mse.joint$model <- 'joint'
df.mse <- rbind(df.mse.single, df.mse.joint)
df.mse

df.mse$MSE_CNT <- (df.mse$Estimated_CNT - df.mse$y)^2
df.mse$MSE_Log_BA <- (df.mse$Estimated_BA - df.mse$log_ba)^2
# avg.log_ba.train <- mean(df.mse[df.mse$year<=2022&df.mse$y>0,'log_ba'])
# df.mse$log_ba1 <- df.mse$log_ba
# df.mse[df.mse$y==0,'log_ba1'] <- avg.log_ba.train
# df.mse$Exp_Z <- df.mse$Estimated_Z*df.mse$log_ba1
# df.mse$MSE_Z <- (df.mse$Exp_Z-df.mse$log_ba)^2

df.mse%>%group_by(model,label)%>%summarise(MSE_CNT=mean(MSE_CNT),MSE_Log_BA=mean(MSE_Log_BA))

df.mse[df.mse$y>3,]%>%group_by(model,label)%>%summarise(MSE_CNT=mean(MSE_CNT))
df.mse[df.mse$log_ba>log(50),]%>%group_by(model,label)%>%summarise(MSE_CNT=mean(MSE_CNT),MSE_Log_BA=mean(MSE_Log_BA))

ggplot(df.mse[(df.mse$label=='Train'),], 
       aes(x = NAME_1, y = MSE_CNT, color=model)) +
  geom_boxplot()+
  ggtitle('MSE of Count in Training set')+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        axis.title.x = element_blank(),
        plot.title = element_text(hjust = 0.5)) 

ggplot(df.mse[(df.mse$label=='Test'),], 
       aes(x = NAME_1, y = MSE_CNT, color=model)) +
  geom_boxplot()+
  ggtitle('MSE of Count in Test set')+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        axis.title.x = element_blank(),
        plot.title = element_text(hjust = 0.5)) 

ggplot(df.mse[(df.mse$label=='Train'&df.mse$y>3),], 
       aes(x = NAME_1, y = MSE_CNT, color=model)) +
  geom_boxplot()+
  ggtitle('MSE of Count > 3 in Training set')+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        axis.title.x = element_blank(),
        plot.title = element_text(hjust = 0.5)) 

ggplot(df.mse[(df.mse$label=='Test'&df.mse$y>3),], 
       aes(x = NAME_1, y = MSE_CNT, color=model)) +
  geom_boxplot()+
  ggtitle('MSE of Count > 3 in Test set')+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        axis.title.x = element_blank(),
        plot.title = element_text(hjust = 0.5)) 

ggplot(df.mse[(df.mse$label=='Train'),], 
       aes(x = NAME_1, y = MSE_Log_BA, color=model)) +
  geom_boxplot()+
  ggtitle('MSE of log Burn Area in Training set')+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        axis.title.x = element_blank(),
        plot.title = element_text(hjust = 0.5)) 


ggplot(df.mse[(df.mse$label=='Test'),], 
       aes(x = NAME_1, y = MSE_Log_BA, color=model)) +
  geom_boxplot()+
  ggtitle('MSE of log Burn Area in Test set')+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        axis.title.x = element_blank(),
        plot.title = element_text(hjust = 0.5)) 

ggplot(df.mse[(df.mse$label=='Train'&df.mse$log_ba>log(50)),], 
       aes(x = NAME_1, y = MSE_Log_BA, color=model)) +
  geom_boxplot()+
  ggtitle('MSE of log Burn Area > log(50) in Train set')+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        axis.title.x = element_blank(),
        plot.title = element_text(hjust = 0.5)) 

ggplot(df.mse[(df.mse$label=='Test'&df.mse$log_ba>log(50)),], 
       aes(x = NAME_1, y = MSE_Log_BA, color=model)) +
  geom_boxplot()+
  ggtitle('MSE of log Burn Area > log(50) in Test set')+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        axis.title.x = element_blank(),
        plot.title = element_text(hjust = 0.5)) 


ggplot(df.mse, 
       aes(x = NAME_1, y = y, color=label)) +
  geom_boxplot()+
  ggtitle('CNT')+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        axis.title.x = element_blank(),
        plot.title = element_text(hjust = 0.5)) 

ggplot(df.mse[df.mse$y>0,], 
       aes(x = NAME_1, y = y, color=label)) +
  geom_boxplot()+
  ggtitle('CNT | CNT> 0')+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        axis.title.x = element_blank(),
        plot.title = element_text(hjust = 0.5)) 

ggplot(df.mse, 
       aes(x = NAME_1, y = log_ba, color = label)) +
  geom_boxplot()+
  ggtitle('Log BA')+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        axis.title.x = element_blank(),
        plot.title = element_text(hjust = 0.5)) 

ggplot(df.mse[df.mse$y>0,], 
       aes(x = NAME_1, y = log_ba, color = label)) +
  geom_boxplot()+
  ggtitle('Log BA | Log BA >0')+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        axis.title.x = element_blank(),
        plot.title = element_text(hjust = 0.5)) 
