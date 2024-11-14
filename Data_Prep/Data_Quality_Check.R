library(ggplot2)
library(dplyr)
library(maps)
library(dbscan)
library(factoextra)
library(RgoogleMaps)
library(maps) 
library(mapdata) 
library(rgdal)
library(viridis)
library(RColorBrewer)
library(tidygeocoder)
library(dplyr)
library(raster)
library(tmap)
library(rgeos)
library(sf)
library(spatstat)
library(geosphere)
library(readxl)
library(deldir)
library(OpenStreetMap)
library(hms)
library(janitor)
library(stringr)
library(stringi)
library(here)
library(stargazer)
library(writexl)
library(formatR)
library(xlsx)
library(knitr)
library(lubridate)
library(validate)

dir.data <- "/home/pgrad2/2448355h/My_PhD_Project/00_Dataset/Porgual_Wildfire"
source("/home/pgrad2/2448355h/My_PhD_Project/Portugal_Wildfire/Data_Prep/Source_Functions.R")

raw1 <- read_xlsx(file.path(dir.data, "Registos_Incendios_SGIF_2011_2020.xlsx"))
raw2 <- read_xlsx(file.path(dir.data, "Registos_Incendios_SGIF_2021_2023.xlsx"))

raw <- rbind(raw1,raw2)


###########################################################
#unzip shaple files

# Districts

dist<-shapefile(file.path(dir.data,'shapefile', "distritos.shp"))
# Municipalities
conc<-shapefile(file.path(dir.data, "shapefile","concelhos.shp"))
# Parishes
freg<-shapefile(file.path(dir.data, "shapefile","freguesias.shp"))


# Districts
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
dist=dist[dist$NAME_1!="A莽ores",]
dist=dist[dist$NAME_1!="Madeira",]
dist$NAME_1<-as.factor(droplevels(dist$NAME_1))


# Municipalities
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


# Parishes
freg$NAME_1 <- as.factor(iconv(as.character(freg$NAME_1), "UTF-8"))
freg$NAME_2 <- as.factor(iconv(as.character(freg$NAME_2), "UTF-8"))
freg$NAME_3 <- as.factor(iconv(as.character(freg$NAME_3), "UTF-8"))
freg=freg[freg$NAME_1!="Azores",]
freg=freg[freg$NAME_1!="Madeira",]
freg$NAME_1<-as.factor(droplevels(freg$NAME_1))
freg$NAME_2<-as.factor(droplevels(freg$NAME_2))
freg$NAME_3<-as.factor(droplevels(freg$NAME_3))


n0=dim(raw);n0

data<-raw[!is.na(raw$AreaTotal_ha),]


# choose only closed record
table(data$Duracao_Horas)
data <- data[!is.na(data$Duracao_Horas), ]
n2 <- dim(data); n2


# Write consistently all names, using lowercase letters, 
# followed by capitalisation, both in the dataset and shapefiles
data$district<-factor(stringr::str_to_title(data$Distrito))
data$municipality<-factor(stringr::str_to_title(data$Concelho))
data$parish<-factor(stringr::str_to_title(data$Freguesia))

#shapefiles
dist$NAME_1<-stringr::str_to_title(dist$NAME_1)
conc$NAME_2<-stringr::str_to_title(conc$NAME_2)
freg$NAME_1<-stringr::str_to_title(freg$NAME_1)
freg$NAME_2<-stringr::str_to_title(freg$NAME_2)
freg$NAME_3<-stringr::str_to_title(freg$NAME_3)




dist.list<-as.character(sort(unique(data$district)))
mun.list<-as.character(sort(unique(data$municipality)))
par.list<-as.character(sort(unique(data$parish)))


select.columns <- c('AreaPov_ha','AreaMato_ha','AreaAgric_ha','AreaTotal_ha','DataHoraAlerta','DataHora_Extincao','Duracao_Horas',
                    'Longitude','Latitude','district','municipality','parish')

data <- data[,select.columns]
colnames(data) <- c('AreaPov_ha','AreaMato_ha','AreaAgric_ha','AreaTotal_ha','open','close','length',
                    'lon','lat','district','municipality','parish')

rules <- validator(
  # r1=is.numeric(data$id),
  # r2=is_unique(data$id),
  # r3=is.factor(data$importance),
  # r4=data$importance %in% c("low","medium","high","unknown"),
  # r5= data$importance[data$deaths==0 & (data$major+data$minor)<5 & data$length<30]=="low",
  # r6= data$importance[data$deaths>=1 | (data$major+data$minor)>=10 | data$length>=60]=="high",
  # r7= data$importance[data$deaths==0 & (data$major+data$minor)<10 & data$length>=30 & data$length<60]=="medium",
  r4=data$AreaPov_ha>=0,
  r5=data$AreaMato_ha>=0,
  r6=data$AreaAgric_ha>=0,
  r7=abs(data$AreaPov_ha+data$AreaMato_ha+data$AreaAgric_ha-data$AreaTotal_ha)<10^-6,
  r8=is.POSIXct(data$open),
  r9=in_range(data$open,min="2011-01-01 00:00:00",max="2023-12-31 23:59:00"),
  r10=is.POSIXct(data$close),
  r11=in_range(data$close,min="2011-01-01 00:00:00",max="2023-12-31 23:59:00"),
  r12=data$close > data$open,
  r13=is.numeric(data$length),
  r14=data$length>0,
  r15=abs(data$length - difftime(data$close,data$open,units = "hours"))<10^-6,
  # r16=is.factor(data$type),
  # r17=length(levels(data$type))<=length(structure.list),
  # r18=data$type %in% structure.list,
  # r19=is.factor(data$region),
  # r20=length(levels(data$region))<=length(region.list),
  # r21=data$region %in% region.list,
  # r22=is.factor(data$subregion),
  # r23=length(levels(data$subregion))<=length(subregion.list),
  # r24=data$subregion %in% subregion.list,
  r25=is.factor(data$district),
  r26=length(levels(data$district)) <= length(levels(as.factor(dist$NAME_1))),
  r27=data$district %in% dist$NAME_1,
  r28=listunmatch(dist.list,matchdist)==0,
  r29=is.factor(data$municipality),
  r30=length(levels(data$municipality)) <= length(levels(as.factor(conc$NAME_2))),  
  r31=data$municipality %in% conc$NAME_2, 
  r32=listunmatch(mun.list,matchmun)==0,
  r33=is.factor(data$parish),
  r34=length(levels(data$parish)) <= length(levels(as.factor(freg$NAME_3))),
  r35=data$parish %in% freg$NAME_3,
  r36=listunmatch(par.list,matchpar)==0,
  r37=is.numeric(data$lat),
  r38=in_range(data$lat, min=dist@bbox[2,1], max=dist@bbox[2,2]),
  r39=is.numeric(data$lon),
  r40=in_range(data$lon, min=dist@bbox[1,1], max=dist@bbox[1,2])
  # r41=is.factor(data$pfd),
  # r42=is.numeric(data$groundhr),
  # r43=is.numeric(data$groundtr),
  # r44=is.numeric(data$airhr),
  # r45=is.numeric(data$airtr),
  # r46=is.numeric(data$deaths),
  # r47=is.numeric(data$major),
  # r48=is.numeric(data$minor),
  # r49=is.numeric(data$assist),
  # r50=is.numeric(data$other),
  # r51=is.numeric(data$apc),
  # r52=is.numeric(data$otherv),
  # r53=data$groundhr >= 0,
  # r54=data$groundtr >= 0,
  # r55=data$airhr >= 0,
  # r56=airtr >= 0,
  # r57=data$deaths >= 0,
  # r58=data$major >= 0,
  # r59=data$minor >= 0,
  # r60=data$assist >= 0,
  # r61=data$other >= 0,
  # r62=data$apc >= 0,
  # r63=data$otherv >= 0,
  # r64=(data$deaths+data$major+data$minor+data$assist+data$other)==(data$apc+data$otherv),
  # r65=is.factor(data$reignit),
  # r66=length(levels(data$reignit))==2,
  # r67=data$reignit %in% c("no","yes")
)



summary(confront(data, rules))[,-c(6:8)]



data = data %>%  filter(data$length>0)
n3=dim(data);n3


# remove points that are outside the Portugal border
bound=dist[dist$NAME_0=="Portugal",]


spdf <- SpatialPointsDataFrame(coords = data[, c("lon", "lat")], data = data, proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs"))

pts<-data[is.na(over(spdf, as(bound, "SpatialPolygons"))),c("lon","lat")]


map=dist[dist$NAME_0=="Portugal",]
pts_sf <- st_as_sf(pts, coords = c("lon", "lat"), crs = 4326)
tmap_options(check.and.fix = TRUE)
tm_shape(map) + 
  tm_borders() +
  tm_text("NAME_1",size = 0.5, fontfamily="sans", ymod = 0.4)+
  tm_shape(pts_sf) + tm_dots(size=0.1,shape=4,col="red")  


data = data %>%  filter(!is.na(over(spdf, as(bound, "SpatialPolygons"))))
n4=dim(data);n4

NameUnmatch=function(df,name1,name2){
  # df=conc
  # name1="municipality"
  # name2="NAME_2"
  comp=as.character(sort(unique(pull(data,name1))))==as.character(sort(unique(df@data[,name2])))
  v1=sort(unique(conc@data[,name2]))[which(comp==FALSE)[1]]
  v2=as.character(sort(unique(pull(data,name1)))[which(comp==FALSE)[1]])
  return(list(correct=v1,incorrect=v2))
}

library(stringdist)

v1 <- c("apple", "ornge", "bannna") # Example vector with typos
v2 <- c("apple", "orange", "banana", "grape") # Correct labels
correct_labels <- function(mislabeled_vec, correct_vec) {
  sapply(mislabeled_vec, function(item) {
    distances <- stringdist::stringdist(item, correct_vec)
    corrected <- correct_vec[which.min(distances)]
    return(corrected)
  })
}

NameUnmatch1=function(df,name1,name2){
  v1 <- as.character(sort(unique(pull(data,name1))))
  v2 <- as.character(sort(unique(df@data[,name2])))
  v1_corrected <-  correct_labels(v1, v2)
  mismatch.ind <- which(names(v1_corrected)!=v1_corrected)
  if (length(mismatch.ind)>0){
    return(list(correct=as.character(v1_corrected[mismatch.ind]),incorrect=names(v1_corrected[mismatch.ind])))
  }else{
    print('All correct')
    return('Empty')
  }
}


ifelse(length(unique((data$district)))==length(unique((dist$NAME_1))),"Same number of district levels","Different number of district levels")
udnames=NameUnmatch1(dist,"district","NAME_1");udnames


ifelse(length(unique(data$municipality))==length(unique((conc$NAME_2))),"Same number of municipality levels","Different number of municipality levels")
umnames=NameUnmatch1(conc,"municipality","NAME_2");umnames

levels(data$municipality)[levels(data$municipality)==umnames$incorrect]<-umnames$correct


#############################################
#correct district label

namesd=sort(unique(data$district))
listoutd<-list()
for (i in 1:length(namesd)) {
  listoutd[[i]]=outlimitdist(namesd[i])
}
names(listoutd)<-namesd
dfs<-Filter(is.data.frame, listoutd)
df<-bind_rows(dfs)
rownames(df)<-NULL
df

# tif_sf <- st_as_sf(df, coords = c("lon", "lat"), crs = 4326)
# tif_sf$new_dist_error=paste0("  ",rownames(tif_sf),".",tif_sf$dist_error)
# 
# tmap_options(check.and.fix = TRUE)
# tm_shape(map) + 
#   tm_borders() +
#   tm_text("NAME_1",size = 0.5, fontfamily="sans", ymod = 0.4)+
#   tm_shape(tif_sf) + tm_dots(size=0.05,shape=4,col="red") +
#   tm_text("new_dist_error",size = 0.4,col="red",just="left") 


# for (i in 2:length(dfs)) {
#   print(i)
#   df=dfs[[i]]
#   spdf <- SpatialPointsDataFrame(coords = df[,c("lon", "lat")], data = df, proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs"))
#   newloc<-as.character(dist$NAME_1[over(spdf, as(dist, "SpatialPolygons"))])
#   data[data$lon %in% df$lon & data$lat %in% df$lat & data$district %in% df$dist_error,"district"]<-newloc
# }
spdf <- SpatialPointsDataFrame(coords = data[,c("lon", "lat")], data = data, proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs"))
data$district <- as.character(dist$NAME_1[over(spdf, as(dist, "SpatialPolygons"))])
clistoutd<-list()
for (i in 1:length(namesd)) {
  clistoutd[[i]]=outlimitdist(namesd[i])
}

n5=dim(data);n5

#############################################
#correct Municipalities label

listoutm<-list()
for (i in 1:length(unique(data$municipality))) {
  listoutm[[i]]=outlimitmun(unique(data$municipality)[i])
}
dfs_m<-listoutm[which(unlist(lapply(X = listoutm,FUN = is.data.frame)))]
df_m<-bind_rows(dfs_m)
rownames(df_m)<-NULL


df<-data.frame(data[,c("municipality","lon","lat")])
coordinates(df)=~lon+lat
proj4string(df) = proj4string(conc)
correctmun<-over(df, conc)

correctmun1 <- over(df, conc)
correctmun2 <- over(df, freg)

municipality1 <- factor(correctmun1$NAME_2)
municipality2 <- factor(correctmun2$NAME_2)
#correction
# if using municipality2, there would exsit inconsistance in Alcoutim and Tavira
data$municipality<- municipality1

clistoutm<-list()
for (i in 1:length(unique(data$municipality))) {#length=278 - certo
  clistoutm[[i]]=outlimitmun(unique(data$municipality)[i])
}
correctedm=length(clistoutm[which(unlist(lapply(X = clistoutm,FUN = is.data.frame)))])

n6=dim(data);n6

#############################################
#correct parishes label
clistoutp<-list()
for (i in 1:length(unique(data$parish))) {
  clistoutp[[i]]=outlimitpar(unique(data$parish)[i])
}

#extract dataframes with points outside boundaries
correctedp=length(clistoutp[which(unlist(lapply(X = clistoutp,FUN = is.data.frame)))])
correctedp

#points outside boundaries
df_p<-data.frame(data[,c("parish","lon","lat")])
coordinates(df_p)=~lon+lat
proj4string(df_p) = proj4string(freg)
correctpar<-over(df_p, freg)

#correction
data$parish<-factor(correctpar$NAME_3)

#checking
clistoutp_1<-list()
for (i in 1:length(unique(data$parish))) {
  clistoutp_1[[i]]=outlimitpar(unique(data$parish)[i])
}

# extract dataframes with points outside boundaries
correctedp_1=length(clistoutp_1[which(unlist(lapply(X = clistoutp_1,FUN = is.data.frame)))])

n7=dim(data);n7

confront(data, rules)
summary(confront(data, rules))[,-c(6:8)]

n8=dim(data);n8

# 
# size=as.data.frame(rbind(n1,n2,n3,n4,n5,n6,n7,n8))
# 
# for (i in 2:nrow(size)) {
#   size$dif.cases[1]=size$V1[1]-size$V1[1]
#   size$dif.cases[i]=size$V1[i-1]-size$V1[i]
# }
# 
# for (i in 2:nrow(size)) {
#   size$dif.vars[1]=size$V2[1]-size$V2[1]
#   size$dif.vars[i]=size$V2[i-1]-size$V2[i]
# }
# 
# rownames(size)=c("Start","Removal of variables","Negative length","Portuguese borders","District labels correction","Municipality labels correction","Parish labels correction","Final")
# colnames(size)=c("n","p","delta.n","delta.p")
# 
# size

str(data)
summary(data)

save(data, file=file.path(dir.data,"Wildfire.RData"))
