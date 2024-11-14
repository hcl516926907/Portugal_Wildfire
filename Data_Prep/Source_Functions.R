##Extract the names of municipalities per district
namesmunf=function(name){
  muns=data %>% count(district, municipality ,sort = TRUE)
  vec=c(unique(muns[muns$district==name,"municipality"]))
  return(vec)
}

## Identifies spatial points falling outside the district boundaries
outlimitdist=function(name){
  Dist=dist[dist$NAME_1==name,]
  spdf_error <- SpatialPointsDataFrame(coords = data[data$district==name, c("lon", "lat")], data = data[data$district==name,], proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs"))
  out<-sum(is.na(over(spdf_error, as(Dist, "SpatialPolygons"))))
  result=list() 
  if(out==0) {
    result<-print(paste0(name,": all points inside bounderies"))
  } else {
    points_out <- spdf_error[is.na(over(spdf_error, as(Dist, "SpatialPolygons"))), ]
    dfout=data.frame(dist_error=rep(name,length(points_out$lon)),lon=points_out$lon,lat=points_out$lat)
    spdf_true<-SpatialPointsDataFrame(coords = dfout[,c("lon", "lat")], data = dfout, proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs"))
    dfout$dist_true<-as.character(dist$NAME_1[over(spdf_true, as(dist, "SpatialPolygons"))])
    result<-dfout
  }
  return(result)
}

## Does it match? Give the number of districts that do not match
matchdist=function(name){
  # name="Aveiro"
  Dist=dist[dist$NAME_1==name,]
  spdf_error <- SpatialPointsDataFrame(coords = data[data$district==name, c("lon", "lat")], data = data[data$district==name,], proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs"))
  out<-sum(is.na(over(spdf_error, as(Dist, "SpatialPolygons"))))
  return(out)
}

## Identifies spatial points falling outside the municipality boundaries
outlimitmun=function(name){
  Mun=conc[conc$NAME_2==name,]
  spdf_error <- SpatialPointsDataFrame(coords = data[data$municipality==name, c("lon", "lat")], 
                                       data = data[data$municipality==name,], 
                                       proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs"))
  out<-sum(is.na(over(spdf_error, as(Mun, "SpatialPolygons"))))
  if(out==0) {
    result<-print(paste0(name,": all points inside bounderies"))
  } else {
    points_out <- spdf_error[is.na(over(spdf_error, as(Mun, "SpatialPolygons"))), ]
    dfout=data.frame(mun_error=rep(name,length(points_out$lon)),lon=points_out$lon,lat=points_out$lat)
    spdf_true<-SpatialPointsDataFrame(coords = dfout[,c("lon", "lat")], 
                                      data = dfout, 
                                      proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs"))
    dfout$mun_true<-as.character(conc$NAME_2[over(spdf_true, as(conc, "SpatialPolygons"))])
    result<-dfout
  }
  return(result)
}

## Does it match? Give the number of mun that do not match
matchmun=function(name){
  Mun=conc[conc$NAME_2==name,]
  spdf_error <- SpatialPointsDataFrame(coords = data[data$municipality==name, c("lon", "lat")], 
                                       data = data[data$municipality==name,], 
                                       proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs"))
  out<-sum(is.na(over(spdf_error, as(Mun, "SpatialPolygons"))))
  return(out)
}

## Identifies spatial points falling outside the parish boundaries
outlimitpar=function(name){
  Par=freg[freg$NAME_3==name,]
  spdf_error <- SpatialPointsDataFrame(coords = data[data$parish==name, c("lon", "lat")], 
                                       data = data[data$parish==name,], 
                                       proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs"))
  out<-sum(is.na(over(spdf_error, as(Par, "SpatialPolygons"))))
  if(out==0) {
    result<-print(paste0(name,": all points inside bounderies"))
  } else {
    points_out <- spdf_error[is.na(over(spdf_error, as(Par, "SpatialPolygons"))), ]
    dfout=data.frame(par_error=rep(name,length(points_out$lon)),lon=points_out$lon,lat=points_out$lat)
    spdf_true<-SpatialPointsDataFrame(coords = dfout[,c("lon", "lat")], 
                                      data = dfout, 
                                      proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs"))
    dfout$par_true<-as.character(freg$NAME_3[over(spdf_true, as(freg, "SpatialPolygons"))])
    result<-dfout
  }
  return(result)
}

## Does it match? Give the number of mun that do not match
matchpar=function(name){
  Par=freg[freg$NAME_3==name,]
  spdf_error <- SpatialPointsDataFrame(coords = data[data$parish==name, c("lon", "lat")], 
                                       data = data[data$parish==name,], 
                                       proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs"))
  out<-sum(is.na(over(spdf_error, as(Par, "SpatialPolygons"))))
  return(out)
}

listunmatch=function(vnames,fun){#vector of names given fpr example by namesd=sort(unique(data$district))
  # vnames=sort(unique(data$district))
  listc<-list()
  for (i in 1:length(vnames)) {
    listc[[i]]=fun(vnames[i])
  }
  out<-sum(unlist(listc))
  return(out)
}


## Plots a particular municipality and a set of spatial points
testplot=function(df,name){
  Mun=conc[conc$NAME_2==name,]
  plot(Mun,main=paste(name),cex.main=1)
  points(df$lon,df$lat,col="red",pch=16,cex=1)
}

## Finds the frequency of fires that occur in the same point within each district
nburn<-function(name,year){
  # year=2013
  # name="Coimbra"
  temp=data[as.numeric(format(data$open,format="%Y"))==year,]
  cpairs=table(temp$lon[temp$district==name],temp$lat[temp$district==name])
  tab=as.data.frame(table(as.numeric(cpairs)))[-1,]
  tab$total=as.numeric(as.character(tab$Var1))*tab$Freq
  tab$Var1=as.numeric(as.character(tab$Var1))
  tab$dist=rep(name,nrow(tab))
  tab$year=rep(year,nrow(tab))
  colnames(tab)=c("N. Fires", "N. Locations","n","District","Year")
  rownames(tab)=NULL
  return(tab)
}

# function to produce table with quantiles for quantitative variables
qtab=function(var){
  tab=as.data.frame(quantile(unlist(data[,paste(var)]),c(seq(0,0.9,0.1),0.95,0.99),na.rm=T))
  colnames(tab)=c(stringr::str_to_title(paste(var)))
  return(tab)
}

# function to produce table with counts and plot
ctab=function(var){
  tab1=data %>%
    mutate(x={{var}}) %>%
    count(x, sort = TRUE) %>% 
    mutate(Proportion=round(n/nrow(data),4)) %>% arrange(-Proportion) 
  
  p=ggplot(tab1, aes(x,Proportion)) + 
    geom_bar(stat = "identity")+
    labs(x=" ",y="Frequency")+theme_bw()
  return(list(Distribution.Table=tab1 %>%
                kbl(format = "simple",align = "l"),Distribution.Plot=p))
}

#identifies unmatching names
NameUnmatch=function(df,name1,name2){
  # df=conc
  # name1="municipality"
  # name2="NAME_2"
  comp=as.character(sort(unique(pull(data,name1))))==as.character(sort(unique(df@data[,name2])))
  v1=sort(unique(conc@data[,name2]))[which(comp==FALSE)[1]]
  v2=as.character(sort(unique(pull(data,name1)))[which(comp==FALSE)[1]])
  return(list(correct=v1,incorrect=v2))
}


# function to produce table with quantiles for quantitative variables
qtab=function(var){
  tab=as.data.frame(quantile(unlist(data[,paste(var)]),c(seq(0,0.9,0.1),0.95,0.99)))
  colnames(tab)=c(stringr::str_to_title(paste(var)))
  return(tab)
}

# function to produce table with counts and plot
ctab=function(var){
  tab1=data %>%
    mutate(x={{var}}) %>%
    count(x, sort = TRUE) %>% 
    mutate(Proportion=round(n/nrow(data),4)) %>% arrange(-Proportion) 
  
  p=ggplot(tab1, aes(x,Proportion)) + 
    geom_bar(stat = "identity")+
    labs(x=" ",y="Frequency")+theme_bw()
  return(list(Distribution.Table=tab1,Distribution.Plot=p))
}
