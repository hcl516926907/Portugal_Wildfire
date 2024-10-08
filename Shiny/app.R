library(shiny)
library(bslib)
library(tidyr)
library(raster)
library(ggplot2)
library(sf)
library(cowplot)
library(boot)
library(scales)
library(dplyr)
library(leaflet)

setwd('/home/pgrad2/2448355h/My_PhD_Project/Portugal_Wildfire/Shiny')
dir.out <- '/home/pgrad2/2448355h/My_PhD_Project/01_Output/Urban_Fire_Model'

# load("grid_cell_map.RData")
# load('Final_Model_3.2_pred.sp_1000.RData')
# n.grid <- 192

# load(file.path(dir.out,'data.fit2.score_0.125.RData'))
# load(file.path(dir.out,'grid_cell_map_0.125.RData'))
# load(file.path(dir.out,'Model_weibull_0.125_pred_sp_200.RData'))

load(file.path(dir.out, 'data.fit2.pred_0.0625.RData'))
load(file=file.path(dir.out,'Model_weibull_0.125_pred_0.0625_pred_sp_200.RData'))
load(file.path(dir.out,'grid_cell_map_0.0625.RData'))
# n.grid <- 681
n.grid <- 2554

pred.cnt <- pred.sp$pred.cnt
pred.ba <- pred.sp$pred.ba
pred.z <- pred.sp$pred.z
param.cnt <- pred.sp$param.cnt
param.z <- pred.sp$param.z
param.ba <- pred.sp$param.ba
hyper.ba <- pred.sp$hyper.ba

rm(pred.sp)


dist<-shapefile("distritos.shp")
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
ggplot(sf_conc)+ geom_sf()



grid.cell.coord <- st_as_sf(data.fit2.pred, coords = c("lon.grid", "lat.grid"), crs = 4326)

# merged_sf <- st_join(grid.cell.coord, sf_districts[,'NAME_1'], join = st_within)
merged_sf1 <- st_drop_geometry(st_join(grid.cell.coord, sf_conc[,'NAME_2'], join = st_nearest_feature))


B2_sf <- st_as_sf(B2.merge.pred)


B2_sf <- gather(B2_sf, y.time.idx, y, paste0("y.", 97:108))
B2_sf$time.idx <- as.integer(substring(B2_sf$y.time.idx, 3, 5))

data.fit2.pred[is.na(data.fit2.pred$log_ba),'log_ba'] <- 0
B2_sf$log_ba <- data.fit2.pred$log_ba

B2_sf <- merge(B2_sf,merged_sf1[,c('grid.idx','time.idx','NAME_2')],by=c('grid.idx','time.idx')) %>%
         arrange(time.idx,grid.idx)



ui <- fluidPage(
  titlePanel("Portugal Wildfire Analysis"),
  sidebarLayout(
    
    sidebarPanel(
      selectInput("year", "Select Year:", 
                  choices = 2020, selected = 2020),
      
      # Input for month selection
      selectInput("month", "Select Month:", 
                  choices = month.name, selected = "August"),
      
      # Input for fire type
      radioButtons("type", "Select Type of Fire:",
                   choices = list("Count" = "count", "Burn Area" = "burn_area",'Joint'='joint')),
      
      # Conditional input for count threshold
      conditionalPanel(
        condition = "input.type == 'count' || input.type == 'joint'",
        numericInput("countThreshold", "Enter Count Threshold:", value = 0, min = 1, max=100)
      ),
      
      # Conditional input for burn area threshold
      conditionalPanel(
        condition = "input.type == 'burn_area' || input.type == 'joint'",
        numericInput("burn_areaThreshold", "Enter Burn Area Threshold:", value = 10, min = 1, max=100000)
      ),
      
      actionButton("generate", "Generate Map"),
      h4("Click on a council to zoom in")
    ),
    mainPanel(
      leafletOutput("portugalMap", height = "60vh"),
      plotOutput("zoomedPlot", height = "40vh")
    )
  )
)

server <- function(input, output, session) {
  
  labels <- c('Low','Low_Median','Median','Median_High',  'High')
  breaks <- c(0, 0.05,0.1,0.25,0.5,1)

  low_colors <- colorRampPalette(c("#a1d99b", "#006d2c"))      # Light to dark green
  median_colors <- colorRampPalette(c("#fed976", "#fd8d3c"))   # Light to dark orange
  high_colors <- colorRampPalette(c("#fc9272", "#de2d26"))  
  
  low_median_colors <- colorRampPalette(c("#006d2c", "#fed976"))  # Dark green to light orange
  median_high_colors <- colorRampPalette(c("#fd8d3c", "#fc9272")) 
  
  
  Input.temp <- eventReactive(input$generate,{
    year <- as.integer(input$year)
    month <- switch(input$month,
                    "January"=1,
                    "February"=2,
                    "March"=3,
                    "April"=4,
                    "May"=5,
                    "June"=6,
                    "July"=7,
                    "August"=8,
                    "September"=9,
                    "October"=10,
                    "November"=11,
                    "December"=12)
    idx <- 12*(year-2012) + month
    B2_sf <- B2_sf[B2_sf$time.idx == idx, ]
    
    grid.month.idx <- (1+(month-1)*n.grid ): (month*n.grid)
    
    if (input$type=='count'){
      thres <- input$countThreshold
      
      prob.mark.mat <- sapply(1:200, function(x) ppois(thres,param.cnt[[x]][grid.month.idx], lower.tail=FALSE))
      prob.z.mat <- sapply(1:200, function(x) param.z[[x]][grid.month.idx])
      prob.mat <- prob.z.mat*prob.mark.mat
      
      
      prob <- rowMeans(prob.mat)
      prob.sd <- apply(prob.mat,1,sd)
      prob.sd.std <- (prob.sd - min(prob.sd))/max(prob.sd)
      
      B2_sf$prob <- prob
      B2_sf$neg.prob.sd.std <- 1- prob.sd.std
      
      prob.bin <- cut(prob, breaks = breaks,labels=labels, include.lowest = TRUE)
      B2_sf$prob.bin <- prob.bin

      B2_sf <- B2_sf %>%
        group_by(prob.bin) %>%
        mutate(
          certainty = rescale(neg.prob.sd.std, to = c(1, 100)),
          color = case_when(
            prob.bin == "Low" ~ low_colors(100)[round(certainty)],
            prob.bin == "Median" ~ median_colors(100)[round(certainty)],
            prob.bin == "High" ~ high_colors(100)[round(certainty)],
            prob.bin == 'Low_Median' ~ low_median_colors(100)[round(certainty)],
            prob.bin == 'Median_High' ~ median_high_colors(100)[round(certainty)]
          )
        ) %>%ungroup()
      
      
      df.prob.councils <- cbind(st_drop_geometry(B2_sf[,'NAME_2']),prob.mat)
      df.prob.councils.1 <- aggregate(.~NAME_2, df.prob.councils, mean)
      
      prob <- rowMeans(df.prob.councils.1[,c(-1)])
      prob.sd <- apply(df.prob.councils.1[,c(-1)],1,sd)
      prob.sd.std <- (prob.sd - min(prob.sd))/max(prob.sd)
      
      df.prob.councils.2 <-  data.frame('NAME_2'=df.prob.councils.1[,'NAME_2'])
      df.prob.councils.2$prob <- prob
      df.prob.councils.2$neg.prob.sd.std <- 1- prob.sd.std
      
      prob.bin <- cut(prob, breaks = breaks ,labels=labels, include.lowest = TRUE)
      
      df.prob.councils.2$prob.bin <- prob.bin
      
      df.prob.councils.2 <- df.prob.councils.2 %>%
        group_by(prob.bin) %>%
        mutate(
          certainty = rescale(neg.prob.sd.std, to = c(1, 100)),
          color = case_when(
            prob.bin == "Low" ~ low_colors(100)[round(certainty)],
            prob.bin == "Median" ~ median_colors(100)[round(certainty)],
            prob.bin == "High" ~ high_colors(100)[round(certainty)],
            prob.bin == 'Low_Median' ~ low_median_colors(100)[round(certainty)],
            prob.bin == 'Median_High' ~ median_high_colors(100)[round(certainty)]
          )
        ) %>%ungroup()
      
      councils <- merge(sf_conc, df.prob.councils.2,by='NAME_2')
    }
    
    
    if (input$type=='burn_area'){
      thres <- log(max(input$burn_areaThreshold, 1))
      
      prob.mark.mat <- sapply(1:200, function(x) ppois(thres,param.ba[[x]][grid.month.idx], lower.tail=FALSE))
      prob.z.mat <- sapply(1:200, function(x) param.z[[x]][grid.month.idx])
      prob.mat <- prob.z.mat*prob.mark.mat
      
      
      prob <- rowMeans(prob.mat)
      prob.sd <- apply(prob.mat,1,sd)
      prob.sd.std <- (prob.sd - min(prob.sd))/max(prob.sd)
      
      B2_sf$prob <- prob
      B2_sf$neg.prob.sd.std <- 1- prob.sd.std
      
      prob.bin <- cut(prob, breaks = breaks ,labels=labels, include.lowest = TRUE)
      B2_sf$prob.bin <- prob.bin
      
      
      B2_sf <- B2_sf %>%
        group_by(prob.bin) %>%
        mutate(
          certainty = rescale(neg.prob.sd.std, to = c(1, 100)),
          color = case_when(
            prob.bin == "Low" ~ low_colors(100)[round(certainty)],
            prob.bin == "Median" ~ median_colors(100)[round(certainty)],
            prob.bin == "High" ~ high_colors(100)[round(certainty)],
            prob.bin == 'Low_Median' ~ low_median_colors(100)[round(certainty)],
            prob.bin == 'Median_High' ~ median_high_colors(100)[round(certainty)]
          )
        ) %>%ungroup()
      
      
      df.prob.councils <- cbind(st_drop_geometry(B2_sf[,'NAME_2']),prob.mat)
      df.prob.councils.1 <- aggregate(.~NAME_2, df.prob.councils, mean)
      
      prob <- rowMeans(df.prob.councils.1[,c(-1)])
      prob.sd <- apply(df.prob.councils.1[,c(-1)],1,sd)
      prob.sd.std <- (prob.sd - min(prob.sd))/max(prob.sd)
      
      df.prob.councils.2 <-  data.frame('NAME_2'=df.prob.councils.1[,'NAME_2'])
      df.prob.councils.2$prob <- prob
      df.prob.councils.2$neg.prob.sd.std <- 1- prob.sd.std
      
      prob.bin <- cut(prob, breaks = breaks ,labels=labels, include.lowest = TRUE)
      
      df.prob.councils.2$prob.bin <- prob.bin
      
      df.prob.councils.2 <- df.prob.councils.2 %>%
        group_by(prob.bin) %>%
        mutate(
          certainty = rescale(neg.prob.sd.std, to = c(1, 100)),
          color = case_when(
            prob.bin == "Low" ~ low_colors(100)[round(certainty)],
            prob.bin == "Median" ~ median_colors(100)[round(certainty)],
            prob.bin == "High" ~ high_colors(100)[round(certainty)],
            prob.bin == 'Low_Median' ~ low_median_colors(100)[round(certainty)],
            prob.bin == 'Median_High' ~ median_high_colors(100)[round(certainty)]
          )
        ) %>%ungroup()
      
      councils <- merge(sf_conc, df.prob.councils.2,by='NAME_2')
      
    }
    if (input$type=='joint'){
      thres.ba <- log(max(input$burn_areaThreshold, 1))
      thres.cnt <- input$countThreshold
      
      prob.mark.mat1 <- sapply(1:200, function(x) ppois(thres.ba, param.cnt[[x]][grid.month.idx], lower.tail=FALSE))
      prob.mark.mat2 <- sapply(1:200, function(x) ppois(thres.cnt, param.ba[[x]][grid.month.idx], lower.tail=FALSE))
      prob.z.mat <- sapply(1:200, function(x) param.z[[x]][grid.month.idx])
      prob.mat <- prob.z.mat*prob.mark.mat1*prob.mark.mat2
      
      prob <- rowMeans(prob.mat)
      prob.sd <- apply(prob.mat,1,sd)
      prob.sd.std <- (prob.sd - min(prob.sd))/max(prob.sd)
      
      B2_sf$prob <- prob
      B2_sf$neg.prob.sd.std <- 1- prob.sd.std
      
      prob.bin <- cut(prob, breaks = breaks,labels=labels, include.lowest = TRUE)
      B2_sf$prob.bin <- prob.bin

      
      B2_sf <- B2_sf %>%
        group_by(prob.bin) %>%
        mutate(
          certainty = rescale(neg.prob.sd.std, to = c(1, 100)),
          color = case_when(
            prob.bin == "Low" ~ low_colors(100)[round(certainty)],
            prob.bin == "Median" ~ median_colors(100)[round(certainty)],
            prob.bin == "High" ~ high_colors(100)[round(certainty)],
            prob.bin == 'Low_Median' ~ low_median_colors(100)[round(certainty)],
            prob.bin == 'Median_High' ~ median_high_colors(100)[round(certainty)]
          )
        ) %>%ungroup()
      
      df.prob.councils <- cbind(st_drop_geometry(B2_sf[,'NAME_2']),prob.mat)
      df.prob.councils.1 <- aggregate(.~NAME_2, df.prob.councils, mean)
      
      prob <- rowMeans(df.prob.councils.1[,c(-1)])
      prob.sd <- apply(df.prob.councils.1[,c(-1)],1,sd)
      prob.sd.std <- (prob.sd - min(prob.sd))/max(prob.sd)
      
      df.prob.councils.2 <-  data.frame('NAME_2'=df.prob.councils.1[,'NAME_2'])
      df.prob.councils.2$prob <- prob
      df.prob.councils.2$neg.prob.sd.std <- 1- prob.sd.std
      
      prob.bin <- cut(prob, breaks = breaks ,labels=labels, include.lowest = TRUE)
      
      df.prob.councils.2$prob.bin <- prob.bin
      
      df.prob.councils.2 <- df.prob.councils.2 %>%
        group_by(prob.bin) %>%
        mutate(
          certainty = rescale(neg.prob.sd.std, to = c(1, 100)),
          color = case_when(
            prob.bin == "Low" ~ low_colors(100)[round(certainty)],
            prob.bin == "Median" ~ median_colors(100)[round(certainty)],
            prob.bin == "High" ~ high_colors(100)[round(certainty)],
            prob.bin == 'Low_Median' ~ low_median_colors(100)[round(certainty)],
            prob.bin == 'Median_High' ~ median_high_colors(100)[round(certainty)]
          )
        ) %>%ungroup()
      
      councils <- merge(sf_conc, df.prob.councils.2,by='NAME_2')
      
    }
    
    list('B2_sf'= B2_sf,
         'type'= input$type,
         'councils'=councils)
    
  })
  
  
  # Render the initial map with council boundaries
  output$portugalMap <- renderLeaflet({
    leaflet(Input.temp()$councils) %>%
      addTiles() %>%
      addPolygons(
        layerId = ~NAME_2,  # Replace 'NAME' with the column name containing council names
        fillColor = ~color,
        weight = 1,
        opacity = 1,
        color = "white",
        fillOpacity = 0.5
      )
  })
  
  # Observe click events on the map
  observeEvent(input$portugalMap_shape_click, {
    click <- input$portugalMap_shape_click
    if (is.null(click)) return()
    
    # Get the name of the clicked council
    selected_council <- click$id
    
    # Filter spatial data for the selected council
    council_data <- Input.temp()$councils[Input.temp()$councils$NAME_2 == selected_council, ]  # Replace 'NAME' with appropriate column
    
    # Generate zoomed-in plot for the selected council
    output$zoomedPlot <- renderPlot({
      # ggplot() +
      #   geom_sf(data = council_data, fill = "red", color = "black", size = 0.5) +
      #   ggtitle(paste("Zoomed-in Plot for", selected_council)) +
      #   theme_minimal()
      
      B2_sf_input <- Input.temp()$B2_sf
      ggplot(council_data)+
        geom_sf(data=B2_sf_input[B2_sf_input$NAME_2==selected_council,], aes(fill = color),lwd = 0.01) +
        scale_fill_identity()+
        geom_sf(alpha=0) +
        theme(legend.position = "right",
              axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)
              # axis.text.y = element_blank(),
              # axis.ticks = element_blank()
        )+ 
        labs(fill  = "Prob",title=paste('Probability in',selected_council ))
    })
  })
}



shinyApp(ui = ui, server = server)
