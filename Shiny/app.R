rm(list = ls())
library(shiny)
library(bslib)
library(tidyr)
library(raster)
library(ggplot2)
library(sf)
library(cowplot)

setwd('/home/pgrad2/2448355h/My_PhD_Project/Portugal_Wildfire/Shiny')
# dir.out <- '/home/pgrad2/2448355h/My_PhD_Project/01_Output/Urban_Fire_Model'

load("grid_cell_map.RData")
load('Final_Model_3.2_pred.sp_1000.RData')

pred.cnt <- pred.sp$pred.cnt
pred.ba <- pred.sp$pred.ba
# pred.z <- pred.sp$pred.z
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




gt.u <- function(x, u){
  return(sum(x>u))
}

gt.u.1 <- function(x, u){
  return(x>u)
}
gt.u.2 <- function(x1, u1, x2, u2){
  ind <- x1>u1 & x2>u2
  return(sum(ind))
}


sapply(pred.cnt,gt.u,2) / sapply(pred.cnt,gt.u,1) 
sapply(pred.cnt,gt.u,2) /200



B2_sf <- st_as_sf(B2.merge)
rm(B2.merge)

B2_sf <- gather(B2_sf, y.time.idx, y, paste0("y.", 1:108))
B2_sf$time.idx <- as.integer(substring(B2_sf$y.time.idx, 3, 5))



ui <- fluidPage(

  
  titlePanel("Forecasted Survival Probability"),
  
  sidebarLayout(
    sidebarPanel(
      # Input for year selection
      selectInput("year", "Select Year:", 
                  choices = 2012:2020, selected = 2019),
      
      # Input for month selection
      selectInput("month", "Select Month:", 
                  choices = month.name, selected = "August"),
      
      # Input for fire type
      radioButtons("type", "Select Type of Fire:",
                   choices = list("Count" = "count", "Burn Area" = "burn_area",'Joint'='joint')),
      
      # Conditional input for count threshold
      conditionalPanel(
        condition = "input.type == 'count' || input.type == 'joint'",
        numericInput("countThreshold", "Enter Count Threshold:", value = 1, min = 1, max=10)
      ),
      
      # Conditional input for burn area threshold
      conditionalPanel(
        condition = "input.type == 'burn_area' || input.type == 'joint'",
        numericInput("burn_areaThreshold", "Enter Burn Area Threshold:", value = 10, min = 1, max=1000)
      ),
      
      actionButton("generate", "Generate Map")
    ),
  
  
    mainPanel(
      # Output map
      card(plotOutput("map", height='90vh'))
      # card(textOutput( "memory_usage "))
    ),
    


  )
)

# Server logic ----
server <- function(input, output) {
  
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
    B2_sf[B2_sf$time.idx == idx, ]
    
    if (input$type=='count'){
      thres <- input$countThreshold
      B2_sf$prob <- sapply(pred.cnt[(1:192)+(idx-1)*192], gt.u, thres) /length(pred.cnt[[1]])
    }
    if (input$type=='burn_area'){
      thres <- log(max(input$burn_areaThreshold, 1))
      B2_sf$prob <- sapply(pred.ba[(1:192)+(idx-1)*192], gt.u, thres) /length(pred.ba[[1]])
    }
    if (input$type=='joint'){
      thres.ba <- log(max(input$burn_areaThreshold, 1))
      thres.cnt <- input$countThreshold
      excd.ba <- sapply(pred.ba[(1:192)+(idx-1)*192], gt.u.1, thres.ba) 
      excd.cnt <- sapply(pred.cnt[(1:192)+(idx-1)*192], gt.u.1, thres.cnt) 
      B2_sf$prob <- colMeans(excd.ba & excd.cnt) 
    }
    
    list('B2_sf'= B2_sf,
         'type'= input$type)

  })
  


  output$map <- renderPlot({
    req(Input.temp())

    csc.scale.fix <- scale_fill_gradient( low = "#F7F7F7", high = "#E4003A", limits = c(0,1))
    
    csc.scale.dym <- scale_fill_gradient( low = "#F7F7F7", high = "#E4003A", limits = c(min(Input.temp()$B2_sf$prob), max(Input.temp()$B2_sf$prob)))
    
    p1 <- ggplot(data=sf_districts) + 
      geom_sf(data=Input.temp()$B2_sf, aes(fill = prob),lwd = 0.01) +
      geom_sf(alpha=0) +
      theme(legend.position = "right",
            axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)
            # axis.text.y = element_blank(),
            # axis.ticks = element_blank()
      )+ csc.scale.fix +
      labs(fill  = "Prob",title='Prob on [0,1]' )
    
    p2 <- ggplot(data=sf_districts) + 
      geom_sf(data=Input.temp()$B2_sf, aes(fill = prob),lwd = 0.01) +
      geom_sf(alpha=0) +
      theme(legend.position = "right",
            axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)
            # axis.text.y = element_blank(),
            # axis.ticks = element_blank()
      )+ csc.scale.dym + 
      labs(fill  = "Prob",title='Prob rescaled' )
    
    p <- plot_grid(p1, p2, align = "hv", nrow = 1)
    
    if (Input.temp()$type=='count'){
      text <- 'P(Count > thres.cnt)'
    }else if (Input.temp()$type=='burn_area'){
      text <- 'P(Burn Area > thres.ba)'
    }else{
      text <- 'P(Count > thres.cnt, Burn Area > thres.ba)'
    }
    title <- ggdraw() + draw_label(text, fontface='bold')
    plot_grid(title, p,  align = "v", nrow = 2, rel_heights = c(0.1, 1))
  })

}

# Run app ----
shinyApp(ui, server)
