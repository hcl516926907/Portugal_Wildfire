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

setwd('/home/pgrad2/2448355h/My_PhD_Project/Portugal_Wildfire/Shiny')
dir.out <- '/home/pgrad2/2448355h/My_PhD_Project/01_Output/Urban_Fire_Model'

# load("grid_cell_map.RData")
# load('Final_Model_3.2_pred.sp_1000.RData')
# n.grid <- 192

load(file.path(dir.out,'data.fit2.score_0.125.RData'))
load(file.path(dir.out,'grid_cell_map_0.125.RData'))
load(file.path(dir.out,'Model_weibull_0.125_pred_sp_200.RData'))
n.grid <- 681

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

gt.u <- function(x, u){
  return(mean(x>u))
}

gt.u.1 <- function(x, u){
  return(x>u)
}
gt.u.2 <- function(x1, u1, x2, u2){
  ind <- x1>u1 & x2>u2
  return(sum(ind))
}


thres.exceed <- function(data, indices, thres) {
  return(gt.u(data[indices],thres))
}



B2_sf <- st_as_sf(B2.merge)


B2_sf <- gather(B2_sf, y.time.idx, y, paste0("y.", 1:108))
B2_sf$time.idx <- as.integer(substring(B2_sf$y.time.idx, 3, 5))

B2_sf$log_ba <- data.fit2$log_ba




ui <- fluidPage(
  
  
  titlePanel("Forecasted Survival Probability"),
  
  sidebarLayout(
    sidebarPanel(
      # Input for year selection
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
    # idx <- 12*(year-2012) + month
    idx <- month
    B2_sf <- B2_sf[B2_sf$time.idx == idx, ]
    
    pred.cnt.subset <- pred.cnt[(1:n.grid)+(idx-1)*n.grid]
    pred.ba.subset <- pred.ba[(1:n.grid)+(idx-1)*n.grid]
    

    
    if (input$type=='count'){
      thres <- input$countThreshold
      
      bootstrap_list <- lapply(pred.cnt.subset, function(sample) {
        boot(data = sample, statistic = thres.exceed, R = length(pred.cnt.subset), thres=thres)
      })
      
      prob = numeric(length(pred.cnt.subset))
      prob.sd = numeric(length(pred.cnt.subset))
      
      for (i in seq_along(bootstrap_list)) {
        prob[i] <- bootstrap_list[[i]]$t0  
        prob.sd[i] <- sd(bootstrap_list[[i]]$t)      
      }
      
      B2_sf$prob <- prob
      
      prob.sd.std <- (prob.sd - min(prob.sd))/max(prob.sd)
      B2_sf$neg.prob.sd.std <- 1- prob.sd.std
      
      labels <- c('Low','Median','High')
      prob.bin <- cut(prob, breaks = c(0, 0.05,0.2,1),labels=labels, include.lowest = TRUE)
      B2_sf$prob.bin <- prob.bin
      
      low_colors <- colorRampPalette(c("#a1d99b", "#006d2c"))      # Light to dark green
      median_colors <- colorRampPalette(c("#fed976", "#fd8d3c"))   # Light to dark orange
      high_colors <- colorRampPalette(c("#fc9272", "#de2d26"))  
      
      B2_sf <- B2_sf %>%
        group_by(prob.bin) %>%
        mutate(
          certainty = rescale(neg.prob.sd.std, to = c(1, 100)),
          color = case_when(
            prob.bin == "Low" ~ low_colors(100)[round(certainty)],
            prob.bin == "Median" ~ median_colors(100)[round(certainty)],
            prob.bin == "High" ~ high_colors(100)[round(certainty)]
          )
        ) %>%ungroup()
    }
    
    
    if (input$type=='burn_area'){
      thres <- log(max(input$burn_areaThreshold, 1))
      
      bootstrap_list <- lapply(pred.ba.subset, function(sample) {
        boot(data = sample, statistic = thres.exceed, R = length(pred.ba.subset), thres=thres)
      })
      
      prob = numeric(length(pred.ba.subset))
      prob.sd = numeric(length(pred.ba.subset))
      
      for (i in seq_along(bootstrap_list)) {
        prob[i] <- bootstrap_list[[i]]$t0  
        prob.sd[i] <- sd(bootstrap_list[[i]]$t)      
      }
      
      B2_sf$prob <- prob
      
      prob.sd.std <- (prob.sd - min(prob.sd))/max(prob.sd)
      B2_sf$neg.prob.sd.std <- 1- prob.sd.std
      
      labels <- c('Low','Median','High')
      prob.bin <- cut(prob, breaks = c(0, 0.1,0.2,1),labels=labels, include.lowest = TRUE)
      B2_sf$prob.bin <- prob.bin
      
      low_colors <- colorRampPalette(c("#a1d99b", "#006d2c"))      # Light to dark green
      median_colors <- colorRampPalette(c("#fed976", "#fd8d3c"))   # Light to dark orange
      high_colors <- colorRampPalette(c("#fc9272", "#de2d26"))  
      
      B2_sf <- B2_sf %>%
        group_by(prob.bin) %>%
        mutate(
          certainty = rescale(neg.prob.sd.std, to = c(1, 100)),
          color = case_when(
            prob.bin == "Low" ~ low_colors(100)[round(certainty)],
            prob.bin == "Median" ~ median_colors(100)[round(certainty)],
            prob.bin == "High" ~ high_colors(100)[round(certainty)]
          )
        ) %>% ungroup()
      
    }
    if (input$type=='joint'){
      thres.ba <- log(max(input$burn_areaThreshold, 1))
      thres.cnt <- input$countThreshold
      excd.ba <- sapply(pred.ba.subset, gt.u.1, thres.ba) 
      excd.cnt <- sapply(pred.ba.subset, gt.u.1, thres.cnt) 
      B2_sf$prob <- colMeans(excd.ba & excd.cnt) 
    }
    
    list('B2_sf'= B2_sf,
         'type'= input$type)

  })
  


  output$map <- renderPlot({
    req(Input.temp())

    csc.scale.fix <- scale_fill_gradient( low = "#F7F7F7", high = "#E4003A", limits = c(0,1))
    
    csc.scale.dym <- scale_fill_gradient( low = "#F7F7F7", high = "#E4003A", limits = c(min(Input.temp()$B2_sf$prob), max(Input.temp()$B2_sf$prob)))
    
   # Light to dark red
    
    
    p1 <- ggplot(data=sf_conc) + 
      geom_sf(data=Input.temp()$B2_sf, aes(fill = prob),lwd = 0.01) +
      geom_sf(alpha=0) +
      theme(legend.position = "right",
            axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)
            # axis.text.y = element_blank(),
            # axis.ticks = element_blank()
      )+ 
      csc.scale.fix +
      labs(fill  = "Prob",title='Prob on [0,1]' )
    
    # p2 <- ggplot(data=sf_districts) + 
    #   geom_sf(data=Input.temp()$B2_sf, aes(fill = prob),lwd = 0.01) +
    #   geom_sf(alpha=0) +
    #   theme(legend.position = "right",
    #         axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)
    #         # axis.text.y = element_blank(),
    #         # axis.ticks = element_blank()
    #   )+ csc.scale.dym + 
    #   labs(fill  = "Prob",title='Prob rescaled' )
    # legend_df <- expand.grid(
    #   prob = c("Low", "Median", "High"),
    #   inv.std = seq(0, 1, length.out = 100)
    # ) %>% group_by(prob) %>% mutate(
    #     color = case_when(
    #       prob == "Low" ~ low_colors(100),
    #       prob == "Median" ~ median_colors(100),
    #       prob == "High" ~ high_colors(100)
    #     )
    #   )%>% ungroup()
    # 
    p2 <- ggplot(sf_conc) +
      geom_sf(data=Input.temp()$B2_sf, aes(fill = color),lwd = 0.01) +
      scale_fill_identity() +
      geom_sf(alpha=0) +
      
      # geom_tile(data = legend_df, aes(x = prob, y = inv.std, fill = color)) +
      # facet_grid(~ prob, scales = "free", space = "free") +
      
      theme(legend.position = "right",
            axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)
            # axis.text.y = element_blank(),
            # axis.ticks = element_blank()
      )+ 
      labs(fill  = "Prob",title='Prob with Certainty' )
    
    p <- plot_grid(p1, p2, align = "h", nrow = 1)
    
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
