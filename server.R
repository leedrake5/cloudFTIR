library(shiny)
library(ggplot2)
library(pbapply)
library(reshape2)
library(dplyr)
library(data.table)
library(DT)
library(gridExtra)
library(rhandsontable)
library(Cairo)
library(broom)
library(shinyjs)
library(formattable)
library(markdown)
library(rmarkdown)
library(XML)
library(corrplot)
library(scales)
library(TTR)


shinyServer(function(input, output, session) {
    
    output$filegrab <- renderUI({
        
        if(input$filetype=="CSV") {
            fileInput('file1', 'Choose CSV', multiple=TRUE,
            accept=c(".csv"))
        } else if(input$filetype=="DPT") {
            fileInput('file1', 'Choose DPT', multiple=TRUE,
            accept=c(".dpt"))
        }
        
    })
    
    
    output$gainshiftui <- renderUI({
        
        if(input$advanced==TRUE){
           numericInput('gainshift', "Gain Shift (nm)", min=-1, max=1, value=0)
        } else {
            p()
        }
        
    })
    
    gainshiftHold <- reactive({
        
        if(input$advanced==TRUE){
            input$gainshift
        } else if(input$advanced==FALSE){
            0
        }
        
    })
    
    
    
    
    readDPT <- reactive({
        
            inFile <- input$file1
            if (is.null(inFile)) return(NULL)
            
            n <- length(inFile$datapath)
            names <- inFile$name
            
            myfiles.frame <- as.data.frame(do.call(rbind, lapply(seq(1, n, 1), function(x) readDPTData(filepath=inFile$datapath[x], filename=inFile$name[x]))))
        
        myfiles.frame$Wavelength <- myfiles.frame$Wavelength + gainshiftHold()
        
        myfiles.frame

    })
    
    
    readCSV <- reactive({
        
            inFile <- input$file1
            if (is.null(inFile)) return(NULL)
            
            n <- length(inFile$datapath)
            names <- inFile$name
            
            myfiles.frame <- as.data.frame(do.call(rbind, lapply(seq(1, n, 1), function(x) readCSVData(filepath=inFile$datapath[x], filename=inFile$name[x]))))
            
        myfiles.frame$Wavelength <- myfiles.frame$Wavelength + gainshiftHold()
        
        myfiles.frame
     
    })
    
    
    observeEvent(input$actionprocess, {
        

        myData <- reactive({
            
            data <- if(input$filetype=="DPT"){
                readDPT()
            } else if(input$filetype=="CSV"){
                readCSV()
            }
            
                data
        })
        
  

        dataHold <- reactive({
            data <- myData()
    
            data <- data[order(as.character(data$Spectrum)),]
    
            data$Spectrum <- gsub(".dpt", "", data$Spectrum)
            data$Spectrum <- gsub(".csv", "", data$Spectrum)
            data$Spectrum <- gsub(".CSV", "", data$Spectrum)

            data
    
        })
        
        
        dataCombine <- reactive({
            
            data <- dataHold()
            
            data <- data[order(data$Wavelength),]
            data$Intensity <- SMA(data$Intensity, 10)
            data$Spectrum <- rep("Combined", length(data$Intensity))
            
            data

        })
        
        dataBackgroundSubtract <- reactive({
            
            data <- if(input$combine==FALSE){
                dataHold()
            } else if(input$combine==TRUE){
                dataCombine()
            }
            
            data$Intensity <- Hodder.v(data$Intensity)
            
            data

            
        })
        
        
        dataManipulate <- reactive({
            
            if(input$combine==FALSE && input$backgroundsubtract==FALSE){
                dataHold()
            } else if(input$combine==TRUE && input$backgroundsubtract==FALSE){
                dataCombine()
            }  else if(input$combine==TRUE && input$backgroundsubtract==TRUE){
                dataBackgroundSubtract()
            } else if(input$combine==FALSE && input$backgroundsubtract==TRUE){
                dataBackgroundSubtract()
            }
            
            
        })


        dataCount <- reactive({
            inFile <- input$file1
    
            length(inFile$datapath)
    
        })


        observeEvent(input$actionplot, {
        
    
        ranges <- reactiveValues(x = NULL, y = NULL)
        
        
        
         plotInput <- reactive({
             
             data <- dataManipulate()

             n <- length(data$Wavelength)
             
            normal <- ggplot(data, aes(Wavelength, Intensity, colour=Spectrum)) +
             geom_line() +
             theme_light()+
             theme(legend.position="bottom") +
             scale_colour_discrete("Spectrum") +
             scale_x_reverse("Wavelength (nm)", breaks=seq(0, 4000, 250)) +
             scale_y_continuous("Intensity") +
             coord_cartesian(xlim = ranges$x, ylim = ranges$y)
             
             combine <- ggplot(data, aes(Wavelength, Intensity)) +
             geom_line() +
             theme_light()+
             theme(legend.position="bottom") +
             scale_x_reverse("Wavelength (nm)", breaks=seq(0, 4000, 250)) +
             scale_y_continuous("Intensity") +
             coord_cartesian(xlim = ranges$x, ylim = ranges$y)
             
             normal.invert <- ggplot(data, aes(Wavelength, Intensity, colour=Spectrum)) +
             geom_line() +
             theme_light()+
             theme(legend.position="bottom") +
             scale_colour_discrete("Spectrum") +
             scale_x_reverse("Wavelength (nm)", breaks=seq(0, 4000, 250)) +
             scale_y_reverse("Intensity") +
             coord_cartesian(xlim = ranges$x, ylim = ranges$y)
             
             combine.invert <- ggplot(data, aes(Wavelength, Intensity)) +
             geom_line() +
             theme_light()+
             theme(legend.position="bottom") +
             scale_x_reverse("Wavelength (nm)", breaks=seq(0, 4000, 250)) +
             scale_y_reverse("Intensity") +
             coord_cartesian(xlim = ranges$x, ylim = ranges$y)
             

             
             if(input$combine==FALSE && input$invert==FALSE){
                 normal
             } else if(input$combine==TRUE && input$invert==FALSE){
                 combine
             } else if(input$combine==FALSE && input$invert==TRUE){
                 normal.invert
             } else if(input$combine==TRUE && input$invert==TRUE){
                 combine.invert
             }
             

             

         })


        output$distPlot <- renderPlot({

            plotInput()

        })
        
        # When a double-click happens, check if there's a brush on the plot.
        # If so, zoom to the brush bounds; if not, reset the zoom.
        observeEvent(input$plot1_dblclick, {
            data <- dataHold()
            brush <- input$plot1_brush
            if (!is.null(brush)) {
                ranges$x <- c(brush$xmin, brush$xmax)
                ranges$y <- c(brush$ymin, brush$ymax)
                
            } else {
                ranges$x <- NULL
                ranges$y <- NULL
            }

        })
        
        output$hover_info_spectrum <- renderUI({
            
            point.table <- dataManipulate()

            hover <- input$plot_hover_spectrum
            point <- nearPoints(point.table,  coordinfo=hover,   threshold = 5, maxpoints = 1, addDist = TRUE)
            if (nrow(point) == 0) return(NULL)
            
            # calculate point position INSIDE the image as percent of total dimensions
            # from left (horizontal) and from top (vertical)
            left_pct <- (hover$x - hover$domain$left) / (hover$domain$right - hover$domain$left)
            top_pct <- (hover$domain$top - hover$y) / (hover$domain$top - hover$domain$bottom)
            
            # calculate distance from left and bottom side of the picture in pixels
            left_px <- hover$range$left + left_pct * (hover$range$right - hover$range$left)
            top_px <- hover$range$top + top_pct * (hover$range$bottom - hover$range$top)
            
            
            # create style property fot tooltip
            # background color is set so tooltip is a bit transparent
            # z-index is set so we are sure are tooltip will be on top
            style <- paste0("position:absolute; z-index:100; background-color: rgba(245, 245, 245, 0.85); ",
            "left:", left_px + 2, "px; top:", top_px + 2, "px;")
            
            # actual tooltip created as wellPanel
            wellPanel(
            style = style,
            p(HTML(paste0("Spectrum:", " ", point$Spectrum))),
            p(HTML(paste0("Wavelength:", " ", round(point$Wavelength, 0)))),
            p(HTML(paste0("Intensity:", " ", round(point$Intensity, 2))))
            )
        })
        
        output$downloadPlot <- downloadHandler(
        filename = function() { paste(input$dataset, '.png', sep='') },
        content = function(file) {
            ggsave(file,plotInput(), width=10, height=10)
        }
        )
        
         })
         
         

 })
 })




