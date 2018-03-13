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
library(peakPick)
library(soil.spec)
library(parallel)

shinyServer(function(input, output, session) {
    
    output$filegrab <- renderUI({
        
        if(input$filetype=="CSV") {
            fileInput('file1', 'Choose CSV', multiple=TRUE,
            accept=c(".csv"))
        } else if(input$filetype=="DPT") {
            fileInput('file1', 'Choose DPT', multiple=TRUE,
            accept=c(".dpt"))
        } else if(input$filetype=="Opus") {
            fileInput('file1', 'Choose Opus File', multiple=TRUE,
            accept=NULL)
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
            
            myfiles.frame <- as.data.frame(do.call(rbind, pblapply(seq(1, n, 1), function(x) readDPTData(filepath=inFile$datapath[x], filename=inFile$name[x]))))
        
        myfiles.frame$Wavelength <- myfiles.frame$Wavelength + gainshiftHold()
        
        myfiles.frame

    })
    
    
    readCSV <- reactive({
        
            inFile <- input$file1
            if (is.null(inFile)) return(NULL)
            
            n <- length(inFile$datapath)
            names <- inFile$name
            
            myfiles.frame <- as.data.frame(do.call(rbind, pblapply(seq(1, n, 1), function(x) readCSVData(filepath=inFile$datapath[x], filename=inFile$name[x]))))
            
        myfiles.frame$Wavelength <- myfiles.frame$Wavelength + gainshiftHold()
        
        myfiles.frame
     
    })
    
    
    readOpus <- reactive({
        
        inFile <- input$file1
        if (is.null(inFile)) return(NULL)
        
        n <- length(inFile$datapath)
        names <- inFile$name
        
        myfiles.frame <- as.data.frame(do.call(rbind, pblapply(seq(1, n, 1), function(x) readOpusData(filepath=inFile$datapath[x], filename=inFile$name[x]))))
        
        myfiles.frame$Wavelength <- myfiles.frame$Wavelength + gainshiftHold()
        
        myfiles.frame
        
        
    })
    
    
    observeEvent(input$actionprocess, {
        

        myData <- reactive({
            
            data <- if(input$filetype=="DPT"){
                readDPT()
            } else if(input$filetype=="CSV"){
                readCSV()
            } else if(input$filetype=="Opus"){
                readOpus()
            }
            
                data
        })
        
  

        dataHold <- reactive({
            data <- myData()
    
            data <- data[order(as.character(data$Spectrum)),]
    
            data$Spectrum <- gsub(".dpt", "", data$Spectrum)
            data$Spectrum <- gsub(".csv", "", data$Spectrum)
            data$Spectrum <- gsub(".CSV", "", data$Spectrum)
            data$Spectrum <- gsub(".0", "", data$Spectrum)

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
        
        findPeaks <- reactive({
            
            data <- dataManipulate()
            
            
            
            
            
            if(input$showpeaks==TRUE){
                as.vector(peakpick(matrix(data[,3], ncol=1), neighlim=input$spikesensitivity, peak.min.sd=input$spikeheight)[,1])
            } else if(input$showpeaks==FALSE){
                NULL
            }
            
        })
        
        
        peakTable <- reactive({
            
            data <- dataManipulate()
            data$findpeaks <- findPeaks()
            
            if(input$showpeaks==FALSE){
                data.frame(Spectrum=c("test"), Wavelength=c(-200), Intensity=c(0))
            } else if(input$showpeaks==TRUE){
                newdata <- subset(data, data$Intensity > input$spikeheight)
                newdata[newdata$findpeaks, ]
            }
            
            
        })
        
        output$peaktable <- renderDataTable({
            
            peakTable()
            
        })
        
        
        output$downloadPeakTable <- downloadHandler(
        filename = function() { paste("FTIR Results", ".csv") },
        content = function(file
        ) {
            write.csv(peakTable(), file)
        }
        )
        
        intensityDescriptive <- reactive({
            
            min <- min(dataManipulate()[,3])
            max <- max(dataManipulate()[,3])
            mean <- mean(dataManipulate()[,3])
            
            
                c(min, max, mean)

        })
        
        output$uispikeheight <- renderUI({
            
            if(input$showpeaks==FALSE){
                sliderInput('spikeheight', "Spike Height", min=0, max=1, value=.1)
            }else if(input$showpeaks==TRUE){
                sliderInput('spikeheight', "Spike Height", min=intensityDescriptive()[1], max=intensityDescriptive()[2], value=intensityDescriptive()[3])

            }
            
        })
        


    
        ranges <- reactiveValues(x = NULL, y = NULL)
        
        
        
         plotInput <- reactive({
             
             data <- dataManipulate()

             n <- length(data$Wavelength)
             
            normal <- ggplot(data) +
             geom_line(aes(Wavelength, Intensity, colour=Spectrum)) +
             theme_light()+
             theme(legend.position="bottom") +
             scale_colour_discrete("Spectrum") +
             scale_x_reverse("Wavelength (nm)", limits=c(max(data[,2]), min(data[,2])), breaks=seq(0, 4000, 250)) +
             scale_y_continuous("Intensity") +
             coord_cartesian(xlim = ranges$x, ylim = ranges$y) +
             geom_point(data=peakTable(), aes(Wavelength, Intensity), shape=1, size=3)
             
             combine <- ggplot(data) +
             geom_line(aes(Wavelength, Intensity)) +
             theme_light()+
             theme(legend.position="bottom") +
             scale_x_reverse("Wavelength (nm)", breaks=seq(0, 4000, 250)) +
             scale_y_continuous("Intensity") +
             coord_cartesian(xlim = ranges$x, ylim = ranges$y) +
             geom_point(data=peakTable(), aes(Wavelength, Intensity), shape=1, size=3)
             
             normal.invert <- ggplot(data) +
             geom_line(aes(Wavelength, Intensity, colour=Spectrum)) +
             theme_light()+
             theme(legend.position="bottom") +
             scale_colour_discrete("Spectrum") +
             scale_x_reverse("Wavelength (nm)", breaks=seq(0, 4000, 250)) +
             scale_y_reverse("Intensity") +
             coord_cartesian(xlim = ranges$x, ylim = ranges$y) +
             geom_point(data=peakTable(), aes(Wavelength, Intensity), shape=1, size=3)
             
             combine.invert <- ggplot(data) +
             geom_line(aes(Wavelength, Intensity)) +
             theme_light()+
             theme(legend.position="bottom") +
             scale_x_reverse("Wavelength (nm)", breaks=seq(0, 4000, 250)) +
             scale_y_reverse("Intensity") +
             coord_cartesian(xlim = ranges$x, ylim = ranges$y) +
             geom_point(data=peakTable(), aes(Wavelength, Intensity), shape=1, size=3)
             

             
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
    
    
    
    
    
    hotableInputBlank <- reactive({
        
        #elements <- elementallinestouse()
        
        
        
        
        spectra.line.table <- if(input$filetype=="Spectra"){
            spectraData()
        } else if(input$filetype=="Elio"){
            spectraData()
        }  else if(input$filetype=="MCA"){
            spectraData()
        }  else if(input$filetype=="SPX"){
            spectraData()
        }  else if(input$filetype=="PDZ"){
            spectraData()
        } else if(input$filetype=="Net"){
            dataHold()
        }
        
        empty.line.table <- spectra.line.table[,elements] * 0.0000
        
        #empty.line.table$Spectrum <- spectra.line.table$Spectrum
        
        hold.frame <- data.frame(spectra.line.table$Spectrum, empty.line.table)
        colnames(hold.frame) <- c("Spectrum", elements)
        
        hold.frame <- as.data.frame(hold.frame)
        
        
        
        hold.frame
        
        
    })
    
    hotableInputCal <- reactive({
        
        #elements <- elementallinestouse()
        
        
        
        
        spectra.line.table <- dataManipulate()
        
        
        
        
        
        empty.line.table <- spectra.line.table[,elements] * 0.0000
        
        #empty.line.table$Spectrum <- spectra.line.table$Spectrum
        
        hold.frame <- data.frame(spectra.line.table$Spectrum, empty.line.table)
        colnames(hold.frame) <- c("Spectrum", elements)
        
        hold.frame <- as.data.frame(hold.frame)
        
        
        value.frame <- calFileContents()$Values
        
        #anna <- rbind(hold.frame, value.frame)
        
        #temp.table <- data.table(anna)[,list(result = sum(result)), elements]
        
        #as.data.frame(temp.table)
        
        #element.matches <- elements[elements %in% ls(value.frame)]
        
        #merge_Sum(.df1=hold.frame, .df2=value.frame, .id_Columns="Spectrum",  .match_Columns=element.matches)
        
        
        #data.frame(calFileContents()$Values, hold.frame[,! names(hold.frame) %in% names(calFileContents()$Values)])
        
        hold.frame.reduced <- hold.frame[2:length(hold.frame)]
        value.frame.reduced <- if(colnames(calFileContents()$Values)[1]=="Spectrum"){
            value.frame[2:length(value.frame)]
        }else if(colnames(calFileContents()$Values)[1]=="Include"){
            value.frame[3:length(value.frame)]
        }
        
        rownames(hold.frame.reduced) <- hold.frame$Spectrum
        rownames(value.frame.reduced) <- value.frame$Spectrum
        
        
        hotable.new = hold.frame.reduced %>% add_rownames %>%
        full_join(value.frame.reduced %>% add_rownames) %>%
        group_by(rowname) %>%
        summarise_all(funs(sum(., na.rm = FALSE)))
        
        colnames(hotable.new)[1] <- "Spectrum"
        
        hotable.new$Spectrum <- gsub(".pdz", "", hotable.new$Spectrum)
        hotable.new$Spectrum <- gsub(".csv", "", hotable.new$Spectrum)
        hotable.new$Spectrum <- gsub(".CSV", "", hotable.new$Spectrum)
        hotable.new$Spectrum <- gsub(".spt", "", hotable.new$Spectrum)
        hotable.new$Spectrum <- gsub(".mca", "", hotable.new$Spectrum)
        hotable.new$Spectrum <- gsub(".spx", "", hotable.new$Spectrum)
        
        hotable.new
        
        
    })
    
    hotableInput <- reactive({
        
        
        hotable.data <- if(input$usecalfile==FALSE){
            hotableInputBlank()
        }else if(input$usecalfile==TRUE){
            hotableInputCal()
        }
        
        
        
        hotable.new <- if(input$usecalfile==FALSE){
            data.frame(Include=rep(TRUE, length(hotable.data$Spectrum)), hotable.data)
        }else if(input$usecalfile==TRUE && colnames(calFileContents()$Values)[1]=="Spectrum"){
            data.frame(Include=rep(TRUE, length(hotable.data$Spectrum)), hotable.data)
        }else if(input$usecalfile==TRUE && colnames(calFileContents()$Values)[1]=="Include"){
            data.frame(Include=calFileContents()$Values[,1], hotable.data)
        }
        
        
        
    })
    
    
    
    
    
    values <- reactiveValues()
    
    
    
    
    
    #observe({
    #    if (!is.null(input$hot)) {
    #        DF <- hot_to_r(input$hot)
    #    } else {
    #        if (input$linecommit)
    #        DF <- hotableInput()
    #        else
    #        DF <- values[["DF"]]
    #    }
    #    values[["DF"]] <- DF
    #})
    
    eventReactive(input$linecommit,{
        
        values[["DF"]] <- hotableInput()
        
    })
    
    
    ## Handsontable
    
    output$hot <- renderRHandsontable({
        
        DF <- values[["DF"]]
        
        DF <- DF[order(as.character(DF$Spectrum)),]
        
        
        
        if (!is.null(DF))
        rhandsontable(DF) %>% hot_col(2:length(DF), renderer=htmlwidgets::JS("safeHtmlRenderer"))
        
        
    })
    
    
    observeEvent(input$resethotable, {
        
        values[["DF"]] <- NULL
        
        values[["DF"]] <- hotableInput()
        
        
    })
    
    output$covarianceplotvalues <- renderPlot({
        
        data.table <- values[["DF"]]
        correlations <- cor(data.table[,3:length(data.table)], use="pairwise.complete.obs")
        corrplot(correlations, method="circle")
        
    })
    
    # randomInterList <- reactive({
    #   if (is.null(input$intercept_vars))
    #   paste(,2)
    #   else
    #   input$intercept_vars
    #})
    
    
    #randomSlopeList <- reactive({
    #   if (is.null(input$intercept_vars))
    #   paste(,2)
    #   else
    #   input$slope_vars
    #})
    
    #output$nullintercept <- randomInterList()
    
    #output$nullslope <- randomSlopeList()

         
         

 })




