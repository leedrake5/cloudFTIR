library(shiny)
library(ggplot2)
library(pbapply)
library(reshape2)
library(dplyr)
library(plyr)
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

options(warn=-1)
assign("last.warning", NULL, envir = baseenv())

shinyServer(function(input, output, session) {
    
    options(warn=-1)

    
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
    
    
    observeEvent(!is.null(input$file1), {
        

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
        
        
        
        dataSplit <- reactive({
            
            data <- dataManipulate()
            
            index <- as.vector(unique(data$Spectrum))
            
            data.list <- lapply(index, function(x) subset(data, data$Spectrum==x))
            names(data.list) <- index
            
            data.list
            
        })
        
        findPeaks <- reactive({
            
            data.list <- dataSplit()
            
            index <- names(data.list)


            
            
            if(input$showpeaks==TRUE){
                data <- lapply(index, function(x) as.vector(peakpick(matrix(data.list[[x]][,3], ncol=1), neighlim=input$spikesensitivity, peak.min.sd=input$spikeheight*max(data.list[[x]][,3]))[,1]))
                names(data) <- index
                data
            } else if(input$showpeaks==FALSE){
                NULL
            }
            
        })
        
        
        peakTable <- reactive({
            
            data <- dataSplit()
            index <- names(data)

            for(i in 1:length(index)){
                data[[i]]$findpeaks <- findPeaks()[[i]]
                    
            }
            
            
            if(input$showpeaks==FALSE){
                data.frame(Spectrum=c("test"), Wavelength=c(-200), Intensity=c(0))
            } else if(input$showpeaks==TRUE){
                newdata <- lapply(names(data), function(x) subset(data[[x]], data[[x]]$Intensity > (input$spikeheight*max(data[[x]]$Intensity))))
                names(newdata) <- names(data)
                newdata <- lapply(names(newdata), function(x) as.data.frame(newdata[[x]][newdata[[x]]$findpeaks, ]))
                final.data <- as.data.frame(do.call(rbind, newdata))
                final.data
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
            
            
            data.list <- dataSplit()
            
            index <- names(data.list)
            
            if(input$showpeaks==FALSE){
                sliderInput('spikeheight', "Spike Height", min=0, max=1, value=.1)
            }else if(input$showpeaks==TRUE){
                sliderInput('spikeheight', "Spike Height", min=intensityDescriptive()[1], max=intensityDescriptive()[2], value=intensityDescriptive()[3])

            }
            
        })
        
        
        
        peakID <- reactive({
            
            pritable[pritable.range.index, ]$Min <- pritable[pritable.range.index, ]$Mid-input$wavethreshold
            pritable[pritable.range.index, ]$Max <- pritable[pritable.range.index, ]$Mid+input$wavethreshold
            
            n <- seq(from=1, to=length(peakTable()[,1]), by=1)
            
            table.list <- pblapply(n, function(x) in_range(spectrum=peakTable()$Spectrum[x],peak=peakTable()$Wavelength[x], pritable=pritable))
            
            data <- do.call(rbind, table.list)
            
            data[,c("Spectrum", "General", "Type", "Peak", "Max", "Mid", "Min")]

            
        })
        
        summaryID <- reactive({
            
            peak.table <- peakID()
            
            peak.table <- peak.table[,c("Spectrum", "Type", "Peak", "Min", "Max")]
            
           peak.summary <-  as.data.frame(peak.table %>%
            group_by(Type) %>%
            summarise_all(funs(toString)))
            
            peak.summary <-  as.data.frame(peak.table %>%
            group_by(Type, Spectrum) %>%
            summarise_all(funs(toString)))
            
            peak.summary$Min <- as.vector(sapply(peak.summary$Min, function(x) keep_singles(as.vector(unlist(strsplit(x, split=","))))))
            peak.summary$Max <- as.vector(sapply(peak.summary$Max, function(x) keep_singles(as.vector(unlist(strsplit(x, split=","))))))
            
            peak.summary

            
        })
        
        
        output$peaktableid <- renderDataTable({
            
            summaryID()
            
        })
        
        
        output$downloadPeakTableID <- downloadHandler(
        filename = function() { paste("FTIR Results", ".csv") },
        content = function(file
        ) {
            write.csv(summaryID(), file)
        }
        )
        
        
        summaryTable <- reactive({
            
            index <- as.vector(unique(peakID()$Spectrum))
            
            table.list <- lapply(index, function(x) subset(peakID(), peakID()$Spectrum==x))
            names(table.list) <- index


            n <- seq(from=1, to=length(table.list), by=1)

            
            table.tables <- lapply(n, function(x) table(table.list[[x]]$General))
            table.tables <- lapply(table.tables, as.data.frame)
            
            for(i in 1:length(index)){
                colnames(table.tables[[i]]) <- c("General", index[i])
            }
            
            #names(table.tables) <- index

            if(length(table.tables)>1){
                frame <- do.call(cbind, table.tables)[,c("General", index)]
                frame[ rowSums(frame[,-1])!=0, ]
            } else if(length(table.tables)==1) {
                table.tables[[1]][table.tables[[1]][,2]!=0, ]
            }
            
        })
        
        
        
        output$peaktablesummary <- renderDataTable({
            
            summaryTable()
            
        })
        
        
        output$downloadPeakTableSummary <- downloadHandler(
        filename = function() { paste("Summary Results", ".csv") },
        content = function(file
        ) {
            write.csv(summaryTable(), file)
        }
        )
        
        
        summaryPlot <- reactive({
            
            datamelt <- melt(summaryTable(), id="General")
            
            ggplot(datamelt, aes(General, value)) +
            geom_bar(stat = "identity", aes(fill = General), position = "dodge") +
            scale_y_continuous("Counts") +
            facet_grid(variable~.) +
            theme_light() +
            theme(strip.text.y = element_text(size=15)) +
            theme(axis.text.x = element_text(angle = 90, hjust = 1))  +
            theme(axis.text.x = element_text(size=15)) +
            theme(axis.text.y = element_text(size=15)) +
            theme(axis.title.x = element_text(size=15)) +
            theme(axis.title.y = element_text(size=15, angle=90)) +
            theme(plot.title=element_text(size=20)) +
            theme(legend.title=element_text(size=15)) +
            theme(legend.text=element_text(size=15))
            
        })
        
        output$summaryplot <- renderPlot({
            
            summaryPlot()
            
        })
        
        
        output$downloadsummaryplot <- downloadHandler(
        filename = function() { paste("Summary Plot") },
        content = function(file
        ) {
            ggsave(file,summaryPlot(), width=10, height=10)
        }
        )
        


    
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
             geom_point(data=peakTable(), aes(Wavelength, Intensity), shape=1, size=3) +
             theme(axis.text.x = element_text(size=15)) +
             theme(axis.text.y = element_text(size=15)) +
             theme(axis.title.x = element_text(size=15)) +
             theme(axis.title.y = element_text(size=15, angle=90)) +
             theme(plot.title=element_text(size=20)) +
             theme(legend.title=element_text(size=15)) +
             theme(legend.text=element_text(size=15))
             
             combine <- ggplot(data) +
             geom_line(aes(Wavelength, Intensity)) +
             theme_light()+
             theme(legend.position="bottom") +
             scale_x_reverse("Wavelength (nm)", breaks=seq(0, 4000, 250)) +
             scale_y_continuous("Intensity") +
             coord_cartesian(xlim = ranges$x, ylim = ranges$y) +
             geom_point(data=peakTable(), aes(Wavelength, Intensity), shape=1, size=3) +
             theme(axis.text.x = element_text(size=15)) +
             theme(axis.text.y = element_text(size=15)) +
             theme(axis.title.x = element_text(size=15)) +
             theme(axis.title.y = element_text(size=15, angle=90)) +
             theme(plot.title=element_text(size=20)) +
             theme(legend.title=element_text(size=15)) +
             theme(legend.text=element_text(size=15))
             
             normal.invert <- ggplot(data) +
             geom_line(aes(Wavelength, Intensity, colour=Spectrum)) +
             theme_light()+
             theme(legend.position="bottom") +
             scale_colour_discrete("Spectrum") +
             scale_x_reverse("Wavelength (nm)", breaks=seq(0, 4000, 250)) +
             scale_y_reverse("Intensity") +
             coord_cartesian(xlim = ranges$x, ylim = ranges$y) +
             geom_point(data=peakTable(), aes(Wavelength, Intensity), shape=1, size=3) +
             theme(axis.text.x = element_text(size=15)) +
             theme(axis.text.y = element_text(size=15)) +
             theme(axis.title.x = element_text(size=15)) +
             theme(axis.title.y = element_text(size=15, angle=90)) +
             theme(plot.title=element_text(size=20)) +
             theme(legend.title=element_text(size=15)) +
             theme(legend.text=element_text(size=15))
             
             combine.invert <- ggplot(data) +
             geom_line(aes(Wavelength, Intensity)) +
             theme_light()+
             theme(legend.position="bottom") +
             scale_x_reverse("Wavelength (nm)", breaks=seq(0, 4000, 250)) +
             scale_y_reverse("Intensity") +
             coord_cartesian(xlim = ranges$x, ylim = ranges$y) +
             geom_point(data=peakTable(), aes(Wavelength, Intensity), shape=1, size=3) +
             theme(axis.text.x = element_text(size=15)) +
             theme(axis.text.y = element_text(size=15)) +
             theme(axis.title.x = element_text(size=15)) +
             theme(axis.title.y = element_text(size=15, angle=90)) +
             theme(plot.title=element_text(size=20)) +
             theme(legend.title=element_text(size=15)) +
             theme(legend.text=element_text(size=15))
             

             
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
    
    
    
    waveInput <- reactive({
        
        blank.frame <- data.frame(
        Name=rep("", 25),
        WaveMin=rep("", 25),
        WaveMax=rep("", 25)
        )
        
        blank.frame$Name <- as.character(blank.frame$Name)
        blank.frame$WaveMin <- as.numeric(blank.frame$WaveMin)
        blank.frame$WaveMax <- as.numeric(blank.frame$WaveMax)
        
        blank.frame

        
    })
    
    
    
    wavevalues <- reactiveValues()
    
    
    
    

    
    

    
    waveInputCal <- reactive({
        
        #elements <- elementallinestouse()
        
        
        
        
        spectra.line.table <- dataManipulate()
        
        
        hold.frame <- waveInput()
        
        
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
    
    waveTableInput <- reactive({
        
        
        #if(input$usecalfile==FALSE){
        #    waveInput()
        #}else if(input$usecalfile==TRUE){
        #    waveInputCal()
        #}
        
        waveInput()
        
        
    })
    
    
    
    
    
    
    
    
    
    
    observe({
        if (!is.null(input$hotwave)) {
            DF <- hot_to_r(input$hotwave)
        } else {
            if (input$linecommitwave)
            DF <- waveTableInput()
            else
            DF <- wavevalues[["DF"]]
        }
        wavevalues[["DF"]] <- DF
    })
    
    eventReactive(input$linecommitwave,{
        
        wavevalues[["DF"]] <- waveTableInput()
        
    })
    
    

    
    
    ## Handsontable
    
    output$hotwave <- renderRHandsontable({
        
        DF <- wavevalues[["DF"]]
        
        DF <- DF[order(as.character(DF$Name)),]
        
        
        
        if (!is.null(DF))
        rhandsontable(DF) %>% hot_col(2:length(DF), renderer=htmlwidgets::JS("safeHtmlRenderer"))
        
        
    })
    
    
    observeEvent(input$resethotablewave, {
        
        wavevalues[["DF"]] <- NULL
        
        wavevalues[["DF"]] <- waveTableInput()
        
        
    })
    
    output$covarianceplotvalueswave <- renderPlot({
        
        data.table <- wavevalues[["DF"]]
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




