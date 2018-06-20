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
library(caret)
library(randomForest)
pdf(NULL)

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
        
        myfiles.frame$Wavenumber <- myfiles.frame$Wavenumber + gainshiftHold()
        
        myfiles.frame

    })
    
    
    readCSV <- reactive({
        
            inFile <- input$file1
            if (is.null(inFile)) return(NULL)
            
            n <- length(inFile$datapath)
            names <- inFile$name
            
            myfiles.frame <- as.data.frame(do.call(rbind, pblapply(seq(1, n, 1), function(x) readCSVData(filepath=inFile$datapath[x], filename=inFile$name[x]))))
            
        myfiles.frame$Wavenumber <- myfiles.frame$Wavenumber + gainshiftHold()
        
        myfiles.frame
     
    })
    
    
    readOpus <- reactive({
        
        inFile <- input$file1
        if (is.null(inFile)) return(NULL)
        
        n <- length(inFile$datapath)
        names <- inFile$name
        
        myfiles.frame <- as.data.frame(do.call(rbind, pblapply(seq(1, n, 1), function(x) readOpusData(filepath=inFile$datapath[x], filename=inFile$name[x]))))
        
        myfiles.frame$Wavenumber <- myfiles.frame$Wavenumber + gainshiftHold()
        
        myfiles.frame
        
        
    })
    
    
    
    calFileContents <- reactive({
        
        existingCalFile <- input$calfileinput
        
        if (is.null(existingCalFile)) return(NULL)
        
        
        Calibration <- readRDS(existingCalFile$datapath)
        
        
        
        Calibration[["Values"]]$Spectrum <- gsub(".spx", "", Calibration[["Values"]]$Spectrum)
        Calibration[["Values"]]$Spectrum <- gsub(".pdz", "", Calibration[["Values"]]$Spectrum)
        Calibration[["Values"]]$Spectrum <- gsub(".CSV", "", Calibration[["Values"]]$Spectrum)
        Calibration[["Values"]]$Spectrum <- gsub(".csv", "", Calibration[["Values"]]$Spectrum)
        Calibration[["Values"]]$Spectrum <- gsub(".spt", "", Calibration[["Values"]]$Spectrum)
        Calibration[["Values"]]$Spectrum <- gsub(".mca", "", Calibration[["Values"]]$Spectrum)
        
        Calibration[["Spectra"]]$Spectrum <- gsub(".spx", "", Calibration[["Spectra"]]$Spectrum)
        Calibration[["Spectra"]]$Spectrum <- gsub(".pdz", "", Calibration[["Spectra"]]$Spectrum)
        Calibration[["Spectra"]]$Spectrum <- gsub(".CSV", "", Calibration[["Spectra"]]$Spectrum)
        Calibration[["Spectra"]]$Spectrum <- gsub(".csv", "", Calibration[["Spectra"]]$Spectrum)
        Calibration[["Spectra"]]$Spectrum <- gsub(".spt", "", Calibration[["Spectra"]]$Spectrum)
        Calibration[["Spectra"]]$Spectrum <- gsub(".mca", "", Calibration[["Spectra"]]$Spectrum)
        
        Calibration
        
    })
    
    
    observeEvent(!is.null(input$file1) | !is.null(input$calfileinput), {
        

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
            data <- if(input$usecalfile==FALSE){
                myData()
            } else if(input$usecalfile==TRUE){
                calFileContents()[["Spectra"]]
            }
            
    
            data <- data[order(as.character(data$Spectrum)),]
    
            data$Spectrum <- gsub(".dpt", "", data$Spectrum)
            data$Spectrum <- gsub(".csv", "", data$Spectrum)
            data$Spectrum <- gsub(".CSV", "", data$Spectrum)
            data$Spectrum <- gsub(".0", "", data$Spectrum)

            data
    
        })
        
        
        dataCombine <- reactive({
            
            data <- dataHold()
            
            data <- data[order(data$Wavenumber),]
            data$Amplitude <- SMA(data$Amplitude, 10)
            data$Spectrum <- rep("Combined", length(data$Amplitude))
            
            data

        })
        
        dataBackgroundSubtract <- reactive({
            
            data <- if(input$combine==FALSE){
                dataHold()
            } else if(input$combine==TRUE){
                dataCombine()
            }
            
            data$Amplitude <- Hodder.v(data$Amplitude)
            
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
        
        
        output$downloadData <- downloadHandler(
        filename = function() { paste0(input$projectname, "_Raw", ".csv") },
        content = function(file
        ) {
            write.csv(dataManipulate(), file)
        }
        )


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
        
        
        idPeaks <- reactive({
            
            data.list <- dataSplit()
            
            index <- names(data.list)

            
                data <- lapply(index, function(x) as.vector(peakpick(matrix(data.list[[x]][,3], ncol=1), neighlim=input$spikesensitivity, peak.min.sd=input$spikeheight*max(data.list[[x]][,3]))[,1]))
                names(data) <- index
                
                data
            
        })
        
        idNull <- reactive({
            
            data.list <- dataSplit()
            
            index <- names(data.list)
            
            data <- lapply(index, function(x) rep(FALSE, length(data.list[[x]][,1])))
            names(data) <- index
            data
            
        })
        
        findPeaks <- reactive({
           
           if(input$showpeaks==TRUE){
               idPeaks()
            } else if(input$showpeaks==FALSE){
                idNull()
            }
           
            
        })
        
        output$downloadfp <- downloadHandler(
        filename <- function(){
            paste(input$projectname, "fp", sep=".")
        },
        
        content = function(file) {
            saveRDS(findPeaks(), file = file, compress="xz")
        }
        )
        
        
        peaks <- reactiveValues()
        peaks$findpeaks <- findPeaks()
        
        observe({
            if (input$showpeaks == FALSE){
            peaks$findpeaks <- idNull()
            }
        })
        
        observe({
            if (input$showpeaks == TRUE){
                peaks$findpeaks <- idPeaks()
            }
        })

        

        
        
        peakTable <- reactive({
            
            data <- dataSplit()
            index <- names(data)

            for(i in 1:length(index)){
                data[[i]]$findpeaks <- peaks$findpeaks[[i]]
                    
            }
            
            
                newdata <- lapply(names(data), function(x) subset(data[[x]], data[[x]]$Amplitude > (input$spikeheight*max(data[[x]]$Amplitude))))
                names(newdata) <- names(data)
                newdata <- lapply(names(newdata), function(x) as.data.frame(newdata[[x]][newdata[[x]]$findpeaks, ]))
                final.data <- as.data.frame(do.call(rbind, newdata))
                final.data[,c("Spectrum", "Wavenumber", "Amplitude")]

            
            
            
        })
        
        output$peaktable <- renderDataTable({
            
            peakTable()
            
        })
        
        
        output$downloadPeakTable <- downloadHandler(
        filename = function() { paste0(input$projectname, "_Peaks", ".csv") },
        content = function(file
        ) {
            write.csv(peakTable(), file)
        }
        )
        
        
        AmplitudeDescriptive <- reactive({
            
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
                sliderInput('spikeheight', "Spike Height", min=AmplitudeDescriptive()[1], max=AmplitudeDescriptive()[2], value=AmplitudeDescriptive()[3])

            }
            
        })
        
        
        peakID <- reactive({
            
            pritable[pritable.range.index, ]$Min <- pritable[pritable.range.index, ]$Mid-input$wavethreshold
            pritable[pritable.range.index, ]$Max <- pritable[pritable.range.index, ]$Mid+input$wavethreshold
            
            n <- seq(from=1, to=length(peakTable()[,1]), by=1)
            
            table.list <- pblapply(n, function(x) in_range(spectrum=peakTable()$Spectrum[x],peak=peakTable()$Wavenumber[x], pritable=pritable))
            
            data <- do.call(rbind, table.list)
            
            data[,c("Spectrum", "General", "Group", "Type", "Peak", "Max", "Mid", "Min", "Peak.Range")]
            
            
        })
        
        
        
        peakIDDisplay <- reactive({
            

            
            new.data <- peakID()
            new.data <- new.data[,c("Spectrum", "Group", "Type", "Peak", "Max", "Mid", "Min", "Peak.Range")]

            new.data <- new.data[!duplicated(new.data), ]
            
            new.data[order(new.data$Group, new.data$Type, -new.data$Peak),]

            
        })
        
        summaryID <- reactive({
            

            peak.table <- peakIDDisplay()
            
            peak.table <- peak.table[,c("Group", "Peak.Range", "Spectrum", "Type", "Peak")]
            
            
            peak.summary <-  as.data.frame(peak.table %>%
            group_by(Group, Peak.Range, Type, Spectrum) %>%
            summarise_all(funs(toString)))
            
            
            
            #peak.summary$Min <- as.vector(sapply(peak.summary$Min, function(x) keep_singles(strsplit(x, split=","))))
            #peak.summary$Max <- as.vector(sapply(peak.summary$Max, function(x) keep_singles(strsplit(x, split=","))))
            
            peak.result <- dcast(peak.summary, Group+Peak.Range+Type ~ Spectrum)
            peak.result <- peak.result[order(peak.result$Group, peak.result$Type),]
            peak.result[is.na(peak.result)]   <- " "
            peak.result
            
        })
        
        
        output$peaktableid <- renderDataTable({
            
            summaryID()
            
        })
        
        
        output$downloadPeakTableID <- downloadHandler(
        filename = function() { paste0(input$projectname, "_PeakID", ".csv") },
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
        filename = function() { paste0(input$projectname, "_Summary", ".csv") },
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
        filename = function() { paste0(input$projectname, "_Summary", ".jpg") },
        content = function(file
        ) {
            ggsave(file,summaryPlot(), width=15, height=10, device="jpeg")
        }
        )
        


    
        ranges <- reactiveValues(x = NULL, y = NULL)
        
        
        
         plotInput <- reactive({
             
             data <- dataManipulate()

             n <- length(data$Wavenumber)
             
            normal <- ggplot(data) +
             geom_line(aes(Wavenumber, Amplitude, colour=Spectrum)) +
             theme_light()+
             theme(legend.position="bottom") +
             scale_colour_discrete("Spectrum") +
             scale_x_reverse(expression(paste("Wavenumber (cm"^"-1"*")")), limits=c(max(data[,2]), min(data[,2])), breaks=seq(0, 4000, 250)) +
             scale_y_continuous("Amplitude") +
             coord_cartesian(xlim = ranges$x, ylim = ranges$y) +
             geom_point(data=peakTable(), aes(Wavenumber, Amplitude), shape=1, size=3) +
             theme(axis.text.x = element_text(size=15)) +
             theme(axis.text.y = element_text(size=15)) +
             theme(axis.title.x = element_text(size=15)) +
             theme(axis.title.y = element_text(size=15, angle=90)) +
             theme(plot.title=element_text(size=20)) +
             theme(legend.title=element_text(size=15)) +
             theme(legend.text=element_text(size=15))
             
             combine <- ggplot(data) +
             geom_line(aes(Wavenumber, Amplitude)) +
             theme_light()+
             theme(legend.position="bottom") +
             scale_x_reverse(expression(paste("Wavenumber (cm"^"-1"*")")), breaks=seq(0, 4000, 250)) +
             scale_y_continuous("Amplitude") +
             coord_cartesian(xlim = ranges$x, ylim = ranges$y) +
             geom_point(data=peakTable(), aes(Wavenumber, Amplitude), shape=1, size=3) +
             theme(axis.text.x = element_text(size=15)) +
             theme(axis.text.y = element_text(size=15)) +
             theme(axis.title.x = element_text(size=15)) +
             theme(axis.title.y = element_text(size=15, angle=90)) +
             theme(plot.title=element_text(size=20)) +
             theme(legend.title=element_text(size=15)) +
             theme(legend.text=element_text(size=15))
             
             normal.invert <- ggplot(data) +
             geom_line(aes(Wavenumber, Amplitude, colour=Spectrum)) +
             theme_light()+
             theme(legend.position="bottom") +
             scale_colour_discrete("Spectrum") +
             scale_x_reverse(expression(paste("Wavenumber (cm"^"-1"*")")), breaks=seq(0, 4000, 250)) +
             scale_y_reverse("Amplitude") +
             coord_cartesian(xlim = ranges$x, ylim = ranges$y) +
             geom_point(data=peakTable(), aes(Wavenumber, Amplitude), shape=1, size=3) +
             theme(axis.text.x = element_text(size=15)) +
             theme(axis.text.y = element_text(size=15)) +
             theme(axis.title.x = element_text(size=15)) +
             theme(axis.title.y = element_text(size=15, angle=90)) +
             theme(plot.title=element_text(size=20)) +
             theme(legend.title=element_text(size=15)) +
             theme(legend.text=element_text(size=15))
             
             combine.invert <- ggplot(data) +
             geom_line(aes(Wavenumber, Amplitude)) +
             theme_light()+
             theme(legend.position="bottom") +
             scale_x_reverse(expression(paste("Wavenumber (cm"^"-1"*")")), breaks=seq(0, 4000, 250)) +
             scale_y_reverse("Amplitude") +
             coord_cartesian(xlim = ranges$x, ylim = ranges$y) +
             geom_point(data=peakTable(), aes(Wavenumber, Amplitude), shape=1, size=3) +
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
            p(HTML(paste0("Wavenumber:", " ", round(point$Wavenumber, 0)))),
            p(HTML(paste0("Amplitude:", " ", round(point$Amplitude, 2))))
            )
        })
        
        output$downloadPlot <- downloadHandler(
        filename = function() { paste0(input$projectname, '.jpg', sep='') },
        content = function(file) {
            ggsave(file,plotInput(), width=14, height=8, device="jpeg")
        }
        )
        
        # Toggle points that are clicked
        observeEvent(input$plot1_click, {
            
            predict.frame <- dataManipulate()
            
            res <- nearPoints(predict.frame, input$plot1_click, allRows = TRUE, maxpoints=1)
            
            temprows <- xor(unlist(peaks$findpeaks), res$selected_)
            
            peaks$findpeaks <- relist(flesh=temprows, skeleton=peaks$findpeaks)
        })
        
        
        
        # Toggle points that are brushed, when button is clicked
        observeEvent(input$exclude_toggle, {
            
            predict.frame <- dataManipulate()
            res <- brushedPoints(predict.frame, input$plot1_brush, allRows = TRUE)
            
            temprows <- xor(unlist(peaks$findpeaks), res$selected_)
            
            peaks$findpeaks <- relist(flesh=temprows, skeleton=peaks$findpeaks)
            
        })
        
        
        
        # Reset all points
        observeEvent(input$exclude_reset, {
            
            predict.frame <- dataManipulate()

            peaks$findpeaks <- findPeaks()
        })
        
    
    
    
    waveInput <- reactive({
        
        blank.frame <- data.frame(
        Name=as.vector(as.character(rep("", 25))),
        WaveMin=as.numeric(rep("", 25)),
        WaveMax=as.numeric(rep("", 25)),
        stringsAsFactors = FALSE
        )
        
        blank.frame
        
    })
    
    
    
    wavevalues <- reactiveValues()
    

    
    waveInputCal <- reactive({

        
       calFileContents()[["Definitions"]]
       
        
    })
    
    waveTableInput <- reactive({
        
        
        if(input$usecalfile==FALSE){
            waveInput()
        }else if(input$usecalfile==TRUE){
            waveInputCal()
        }

    })
    
    
    
    
    
    
    
    
    
    
    observe({
        if (!is.null(input$hotwave)) {
            DF <- hot_to_r(input$hotwave)
        } else {
            DF <- waveTableInput()
        }
        wavevalues[["DF"]] <- DF
    })
    
    eventReactive(input$linecommitwave,{
        
        wavevalues[["DF"]] <- waveTableInput()
        
    })
    
    

    
    
    ## Handsontable
    
    output$hotwave <- renderRHandsontable({
        
        DF <- wavevalues[["DF"]]
        
        
        
        
        rhandsontable(DF) %>% hot_col(1:length(DF), renderer=htmlwidgets::JS("safeHtmlRenderer"))
        
        
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
    
    waveSubset <- reactive({
        
        ftir_parse(range.table = wavevalues[["DF"]], data=dataManipulate())

    })
    
    
    output$LineValues <- renderDataTable({
        waveSubset()
    })
    
    
    
    hotableInputBlank <- reactive({

        
        spectra.line.table <- waveSubset()
      
        empty <- spectra.line.table[,-1] * 0.0000
        empty.line.table <- data.frame(Spectrum=spectra.line.table$Spectrum, empty)

        empty.line.table
        
        
    })
    
    output$testtable <- renderDataTable({
        hotableInputBlank()
    })
    
    hotableInputCal <- reactive({
        
        calFileContents()[["Values"]]
        
        
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
    
    
    
    
    
    observe({
        if (!is.null(input$hot)) {
            DF <- hot_to_r(input$hot)
        } else {
            if (input$linecommitwave)
            DF <- hotableInput()
            else
            DF <- values[["DF"]]
        }
        values[["DF"]] <- DF
    })
    
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
    
    
    Wavenumberlinestouse <- reactive({
        
        table <- wavevalues[["DF"]]
        table <- table[complete.cases(table),]
        
        as.vector(table$Name)
        
    })
    
    
    
    
    
    outVar <- reactive({
        input$hotableprocess2
        
        myelements <- Wavenumberlinestouse()
        
        result <- if(is.null(myelements)){
            NULL
        }else{
            myelements
        }
        
        result
        
        
    })
    
    outVaralt <- reactive({
        input$hotableprocess2
        
        
        myelements <- c(Wavenumberlinestouse())
        
        
        if(is.null(myelements)){
            NULL
        }else{
            myelements
        }
        
    })
    
    outVaralt2 <- reactive({
        input$hotableprocess2
        
        
        mylines <- c(Wavenumberlinestouse())
        
        
        if(is.null(mylines)){
            NULL
        }else{
            mylines[! mylines %in% c(input$calcurveline)]
        }
        
    })
    
    output$inVar2 <- renderUI({
        selectInput(inputId = "calcurveline", label = h4("Line"), choices =  outVar())
    })
    
    inVar3Selectedpre <- reactive({
        
        hold <- values[["DF"]]
        
        optionhold <- if(is.null(input$calcurveline)){
            ls(hold)[2]
        }else{
            input$calcurveline
        }
        
        if(input$usecalfile==FALSE && is.null(calList[[optionhold]])==TRUE){
            optionhold
        } else if(input$usecalfile==TRUE && is.null(calList[[optionhold]])==TRUE  && is.null(calFileContents()$calList[[optionhold]])==FALSE){
            calFileContents()$calList[[optionhold]][[1]]$Intercept
        } else if(input$usecalfile==FALSE && is.null(calList[[optionhold]])==FALSE){
            calList[[optionhold]][[1]]$Intercept
        } else if(input$usecalfile==TRUE && is.null(calList[[optionhold]])==FALSE){
            calList[[optionhold]][[1]]$Intercept
        } else if(input$usecalfile==TRUE && is.null(calList[[optionhold]])==TRUE && is.null(calFileContents()$calList[[optionhold]])==TRUE){
            optionhold
        }
        
        
    })
    
    
    
    
    inVar4Selectedpre <- reactive({
        
        hold <- values[["DF"]]
        
        optionhold <- if(is.null(input$calcurveline)){
            ls(hold)[2]
        }else{
            input$calcurveline
        }
        
        
        if(input$usecalfile==FALSE && is.null(calList[[optionhold]])==TRUE){
            optionhold
        } else if(input$usecalfile==TRUE && is.null(calList[[optionhold]])==TRUE && is.null(calFileContents()$calList[[optionhold]])==FALSE){
            calFileContents()$calList[[optionhold]][[1]]$Slope
        } else if(input$usecalfile==FALSE && is.null(calList[[optionhold]])==FALSE){
            calList[[optionhold]][[1]]$Slope
        } else if(input$usecalfile==TRUE && is.null(calList[[optionhold]])==FALSE){
            calList[[optionhold]][[1]]$Slope
        }  else if(input$usecalfile==TRUE && is.null(calList[[optionhold]])==TRUE && is.null(calFileContents()$calList[[optionhold]])==TRUE){
            optionhold
        }
    })
    
    
    dataNorm <- reactive({
        
        dataManipulate()
        
    })
    
    
    #####Set Defaults
    
    
    calConditons <- reactiveValues()
    calList <- reactiveValues()
    calList <- NULL
    
    observeEvent(input$hotableprocess2, {
        
        cal.condition <- 3
        norm.condition <- 1
        
        norm.min <- 2000
        norm.max <- 2250
        
        cal.table <- data.frame(cal.condition, norm.condition, norm.min, norm.max)
        colnames(cal.table) <- c("CalType", "NormType", "Min", "Max")
        
        slope.corrections <- NULL
        intercept.corrections <- NULL
        
        standards.used <- vals$keeprows
        
        cal.mode.list <- list(cal.table, slope.corrections, intercept.corrections, standards.used)
        names(cal.mode.list) <- c("CalTable", "Slope", "Intercept", "StandardsUsed")
        
        calConditons <<- cal.mode.list
        
    })
    
    ########Machine Learning: Normalization
    
    normMinPre <- reactive({
        
        hold <- values[["DF"]]
        
        optionhold <- if(is.null(input$calcurveline)){
            ls(hold)[2]
        }else{
            input$calcurveline
        }
        
        if(input$usecalfile==FALSE && is.null(calList[[optionhold]])==TRUE){
            calConditons[["CalTable"]][["Min"]]
        }else if(input$usecalfile==TRUE && is.null(calFileContents()$calList[[optionhold]])==FALSE && is.null(calList[[optionhold]])==TRUE){
            calFileContents()$calList[[optionhold]][[1]]$CalTable$Min
        } else if(input$usecalfile==FALSE && is.null(calList[[optionhold]])==FALSE){
            calList[[optionhold]][[1]]$CalTable$Min
        } else if(input$usecalfile==TRUE && is.null(calList[[optionhold]])==FALSE){
            calList[[optionhold]][[1]]$CalTable$Min
        } else if(input$usecalfile==TRUE && is.null(calList[[optionhold]])==TRUE && is.null(calFileContents()$calList[[optionhold]])==FALSE){
            calConditons[["CalTable"]][["Min"]]
        }
        
    })
    
    normMaxPre <- reactive({
        
        hold <- values[["DF"]]
        
        optionhold <- if(is.null(input$calcurveline)){
            ls(hold)[2]
        }else{
            input$calcurveline
        }
        
        if(input$usecalfile==FALSE && is.null(calList[[optionhold]])==TRUE){
            calConditons[["CalTable"]][["Max"]]
        }else if(input$usecalfile==TRUE && is.null(calList[[optionhold]])==TRUE && is.null(calFileContents()$calList[[optionhold]])==FALSE){
            calFileContents()$calList[[optionhold]][[1]]$CalTable$Max
        } else if(input$usecalfile==FALSE && is.null(calList[[optionhold]])==FALSE){
            calList[[optionhold]][[1]]$CalTable$Max
        } else if(input$usecalfile==TRUE && is.null(calList[[optionhold]])==FALSE){
            calList[[optionhold]][[1]]$CalTable$Max
        } else if(input$usecalfile==TRUE && is.null(calList[[optionhold]])==TRUE && is.null(calFileContents()$calList[[optionhold]])==TRUE){
            calConditons[["CalTable"]][["Max"]]
        }
    })
    
    
    
    planktonVector <- reactive({
        
        #0.7, 0.9
        #2.5, 2.8
        #11.0, 11.2
        #18.4, 19.4
        #19.5, 22
        #21, 22
        #30, 35
        #35, 40
        
        mins <- c(3500, 3000, 2500, 2250, 2000, 1750)
        maxs <- c(4000, 3500, 3000, 2500, 2250, 2000)
        
        norm.list <- list(mins, maxs)
        names(norm.list) <- c("Min", "Max")
        norm.list
        
        
    })
    
    bestNormVars <- reactive({
        
        norm.list <- planktonVector()
        
        line <- input$calcurveline
        
        choices <- Wavenumberlinestouse()
        spectra.line.table <- spectraLineTable()
        data <- dataNorm()
        concentration.table <- concentrationTable()
        
        index <- seq(1, length(norm.list[[1]]), 1)
        
        
        concentration.table <- concentration.table[complete.cases(concentration.table[,input$calcurveline]),]
        
        
        
        spectra.line.table <- spectraLineTable()[spectraLineTable()$Spectrum %in% holdFrame()$Spectrum, ]
        
        #spectra.line.table <- spectra.line.table[spectra.line.table$Spectrum %in% concentration.table$Spectrum, ]
        
        spectra.line.table <- spectra.line.table[complete.cases(concentration.table[, line]),]
        
        data <- data[data$Spectrum %in% concentration.table$Spectrum, ]
        
        
        time.bic <- if(dataType()=="Spectra"){
            extractAIC(lm(concentration.table[, input$calcurveline]~general.prep(spectra.line.table, input$calcurveline)$Amplitude, na.action=na.exclude), k=log(length(1)))[2]
        } else if(dataType()=="Net"){
            extractAIC(lm(concentration.table[, input$calcurveline]~general.prep.net(spectra.line.table, input$calcurveline)$Amplitude, na.action=na.exclude), k=log(length(1)))[2]
        }
        
        tc.bic <- if(dataType()=="Spectra"){
            extractAIC(lm(concentration.table[, input$calcurveline]~simple.tc.prep(data, spectra.line.table, input$calcurveline)$Amplitude, na.action=na.exclude), k=log(length(1)))[2]
        } else if(dataType()=="Net"){
            extractAIC(lm(concentration.table[, input$calcurveline]~simple.tc.prep.net(data, spectra.line.table, input$calcurveline)$Amplitude, na.action=na.exclude), k=log(length(1)))[2]
        }
        
        comp.bic <- if(dataType()=="Spectra"){
            optimal_norm_chain(data=data, element=line, spectra.line.table=spectra.line.table, values=concentration.table, possible.mins=norm.list[["Min"]], possible.maxs=norm.list[["Max"]])
        } else if(dataType()=="Net"){
            time.bic
        }
        
        norm.chain <- c(time.bic, tc.bic, comp.bic)
        type.chain <- c(1, 2, 3)
        
        best <- index[[which.min(unlist(norm.chain))]]
        best.comp <- c(planktonVector()[["Min"]][best], planktonVector()[["Max"]][best])
        best.type <- type.chain[which.min(unlist(norm.chain))]
        result.list <- list(best.type, best.comp)
        names(result.list) <- c("Type", "Compton")
        result.list
    })
    
    
    normhold <- reactiveValues()
    
    observeEvent(input$hotableprocess2, {
        normhold$norms <- c(normMinPre(), normMaxPre())
        normhold$normtype <- calNormSelectionpre()
    })
    
    
    observeEvent(input$trainslopes, {
        
        isolate(normhold$norms[1] <- bestNormVars()[["Compton"]][1])
        isolate(normhold$norms[2] <- bestNormVars()[["Compton"]][2])
        isolate(normhold$normtype <- bestNormVars()[["Type"]])
        
    })
    
    calNormSelection <- reactive({
        normhold$normtype
    })
    
    normMinSelection <- reactive({
        normhold$norms[1]
    })
    
    normMaxSelection <- reactive({
        normhold$norms[2]
    })
    
    
    
    
    calNormSelectionpre <- reactive({
        
        hold <- values[["DF"]]
        
        optionhold <- if(is.null(input$calcurveline)){
            ls(hold)[2]
        }else{
            input$calcurveline
        }
        
        if(input$usecalfile==FALSE && is.null(calList[[optionhold]])==TRUE){
            calConditons[["CalTable"]][["NormType"]]
        }else if(input$usecalfile==TRUE && is.null(calList[[optionhold]])==TRUE && is.null(calFileContents()$calList[[optionhold]])==FALSE){
            calFileContents()$calList[[optionhold]][[1]]$CalTable$NormType
        } else if(input$usecalfile==FALSE && is.null(calList[[optionhold]])==FALSE){
            calList[[optionhold]][[1]]$CalTable$NormType
        } else if(input$usecalfile==TRUE && is.null(calList[[optionhold]])==FALSE){
            calList[[optionhold]][[1]]$CalTable$NormType
        } else if(input$usecalfile==TRUE && is.null(calList[[optionhold]])==TRUE && is.null(calFileContents()$calList[[optionhold]])==TRUE){
            calConditons[["CalTable"]][["NormType"]]
        }
        
    })
    
    
    output$normTypeInput <- renderUI({
        
        selectInput("normcal", label = "Normalization",
        choices = list("Time" = 1, "Total Counts" = 2, "Compton" = 3),
        selected = calNormSelection())
        
        
    })
    
    
    output$comptonMinInput <- renderUI({
        
        numericInput('comptonmin', label=h6("Min"), step=1, value=normMinSelection(), min=400, max=4000, width='30%')
        
    })
    
    output$comptonMaxInput <- renderUI({
        
        numericInput('comptonmax', label=h6("Max"), step=1, value=normMaxSelection(), min=400, max=4000, width='30%')
        
    })
    
    
    #####Machine Learning: Intercepts
    
    
    
    cephlopodVector <- reactive({
        
        combos_mod <- function(a.vector){
            
            so <- seq(from=1, to=length(a.vector), by=1)
            
            long <- pblapply(so, function(x) gRbase::combnPrim(x=a.vector, m=x), cl=6L)
            and <- pblapply(long, function(x) plyr::alply(x, 2), cl=6L)
            thanks.for.all.the.fish <- do.call(list, unlist(and, recursive=FALSE))
            
            thanks.for.all.the.fish
            
        }
        
        if(!is.null(likely_intercepts(input$calcurveline))){
            combos_mod(likely_intercepts(input$calcurveline))
        } else if(is.null(likely_intercepts(input$calcurveline))){
            c("Rh.K.alpha", "Rh.L.alpha")
        }
        
    })
    
    
    bestInterceptVars <- reactive({
        
        line <- input$calcurveline
        
        choices <- Wavenumberlinestouse()
        
        spectra.line.table <- if(all(cephlopodVector() %in% colnames(spectraLineTable()))==TRUE){
            spectraLineTable()
        } else if(all(cephlopodVector() %in% colnames(spectraLineTable()))==FALSE){
            merge(spectraLineTable(), elementFrame(data=dataHold(), elements=cephlopodVector()[cephlopodVector() %in% colnames(spectraLineTable())]))
        }
        
        data <- dataNorm()
        concentration.table <- concentrationTable()
        
        
        spectra.line.table <- spectra.line.table[spectra.line.table$Spectrum %in% holdFrame()$Spectrum, ]
        
        
        predict.amplitude.list <- if(input$normcal==1){
            pblapply(cephlopodVector(), function(x) lucas.simp.prep(spectra.line.table=spectra.line.table, element.line=line, slope.element.lines=line, intercept.element.lines=c(line, x)))
        } else if(input$normcal==2){
            pblapply(cephlopodVector(), function(x) lucas.tc.prep(data=data, spectra.line.table=spectra.line.table, element.line=line, slope.element.lines=line, intercept.element.lines=c(line, x)))
        } else if(input$normcal==3){
            pblapply(cephlopodVector(), function(x) lucas.comp.prep(data=data, spectra.line.table=spectra.line.table, element.line=line, slope.element.lines=line, intercept.element.lines=c(line, x), norm.min=input$comptonmin, norm.max=input$comptonmax))
        }
        
        optimal_intercept_chain(element=line, intensities=predict.amplitude.list, values=concentration.table, keep=vals$keeprows)
        
        
    })
    
    
    intercepthold <- reactiveValues()
    intercepthold$intercepts <- NULL
    
    observeEvent(input$calcurveline, {
        
        isolate(intercepthold$intercepts <- inVar3Selectedpre())
        
    })
    
    
    #observeEvent(input$trainslopes, {
    
    #    isolate(intercepthold$intercepts <- bestInterceptVars())
    
    #})
    
    
    inVar3Selected <- reactive({
        
        intercepthold$intercepts
        
        
    })
    
    output$inVar3 <- renderUI({
        
        selectInput(inputId = "intercept_vars", label = h4("Intercept"), choices =  outVaralt2(), selected=inVar3Selected(), multiple=TRUE)
    })
    
    
    
    ####Machine Learning: Slopes
    
    caretSlope <- reactive({
        
        line <- input$calcurveline
        
        # prepare simple test suite
        control <- trainControl(method="cv", number=5)
        seed <- 7
        metric <- "RMSE"
        set.seed(seed)
        
        data <- dataNorm()
        concentration.table <- concentrationTable()
        
        concentration.table <- concentration.table[complete.cases(concentration.table[,input$calcurveline]),]
        
        choices <- Wavenumberlinestouse()
        
        
        spectra.line.table <- spectraLineTable()[spectraLineTable()$Spectrum %in% holdFrame()$Spectrum, ]
        
        #spectra.line.table <- spectra.line.table[spectra.line.table$Spectrum %in% concentration.table$Spectrum, ]
        
        spectra.line.table <- spectra.line.table[complete.cases(concentration.table[, line]),]
        
        data <- data[data$Spectrum %in% concentration.table$Spectrum, ]
        
        
        cal.table <- if(dataType()=="Spectra"){
            if(input$normcal==1){
                lucas.simp.prep(spectra.line.table=spectra.line.table, element.line=input$calcurveelemenet,slope.element.lines=choices, intercept.element.lines=input$intercept_vars)
            } else if(input$normcal==2){
                lucas.tc.prep(data=data, spectra.line.table=spectra.line.table, element.line=input$calcurveline, slope.element.lines=choices, intercept.element.lines=input$intercept_vars)
            } else if(input$normcal==3){
                lucas.comp.prep(data=data, spectra.line.table=spectra.line.table, element.line=input$calcurveline, slope.element.lines=choices, intercept.element.lines=input$intercept_vars, norm.min=input$comptonmin, norm.max=input$comptonmax)
            }
        } else if(dataType()=="Net"){
            if(input$normcal==1){
                lucas.simp.prep.net(spectra.line.table=spectra.line.table, element.line=input$calcurveelemenet,slope.element.lines=choices, intercept.element.lines=input$intercept_vars)
            } else if(input$normcal==2){
                lucas.tc.prep.net(data=data, spectra.line.table=spectra.line.table, element.line=input$calcurveline, slope.element.lines=choices, intercept.element.lines=input$intercept_vars)
            } else if(input$normcal==3){
                lucas.comp.prep.net(data=data, spectra.line.table=spectra.line.table, element.line=input$calcurveline, slope.element.lines=choices, intercept.element.lines=input$intercept_vars, norm.min=input$comptonmin, norm.max=input$comptonmax)
            }
        }
        
        #cal.table <- cal.table[,!colnames(cal.table) %in% "Amplitude"]
        cal.table$Concentration <- concentration.table[,input$calcurveline]
        
        
        train(Concentration~., data=cal.table[,-1], method="lm", metric=metric, preProc=c("center", "scale"), trControl=control)
        
        
    })
    
    
    
    slopeImportance <- reactive({
        
        varImp(caretSlope(), scale=FALSE)
        
    })
    
    output$importanceplot <- renderPlot({
        
        plot(slopeImportance())
        
    })
    
    
    
    fishVector <- reactive({
        
        combos_mod <- function(a.vector){
            
            so <- seq(from=2, to=input$nvariables, by=1)
            
            long <- pblapply(so, function(x) gRbase::combnPrim(x=a.vector, m=x), cl=6L)
            and <- pblapply(long, function(x) plyr::alply(x, 2), cl=6L)
            thanks.for.all.the.fish <- do.call(list, unlist(and, recursive=FALSE))
            thanks.for.all.the.fish <- pblapply(thanks.for.all.the.fish, function(x) c(input$calcurveline, x))
            
            thanks.for.all.the.fish
            
        }
        
        fit.lm <- caretSlope()
        
        first.combos <- c(Wavenumberlinestouse()[!Wavenumberlinestouse() %in% input$calcurveline])
        
        coef.frame <- as.data.frame(summary(fit.lm)$coefficients)
        sig.frame <- subset(coef.frame, coef.frame[,4] < 0.05)
        
        second.combos <- first.combos[c(first.combos %in% rownames(sig.frame)[c(!rownames(sig.frame) %in% "(Intercept)")])]
        
        
        
        combos_mod(second.combos)
        
        
    })
    
    
    bestSlopeVars <- reactive({
        
        line <- input$calcurveline
        
        choices <- Wavenumberlinestouse()
        spectra.line.table <- spectraLineTable()
        data <- dataNorm()
        concentration.table <- concentrationTable()
        
        #concentration.table[complete.cases(concentration.table[,input$calcurveline]),]
        
        #index <- complete.cases(concentration.table[,input$calcurveline])
        
        
        spectra.line.table <- spectraLineTable()[spectraLineTable()$Spectrum %in% holdFrame()$Spectrum, ]
        
        #spectra.line.table <- spectra.line.table[spectra.line.table$Spectrum %in% concentration.table$Spectrum, ]
        
        spectra.line.table <- spectra.line.table[complete.cases(concentration.table[, line]),]
        
        data <- data[data$Spectrum %in% concentration.table$Spectrum, ]
        
        
        
        
        predict.amplitude <- if(input$normcal==1){
            if(dataType()=="Spectra"){
                lucas.simp.prep(spectra.line.table=spectra.line.table, element.line=line, slope.element.lines=choices, intercept.element.lines=input$intercept_vars)
            } else if(dataType()=="Net"){
                lucas.simp.prep.net(spectra.line.table=spectra.line.table, element.line=line, slope.element.lines=choices, intercept.element.lines=input$intercept_vars)
            }
        } else if(input$normcal==2){
            predict.amplitude <- if(dataType()=="Spectra"){
                lucas.tc.prep(data=data, spectra.line.table=spectra.line.table, element.line=line, slope.element.lines=choices, intercept.element.lines=input$intercept_vars)
            } else if(dataType()=="Net"){
                lucas.tc.prep.net(data=data, spectra.line.table=spectra.line.table, element.line=line, slope.element.lines=choices, intercept.element.lines=input$intercept_vars)
            }
        } else if(input$normcal==3){
            predict.amplitude <- if(dataType()=="Spectra"){
                lucas.comp.prep(data=data, spectra.line.table=spectra.line.table, element.line=line, slope.element.lines=choices, intercept.element.lines=input$intercept_vars, norm.min=input$comptonmin, norm.max=input$comptonmax)
            } else if(dataType()=="Net"){
                lucas.comp.prep.net(data=data, spectra.line.table=spectra.line.table, element.line=line, slope.element.lines=choices, intercept.element.lines=input$intercept_vars, norm.min=input$comptonmin, norm.max=input$comptonmax)
            }
        }
        
        
        
        #optimal_r_chain(element=element, intensities=predict.amplitude, values= concentration.table, possible.slopes=fishVector(), keep=vals$keeprows)
        
        results <- variable_select_short(slopeImportance())
        
        c(input$calcurveline, results[!results %in% input$calcurveline])
        
        
    })
    
    slopehold <- reactiveValues()
    slopehold$slopes <- NULL
    
    observeEvent(input$calcurveline, {
        
        isolate(slopehold$slopes <- inVar4Selectedpre())
        
    })
    observeEvent(input$trainslopes, {
        
        isolate(slopehold$slopes <- bestSlopeVars())
        
    })
    
    
    inVar4Selected <- reactive({
        
        slopehold$slopes
        
        
    })
    
    
    
    output$inVar4 <- renderUI({
        
        selectInput(inputId = "slope_vars", label = h4("Slope"), choices =  outVaralt(), selected=inVar4Selected(), multiple=TRUE)
    })
    
    
    
    
    #####Machine Learning: Cal Type
    
    
    calTypeSelectionPre <- reactive({
        
        hold <- values[["DF"]]
        
        optionhold <- if(is.null(input$calcurveline)){
            ls(hold)[2]
        }else{
            input$calcurveline
        }
        
        if(input$usecalfile==FALSE && is.null(calList[[optionhold]])==TRUE){
            calConditons[["CalTable"]][["CalType"]]
        } else if(input$usecalfile==TRUE && is.null(calList[[optionhold]])==TRUE && is.null(calFileContents()$calList[[optionhold]])==FALSE){
            calFileContents()$calList[[optionhold]][[1]]$CalTable$CalType
        } else if(input$usecalfile==FALSE && is.null(calList[[optionhold]])==FALSE){
            calList[[optionhold]][[1]]$CalTable$CalType
        } else if(input$usecalfile==TRUE && is.null(calList[[optionhold]])==FALSE){
            calList[[optionhold]][[1]]$CalTable$CalType
        } else if(input$usecalfile==TRUE && is.null(calList[[optionhold]])==TRUE && is.null(calFileContents()$calList[[optionhold]])==TRUE){
            calConditons[["CalTable"]][["CalType"]]
        }
        
    })
    
    
    
    
    
    
    bestCalType <- reactive({
        
        concentration.table <- concentrationTable()
        data <- dataNorm()
        spectra.line.table <- spectraLineTable()
        
        
        #concentration.table <- concentration.table[complete.cases(concentration.table[, input$calcurveline]),]
        
        #spectra.line.table <- spectra.line.table[complete.cases(concentration.table[, input$calcurveline]),]
        #data <- data[data$Spectrum %in% concentration.table$Spectrum, ]
        
        predict.amplitude <- if(input$normcal==1){
            if(dataType()=="Spectra"){
                lucas.simp.prep(spectra.line.table=spectra.line.table, element.line=input$calcurveline, slope.element.lines=Wavenumberlinestouse(), intercept.element.lines=input$intercept_vars)
            } else if(dataType()=="Net"){
                lucas.simp.prep.net(spectra.line.table=spectra.line.table, element.line=input$calcurveline, slope.element.lines=Wavenumberlinestouse(), intercept.element.lines=input$intercept_vars)
            }
        } else if(input$normcal==2){
            predict.amplitude <- if(dataType()=="Spectra"){
                lucas.tc.prep(data=data, spectra.line.table=spectra.line.table, element.line=input$calcurveline, slope.element.lines=Wavenumberlinestouse(), intercept.element.lines=input$intercept_vars)
            } else if(dataType()=="Net"){
                lucas.tc.prep.net(data=data, spectra.line.table=spectra.line.table, element.line=input$calcurveline, slope.element.lines=Wavenumberlinestouse(), intercept.element.lines=input$intercept_vars)
            }
        } else if(input$normcal==3){
            predict.amplitude <- if(dataType()=="Spectra"){
                lucas.comp.prep(data=data, spectra.line.table=spectra.line.table, element.line=input$calcurveline, slope.element.lines=Wavenumberlinestouse(), intercept.element.lines=input$intercept_vars, norm.min=input$comptonmin, norm.max=input$comptonmax)
            } else if(dataType()=="Net"){
                lucas.comp.prep.net(data=data, spectra.line.table=spectra.line.table, element.line=input$calcurveline, slope.element.lines=Wavenumberlinestouse(), intercept.element.lines=input$intercept_vars, norm.min=input$comptonmin, norm.max=input$comptonmax)
            }
        }
        
        
        
        predict.frame <- data.frame(predict.amplitude, concentration.table[,input$calcurveline])
        predict.frame <- predict.frame[complete.cases(predict.frame),]
        predict.frame <- predict.frame[vals$keeprows,]
        colnames(predict.frame) <- c(names(predict.amplitude), "Concentration")
        predict.frame <- predict.frame[complete.cases(predict.frame$Concentration),]
        
        
        predict.frame.simp <- predict.frame[,c("Concentration", "Amplitude")]
        predict.frame.luc <- predict.frame[, c("Concentration", "Amplitude", input$slope_vars)]
        predict.frame.forest <- predict.frame
        
        
        cal.lm.simp <- lm(Concentration~Amplitude, data=predict.frame.simp, na.action=na.exclude)
        
        cal.lm.two <- lm(Concentration~Amplitude + I(Amplitude^2), data=predict.frame.simp, na.action=na.exclude)
        
        cal.lm.luc <- lm(Concentration~., data=predict.frame.luc, na.action=na.exclude)
        
        cal.lm.forest <- randomForest(Concentration~., data=predict.frame.forest, na.action=na.omit)
        
        forest.predict <- predict(cal.lm.forest, new.data=predict.frame.forest, proximity=FALSE)
        forest.sum <- lm(predict.frame$Concentration~forest.predict, na.action=na.exclude)
        
        
        r2.vector <- c(summary(cal.lm.simp)$adj.r.squared, summary(cal.lm.two)$adj.r.squared-.5, summary(cal.lm.luc)$adj.r.squared, summary(forest.sum)$adj.r.squared)
        which.max(r2.vector)
        
        
    })
    
    testing2 <- reactive({
        
        predict.frame <- predictFrame()
        predict.frame <- predict.frame[complete.cases(predict.frame$Concentration),]
        
        predict.frame
        
    })
    
    
    testing <- reactive({
        predict.frame <- predictFrame()[complete.cases(predictFrame()$Concentration),]
        cal.lm.forest <- randomForest(Concentration~., data=predict.frame, na.action=na.omit)
        
        forest.predict <- predict(cal.lm.forest, new.data=predict.frame, proximity=FALSE)
        as.data.frame(forest.predict)
        
    })
    
    output$testtable <- renderDataTable({
        
        testing()
        
    })
    
    output$testtable2 <- renderDataTable({
        
        testing2()
        
    })
    
    
    calhold <- reactiveValues()
    
    observeEvent(input$hotableprocess2, {
        calhold$caltype <- calTypeSelectionPre()
    })
    
    
    observeEvent(input$trainslopes, {
        
        isolate(calhold$caltype <- bestCalType())
        
    })
    
    
    
    calTypeSelection <- reactive({
        calhold$caltype
    })
    
    
    output$calTypeInput <- renderUI({
        
        selectInput("radiocal", label = "Calibration Curve",
        choices = list("Linear" = 1, "Non-Linear" = 2, "Lucas-Tooth" = 3, "Forest" = 4),
        selected = calTypeSelection())
        
        
    })
    
    
    

    
    
    
    
    
    
    
    lineHold <- reactive({
        
        if(is.null(input$calcurveline)==TRUE){
            ls(dataHold())[1]
        } else{
            input$calcurveline
        }
        
    })
    
    
    
    
    calFileStandards <- reactive({
        
        
        
        if(input$usecalfile==TRUE && is.null(calList[[lineHold()]])==TRUE && is.null(calFileContents()$calList[[lineHold()]])==FALSE){
            calFileContents()$calList[[lineHold()]][[1]][[4]]
        } else if(input$usecalfile==FALSE && is.null(calList[[lineHold()]])==TRUE && is.null(calFileContents()$calList[[lineHold()]])==TRUE){
            rep(TRUE, dataCount())
        } else if(input$usecalfile==TRUE && is.null(calList[[lineHold()]])==FALSE && is.null(calFileContents()$calList[[lineHold()]])==TRUE){
            calList[[lineHold()]][[1]][[4]]
        } else if(input$usecalfile==TRUE && is.null(calList[[lineHold()]])==FALSE && is.null(calFileContents()$calList[[lineHold()]])==FALSE){
            calList[[lineHold()]][[1]][[4]]
        } else if(input$usecalfile==FALSE && is.null(calList[[lineHold()]])==FALSE && is.null(calFileContents()$calList[[lineHold()]])==TRUE){
            calList[[lineHold()]][[1]][[4]]
        } else if(input$usecalfile==FALSE && is.null(calList[[lineHold()]])==FALSE && is.null(calFileContents()$calList[[lineHold()]])==FALSE){
            calList[[lineHold()]][[1]][[4]]
        } else if(input$usecalfile==TRUE && is.null(calList[[lineHold()]])==TRUE && is.null(calFileContents()$calList[[lineHold()]])==TRUE){
            rep(TRUE, dataCount())
        }
        
        
        
        
    })
    
    
    
    
    vals <- reactiveValues()
    
    
    vals$keeprows <- if(input$usecalfile==TRUE){
        calFileStandards()
    }else{
        rep(TRUE, dataCount())
    }
    
    keepRowsFrame <- reactive({
        
        spectra.stuff <- values[["DF"]]
        rows <- vals$keeprows
        
        the.frame <- data.frame(Spectrum=spectra.stuff$Spectrum, Standards=rows)
        the.frame
        
    })
    
    
    output$whichrowstokeep <- renderRHandsontable({
        
        DF <- keepRowsFrame()
        
        DF <- DF[order(as.character(DF$Spectrum)),]
        
        
        
        if (!is.null(DF))
        rhandsontable(DF) %>% hot_col(2:length(DF), renderer=htmlwidgets::JS("safeHtmlRenderer"))
        
        
    })
    
    
    
    observe({
        if (!is.null(input$hot)) {
            DF <- hot_to_r(input$hot)
        } else {
            if (input$linecommitwave)
            DF <- hotableInput()
            else
            DF <- values[["DF"]]
        }
        values[["DF"]] <- DF
    })
    
    
    
    #if(input$hotableprocess2){vals$keeprows <- vals$keeprows[dropStandard()]}
    
    
    output$temp <- renderTable({
        
        as.data.frame(vals$keeprows)
        
    })
    
    
    dataType <- reactive({
        if(input$filetype=="DPT"){
            "Spectra"
        } else if(input$filetype=="Opus"){
            "Spectra"
        } else if(input$filetype=="CSV"){
            "Spectra"
        }
        
    })
    
    
    
    concentrationTable <- reactive({
        
        concentration.table <- as.data.frame(values[["DF"]], stringsAsFactors=FALSE)
        concentration.table[concentration.table==""] <- NA
        concentration.table[values[["DF"]]$Include,]
        
    })
    
    spectraLineTable <- reactive({
        
        spectra.line.table <- if(dataType()=="Spectra"){
            waveSubset()[values[["DF"]]$Include,]
        }else if(dataType()=="Net"){
            dataHold()[values[["DF"]]$Include,]
        }
        
        
        spectra.line.table <- spectra.line.table[order(as.character(spectra.line.table$Spectrum)),]
        spectra.line.table
        
        
        
    })
    
    
    holdFrame <- reactive({
        
        spectra.line.table <- spectraLineTable()
        
        concentration <- as.vector(as.numeric(unlist(concentrationTable()[,input$calcurveline])))
        
        
        
        Amplitude <- as.vector(as.numeric(unlist(spectraLineTable()[,input$calcurveline])))
        
        spectra.names <- spectra.line.table$Spectrum
        
        hold.frame <- data.frame(spectra.names, concentration, Amplitude)
        colnames(hold.frame) <- c("Spectrum", "Concentration", "Amplitude")
        hold.frame <- na.omit(hold.frame)
        
        hold.frame <- hold.frame[order(as.character(hold.frame$Spectrum)),]
        
        
        hold.frame[complete.cases(hold.frame),]
        
        
    })
    
    dataNorm <- reactive({
        
        data <- dataHold()
        data[data$Spectrum %in% holdFrame()$Spectrum, ]
        
        
    })
    
    
    predictFramePre <- reactive({
        
        Amplitude <- holdFrame()$Amplitude
        
        concentration <- holdFrame()$Concentration
        
        predict.frame <- data.frame(concentration, Amplitude)
        colnames(predict.frame) <- c("Concentration", "Amplitude")
        
        
        predict.frame
        
        
    })
    
    
    
    predictAmplitude <- reactive({
        
        spectra.line.table <- spectraLineTable()
        data <- dataNorm()
        
        
        spectra.line.table <- spectraLineTable()[spectraLineTable()$Spectrum %in% holdFrame()$Spectrum, ]
        
        lines <- names(spectra.line.table[,-1])
        
        
        if (input$radiocal==1){
            
            if(input$normcal==1){
                predict.amplitude <- if(dataType()=="Spectra"){
                    general.prep(spectra.line.table=spectra.line.table, element.line=input$calcurveline)
                } else if(dataType()=="Net"){
                    general.prep.net(spectra.line.table=spectra.line.table, element.line=input$calcurveline)
                }
            }
            
            if(input$normcal==2){
                predict.amplitude <- if(dataType()=="Spectra"){
                    simple.tc.prep(data=data, spectra.line.table=spectra.line.table, element.line=input$calcurveline)
                } else if(dataType()=="Net"){
                    simple.tc.prep.net(data=data, spectra.line.table=spectra.line.table, element.line=input$calcurveline)
                }
            }
            
            if(input$normcal==3){
                predict.amplitude <- if(dataType()=="Spectra"){
                    simple.comp.prep(data=data, spectra.line.table=spectra.line.table, element.line=input$calcurveline, norm.min=input$comptonmin, norm.max=input$comptonmax)
                } else if(dataType()=="Net"){
                    simple.comp.prep.net(data=data, spectra.line.table=spectra.line.table, element.line=input$calcurveline, norm.min=input$comptonmin, norm.max=input$comptonmax)
                }
            }
            
        }
        
        if (input$radiocal==2){
            
            if(input$normcal==1){
                predict.amplitude <- if(dataType()=="Spectra"){
                    general.prep(spectra.line.table=spectra.line.table, element.line=input$calcurveline)
                } else if(dataType()=="Net"){
                    general.prep.net(spectra.line.table=spectra.line.table, element.line=input$calcurveline)
                }
            }
            
            if(input$normcal==2){
                predict.amplitude <- if(dataType()=="Spectra"){
                    simple.tc.prep(data=data, spectra.line.table=spectra.line.table, element.line=input$calcurveline)
                } else if(dataType()=="Net"){
                    simple.tc.prep.net(data=data, spectra.line.table=spectra.line.table, element.line=input$calcurveline)
                }
            }
            
            if(input$normcal==3){
                predict.amplitude <- if(dataType()=="Spectra"){
                    simple.comp.prep(data=data, spectra.line.table=spectra.line.table, element.line=input$calcurveline, norm.min=input$comptonmin, norm.max=input$comptonmax)
                } else if(dataType()=="Net"){
                    simple.comp.prep.net(data=data, spectra.line.table=spectra.line.table, element.line=input$calcurveline, norm.min=input$comptonmin, norm.max=input$comptonmax)
                }
            }
            
        }
        
        
        if (input$radiocal==3){
            
            if(input$normcal==1){
                predict.amplitude <- if(dataType()=="Spectra"){
                    lucas.simp.prep(spectra.line.table=spectra.line.table, element.line=input$calcurveline, slope.element.lines=input$slope_vars, intercept.element.lines=input$intercept_vars)
                } else if(dataType()=="Net"){
                    lucas.simp.prep.net(spectra.line.table=spectra.line.table, element.line=input$calcurveline, slope.element.lines=input$slope_vars, intercept.element.lines=input$intercept_vars)
                }
            }
            
            if(input$normcal==2){
                predict.amplitude <- if(dataType()=="Spectra"){
                    lucas.tc.prep(data=data, spectra.line.table=spectra.line.table, element.line=input$calcurveline, slope.element.lines=input$slope_vars, intercept.element.lines=input$intercept_vars)
                } else if(dataType()=="Net"){
                    lucas.tc.prep.net(data=data, spectra.line.table=spectra.line.table, element.line=input$calcurveline, slope.element.lines=input$slope_vars, intercept.element.lines=input$intercept_vars)
                }
            }
            
            if(input$normcal==3){
                predict.amplitude <- if(dataType()=="Spectra"){
                    lucas.comp.prep(data=data, spectra.line.table=spectra.line.table, element.line=input$calcurveline, slope.element.lines=input$slope_vars, intercept.element.lines=input$intercept_vars, norm.min=input$comptonmin, norm.max=input$comptonmax)
                } else if(dataType()=="Net"){
                    lucas.comp.prep.net(data=data, spectra.line.table=spectra.line.table, element.line=input$calcurveline, slope.element.lines=input$slope_vars, intercept.element.lines=input$intercept_vars, norm.min=input$comptonmin, norm.max=input$comptonmax)
                }
            }
            
        }
        
        if (input$radiocal==4){
            concentration.table <- concentrationTable()
            concentration.table <- concentration.table[complete.cases(concentration.table[, input$calcurveline]),]
            spectra.line.table <- spectra.line.table[complete.cases(concentration.table[, input$calcurveline]),]
            data <- data[data$Spectrum %in% concentration.table$Spectrum, ]
            
            predict.amplitude <- if(input$normcal==1){
                if(dataType()=="Spectra"){
                    lucas.simp.prep(spectra.line.table=spectra.line.table, element.line=input$calcurveline, slope.element.lines=lines, intercept.element.lines=input$intercept_vars)
                } else if(dataType()=="Net"){
                    lucas.simp.prep.net(spectra.line.table=spectra.line.table, element.line=input$calcurveline, slope.element.lines=lines, intercept.element.lines=input$intercept_vars)
                }
            } else if(input$normcal==2){
                predict.amplitude <- if(dataType()=="Spectra"){
                    lucas.tc.prep(data=data, spectra.line.table=spectra.line.table, element.line=input$calcurveline, slope.element.lines=lines, intercept.element.lines=input$intercept_vars)
                } else if(dataType()=="Net"){
                    lucas.tc.prep.net(data=data, spectra.line.table=spectra.line.table, element.line=input$calcurveline, slope.element.lines=lines, intercept.element.lines=input$intercept_vars)
                }
            } else if(input$normcal==3){
                predict.amplitude <- if(dataType()=="Spectra"){
                    lucas.comp.prep(data=data, spectra.line.table=spectra.line.table, element.line=input$calcurveline, slope.element.lines=lines, intercept.element.lines=input$intercept_vars, norm.min=input$comptonmin, norm.max=input$comptonmax)
                } else if(dataType()=="Net"){
                    lucas.comp.prep.net(data=data, spectra.line.table=spectra.line.table, element.line=input$calcurveline, slope.element.lines=lines, intercept.element.lines=input$intercept_vars, norm.min=input$comptonmin, norm.max=input$comptonmax)
                }
            }
            
        }
        
        
        predict.amplitude
        
        
    })
    
    
    predictFrame <- reactive({
        
        predict.frame <- predictFramePre()
        predict.amplitude <- predictAmplitude()
        
        
        
        predict.frame <- data.frame(predict.amplitude, predict.frame$Concentration)
        colnames(predict.frame) <- c(names(predict.amplitude), "Concentration")
        
        
        
        predict.frame
        
        
    })
    
    
    predictFrameName <- reactive({
        
        predict.frame <- predictFrame()[ vals$keeprows, , drop = FALSE]
        spectra.line.table <- spectraLineTable()[ vals$keeprows, , drop = FALSE]
        
        predict.frame.name <- data.frame(spectra.line.table$Spectrum, predict.frame)
        colnames(predict.frame.name) <- c("Spectrum", names(predict.frame))
        predict.frame.name
        
    })
    
    
    
    calCurveFrame <- reactive({
        
        predictFrame()
        
    })
    
    
    lineModel <- reactive({
        
        predict.frame <- predictFrame()
        
        
        if (input$radiocal==1){
            cal.lm <- lm(Concentration~Amplitude, data=predict.frame[ vals$keeprows, , drop = FALSE])
        }
        
        
        if (input$radiocal==2){
            cal.lm <- lm(Concentration~Amplitude + I(Amplitude^2), data=predict.frame[ vals$keeprows, , drop = FALSE])
        }
        
        if (input$radiocal==3){
            cal.lm <- lm(Concentration~., data=predict.frame[ vals$keeprows, , drop = FALSE])
        }
        
        if (input$radiocal==4){
            cal.lm <- randomForest(Concentration~., data=predict.frame[ vals$keeprows, , drop = FALSE], na.action=na.omit)
        }
        
        cal.lm
        
    })
    
    
    valFrame <- reactive({
        
        predict.amplitude <- predictAmplitude()
        predict.frame <- predictFrame()
        line.model <- lineModel()
        
        
        if (input$radiocal==1){
            cal.est.conc.pred <- predict(object=line.model, newdata=predict.amplitude, interval='confidence')
            cal.est.conc.tab <- data.frame(cal.est.conc.pred)
            cal.est.conc <- cal.est.conc.tab$fit
            
            val.frame <- data.frame(predict.frame$Concentration, cal.est.conc)
            colnames(val.frame) <- c("Concentration", "Prediction")
        }
        
        if (input$radiocal==2){
            cal.est.conc.pred <- predict(object=line.model, newdata=predict.amplitude, interval='confidence')
            cal.est.conc.tab <- data.frame(cal.est.conc.pred)
            cal.est.conc <- cal.est.conc.tab$fit
            
            val.frame <- data.frame(predict.frame$Concentration, cal.est.conc)
            colnames(val.frame) <- c("Concentration", "Prediction")
        }
        
        if (input$radiocal==3){
            
            xmin = 0; xmax=10
            N = length(predict.frame$Concentration)
            means = colMeans(predict.frame)
            dummyDF = t(as.data.frame(means))
            for(i in 2:N){dummyDF=rbind(dummyDF,means)}
            xv=seq(xmin,xmax, length.out=N)
            dummyDF$Concentration = xv
            yv=predict(line.model, newdata=predict.amplitude)
            
            
            lucas.x <- yv
            
            cal.est.conc.pred.luc <- predict(object=line.model , newdata=predict.amplitude, interval='confidence')
            cal.est.conc.tab <- data.frame(cal.est.conc.pred.luc)
            cal.est.conc.luc <- cal.est.conc.tab$fit
            cal.est.conc.luc.up <- cal.est.conc.tab$upr
            cal.est.conc.luc.low <- cal.est.conc.tab$lwr
            
            
            val.frame <- data.frame(predict.frame$Concentration, predict.amplitude$Amplitude, lucas.x, cal.est.conc.luc, cal.est.conc.luc.up, cal.est.conc.luc.low)
            colnames(val.frame) <- c("Concentration", "Amplitude", "AmplitudeNorm", "Prediction", "Upper", "Lower")
        }
        
        if (input$radiocal==4){
            
            xmin = 0; xmax=10
            N = length(predict.frame$Concentration)
            means = colMeans(predict.frame)
            dummyDF = t(as.data.frame(means))
            for(i in 2:N){dummyDF=rbind(dummyDF,means)}
            xv=seq(xmin,xmax, length.out=N)
            dummyDF$Concentration = xv
            yv=predict(line.model, newdata=predict.amplitude)
            
            
            lucas.x <- yv
            
            cal.est.conc.pred.luc <- predict(object=line.model , newdata=predict.amplitude, interval='confidence')
            #cal.est.conc.tab <- data.frame(cal.est.conc.pred.luc)
            #cal.est.conc.luc <- cal.est.conc.tab$fit
            #cal.est.conc.luc.up <- cal.est.conc.tab$upr
            #cal.est.conc.luc.low <- cal.est.conc.tab$lwr
            
            
            val.frame <- data.frame(predict.frame$Concentration, predict.amplitude$Amplitude, lucas.x, as.vector(cal.est.conc.pred.luc))
            colnames(val.frame) <- c("Concentration", "Amplitude", "AmplitudeNorm", "Prediction")
        }
        
        
        
        
        val.frame
        
    })
    
    calValFrame <- reactive({
        
        valFrame()
        
    })
    
    
    calType <- reactive({
        
        
        if(input$radiocal==1){
            1
        } else if(input$radiocal==2){
            2
        } else if(input$radiocal==3){
            3
        } else if(input$radiocal==4){
            3
        }
        
    })
    
    rangescalcurve <- reactiveValues(x = NULL, y = NULL)
    
    
    
    calCurvePlot <- reactive({
        
        predict.frame <- predictFrame()
        line.model <- lineModel()
        val.frame <- valFrame()
        
        element.name <- gsub("[.]", "", substr(input$calcurveline, 1, 2))
        intens <- " Counts per Second"
        norma <- " Normalized"
        norma.comp <- " Compton Normalized"
        norma.tc <- " Valid Counts Normalized"
        conen <- " (%)"
        predi <- " Estimate (%)"
        log <- "Log "
        
        
        amplitude.name <- c(element.name, intens)
        concentration.name <- c(element.name, conen)
        prediction.name <- c(element.name, predi)
        
        
        if(input$radiocal==1){
            calcurve.plot <- ggplot(data=predict.frame[ vals$keeprows, , drop = FALSE], aes(Amplitude, Concentration)) +
            theme_light() +
            annotate("text", label=lm_eqn(lm(Concentration~Amplitude, predict.frame[ vals$keeprows, , drop = FALSE])), x=0, y=Inf, hjust=0, vjust=1, parse=TRUE)+
            geom_point() +
            geom_point(data = predict.frame[!vals$keeprows, , drop = FALSE], shape = 21, fill = "red", color = "black", alpha = 0.25) +
            stat_smooth(method="lm", fullrange = TRUE) +
            scale_x_continuous(paste(element.name, intens)) +
            scale_y_continuous(paste(element.name, conen)) +
            coord_cartesian(xlim = rangescalcurve$x, ylim = rangescalcurve$y, expand = TRUE)
            
        }
        
        if(input$radiocal==2){
            calcurve.plot <- ggplot(data=predict.frame[ vals$keeprows, , drop = FALSE], aes(Amplitude, Concentration)) +
            theme_light() +
            annotate("text", label=lm_eqn_poly(lm(Concentration~Amplitude + I(Amplitude^2), predict.frame[ vals$keeprows, , drop = FALSE])), x=0, y=Inf, hjust=0, vjust=1, parse=TRUE)+
            geom_point() +
            geom_point(data = predict.frame[!vals$keeprows, , drop = FALSE], shape = 21, fill = "red", color = "black", alpha = 0.25) +
            stat_smooth(method="lm", formula=y~poly(x,2)) +
            scale_x_continuous(paste(element.name, intens)) +
            scale_y_continuous(paste(element.name, conen)) +
            coord_cartesian(xlim = rangescalcurve$x, ylim = rangescalcurve$y, expand = TRUE)
            
        }
        
        if(input$radiocal==3){
            calcurve.plot <- ggplot(data=val.frame[ vals$keeprows, , drop = FALSE], aes(AmplitudeNorm, Concentration)) +
            theme_light() +
            annotate("text", label=lm_eqn(lm(Concentration~., val.frame[ vals$keeprows, , drop = FALSE])), x=0, y=Inf, hjust=0, vjust=1, parse=TRUE)+
            geom_point() +
            geom_point(aes(AmplitudeNorm, Concentration), data = val.frame[!vals$keeprows, , drop = FALSE], shape = 21, fill = "red", color = "black", alpha = 0.25) +
            geom_smooth(aes(x=AmplitudeNorm, y=Concentration, ymin = Lower, ymax = Upper)) +
            scale_x_continuous(paste(element.name, norma)) +
            scale_y_continuous(paste(element.name, conen)) +
            coord_cartesian(xlim = rangescalcurve$x, ylim = rangescalcurve$y, expand = TRUE)
            
        }
        
        if(input$radiocal==4){
            calcurve.plot <- ggplot(data=val.frame[ vals$keeprows, , drop = FALSE], aes(AmplitudeNorm, Concentration)) +
            theme_light() +
            annotate("text", label=lm_eqn(lm(Concentration~., val.frame[ vals$keeprows, , drop = FALSE])), x=0, y=Inf, hjust=0, vjust=1, parse=TRUE)+
            geom_point() +
            geom_point(aes(AmplitudeNorm, Concentration), data = val.frame[!vals$keeprows, , drop = FALSE], shape = 21, fill = "red", color = "black", alpha = 0.25) +
            geom_smooth() +
            scale_x_continuous(paste(element.name, norma)) +
            scale_y_continuous(paste(element.name, conen)) +
            coord_cartesian(xlim = rangescalcurve$x, ylim = rangescalcurve$y, expand = TRUE)
            
        }
        
        
        
        calcurve.plot
        
        
    })
    
    observeEvent(input$cropcal, {
        brush <- input$plot_cal_brush
        if (!is.null(brush)) {
            rangescalcurve$x <- c(brush$xmin, brush$xmax)
            rangescalcurve$y <- c(brush$ymin, brush$ymax)
            
        } else {
            rangescalcurve$x <- NULL
            rangescalcurve$y <- NULL
        }
    })
    
    output$calcurveplots <- renderPlot({
        calCurvePlot()
    })
    
    
    rangesvalcurve <- reactiveValues(x = NULL, y = NULL)
    
    
    valCurvePlot <- reactive({
        
        predict.frame <- predictFrame()
        line.model <- lineModel()
        
        
        
        element.name <- gsub("[.]", "", substr(input$calcurveline, 1, 2))
        intens <- " Counts per Second"
        norma <- " Normalized"
        norma.comp <- " Compton Normalized"
        norma.tc <- " Valid Counts Normalized"
        conen <- " (%)"
        predi <- " Estimate (%)"
        log <- "Log "
        
        amplitude.name <- c(element.name, intens)
        concentration.name <- c(element.name, conen)
        prediction.name <- c(element.name, predi)
        val.frame <- valFrame()
        
        
        valcurve.plot <- ggplot(data=val.frame[ vals$keeprows, , drop = FALSE], aes(Prediction, Concentration)) +
        theme_bw() +
        annotate("text", label=lm_eqn_val(lm(Concentration~Prediction, val.frame[ vals$keeprows, , drop = FALSE])), x=0, y=Inf, hjust=0, vjust=1, parse=TRUE)+
        geom_abline(intercept=0, slope=1, lty=2) +
        stat_smooth(method="lm") +
        geom_point() +
        geom_point(aes(Prediction, Concentration),  data = val.frame[!vals$keeprows, , drop = FALSE], shape = 21, fill = "red", color = "black", alpha = 0.25) +
        scale_x_continuous(paste(element.name, predi)) +
        scale_y_continuous(paste(element.name, conen)) +
        coord_cartesian(xlim = rangesvalcurve$x, ylim = rangesvalcurve$y, expand = TRUE)
        
        
        
        
        
        valcurve.plot
        
    })
    
    observeEvent(input$cropval, {
        brush <- input$plot_val_brush
        if (!is.null(brush)) {
            rangesvalcurve$x <- c(brush$xmin, brush$xmax)
            rangesvalcurve$y <- c(brush$ymin, brush$ymax)
            
        } else {
            rangesvalcurve$x <- NULL
            rangesvalcurve$y <- NULL
        }
    })
    
    
    output$valcurveplots <- renderPlot({
        valCurvePlot()
    })
    
    
    
    calPlotDownload <- reactive({
        
        grid.arrange(calCurvePlot(), valCurvePlot(), ncol=2)
        
    })
    
    
    output$downloadcloudplot <- downloadHandler(
    filename = function() { paste(paste(c(input$calname, "_", input$calcurveline), collapse=''), '.tiff',  sep='') },
    content = function(file) {
        ggsave(file,calPlotDownload(), device="tiff", compression="lzw", type="cairo", dpi=300, width=12, height=7)
    }
    )
    
    
    calValTable <- reactive({
        
        standard.table <- valFrame()
        hold.frame <- holdFrame()
        
        standard.table.summary <- data.frame(hold.frame$Spectrum, standard.table$Concentration, standard.table$Prediction, standard.table$Concentration-standard.table$Prediction, ((standard.table$Concentration-standard.table$Prediction)/standard.table$Concentration))
        colnames(standard.table.summary) <- c("Standard", "Concentration", "Prediction", "Difference", "Relative")
        
        standard.table.summary[,-1] <-round(standard.table.summary[,-1],4)
        standard.table.summary[,5] <- as.character(percent(standard.table.summary[,5]))
        
        this.table <- standard.table.summary
        this.table
        
    })
    
    
    output$standardsperformance <- DT::renderDataTable({
        
        
        standard.table <- calValTable()
        standard.table
        
    }, options =list(aoColumnDefs = list(list(sClass="alignRight",aTargets=c(list(2), list(3),list(4),list(5))))  ))
    
    
    randomizeData <- reactive({
        
        cal.frame <- concentrationTable()
        cal.frame <- cal.frame[ vals$keeprows, , drop = FALSE]
        total.number <- length(cal.frame[,1])
        sample.number <- total.number-round(input$percentrandom*total.number, 0)
        
        hold <- cal.frame[sample(nrow(cal.frame), sample.number),]
        cal.frame$Spectrum %in% hold$Spectrum
        
    })
    
    
    
    calCurveFrameRandomized <- reactive({
        
        predict.frame <- predictFrame()
        predict.frame <- predict.frame[ vals$keeprows, , drop = FALSE]
        
        predict.frame[randomizeData(),]
        
    })
    
    
    lineModelRandom <- reactive({
        
        predict.frame <- calCurveFrameRandomized()
        
        
        if (input$radiocal==1){
            cal.lm <- lm(Concentration~Amplitude, data=predict.frame)
        }
        
        
        if (input$radiocal==2){
            cal.lm <- lm(Concentration~Amplitude + I(Amplitude^2), data=predict.frame)
        }
        
        if (input$radiocal==3){
            cal.lm <- lm(Concentration~., data=predict.frame)
        }
        
        if (input$radiocal==4){
            cal.lm <- randomForest(Concentration~., data=predict.frame, na.action=na.omit)
        }
        
        cal.lm
        
    })
    
    
    valFrameRandomized <- reactive({
        
        predict.amplitude <- predictAmplitude()[ vals$keeprows, , drop = FALSE]
        predict.frame <- predictFrame()[ vals$keeprows, , drop = FALSE]
        
        predict.amplitude <- predict.amplitude[!(randomizeData()), , drop = FALSE]
        predict.frame <- predict.frame[!(randomizeData()), , drop = FALSE]
        line.model <- lineModelRandom()
        
        
        
        if (input$radiocal==1){
            cal.est.conc.pred <- predict(object=line.model, newdata=predict.amplitude, interval='confidence')
            cal.est.conc.tab <- data.frame(cal.est.conc.pred)
            cal.est.conc <- cal.est.conc.tab$fit
            
            val.frame <- data.frame(predict.frame$Concentration, cal.est.conc)
            colnames(val.frame) <- c("Concentration", "Prediction")
        }
        
        if (input$radiocal==2){
            cal.est.conc.pred <- predict(object=line.model, newdata=predict.amplitude, interval='confidence')
            cal.est.conc.tab <- data.frame(cal.est.conc.pred)
            cal.est.conc <- cal.est.conc.tab$fit
            
            val.frame <- data.frame(predict.frame$Concentration, cal.est.conc)
            colnames(val.frame) <- c("Concentration", "Prediction")
        }
        
        if (input$radiocal==3){
            
            xmin = 0; xmax=10
            N = length(predict.frame$Concentration)
            means = colMeans(predict.frame)
            dummyDF = t(as.data.frame(means))
            for(i in 2:N){dummyDF=rbind(dummyDF,means)}
            xv=seq(xmin,xmax, length.out=N)
            dummyDF$Concentration = xv
            yv=predict(line.model, newdata=predict.amplitude)
            
            
            lucas.x <- yv
            
            cal.est.conc.pred.luc <- predict(object=line.model , newdata=predict.amplitude, interval='confidence')
            cal.est.conc.tab <- data.frame(cal.est.conc.pred.luc)
            cal.est.conc.luc <- cal.est.conc.tab$fit
            cal.est.conc.luc.up <- cal.est.conc.tab$upr
            cal.est.conc.luc.low <- cal.est.conc.tab$lwr
            
            
            val.frame <- data.frame(predict.frame$Concentration, predict.amplitude$Amplitude, lucas.x, cal.est.conc.luc, cal.est.conc.luc.up, cal.est.conc.luc.low)
            colnames(val.frame) <- c("Concentration", "Amplitude", "AmplitudeNorm", "Prediction", "Upper", "Lower")
        }
        
        if (input$radiocal==4){
            
            xmin = 0; xmax=10
            N = length(predict.frame$Concentration)
            means = colMeans(predict.frame)
            dummyDF = t(as.data.frame(means))
            for(i in 2:N){dummyDF=rbind(dummyDF,means)}
            xv=seq(xmin,xmax, length.out=N)
            dummyDF$Concentration = xv
            yv=predict(line.model, newdata=predict.amplitude)
            
            
            lucas.x <- yv
            
            cal.est.conc.pred.luc <- predict(object=line.model , newdata=predict.amplitude, interval='confidence')
            #cal.est.conc.tab <- data.frame(cal.est.conc.pred.luc)
            #cal.est.conc.luc <- cal.est.conc.tab$fit
            #cal.est.conc.luc.up <- cal.est.conc.tab$upr
            #cal.est.conc.luc.low <- cal.est.conc.tab$lwr
            
            
            val.frame <- data.frame(predict.frame$Concentration, predict.amplitude$Amplitude, lucas.x, as.vector(cal.est.conc.pred.luc))
            colnames(val.frame) <- c("Concentration", "Amplitude", "AmplitudeNorm", "Prediction")
        }
        
        
        
        
        val.frame
        
    })
    
    
    valFrameRandomizedRev <- reactive({
        
        predict.amplitude <- predictAmplitude()[ vals$keeprows, , drop = FALSE]
        predict.frame <- predictFrame()[ vals$keeprows, , drop = FALSE]
        
        predict.amplitude <- predict.amplitude[(randomizeData()), ]
        predict.frame <- predict.frame[(randomizeData()), ]
        line.model <- lineModelRandom()
        
        
        
        if (input$radiocal==1){
            cal.est.conc.pred <- predict(object=line.model, newdata=predict.amplitude, interval='confidence')
            cal.est.conc.tab <- data.frame(cal.est.conc.pred)
            cal.est.conc <- cal.est.conc.tab$fit
            
            val.frame <- data.frame(predict.frame$Concentration, cal.est.conc)
            colnames(val.frame) <- c("Concentration", "Prediction")
        }
        
        if (input$radiocal==2){
            cal.est.conc.pred <- predict(object=line.model, newdata=predict.amplitude, interval='confidence')
            cal.est.conc.tab <- data.frame(cal.est.conc.pred)
            cal.est.conc <- cal.est.conc.tab$fit
            
            val.frame <- data.frame(predict.frame$Concentration, cal.est.conc)
            colnames(val.frame) <- c("Concentration", "Prediction")
        }
        
        if (input$radiocal==3){
            
            xmin = 0; xmax=10
            N = length(predict.frame$Concentration)
            means = colMeans(predict.frame)
            dummyDF = t(as.data.frame(means))
            for(i in 2:N){dummyDF=rbind(dummyDF,means)}
            xv=seq(xmin,xmax, length.out=N)
            dummyDF$Concentration = xv
            yv=predict(line.model, newdata=predict.amplitude)
            
            
            lucas.x <- yv
            
            cal.est.conc.pred.luc <- predict(object=line.model , newdata=predict.amplitude, interval='confidence')
            cal.est.conc.tab <- data.frame(cal.est.conc.pred.luc)
            cal.est.conc.luc <- cal.est.conc.tab$fit
            cal.est.conc.luc.up <- cal.est.conc.tab$upr
            cal.est.conc.luc.low <- cal.est.conc.tab$lwr
            
            
            val.frame <- data.frame(predict.frame$Concentration, predict.amplitude$Amplitude, lucas.x, cal.est.conc.luc, cal.est.conc.luc.up, cal.est.conc.luc.low)
            colnames(val.frame) <- c("Concentration", "Amplitude", "AmplitudeNorm", "Prediction", "Upper", "Lower")
        }
        
        if (input$radiocal==4){
            
            xmin = 0; xmax=10
            N = length(predict.frame$Concentration)
            means = colMeans(predict.frame)
            dummyDF = t(as.data.frame(means))
            for(i in 2:N){dummyDF=rbind(dummyDF,means)}
            xv=seq(xmin,xmax, length.out=N)
            dummyDF$Concentration = xv
            yv=predict(line.model, newdata=predict.amplitude)
            
            
            lucas.x <- yv
            
            cal.est.conc.pred.luc <- predict(object=line.model , newdata=predict.amplitude, interval='confidence')
            
            
            
            val.frame <- data.frame(predict.frame$Concentration, predict.amplitude$Amplitude, lucas.x, as.vector(cal.est.conc.pred.luc))
            colnames(val.frame) <- c("Concentration", "Amplitude", "AmplitudeNorm", "Prediction")
        }
        
        
        
        
        val.frame
        
    })
    
    rangescalcurverandom <- reactiveValues(x = NULL, y = NULL)
    
    
    calCurvePlotRandom <- reactive({
        
        predict.frame <- calCurveFrameRandomized()
        line.model <- lineModelRandom()
        
        
        element.name <- gsub("[.]", "", substr(input$calcurveline, 1, 2))
        intens <- " Counts per Second"
        norma <- " Normalized"
        norma.comp <- " Compton Normalized"
        norma.tc <- " Valid Counts Normalized"
        conen <- " (%)"
        predi <- " Estimate (%)"
        
        amplitude.name <- c(element.name, intens)
        concentration.name <- c(element.name, conen)
        prediction.name <- c(element.name, predi)
        
        
        if(input$radiocal==1){
            calcurve.plot <- ggplot(data=predict.frame, aes(Amplitude, Concentration)) +
            theme_light() +
            annotate("text", label=lm_eqn(line.model), x=0, y=Inf, hjust=0, vjust=1, parse=TRUE)+
            geom_point() +
            geom_point(data = predict.frame, shape = 21, fill = "red", color = "black", alpha = 0.25) +
            stat_smooth(method="lm", fullrange = TRUE) +
            scale_x_continuous(paste(element.name, intens)) +
            scale_y_continuous(paste(element.name, conen)) +
            coord_cartesian(xlim = rangescalcurverandom$x, ylim = rangescalcurverandom$y, expand = TRUE)
            
        }
        
        if(input$radiocal==2){
            
            calcurve.plot <- ggplot(data=predict.frame, aes(Amplitude, Concentration)) +
            theme_light() +
            annotate("text", label=lm_eqn_poly(line.model), x=0, y=Inf, hjust=0, vjust=1, parse=TRUE)+
            geom_point() +
            geom_point(data = predict.frame, shape = 21, fill = "red", color = "black", alpha = 0.25) +
            stat_smooth(method="lm", formula=y~poly(x,2)) +
            scale_x_continuous(paste(element.name, intens)) +
            scale_y_continuous(paste(element.name, conen)) +
            coord_cartesian(xlim = rangescalcurverandom$x, ylim = rangescalcurverandom$y, expand = TRUE)
        }
        
        if(input$radiocal==3){
            val.frame <- valFrameRandomizedRev()
            
            calcurve.plot <- ggplot(data=val.frame, aes(AmplitudeNorm, Concentration)) +
            theme_light() +
            annotate("text", label=lm_eqn(line.model), x=0, y=Inf, hjust=0, vjust=1, parse=TRUE)+
            geom_point() +
            geom_point(aes(AmplitudeNorm, Concentration), data = val.frame, shape = 21, fill = "red", color = "black", alpha = 0.25) +
            geom_smooth(aes(x=AmplitudeNorm, y=Concentration, ymin = Lower, ymax = Upper)) +
            scale_x_continuous(paste(element.name, norma)) +
            scale_y_continuous(paste(element.name, conen)) +
            coord_cartesian(xlim = rangescalcurverandom$x, ylim = rangescalcurverandom$y, expand = TRUE)
        }
        
        if(input$radiocal==4){
            val.frame <- valFrameRandomizedRev()
            
            calcurve.plot <- ggplot(data=val.frame, aes(AmplitudeNorm, Concentration)) +
            theme_light() +
            annotate("text", label=lm_eqn(lm(Concentration~., val.frame)), x=0, y=Inf, hjust=0, vjust=1, parse=TRUE)+
            geom_point() +
            geom_point(aes(AmplitudeNorm, Concentration), data = val.frame, shape = 21, fill = "red", color = "black", alpha = 0.25) +
            geom_smooth() +
            scale_x_continuous(paste(element.name, norma)) +
            scale_y_continuous(paste(element.name, conen)) +
            coord_cartesian(xlim = rangescalcurverandom$x, ylim = rangescalcurverandom$y, expand = TRUE)
        }
        
        calcurve.plot
        
        
    })
    
    observeEvent(input$cropcalrandom, {
        brush <- input$plot_cal_brush_random
        if (!is.null(brush)) {
            rangescalcurverandom$x <- c(brush$xmin, brush$xmax)
            rangescalcurverandom$y <- c(brush$ymin, brush$ymax)
            
        } else {
            rangescalcurverandom$x <- NULL
            rangescalcurverandom$y <- NULL
        }
    })
    
    
    
    output$calcurveplotsrandom <- renderPlot({
        calCurvePlotRandom()
    })
    
    
    rangesvalcurverandom <- reactiveValues(x = NULL, y = NULL)
    
    valCurvePlotRandom <- reactive({
        
        
        
        element.name <- gsub("[.]", "", substr(input$calcurveline, 1, 2))
        intens <- " Counts per Second"
        norma <- " Normalized"
        norma.comp <- " Compton Normalized"
        norma.tc <- " Valid Counts Normalized"
        conen <- " (%)"
        predi <- " Estimate (%)"
        
        amplitude.name <- c(element.name, intens)
        concentration.name <- c(element.name, conen)
        prediction.name <- c(element.name, predi)
        
        val.frame <- valFrameRandomized()
        
        valcurve.plot <- ggplot(data=val.frame, aes(Prediction, Concentration)) +
        theme_bw() +
        annotate("text", label=lm_eqn_val(lm(Concentration~Prediction, val.frame)), x=0, y=Inf, hjust=0, vjust=1, parse=TRUE)+
        geom_abline(intercept=0, slope=1, lty=2) +
        stat_smooth(method="lm") +
        geom_point() +
        geom_point(aes(Prediction, Concentration),  data = val.frame, shape = 21, fill = "red", color = "black", alpha = 0.25) +
        scale_x_continuous(paste(element.name, predi)) +
        scale_y_continuous(paste(element.name, conen)) +
        coord_cartesian(xlim = rangesvalcurverandom$x, ylim = rangesvalcurverandom$y, expand = TRUE)
        
        valcurve.plot
        
    })
    
    observeEvent(input$cropvalrandom, {
        brush <- input$plot_val_brush_random
        if (!is.null(brush)) {
            rangesvalcurverandom$x <- c(brush$xmin, brush$xmax)
            rangesvalcurverandom$y <- c(brush$ymin, brush$ymax)
            
        } else {
            rangesvalcurverandom$x <- NULL
            rangesvalcurverandom$y <- NULL
        }
    })
    
    output$valcurveplotsrandom <- renderPlot({
        valCurvePlotRandom()
    })
    
    
    ####CalCurves
    
    # Float over info
    output$hover_infocal <- renderUI({
        
        point.table <- if(calType()!=3){
            calCurveFrame()
        } else if(calType()==3) {
            calValFrame()
        }
        
        concentration.table <- concentrationTable()
        hold.table <- concentration.table[,c("Spectrum", input$calcurveline)]
        colnames(hold.table) <- c("Spectrum", "Selection")
        hold.table$Selection[hold.table$Selection==""] <- NA
        hold.table <- hold.table[complete.cases(hold.table), ]
        
        point.table$Spectrum <- hold.table["Spectrum"]
        
        hover <- input$plot_hovercal
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
        p(HTML(paste0(point$Spectrum
        
        )))
        )
    })
    
    
    # Float over info
    output$hover_infocal_random <- renderUI({
        
        point.table <- if(calType()!=3){
            calCurveFrame()
        } else if(calType()==3) {
            valFrame()
        }
        
        randomized <- randomizeData()
        
        
        point.table <- point.table[ vals$keeprows, , drop = FALSE]
        point.table <- point.table[randomized,]
        
        
        concentration.table <- concentrationTable()
        
        concentration.table <- concentration.table[ vals$keeprows, , drop = FALSE]
        concentration.table <- concentration.table[randomized,]
        
        hold.table <- concentration.table[,c("Spectrum", input$calcurveline)]
        colnames(hold.table) <- c("Spectrum", "Selection")
        
        
        point.table$Spectrum <- hold.table["Spectrum"]
        
        hover <- input$plot_hovercal_random
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
        p(HTML(paste0(point$Spectrum
        
        )))
        )
    })
    
    
    # Toggle points that are clicked
    observeEvent(input$plot_cal_click, {
        
        predict.frame <- if(calType()!=3){
            calCurveFrame()
        } else if(calType()==3) {
            calValFrame()
        }
        
        res <- nearPoints(predict.frame, input$plot_cal_click, allRows = TRUE)
        
        vals$keeprows <- xor(vals$keeprows, res$selected_)
    })
    
    
    # Toggle points that are brushed, when button is clicked
    observeEvent(input$exclude_toggle, {
        
        predict.frame <- if(calType()!=3){
            calCurveFrame()
        } else if(calType()==3) {
            calValFrame()
        }
        res <- brushedPoints(predict.frame, input$plot_cal_brush, allRows = TRUE)
        
        vals$keeprows <- xor(vals$keeprows, res$selected_)
    })
    
    # Reset all points
    observeEvent(input$exclude_reset, {
        
        predict.frame <- if(calType()!=3){
            calCurveFrame()
        } else if(calType()==3) {
            calValFrame()
        }
        vals$keeprows <- rep(TRUE, nrow(predict.frame))
    })
    
    # Reset all points on element change
    observeEvent(input$calcurveline, {
        
        
        
        vals$keeprows <- calFileStandards()
        
        
        
    })
    
    
    
    
    
    ####ValCurves
    
    
    
    
    # Float over info
    output$hover_infoval <- renderUI({
        
        point.table <- calValFrame()
        concentration.table <- concentrationTable()
        hold.table <- concentration.table[,c("Spectrum", input$calcurveline)]
        colnames(hold.table) <- c("Spectrum", "Selection")
        hold.table$Selection[hold.table$Selection==""] <- NA
        hold.table <- hold.table[complete.cases(hold.table), ]
        
        point.table$Spectrum <- hold.table["Spectrum"]
        
        
        hover <- input$plot_hoverval
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
        p(HTML(paste0(point$Spectrum
        
        )))
        )
    })
    
    
    output$hover_infoval_random <- renderUI({
        
        point.table <- calValFrame()
        
        randomized <- randomizeData()
        
        
        point.table <- point.table[ vals$keeprows, , drop = FALSE]
        point.table <- point.table[!(randomized),]
        
        concentration.table <- concentrationTable()
        
        concentration.table <- concentration.table[ vals$keeprows, , drop = FALSE]
        concentration.table <- concentration.table[!(randomized),]
        concentration.table.rev <- concentration.table[(randomized),]
        
        hold.table <- concentration.table[,c("Spectrum", input$calcurveline)]
        colnames(hold.table) <- c("Spectrum", "Selection")
        
        
        point.table$Spectrum <- hold.table["Spectrum"]
        
        
        point.table <- point.table[point.table$Concentration > min(concentration.table.rev[,input$calcurveline], na.rm = TRUE) & point.table$Concentration < max(concentration.table.rev[,input$calcurveline], na.rm = TRUE), ]
        
        
        
        hover <- input$plot_hoverval_random
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
        p(HTML(paste0(point$Spectrum
        
        )))
        )
    })
    
    # Toggle points that are clicked
    observeEvent(input$plot_val_click, {
        
        predict.frame <- calValFrame()
        
        res <- nearPoints(predict.frame, input$plot_val_click, allRows = TRUE)
        
        vals$keeprows <- xor(vals$keeprows, res$selected_)
    })
    
    
    # Toggle points that are brushed, when button is clicked
    observeEvent(input$exclude_toggle, {
        predict.frame <- calValFrame()
        
        res <- brushedPoints(predict.frame, input$plot_val_brush, allRows = TRUE)
        
        vals$keeprows <- xor(vals$keeprows, res$selected_)
    })
    
    # Reset all points
    observeEvent(input$exclude_reset, {
        predict.frame <- calValFrame()
        
        vals$keeprows <- rep(TRUE, nrow(predict.frame))
    })
    
    
    normalLM <- reactive({
        
        
        model <- lineModel()
        
        model.frame <- as.data.frame(augment(model))
        
        model.frame$qq <- qqnorm(model.frame$.std.resid)[[1]]
        
        model.frame$sqrt.std.resid <- sqrt(abs(model.frame$.std.resid))
        
        model.frame$seq.cooksd <- seq_along(model.frame$.cooksd)
        
        #model.frame$Spectrum <- predictFrameName()$Spectrum
        
        
        
        model.frame
        
    })
    
    
    forestLM <- reactive({
        
        
        model <- lm(Concentration~Prediction, data=as.data.frame(calValTable()))
        
        model.frame <- as.data.frame(augment(model))
        
        model.frame$qq <- qqnorm(model.frame$.std.resid)[[1]]
        
        model.frame$sqrt.std.resid <- sqrt(abs(model.frame$.std.resid))
        
        model.frame$seq.cooksd <- seq_along(model.frame$.cooksd)
        
        #model.frame$Spectrum <- predictFrameName()$Spectrum
        
        
        
        model.frame
        
        
    })
    
    
    modelFrame <- reactive({
        
        if(input$radiocal!=4){
            normalLM()
        } else if(input$radiocal==4){
            forestLM()
        }
        
        
        
    })
    
    
    
    
    
    diagResidualsFitted <- reactive({
        
        model <- modelFrame()
        
        p1 <- ggplot(model[ vals$keeprows, , drop = FALSE], aes(.fitted, .resid)) +
        stat_smooth(method="loess") +
        geom_hline(yintercept=0, col="red", linetype="dashed") +
        xlab("Fitted values") +
        ylab("Residuals") +
        ggtitle("Residual vs Fitted Plot") +
        theme_light() +
        geom_point() +
        geom_point(data=model[ !vals$keeprows, , drop = FALSE], aes(.fitted, .resid), shape = 21, fill = "red", color = "black", alpha = 0.25)
        
        p1
        
    })
    
    output$residualsfitted <- renderPlot({
        
        diagResidualsFitted()
        
    })
    
    
    diagQQ <- reactive({
        
        model <- modelFrame()
        
        p2 <- ggplot(model[ vals$keeprows, , drop = FALSE], aes(qq, .std.resid))+geom_point(na.rm = TRUE) +
        geom_abline() +
        xlab("Theoretical Quantiles") +
        ylab("Standardized Residuals") +
        ggtitle("Normal Q-Q") +
        theme_light() +
        geom_point(data=model[ !vals$keeprows, , drop = FALSE], aes(qq, .std.resid), shape = 21, fill = "red", color = "black", alpha = 0.25)
        
        
        p2
        
    })
    
    
    output$qq <- renderPlot({
        
        diagQQ()
        
    })
    
    diagScaleLocation <- reactive({
        
        model <- modelFrame()
        
        
        p3 <- ggplot(model[ vals$keeprows, , drop = FALSE], aes(.fitted, sqrt.std.resid)) +
        stat_smooth(method="loess", na.rm = TRUE) +
        xlab("Fitted Value") +
        ylab(expression(sqrt("|Standardized residuals|"))) +
        ggtitle("Scale-Location") +
        theme_light() +
        geom_point(na.rm=TRUE) +
        geom_point(data=model[ !vals$keeprows, , drop = FALSE], aes(.fitted, sqrt.std.resid), shape = 21, fill = "red", color = "black", alpha = 0.25)
        
        
        p3
        
        
    })
    
    
    output$scalelocation <- renderPlot({
        
        diagScaleLocation()
        
    })
    
    
    diagCooksDistance <- reactive({
        
        model <- modelFrame()
        
        p4 <- ggplot(model, aes(seq.cooksd, .cooksd)) +
        geom_bar(stat="identity", position="identity") +
        xlab("Obs. Number") +
        ylab("Cook's distance") +
        ggtitle("Cook's distance") +
        theme_light()
        
        p4
        
    })
    
    output$cooksdistance <- renderPlot({
        
        diagCooksDistance()
        
    })
    
    
    diagResidualLeverage <- reactive({
        
        model <- modelFrame()
        
        
        p5<-ggplot(model[ vals$keeprows, , drop = FALSE], aes(.hat, .std.resid))+
        geom_point(aes(size=.cooksd), na.rm=TRUE) +
        geom_point(data=model[ !vals$keeprows, , drop = FALSE], aes(.hat, .std.resid), shape = 21, fill = "red", color = "black", alpha = 0.25) +
        stat_smooth(method="loess", na.rm=TRUE) +
        xlab("Leverage") +
        ylab("Standardized Residuals") +
        ggtitle("Residual vs Leverage Plot") +
        scale_size_continuous("Cook's Distance", range=c(1,5)) +
        theme_light() +
        theme(legend.position="bottom")
        
        p5
        
    })
    
    output$residualleverage <- renderPlot({
        
        diagResidualLeverage()
        
    })
    
    
    diagCooksLeverage <- reactive({
        
        model <- modelFrame()
        
        
        p6 <- ggplot(model[ vals$keeprows, , drop = FALSE], aes(.hat, .cooksd)) +
        stat_smooth(method="loess", na.rm=TRUE) +
        xlab("Leverage hii") +
        ylab("Cook's Distance") +
        ggtitle("Cook's dist vs Leverage hii/(1-hii)") +
        geom_abline(slope=seq(0,3,0.5), color="gray", linetype="dashed") +
        theme_light() +
        geom_point(na.rm=TRUE) +
        geom_point(data=model[ vals$keeprows, , drop = FALSE], aes(.hat, .cooksd), shape = 21, fill = "red", color = "black", alpha = 0.25)
        
        p6
        
    })
    
    
    output$cooksleverage <- renderPlot({
        
        diagCooksLeverage()
        
    })
    
    
    diagPlotDownload <- reactive({
        
        grid.arrange(diagResidualsFitted(), diagQQ(),
        diagScaleLocation(), diagCooksDistance(),
        diagResidualLeverage(), diagCooksLeverage(),
        ncol=2, nrow=3)
        
    })
    
    #########Diagnostic Plot Controls#######
    ####Residuals Fitted
    # Float over info
    output$hover_inforesidualsfitted <- renderUI({
        
        point.table <- modelFrame()
        
        
        hover <- input$plot_hoverresidualsfitted
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
        p(HTML(paste0(point$Spectrum
        
        )))
        )
    })
    
    # Toggle points that are clicked
    observeEvent(input$plot_residualsfitted_click, {
        
        predict.frame <- modelFrame()
        
        res <- nearPoints(predict.frame, input$plot_residualsfitted_click, allRows = TRUE)
        
        vals$keeprows <- xor(vals$keeprows, res$selected_)
    })
    
    
    # Toggle points that are brushed, when button is clicked
    observeEvent(input$exclude_toggle_diag, {
        
        predict.frame <- modelFrame()
        
        res <- brushedPoints(predict.frame, input$plot_residualsfitted_brush, allRows = TRUE)
        
        vals$keeprows <- xor(vals$keeprows, res$selected_)
    })
    
    
    ####QQ Norm
    # Float over info
    output$hover_infoqq <- renderUI({
        
        point.table <- modelFrame()
        
        
        hover <- input$plot_hoverqq
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
        p(HTML(paste0(point$Spectrum
        
        )))
        )
    })
    
    # Toggle points that are clicked
    observeEvent(input$plot_qq_click, {
        
        predict.frame <- modelFrame()
        
        res <- nearPoints(predict.frame, input$plot_qq_click, allRows = TRUE)
        
        vals$keeprows <- xor(vals$keeprows, res$selected_)
    })
    
    
    # Toggle points that are brushed, when button is clicked
    observeEvent(input$exclude_toggle_diag, {
        
        predict.frame <- modelFrame()
        
        res <- brushedPoints(predict.frame, input$plot_qq_brush, allRows = TRUE)
        
        vals$keeprows <- xor(vals$keeprows, res$selected_)
    })
    
    
    ####Scaled Residuals
    # Float over info
    output$hover_infoscalelocation <- renderUI({
        
        point.table <- modelFrame()
        
        
        hover <- input$plot_hoverscalelocation
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
        p(HTML(paste0(point$Spectrum
        
        )))
        )
    })
    
    # Toggle points that are clicked
    observeEvent(input$plot_scalelocation_click, {
        
        predict.frame <- modelFrame()
        
        res <- nearPoints(predict.frame, input$plot_scalelocation_click, allRows = TRUE)
        
        vals$keeprows <- xor(vals$keeprows, res$selected_)
    })
    
    
    # Toggle points that are brushed, when button is clicked
    observeEvent(input$exclude_toggle_diag, {
        
        predict.frame <- modelFrame()
        
        res <- brushedPoints(predict.frame, input$plot_scalelocation_brush, allRows = TRUE)
        
        vals$keeprows <- xor(vals$keeprows, res$selected_)
    })
    
    
    ####Residuals Leverage
    # Float over info
    output$hover_inforesidualleverage <- renderUI({
        
        point.table <- modelFrame()
        
        
        hover <- input$plot_hoverresidualleverage
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
        p(HTML(paste0(point$Spectrum
        
        )))
        )
    })
    
    # Toggle points that are clicked
    observeEvent(input$plot_residualleverage_click, {
        
        predict.frame <- modelFrame()
        
        res <- nearPoints(predict.frame, input$plot_residualleverage_click, allRows = TRUE)
        
        vals$keeprows <- xor(vals$keeprows, res$selected_)
    })
    
    
    # Toggle points that are brushed, when button is clicked
    observeEvent(input$exclude_toggle_diag, {
        
        predict.frame <- modelFrame()
        
        res <- brushedPoints(predict.frame, input$plot_residualleverage_brush, allRows = TRUE)
        
        vals$keeprows <- xor(vals$keeprows, res$selected_)
    })
    
    
    ####Cooks Leverage
    # Float over info
    output$hover_infocooksleverage <- renderUI({
        
        point.table <- modelFrame()
        
        
        hover <- input$plot_hovercooksleverage
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
        p(HTML(paste0(point$Spectrum
        
        )))
        )
    })
    
    # Toggle points that are clicked
    observeEvent(input$plot_cooksleverage_click, {
        
        predict.frame <- modelFrame()
        
        res <- nearPoints(predict.frame, input$plot_cooksleverage_click, allRows = TRUE)
        
        vals$keeprows <- xor(vals$keeprows, res$selected_)
    })
    
    
    # Toggle points that are brushed, when button is clicked
    observeEvent(input$exclude_toggle_diag, {
        
        predict.frame <- modelFrame()
        
        res <- brushedPoints(predict.frame, input$plot_cooksleverage_brush, allRows = TRUE)
        
        vals$keeprows <- xor(vals$keeprows, res$selected_)
    })
    
    
    
    # Reset all points
    observeEvent(input$exclude_reset_diag, {
        
        predict.frame <- modelFrame()
        
        vals$keeprows <- rep(TRUE, nrow(predict.frame))
    })
    
    
    
    
    #output$downloadcal <- downloadHandler(
    #filename = function() { paste(input$dataset, '.csv', sep=',') },
    #content = function(file
    #) {
    #  write.csv(metadataForm(), file)
    # }
    #)
    
    
    
    ####alternative take
    
    nullList <- reactive({
        
        spectra.line.table <- waveSubset()
        
        cal.vector <- Wavenumberlinestouse()
        cal.vector2 <- cal.vector[2:length(cal.vector)]
        cal.list <- as.list(cal.vector2)
        setNames(cal.list, cal.vector2)
        cal.list <- pblapply(cal.list, function(x) return(NULL))
        nullList <- cal.list
        
        
    })
    
    
    
    
    #rf2 <- reactiveValues()
    #observe({
    #    if(input$createcalelement > 0){
    #    calList[[input$calcurveline]] <- lineModel()
    #    }
    #    rf2 <<- calList
    #})
    
    emptyList <- reactive({
        a.list <- list()
        a.list
    })
    
    
    
    
    
    calList <- reactiveValues()
    
    observeEvent(!is.null(input$file1), {
        isolate(calList <- emptyList())
        calList <<- calList
    })
    
    
    
    
    observeEvent(input$createcalelement, {
        
        cal.condition <- input$radiocal
        norm.condition <- input$normcal
        
        norm.min <- input$comptonmin
        norm.max <- input$comptonmax
        
        cal.table <- data.frame(cal.condition, norm.condition, norm.min, norm.max)
        colnames(cal.table) <- c("CalType", "NormType", "Min", "Max")
        
        slope.corrections <- input$slope_vars
        intercept.corrections <- input$intercept_vars
        
        standards.used <- vals$keeprows
        
        cal.mode.list <- list(cal.table, slope.corrections, intercept.corrections, standards.used)
        names(cal.mode.list) <- c("CalTable", "Slope", "Intercept", "StandardsUsed")
        
        calConditons <<- cal.mode.list
        
    })
    
    
    observeEvent(input$createcalelement, {
        
        calList[[input$calcurveline]] <- list(isolate(calConditons), isolate(strip_glm(lineModel())))
        
        calList <<- calList
        
    })
    
    
    
    calPlotList <- reactiveValues()
    calPlotList <- emptyList()
    observeEvent(input$createcalelement, {
        
        
        calPlotList[[input$calcurveline]] <- isolate(calPlotDownload())
        
        calPlotList <<- calPlotList
        
    })
    
    diagPlotList <- reactiveValues()
    diagPlotList <- emptyList()
    #observeEvent(input$createcalelement, {
    
    
    #diagPlotList[[input$calcurveline]] <- isolate(diagPlotDownload())
    
    #diagPlotList <<- diagPlotList
    
    #})
    
    Calibration <- reactiveValues()
    observeEvent(input$createcal, {
        
        
        spectra.line.table <- if(dataType()=="Spectra"){
            waveSubset()
        } else if(dataType()=="Net"){
            waveSubset()
        }
        cal.intensities <- spectra.line.table
        cal.values <- values[["DF"]]
        cal.data <- if(dataType()=="Spectra"){
            dataManipulate()
        } else if(dataType()=="Net"){
            dataManipulate()
        }
        
        cal.defs <- wavevalues[["DF"]]
        
        
        calibrationList <- NULL
        calibrationList <- list(input$filetype, input$calunits, cal.data, cal.intensities, cal.values, cal.defs, calList)
        names(calibrationList) <- c("FileType", "Units", "Spectra", "Intensities", "Values", "Definitions", "calList")
        
        Calibration <<- calibrationList
        
        
    })
    
    CalibrationPlots <- reactiveValues()
    observeEvent(input$createcal, {
        
        CalibrationPlots$calCurves <<- calPlotList
        
        
    })
    
    #observeEvent(input$createcal, {
    
    #CalibrationPlots$diagPlots <<- diagPlotList
    
    
    #})
    
    
    
    output$downloadModel <- downloadHandler(
    filename <- function(){
        paste(input$projectname, "quant", sep=".")
    },
    
    content = function(file) {
        saveRDS(Calibration, file = file, compress="xz")
    }
    )
    
    
    output$downloadReport <- downloadHandler(
    function() { paste(paste(c(input$projectname), collapse=''), '.pdf',  sep='') },
    content = function(file){
        ml = marrangeGrob(grobs=CalibrationPlots$calCurves, nrow=1, ncol=1)
        ggsave(file, ml, device="pdf", dpi=300, width=12, height=7)
        
        dev.off()
    })
    
    
    })
    
    
    output$filevalgrab <- renderUI({
        
        if(input$filetype=="CSV") {
            fileInput('loadvaldata', 'Choose CSV', multiple=TRUE,
            accept=c(".csv"))
        } else if(input$filetype=="DPT") {
            fileInput('loadvaldata', 'Choose DPT', multiple=TRUE,
            accept=c(".dpt"))
        } else if(input$filetype=="Opus") {
            fileInput('loadvaldata', 'Choose Opus File', multiple=TRUE,
            accept=NULL)
        }
        
    })
    
    
    
    observeEvent(input$processvalspectra, {
        
        readValDPT <- reactive({
            
            inFile <- input$loadvaldata
            if (is.null(inFile)) return(NULL)
            
            n <- length(inFile$datapath)
            names <- inFile$name
            
            myfiles.frame <- as.data.frame(do.call(rbind, pblapply(seq(1, n, 1), function(x) readDPTData(filepath=inFile$datapath[x], filename=inFile$name[x]))))
            
            myfiles.frame$Wavenumber <- myfiles.frame$Wavenumber + gainshiftHold()
            
            myfiles.frame
            
        })
        
        
        readValCSV <- reactive({
            
            inFile <- input$loadvaldata
            if (is.null(inFile)) return(NULL)
            
            n <- length(inFile$datapath)
            names <- inFile$name
            
            myfiles.frame <- as.data.frame(do.call(rbind, pblapply(seq(1, n, 1), function(x) readCSVData(filepath=inFile$datapath[x], filename=inFile$name[x]))))
            
            myfiles.frame$Wavenumber <- myfiles.frame$Wavenumber + gainshiftHold()
            
            myfiles.frame
            
        })
        
        
        readValOpus <- reactive({
            
            inFile <- input$loadvaldata
            if (is.null(inFile)) return(NULL)
            
            n <- length(inFile$datapath)
            names <- inFile$name
            
            myfiles.frame <- as.data.frame(do.call(rbind, pblapply(seq(1, n, 1), function(x) readOpusData(filepath=inFile$datapath[x], filename=inFile$name[x]))))
            
            myfiles.frame$Wavenumber <- myfiles.frame$Wavenumber + gainshiftHold()
            
            myfiles.frame
            
            
        })
        
        
        
        
        
        myValData <- reactive({
            
            data <- if(input$valfiletype=="Opus"){
                readValOpus()
            } else if(input$valfiletype=="DPT"){
                readValDPT()
            } else if(input$valfiletype=="CSV"){
                readValCSV()
            } 
            
            data
            
        })
        
        
        
        calFileContents2 <- reactive({
            
            existingCalFile <- input$calfileinput2
            
            if (is.null(existingCalFile)) return(NULL)
            
            
            Calibration <- readRDS(existingCalFile$datapath)
            
            Calibration
            
        })
        
        valdata <- myValData()
        
        
        
        output$contents2 <- renderTable({
            
            
            
            myValData()
            
        })
        
        calDefinitions <- reactive({
            
            calFileContents2()[["Definitions"]]
            
        })
        
        calValHold <- reactive({
            
            
            calFileContents2()[["calList"]]
 
            
        })
        
        calVariables <- reactive({
            
            
            calFileContents2()$Intensities
            
            
            
        })
        
        calValElements <- reactive({
            
            as.vector(calDefinitions()$Name)
          
        })
        
        calVariableElements <- reactive({
            variables <- calVariables()
            variableelements <- ls(variables)
            
            #variableelements.simp <- gsub(".K.alpha", "", variableelements)
            #variableelements.simp <- gsub(".K.beta", "", variableelements)
            #variableelements.simp <- gsub(".L.alpha", "", variableelements)
            #variableelements.simp <- gsub(".L.beta", "", variableelements)
            #variableelements.simp <- gsub(".M.line", "", variableelements)
            
            #variableelements <- as.vector(as.character(na.omit(variableelements[match(as.character(fluorescence.lines$Symbol), variableelements.simp)])))
            
            variableelements
        })
        
        valDataType <- reactive({
            
            if(input$valfiletype=="Opus"){
                "Spectra"
            } else if(input$valfiletype=="DPT"){
                "Spectra"
            } else if(input$valfiletype=="CSV"){
                "Spectra"
            }
            
        })
        
        
        
        
        
        
        tableInputValCounts <- reactive({
            
            ftir_parse(range.table = calDefinitions(), data=myValData())
            
        })
        
        
        
        
        fullInputValCounts <- reactive({
            
            tableInputValCounts()
            
        })
        
        
        
        output$myvaltable1 <- renderDataTable({
            
            fullInputValCounts()
            
        })
        
        
        
        
        tableInputValQuant <- reactive({
            
            
            
            count.table <- data.frame(fullInputValCounts())
            the.cal <- calValHold()
            elements.cal <- calValElements()
            elements <- elements.cal[!is.na(match(elements.cal, ls(count.table)))]
            variables <- calVariableElements()
            valdata <- myValData()
            
            #elements <- fluorescence.lines$Symbol[sort(order(fluorescence.lines$Symbol)[elements])]
            
            cal_type <- function(element){
                
                
                if(the.cal[[element]][[1]]$CalTable$CalType==1){
                    1
                } else if(the.cal[[element]][[1]]$CalTable$CalType==2){
                    1
                } else if(the.cal[[element]][[1]]$CalTable$CalType==3){
                    3
                } else if(the.cal[[element]][[1]]$CalTable$CalType==4){
                    4
                }
                
            }
            
            
            
            predicted.list <- lapply(elements, function (x)
            if(valDataType()=="Spectra" && cal_type(x)==1 && the.cal[[x]][[1]]$CalTable$NormType==1){
                predict(
                object=the.cal[[x]][[2]],
                newdata=general.prep(
                spectra.line.table=as.data.frame(
                count.table
                ),
                element.line=x)
                )
            } else if(valDataType()=="Spectra" && cal_type(x)==1 && the.cal[[x]][[1]]$CalTable$NormType==2) {
                predict(
                object=the.cal[[x]][[2]],
                newdata=simple.tc.prep(
                data=valdata,
                spectra.line.table=as.data.frame(
                count.table
                ),
                element.line=x
                )
                )
            } else if(valDataType()=="Spectra" && cal_type(x)==1 && the.cal[[x]][[1]]$CalTable$NormType==3) {
                predict(
                object=the.cal[[x]][[2]],
                newdata=simple.comp.prep(
                data=valdata,
                spectra.line.table=as.data.frame(
                count.table
                ),
                element.line=x,
                norm.min=the.cal[[x]][[1]][1]$CalTable$Min,
                norm.max=the.cal[[x]][[1]][1]$CalTable$Max
                )
                )
            } else if(valDataType()=="Spectra" && cal_type(x)==3 && the.cal[[x]][[1]]$CalTable$NormType==1){
                predict(
                object=the.cal[[x]][[2]],
                newdata=lucas.simp.prep(
                spectra.line.table=as.data.frame(
                count.table
                ),
                element.line=x,
                slope.element.lines=the.cal[[x]][[1]][2]$Slope,
                intercept.element.lines=the.cal[[x]][[1]][3]$Intercept
                )
                )
            } else if(valDataType()=="Spectra" && cal_type(x)==3 && the.cal[[x]][[1]]$CalTable$NormType==2){
                predict(
                object=the.cal[[x]][[2]],
                newdata=lucas.tc.prep(
                data=valdata,
                spectra.line.table=as.data.frame(
                count.table
                ),
                element.line=x,
                slope.element.lines=the.cal[[x]][[1]][2]$Slope,
                intercept.element.lines=the.cal[[x]][[1]][3]$Intercept
                )
                )
            } else if(valDataType()=="Spectra" && cal_type(x)==3 && the.cal[[x]][[1]]$CalTable$NormType==3){
                predict(
                object=the.cal[[x]][[2]],
                newdata=lucas.comp.prep(
                data=valdata,
                spectra.line.table=as.data.frame(
                count.table
                ),
                element.line=x,
                slope.element.lines=the.cal[[x]][[1]][2]$Slope,
                intercept.element.lines=the.cal[[x]][[1]][3]$Intercept,
                norm.min=the.cal[[x]][[1]][1]$CalTable$Min,
                norm.max=the.cal[[x]][[1]][1]$CalTable$Max
                )
                )
            } else if(valDataType()=="Spectra" && cal_type(x)==4 && the.cal[[x]][[1]]$CalTable$NormType==1){
                predict(
                object=the.cal[[x]][[2]],
                newdata=lucas.simp.prep(
                spectra.line.table=as.data.frame(
                count.table
                ),
                element.line=x,
                slope.element.lines=variables,
                intercept.element.lines=the.cal[[x]][[1]][3]$Intercept
                )
                )
            } else if(valDataType()=="Spectra" && cal_type(x)==4 && the.cal[[x]][[1]]$CalTable$NormType==2){
                predict(
                object=the.cal[[x]][[2]],
                newdata=lucas.tc.prep(
                data=valdata,
                spectra.line.table=as.data.frame(
                count.table
                ),
                element.line=x,
                slope.element.lines=variables,
                intercept.element.lines=the.cal[[x]][[1]][3]$Intercept
                )
                )
            } else if(valDataType()=="Spectra" && cal_type(x)==4 && the.cal[[x]][[1]]$CalTable$NormType==3){
                predict(
                object=the.cal[[x]][[2]],
                newdata=lucas.comp.prep(
                data=valdata,
                spectra.line.table=as.data.frame(
                count.table
                ),
                element.line=x,
                slope.element.lines=variables,
                intercept.element.lines=the.cal[[x]][[1]][3]$Intercept,
                norm.min=the.cal[[x]][[1]][1]$CalTable$Min,
                norm.max=the.cal[[x]][[1]][1]$CalTable$Max
                )
                )
            } else if(valDataType()=="Net" && cal_type(x)==1 && the.cal[[x]][[1]]$CalTable$NormType==1){
                predict(
                object=the.cal[[x]][[2]],
                newdata=general.prep.net(
                spectra.line.table=as.data.frame(
                count.table
                ),
                element.line=x)
                )
            } else if(valDataType()=="Net" && cal_type(x)==1 && the.cal[[x]][[1]]$CalTable$NormType==2) {
                predict(
                object=the.cal[[x]][[2]],
                newdata=simple.tc.prep.net(
                data=valdata,
                spectra.line.table=as.data.frame(
                count.table
                ),
                element.line=x
                )
                )
            } else if(valDataType()=="Net" && cal_type(x)==1 && the.cal[[x]][[1]]$CalTable$NormType==3) {
                predict(
                object=the.cal[[x]][[2]],
                newdata=simple.comp.prep.net(
                data=valdata,
                spectra.line.table=as.data.frame(
                count.table
                ),
                element.line=x,
                norm.min=the.cal[[x]][[1]][1]$CalTable$Min,
                norm.max=the.cal[[x]][[1]][1]$CalTable$Max
                )
                )
            } else if(valDataType()=="Net" && cal_type(x)==3 && the.cal[[x]][[1]]$CalTable$NormType==1){
                predict(
                object=the.cal[[x]][[2]],
                newdata=lucas.simp.prep.net(
                spectra.line.table=as.data.frame(
                count.table
                ),
                element.line=x,
                slope.element.lines=the.cal[[x]][[1]][2]$Slope,
                intercept.element.lines=the.cal[[x]][[1]][3]$Intercept
                )
                )
            } else if(valDataType()=="Net" && cal_type(x)==3 && the.cal[[x]][[1]]$CalTable$NormType==2){
                predict(
                object=the.cal[[x]][[2]],
                newdata=lucas.tc.prep.net(
                data=valdata,
                spectra.line.table=as.data.frame(
                count.table
                ),
                element.line=x,
                slope.element.lines=the.cal[[x]][[1]][2]$Slope,
                intercept.element.lines=the.cal[[x]][[1]][3]$Intercept
                )
                )
            } else if(valDataType()=="Net" && cal_type(x)==3 && the.cal[[x]][[1]]$CalTable$NormType==3){
                predict(
                object=the.cal[[x]][[2]],
                newdata=lucas.comp.prep.net(
                data=valdata,
                spectra.line.table=as.data.frame(
                count.table
                ),
                element.line=x,
                slope.element.lines=the.cal[[x]][[1]][2]$Slope,
                intercept.element.lines=the.cal[[x]][[1]][3]$Intercept,
                norm.min=the.cal[[x]][[1]][1]$CalTable$Min,
                norm.max=the.cal[[x]][[1]][1]$CalTable$Max
                )
                )
            } else if(valDataType()=="Net" && cal_type(x)==4 && the.cal[[x]][[1]]$CalTable$NormType==1){
                predict(
                object=the.cal[[x]][[2]],
                newdata=lucas.simp.prep.net(
                spectra.line.table=as.data.frame(
                count.table
                ),
                element.line=x,
                slope.element.lines=variables,
                intercept.element.lines=the.cal[[x]][[1]][3]$Intercept
                )
                )
            } else if(valDataType()=="Net" && cal_type(x)==4 && the.cal[[x]][[1]]$CalTable$NormType==2){
                predict(
                object=the.cal[[x]][[2]],
                newdata=lucas.tc.prep.net(
                data=valdata,
                spectra.line.table=as.data.frame(
                count.table
                ),
                element.line=x,
                slope.element.lines=variables,
                intercept.element.lines=the.cal[[x]][[1]][3]$Intercept
                )
                )
            } else if(valDataType()=="Net" && cal_type(x)==4 && the.cal[[x]][[1]]$CalTable$NormType==3){
                predict(
                object=the.cal[[x]][[2]],
                newdata=lucas.comp.prep.net(
                data=valdata,
                spectra.line.table=as.data.frame(
                count.table
                ),
                element.line=x,
                slope.element.lines=variables,
                intercept.element.lines=the.cal[[x]][[1]][3]$Intercept,
                norm.min=the.cal[[x]][[1]][1]$CalTable$Min,
                norm.max=the.cal[[x]][[1]][1]$CalTable$Max
                )
                )
            }
            )
            
            predicted.vector <- unlist(predicted.list)
            
            dim(predicted.vector) <- c(length(count.table$Spectrum), length(elements))
            
            predicted.frame <- data.frame(count.table$Spectrum, predicted.vector)
            
            colnames(predicted.frame) <- c("Spectrum", elements)
            #elements <- elements[order(match(fluorescence.lines$Symbol, elements))]
            
            
            
            predicted.data.table <- predicted.frame
            
            #predicted.values <- t(predicted.values)
            predicted.data.table
            
            
        })
        
        output$myvaltable2 <- renderDataTable({
            
            tableInputValQuant()
            
        })
        
        
        # valtest <- lapply(valelements, function(x) predict(calsList[[x]], as.data.frame(val.line.table[x])))
        
        output$downloadValData <- downloadHandler(
        filename = function() { paste(input$calname, "_ValData", '.csv', sep='', collapse='') },
        content = function(file
        ) {
            write.csv(tableInputValQuant(), file)
        }
        )
        
        
        
        
        
    })

         
         

 })




