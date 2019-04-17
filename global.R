list.of.packages <- c("pbapply", "reshape2", "TTR", "dplyr", "ggtern", "ggplot2", "shiny", "rhandsontable", "random", "data.table", "DT", "shinythemes", "Cairo", "broom", "shinyjs", "gridExtra", "dtplyr", "formattable", "XML", "corrplot", "scales", "rmarkdown", "markdown", "peakPick", "shinyWidgets", "data.table", "baseline", "pls", "prospectr", "doSNOW", "parallel", "caret", "wavelets", "hexView", "GSIF", "gstat", "nnet", "neuralnet", "grid", "gridExtra", "compiler")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

if ("soil.spec" %in% installed.packages()[,"Package"]==FALSE){
    install.packages("soil.spec_2.1.4.tar", repos=NULL, type="source")
}

library(soil.spec)

#options(warn=-1)
#assign("last.warning", NULL, envir = baseenv())

library(shiny)
library(ggplot2)
library(pbapply)
library(reshape2)
library(dplyr)
library(DT)
library(XML)
library(soil.spec)
library(parallel)
library(caret)
library(data.table)
library(compiler)
library(nnet)
library(neuralnet)
library(gridExtra)
library(grid)

Sys.setenv(R_MAX_VSIZE = 16e9)

get_os <- function(){
    sysinf <- Sys.info()
    if (!is.null(sysinf)){
        os <- sysinf['sysname']
        if (os == 'Darwin')
        os <- "osx"
    } else { ## mystery machine
        os <- .Platform$OS.type
        if (grepl("^darwin", R.version$os))
        os <- "osx"
        if (grepl("linux-gnu", R.version$os))
        os <- "linux"
    }
    tolower(os)
}


my.cores <- if(parallel::detectCores()>=3){
                paste0(parallel::detectCores()-2)
            } else if(parallel::detectCores()<=2){
                "1"
            }


my_combine <- function (...)
{
    pad0 <- function(x, len) c(x, rep(0, len - length(x)))
    padm0 <- function(x, len) rbind(x, matrix(0, nrow = len -
    nrow(x), ncol = ncol(x)))
    rflist <- list(...)
    areForest <- sapply(rflist, function(x) inherits(x, "randomForest"))
    if (any(!areForest))
    stop("Argument must be a list of randomForest objects")
    rf <- rflist[[1]]
    classRF <- rf$type == "classification"
    trees <- sapply(rflist, function(x) x$ntree)
    ntree <- sum(trees)
    rf$ntree <- ntree
    nforest <- length(rflist)
    haveTest <- !any(sapply(rflist, function(x) is.null(x$test)))
    vlist <- lapply(rflist, function(x) rownames(importance(x)))
    numvars <- sapply(vlist, length)
    if (!all(numvars[1] == numvars[-1]))
    stop("Unequal number of predictor variables in the randomForest objects.")
    for (i in seq_along(vlist)) {
        if (!all(vlist[[i]] == vlist[[1]]))
        stop("Predictor variables are different in the randomForest objects.")
    }
    haveForest <- sapply(rflist, function(x) !is.null(x$forest))
    if (all(haveForest)) {
        nrnodes <- max(sapply(rflist, function(x) x$forest$nrnodes))
        rf$forest$nrnodes <- nrnodes
        rf$forest$ndbigtree <- unlist(sapply(rflist, function(x) x$forest$ndbigtree))
        rf$forest$nodestatus <- do.call("cbind", lapply(rflist,
        function(x) padm0(x$forest$nodestatus, nrnodes)))
        rf$forest$bestvar <- do.call("cbind", lapply(rflist,
        function(x) padm0(x$forest$bestvar, nrnodes)))
        rf$forest$xbestsplit <- do.call("cbind", lapply(rflist,
        function(x) padm0(x$forest$xbestsplit, nrnodes)))
        rf$forest$nodepred <- do.call("cbind", lapply(rflist,
        function(x) padm0(x$forest$nodepred, nrnodes)))
        tree.dim <- dim(rf$forest$treemap)
        if (classRF) {
            rf$forest$treemap <- array(unlist(lapply(rflist,
            function(x) apply(x$forest$treemap, 2:3, pad0,
            nrnodes))), c(nrnodes, 2, ntree))
        }
        else {
            rf$forest$leftDaughter <- do.call("cbind", lapply(rflist,
            function(x) padm0(x$forest$leftDaughter, nrnodes)))
            rf$forest$rightDaughter <- do.call("cbind", lapply(rflist,
            function(x) padm0(x$forest$rightDaughter, nrnodes)))
        }
        rf$forest$ntree <- ntree
        if (classRF)
        rf$forest$cutoff <- rflist[[1]]$forest$cutoff
    }
    else {
        rf$forest <- NULL
    }
    #
    #Tons of stuff removed here...
    #
    if (classRF) {
        rf$confusion <- NULL
        rf$err.rate <- NULL
        if (haveTest) {
            rf$test$confusion <- NULL
            rf$err.rate <- NULL
        }
    }
    else {
        rf$mse <- rf$rsq <- NULL
        if (haveTest)
        rf$test$mse <- rf$test$rsq <- NULL
    }
    rf
}

#registerDoSNOW(makeCluster(as.numeric(my.cores), type="SOCK"))


options(digits=4)
options(warn=-1)
assign("last.warning", NULL, envir = baseenv())

keep_singles <- function(v){
    v[!(v %in% v[duplicated(v)])]
}

Hodder.v <- function(y)
{
    
    n <-length(y)
    
    for(i in 1:(n-1)) {
        y[i] <- y[i+1] - y[i]
        y[1:(n-1)]
        y <- abs(y)
    }
    y <- c(0, y[1:(n-1)])
    
    return(y)
}

Shepherd.v <- function(y) {
    n <-length(y)
    mean <- mean(y)
    sd <- sd(y)
    
    for(i in 1:n) {
        y[i] <- (y[i]/mean)*sd
    }
    
    return(y)

}

cal.lmsummary <-function(lm.object){
    res<-c(paste(as.character(summary(lm.object)$call),collapse=" "),
    length(lm.object$model),
    summary(lm.object)$r.squared,
    summary(lm.object)$adj.r.squared,
    summary(lm.object)$fstatistic,
    pf(summary(lm.object)$fstatistic[1],summary(lm.object)$fstatistic[2],summary(lm.object)$fstatistic[3],lower.tail=FALSE))
    names(res)<-c("Call","n", "R2","Adj. R2",
    "F-statistic","numdf","dendf","p-value")
    return(res)}


val.lmsummary <-function(lm.object){
    res<-c(paste(as.character(summary(lm.object)$call),collapse=" "),
    lm.object$coefficients[1],
    lm.object$coefficients[2],
    length(lm.object$model),
    summary(lm.object)$coefficients[2,2],
    summary(lm.object)$r.squared,
    summary(lm.object)$adj.r.squared,
    summary(lm.object)$fstatistic,
    pf(summary(lm.object)$fstatistic[1],summary(lm.object)$fstatistic[2],summary(lm.object)$fstatistic[3],lower.tail=FALSE))
    names(res)<-c("Call","Intercept","Slope","n","Slope SE","R2","Adj. R2",
    "F-statistic","numdf","dendf","p-value")
    return(res)}



read_csv_filename_x <- function(filename){
    ret <- read.csv(file=filename, sep=",", header=FALSE)
    return.res <- as.numeric(as.vector(ret$V2[18]))/1000
    return.chan.counts <-as.numeric(as.vector(ret$V1[22:2069]))
    return.Wavenumber <- return.chan.counts*return.res
    return(return.Wavenumber)
}

read_csv_filename_y <- function(filename){
    ret <- read.csv(file=filename, sep=",", header=FALSE)
    return.live.time <- as.numeric(as.vector(ret$V2[10]))
    return.counts <- as.numeric(as.vector(ret$V2[22:2069]))
    return.cps <- return.counts/return.live.time
    return(return.cps)
}


read_csv_net <- function(filepath) {
    
    ret <- read.csv(file=filepath, sep=",", header=TRUE)
    element <- ret$Element
    line <- ret$Line
    net <- ret$Net
    background <- ret$Backgr.
    eline <- paste(element, line, sep="-")
    
    simple.table <- data.frame(net)
    colnames(simple.table) <- NULL
    simple.transpose <- as.data.frame(t(simple.table))
    colnames(simple.transpose) <- eline
    
    simple.transpose
    
}



readDPTData <- function(filepath, filename){
    
    filename <- gsub(".dpt", "", filename)
    
    file <- read.csv(filepath)
    
    file.frame <- data.frame(
    Spectrum=rep(filename, length(file[,1])),
    Wavenumber=file[,1],
    Amplitude=file[,2]
    )
    
    return(file.frame)
    
}

readCSVData <- function(filepath, filename){
    
    filename <- gsub(".csv", "", filename)
    
    file <- read.csv(filepath)
    
    file.frame <- data.frame(
    Spectrum=rep(filename, length(file[,1])),
    Wavenumber=file[,1],
    Amplitude=file[,2]
    )
    
    return(file.frame)
    
}

readOpusData <- function(filepath, filename){
    
    filename <- make.names(filename)
    
    alldata <- soil.spec::read.opus(file=filepath)
    
    thedata <- t(slot(slot(alldata, "data"), "ab"))
    
    file.frame <- data.frame(
    Spectrum=rep(filename, length(thedata[,1])-1),
    Wavenumber=as.numeric(as.vector(substring(rownames(thedata), 2)[2:length(thedata[,1])])),
    Amplitude=as.numeric(as.vector(thedata[,1][2:length(thedata[,1])]))
    )
    
    return(file.frame)
    
}



    
    
    
    


range_subset_ftir <- function(range.frame, data){
    
    new.data <- subset(data, Wavenumber >= range.frame$WaveMin & Wavenumber <= range.frame$WaveMax, drop=TRUE)
    newer.data <- aggregate(new.data, by=list(new.data$Spectrum), FUN=mean, na.rm=TRUE)[,c("Group.1", "Amplitude")]
    colnames(newer.data) <- c("Spectrum", as.character(range.frame$Name))
    newer.data
}

ftir_parse <- function(range.table, data){
    
    choice.lines <- range.table[complete.cases(range.table),]
    
    choice.list <- split(choice.lines, f=choice.lines$Name)
    names(choice.list) <- choice.lines[,"Name"]
    
    index <- choice.lines[,"Name"]
    
    selected.list <- lapply(index, function(x) range_subset_ftir(range.frame=choice.list[[x]], data=data))
    
    Reduce(function(...) merge(..., all=T), selected.list)
}


spectra_frame <- function(spectra){
    
    data <- reshape2::dcast(spectra, Spectrum~Wavenumber)
    
    #test <- apply(test, 2, as.numeric)
    colnames(data) <- make.names(colnames(data))
    data[complete.cases(data),]
}


spectra_table <- function(spectra, concentration){
    
    data <- reshape2::dcast(spectra, Spectrum~Wavenumber)
    data$Concentration <- concentration
    
    #test <- apply(test, 2, as.numeric)
    colnames(data) <- make.names(colnames(data))
    data[complete.cases(data),]
}

###Train Functions

pull_test <- function(a.vector, a.value.position){
    
    scaled <- scale(a.vector)[,1]
    
    value <- scaled[a.value.position]
    scale.vector <- scaled[-a.value.position]
    
    ZScore <- (value-mean(scale.vector))/sd(scale.vector)
    pvalue <- pnorm(-abs(ZScore))
    is.sig <- pvalue < 0.05
    
    data.frame(Value=a.vector[a.value.position], ZScore=ZScore, pvalue=pvalue, Sig=is.sig)
}


Z_frame <- function(a.vector){
    
    do.call("rbind", lapply(seq(1, length(a.vector), 1), function(x) pull_test(a.vector, x)))
}


Z_choose <- function(a.vector){
    
    full <- Z_frame(a.vector)
    full[full$Sig,]
    
}

variable_select <- function(amplitudes, values, analyte){
    
    control <- trainControl(method="cv", number=5)
    seed <- 7
    metric <- "RMSE"
    set.seed(seed)
    
    cal.table <- data.frame(amplitudes, Concentration=values[,analyte])
    fit.lm <- train(Concentration~., data=cal.table, method="lm", metric=metric, preProc=c("center", "scale"), trControl=control)
    importance <- varImp(fit.lm, scale=FALSE)
    importance.frame <- as.data.frame(importance$importance)
    elements <- rownames(importance$importance)
    elements[as.numeric(rownames(Z_choose(importance.frame$Overall)))]
    
}


variable_select_short <- function(importance){
    importance.frame <- as.data.frame(importance$importance)
    lines <- rownames(importance$importance)
    lines[as.numeric(rownames(Z_choose(importance.frame$Overall)))]
}




file.0 <- function(file) {
    if (length(file) > 0)
    {
    return(file)
    }else{
        return(levels(file))
    }
}

is.0 <- function(amplitude, file) {
    file.0 <- function(file) {
        if (length(file) > 0)
        {
            return(file)
        }else{
            return(levels(file))
        }
    }
    if (length(amplitude) > 0)
    {
        hope <-data.frame(amplitude, file.0(file))
        return(hope)
    } else {
        empty <- rep(0, length(file.0(file)))
        framed <- data.frame(empty, file.0(file))
        return(framed)
    }
}

dt_options <- reactive({
    # dynamically create options for `aoColumns` depending on how many columns are selected.
    toggles <- lapply(1:length(input$show_vars), function(x) list(bSearchable = F))
    # for `species` columns
    toggles[[length(toggles) + 1]] <- list(bSearchable = T)
    
    list(
    aoColumns = toggles,
    bFilter = 1, bSortClasses = 1,
    aLengthMenu = list(c(10,25,50, -1), list('10','25', '50', 'Todas')),
    iDisplayLength = 10
    )
})

ifrm <- function(obj, env = globalenv()) {
    obj <- deparse(substitute(obj))
    if(exists(obj, envir = env)) {
        rm(list = obj, envir = env)
    }
}

lmp <- function (modelobject) {
    if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
    f <- summary(modelobject)$fstatistic
    p <- pf(f[1],f[2],f[3],lower.tail=F)
    attributes(p) <- NULL
    return(p)
}

lm_eqn.old <- function(df){
    m <- lm(y ~ x, df);
    eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2,
    list(a = format(coef(m)[1], digits = 2),
    b = format(coef(m)[2], digits = 2),
    r2 = format(summary(m)$r.squared, digits = 3)))
    as.character(as.expression(eq));
}


lm_eqn = function(m) {
    
    l <- list(a = as.numeric(format(coef(m)[1], digits = 2)),
    b = as.numeric(format(abs(coef(m)[2]), digits = 2)),
    r2 = format(summary(m)$r.squared, digits = 3));
    
    eq <- substitute(italic(C)[i] == a + b %.% italic(I)[i]*","~~italic(r)^2~"="~r2,l)
    
    
    as.character(as.expression(eq));
}

lm_eqn_poly = function(m) {
    
    l <- list(a = as.numeric(format(coef(m)[1], digits = 2)),
    b = as.numeric(format(abs(coef(m)[2]), digits = 2)),
    c = as.numeric(format(abs(coef(m)[3]), digits = 2)),
    r2 = format(summary(m)$r.squared, digits = 3));
    
    eq <- substitute(italic(C)[i] == a + c %.% italic(I)[i]^2 + b %.% italic(I)[i]*","~~italic(r)^2~"="~r2,l)
    
    
    as.character(as.expression(eq));
}

lm_eqn_val = function(m) {
    
    l <- list(a = as.numeric(format(coef(m)[1], digits = 2)),
    b = as.numeric(format(abs(coef(m)[2]), digits = 2)),
    r2 = format(summary(m)$r.squared, digits = 3));
    
    eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2,l)
    
    
    as.character(as.expression(eq));
}

numericInput2<-function (inputId, label, value = "",...)
{
    div(style="display:inline-block",
    tags$label(label, `for` = inputId),
    tags$input(id = inputId, type = "text", value = value,...))
}

numericInputRow<-function (inputId, label, min, max,  value = "")
{
    div(style="display:inline-block",
    tags$label(label, `for` = inputId),
    tags$input(id = inputId, type = "text", value = value, class="input-mini", width='20%'))
}



diagPlot<-function(model){
    p1<-ggplot(model, aes(as.vector(.fitted), as.vector(.resid)))+geom_point()
    p1<-p1+stat_smooth(method="loess")+geom_hline(yintercept=0, col="red", linetype="dashed")
    p1<-p1+xlab("Fitted values")+ylab("Residuals")
    p1<-p1+ggtitle("Residual vs Fitted Plot")+theme_light()
    
    p2 <- ggplot(model, aes(qqnorm(.stdresid)[[1]], .stdresid))+geom_point(na.rm = TRUE)
    p2 <- p2+geom_abline()+xlab("Theoretical Quantiles")+ylab("Standardized Residuals")
    p2 <- p2+ggtitle("Normal Q-Q")+theme_bw()
    p2
    
    p3<-ggplot(model, aes(as.vector(.fitted), sqrt(abs(as.vector(.stdresid)))))+geom_point(na.rm=TRUE)
    p3<-p3+stat_smooth(method="loess", na.rm = TRUE)+xlab("Fitted Value")
    p3<-p3+ylab(expression(sqrt("|Standardized residuals|")))
    p3<-p3+ggtitle("Scale-Location")+theme_light()
    
    p4<-ggplot(model, aes(seq_along(as.vector(.cooksd)), as.vector(.cooksd)))+geom_bar(stat="identity", position="identity")
    p4<-p4+xlab("Obs. Number")+ylab("Cook's distance")
    p4<-p4+ggtitle("Cook's distance")+theme_light()
    
    p5<-ggplot(model, aes(as.vector(.hat), as.vector(.stdresid)))+geom_point(aes(size=as.vector(.cooksd)), na.rm=TRUE)
    p5<-p5+stat_smooth(method="loess", na.rm=TRUE)
    p5<-p5+xlab("Leverage")+ylab("Standardized Residuals")
    p5<-p5+ggtitle("Residual vs Leverage Plot")
    p5<-p5+scale_size_continuous("Cook's Distance", range=c(1,5))
    p5<-p5+theme_light()+theme(legend.position="bottom")
    
    p6<-ggplot(model, aes(as.vector(.hat), as.vector(.cooksd)))+geom_point(na.rm=TRUE)+stat_smooth(method="loess", na.rm=TRUE)
    p6<-p6+xlab("Leverage hii")+ylab("Cook's Distance")
    p6<-p6+ggtitle("Cook's dist vs Leverage hii/(1-hii)")
    p6<-p6+geom_abline(slope=seq(0,3,0.5), color="gray", linetype="dashed")
    p6<-p6+theme_light()
    
    return(list(rvfPlot=p1, qqPlot=p2, sclLocPlot=p3, cdPlot=p4, rvlevPlot=p5, cvlPlot=p6))
}

rbind.match.columns <- function(input1, input2) {
    n.input1 <- ncol(input1)
    n.input2 <- ncol(input2)
    
    if (n.input2 < n.input1) {
        TF.names <- which(names(input2) %in% names(input1))
        column.names <- names(input2[, TF.names])
    } else {
        TF.names <- which(names(input1) %in% names(input2))
        column.names <- names(input1[, TF.names])
    }
    
    return(rbind(input1[, column.names], input2[, column.names]))
}


strip_glm <- function(cm) {
    cm$y = c()
    cm$model = c()
    
    cm$residuals = c()
    cm$fitted.values = c()
    cm$effects = c()
    cm$qr$qr = c()
    cm$linear.predictors = c()
    cm$weights = c()
    cm$prior.weights = c()
    cm$data = c()
    
    
    cm$family$variance = c()
    cm$family$dev.resids = c()
    cm$family$aic = c()
    cm$family$validmu = c()
    cm$family$simulate = c()
    attr(cm$terms,".Environment") = c()
    attr(cm$formula,".Environment") = c()
    
    cm
}

merge_Sum <- function(.df1, .df2, .id_Columns, .match_Columns){
    merged_Columns <- unique(c(names(.df1),names(.df2)))
    merged_df1 <- data.frame(matrix(nrow=nrow(.df1), ncol=length(merged_Columns)))
    names(merged_df1) <- merged_Columns
    for (column in merged_Columns){
        if(column %in% .id_Columns | !column %in% names(.df2)){
            merged_df1[, column] <- .df1[, column]
        } else if (!column %in% names(.df1)){
            merged_df1[, column] <- .df2[match(.df1[, .match_Columns],.df2[, .match_Columns]), column]
        } else {
            df1_Values=.df1[, column]
            df2_Values=.df2[match(.df1[, .match_Columns],.df2[, .match_Columns]), column]
            df2_Values[is.na(df2_Values)] <- 0
            merged_df1[, column] <- df1_Values + df2_Values
        }
    }
    return(merged_df1)
}


GG_save_pdf = function(list, filename) {
    #start pdf
    pdf(filename)
    
    #loop
    for (p in list) {
        print(p)
    }
    
    #end pdf
    dev.off()
    
    invisible(NULL)
}



####Cal Models

linear_simp_ftir <- function(concentration.table, spectra.line.table, element.line) {
    
    
    
    concentration <- na.omit(as.vector(as.numeric(unlist(concentration.table[element.line]))))
    
    
    Amplitude <- na.omit(as.vector(as.numeric(unlist(spectra.line.table[element.line]))))
    
    
    predict.frame <- data.frame(concentration, Amplitude)
    colnames(predict.frame) <- c("Concentration", "Amplitude")
    
    
    
    predict.amplitude <- data.frame(predict.frame$Amplitude)
    colnames(predict.amplitude) <- c("Amplitude")
    
    cal.lm <- lm(predict.frame$Concentration~predict.frame$Amplitude)
    
    cal.lm
    
}

poly_simp_ftir <- function(concentration.table, spectra.line.table, element.line) {
    
    
    
    concentration <- na.omit(as.vector(as.numeric(unlist(concentration.table[element.line]))))
    
    
    Amplitude <- na.omit(as.vector(as.numeric(unlist(spectra.line.table[element.line]))))
    
    
    predict.frame <- data.frame(concentration, Amplitude)
    colnames(predict.frame) <- c("Concentration", "Amplitude")
    
    
    
    predict.amplitude <- data.frame(predict.frame$Amplitude)
    colnames(predict.amplitude) <- c("Amplitude")
    
    cal.lm.poly <- lm(predict.frame$Concentration~poly(predict.frame$Amplitude, 2))
    
    cal.lm.poly
    
}

lucas_simp_ftir <- function(concentration.table, spectra.line.table, element.line, slope.element.lines, intercept.element.lines) {
    
    
    concentration <- na.omit(as.vector(as.numeric(unlist(concentration.table[element.line]))))
    
    
    Amplitude <- na.omit(as.vector(as.numeric(unlist(spectra.line.table[element.line]))))
    
    lucas.intercept.table <- data.frame(rowSums(lucas.intercept.table.x[intercept.element.lines]))
    colnames(lucas.intercept.table) <- c("first")
    
    
    
    lucas.intercept <- lucas.intercept.table$first
    lucas.slope <- data.frame(lucas.slope.table[slope.element.lines])
    
    
    
    predict.frame.luk <- data.frame(concentration, ((1+Amplitude/(Amplitude+lucas.intercept))-lucas.intercept/(Amplitude+lucas.intercept)),lucas.slope)
    colnames(predict.frame.luk) <- c("Concentration", "Amplitude", names(lucas.slope))
    
    
    
    predict.amplitude.luk <- data.frame(predict.frame.luk$Amplitude, lucas.slope)
    colnames(predict.amplitude.luk) <- c("Amplitude", names(lucas.slope))
    
    lucas.lm <- lm(Concentration~., data=predict.frame.luk)
    
    lucas.lm
    
    
}


linear_tc_ftir <- function(concentration.table, spectra.line.table, element.line) {
    
    
    
    concentration <- na.omit(as.vector(as.numeric(unlist(concentration.table[element.line]))))
    
    
    Amplitude <- na.omit(as.vector(as.numeric(unlist(spectra.line.table[element.line]))))
    
    
    total.counts <- aggregate(Amplitude~Spectrum, data=data, sum)
    colnames(total.counts) <- c("Spectrum", "Amplitude")
    
    
    
    predict.frame.tc <- data.frame(concentration, Amplitude/total.counts$Amplitude)
    colnames(predict.frame.tc) <- c("Concentration", "Amplitude")
    
    
    
    predict.amplitude.tc <- data.frame(predict.frame.tc$Amplitude)
    colnames(predict.amplitude.tc) <- c("Amplitude")
    
    cal.lm.tc <- lm(predict.frame.tc$Concentration~predict.frame.tc$Amplitude)
    
    cal.lm.tc
    
}

poly_tc_ftir <- function(concentration.table, spectra.line.table, element.line) {
    
    
    
    concentration <- na.omit(as.vector(as.numeric(unlist(concentration.table[element.line]))))
    
    
    Amplitude <- na.omit(as.vector(as.numeric(unlist(spectra.line.table[element.line]))))
    
    
    
    
    total.counts <- aggregate(Amplitude~Spectrum, data=data, sum)
    colnames(total.counts) <- c("Spectrum", "Amplitude")
    
    
    
    predict.frame.tc <- data.frame(concentration, Amplitude/total.counts$Amplitude)
    colnames(predict.frame.tc) <- c("Concentration", "Amplitude")
    
    
    
    predict.amplitude.tc <- data.frame(predict.frame.tc$Amplitude)
    colnames(predict.amplitude.tc) <- c("Amplitude")
    
    cal.lm.poly.tc <- lm(predict.frame.tc$Concentration~poly(predict.frame.tc$Amplitude, 2))
    
    cal.lm.poly.tc
    
    
    
}




lucas_tc_ftir <- function(concentration.table, spectra.line.table, element.line, slope.element.lines, intercept.element.lines) {
    
    
    concentration <- na.omit(as.vector(as.numeric(unlist(concentration.table[element.line]))))
    
    
    Amplitude <- na.omit(as.vector(as.numeric(unlist(spectra.line.table[element.line]))))
    
    lucas.intercept.table.tc <- data.frame(rowSums(lucas.intercept.table.x[intercept.element.lines]))/total.counts$Amplitude
    colnames(lucas.intercept.table.tc) <- c("first")
    
    
    
    lucas.intercept.tc <- lucas.intercept.table.tc$first
    lucas.slope.tc <- data.frame(lucas.slope.table[slope.element.lines])/total.counts$Amplitude
    
    
    
    predict.frame.luc.tc <- data.frame(concentration, ((Amplitude/total.counts$Amplitude-lucas.intercept.tc)/(Amplitude/total.counts$Amplitude+lucas.intercept.tc)),lucas.slope.tc)
    colnames(predict.frame.luc.tc) <- c("Concentration", "Amplitude", names(lucas.slope.tc))
    
    
    
    predict.amplitude.luc.tc <- data.frame(predict.frame.luc.tc$Amplitude, lucas.slope.tc)
    colnames(predict.amplitude.luc.tc) <- c("Amplitude", names(lucas.slope.tc))
    
    lucas.lm.tc <- lm(Concentration~., data=predict.frame.luc.tc)
    
    lucas.lm.tc
    
    
}

linear_comp_ftir <- function(data, concentration.table, spectra.line.table, element.line) {
    
    
    
    concentration <- na.omit(as.vector(as.numeric(unlist(concentration.table[element.line]))))
    
    
    Amplitude <- na.omit(as.vector(as.numeric(unlist(spectra.line.table[element.line]))))
    
    
    compton.norm <- subset(data$Amplitude, !(data$Wavenumber < input$comptonmin | data$Wavenumber > input$comptonmax))
    compton.file <- subset(data$Spectrum, !(data$Wavenumber < input$comptonmin | data$Wavenumber > input$comptonmax))
    compton.frame <- data.frame(is.0(compton.norm, compton.file))
    colnames(compton.frame) <- c("Compton", "Spectrum")
    compton.frame.ag <- aggregate(list(compton.frame$Compton), by=list(compton.frame$Spectrum), FUN="sum")
    colnames(compton.frame.ag) <- c("Spectrum", "Compton")
    
    predict.frame.comp <- data.frame(concentration, Amplitude/compton.frame.ag$Compton)
    colnames(predict.frame.comp) <- c("Concentration", "Amplitude")
    
    
    
    predict.amplitude.comp <- data.frame(predict.frame.comp$Amplitude)
    colnames(predict.amplitude.comp) <- c("Amplitude")
    
    cal.lm.comp <- lm(predict.frame.comp$Concentration~predict.frame.comp$Amplitude)
    
    cal.lm.comp
    
}

poly_comp_ftir <- function(data, concentration.table, spectra.line.table, element.line) {
    
    
    
    concentration <- na.omit(as.vector(as.numeric(unlist(concentration.table[element.line]))))
    
    
    Amplitude <- na.omit(as.vector(as.numeric(unlist(spectra.line.table[element.line]))))
    
    
    compton.norm <- subset(data$Amplitude, !(data$Wavenumber < input$comptonmin | data$Wavenumber > input$comptonmax))
    compton.file <- subset(data$Spectrum, !(data$Wavenumber < input$comptonmin | data$Wavenumber > input$comptonmax))
    compton.frame <- data.frame(is.0(compton.norm, compton.file))
    colnames(compton.frame) <- c("Compton", "Spectrum")
    compton.frame.ag <- aggregate(list(compton.frame$Compton), by=list(compton.frame$Spectrum), FUN="sum")
    colnames(compton.frame.ag) <- c("Spectrum", "Compton")
    
    predict.frame.comp <- data.frame(concentration, Amplitude/compton.frame.ag$Compton)
    colnames(predict.frame.comp) <- c("Concentration", "Amplitude")
    
    
    
    predict.amplitude.comp <- data.frame(predict.frame.comp$Amplitude)
    colnames(predict.amplitude.comp) <- c("Amplitude")
    
    cal.lm.poly.comp <- lm(predict.frame.comp$Concentration~poly(predict.frame.comp$Amplitude, 2))
    
    cal.lm.poly.comp
    
}

lucas_comp_ftir <- function(data, concentration.table, spectra.line.table, element.line, slope.element.lines, intercept.element.lines) {
    
    
    concentration <- na.omit(as.vector(as.numeric(unlist(concentration.table[element.line]))))
    
    
    Amplitude <- na.omit(as.vector(as.numeric(unlist(spectra.line.table[element.line]))))
    
    
    compton.norm <- subset(data$Amplitude, !(data$Wavenumber < input$comptonmin | data$Wavenumber > input$comptonmax))
    compton.file <- subset(data$Spectrum, !(data$Wavenumber < input$comptonmin | data$Wavenumber > input$comptonmax))
    compton.frame <- data.frame(is.0(compton.norm, compton.file))
    colnames(compton.frame) <- c("Compton", "Spectrum")
    compton.frame.ag <- aggregate(list(compton.frame$Compton), by=list(compton.frame$Spectrum), FUN="sum")
    colnames(compton.frame.ag) <- c("Spectrum", "Compton")
    
    
    
    lucas.intercept.table.comp <- data.frame(rowSums(lucas.intercept.table.x[intercept.element.lines]))/compton.frame.ag$Compton
    colnames(lucas.intercept.table.comp) <- c("first")
    
    
    
    lucas.intercept.comp <- lucas.intercept.table.comp$first
    lucas.slope.comp <- data.frame(lucas.slope.table[slope.element.lines])/compton.frame.ag$Compton
    
    
    
    
    predict.frame.luc.comp <- data.frame(concentration, ((1+Amplitude/compton.frame.ag$Compton)/(Amplitude/compton.frame.ag$Compton+lucas.intercept.comp)-lucas.intercept.comp/(Amplitude/compton.frame.ag$Compton+lucas.intercept.comp)),lucas.slope.comp)
    colnames(predict.frame.luc.comp) <- c("Concentration", "Amplitude", names(lucas.slope.comp))
    
    
    
    predict.amplitude.luc.comp <- data.frame(predict.frame.luc.comp$Amplitude, lucas.slope.comp)
    colnames(predict.amplitude.luc.comp) <- c("Amplitude", names(lucas.slope.comp))
    
    lucas.lm.comp <- lm(Concentration~., data=predict.frame.luc.comp)
    
    lucas.lm.comp
    
}



###############
###Prep Data###
###############


###############
###Full Spectra##
###############


spectra_frame_ftir <- function(spectra){
    
    data <- reshape2::dcast(spectra, Spectrum~Wavenumber)
    
    #test <- apply(test, 2, as.numeric)
    colnames(data) <- make.names(colnames(data))
    data[complete.cases(data),]
}


spectra_table_ftir <- function(spectra, concentration){
    
    data <- reshape2::dcast(spectra, Spectrum~Wavenumber)
    data$Concentration <- concentration
    
    #test <- apply(test, 2, as.numeric)
    colnames(data) <- make.names(colnames(data))
    data[complete.cases(data),]
}




spectra_simp_prep_ftir <- function(spectra, compression=0){
    
    spectra$Wavenumber <- round(spectra$Wavenumber, compression)
    spectra <- data.table(spectra)
    spectra.aggregate <- spectra[, list(Amplitude=mean(Amplitude, na.rm = TRUE)), by = list(Spectrum,Wavenumber)]
    
    data <- as.data.frame(dcast.data.table(spectra.aggregate, Spectrum~Wavenumber, value.var="Amplitude"))
    
    #test <- apply(test, 2, as.numeric)
    colnames(data) <- make.names(colnames(data))
    do.call(data.frame,lapply(data, function(x) replace(x, is.infinite(x),0)))
    
}

spectra_tc_prep_ftir <- function(spectra, compression=0){
    
    spectra$Wavenumber <- round(spectra$Wavenumber, compression)
    
    spectra <- data.table(spectra)
    spectra.aggregate <- spectra[, list(Amplitude=mean(Amplitude, na.rm = TRUE)), by = list(Spectrum,Wavenumber)]
    
    data <- as.data.frame(dcast.data.table(spectra.aggregate, Spectrum~Wavenumber, value.var="Amplitude"))
    
    #test <- apply(test, 2, as.numeric)
    colnames(data) <- make.names(colnames(data))
    data <- data[complete.cases(data),]
    
    total.counts <- rowSums(data[,-1], na.rm=TRUE)
    
    data <- data.frame(Spectrum=data$Spectrum, data[,-1]/total.counts)
    do.call(data.frame,lapply(data, function(x) replace(x, is.infinite(x),0)))
    
    
    
}

spectra_comp_prep_ftir <- function(spectra, compression=0, norm.min, norm.max){
    
    compton.norm <- subset(spectra$Amplitude, !(spectra$Wavenumber < norm.min | spectra$Wavenumber > norm.max))
    compton.file <- subset(spectra$Spectrum, !(spectra$Wavenumber < norm.min | spectra$Wavenumber > norm.max))
    compton.frame <- data.frame(is.0(compton.norm, compton.file))
    colnames(compton.frame) <- c("Compton", "Spectrum")
    compton.frame.ag <- aggregate(list(compton.frame$Compton), by=list(compton.frame$Spectrum), FUN="sum")
    colnames(compton.frame.ag) <- c("Spectrum", "Compton")
    
    spectra$Wavenumber <- round(spectra$Wavenumber, compression)
    
    spectra <- data.table(spectra)
    spectra.aggregate <- spectra[, list(Amplitude=mean(Amplitude, na.rm = TRUE)), by = list(Spectrum,Wavenumber)]
    
    data <- as.data.frame(dcast.data.table(spectra.aggregate, Spectrum~Wavenumber, value.var="Amplitude"))
    #test <- apply(test, 2, as.numeric)
    colnames(data) <- make.names(colnames(data))
    
    data <- data.frame(Spectrum=data$Spectrum, data[,-1]/compton.frame.ag$Compton)
    do.call(data.frame,lapply(data, function(x) replace(x, is.infinite(x),0)))
    
}



###############
###Prep Data###
###############


###############
###Raw Spectra##
###############



general_prep_ftir <- function(spectra.line.table, element.line) {
    
    Amplitude <- spectra.line.table[,element.line]
    
    
    data.frame(Amplitude=spectra.line.table[,element.line])

}

simple_tc_prep_ftir <- function(data,spectra.line.table, element.line) {
    
    Amplitude <- spectra.line.table[,element.line]
    
    
    total.counts <- aggregate(Amplitude~Spectrum, data=data, sum)
    colnames(total.counts) <- c("Spectrum", "Amplitude")
    
    
    
    predict.frame.tc <- data.frame(Amplitude/total.counts$Amplitude)
    colnames(predict.frame.tc) <- c("Amplitude")
    
    
    
    predict.amplitude.tc <- data.frame(predict.frame.tc$Amplitude)
    colnames(predict.amplitude.tc) <- c("Amplitude")
    
    predict.amplitude.tc
}


simple_comp_prep_ftir <- function(data, spectra.line.table, element.line, norm.min, norm.max) {
    
    
    
    Amplitude <- spectra.line.table[,element.line]
    
    
    compton.norm <- subset(data$Amplitude, !(data$Wavenumber < norm.min | data$Wavenumber > norm.max))
    compton.file <- subset(data$Spectrum, !(data$Wavenumber < norm.min | data$Wavenumber > norm.max))
    compton.frame <- data.frame(is.0(compton.norm, compton.file))
    colnames(compton.frame) <- c("Compton", "Spectrum")
    compton.frame.ag <- aggregate(list(compton.frame$Compton), by=list(compton.frame$Spectrum), FUN="sum")
    colnames(compton.frame.ag) <- c("Spectrum", "Compton")
    
    predict.frame.comp <- data.frame( Amplitude/compton.frame.ag$Compton)
    colnames(predict.frame.comp) <- c("Amplitude")
    
    
    
    predict.amplitude.comp <- data.frame(predict.frame.comp$Amplitude)
    colnames(predict.amplitude.comp) <- c("Amplitude")
    
    predict.amplitude.comp
    
}





###Prep Data



lucas_simp_prep_ftir <- function(spectra.line.table, element.line, slope.element.lines, intercept.element.lines) {
    
    
    Amplitude <- spectra.line.table[,element.line]
    
    
    intercept.none <- rep(0, length(spectra.line.table[,1]))
    lucas.intercept.table.x <- data.frame(spectra.line.table, intercept.none, intercept.none)
    colnames(lucas.intercept.table.x) <- c(names(spectra.line.table), "None", "NoneNull")
    
    
    
    
    slope.none <- rep(1, length(spectra.line.table[,1]))
    lucas.slope.table <- data.frame(spectra.line.table, slope.none)
    colnames(lucas.slope.table) <- c(names(spectra.line.table), "None")
    
    
    lucas.intercept.table <- data.frame(rowSums(lucas.intercept.table.x[,c(intercept.element.lines, "None", "NoneNull")]))
    colnames(lucas.intercept.table) <- c("first")
    
    
    
    lucas.intercept <- lucas.intercept.table$first
    lucas.slope <- data.frame(lucas.slope.table[,slope.element.lines])
    colnames(lucas.slope) <- slope.element.lines
    
    
    
    predict.frame.luk <- data.frame(Amplitude=((1+Amplitude/(Amplitude+lucas.intercept))-lucas.intercept/(Amplitude+lucas.intercept)),lucas.slope)

    
    predict.frame.luk
    
    
}



lucas_tc_prep_ftir <- function(data, spectra.line.table, element.line, slope.element.lines, intercept.element.lines) {
    
    
    Amplitude <- spectra.line.table[,element.line]
    
    
    total.counts <- aggregate(Amplitude~Spectrum, data=data, sum)
    colnames(total.counts) <- c("Spectrum", "Amplitude")
    
    
    intercept.none <- rep(0, length(spectra.line.table[,1]))
    lucas.intercept.table.x <- data.frame(spectra.line.table, intercept.none, intercept.none)
    colnames(lucas.intercept.table.x) <- c(names(spectra.line.table), "None", "NoneNull")
    
    
    
    
    slope.none <- rep(1, length(spectra.line.table[,1]))
    lucas.slope.table <- data.frame(spectra.line.table, slope.none)
    colnames(lucas.slope.table) <- c(names(spectra.line.table), "None")
    
    
    lucas.intercept.table.tc <- data.frame(rowSums(lucas.intercept.table.x[,c(intercept.element.lines, "None", "NoneNull")]))/total.counts$Amplitude
    colnames(lucas.intercept.table.tc) <- c("first")
    
    
    
    lucas.intercept.tc <- lucas.intercept.table.tc$first
    lucas.slope.tc <- data.frame(lucas.slope.table[,slope.element.lines])/total.counts$Amplitude
    colnames(lucas.slope.tc) <- slope.element.lines
    
    
    
    predict.amplitude.luc.tc <- data.frame(Amplitude=((1+Amplitude/(Amplitude+lucas.intercept.tc)-lucas.intercept.tc/(Amplitude+lucas.intercept.tc))),lucas.slope.tc)
    
    
    predict.amplitude.luc.tc
}





lucas_comp_prep_ftir <- function(data, spectra.line.table, element.line, slope.element.lines, intercept.element.lines, norm.min, norm.max) {
    
    
    Amplitude <- spectra.line.table[,element.line]
    
    
    
    compton.norm <- subset(data$Amplitude, !(data$Wavenumber < norm.min | data$Wavenumber > norm.max))
    compton.file <- subset(data$Spectrum, !(data$Wavenumber < norm.min | data$Wavenumber > norm.max))
    compton.frame <- data.frame(is.0(compton.norm, compton.file))
    colnames(compton.frame) <- c("Compton", "Spectrum")
    compton.frame.ag <- aggregate(list(compton.frame$Compton), by=list(compton.frame$Spectrum), FUN="sum")
    colnames(compton.frame.ag) <- c("Spectrum", "Compton")
    
    
    intercept.none <- rep(0, length(spectra.line.table[,1]))
    lucas.intercept.table.x <- data.frame(spectra.line.table, intercept.none, intercept.none)
    colnames(lucas.intercept.table.x) <- c(names(spectra.line.table), "None", "NoneNull")
    
    
    
    
    slope.none <- rep(1, length(spectra.line.table[,1]))
    lucas.slope.table <- data.frame(spectra.line.table, slope.none)
    colnames(lucas.slope.table) <- c(names(spectra.line.table), "None")
    
    
    
    lucas.intercept.table.comp <- data.frame(rowSums(lucas.intercept.table.x[,c(intercept.element.lines, "None", "NoneNull")])/compton.frame.ag$Compton)
    colnames(lucas.intercept.table.comp) <- c("first")
    
    
    
    lucas.intercept.comp <- lucas.intercept.table.comp$first
    lucas.slope.comp <- data.frame(lucas.slope.table[,slope.element.lines]/compton.frame.ag$Compton)
    colnames(lucas.slope.comp) <- slope.element.lines
    
    
    predict.frame.luc.comp <- data.frame(Amplitude=((1+Amplitude/compton.frame.ag$Compton)/(Amplitude/compton.frame.ag$Compton+lucas.intercept.comp)-lucas.intercept.comp/(Amplitude/compton.frame.ag$Compton+lucas.intercept.comp)),lucas.slope.comp)
    
    
    predict.frame.luc.comp
}





###############
###Prep Data###
###############


###############
###Net Counts##
###############


general_prep_ftir_net <- function(spectra.line.table, element.line) {
    
    Amplitude <- spectra.line.table[,element.line]
    
    
    predict.frame <- data.frame(Amplitude)
    colnames(predict.frame) <- c("Amplitude")
    
    
    
    predict.amplitude <- data.frame(predict.frame$Amplitude)
    colnames(predict.amplitude) <- c("Amplitude")
    
    predict.amplitude
}

simple_tc_prep_ftir_net <- function(data,spectra.line.table, element.line) {
    
    Amplitude <- spectra.line.table[,element.line]
    
    total.counts.net <- rowSums(spectra.line.table[,-1])
    total.counts <- data.frame(data$Spectrum, total.counts.net)
    colnames(total.counts) <- c("Spectrum", "Amplitude")
    
    
    
    predict.frame.tc <- data.frame(Amplitude/total.counts$Amplitude)
    colnames(predict.frame.tc) <- c("Amplitude")
    
    
    
    predict.amplitude.tc <- data.frame(predict.frame.tc$Amplitude)
    colnames(predict.amplitude.tc) <- c("Amplitude")
    
    predict.amplitude.tc
}


simple_comp_prep_ftir_net <- function(data, spectra.line.table, element.line, norm.min, norm.max) {
    
    
    
    Amplitude <- spectra.line.table[,element.line]
    
    
    compton.ag.fake.Spectrum <- data$Spectrum
    compton.ag.fake.Compton <- rep(1, length(data$Spectrum))
    compton.ag.fake <- data.frame(compton.ag.fake.Spectrum,compton.ag.fake.Compton)
    colnames(compton.ag.fake) <- c("Spectrum", "Compton")
    
    predict.frame.comp <- data.frame( Amplitude/compton.ag.fake$Compton)
    colnames(predict.frame.comp) <- c("Amplitude")
    
    
    
    predict.amplitude.comp <- data.frame(predict.frame.comp$Amplitude)
    colnames(predict.amplitude.comp) <- c("Amplitude")
    
    predict.amplitude.comp
    
}



###Prep Data



lucas_simp_prep_ftir_net <- function(spectra.line.table, element.line, slope.element.lines, intercept.element.lines) {
    
    
    Amplitude <- spectra.line.table[,element.line]
    
    intercept.none <- rep(0, length(spectra.line.table[,1]))
    lucas.intercept.table.x <- data.frame(spectra.line.table, intercept.none, intercept.none)
    colnames(lucas.intercept.table.x) <- c(names(spectra.line.table), "None", "NoneNull")
    
    
    
    
    slope.none <- rep(1, length(spectra.line.table[,1]))
    lucas.slope.table <- data.frame(spectra.line.table, slope.none)
    colnames(lucas.slope.table) <- c(names(spectra.line.table), "None")
    
    
    lucas.intercept.table <- data.frame(rowSums(lucas.intercept.table.x[,c(intercept.element.lines, "None", "NoneNull")]))
    colnames(lucas.intercept.table) <- c("first")
    
    
    
    lucas.intercept <- lucas.intercept.table$first
    lucas.slope <- data.frame(lucas.slope.table[,slope.element.lines])
    colnames(lucas.slope) <- slope.element.lines
    
    
    
    predict.frame.luk <- data.frame(((1+Amplitude/(Amplitude+lucas.intercept))-lucas.intercept/(Amplitude+lucas.intercept)),lucas.slope)
    colnames(predict.frame.luk) <- c("Amplitude", names(lucas.slope))
    
    
    
    predict.amplitude.luk <- data.frame(predict.frame.luk$Amplitude, lucas.slope)
    colnames(predict.amplitude.luk) <- c("Amplitude", names(lucas.slope))
    
    predict.amplitude.luk
    
    
}



lucas_tc_prep_ftir_net <- function(data, spectra.line.table, element.line, slope.element.lines, intercept.element.lines) {
    
    
    Amplitude <- spectra.line.table[,element.line]
    
    
    total.counts.net <- rowSums(spectra.line.table[,-1])
    total.counts <- data.frame(data$Spectrum, total.counts.net)
    colnames(total.counts) <- c("Spectrum", "Amplitude")
    
    
    
    
    intercept.none <- rep(0, length(spectra.line.table[,1]))
    lucas.intercept.table.x <- data.frame(spectra.line.table, intercept.none, intercept.none)
    colnames(lucas.intercept.table.x) <- c(names(spectra.line.table), "None", "NoneNull")
    
    
    
    
    slope.none <- rep(1, length(spectra.line.table[,1]))
    lucas.slope.table <- data.frame(spectra.line.table, slope.none)
    colnames(lucas.slope.table) <- c(names(spectra.line.table), "None")
    
    
    
    lucas.intercept.table.tc <- data.frame(rowSums(lucas.intercept.table.x[,c(intercept.element.lines, "None", "NoneNull")]))/total.counts$Amplitude
    colnames(lucas.intercept.table.tc) <- c("first")
    
    
    
    
    lucas.intercept.tc <- lucas.intercept.table.tc$first
    lucas.slope.tc <- data.frame(lucas.slope.table[,slope.element.lines])/total.counts$Amplitude
    colnames(lucas.slope.tc) <- slope.element.lines
    
    
    predict.amplitude.luc.tc <- data.frame(((1+Amplitude/(Amplitude+lucas.intercept.tc)-lucas.intercept.tc/(Amplitude+lucas.intercept.tc))),lucas.slope.tc)
    colnames(predict.amplitude.luc.tc) <- c("Amplitude", names(lucas.slope.tc))
    
    
    predict.amplitude.luc.tc
}


lucas_comp_prep_ftir_net <- function(data, spectra.line.table, element.line, slope.element.lines, intercept.element.lines, norm.min, norm.max) {
    
    
    Amplitude <- spectra.line.table[,element.line]
    
    
    
    
    compton.ag.fake.Spectrum <- data$Spectrum
    compton.ag.fake.Compton <- rep(1, length(data$Spectrum))
    compton.ag.fake <- data.frame(compton.ag.fake.Spectrum,compton.ag.fake.Compton)
    colnames(compton.ag.fake) <- c("Spectrum", "Compton")
    
    
    intercept.none <- rep(0, length(spectra.line.table[,1]))
    lucas.intercept.table.x <- data.frame(spectra.line.table, intercept.none, intercept.none)
    colnames(lucas.intercept.table.x) <- c(names(spectra.line.table), "None", "NoneNull")
    
    
    
    
    slope.none <- rep(1, length(spectra.line.table[,1]))
    lucas.slope.table <- data.frame(spectra.line.table, slope.none)
    colnames(lucas.slope.table) <- c(names(spectra.line.table), "None")
    
    
    
    lucas.intercept.table.comp <- data.frame(rowSums(lucas.intercept.table.x[,c(intercept.element.lines, "None", "NoneNull")]))/compton.ag.fake$Compton
    colnames(lucas.intercept.table.comp) <- c("first")
    
    
    
    
    lucas.intercept.comp <- lucas.intercept.table.comp$first
    lucas.slope.comp <- data.frame(lucas.slope.table[,slope.element.lines])/compton.ag.fake$Compton
    colnames(lucas.slope.comp) <- slope.element.lines
    
    
    
    predict.frame.luc.comp <- data.frame(((1+predict.frame.comp$Amplitude/(predict.frame.comp$Amplitude+lucas.intercept.comp)-lucas.intercept.comp/(predict.frame.comp$Amplitude+lucas.intercept.comp))),lucas.slope.comp)
    colnames(predict.frame.luc.comp) <- c("Amplitude", names(lucas.slope.comp))
    
    
    
    predict.amplitude.luc.comp <- data.frame(predict.frame.luc.comp$Amplitude, lucas.slope.comp)
    colnames(predict.amplitude.luc.comp) <- c("Amplitude", names(lucas.slope.comp))
    
    
    predict.amplitude.luc.comp
}


in_range_ftir <- function(spectrum, peak, pritable){
    
    temp <- pritable[(pritable$Max>=peak & pritable$Min<=peak),]
    temp <- temp[complete.cases(temp[,c("id", "Peak.Range", "Wavenumbers", "Max", "Min", "Type", "General")]),]
    temp$Spectrum <- rep(spectrum, length(temp[,1]))
    temp$Peak <- rep(peak, length(temp[,1]))
    
    temp
    
}



####Machine Learning Functions


create_frame_slopes_ftir <- function(element, slopes, values, amplitudes){
    values <- values[complete.cases(values[,element]),]
    amplitudes <- amplitudes[complete.cases(values[,element]),]
    
    data.frame(Value=values[,element],
    Amplitude=amplitudes[,"Amplitude"],
    amplitudes[,slopes])
    
}

create_frame_intercepts_ftir <- function(element, slopes, values, amplitudes){
    
    data.frame(Value=values[,element],
    Amplitude=amplitudes[,"Amplitude"],
    amplitudes[,slopes])
    
}



optimal_r_chain_ftir <- function(element, amplitudes, values, possible.slopes, keep){
    
    values <- values[complete.cases(values[,element]),]
    amplitudes <- amplitudes[complete.cases(values[,element]),]
    index <- seq(1, length(possible.slopes), 1)
    
    chain.lm <- pbapply::pblapply(possible.slopes, function(x) lm(Value~Amplitude+., data=create.frame.slopes(element=element, slopes=x, values=values[keep,], amplitudes=amplitudes)[keep,]))
    
    #chain.predict <- pblapply(index, function(x) predict(object=chain.lm[[x]], newdata=create.frame.slopes(element=element, slopes=possible.slopes[[x]], values=values[keep,], intensities=intensities)[keep,], interval='confidence'))
    #chain.fits <- pblapply(chain.predict, function(x) data.frame(x)$fit)
    #val.lm <- pblapply(chain.fits, function(x) lm(values[,element]~x))
    
    aic <- lapply(chain.lm, function(x) extractAIC(x, k=log(length(possible.slopes)))[2])
    best <- chain.lm[[which.min(unlist(aic))]]
    best.aic <- unlist(aic)[which.min(unlist(aic))]
    #r.adj <- lapply(chain.lm, function(x) summary(x)$adj.r.squared)
    #best <- chain.lm[[which.max(unlist(r.adj))]]
    coef <- data.frame(best$coefficients)
    best.var <- rownames(coef)[3:length(rownames(coef))]
    
    simple_lm_ftir <- lm(Value~Amplitude, data=create_frame_slopes_ftir(element=element, slopes=element, values=values, amplitudes=amplitudes)[keep,])
    #simple.predict <- as.data.frame(predict(simple.lm, newdata=create.frame.slopes(element=element, slopes=element, values=values[keep,], intensities=intensities)[keep,], interval='confidence'), interval='confidence')$fit
    #simple.val <- lm(values[,element]~simple.predict)
    simple.aic <- extractAIC(simple.lm, k=log(length(1)))[2]
    
    if(simple.aic <= best.aic){
        element
    } else if(best.aic < simple.aic){
        best.var
    }
    
    #best.var
}


optimal_norm_chain <- function(data, element, spectra.line.table, values, possible.mins, possible.maxs){
    
    index <- seq(1, length(possible.mins), 1)
    
    chain.lm <- pbapply::pblapply(index, function(x) lm(values[,element]~simple_comp_prep_ftir(data=data, spectra.line.table=spectra.line.table, element.line=element, norm.min=possible.mins[x], norm.max=possible.maxs[x])$Amplitude, na.action=na.exclude))
    aic <- lapply(chain.lm, function(x) extractAIC(x, k=log(length(1)))[2])
    best <- index[[which.min(unlist(aic))]]
    
    
    best
    
}


optimal_intercept_chain_ftir <- function(element, amplitudes, values, keep){
    
    
    chain.lm <- pbapply::pblapply(amplitudes, function(x) lm(values[,element]~Amplitude, data=x[keep,]))
    aic <- lapply(chain.lm, function(x) extractAIC(x, k=log(1))[2])
    best <- chain.lm[[which.min(unlist(aic))]]
    coef <- data.frame(best$coefficients)
    best.var <- rownames(coef)[3:length(rownames(coef))]
    
    best.var
    
}



####Data
pritable <- read.csv(file="data/FTIR Master Table - Sheet1.csv")

pritable.range.index <- as.vector(is.na(pritable$Max))

###Machine Learning UI

forestTryUI <- function(radiocal=3, neuralhiddenlayers=1, selection, maxsample){
    if(radiocal==1){
        NULL
    } else if(radiocal==2){
        NULL
    } else if(radiocal==3){
        NULL
    } else if(radiocal==4){
        sliderInput("foresttry", label="Sampling", min=2, max=maxsample-2, value=selection)
    }  else if(radiocal==5){
        sliderInput("foresttry", label="Sampling", min=2, max=maxsample-2, value=selection)
    } else if(radiocal==6 && neuralhiddenlayers == 1){
        NULL
    } else if(radiocal==6 && neuralhiddenlayers > 1){
        sliderInput("foresttry", label="Sampling", min=2, max=maxsample-2, value=selection)
    } else if(radiocal==7 && neuralhiddenlayers == 1){
        NULL
    } else if(radiocal==7 && neuralhiddenlayers > 1){
        sliderInput("foresttry", label="Sampling", min=2, max=maxsample-2, value=selection)
    }
}

forestMetricUI <- function(radiocal, selection){
    if(radiocal==1){
        NULL
    } else if(radiocal==2){
        NULL
    } else if(radiocal==3){
        NULL
    } else if(radiocal==4){
        selectInput("forestmetric", label="Metric", choices=c("Root Mean Square Error"="RMSE", "R2"="Rsquared", "Mean Absolute Error"="MAE", "Kappa"="Kappa", "Logarithmic Loss"="logLoss"), selected=selection)
    } else if(radiocal==5){
        selectInput("forestmetric", label="Metric", choices=c("Root Mean Square Error"="RMSE", "R2"="Rsquared", "Logarithmic Loss"="logLoss"), selected=selection)
    } else if(radiocal==6){
        selectInput("forestmetric", label="Metric", choices=c("Root Mean Square Error"="RMSE", "R2"="Rsquared", "Logarithmic Loss"="logLoss"), selected=selection)
    } else if(radiocal==7){
        selectInput("forestmetric", label="Metric", choices=c("Root Mean Square Error"="RMSE", "R2"="Rsquared", "Logarithmic Loss"="logLoss"), selected=selection)
    }
}

forestTrainUI <- function(radiocal, selection){
    if(radiocal==1){
        NULL
    } else if(radiocal==2){
        NULL
    } else if(radiocal==3){
        NULL
    } else if(radiocal==4){
        selectInput("foresttrain", label="Train Control", choices=c("k-fold Cross Validation"="cv", "Bootstrap"="boot", "0.632 Bootstrap"="boot632", "Optimism Bootstrap"="optimism_boot", "Repeated k-fold Cross Validation"="repeatedcv", "Leave One Out Cross Validation"="LOOCV", "Out of Bag Estimation"="oob"), selected=selection)
    }  else if(radiocal==5){
        selectInput("foresttrain", label="Train Control", choices=c("k-fold Cross Validation"="cv", "Bootstrap"="boot", "0.632 Bootstrap"="boot632", "Optimism Bootstrap"="optimism_boot", "Repeated k-fold Cross Validation"="repeatedcv", "Leave One Out Cross Validation"="LOOCV", "Out of Bag Estimation"="oob"), selected=selection)
    } else if(radiocal==6){
        selectInput("foresttrain", label="Train Control", choices=c("k-fold Cross Validation"="cv", "Bootstrap"="boot", "0.632 Bootstrap"="boot632", "Optimism Bootstrap"="optimism_boot", "Repeated k-fold Cross Validation"="repeatedcv", "Leave One Out Cross Validation"="LOOCV"), selected=selection)
    } else if(radiocal==7){
        selectInput("foresttrain", label="Train Control", choices=c("k-fold Cross Validation"="cv", "Bootstrap"="boot", "0.632 Bootstrap"="boot632", "Optimism Bootstrap"="optimism_boot", "Repeated k-fold Cross Validation"="repeatedcv", "Leave One Out Cross Validation"="LOOCV"), selected=selection)
    }
}

forestNumberUI <- function(radiocal, selection){
    if(radiocal==1){
        NULL
    } else if(radiocal==2){
        NULL
    } else if(radiocal==3){
        NULL
    } else if(radiocal==4){
        sliderInput("forestnumber", label="Iterations", min=5, max=500, value=selection)
    }  else if(radiocal==5){
        sliderInput("forestnumber", label="Iterations", min=5, max=500, value=selection)
    } else if(radiocal==6){
        sliderInput("forestnumber", label="Iterations", min=5, max=500, value=selection)
    } else if(radiocal==7){
        sliderInput("forestnumber", label="Iterations", min=5, max=500, value=selection)
    }
}

forestTreesUI <- function(radiocal, selection){
    if(radiocal==1){
        NULL
    } else if(radiocal==2){
        NULL
    } else if(radiocal==3){
        NULL
    } else if(radiocal==4){
        sliderInput("foresttrees", label="Trees", min=50, max=2000, value=selection)
    } else if(radiocal==5){
        sliderInput("foresttrees", label="Trees", min=50, max=2000, value=selection)
    } else if(radiocal==6){
        NULL
    } else if(radiocal==7){
        NULL
    }
}

neuralHiddenLayersUI <- function(radiocal, selection){
    if(radiocal==1){
        NULL
    } else if(radiocal==2){
        NULL
    } else if(radiocal==3){
        NULL
    } else if(radiocal==4){
        NULL
    }  else if(radiocal==5){
        NULL
    } else if(radiocal==6){
        sliderInput("neuralhiddenlayers", label="Hidden Layers", min=1, max=3, value=selection)
    } else if(radiocal==7){
        sliderInput("neuralhiddenlayers", label="Hidden Layers", min=1, max=3, value=selection)
    }
}

neuralHiddenUnitsUi <- function(radiocal, selection){
    if(radiocal==1){
        NULL
    } else if(radiocal==2){
        NULL
    } else if(radiocal==3){
        NULL
    } else if(radiocal==4){
        NULL
    }  else if(radiocal==5){
        NULL
    } else if(radiocal==6){
        sliderInput("neuralhiddenunits", label="Hidden Units", min=1, max=10, value=selection)
    } else if(radiocal==7){
        sliderInput("neuralhiddenunits", label="Hidden Units", min=1, max=10, value=selection)
    }
}

neuralWeightDecayUI <- function(radiocal, selection, neuralhiddenlayers=1){
    if(radiocal==1){
        NULL
    } else if(radiocal==2){
        NULL
    } else if(radiocal==3){
        NULL
    } else if(radiocal==4){
        NULL
    }  else if(radiocal==5){
        NULL
    } else if(radiocal==6 && neuralhiddenlayers == 1){
        sliderInput("neuralweightdecay", label="Weight Decay", min=0.1, max=0.7, step=0.1, value=selection)
    } else if(radiocal==6 && neuralhiddenlayers > 1){
        NULL
    } else if(radiocal==7 && neuralhiddenlayers == 1){
        sliderInput("neuralweightdecay", label="Weight Decay", min=0.1, max=0.7, step=0.1, value=selection)
    } else if(radiocal==7 && neuralhiddenlayers > 1){
        NULL
    }
}

neuralMaxIterationsUI <- function(radiocal, selection, neuralhiddenlayers=1){
    if(radiocal==1){
        NULL
    } else if(radiocal==2){
        NULL
    } else if(radiocal==3){
        NULL
    } else if(radiocal==4){
        NULL
    }  else if(radiocal==5){
        NULL
    } else if(radiocal==6 && neuralhiddenlayers == 1){
        sliderInput("neuralmaxiterations", label="Max Iterations", min=50, max=2000, value=selection)
    } else if(radiocal==6 && neuralhiddenlayers > 1){
        NULL
    } else if(radiocal==7 && neuralhiddenlayers == 1){
        sliderInput("neuralmaxiterations", label="Max Iterations", min=50, max=2000, value=selection)
    } else if(radiocal==7 && neuralhiddenlayers > 1){
        NULL
    }
}


###Neural Net Plots
plot.nnet<-function(mod.in,nid=T,all.out=T,all.in=T,bias=T,wts.only=F,rel.rsc=5,
circle.cex=5,node.labs=T,var.labs=T,x.lab=NULL,y.lab=NULL,
line.stag=NULL,struct=NULL,cex.val=1,alpha.val=1,
circle.col='lightblue',pos.col='black',neg.col='grey',
bord.col='lightblue', max.sp = F,...){
    
    require(scales)
    
    #sanity checks
    if('mlp' %in% class(mod.in)) warning('Bias layer not applicable for rsnns object')
    if('numeric' %in% class(mod.in)){
        if(is.null(struct)) stop('Three-element vector required for struct')
        if(length(mod.in) != ((struct[1]*struct[2]+struct[2]*struct[3])+(struct[3]+struct[2])))
        stop('Incorrect length of weight matrix for given network structure')
    }
    if('train' %in% class(mod.in)){
        if('nnet' %in% class(mod.in$finalModel)){
            mod.in<-mod.in$finalModel
            warning('Using best nnet model from train output')
        } else if('nn' %in% class(mod.in$finalModel)){
            mod.o <- mod.in
            mod.in<-mod.in$finalModel
            warning('Using best nn model from train output')
        }
        else stop('Only nnet method can be used with train object')
    }
    
    #gets weights for neural network, output is list
    #if rescaled argument is true, weights are returned but rescaled based on abs value
    nnet.vals<-function(mod.in,nid,rel.rsc,struct.out=struct){
        
        require(scales)
        require(reshape)
        
        if('numeric' %in% class(mod.in)){
            struct.out<-struct
            wts<-mod.in
        }
        
        #neuralnet package
        if('nn' %in% class(mod.in)){
            struct.out<-unlist(lapply(mod.in$weights[[1]],ncol))
            struct.out<-struct.out[-length(struct.out)]
            struct.out<-c(
            length(mod.in$model.list$variables),
            struct.out,
            length(mod.in$model.list$response)
            )
            wts<-unlist(mod.in$weights[[1]])
        }
        
        #nnet package
        if('nnet' %in% class(mod.in)){
            struct.out<-mod.in$n
            wts<-mod.in$wts
        }
        
        #RSNNS package
        if('mlp' %in% class(mod.in)){
            struct.out<-c(mod.in$nInputs,mod.in$archParams$size,mod.in$nOutputs)
            hid.num<-length(struct.out)-2
            wts<-mod.in$snnsObject$getCompleteWeightMatrix()
            
            #get all input-hidden and hidden-hidden wts
            inps<-wts[grep('Input',row.names(wts)),grep('Hidden_2',colnames(wts)),drop=F]
            inps<-melt(rbind(rep(NA,ncol(inps)),inps))$value
            uni.hids<-paste0('Hidden_',1+seq(1,hid.num))
            for(i in 1:length(uni.hids)){
                if(is.na(uni.hids[i+1])) break
                tmp<-wts[grep(uni.hids[i],rownames(wts)),grep(uni.hids[i+1],colnames(wts)),drop=F]
                inps<-c(inps,melt(rbind(rep(NA,ncol(tmp)),tmp))$value)
            }
            
            #get connections from last hidden to output layers
            outs<-wts[grep(paste0('Hidden_',hid.num+1),row.names(wts)),grep('Output',colnames(wts)),drop=F]
            outs<-rbind(rep(NA,ncol(outs)),outs)
            
            #weight vector for all
            wts<-c(inps,melt(outs)$value)
            assign('bias',F,envir=environment(nnet.vals))
        }
        
        if(nid) wts<-rescale(abs(wts),c(1,rel.rsc))
        
        #convert wts to list with appropriate names
        hid.struct<-struct.out[-c(length(struct.out))]
        row.nms<-NULL
        for(i in 1:length(hid.struct)){
            if(is.na(hid.struct[i+1])) break
            row.nms<-c(row.nms,rep(paste('hidden',i,seq(1:hid.struct[i+1])),each=1+hid.struct[i]))
        }
        row.nms<-c(
        row.nms,
        rep(paste('out',seq(1:struct.out[length(struct.out)])),each=1+struct.out[length(struct.out)-1])
        )
        out.ls<-data.frame(wts,row.nms, stringsAsFactors=FALSE)
        out.ls$row.nms<-factor(row.nms,levels=unique(row.nms),labels=unique(row.nms))
        out.ls<-split(out.ls$wts,f=out.ls$row.nms)
        
        assign('struct',struct.out,envir=environment(nnet.vals))
        
        out.ls
        
    }
    
    wts<-nnet.vals(mod.in,nid=F)
    
    if(wts.only) return(wts)
    
    #circle colors for input, if desired, must be two-vector list, first vector is for input layer
    if(is.list(circle.col)){
        circle.col.inp<-circle.col[[1]]
        circle.col<-circle.col[[2]]
    } else circle.col.inp<-circle.col
    
    #initiate plotting
    x.range<-c(0,100)
    y.range<-c(0,100)
    #these are all proportions from 0-1
    if(is.null(line.stag)) line.stag<-0.011*circle.cex/2
    layer.x<-seq(0.17,0.9,length=length(struct))
    bias.x<-layer.x[-length(layer.x)]+diff(layer.x)/2
    bias.y<-0.95
    circle.cex<-circle.cex
    
    #get variable names from mod.in object
    #change to user input if supplied
    if('numeric' %in% class(mod.in)){
        x.names<-paste0(rep('X',struct[1]),seq(1:struct[1]))
        y.names<-paste0(rep('Y',struct[3]),seq(1:struct[3]))
    }
    if('mlp' %in% class(mod.in)){
        all.names<-mod.in$snnsObject$getUnitDefinitions()
        x.names<-all.names[grep('Input',all.names$unitName),'unitName']
        y.names<-all.names[grep('Output',all.names$unitName),'unitName']
    }
    if('nn' %in% class(mod.in)){
        x.names<-mod.in$model.list$variables
        y.names<-mod.in$model.list$respons
    }
    if('xNames' %in% names(mod.in)){
        x.names<-mod.in$xNames
        y.names<-if('nn' %in% class(mod.in)){
            attr(terms(mod.o),'factor')
        } else {
            attr(terms(mod.in),'factor')
        }
        
        y.names<-row.names(y.names)[!row.names(y.names) %in% x.names]
    }
    if(!'xNames' %in% names(mod.in) & 'nnet' %in% class(mod.in)){
        if(is.null(mod.in$call$formula)){
            x.names<-colnames(eval(mod.in$call$x))
            y.names<-colnames(eval(mod.in$call$y))
        }
        else{
            forms<-eval(mod.in$call$formula)
            x.names<-mod.in$coefnames
            facts<-attr(terms(mod.in),'factors')
            y.check<-mod.in$fitted
            if(ncol(y.check)>1) y.names<-colnames(y.check)
            else y.names<-as.character(forms)[2]
        }
    }
    #change variables names to user sub
    if(!is.null(x.lab)){
        if(length(x.names) != length(x.lab)) stop('x.lab length not equal to number of input variables')
        else x.names<-x.lab
    }
    if(!is.null(y.lab)){
        if(length(y.names) != length(y.lab)) stop('y.lab length not equal to number of output variables')
        else y.names<-y.lab
    }
    
    #initiate plot
    plot(x.range,y.range,type='n',axes=F,ylab='',xlab='',...)
    
    #function for getting y locations for input, hidden, output layers
    #input is integer value from 'struct'
    get.ys<-function(lyr, max_space = max.sp){
        if(max_space){
            spacing <- diff(c(0*diff(y.range),0.9*diff(y.range)))/lyr
        } else {
            spacing<-diff(c(0*diff(y.range),0.9*diff(y.range)))/max(struct)
        }
        
        seq(0.5*(diff(y.range)+spacing*(lyr-1)),0.5*(diff(y.range)-spacing*(lyr-1)),
        length=lyr)
    }
    
    #function for plotting nodes
    #'layer' specifies which layer, integer from 'struct'
    #'x.loc' indicates x location for layer, integer from 'layer.x'
    #'layer.name' is string indicating text to put in node
    layer.points<-function(layer,x.loc,layer.name,cex=cex.val){
        x<-rep(x.loc*diff(x.range),layer)
        y<-get.ys(layer)
        points(x,y,pch=21,cex=circle.cex,col=bord.col,bg=in.col)
        if(node.labs) text(x,y,paste(layer.name,1:layer,sep=''),cex=cex.val)
        if(layer.name=='I' & var.labs) text(x-line.stag*diff(x.range),y,x.names,pos=2,cex=cex.val)
        if(layer.name=='O' & var.labs) text(x+line.stag*diff(x.range),y,y.names,pos=4,cex=cex.val)
    }
    
    #function for plotting bias points
    #'bias.x' is vector of values for x locations
    #'bias.y' is vector for y location
    #'layer.name' is  string indicating text to put in node
    bias.points<-function(bias.x,bias.y,layer.name,cex,...){
        for(val in 1:length(bias.x)){
            points(
            diff(x.range)*bias.x[val],
            bias.y*diff(y.range),
            pch=21,col=bord.col,bg=in.col,cex=circle.cex
            )
            if(node.labs)
            text(
            diff(x.range)*bias.x[val],
            bias.y*diff(y.range),
            paste(layer.name,val,sep=''),
            cex=cex.val
            )
        }
    }
    
    #function creates lines colored by direction and width as proportion of magnitude
    #use 'all.in' argument if you want to plot connection lines for only a single input node
    layer.lines<-function(mod.in,h.layer,layer1=1,layer2=2,out.layer=F,nid,rel.rsc,all.in,pos.col,
    neg.col,...){
        
        x0<-rep(layer.x[layer1]*diff(x.range)+line.stag*diff(x.range),struct[layer1])
        x1<-rep(layer.x[layer2]*diff(x.range)-line.stag*diff(x.range),struct[layer1])
        
        if(out.layer==T){
            
            y0<-get.ys(struct[layer1])
            y1<-rep(get.ys(struct[layer2])[h.layer],struct[layer1])
            src.str<-paste('out',h.layer)
            
            wts<-nnet.vals(mod.in,nid=F,rel.rsc)
            wts<-wts[grep(src.str,names(wts))][[1]][-1]
            wts.rs<-nnet.vals(mod.in,nid=T,rel.rsc)
            wts.rs<-wts.rs[grep(src.str,names(wts.rs))][[1]][-1]
            
            cols<-rep(pos.col,struct[layer1])
            cols[wts<0]<-neg.col
            
            if(nid) segments(x0,y0,x1,y1,col=cols,lwd=wts.rs)
            else segments(x0,y0,x1,y1)
            
        }
        
        else{
            
            if(is.logical(all.in)) all.in<-h.layer
            else all.in<-which(x.names==all.in)
            
            y0<-rep(get.ys(struct[layer1])[all.in],struct[2])
            y1<-get.ys(struct[layer2])
            src.str<-paste('hidden',layer1)
            
            wts<-nnet.vals(mod.in,nid=F,rel.rsc)
            wts<-unlist(lapply(wts[grep(src.str,names(wts))],function(x) x[all.in+1]))
            wts.rs<-nnet.vals(mod.in,nid=T,rel.rsc)
            wts.rs<-unlist(lapply(wts.rs[grep(src.str,names(wts.rs))],function(x) x[all.in+1]))
            
            cols<-rep(pos.col,struct[layer2])
            cols[wts<0]<-neg.col
            
            if(nid) segments(x0,y0,x1,y1,col=cols,lwd=wts.rs)
            else segments(x0,y0,x1,y1)
            
        }
        
    }
    
    bias.lines<-function(bias.x,mod.in,nid,rel.rsc,all.out,pos.col,neg.col,...){
        
        if(is.logical(all.out)) all.out<-1:struct[length(struct)]
        else all.out<-which(y.names==all.out)
        
        for(val in 1:length(bias.x)){
            
            wts<-nnet.vals(mod.in,nid=F,rel.rsc)
            wts.rs<-nnet.vals(mod.in,nid=T,rel.rsc)
            
            if(val != length(bias.x)){
                wts<-wts[grep('out',names(wts),invert=T)]
                wts.rs<-wts.rs[grep('out',names(wts.rs),invert=T)]
                sel.val<-grep(val,substr(names(wts.rs),8,8))
                wts<-wts[sel.val]
                wts.rs<-wts.rs[sel.val]
            }
            
            else{
                wts<-wts[grep('out',names(wts))]
                wts.rs<-wts.rs[grep('out',names(wts.rs))]
            }
            
            cols<-rep(pos.col,length(wts))
            cols[unlist(lapply(wts,function(x) x[1]))<0]<-neg.col
            wts.rs<-unlist(lapply(wts.rs,function(x) x[1]))
            
            if(nid==F){
                wts.rs<-rep(1,struct[val+1])
                cols<-rep('black',struct[val+1])
            }
            
            if(val != length(bias.x)){
                segments(
                rep(diff(x.range)*bias.x[val]+diff(x.range)*line.stag,struct[val+1]),
                rep(bias.y*diff(y.range),struct[val+1]),
                rep(diff(x.range)*layer.x[val+1]-diff(x.range)*line.stag,struct[val+1]),
                get.ys(struct[val+1]),
                lwd=wts.rs,
                col=cols
                )
            }
            
            else{
                segments(
                rep(diff(x.range)*bias.x[val]+diff(x.range)*line.stag,struct[val+1]),
                rep(bias.y*diff(y.range),struct[val+1]),
                rep(diff(x.range)*layer.x[val+1]-diff(x.range)*line.stag,struct[val+1]),
                get.ys(struct[val+1])[all.out],
                lwd=wts.rs[all.out],
                col=cols[all.out]
                )
            }
            
        }
    }
    
    #use functions to plot connections between layers
    #bias lines
    if(bias) bias.lines(bias.x,mod.in,nid=nid,rel.rsc=rel.rsc,all.out=all.out,pos.col=alpha(pos.col,alpha.val),
    neg.col=alpha(neg.col,alpha.val))
    
    #layer lines, makes use of arguments to plot all or for individual layers
    #starts with input-hidden
    #uses 'all.in' argument to plot connection lines for all input nodes or a single node
    if(is.logical(all.in)){
        mapply(
        function(x) layer.lines(mod.in,x,layer1=1,layer2=2,nid=nid,rel.rsc=rel.rsc,
        all.in=all.in,pos.col=alpha(pos.col,alpha.val),neg.col=alpha(neg.col,alpha.val)),
        1:struct[1]
        )
    }
    else{
        node.in<-which(x.names==all.in)
        layer.lines(mod.in,node.in,layer1=1,layer2=2,nid=nid,rel.rsc=rel.rsc,all.in=all.in,
        pos.col=alpha(pos.col,alpha.val),neg.col=alpha(neg.col,alpha.val))
    }
    #connections between hidden layers
    lays<-split(c(1,rep(2:(length(struct)-1),each=2),length(struct)),
    f=rep(1:(length(struct)-1),each=2))
    lays<-lays[-c(1,(length(struct)-1))]
    for(lay in lays){
        for(node in 1:struct[lay[1]]){
            layer.lines(mod.in,node,layer1=lay[1],layer2=lay[2],nid=nid,rel.rsc=rel.rsc,all.in=T,
            pos.col=alpha(pos.col,alpha.val),neg.col=alpha(neg.col,alpha.val))
        }
    }
    #lines for hidden-output
    #uses 'all.out' argument to plot connection lines for all output nodes or a single node
    if(is.logical(all.out))
    mapply(
    function(x) layer.lines(mod.in,x,layer1=length(struct)-1,layer2=length(struct),out.layer=T,nid=nid,rel.rsc=rel.rsc,
    all.in=all.in,pos.col=alpha(pos.col,alpha.val),neg.col=alpha(neg.col,alpha.val)),
    1:struct[length(struct)]
    )
    else{
        node.in<-which(y.names==all.out)
        layer.lines(mod.in,node.in,layer1=length(struct)-1,layer2=length(struct),out.layer=T,nid=nid,rel.rsc=rel.rsc,
        pos.col=pos.col,neg.col=neg.col,all.out=all.out)
    }
    
    #use functions to plot nodes
    for(i in 1:length(struct)){
        in.col<-circle.col
        layer.name<-'H'
        if(i==1) { layer.name<-'I'; in.col<-circle.col.inp}
        if(i==length(struct)) layer.name<-'O'
        layer.points(struct[i],layer.x[i],layer.name)
    }
    
    if(bias) bias.points(bias.x,bias.y,'B')
    
}
plot.nnet <- cmpfun(plot.nnet)








