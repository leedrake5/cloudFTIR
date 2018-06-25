list.of.packages <- c("pbapply", "reshape2", "TTR", "dplyr", "ggtern", "ggplot2", "shiny", "rhandsontable", "random", "data.table", "DT", "shinythemes", "Cairo", "broom", "shinyjs", "gridExtra", "dtplyr", "formattable", "XML", "corrplot", "scales", "rmarkdown", "markdown", "peakPick", "shinyWidgets", "soil.spec", "data.table")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
#if(length(new.packages)) install.packages(new.packages)

options(warn=-1)
assign("last.warning", NULL, envir = baseenv())

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

my.cores <- if(parallel::detectCores()>=3){
                paste0(parallel::detectCores()-2)
            } else if(parallel::detectCores()<=2){
                "1"
            }

options(digits=4)
options(warn=-1)
assign("last.warning", NULL, envir = baseenv())

keep_singles <- function(v){
    v[!(v %in% v[duplicated(v)])]
}

Hodder.v <- function(y)
{
    
    n<-length(y)
    
    for(i in 1:(n-1)) {
        y[i] <- y[i+1] - y[i]
        y[1:(n-1)]
        y <- abs(y)
    }
    y <- c(0, y[1:(n-1)])
    
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

range_subset <- function(range.frame, data){
    
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
    
    selected.list <- lapply(index, function(x) range_subset(range.frame=choice.list[[x]], data=data))
    
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
    
    l <- list(a = format(coef(m)[1], digits = 2),
    b = format(abs(coef(m)[2]), digits = 2),
    r2 = format(summary(m)$r.squared, digits = 3));
    
        eq <- substitute(italic(C)[i] == a + b %.% italic(I)[i]*","~~italic(r)^2~"="~r2,l)
  
    
    as.character(as.expression(eq));
}

lm_eqn_poly = function(m) {
    
    l <- list(a = format(coef(m)[1], digits = 2),
    b = format(abs(coef(m)[2]), digits = 2),
    c = format(abs(coef(m)[3]), digits = 2),
    r2 = format(summary(m)$r.squared, digits = 3));
    
        eq <- substitute(italic(C)[i] == a + c %.% italic(I)[i]^2 + b %.% italic(I)[i]*","~~italic(r)^2~"="~r2,l)
   
    
    as.character(as.expression(eq));
}

lm_eqn_val = function(m) {
    
    l <- list(a = format(coef(m)[1], digits = 2),
    b = format(abs(coef(m)[2]), digits = 2),
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

pblapply <- function (X, FUN, ..., cl = my.cores)
{
    FUN <- match.fun(FUN)
    if (!is.vector(X) || is.object(X))
    X <- as.list(X)
    if (!is.null(cl)) {
        if (.Platform$OS.type == "windows") {
            if (!inherits(cl, "cluster"))
            cl <- NULL
        }
        else {
            if (inherits(cl, "cluster")) {
                if (length(cl) < 2L)
                cl <- NULL
            }
            else {
                if (cl < 2)
                cl <- NULL
            }
        }
    }
    nout <- as.integer(getOption("pboptions")$nout)
    if (is.null(cl)) {
        if (!dopb())
        return(lapply(X, FUN, ...))
        Split <- splitpb(length(X), 1L, nout = nout)
        B <- length(Split)
        pb <- startpb(0, B)
        on.exit(closepb(pb), add = TRUE)
        rval <- vector("list", B)
        for (i in seq_len(B)) {
            rval[i] <- list(lapply(X[Split[[i]]], FUN, ...))
            setpb(pb, i)
        }
    }
    else {
        if (inherits(cl, "cluster")) {
            PAR_FUN <- if (isTRUE(getOption("pboptions")$use_lb))
            parallel::parLapplyLB
            else parallel::parLapply
            if (!dopb())
            return(PAR_FUN(my.cores, X, FUN, ...))
            Split <- splitpb(length(X), length(my.cores), nout = nout)
            B <- length(Split)
            pb <- startpb(0, B)
            on.exit(closepb(pb), add = TRUE)
            rval <- vector("list", B)
            for (i in seq_len(B)) {
                rval[i] <- list(PAR_FUN(cl, X[Split[[i]]], FUN,
                ...))
                setpb(pb, i)
            }
        }
        else {
            if (!dopb())
            return(parallel::mclapply(X, FUN, ..., mc.cores = as.integer(my.cores)))
            Split <- splitpb(length(X), as.integer(my.cores), nout = nout)
            B <- length(Split)
            pb <- startpb(0, B)
            on.exit(closepb(pb), add = TRUE)
            rval <- vector("list", B)
            for (i in seq_len(B)) {
                rval[i] <- list(parallel::mclapply(X[Split[[i]]],
                FUN, ..., mc.cores = as.integer(my.cores)))
                setpb(pb, i)
            }
        }
    }
    rval <- do.call(c, rval, quote = TRUE)
    names(rval) <- names(X)
    rval
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




spectra_simp_prep_ftir <- function(spectra){
    
    spectra$Energy <- round(spectra$Wavenumber, 0)
    spectra <- data.table(spectra)
    spectra.aggregate <- spectra[, list(Amplitude=mean(Amplitude, na.rm = TRUE)), by = list(Spectrum,Wavenumber)]
    
    data <- as.data.frame(dcast.data.table(spectra.aggregate, Spectrum~Wavenumber, value.var="Amplitude"))
    
    #test <- apply(test, 2, as.numeric)
    colnames(data) <- make.names(colnames(data))
    do.call(data.frame,lapply(data, function(x) replace(x, is.infinite(x),0)))
    
}

spectra_tc_prep_ftir <- function(spectra){
    
    spectra$Energy <- round(spectra$Wavenumber, 0)
    
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

spectra_comp_prep_ftir <- function(spectra, norm.min, norm.max){
    
    compton.norm <- subset(spectra$Amplitude, !(spectra$Wavenumber < norm.min | spectra$Wavenumber > norm.max))
    compton.file <- subset(spectra$Spectrum, !(spectra$Wavenumber < norm.min | spectra$Wavenumber > norm.max))
    compton.frame <- data.frame(is.0(compton.norm, compton.file))
    colnames(compton.frame) <- c("Compton", "Spectrum")
    compton.frame.ag <- aggregate(list(compton.frame$Compton), by=list(compton.frame$Spectrum), FUN="sum")
    colnames(compton.frame.ag) <- c("Spectrum", "Compton")
    
    spectra$Wavenumber <- round(spectra$Wavenumber, 0)
    
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











