library(shiny)
library(DT)
library(dplyr)
library(shinythemes)
library(data.table)
library(dtplyr)
library(rhandsontable)
library(Cairo)
library(shinyWidgets)

options(warn=-1)
assign("last.warning", NULL, envir = baseenv())




ui=list(
tagList(
header=tags$head(tags$style(".table .alignRight {color: black; text-align:right;}"))),

tags$head(tags$script("$(function() {$.fn.dataTableExt.errMode = 'throw';});")),

shinyUI(navbarPage("cloudFTIR", id="nav", theme = shinytheme("flatly"),



tabPanel("Spectrum",
div(class="outer",
headerPanel("FTIR Spectrum Viewer"),
sidebarLayout(
sidebarPanel(

tags$style(type="text/css",
".shiny-output-error { visibility: hidden; }",
".shiny-output-error:before { visibility: hidden; }"
),


textInput('projectname', label = "Project Name", value="myFTIR"),

tags$hr(),

#actionButton('actionprocess', label = "Process Data"),
downloadButton('downloadPlot', "Plot"),
downloadButton('downloadPeakTableID', "Peak ID"),
downloadButton('downloadPeakTableSummary', "Summary"),
downloadButton('downloadsummaryplot', "Summary Plot"),
downloadButton('downloadPeakTable', "Raw"),
#downloadButton('downloadfp', "FP"),



dropdownButton(
tags$h3("Manual Changes"), icon = icon("gear"),
checkboxInput('advanced', "Advanced", value=FALSE),
uiOutput('gainshiftui'),
checkboxInput('backgroundsubtract', "Background Subtract", value=FALSE),
checkboxInput('combine', "Combine", value=FALSE),
checkboxInput('invert', "Invert", value=FALSE),
tooltip = tooltipOptions(title = "Click for manual options")
),

checkboxInput('showpeaks', "Show Peaks", value=TRUE),


sliderInput('spikesensitivity', "Spike Proximity", min=0.1, max=100, value=20),
#uiOutput('uispikeheight'),
sliderInput('spikeheight', "Spike Height", min=0, max=1, value=.2),
sliderInput('wavethreshold', "Accepted Wavenumber Buffer", min=0, max=20, value=3),






tags$hr(),

uiOutput('filegrab'),

selectInput('filetype', label="Filetype", c("DPT", "CSV", "Opus"), selected="DPT")
),


mainPanel(

tags$style(type="text/css",
".shiny-output-error { visibility: hidden; }",
".shiny-output-error:before { visibility: hidden; }"
),
tabsetPanel(
tabPanel("Plot",
    div(
        style = "position:relative",
        plotOutput('distPlot', height = 685, click='plot1_click',
            dblclick = 'plot1_dblclick',
            brush = brushOpts(id = 'plot1_brush', resetOnNew = TRUE),
            hover = hoverOpts('plot_hover_spectrum', delay = 100, delayType = "debounce")),
        uiOutput('hover_info_spectrum')),
actionButton("exclude_toggle", "Toggle points"),
actionButton("exclude_reset", "Reset")
),
tabPanel("Peak ID", dataTableOutput('peaktableid')),
tabPanel("Summary", dataTableOutput('peaktablesummary')),
tabPanel("Summary Plot",
    div(
        style = "position:relative",
        plotOutput('summaryplot', height = 685,
            dblclick = 'sumplot1_dblclick',
            brush = brushOpts(id = 'sumplot1_brush', resetOnNew = TRUE),
        hover = hoverOpts('sumplot_hover_spectrum', delay = 100, delayType = "debounce"))

)),
tabPanel("Raw Data", dataTableOutput('peaktable'))

))
)
)),



tabPanel("Add Wavenumbers",
div(class="outer",

fluidRow(
sidebarLayout(
sidebarPanel(

actionButton('resethotablewave', "Reset"),
actionButton('linecommitwave', "Add Wavenumbers")



),


mainPanel(
tabsetPanel(
tabPanel('Enter Concentrations', rHandsontableOutput('hotwave')),
tabPanel('Covariance', plotOutput('covarianceplotvalueswave'))
))
))
)),

tabPanel("Add Concentrations",
div(class="outer",

fluidRow(
sidebarLayout(
sidebarPanel(

actionButton('resethotable', "Reset"),
actionButton('hotableprocess2', "Enter Values"),
tags$hr(),
textInput("calunits", label = "Units", value="Weight %")


),


mainPanel(
tabsetPanel(
id = 'dataset',
tabPanel('Enter Concentrations', rHandsontableOutput('hot')),
tabPanel('Covariance', plotOutput('covarianceplotvalues'))

))
))
))

))
)
