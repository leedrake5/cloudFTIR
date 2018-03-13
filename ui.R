library(shiny)
library(DT)
library(dplyr)
library(shinythemes)
library(data.table)
library(dtplyr)
library(rhandsontable)
library(Cairo)
library(shinyWidgets)


ui=list(
tagList(
header=tags$head(tags$style(".table .alignRight {color: black; text-align:right;}"))),


shinyUI(navbarPage("cloudFTIR", id="nav", theme = shinytheme("flatly"),

tabPanel("Spectrum",
div(class="outer",
headerPanel("FTIR Spectrum Viewer"),
sidebarLayout(
sidebarPanel(


textInput('projectname', label = "Project Name", value="myFTIR"),

tags$hr(),

actionButton('actionprocess', label = "Process Data"),
downloadButton('downloadPlot', "Plot"),
downloadButton('downloadPeakTable', "Table"),




dropdownButton(
tags$h3("Manual Changes"), icon = icon("gear"),
checkboxInput('advanced', "Advanced", value=FALSE),
uiOutput('gainshiftui'),
checkboxInput('backgroundsubtract', "Background Subtract", value=FALSE),
checkboxInput('combine', "Combine", value=FALSE),
checkboxInput('invert', "Invert", value=FALSE),
tooltip = tooltipOptions(title = "Click for manual options")
),

dropdownButton(
tags$h3("Analytics"), icon = icon("code"),
checkboxInput('showpeaks', "Show Peaks", value=FALSE),
sliderInput('spikesensitivity', "Spike Sensitivity", min=0.1, max=100, value=20),
uiOutput('uispikeheight'),
tooltip = tooltipOptions(title = "Click for analytics")
),





tags$hr(),

uiOutput('filegrab'),

selectInput('filetype', label="Filetype", c("DPT", "CSV", "Opus"), selected="DPT")
),


mainPanel(
tabsetPanel(
tabPanel("Plot",

div(
style = "position:relative",
plotOutput('distPlot', height = 685,
dblclick = 'plot1_dblclick',
brush = brushOpts(id = 'plot1_brush', resetOnNew = TRUE),
hover = hoverOpts('plot_hover_spectrum', delay = 100, delayType = "debounce")),
uiOutput('hover_info_spectrum'))
),
tabPanel("Table", dataTableOutput('peaktable'))
))
)
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
