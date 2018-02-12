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


dropdownButton(
tags$h3("Manual Changes"), icon = icon("gear"),
checkboxInput('advanced', "Advanced", value=FALSE),
checkboxInput('backgroundsubtract', "Background Subtract", value=FALSE),
checkboxInput('combine', "Combine", value=FALSE),
checkboxInput('invert', "Invert", value=FALSE),
tooltip = tooltipOptions(title = "Click for manual options")
),




uiOutput('gainshiftui'),
             
             tags$hr(),

actionButton('actionprocess', label = "Process Data"),
actionButton('actionplot', label = "Plot Spectrum"),
downloadButton('downloadPlot', "Plot"),


tags$hr(),

uiOutput('filegrab'),

selectInput('filetype', label="Filetype", c("DPT", "CSV"), selected="DPT")
),


mainPanel(

div(
style = "position:relative",
plotOutput('distPlot', height = 685,
dblclick = 'plot1_dblclick',
brush = brushOpts(id = 'plot1_brush', resetOnNew = TRUE),
hover = hoverOpts('plot_hover_spectrum', delay = 100, delayType = "debounce")),
uiOutput('hover_info_spectrum'))
)
))
)


))
)
