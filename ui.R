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

shinyUI(navbarPage("cloudFTIR", id="nav", theme = shinytheme("yeti"),



tabPanel("Spectrum",
div(class="outer",
headerPanel("FTIR Spectrum Viewer"),
sidebarLayout(
sidebarPanel(

#tags$style(type="text/css",
#".shiny-output-error { visibility: hidden; }",
#".shiny-output-error:before { visibility: hidden; }"
#),


textInput('projectname', label = "Project Name", value="myFTIR"),

tags$hr(),

#actionButton('actionprocess', label = "Process Data"),
#downloadButton('downloadfp', "FP"),

checkboxInput('showspectralplot', "Show Spectra", value=TRUE),

tags$hr(),

dropdownButton(
tags$h3("Manual Changes"), icon = icon("stats", lib="glyphicon"),
uiOutput('datatransformationsui'),
checkboxInput('advanced', "Advanced", value=FALSE),
uiOutput('gainshiftui'),
checkboxInput('combine', "Combine", value=FALSE),
checkboxInput('invert', "Invert", value=FALSE),
tooltip = tooltipOptions(title = "Click for transformation options")
),

tags$hr(),

dropdownButton(
tags$h3("Peaks"), icon = icon("screenshot", lib="glyphicon"),
checkboxInput('showpeaks', "Show Peaks", value=FALSE),
sliderInput('spikesensitivity', "Spike Proximity", min=0.1, max=100, value=20),
#uiOutput('uispikeheight'),
sliderInput('spikeheight', "Spike Height", min=0, max=1, value=.2),
sliderInput('wavethreshold', "Accepted Wavenumber Buffer", min=0, max=20, value=3),
tooltip = tooltipOptions(title = "Click for peak options")
),



tags$hr(),

uiOutput('filegrab'),

selectInput('filetype', label="Filetype", c("DPT", "CSV", "Opus"), selected="Opus"),

tags$hr(),

fileInput('calfileinput', 'Load Cal File', accept=".quant", multiple=FALSE),
checkboxInput('usecalfile', "Use Cal File")
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
tags$hr(),
actionButton("exclude_toggle", "Toggle points"),
actionButton("exclude_reset", "Reset"),
downloadButton('downloadPlot', "Plot"),
checkboxInput('showlegend', "Show Legend", value=TRUE)
),
tabPanel("Peak ID", dataTableOutput('peaktableid'),
    tags$hr(),
    downloadButton('downloadPeakTableID', "Peak ID")
),
tabPanel("Summary", dataTableOutput('peaktablesummary'),
    tags$hr(),
    downloadButton('downloadPeakTableSummary', "Summary")),
tabPanel("Summary Plot",
    div(
        style = "position:relative",
        plotOutput('summaryplot', height = 685,
            dblclick = 'sumplot1_dblclick',
            brush = brushOpts(id = 'sumplot1_brush', resetOnNew = TRUE),
        hover = hoverOpts('sumplot_hover_spectrum', delay = 100, delayType = "debounce"))
    ),
        tags$hr(),
        downloadButton('downloadsummaryplot', "Summary Plot")),
tabPanel("Raw Data", dataTableOutput('peaktable'),
    downloadButton('downloadData', "Raw")

    )

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
tabPanel('Line Counts', dataTableOutput('LineValues')),
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
tabPanel('Enter Concentrations', rHandsontableOutput('hot'),
tags$script(
'
setTimeout(
function() {
    HTMLWidgets.find("#hot").hot.addHook(
    "afterColumnSort",
    function(){
        console.log("sort",this);
        Shiny.onInputChange(
        "hot_sort",
        {
            data: this.getData()
        }
        )
    }
    )
},
1000
)'
)),
tabPanel('Covariance', plotOutput('covarianceplotvalues'))

))
))
)),


tabPanel("Cal Curves",
div(class="outer",

fluidRow(
sidebarLayout(
sidebarPanel(width=3,


actionButton('createcalelement', "Update"),
actionButton('createcal', "Save"),
tags$hr(),
downloadButton('downloadModel', "Model"),
downloadButton('downloadReport', "Report"),


tags$hr(),

actionButton('trainslopes', "Train"),

tags$hr(),

#uiOutput('testing'),


uiOutput('inVar2'),

uiOutput('calTypeInput'),

uiOutput('normTypeInput'),

uiOutput('comptonMinInput'),

uiOutput('comptonMaxInput'),

uiOutput('inVar3'),
uiOutput('inVar4')
#sliderInput("nvariables", label = "# Elements", min=1, max=7, value=2)

),

mainPanel(
tabsetPanel(
#tabPanel("testing", dataTableOutput('testtable2')),
tabPanel("Cal Curves",
splitLayout(cellWidths = c("50%", "50%"),
column(width=12,
div(
style = "position:relative",
plotOutput("calcurveplots", height = 455, click = "plot_cal_click",
dblclick = dblclickOpts(id="plot_cal_dblclick"),
brush = brushOpts(id = "plot_cal_brush", resetOnNew = TRUE),
hover = hoverOpts("plot_hovercal", delay = 100, delayType = "debounce")),
uiOutput("hover_infocal")),
actionButton("cropcal", "Zoom")),
column(width=12,
div(
style = "position:relative",
plotOutput("valcurveplots", height = 455, click = "plot_val_click",
dblclick = "plot_val_dblclick",
brush = brushOpts(id = "plot_val_brush", resetOnNew = TRUE),
hover = hoverOpts("plot_hoverval", delay = 100, delayType = "debounce")),
uiOutput("hover_infoval")),
actionButton("cropval", "Zoom")
)
),
tags$hr(),
actionButton("exclude_toggle", "Toggle points"),
actionButton("exclude_reset", "Reset"),
downloadButton('downloadcloudplot', "Plot")

),

tabPanel("Cross Validation",
splitLayout(cellWidths = c("50%", "50%"),
column(width=12,
div(
style = "position:relative",
plotOutput("calcurveplotsrandom", height = 455,  click = "plot_cal_click_random",
dblclick = "plot_cal_dblclick_random",
brush = brushOpts(id = "plot_cal_brush_random", resetOnNew = TRUE),
hover = hoverOpts("plot_hovercal_random", delay = 100, delayType = "debounce")),
uiOutput("hover_infocal_random")),
actionButton("cropcalrandom", "Zoom")
),
column(width=12,
div(
style = "position:relative",
plotOutput("valcurveplotsrandom", height = 455, click = "plot_val_click_random",
dblclick = "plot_val_dblclick_random",
brush = brushOpts(id = "plot_val_brush_random", resetOnNew = TRUE),
hover = hoverOpts("plot_hoverval_random", delay = 100, delayType = "debounce")),
uiOutput("hover_infoval_random")),
actionButton("cropvalrandom", "Zoom")
)),
tags$hr(),
sliderInput('percentrandom', "Randomize", min=.01, max=.99, value=.33),
downloadButton('downloadcloudplotrandomized', "Plot")


),

tabPanel("Models", dataTableOutput("models")),


tabPanel("Diagnostics",
splitLayout(cellWidths = c("50%", "50%"),
div(
style = "position:relative",
plotOutput("residualsfitted", height=250, click="plot_residualsfitted_click", brush=brushOpts(id="plot_residualsfitted_brush", resetOnNew = TRUE),
hover = hoverOpts("plot_hoverresidualsfitted", delay = 100, delayType = "debounce")),
uiOutput("hover_inforesidualsfitted")),
div(
style = "position:relative",
plotOutput("qq", height=250, click="plot_qq_click",
brush=brushOpts(id="plot_qq_brush", resetOnNew = TRUE),
hover = hoverOpts("plot_hoverqq", delay = 100, delayType = "debounce")),
uiOutput("hover_infoqq"))
),
splitLayout(cellWidths = c("50%", "50%"),
div(
style = "position:relative",
plotOutput("scalelocation", height=250, click="plot_scalelocation_click", brush=brushOpts(id="plot_scalelocation_brush", resetOnNew = TRUE),
hover = hoverOpts("plot_hoverscalelocation", delay = 100, delayType = "debounce")),
uiOutput("hover_infoscalelocation")),
plotOutput("cooksdistance", height=250, click="plot_cooksdistance_click", brush=brushOpts(id="plot_cooksdistance_brush", resetOnNew = TRUE))
),
splitLayout(cellWidths = c("50%", "50%"),
div(
style = "position:relative",
plotOutput("residualleverage", height=250, click="plot_residualleverage_click", brush=brushOpts(id="plot_residualleverage_brush", resetOnNew = TRUE),
hover = hoverOpts("plot_hoverresidualleverage", delay = 100, delayType = "debounce")),
uiOutput("hover_inforesidualleverage")),
div(
style = "position:relative",
plotOutput("cooksleverage", height=250, click="plot_cooksleverage_click", brush=brushOpts(id="plot_cooksleverage_brush", resetOnNew = TRUE),
hover = hoverOpts("plot_hovercooksleverage", delay = 100, delayType = "debounce")),
uiOutput("hover_infocooksleverage"))
),
actionButton("exclude_toggle_diag", "Toggle points"),
actionButton("exclude_reset_diag", "Reset")),

tabPanel("Variables",
div(
style = "position:relative",
plotOutput('importanceplot', hover = hoverOpts('plot_hover_variable', delay = 100, delayType = "debounce")),
uiOutput('hover_info_variable')),
tags$hr(),
downloadButton('variablePlot', "Plot")
),

#tabPanel("Testing", dataTableOutput('testtable')),
#tabPanel("Testing2", dataTableOutput('testtable2')),



tabPanel("Standards",
tabsetPanel(
tabPanel("Validation", dataTableOutput("standardsperformance")),
tabPanel("Used", rHandsontableOutput("whichrowstokeep"))))

))


))

)),

tabPanel("Apply Calibration",
div(class="outer",

fluidRow(
sidebarLayout(
sidebarPanel(width=3,

actionButton('processvalspectra', "Quantify"),


tags$hr(),

uiOutput('filevalgrab'),

selectInput("valfiletype", label="Filetype", c("DPT", "CSV", "Opus"), selected="DPT"),


tags$hr(),
tags$hr(),
tags$hr(),

fileInput('calfileinput2', 'Load Cal File', accept=".quant", multiple=FALSE),

downloadButton('downloadValData', "Results")

),


mainPanel(
tabsetPanel(
id = 'dataset2',
tabPanel('Validation', dataTableOutput('myvaltable2')),
tabPanel('Counts', dataTableOutput('myvaltable1'))

))
))
))

))
)
