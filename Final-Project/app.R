#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)

# Define UI for application that draws a histogram
ui <- fluidPage(

    # Application title
    titlePanel("BF591 Final Project: Bioinformatics Analysis Webapp"),
    p("paragraph describing what app does"),
    
    
    tabsetPanel(
      tabPanel("Sample", h3("Sample Information"), p("Load dataset here and explore pre-processed data"),
               sidebarLayout(
                 sidebarPanel(
                   fileInput(inputId = "file","Load gene expression dataset:", 
                             placeholder = "deseq_res.csv", accept = ".csv"),
                 ),
                 mainPanel(
                   tabsetPanel(
                     tabPanel("Summary", "This is sub-tab 1A content."),
                     tabPanel("Data", "This is sub-tab 1B content."),
                     tabPanel("Plots", "This is sub-tab 1B content.")))),
               ),
      tabPanel("Counts",
               sidebarLayout(
                 sidebarPanel("Input data here"),
                 mainPanel(
                   tabsetPanel(
                     tabPanel("Sub-tab 1A", "This is sub-tab 1A content."),
                     tabPanel("Sub-tab 1B", "This is sub-tab 1B content.")))),
      ),
      tabPanel("DE",
               sidebarLayout(
                 sidebarPanel("Input data here"),
                 mainPanel(
                   tabsetPanel(
                     tabPanel("Sub-tab 1A", "This is sub-tab 1A content."),
                     tabPanel("Sub-tab 1B", "This is sub-tab 1B content.")))),
      ),
      tabPanel("Tab4",
               sidebarLayout(
                 sidebarPanel("Input data here"),
                 mainPanel(
                   tabsetPanel(
                     tabPanel("Sub-tab 1A", "This is sub-tab 1A content."),
                     tabPanel("Sub-tab 1B", "This is sub-tab 1B content.")))),
      ),
    )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
    
    # Function to load data 
    load_data <- reactive({
      req(input$file)
      dataf <- read.csv(input$file$datapath)
      return(dataf)
    })
    
    # Function to create summary table 
    summ_table <- 

    # output$distPlot <- renderPlot({
    #     # generate bins based on input$bins from ui.R
    #     x    <- faithful[, 2]
    #     bins <- seq(min(x), max(x), length.out = input$bins + 1)
    # 
    #     # draw the histogram with the specified number of bins
    #     hist(x, breaks = bins, col = 'darkgray', border = 'white',
    #          xlab = 'Waiting time to next eruption (in mins)',
    #          main = 'Histogram of waiting times')
    # })
}

# Run the app ----
shinyApp(ui = ui, server = server)