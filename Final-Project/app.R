#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(shinythemes)
library(colourpicker) # you might need to install this.
library(tidyverse)
library(glue)
library(DT)
library(ggplot2)
#source("~/BF591-Final-Project-rpragya17/Final-Project/functions.R")

# Define UI for application that draws a histogram
ui <- fluidPage(
    
    theme = shinytheme("united"),
    
    # Application title
    titlePanel("BF591 Final Project: Bioinformatics Analysis Webapp"),
    p("paragraph describing what app does"),
    
    
    tabsetPanel(
      
      # ui for Tab1 - Sample information
      tabPanel("Sample", h3("Sample Information"), p("Load dataset here to explore pre-processed data"),
               sidebarLayout(
                 sidebarPanel(
                   fileInput(inputId = "file1","Load gene expression dataset:", 
                             placeholder = "metadata.csv", accept = ".csv"),
                 ),
                 mainPanel(
                   tabsetPanel(
                     tabPanel("Summary", tableOutput("summary_table")),
                     tabPanel("Data", dataTableOutput("data1")),
                     tabPanel("Plots", plotOutput("density_plots")))),
                 ),
      ),
      # ui for Tab2 - Counts data analysis
      tabPanel("Counts",
               sidebarLayout(
                 sidebarPanel("Input data here"),
                 mainPanel(
                   tabsetPanel(
                     tabPanel("Sub-tab 1A", "This is sub-tab 1A content."),
                     tabPanel("Sub-tab 1B", "This is sub-tab 1B content.")))),
      ),
      
      # ui for Tab3 - Differential Expression Analysis result visualization
      tabPanel("Differential Expression", h3("Differential Expression Analysis"), 
               p("Load differential expression analysis results for visualization."),
               sidebarLayout(
                 sidebarPanel(fileInput(inputId = "file3","Load differential expression results:", placeholder = "deseq_res.csv", accept = ".csv"), 
                              radioButtons(inputId = "x_axis", "Choose the column for x-axis", choices=c("baseMean","log2FoldChange", "lfcSE", 
                                                                                                         "stat", "pvalue", "padj")),
                              radioButtons(inputId = "y_axis", "Choose the column for y-axis", choices=c("baseMean","log2FoldChange", "lfcSE", 
                                                                                                         "stat", "pvalue", "padj")),
                              colourInput(inputId = "base", "Select base point color", value = "midnightblue", showColour = c("both", "text", "background"),
                                          palette = c("square", "limited"), allowedCols = NULL, allowTransparent = FALSE, returnName = FALSE, closeOnClick = FALSE),
                              colourInput(inputId = "highlight", "Highlight point color", value = "springgreen", showColour = c("both", "text", "background"), 
                                          palette = c("square", "limited"), allowedCols = NULL, allowTransparent = FALSE, returnName = FALSE, closeOnClick = FALSE),
                              sliderInput(inputId = "slider", "Select the magnitude of the p adjusted coloring:", min = -100, max = 0,
                                          value = -50),
                              submitButton("plot")
                 ),
                 
                 mainPanel(
                   tabsetPanel(
                     tabPanel("Table", dataTableOutput("deg_table")),
                     tabPanel("Volcano plot", plotOutput("volcano", width = '100%', height = '700px')),
                     tabPanel("Volcano plot results", dataTableOutput("volc_table")))),
                  ),
              ),
      
      # ui for Tab4 
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
server <- function(input, output, session) {
    
    ## Tab 1 
    # Function to load data
    load_data_1 <- reactive({
      if (!is.null(input$file1)){
        file1 <- read_csv(input$file1$datapath)
        return(file1)}
      else{return(NULL)}
    })
    
    # Function to create summary table 
    make_summary <- function(file1){
      Columns <- c(colnames(file1))
      Type <- sapply(file1, class)
      Mean <- sapply(file1, mean, na.rm = TRUE) 
      SD <- sapply(file1, sd, na.rm = TRUE)
      summary_table <- data.frame(Columns, Type, Mean, SD, stringsAsFactors = FALSE)
      return(summary_table)
    }
    
    # Functions to create violin plot for continuous variables in metadata
    draw_density_1 <- function(file1){
      density_plot_1 <- ggplot(file1) +
        geom_density(mapping=aes(x=PMI), fill="salmon") +
        labs(title="Density plot of PMI", 
              x = "PMI", y = "density")+
        theme_classic()
      return(density_plot_1)
    }
    draw_density_2 <- function(file1){
      density_plot_2 <- ggplot(file1) +
        geom_density(mapping=aes(x=Age_of_death), fill="salmon") +
        labs(title="Density plot of Age of Death", 
             x = "Age of Death", y = "density")+
        theme_classic()
      return(density_plot_2)
    }
    draw_density_3 <- function(file1){
      density_plot_3 <- ggplot(file1) +
        geom_density(mapping=aes(x=RIN), fill="salmon") +
        labs(title="Density plot of RIN", 
             x = "RIN", y = "density")+
        theme_classic()
      return(density_plot_3)
    }
    
    # Render output for Tab 1
    output$summary_table <- renderTable({
      make_summary(load_data_1())
    })
    
    # Render output of second subtab as a table 
    output$data1 <- DT::renderDataTable({
      DT::datatable(load_data_1())
    })
    
    # Render density plots for third subtab 
    output$density_plots <- renderPlot({
      density_1 <- draw_density_1(load_data_1())
      density_2 <- draw_density_2(load_data_1())
      density_3 <- draw_density_3(load_data_1())
      
      grid.arrange(plot1, plot2, plot3, ncol = 3)
      #ggarrange(plot1, plot2, plot3, ncol = 3, nrow = 1)
    })
    # output$density1 <- renderPlot({
    #   draw_density_1(load_data_1())
    # })
    # 
    # output$density2 <- renderPlot({
    #   draw_density_2(load_data_1())
    # })
    # 
    # output$density3 <- renderPlot({
    #   draw_density_3(load_data_1())
    # })
    
  
    ## Tab 3
    # Function to load data 
    load_data_3 <- reactive({
      if (!is.null(input$file3)){
        file3 <- read_csv(input$file3$datapath)
        return(file3)}
      else{return(NULL)}
    })
    
    # Function to plot volcano plot using input filters from sliders 
    volcano_plot <- function(file3, x_name, y_name, slider, color1, color2) {
        dataf <- file3 %>% drop_na()
        v_plot <- ggplot(dataf, aes(x= !!sym(x_name), y= -log10(!!sym(y_name)), color= color))+ 
          geom_point(aes(color = if_else(dataf[[y_name]] < 10^slider, TRUE, FALSE))) +   
          scale_color_manual(name = glue("{x_name} < 1*10^{slider}"),
                             values = c("FALSE" = color1,
                                        "TRUE" = color2),
                             labels = c("FALSE", "TRUE")) + 
          labs(title = "Volcano Plot", x = x_name, y = y_name) +
          theme_linedraw() +
          theme(plot.title = element_text(hjust = 0.5), legend.position = "bottom")
        return(v_plot)
    }
    
    # Function to generate table with data for degs filtered in volcano plot 
    draw_table <- function(file3, slider) {
      new_df <- subset(file3, padj < 10^slider)
      # Format p-value columns to display more digits
      new_df$pvalue <- formatC(new_df$pvalue, digits = 5, format = "e")
      new_df$padj <- formatC(new_df$padj, digits = 5, format = "e")
      colnames(new_df)[1] <- "gene"
      
      # Return filtered and formatted data frame
      return(new_df)
    }
    
    ## Tab 3 
    # Render output of first tab as a table 
    output$deg_table <- DT::renderDataTable({
      DT::datatable(load_data_3())
    })
    
    # Render output of second tab as plot
    output$volcano <- renderPlot({
      input$plot
      data3 <- load_data_3()
      volcano_plot(file3 = data3, x_name = input$x_axis, y_name = input$y_axis, 
                   slider = input$slider, color1 = input$base, color2 = input$highlight)
    })
    
    # Render output of third tab as a table 
    output$volc_table <- DT::renderDataTable({
      input$plot
      draw_table(load_data_3(), input$slider)
    })
}

# Run the app ----
shinyApp(ui = ui, server = server)