## BF591 Final RShiny Project script

# Load required libraries
library(shiny)
library(shinythemes)
library(colourpicker) 
library(tidyverse)
library(glue)
library(DT)
library(ggplot2)
library(gridExtra)
library(scales)
library(gplots)
library(ggbeeswarm)
library(stats)


# Define UI for application that draws a histogram
ui <- fluidPage(
    
    theme = shinytheme("united"),
    
    # Application title
    titlePanel("BF591 Final Project: Bioinformatics Analysis Webapp"),
    p("paragraph describing what app does"),
    
    
    tabsetPanel(
      
      # ui for Tab1 - Sample information
      tabPanel("Sample", h3("Sample Information"), p("Load and examine sample information matrix of your dataset here"),
               sidebarLayout(
                 sidebarPanel(
                   fileInput(inputId = "file1","Load metadata:", 
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
      tabPanel("Counts", h3("Counts Analysis"), p("Explore and visualize counts matrices here"),
               sidebarLayout(
                 sidebarPanel(fileInput(inputId = "file2","Load normalized counts data:", 
                                        placeholder = "norm_counts.csv", accept = ".csv"),
                              h4("Select gene filters"),
                              sliderInput(inputId = "min_var", "Minimum percentile of variance:", 
                                          min = 0, max = 100, value = 50),
                              sliderInput(inputId = "min_samp", "Minimum number of non-zero samples:",
                                          min = 0, max = 69, value = 35),
                              submitButton("Submit")),
                
                 mainPanel(
                   tabsetPanel(
                     tabPanel("Summary", tableOutput("counts_summary")),
                     tabPanel("Diagnostic Plots", plotOutput("scatter_plots")),
                     tabPanel("Heatmap and Clustering", plotOutput("heatmap", width = '100%', height = '700px')),
                     tabPanel("PC Analysis", 
                              numericInput(inputId = "top", value = 2,
                                           "Select number of Principal Components (>1) you would like to visualize", width = '100%'),
                              plotOutput("top_pca", width = '100%', height = '700px')))),
      ),),
      
      # ui for Tab3 - Differential Expression Analysis result visualization
      tabPanel("Differential Expression", h3("Differential Expression Analysis"), 
               p("Load differential expression analysis results for visualization."),
               sidebarLayout(
                 sidebarPanel(fileInput(inputId = "file3","Load differential expression results:", placeholder = "deseq_res.csv", accept = ".csv"), 
                              radioButtons(inputId = "x_axis", "Choose the column for x-axis", choices=c("baseMean","log2FoldChange", "lfcSE", 
                                                                                                         "stat", "pvalue", "padj")),
                              radioButtons(inputId = "y_axis", "Choose the column for y-axis", choices=c("baseMean","log2FoldChange", "lfcSE", 
                                                                                                         "stat", "pvalue", "padj")),
                              colourInput(inputId = "base", "Select base point color", value = "#010100", showColour = c("both", "text", "background"),
                                          palette = c("square", "limited"), allowedCols = NULL, allowTransparent = FALSE, returnName = FALSE, closeOnClick = FALSE),
                              colourInput(inputId = "highlight", "Highlight point color", value = "#FC5E70", showColour = c("both", "text", "background"), 
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
      
      # ui for Tab4 - Visualization of Individual Gene Expression
      tabPanel("Individual Gene Expresion", h3("Individual Gene Expression Visualization"), 
               p("Upload the normalized counts and metadate here to visualize individual gene count by a desired sample information variable"),
               sidebarLayout(
                 sidebarPanel(
                   fileInput(inputId = 'file4_counts', label = 'Load normalized counts matrix CSV'),
                   fileInput(inputId = 'file4_meta', label = 'Load sample information matrix CSV'),
                   uiOutput("meta_selector"),
                   textInput("gene", label = "Enter gene to search for", placeholder = "ENSG00000000003.10"),
                   selectInput("plotType", label = "Choose what type of plot to make", choices = c("Bar", "Box", "Violin", "Beeswarm")),
                   submitButton(text='Go')),
                 mainPanel(plotOutput("plot_ige"))
                 )),
              ),
      )


# Define server logic required to draw a histogram
server <- function(input, output, session) {
  options(shiny.maxRequestSize=30*1024^2)
    
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
        geom_density(mapping=aes(x=PMI), fill="#FBBC54") +
        labs(title="Density plot of PMI", 
              x = "PMI", y = "density")+
        theme_linedraw()
      return(density_plot_1)
    }
    draw_density_2 <- function(file1){
      density_plot_2 <- ggplot(file1) +
        geom_density(mapping=aes(x=Age_of_death), fill="#F8DCB0") +
        labs(title="Density plot of Age of Death", 
             x = "Age of Death", y = "density")+
        theme_linedraw()
      return(density_plot_2)
    }
    draw_density_3 <- function(file1){
      density_plot_3 <- ggplot(file1) +
        geom_density(mapping=aes(x=RIN), fill="#FC5E70") +
        labs(title="Density plot of RIN", 
             x = "RIN", y = "density")+
        theme_linedraw()
      return(density_plot_3)
    }
    
    # Render output for first subtab
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
      
      grid.arrange(density_1, density_2, density_3, ncol = 3)
    })
    
    ## Tab 2
    # Function to load data 
    load_data_2 <- reactive({
      if (!is.null(input$file2)){
        file2 <- read_csv(input$file2$datapath)
        return(file2)}
      else{return(NULL)}
    })
    
    # Function to create summary table 
    make_counts_summary <- function(file3, min_var, min_samp) {
      # Filter genes based on variance and number of non-zero samples
      filtered_genes <- file3[rowSums(file3 > 0) >= min_samp & apply(file3, 1, var) >= min_var, ]

      # Get summary for the filtered data
      rows_filtered <- nrow(filtered_genes)
      rows_not_filtered <- nrow(file3) - rows_filtered
      gene_filter <- percent((rows_filtered / nrow(file3)))
      gene_unfilter <- percent((rows_not_filtered) / nrow(file3))

      # Create summary data frame
      summary_df <- data.frame(Description = c("Number of Genes", "Number of Samples", "Number of Genes (Filtered):", 
                                               "% of genes filtered", "Number of Genes (not Filtered:)", "%genes not filtered:"), 
                               Value = c(nrow(file3), ncol(file3), rows_filtered, gene_filter, rows_not_filtered, gene_unfilter))
      return(summary_df)
    }
    
    # Function to draw diagnostic scatter plots 
    draw_scatter1 <- function(file3, min_var){
      # filter data according to parameters
      med_data<-file3 %>% mutate(Median = apply(file3[-1], MARGIN = 1, FUN = median), 
                                  Variance = apply(file3[-1], MARGIN = 1, FUN = var))
      perc_val <- quantile(med_data$Variance, probs = min_var/100)   #calculate percentile
      med_data <- med_data %>% mutate(thresh = case_when(Variance >= perc_val ~ "TRUE", TRUE ~ "FALSE")) #sort genes by percentile
      
      # plot scatter plot
      cols <- c("FALSE" = "#FC5E70", "TRUE" = "#010100")
      scatter1 <- ggplot(med_data, aes(Median, Variance))+
        geom_point(aes(color=thresh), alpha=0.75)+
        scale_color_manual(values = cols)+
        labs(title = 'Plot of Median vs Variance.', subtitle = "Genes filtered out are in red. X and Y axes are log-scaled.")+
        scale_y_log10()+
        scale_x_log10()+
        theme_linedraw()+
        theme(legend.position = 'bottom')
      return(scatter1)
    }
    
    draw_scatter2 <- function(file3, min_samp){
      all_samples <- ncol(file3)-1  
      
      to_plot <- file3 %>% mutate(Median = apply(file3[-1], MARGIN = 1, FUN = median)) %>%
        mutate(Median = as.numeric(Median))   
      to_plot[to_plot == 0] <- NA
      to_plot$no_zeros <- rowSums(is.na(to_plot))  #make new col, with counts.
      to_plot <- to_plot %>% mutate(thresh = case_when(no_zeros <= min_samp ~ "TRUE", TRUE ~ "FALSE"))
      
      #plot scatter plot
      cols <- c("FALSE" = "#FC5E70", "TRUE" = "#010100")
      scatter2 <- ggplot(to_plot, aes(Median, no_zeros))+
        geom_point(aes(color=thresh), alpha=0.75)+
        scale_color_manual(values = cols)+
        scale_x_log10()+
        labs(title = 'Plot of Median vs Number of Non-Zero genes', subtitle = "Genes filtered out are in red.")+
        theme_linedraw()+
        ylab('Number of samples with zero count')+
        theme(legend.position = 'bottom')
      return(scatter2)
      }
    
    # Function to draw clustered heatmap of counts after filtering
    draw_heatmap <- function(file3, min_var, min_samp){
      file3[file3 == 0] <- NA
      #file3 <- na_if(file3, 0)
      file3$no_zeros <- rowSums(is.na(file3))  #make new col, with counts.
      file3 <- filter(file3, no_zeros <= min_samp)
      file3 <- log10(file3[,!colnames(file3) %in% c("gene", "no_zeros")]) #exclude the gene names column and log scale the values  
      
      # plot
      to_plot_heat <- file3 %>% 
        mutate(variance = apply(file3, MARGIN = 1, FUN = var)) #compute variance to filter the data
      min_val <- quantile(to_plot_heat$variance, probs = min_var/100, na.rm = TRUE)   #calculate percentile
      to_plot_heat <- filter(to_plot_heat, variance >= min_val) #filter the tibble
      hmap <- heatmap.2(as.matrix(to_plot_heat[-ncol(to_plot_heat)]), scale = "row")
      return(hmap)
    }
    
    # Function to draw scatterplot of PCA projections
    
    draw_pca <- function(file3, top, min_var){
      filtered_counts <- file3 
      filtered_counts <- filtered_counts %>% 
        mutate(variance = apply(filtered_counts[-1], MARGIN = 1, FUN = var), .after = gene) # calculate variance for filtering
      perc_val <- quantile(filtered_counts$variance, probs = min_var/100, na.rm = TRUE)   # calculate percentile
      filtered_counts <- filter(filtered_counts, variance >= perc_val) # filter the tibble
      pca_res <- prcomp(t(filtered_counts[-c(1,2)]), scale = FALSE) # transpose the data and perform PCA

      # Extract the scores and variances
      scores <- as.data.frame(pca_res$x)
      variances <- pca_res$sdev^2
      
      # Get the top n principal components
      top_pcs <- scores[, 1:top]
      
      # Filter out principal components with variance less than min_var
      top_pcs <- top_pcs[, variances[1:top] >= min_var]
      
      # Reshape the data for plotting with ggbeeswarm
      df <- top_pcs %>% pivot_longer(cols = 1:top, names_to = "PC", values_to = "Score")
      
      # plot top pcs 
      plotpca <- ggplot(df, aes(PC, Score)) +
        geom_quasirandom() +
        labs(x = paste0("PC1-PC", top), y = "Scores") +
        theme_linedraw()
      return(plotpca)
    }
    
    # Render output of summary table 
    output$counts_summary <- renderTable({
      input$Submit
      make_counts_summary(load_data_2())
    })
    
    # Render output of diagnostic plots 
    output$scatter_plots <- renderPlot({
      input$Submit
      scatter1 <- draw_scatter1(load_data_2(), min_var = input$min_var)
      scatter2 <- draw_scatter2(load_data_2(), min_samp = input$min_samp)
      grid.arrange(scatter1, scatter2, ncol=2)
    })
    
    # Render heatmap 
    output$heatmap <- renderPlot({
      input$Submit
      draw_heatmap(load_data_2(), input$min_var, input$min_samp)
    })
    
    # Render pca plot
    output$top_pca <- renderPlot({
      input$Submit
      draw_pca(load_data_2(), input$top, input$min_var)
    })

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
    
    ## Tab 4 
    # function to load normalized counts input file
    load_data_4_counts <- reactive({
      req(input$file4_counts)
      read_csv(input$file4_counts$datapath)
    })
    
    # function to load sample information input file
    load_data_4_meta <- reactive({
      req(input$file4_meta)
      read_csv(input$file4_meta$datapath)
    })
    
    # function to dynamically generate categorical variable options
    output$meta_selector <- renderUI({
      meta <- load_data_4_meta()
      if (!is.null(meta)) {
        categorical_vars <- names(meta)[sapply(meta, is.numeric)]
        selectInput("metachoice", "Select categorical variable:", choices = categorical_vars)
      }
    })
    
    # function to dynamically generate plot type options
    output$plotTypeSelector <- renderUI({
      selectInput("plotType", "Select plot type:", choices = c("Beeswarm", "Violin", "Bar"))
    })
    
    # function to make distribution plots
    plot_ige <- function(counts, meta, categorical_vars, select_gene, plot_type) {
      counts <- column_to_rownames(counts, var = "gene")
      gene_counts <- as.numeric(as.vector(counts[select_gene, ]))
      plot_tib <- tibble(Gene_Counts = gene_counts, meta_value = meta[[categorical_vars]])
      
      if (plot_type == "Beeswarm") {
        plot <- ggplot(plot_tib, aes(x = meta_value, y = Gene_Counts)) +
          geom_beeswarm() +
          theme_linedraw() +
          labs(title = str_c("Plot of gene counts vs ", categorical_vars))
      } else if (plot_type == "Violin") {
        plot <- ggplot(plot_tib, aes(x = meta_value, y = Gene_Counts)) +
          geom_violin() +
          theme_linedraw() +
          labs(title = str_c("Plot of gene counts vs ", categorical_vars))
      } else if (plot_type == "Bar") {
        plot <- ggplot(plot_tib, aes(x = meta_value)) +
          geom_bar() +
          theme_linedraw() +
          labs(title = str_c("Plot of gene counts vs ", categorical_vars))
      } else{
          plot <- ggplot(plot_tib, aes(x = meta_value, y = Gene_Counts)) +
            geom_boxplot() +
            theme_linedraw() +
            labs(title = str_c("Plot of gene counts vs ", categorical_vars))
      }
      
      return(plot)
    }
    
    # Render output plot
    output$plot_ige <- renderPlot({
      input$Go
      plot_ige(load_data_4_counts(), load_data_4_meta(), input$metachoice, input$gene, input$plotType)
    })
}

# Run the app ----
shinyApp(ui = ui, server = server)