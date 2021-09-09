# save variables
# data <- read.table("~/Documents/work/dataCore/shiny/diffEx/data/counts.csv", sep = ",", header = TRUE)
# condition1_name <- "wt"
# condition2_name <- "trt"
# condition1_selected <- c("SRR1039508", "SRR1039509", "SRR1039512")
# condition2_selected <- c("SRR1039521", "SRR1039520", "SRR1039517")
# gene_col <- "gene"
# de_package <- "DESeq2"
# fdr <- TRUE
# pvalue_threshold <- .01
# logfc_threshold <- 2

# SET UP -----
# load libraries
library(shiny)
library(DESeq2)
library(edgeR)
library(tidyverse)
library(DT)
library(scales)

# source scripts
source("./diffEx.R")
source("./plot.R")

# read in data
data <- read.table("~/Documents/work/dataCore/shiny/diffEx/data/counts.csv", sep = ",", header = TRUE)

# check cols
counts_candidate_f <- function(x) {
    if (is.numeric(data[[x]])) {
        return(TRUE)}
    return(FALSE)
}

gene_candidate_f <- function(x) {
    if (is.character(data[[x]])) {
        return(TRUE)}
    return(FALSE)
}

counts_cols <- names(data)[sapply(names(data), counts_candidate_f)]
gene_col <- names(data)[sapply(names(data), gene_candidate_f)]

# UI -----
ui <- fluidPage(
    # Use tabset panel
    tabsetPanel(
        
        # run analysis tab -----
        tabPanel("Run Analysis",
            sidebarLayout(
                # tab 1 sidebar panel
                sidebarPanel(
                             br(),
                             h4("Run Analysis"),
                             textInput("condition1_name",
                                       "Condition label",
                                       value = "control"),
                             uiOutput("condition1_selector"),
                             textInput("condition2_name",
                                       "Condition label",
                                       value = "treatment"),
                             uiOutput("condition2_selector"),
                             selectInput("de_package",
                                         "Select preferred differential analysis package:",
                                         choices = c("DESeq2", "edgeR")),
                             actionButton("apply", "Apply")),
                # tab 1 main panel
                mainPanel(
                    h4("Sample Table"),
                    dataTableOutput("sample_table"))
                ) # end sidebar layout
            ), # end tabPanel 1
        
        # view results tab -----
        tabPanel("View Results",
            sidebarLayout(
                # tab 2 sidebar panel
                sidebarPanel(
                       h4("Filter Dataset"),
                       sliderInput("pvalue_threshold",
                                   "Set significance threshold",
                                   min = 0,
                                   max = 1,
                                   value = .05,
                                   step = .01),
                       checkboxInput("fdr",
                                     "Use False Discovery Rate",
                                     TRUE),
                       sliderInput("logfc_threshold",
                                   "Set log fold change threshold",
                                   min = 0,
                                   max = 10,
                                   value = 2,
                                   step = .5),
                       checkboxInput("fdr",
                                     "Use adjusted P-value (False Discovery Rate)",
                                     TRUE),
                       checkboxInput("de_column",
                                     "Show differentially expressed column",
                                     FALSE),
                       checkboxInput("de_filter",
                                     "Filter table to only show differentially expressed genes",
                                     FALSE),
                       br(),
                       downloadButton("download_de_res",
                                      "Download Results Table"),
                       br(),
                       br(),
                       em("Results table will download exactly as it is displayed")),
                # tab 2 main panel
                mainPanel(
                    h4("Differential Expression Results Table"),
                    dataTableOutput("de_res_table"))
                ) # end tab 2 sidebar layout
        ), # end tabpanel 2
        tabPanel("Plots",
                 fluidRow(
                     column(width = 12,
                            align = "center",
                            h4("MA Plot"),
                            plotOutput("ma_plot"),
                            br(),
                            h4("Volcano Plot"),
                            plotOutput("volcano_plot"))
                 )
                 )
        ) # end tabsetPanel
    ) # end fluidPage

# SERVER -----
server <- function(input, output) {
    
    # select condition colnames -----
    # condition 1
    output$condition1_selector <- renderUI({
        selectInput("condition1_selected",
                    "Select replicates: ",
                    counts_cols,
                    multiple = TRUE,
                    selectize= TRUE)
    })
    # condition 2
    output$condition2_selector <- renderUI({
        selectInput("condition2_selected",
                    "Select replicates: ",
                    counts_cols,
                    multiple = TRUE,
                    selectize= TRUE)
    })
    
    # Differential expression -----
    # create sample table
    reactive_sample_table <- reactive({
        createSampleMatrix(condition1_name = input$condition1_name,
                           condition2_name = input$condition2_name,
                           condition1_selected = input$condition1_selected,
                           condition2_selected = input$condition2_selected)
    })
    
    # render sample table 
    output$sample_table <- renderDataTable(
        reactive_sample_table() %>%
            rownames_to_column(var = "sample") %>%
            select(sample, everything())
    )
    
    # run differential expression analysis
    reactive_de_res <- reactive({
        diffEx(data = data,
              condition1_name = input$condition1_name,
              condition2_name = input$condition2_name,
              condition1_selected = input$condition1_selected,
              condition2_selected = input$condition2_selected,
              gene_col = gene_col,
              de_package = input$de_package)
    })
    
    # filter de results
    de_res_formatted <- reactive({
        formatResults(de_res = reactive_de_res(),
                      pvalue_threshold = input$pvalue_threshold,
                      logfc_threshold = input$logfc_threshold,
                      fdr = input$fdr,
                      de_column = input$de_column,
                      de_filter = input$de_filter,
                      de_package = input$de_package)
    })
    
    # render data table of results on apply button click
    # FIXME: maybe I should have the button run the analysis instead...?
    observeEvent(input$apply, {
        output$de_res_table <- renderDataTable(
            DT::datatable(de_res_formatted(), options = list(pageLength = 25))
        )
    })
    
    # ma plot -----
    reactive_ma <- reactive({
        resultsToMa(results = reactive_de_res(),
                    de_package = input$de_package,
                    pvalue_threshold = input$pvalue_threshold,
                    logfc_threshold = input$logfc_threshold,
                    fdr = input$fdr)
    })
    
    output$ma_plot <- renderPlot(
        reactive_ma()
    )
    
    # volcano plot -----
    reactive_volcano <- reactive({
        resultsToVolcano(result = reactive_de_res(),
                         de_package = input$de_package,
                         pvalue_threshold = input$pvalue_threshold,
                         logfc_threshold = input$logfc_threshold,
                         fdr = input$fdr)
    })
    
    output$volcano_plot <- renderPlot(
        reactive_volcano()
    )
    
    # download handler -----
    
    output$download_de_res <- downloadHandler(
        filename = function() {
            paste0("diffEx-", input$de_package, "-", Sys.Date(), ".csv")
        },
        
        content = function(file) {
            write.csv(de_res_formatted(), file)
        })
}

# Run the application 
shinyApp(ui = ui, server = server)