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
library(reactlog)

reactlog_enable()
# source scripts
source("./diffEx.R")
source("./plot.R")
source("https://raw.githubusercontent.com/FredHutch/interactiveVolcano/master/volcano.R")

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
                    dataTableOutput("sample_table"),
                    hr(),
                    br(),
                    plotOutput("sample_pca")
                    )
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
        tabPanel("MA Plot",
                 sidebarLayout(
                     sidebarPanel(
                         h4("Customize Plot"),
                         textInput("ma_x_label",
                                   "X axis label:",
                                   value = "Mean of normalized counts"),
                         textInput("ma_y_label",
                                   "Y axis label:",
                                   value = "Log fold change"),
                         textInput("ma_legend_title",
                                   "Legend title:",
                                   value = "Differentially expressed")),
                     mainPanel(
                         plotOutput("ma_plot"),
                         downloadButton("download_ma",
                                        "Download MA Plot")))
                 ), # end tab 3
        tabPanel("Volcano Plot",
                 sidebarLayout(
                     sidebarPanel(
                         h4("Customize Plot"),
                         textInput("volcano_x_label",
                                   "X axis label:",
                                   value = "Log fold change"),
                         textInput("volcano_y_label",
                                   "Y axis label:",
                                   value = "Significance (-log10)"),
                         textInput("volcano_legend_title",
                                   "Legend title:",
                                   value = "Differentially expressed")),
                     mainPanel(
                         plotOutput("volcano_plot"),
                         downloadButton("download_volcano",
                                        "Download Volcano Plot")))
        ) # end tab 4
        ) # end tabsetPanel
    ) # end fluidPage

# SERVER -----
server <- function(input, output) {
    
    # reactive values of column names
    # these are updated every time the action button 
    # is clicked based on selected de_package
    col_names <- reactiveValues(logfc_col = NULL, 
                                pvalue_col = NULL,
                                fdr_col = NULL,
                                mean_counts_col = NULL)
    
    # select condition colnames to create sample matrix -----
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
    
    # Show selected sample table -----
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
    
    # Differential expression analysis -----
    # on click: 
    # update stored column names based on de_package selected
    # run differential expression
    de_out <- eventReactive(input$apply, {
        
        if (input$de_package == "DESeq2") {
            col_names$mean_counts_col <- "baseMean"
            col_names$logfc_col <- "log2FoldChange"
            col_names$pvalue_col <- ifelse(input$fdr == TRUE, "padj", "pvalue")
        } else if (input$de_package == "edgeR") {
            col_names$mean_counts_col <- "logCPM"
            col_names$logfc_col <- "logFC"
            col_names$pvalue_col <- ifelse(input$fdr == TRUE, "FDR", "PValue")
        }
        
        diffEx(data = data,
               condition1_name = input$condition1_name,
               condition2_name = input$condition2_name,
               condition1_selected = input$condition1_selected,
               condition2_selected = input$condition2_selected,
               gene_col = gene_col,
               de_package = input$de_package)
    })

    # get de results from de object
    # this reactive element will go into plots 
    de_res_table <- reactive({
        getResults(de_out = de_out(),
                   logfc_col = col_names$logfc_col,
                   pvalue_col = col_names$pvalue_col,
                   fdr_col = col_names$fdr_col,
                   logfc_threshold = input$logfc_threshold,
                   pvalue_threshold = input$pvalue_threshold,
                   fdr = input$fdr,
                   de_package = input$de_package)
    })
    
    # format results based on user inputs
    formatted_res_table <- reactive({
        formatResults(de_res = de_res_table(),
                      de_column = input$de_column,
                      de_filter = input$de_filter)
    })
    
    
    # render de_res_table()
    output$de_res_table <- renderDataTable(
        DT::datatable(formatted_res_table(), options = list(pageLength = 25))
    )
    
    # plots -----
    
    # plot pca using counts from de_out
    # plot pca
    output$sample_pca <- renderPlot({
        countsToPca(de_out = de_out(),
                    sample_matrix = reactive_sample_table(),
                    de_package = input$de_package)
    })

    # ma plot -----
    reactive_ma <- reactive({
        resultsToMa(de_res = de_res_table(),
                    logfc_col = col_names$logfc_col,
                    mean_counts_col = col_names$mean_counts_col,
                    x_label = input$ma_x_label,
                    y_label = input$ma_y_label,
                    legend_title = input$ma_legend_title,
                    de_vec = de_res_table()$isDE)
    })
    
    output$ma_plot <- renderPlot(
        reactive_ma()
    )
    
    # volcano plot -----
    reactive_volcano <- reactive({
        plotVolcano(data = de_res_table(),
                    logfc_col = col_names$logfc_col,
                    pvalue_col = col_names$pvalue_col,
                    gene_col = NULL,
                    pvalue_thresh = input$pvalue_threshold,
                    logfc_thresh = input$logfc_threshold,
                    color_by_de = TRUE,
                    show_logfc_thresh = TRUE,
                    show_pvalue_thresh = TRUE,
                    highlight_genes = NULL,
                    x_label = input$volcano_x_label,
                    y_label = input$volcano_y_label,
                    legend_title = input$volcano_legend_title,
                    xlim = NULL,
                    ylim = NULL,
                    de_vec = de_res_table()$isDE)
    })
    
    output$volcano_plot <- renderPlot(
        reactive_volcano()
    )
    
    # downloads -----
    
    output$download_de_res <- downloadHandler(
        filename = function() {
            paste0("diffEx-", input$de_package, "-results-", Sys.Date(), ".csv")
        },
        
        content = function(file) {
            write.csv(de_res_formatted_table(), file)
        })
    
    output$download_ma <- downloadHandler(
        filename = function() {
            paste0("diffEx-", input$de_package, "-ma-plot-", Sys.Date(), ".csv")
        },
        
        content = function(file) {
            ggsave(file, reactive_ma(), device = "pdf", width = 10, height = 5, units = "in")
        })
    
    output$download_volcano <- downloadHandler(
        filename = function() {
            paste0("diffEx-", input$de_package, "-volcano-", Sys.Date(), ".csv")
        },
        
        content = function(file) {
            ggsave(file, reactive_volcano(), device = "pdf", width = 10, height = 5, units = "in")
        })
}

# Run the application 
shinyApp(ui, server)