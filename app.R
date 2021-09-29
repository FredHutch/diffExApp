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
library(pheatmap)

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
                    dataTableOutput("sample_matrix"),
                    hr(),
                    br(),
                    plotOutput("sample_pca")
                    )
                ) # end sidebar layout
            ), # end tabPanel 1
        
        # view results tab -----
        tabPanel("View Results",
            sidebarLayout(
                # filter results sidebar panel
                sidebarPanel(
                       h4("Filter Results"),
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
                                   step = .1),
                       checkboxInput("fdr",
                                     "Use adjusted P-value (False Discovery Rate)",
                                     TRUE),
                       checkboxInput("de_column",
                                     "Show differentially expressed column",
                                     FALSE),
                       checkboxInput("de_filter",
                                     "Filter table to only show differentially expressed genes",
                                     FALSE),
                       br()),
                # tab 2 main panel
                mainPanel(tabsetPanel(
                    # ...de results subtab -----
                    tabPanel("Differential Expression Results Table",
                             dataTableOutput("de_res_table"),
                             br(),
                             downloadButton("download_de_res",
                                            "Download Results Table"),
                             br(),
                             br(),
                             em("Results table will download exactly as it is displayed")),
                    # ...ma plot subtab -----
                    tabPanel("MA Plot",
                             fluidPage(
                                 # plot
                                 plotOutput("ma_plot"),
                                 # customization
                                 fluidRow(
                                     column(12,
                                            h4("Customize Plot"))),
                                 fluidRow(
                                     column(3,
                                     textInput("ma_x_label",
                                               "X axis label:",
                                               value = "Mean of normalized counts"),
                                     br(),
                                     downloadButton("download_ma",
                                                    "Download MA Plot")),
                                 column(4,
                                        offset = 1,
                                        textInput("ma_y_label",
                                                  "Y axis label:",
                                                  value = "Log fold change")),
                                 column(4,
                                        textInput("ma_legend_title",
                                                  "Legend title:",
                                                  value = "Differentially expressed")))
                                 
                             )),
                    # ...volcano plot subtab -----
                    tabPanel("Volcano Plot",
                             fluidPage(
                                 # plot
                                 plotOutput("volcano_plot"),
                                 # customization
                                 fluidRow(
                                     column(12,
                                            h4("Customize Plot"))),
                                 fluidRow(
                                     column(3,
                                            textInput("volcano_x_label",
                                                      "X axis label:",
                                                      value = "Log fold change"),
                                            br(),
                                            downloadButton("download_volcano",
                                                           "Download Volcano Plot")),
                                     column(4,
                                            offset = 1,
                                            textInput("volcano_y_label",
                                                      "Y axis label:",
                                                      value = "Significance (-log10)")),
                                     column(4,
                                            textInput("volcano_legend_title",
                                                      "Legend title:",
                                                      value = "Differentially expressed")))
                             )),
                    # ...heatmap subtab -----
                    tabPanel("Heatmap of DE Genes",
                             plotOutput("de_heatmap",
                                        height = "2500px"),
                             downloadButton("download_heatmap",
                                            "Download Heatmap"))
                    ) # end tabset panel within mainpanel
                ) # end mainPanel
        ) # end sidebarlayout
        ) # end tabsetPanel
    )) # end fluidPage

# SERVER -----
server <- function(input, output, session) {
    
    # initialize reactive values -----
    # these columns change based on de_package selected and if FDR = TRUE/FALSE
    col_names <- reactiveValues(logfc_col = NULL, 
                                pvalue_col = NULL,
                                mean_counts_col = NULL)
    
    # create sample matrix -----
    # condition 1 selector
    output$condition1_selector <- renderUI({
        selectInput("condition1_selected",
                    "Select replicates: ",
                    counts_cols,
                    multiple = TRUE,
                    selectize= TRUE)
    })
    # condition 2 selector
    output$condition2_selector <- renderUI({
        selectInput("condition2_selected",
                    "Select replicates: ",
                    counts_cols,
                    multiple = TRUE,
                    selectize= TRUE)
    })
    
    # create sample table based on user selections
    sample_matrix <- reactive({
        createSampleMatrix(condition1_name = input$condition1_name,
                           condition2_name = input$condition2_name,
                           condition1_selected = input$condition1_selected,
                           condition2_selected = input$condition2_selected)
    })
    
    # render sample table
    output$sample_matrix <- renderDataTable(
        sample_matrix() %>%
            rownames_to_column(var = "sample") %>%
            select(sample, everything())
    )
    
    # de analysis / update reactive values -----
    # on click (apply) 
    # update stored column names based on de_package selected
    # run differential expression
    de_out <- eventReactive(input$apply, {
        
        if (input$de_package == "DESeq2") {
            col_names$mean_counts_col <- "baseMean"
            col_names$logfc_col <- "log2FoldChange"
        } else if (input$de_package == "edgeR") {
            col_names$mean_counts_col <- "logCPM"
            col_names$logfc_col <- "logFC"
        }
        
        diffEx(data = data,
               condition1_name = input$condition1_name,
               condition2_name = input$condition2_name,
               condition1_selected = input$condition1_selected,
               condition2_selected = input$condition2_selected,
               gene_col = gene_col,
               de_package = input$de_package)
    })
    
    # update reactive values (fdr) -----
    # observe on checkbox click (fdr), colnames based on selected de_package 
    observeEvent(input$fdr, {
        if (input$de_package == "DESeq2") {
            col_names$pvalue_col <- ifelse(input$fdr, "padj", "pvalue")
        } else {
            col_names$pvalue_col <- ifelse(input$fdr, "FDR", "PValue")
        }
    })

    # get de results -----
    de_res <- reactive({
        getResults(de_out = de_out(),
                   de_package = input$de_package)
    })
    
    # get norm counts -----
    counts <- reactive({
        getCounts(de_out = de_out(),
                  de_package = input$de_package,
                  normalized = TRUE)
    })
    
    # format results based on user inputs -----
    # this is specifically for the downloadable results table
    formatted_res_render <- reactive({
        formatResults(de_res = req(de_res()),
                      logfc_threshold = input$logfc_threshold,
                      pvalue_threshold = input$pvalue_threshold,
                      logfc_col = col_names$logfc_col,
                      pvalue_col = col_names$pvalue_col,
                      de_column = input$de_column,
                      de_filter = input$de_filter)
    })
    
    # render results table
    # will re-render baesd on user selections
    output$de_res_table <- renderDataTable(
        DT::datatable(formatted_res_render(), options = list(pageLength = 50))
    )
    
    # format results for plotting -----
    # plots cannot use data that is filtered by user selection and require the de_column!
    # FIXME: maybe there is a more streamlined way to do this?
    formatted_res_plots <- reactive({
        formatResults(de_res = req(de_res()),
                      logfc_threshold = input$logfc_threshold,
                      pvalue_threshold = input$pvalue_threshold,
                      logfc_col = col_names$logfc_col,
                      pvalue_col = col_names$pvalue_col,
                      de_column = TRUE,
                      de_filter = FALSE)
    })

    # pca plot -----    
    # use counts from de_out
    # plot on apply click
    pca <-  eventReactive(input$apply, {
       plotPca(counts = counts(),
               sample_matrix = sample_matrix())
    })
    
    # render pca
    output$sample_pca <- renderPlot({
        pca()
    })
    
    # heatmap -----
    heatmap <- reactive({
        plotHeatmap(counts = counts(),
                    sample_matrix = sample_matrix(),
                    de_vec = formatted_res_plots()$isDE,
                    silent = FALSE)
    })
    
    # output$de_heatmap <- renderPlot({
    #     heatmap()
    # }, height = function() {
    #     session$clientData$output_de_heatmap_width
    # })
    
    output$de_heatmap <- renderPlot({
        heatmap()
    })

    # ma plot -----
    ma <- reactive({
        resultsToMa(de_res = formatted_res_plots(),
                    logfc_col = col_names$logfc_col,
                    mean_counts_col = col_names$mean_counts_col,
                    x_label = input$ma_x_label,
                    y_label = input$ma_y_label,
                    legend_title = input$ma_legend_title,
                    de_vec = formatted_res_plots()$isDE)
    })
    
    output$ma_plot <- renderPlot(
        ma()
    )
    
    # volcano plot -----
    volcano <- reactive({
        plotVolcano(data = formatted_res_plots(),
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
                    de_vec = formatted_res_plots()$isDE)
    })
    
    output$volcano_plot <- renderPlot(
        volcano()
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
            paste0("diffEx-", input$de_package, "-ma-plot-", Sys.Date(), ".pdf")
        },
        
        content = function(file) {
            ggsave(file, ma(), device = "pdf", width = 10, height = 5, units = "in")
        })
    
    output$download_volcano <- downloadHandler(
        filename = function() {
            paste0("diffEx-", input$de_package, "-volcano-", Sys.Date(), ".pdf")
        },
        
        content = function(file) {
            ggsave(file, volcano(), device = "pdf", width = 10, height = 5, units = "in")
        })
    
    output$download_heatmap <- downloadHandler(
        filename = function() {
            paste0("diffEx-", input$de_package, "-de-heatmap-", Sys.Date(), ".pdf")
        },
        
        content = function(file) {
            ggsave(file, heatmap(), device = "pdf", width = 10, height = 5, units = "in")
        })
}

# Run the application 
shinyApp(ui, server)