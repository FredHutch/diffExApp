# save variables
# data <- read.table("~/Documents/work/dataCore/shiny/diffEx/data/counts.csv", sep = ",", header = TRUE)
# condition1_name <- "wt"
# condition2_name <- "trt"
# condition1_selected <- c("SRR1039508", "SRR1039509")
# condition2_selected <- c( "SRR1039512", "SRR1039513")
# gene_col <- "gene"
# de_package <- "edgeR"
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
                sidebarPanel(width = 3,
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
                    plotOutput("sample_pca")))
            ), # end tabPanel 1
        
        # view results tab -----
        tabPanel("View Results",
            sidebarLayout(
                # filter results sidebar panel
                sidebarPanel(width = 3,
                             h4("Set Differential Expression Thresholds"),
                             sliderInput("pvalue_threshold",
                                         "Set significance threshold",
                                         min = 0,
                                         max = 1,
                                         value = .05,
                                         step = .01),
                             checkboxInput("fdr",
                                           "Use adjusted P value (Benjamini-Hochberg)",
                                           TRUE),
                             sliderInput("logfc_threshold",
                                         "Set log fold change threshold",
                                         min = 0,
                                         max = 10,
                                         value = 2,
                                         step = .1),
                             radioButtons("de_filter",
                                          "DE filtering options:",
                                          choices = c("All differentially expressed" = "both",
                                                      "Upregulated only" = "up",
                                                      "Downregulated only" = "down"),
                                          selected = "both"),
                             br()),
                # tab 2 main panel
                mainPanel(
                    tabsetPanel(
                    # ...de results subtab -----
                    tabPanel("Differential Expression Results Table",
                             fluidPage(
                                 # results table
                                 dataTableOutput("de_res_table"),
                                 # customize
                                 fluidRow(h4("Subset Table"),
                                          em("Add/Remove columns, filter, etc")),
                                 fluidRow(
                                     column(2,
                                            br(),
                                            checkboxInput("de_column",
                                                          "Show isDE column",
                                                          value = FALSE)),
                                     column(2,
                                            br(),
                                            checkboxInput("de_subset",
                                                          "Only show DE genes",
                                                          value = FALSE))),
                                 hr(),
                                 # customization
                                 fluidRow(h4("Download")),
                                 fluidRow(
                                     column(1,
                                            style = "margin-top: 25px;",
                                            downloadButton("download_de_res",
                                                           "Download Results Table")),
                                     column(1,
                                            offset = 1,
                                            selectInput("results_file_type",
                                                        "File type",
                                                         choices = c(".tsv", ".csv", ".txt"),
                                                         width = "80px")))
                                 )),
                    # ...ma plot subtab -----
                    tabPanel("MA Plot",
                             fluidPage(
                                 # plot
                                 plotOutput("ma_plot"),
                                 # customization
                                 fluidRow(h4("Customize Plot"),
                                          em("Add custom axes lables and/or legend title")),
                                 fluidRow(
                                     column(3,
                                            br(),
                                     textInput("ma_x_label",
                                               "X axis label:",
                                               value = "Mean of normalized counts")),
                                 column(4,
                                        offset = 1,
                                        br(),
                                        textInput("ma_y_label",
                                                  "Y axis label:",
                                                  value = "Log fold change")),
                                 column(4,
                                        br(),
                                        textInput("ma_legend_title",
                                                  "Legend title:",
                                                  value = "Differentially expressed"))),
                                 hr(),
                                 fluidRow(
                                     h4("Download"),
                                     em("Download plot as a .pdf. Select desired height and width of plot in inches.")),
                                 fluidRow(
                                     column(1,
                                            style = "margin-top: 25px;",
                                            br(),
                                            downloadButton("download_ma",
                                                           "Download MA Plot")),
                                     column(1,
                                            offset = 1,
                                            br(),
                                            numericInput("ma_h",
                                                         "Height (in)",
                                                         value = 5,
                                                         width = "80px")),
                                     column(1,
                                            br(),
                                            numericInput("ma_w",
                                                         "Width (in)",
                                                         value = 10,
                                                         width = "80px"))
                                 )
                             )),
                    # ...volcano plot subtab -----
                    tabPanel("Volcano Plot",
                             fluidPage(
                                 # plot
                                 plotOutput("volcano_plot"),
                                 # customization
                                 fluidRow(h4("Customize Plot"),
                                          em("Add custom axes lables and/or legend title")),
                                 fluidRow(
                                     column(3,
                                            br(),
                                            textInput("volcano_x_label",
                                                      "X axis label:",
                                                      value = "Log fold change")),
                                     column(4,
                                            offset = 1,
                                            br(),
                                            textInput("volcano_y_label",
                                                      "Y axis label:",
                                                      value = "Significance (-log10)")),
                                     column(4,
                                            br(),
                                            textInput("volcano_legend_title",
                                                      "Legend title:",
                                                      value = "Differentially expressed"))),
                                 hr(),
                                 fluidRow(
                                     h4("Download"),
                                     em("Download plot as a .pdf. Select desired height and width of plot in inches.")),
                                 fluidRow(
                                     column(1,
                                            style = "margin-top: 25px;",
                                            br(),
                                            downloadButton("download_volcano",
                                                           "Download Volcano Plot")),
                                     column(1,
                                            offset = 1,
                                            br(),
                                            numericInput("volcano_h",
                                                         "Height (in)",
                                                         value = 5,
                                                         width = "80px")),
                                     column(1,
                                            br(),
                                            numericInput("volcano_w",
                                                         "Width (in)",
                                                         value = 10,
                                                         width = "80px"))
                                 )
                             )),
                    # ...heatmap subtab -----
                    tabPanel("Heatmap of DE Genes",
                             # plot
                             plotOutput("de_heatmap", height = "auto"),
                             # customization
                             fluidRow(h4("Customize Heatmap")),
                             fluidRow(
                                 column(2,
                                        checkboxInput("heatmap_gene_names",
                                                      "Show gene labels",
                                                      value = TRUE)),
                                 column(2,
                                        checkboxInput("heatmap_sample_names",
                                                      "Show sample labels",
                                                      value = TRUE))),
                             hr(),
                             fluidRow(
                                 h4("Download"),
                                 em("Download plot as a .pdf. Select desired height and width of plot in inches.")),
                             fluidRow(
                                 column(1,
                                        style = "margin-top: 25px;",
                                        br(),
                                        downloadButton("download_heatmap",
                                                       "Download Heatmap")),
                                 column(1,
                                        offset = 1,
                                        br(),
                                        numericInput("heatmap_h",
                                                     "Height (in)",
                                                     value = 10,
                                                     width = "80px")),
                                 column(1,
                                        br(),
                                        numericInput("heatmap_w",
                                                     "Width (in)",
                                                     value = 5,
                                                     width = "80px"))
                             ))
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
    # update all cols on apply button click
    observeEvent(input$apply, {
        if (input$de_package == "DESeq2") {
            col_names$mean_counts_col <- "baseMean"
            col_names$logfc_col <- "log2FoldChange"
            col_names$pvalue_col <- ifelse(input$fdr, "padj", "pvalue")
        } else if (input$de_package == "edgeR") {
            col_names$mean_counts_col <- "logCPM"
            col_names$logfc_col <- "logFC"
            col_names$pvalue_col <- ifelse(input$fdr, "FDR", "PValue")
        }
    })
    
    # update pvalue col on fdr checkbox click
    observeEvent(input$fdr, {
        if (input$de_package == "DESeq2") {
            col_names$pvalue_col <- ifelse(input$fdr, "padj", "pvalue")
        } else if (input$de_package == "edgeR") {
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
                      de_filter = input$de_filter,
                      subset = input$de_subset)
    })
    
    # render results table
    # will re-render baesd on user selections
    output$de_res_table <- renderDataTable(
        DT::datatable(formatted_res_render(), options = list(pageLength = 20))
    )
    
    # format results for plotting -----
    # plots cannot use data that is filtered by user selection and require the de_column!
    formatted_res_plots <- reactive({
        formatResults(de_res = req(de_res()),
                      logfc_threshold = input$logfc_threshold,
                      pvalue_threshold = input$pvalue_threshold,
                      logfc_col = col_names$logfc_col,
                      pvalue_col = col_names$pvalue_col,
                      de_column = TRUE,
                      de_filter = input$de_filter,
                      subset = FALSE)
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
                    row_names = input$heatmap_gene_names,
                    col_names = input$heatmap_sample_names)
    })
    
    plot_height <- reactive({
       de_num <- sum(formatted_res_plots()$isDE, na.rm = TRUE)
       de_num * 15 + 300
    })
    
    output$de_heatmap <- renderPlot({
         heatmap()
    }, height = function() {
        plot_height()
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
            paste0("diffEx-", input$de_package, "-results-", Sys.Date(), input$results_file_type)
        },
        
        content = function(file) {
            sep <- switch(input$results_file_type,
                          ".tsv" = "\t",
                          ".csv" = ",",
                          ".txt" = " ")
            write.table(formatted_res_render(), file, sep = sep)
        })
    
    output$download_ma <- downloadHandler(
        filename = function() {
            paste0("diffEx-", input$de_package, "-ma-plot-", Sys.Date(), ".pdf")
        },
        
        content = function(file) {
            ggsave(file, ma(), device = "pdf", width = input$ma_w, height = input$ma_h, units = "in")
        })
    
    output$download_volcano <- downloadHandler(
        filename = function() {
            paste0("diffEx-", input$de_package, "-volcano-", Sys.Date(), ".pdf")
        },
        
        content = function(file) {
            ggsave(file, volcano(), device = "pdf", width = input$volcano_w, height = input$volcano_h, units = "in")
        })
    
    output$download_heatmap <- downloadHandler(
        filename = function() {
            paste0("diffEx-", input$de_package, "-de-heatmap-", Sys.Date(), ".pdf")
        },
        
        content = function(file) {
            ggsave(file, heatmap(), device = "pdf", width = input$heatmap_w, height = input$heatmap_h, units = "in")
        })
}

# Run the application 
shinyApp(ui, server)