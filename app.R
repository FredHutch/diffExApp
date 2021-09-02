# SET UP -----
# load libraries
library(shiny)
library(tidyverse)
library(DESeq2)
library(DT)

# source scripts
source("./helpers.R")

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
gene_cols <- names(data)[sapply(names(data), gene_candidate_f)]

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
                                         choices = c("edgeR", "DESeq2", "limma/voom")),
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
                                   value = 0,
                                   step = .01),
                       checkboxInput("fdr",
                                     "Use False Discovery Rate",
                                     TRUE),
                       sliderInput("logfc_threshold",
                                   "Set log fold change threshold",
                                   min = 0,
                                   max = 10,
                                   value = 0,
                                   step = .5),
                       checkboxInput("fdr",
                                     "Use adjusted P-value (False Discovery Rate)",
                                     TRUE),
                       checkboxInput("de_column",
                                     "Show differentially expressed column",
                                     FALSE),
                       checkboxInput("de_filter",
                                     "Filter table to only show differentially expressed genes",
                                     FALSE)),
                # tab 2 main panel
                mainPanel(
                    dataTableOutput("deseq_res_table"))
                ) # end tab 2 sidebar layout
        ) # end tabpanel 2
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
    
    # Deseq2 differential expression -----
    reactive_sample_table <- reactive({
        createSampleMatrix(condition1_name = input$condition1_name,
                           condition2_name = input$condition2_name,
                           condition1_selected = input$condition1_selected,
                           condition2_selected = input$condition2_selected)
    })
    
    output$sample_table <- renderDataTable(
        reactive_sample_table() %>%
            rownames_to_column(var = "sample") %>%
            select(sample, everything())
    )
    
    de_res_table <- reactive({
        deseq(data = data,
              condition1_name = input$condition1_name,
              condition2_name = input$condition2_name,
              condition1_selected = input$condition1_selected,
              condition2_selected = input$condition2_selected,
              gene_col = gene_col)
    })
    
    de_res_table_filtered <- reactive({
        showDiffEx(de_res = de_res_table(),
                   pvalue_threshold = input$pvalue_threshold,
                   logfc_threshold = input$logfc_threshold,
                   fdr = input$fdr,
                   de_column = input$de_column,
                   de_filter = input$de_filter)
    })
    
    observeEvent(input$apply, {
        output$deseq_res_table <- renderDataTable(
            de_res_table_filtered()
        )
    })
}

# Run the application 
shinyApp(ui = ui, server = server)
