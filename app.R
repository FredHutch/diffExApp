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
    fluidRow(
        column(width = 4,
               h4("Specify which columns align with conditions"),
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
               br(),
               actionButton("apply", "Apply")),
        column(width = 8,
               h4("Sample Table"),
               dataTableOutput("sample_table"),
               br(),
               hr(),
               h4("Differential Expression Results"),
               em("Select replicates, name conditions, select desired analysis method, and click apply"),
               br(),
               dataTableOutput("deseq_res_table"))
    ) # end fluidRow
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
    
    observeEvent(input$apply, {
        output$deseq_res_table <- renderDataTable(
            de_res_table()
        )
    })
}

# Run the application 
shinyApp(ui = ui, server = server)
