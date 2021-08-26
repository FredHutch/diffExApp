# SET UP -----
# load libraries
library(shiny)
library(tidyverse)

# read in data
data <- read.table("~/Documents/work/dataCore/shiny/diffEx/data/counts.csv", sep = ",", header = TRUE)

# check cols
counts_candidate_f <- function(x) {
    if (is.numeric(data[[x]])) {
        return(TRUE)}
    return(FALSE)
}

info_candidate_f <- function(x) {
    if (is.character(data[[x]])) {
        return(TRUE)}
    return(FALSE)
}

counts_cols <- names(data)[sapply(names(data), counts_candidate_f)]
info_cols <- names(data)[sapply(names(data), info_candidate_f)]


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
                           choices = c("edgeR", "DESeq2", "limma/voom"))
               
        )
    ) # end fluidRow
) # end fluidPage

# SERVER -----
server <- function(input, output) {
    
    # select condition colnames -----
    # condition 1
    output$condition1_selector <- renderUI({
        selectInput("select_condition1",
                    "Select replicates: ",
                    counts_cols,
                    multiple = TRUE,
                    selectize= TRUE)
    })
    # condition 2
    output$condition2_selector <- renderUI({
        selectInput("select_condition2",
                    "Select replicates: ",
                    counts_cols,
                    multiple = TRUE,
                    selectize= TRUE)
    })
}

# Run the application 
shinyApp(ui = ui, server = server)
