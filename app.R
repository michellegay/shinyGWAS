
library(shiny)
library(plotly)
library(magrittr)


# Allow file imports up to 1GB in size
options(shiny.maxRequestSize = 1000*1024^2)

########################################################
########################################################
########################################################

ui = fluidPage(
    
    # Multi-tab layout ----
    tabsetPanel(
        
        # Tab 1 ----
        tabPanel(title = "Phenotype & Covariate Data", fluid = TRUE,
                 
                 # Title ----
                 titlePanel("Upload Phenotype & Covariate Data"),
                 
                 # Row 1 ----
                 fluidRow(
                     column(4,
                            
                            # Horizontal line ----
                            tags$hr(),
                            
                            # Input: Select phenotype data file ----
                            fileInput(inputId = "phenoFile", 
                                      label = "Select .CSV file",
                                      multiple = TRUE,
                                      accept = c(".csv")),
                            
                            # Input: Select the id column ----
                            selectInput(inputId = "id",
                                        label = "ID",
                                        character(0),
                                        selectize = FALSE),
                            
                            # Input: Select the phenotype column ----
                            selectInput(inputId = "pheno",
                                        label = "Phenotype",
                                        character(0),
                                        selectize = FALSE),
                            
                            # Input: Select the sex column ----
                            selectInput(inputId = "sex",
                                        label = "Sex",
                                        character(0),
                                        selectize = FALSE),
                            
                            # Input: Select additional covariates ----
                            selectInput(inputId = "covars",
                                        label = "Other covariates",
                                        character(0),
                                        selectize = FALSE,
                                        multiple = TRUE),
                            
                            # Input: Submit button to create scanObj
                            actionButton(inputId = "phenoDo",
                                         label = "Submit"),
                            
                            ),
                     
                     column(8,
                            
                            # Output: show summary table for phenotype data
                            tableOutput(outputId = "phenoSummary")
                            
                            )
                     )
        ),
        
        
        # Tab 2 ----
        tabPanel("Genome Data", fluid = TRUE,
                 
                 # Title ----
                 titlePanel("Upload Genome Data"),
                 
                 # Row 1 ----
                 fluidRow(
                     column(4,
                            
                            # Horizontal line ----
                            tags$hr(),
                            
                            # Input: Select file type(s) for upload
                            radioButtons(inputId = "fileType",
                                         label = "Select file format(s)",
                                         choices = c("GDS" = "gds",
                                                     "PLINK" = "plink"),
                                         selected = "gds"),
                            
                            
                            conditionalPanel(
                                condition = "input.fileType == 'gds'",
                                
                                # Input: Select genome/SNP data file ----
                                fileInput(inputId = "gdsFile", 
                                          label = "Select GDS file",
                                          multiple = FALSE,
                                          accept = c(".gds"))
                            ),
                            
                            conditionalPanel(
                                condition = "input.fileType == 'plink'",
                                
                                # Input: Select genome/SNP data file ----
                                fileInput(inputId = "bedFile", 
                                          label = "Select PLINK .BED file",
                                          multiple = FALSE,
                                          accept = c(".bed")),
                                
                                # Input: Select genome/SNP data file ----
                                fileInput(inputId = "bimFile", 
                                          label = "Select PLINK .BIM file",
                                          multiple = FALSE,
                                          accept = c(".bim")),
                                
                                # Input: Select genome/SNP data file ----
                                fileInput(inputId = "famFile", 
                                          label = "Select PLINK .FAM file",
                                          multiple = FALSE,
                                          accept = c(".fam")),
                                
                            ),
                            
                            
                            # Input: Is genome data imputed? ----
                            checkboxInput(inputId = "isImpute", 
                                          label = "Check box if data has been imputed", 
                                          value = FALSE),
                            
                            # Input: Submit button to create genObj
                            actionButton(inputId = "genoDo",
                                         label = "Submit"),
                            
                            ),
                     
                     column(8,
                            
                            # Output: show summary table for genome data
                            tableOutput(outputId = "genoSummary")
                            
                            )
                 )
        )
    )
)



server <- function(input, output, session) {
    
    source("./shinyGWAS_fns.R")
    
    ##############################
    ## PHENOTYPE/COVARIATE DATA ##
    ##############################
    
    # Save csv data object
    data <- reactive({
        req(input$phenoFile)
        inPheno <- input$phenoFile
        read.csv(inPheno$datapath)
    })
    
    # Request ID variable
    observeEvent(data(), {
        updateSelectInput(session = session, inputId = "id", 
                          choices = c(" " = "", names(data())))
    })
    
    # Request phenotype variable
    observeEvent(data(), {
        updateSelectInput(session = session, inputId = "pheno", 
                          choices = c(" " = "", names(data())))
    })
    
    # Request sex variable
    observeEvent(data(), {
        updateSelectInput(session = session, inputId = "sex", 
                          choices = c(" " = "", names(data())))
    })
    
    # Request additional covariates
    observeEvent(data(), {
        updateSelectInput(session = session, inputId = "covars", 
                          choices = c(" " = "", names(data())), 
                          selected = "")
    })
    
    
    # Save data to ScanAnnotationDataFrame object 
    scanAnnot <- eventReactive(input$phenoDo, {
        # Check required inputs
        req(input$phenoFile)
        validate(
            need(input$id != "", "Select ID variable"),
            need(input$pheno != "", "Select phenotype variable"),
            need(input$sex != "", "Select sex variable"))
        
        scan <- scanObj(data = data(),
                        id = input$id,
                        pheno = input$pheno,
                        sex = input$sex,
                        covars = input$covars)
    })
    
    # Output pheno data summary
    output$phenoSummary <- renderTable({
        pData(scanAnnot()) %>% head
    })
 
    
    #################
    ## GENOME DATA ##
    #################
    
    gdsPath <- reactive({
        req(input$gdsFile)
        input$gdsFile$datapath %>% return
    })
    
    plinkPath <- reactive({
        req(input$bedFile)
        req(input$bimFile)
        req(input$famFile)
        path <- list("bed" = input$bedFile$datapath,
                     "bim" = input$bimFile$datapath,
                     "fam" = input$famFile$datapath)
        plinkToGds(path) %>% return
    })

    # Output genome data summary
    output$genoSummary <- renderTable({
        
        if (input$fileType == "gds"){
            gData <- openGDS(gdsPath)
        } else if (input$fileType == "plink"){
            gData <- openGDS(plinkPath)
        }
        
        data.frame(snpID = getSnpID(gData),
                   chromosome = getChromosome(gData),
                   position = getPosition(gData))
    })
    
    
    
}

# Run the application 
shinyApp(ui = ui, server = server)