
library(shiny)
library(plotly)
library(magrittr)
library(shinybusy)



# Allow file imports up to 1GB in size
options(shiny.maxRequestSize = 1000*1024^2)

ui = fluidPage(
    tabsetPanel(
        
        ################################################
        ## Tab 1 Uploading Phenotype & Covariate data ##
        ################################################
        
        tabPanel(title = "Phenotype & Covariate Data", fluid = TRUE,
                 
                 # Title ----
                 titlePanel("Upload Phenotype & Covariate Data"),
                 
                 fluidRow(
                     column(4,
                            tags$hr(),
                            
                            # Input: Select phenotype data file ----
                            fileInput(inputId = "phenoFile", 
                                      label = "Select .CSV file",
                                      multiple = TRUE,
                                      accept = c(".csv")),
                            
                            tags$hr(),
                            
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
                            
                            # Input: Submit button to create scanObj ----
                            actionButton(inputId = "phenoDo",
                                         label = "Submit")
                            
                            ),
                     
                     column(8,
                            
                            # Output: show summary table for phenotype data ----
                            tableOutput(outputId = "phenoSummary")
                            
                            )
                     )
        ),
        
        ################################################
        ## Tab 2 Uploading Genome data                ##
        ################################################
        
        tabPanel("Genome Data", fluid = TRUE,
                 
                 # Title ----
                 titlePanel("Upload Genome Data"),
                 
                 fluidRow(
                     column(4,
                            tags$hr(),
                            
                            # Input: Select file format(s) for upload ----
                            radioButtons(inputId = "fileType",
                                         label = "Select file format(s)",
                                         choices = c("GDS" = "gds",
                                                     "PLINK" = "plink"),
                                         selected = "gds"),
                            
                            tags$hr(),
                            
                            # Show if file format is GDS ----
                            conditionalPanel(
                                condition = "input.fileType == 'gds'",
                                
                                # Input: Select genome/SNP data file ----
                                fileInput(inputId = "gdsFile", 
                                          label = "Select GDS file",
                                          multiple = FALSE,
                                          accept = c(".gds"))
                            ),
                            
                            # Show if file format is PLINK ----
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
                                          accept = c(".fam"))
                                
                            ),
                            
                            
                            # Input: Is genome data imputed? ----
                            checkboxInput(inputId = "isImpute", 
                                          label = "Check box if data has been imputed", 
                                          value = FALSE),
                            
                            # Input: Submit button to create genObj ----
                            actionButton(inputId = "genoDo",
                                         label = "Submit")
                            
                            ),
                     
                     column(8,
                            
                            # Output: A busy status indicator
                            add_busy_spinner(spin = "fading-circle"),
                            
                            # Output: show summary table for genome data ----
                            tableOutput(outputId = "chromTable")
                            
                            )
                 )
        )
    )
)



server <- function(input, output, session) {
    
    source("./shinyGWAS_fns.R")
    
    ################################################
    ## Tab 1 Uploading Phenotype & Covariate data ##
    ################################################
    
    # Read csv data object
    sData <- reactive({
        req(input$phenoFile)
        inPheno <- input$phenoFile
        read.csv(inPheno$datapath)
    })
    
    # Request ID variable
    observeEvent(sData(), {
        updateSelectInput(session = session, inputId = "id", 
                          choices = c(" " = "", names(sData())))
    })
    
    # Request phenotype variable
    observeEvent(sData(), {
        updateSelectInput(session = session, inputId = "pheno", 
                          choices = c(" " = "", names(sData())))
    })
    
    # Request sex variable
    observeEvent(sData(), {
        updateSelectInput(session = session, inputId = "sex", 
                          choices = c(" " = "", names(sData())))
    })
    
    # Request additional covariates
    observeEvent(sData(), {
        updateSelectInput(session = session, inputId = "covars", 
                          choices = c(" " = "", names(sData())), 
                          selected = "")
    })
    
    
    # Save data to ScanAnnotationDataFrame object on "Submit"
    scanAnnot <- eventReactive(input$phenoDo, {
        # Check required inputs
        req(input$phenoFile)
        validate(
            need(input$id != "", "Select ID variable."),
            need(input$pheno != "", "Select phenotype variable."),
            need(input$sex != "", "Select sex variable.")
            )
        
        scan <- scanObj(data = sData(),
                        id = input$id,
                        pheno = input$pheno,
                        sex = input$sex,
                        covars = input$covars)
    })
    
    # Output pheno data summary
    output$phenoSummary <- renderTable({
        pData(scanAnnot()) %>% head
    })
 
    
    ################################################
    ## Tab 2 Uploading Genome data                ##
    ################################################
    
    # Save file path for gds file
    gData <- eventReactive(input$genoDo, {
        req(input$fileType)
        
        # Close files so instruction can be re-executed
        showfile.gds(closeall = TRUE)
        
        
        if (input$fileType == "gds"){
            req(input$gdsFile)
            path <- input$gdsFile$datapath %>% 
                as.character
            
        } else if (input$fileType == "plink"){
            req(input$bedFile)
            req(input$bimFile)
            req(input$famFile)
            paths <- list("bed" = input$bedFile$datapath,
                         "bim" = input$bimFile$datapath,
                         "fam" = input$famFile$datapath)
            path <- plinkToGds(paths)
        }
        
        validate(
            need(isS4(scanAnnot) == TRUE, "Upload phenotype data."))
        openGDS(path, scanAnnot)
    })
    
    # Output: excerpt of table to show reading in correcctly.
    output$chromTable <- renderTable({
        data.frame(chr = getChromosome(gData())[1:10],
                   id = getSnpID(gData())[1:10],
                   pos = getPosition(gData())[1:10])
    })

}

# Run the application 
shinyApp(ui = ui, server = server)


