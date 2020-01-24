
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
                                        choices = character(0),
                                        selectize = FALSE),
                            
                            # Input: Select the phenotype column ----
                            selectInput(inputId = "pheno",
                                        label = "Phenotype",
                                        choices = character(0),
                                        selectize = FALSE),
                            
                            # Input: Select the sex column ----
                            selectInput(inputId = "sex",
                                        label = "Sex",
                                        choices = character(0),
                                        selectize = FALSE),
                            
                            # Input: Select additional covariates ----
                            selectInput(inputId = "covars",
                                        label = "Other covariates",
                                        choices = character(0),
                                        selectize = FALSE,
                                        multiple = TRUE),
                            
                            tags$hr(),
                            
                            # Input: Submit button to create scanObj ----
                            actionButton(inputId = "phenoDo",
                                         label = "Submit")
                            
                            ), #END column
                     
                     column(8,
                            
                            # Output: A busy status indicator
                            add_busy_spinner(spin = "fading-circle"),
                            
                            # Output: show summary table for phenotype data ----
                            tableOutput(outputId = "phenoSummary")
                            
                            ) #END column
                     ) #END row
        ), #END tab
        
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
                            ), #END panel
                            
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
                                
                            ), #END panel
                            
                            tags$hr(),
                            
                            # Input: Submit button to create genObj ----
                            actionButton(inputId = "genoDo",
                                         label = "Submit")
                            
                            ), #END column
                     
                     column(8,
                            
                            # Output: A busy status indicator
                            add_busy_spinner(spin = "fading-circle"),
                            
                            # Output: show summary table for genome data ----
                            tableOutput(outputId = "chromTable")
                            
                            ) #END column
                 ) #END row
        ), #END tab
        
        ################################################
        ## Tab 3 Configure Association Test           ##
        ################################################
        
        tabPanel("Association", fluid = TRUE,
                 
                 # Title ----
                 titlePanel("Configure Settings for Association Testing"),
                 
                 fluidRow(
                     column(4,
                            tags$hr(),
                            
                            # Input: Is genome data imputed? ----
                            checkboxInput(inputId = "isImpute", 
                                          label = "Check box if data has been imputed", 
                                          value = FALSE),
                            
                            tags$hr(),
                            
                            # Input: Select datatype of Phenotype ----
                            radioButtons(inputId = "phenoDatType",
                                         label = "Select phenotype data-type",
                                         choices = c("Binary" = "bin",
                                                     "Ordinal" = "ord",
                                                     "Continuous" = "cont",
                                                     "Categorical/Nominal" = "cat"),
                                         selected = "bin"),
                            
                            tags$hr(),
                            
                            # INPUT: Select association model ----
                            radioButtons(inputId = "assocTest",
                                         label = "Select association model",
                                         choices = c("Logistic Regression" = "logistic", 
                                                     "Logistic Mixed Model" = "logMM",
                                                     "Linear Regression" = "linear",
                                                     "Linear Mixed Model" = "linMM",
                                                     "Ordinal Logit" = "ordLog",
                                                     "(-_-)_/¯" = "idk")),
                            
                            tags$hr(),
                            tags$hr(),
                            tags$hr(),
                            
                            # Additional inputs/parameter config for linear/logistic mixed models ----
                            conditionalPanel(
                                condition = "input.assocTest == 'linMM' || input.assocTest == 'logMM'",
                                
                                # Input: Select random effects (e.g. kinship matrix) ----
                                fileInput(inputId = "randEffect", 
                                          label = "Select random effects file",
                                          multiple = FALSE,
                                          accept = c(".csv")) #Not sure what to accept here
                                
                            ), #END panel
                            
                            h4("Additional model-specific inputs/parameter configs go here"),
                            
                            tags$hr(),
                            tags$hr(),
                            tags$hr(),
                            
                            # Input: Submit button to run association analysis ----
                            actionButton(inputId = "assocDo",
                                         label = "Run Association")
                            
                            ), #END column
                     
                     column(8,
                            
                            # Output: A busy status indicator
                            add_busy_spinner(spin = "fading-circle"),
                            
                            # Output: show summary table for phenotype data ----
                            tableOutput(outputId = "gwasSummary")
                            
                     ) #END column
                 ) #END row
        )#END tab
     ) #END tabset
) #END ui



server <- function(input, output, session) {
    
    source("./shinyGWAS_fns.R")
    
    ################################################
    ## Tab 1 Uploading Phenotype & Covariate Data ##
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
    ## Tab 2 Uploading Genome Data                ##
    ################################################
    
    # Save genome data to GenotypeData object on "Submit"
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
        
        # validate(
        #     need(isS4(scanAnnot) == TRUE, "Upload phenotype data."))
        openGDS(path, scanAnnot())
    })
    
    # Output: excerpt of table to show reading in correcctly.
    output$chromTable <- renderTable({
        data.frame(chr = getChromosome(gData())[1:10],
                   id = getSnpID(gData())[1:10],
                   pos = getPosition(gData())[1:10])
    })
    
    ################################################
    ## Tab 3 Configuring Association Test Settigs ##
    ################################################
    
    # Update association test options based on data type selection
    observeEvent(input$phenoDatType, {
        req(input$phenoDatType)
        
        if (input$phenoDatType == "bin"){
            choices <- c("Logistic Regression" = "logistic", 
                         "Logistic Mixed Model" = "logMM",
                         "Linear Mixed Model" = "linMM")
            selected <- "logistic"
            
        } else if (input$phenoDatType == "ord"){
            choices <- c("Linear Regression" = "linear",
                         "Linear Mixed Model" = "linMM",
                         "Ordinal Logit" = "ordLog")
            selected <- "linear"
            
        } else if (input$phenoDatType == "cont"){
            choices <- c("Linear Regression" = "linear",
                         "Linear Mixed Model" = "linMM")
            selected <- "linear"
            
        } else if (input$phenoDatType == "cat"){
            choices <- c("(-_-)_/¯" = "idk")
            selected <- "idk"
        }
        
        updateRadioButtons(session,
                           inputId = "assocTest",
                           choices = choices,
                           selected = selected)
    })
    
    # Save genome data to GenotypeData object on "Submit"
    assocData <- eventReactive(input$assocDo, {
        req(gData())

        if (input$assocTest == "logistic" | input$assocTest == "linear"){
            gwasLogLin(gData(), input$assocTest)

        } else if (input$assocTest == "linMM" | input$assocTest == "logMM"){
            gwasMixed(gData())

        } else {
            NULL
        }
    })
    
    # Output pheno data summary
    output$gwasSummary <- renderTable({
        assocData() %>% head
    })
    
} #END server

# Run the application 
shinyApp(ui = ui, server = server)


