
library(shiny)
library(plotly)
library(magrittr)
library(GENESIS)
library(SNPRelate)
library(lme4)
library(GWASTools)
library(GWASdata)
library(qqman)

# Allow file imports up to 1GB in size

options(shiny.maxRequestSize = 1000*1024^2)


########################################################
## Read in and save phenotype and covariate data into ##
## a ScanAnnotationDataFrame object                   ##
########################################################

scanObj <- function(data, id, pheno, gender, covars){
    
    dat <- data.frame(scanID = data[,id],
                      phenotype = data[,pheno],
                      sex = data[,gender],
                      data[,covars])
    
    scanAnnot <- ScanAnnotationDataFrame(dat)
    scanAnnot$sex <- as.factor(scanAnnot$sex)
    scanAnnot %>% return
}

########################################################
## Read in and save genotype data into a              ##
## GenotypeData object                                ##
########################################################

genObj <- function(type, path, scan) {
    if (type == "gds") {
        NULL
        
    } else if (type == "PLINK") {
        bed <- INPUT
        bim <- INPUT
        fam <- INPUT
        
        snpgdsBED2GDS(bed.fn = bed,
                      bim.fn = bim,
                      fam.fn = fam,
                      out.gdsfn = path)
        
    } else if (type == "txt") {
        NULL
    }
    
    geno <- GdsGenotypeReader(filename = path)
    genoData <- GenotypeData(geno, scanAnnot = scan)
}

########################################################
## Perform association test by Linear Mixed Model     ##
########################################################

gwasMixed <- function(genoData, covars) {
    
    # Include all covariates, remove the scanID and outcome var's
    # PROBS BETTER TO HAVE A SELECTION TO SELECT/DESELECT COVARs
    # covars <- (genoData %>% getScanAnnotation %>% getVariableNames)[-c(1,2)]
    
    null.model <- genoData %>% getScanAnnotation %>% 
        fitNullModel(outcome = "phenotype",
                     covars = covars,
                     family = gaussian)
    
    iter <- GenotypeBlockIterator(genoData = genoData,
                                  snpBlock = 10000)
    
    assoc <- assocTestSingle(iter, 
                             null.model = null.model)
    
    # Only returning some of the results from the null model
    list("Null Model" = null.model[c(1,3,5)], 
         "Association" = assoc) %>% return
}

########################################################
## Perform association test by Linear or Logistic     ##
## Regression                                         ##
########################################################

gwasLogLin <- function(genoData, model.type, covars) {
    # model.type = "linear" or "logistic"
    
    # Include all covariates, remove the scanID and outcome var's
    # covars <- (genoData %>% getScanAnnotation %>% getVariableNames)[-c(1,2)]
    
    
    chr <- genoData %>% getChromosome
    
    assoc.list <- unique(chr) %>% lapply(function(x) {
        ## Y chromosome only includes males, cannot have sex as covariate
        ## copied from GWASTools vignette
        ## NEED TO IMPROVE THIS, MAY NOT HAVE SEX AS COVAR IN FIRST PLACE
        covar <- ifelse(x == 25, yes = covars[-1], no = covars)
        start <- which(chr == x)[1]
        end <- start + (which(chr == x) %>% length) -1
        
        assocRegression(genoData = genoData,
                        outcome = "phenotype",
                        model.type = model.type,
                        covar = covar,
                        snpStart = start,
                        snpEnd = end)
        
    })
    
    assoc <- do.call(rbind, assoc.list)
    assoc %>% return
}

########################################################
########################################################
########################################################

ui = fluidPage(
    
    # Multi-tab layout ----
    tabsetPanel(
        
        # Tab 1 ----
        tabPanel(title = "Import Data", fluid = TRUE,
                 
                 # Title ----
                 titlePanel("Input Data"),
                 
                 # Row 1 ----
                 fluidRow(
                     
                     # R1, Col 1 ----
                     column(8,
                            
                            column(6,
                            
                                # Sub-heading ----
                                h3("Phenotype and Covariate Data:"),
                                
                                # Horizontal line ----
                                tags$hr(),
                                
                                # Input: Select phenotype data file ----
                                fileInput(inputId = "phenoFile", 
                                          label = "Select file(s)",
                                          multiple = TRUE,
                                          accept = c(".csv"),
                                          placeholder = ".csv"),
                                
                                column(4,
                                
                                    # Input: Select the id column ----
                                    selectInput(inputId = "idIdx",
                                                label = "ID",
                                                character(0),
                                                selectize = FALSE)
                                ),
                                
                                column(4,
                                
                                    # Input: Select the phenotype column ----
                                    selectInput(inputId = "phenoIdx",
                                                label = "Phenotype",
                                                character(0),
                                                selectize = FALSE)
                                ),
                                
                                column(4,
                                
                                    # Input: Select the sex column ----
                                    selectInput(inputId = "sexIdx",
                                                label = "Sex",
                                                character(0),
                                                selectize = FALSE)
                                ),
                        
                            # Input: Select additional covariates ----
                            selectInput(inputId = "covarsIdx",
                                        label = "Other covariates",
                                        character(0),
                                        selectize = FALSE,
                                        multiple = TRUE),
                            
                            # Action: save data to ScanAnnotationDataFrame object
                           # actionButton()
                            
                            ),
                           
                           column(6,
                                  
                                  # Sub-heading ----
                                  h3("Additional Covariate Data:"),
                                  
                                  # Horizontal line ----
                                  tags$hr(),
                                  
                                  # Input: Select phenotype data file ----
                                  fileInput(inputId = "covarFile", 
                                            label = "Select file(s)",
                                            multiple = TRUE,
                                            accept = c(".csv"),
                                            placeholder = ".csv"),
                                  
                                  
                                  )
                           
                           ),
                     
                     
                     # R1, Col 2 ----
                     column(4,
                             
                             # Sub-heading ----
                             h3("Genome Data"),
                             
                             # Horizontal line ----
                             tags$hr(),
                             
                             # Input: Select genome file type ----
                             radioButtons(inputId = "genomeType",
                                          "Select file type(s)",
                                          choices = c(GDS = "gds",
                                                      CSV = "csv",
                                                      "BED, BIM, FAM" = "PLINK"),
                                          selected = "gds"),
                             
                             # Input: Select genome/SNP data file ----
                             fileInput(inputId = "genomeFile", 
                                       label = "Select file(s)",
                                       multiple = TRUE,
                                       accept = c(".csv",
                                                  ".bed",
                                                  ".bim",
                                                  ".fam",
                                                  ".gds")),
                             
                             # Input: Is genome data imputed? ----
                             checkboxInput(inputId = "isImpute", 
                                           label = "Check box if data has been imputed", 
                                           value = FALSE)
                            ),
                     ),
                 
                 
                 # Row 2 ----
                 fluidRow(
                     
                     # R2, Col 1 ----
                     column(12,
                            
                            # Output text ----
                            textOutput("selected"))
                     ),
                 ),
        
        
        # Tab 2 ----
        tabPanel("Plot", fluid = TRUE)
        
        )
    )





server <- function(input, output, session) {
    
    ##############################
    ## PHENOTYPE/COVARIATE DATA ##
    ##############################
    
    data <- reactive({
        inPheno <- input$phenoFile
        if (is.null(inPheno)) return(NULL)
        read.csv(inPheno$datapath)
    })
    
    observeEvent(data(), {
        updateSelectInput(session = session, inputId = "idIdx", 
                          choices = names(data()))
    })
    
    observeEvent(data(), {
        updateSelectInput(session = session, inputId = "phenoIdx", 
                          choices = names(data()))
    })
    
    observeEvent(data(), {
        updateSelectInput(session = session, inputId = "sexIdx", 
                          choices = c("N/A" = "", names(data())))
    })
    
    observeEvent(data(), {
        updateSelectInput(session = session, inputId = "covarsIdx", 
                          choices = c("N/A" = "", names(data())))
    })
    
    
    
    
    output$selected <- renderText({
        req(data())
        paste0("You've selected ", input$idIdx,
               " as the ID column.")
    })
    
}

# Run the application 
shinyApp(ui = ui, server = server)
