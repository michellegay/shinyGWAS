
library(shiny)
library(shinybusy)
library(magrittr)
library(GENESIS)
library(SNPRelate)
library(lme4)
library(GWASTools)
library(GWASdata)
library(qqman)

library(summarytools)
library(GGally)
library(shinyFiles)
library(DT)
library(plotly)
library(manhattanly)
library(qqman)
library(dplyr)

# Allow file imports up to 3GB in size
options(shiny.maxRequestSize = 3000*1024^2)

ui = fluidPage(
    br(),
    tabsetPanel(
        
        ################################################
        ## Tab 1 Uploading Phenotype & Covariate data ##
        ################################################
        
        tabPanel(title = "Upload Data", fluid = TRUE,
                 
                 titlePanel("Upload Data"),
                 
                 fluidRow( 
                     column(4,
                            tags$hr(),
                            
                            tabsetPanel(
                                tabPanel(title = "Genome & Phenotype", fluid = TRUE,
                                         
                                         br(),
                                         
                                         # Input: Is genome data imputed? ----
                                         checkboxInput(inputId = "isImpute", 
                                                       label = "Check box if data has been imputed", 
                                                       value = FALSE),
                                         
                                         tags$hr(),
                                         
                                         # Input: Select genome file format(s) ----
                                         radioButtons(inputId = "fileType",
                                                      label = "Select file format",
                                                      choices = c("GDS" = "gds",
                                                                  "VCF" = "vcf",
                                                                  "PLINK" = "plink",
                                                                  "IMPUTE2" = "IMPUTE2",
                                                                  "BEAGLE" = "BEAGLE",
                                                                  "MaCH" = "MaCH")),
                                         
                                         tags$hr(),
                                         
                                         # Show only if data is imputed ----
                                         conditionalPanel(
                                             condition = "input.isImpute == 1",
                                             
                                             radioButtons(inputId = "isDosage",
                                                          label = "Indicate how imputed data is represented",
                                                          choices = c("Dosage" = TRUE,
                                                                      "Genotype Probabilities" = FALSE),
                                                          selected = TRUE),
                                             
                                             tags$hr()
                                             
                                         ), #END panel
                                         
                                         # Show only if file format is GDS ----
                                         conditionalPanel(
                                             condition = "input.fileType == 'gds'",
                                             
                                             # Input: Select genome/SNP file ----
                                             fileInput(inputId = "gdsFile", 
                                                       label = "Select GDS file",
                                                       multiple = FALSE,
                                                       accept = c(".gds"))
                                         ), #END panel
                                         
                                         # Show only if file format is VCF ----
                                         conditionalPanel(
                                             condition = "input.fileType == 'vcf'",
                                             
                                             # Input: Select genome/SNP file ----
                                             fileInput(inputId = "vcfFile", 
                                                       label = "Select VCF file(s)",
                                                       multiple = TRUE,
                                                       accept = c(".vcf"))
                                         
                                          ), #END panel
                                         
                                         # Show only if file format is PLINK ----
                                         conditionalPanel(
                                             condition = "input.fileType == 'plink'",
                                             
                                             # Input: Select .bed file ----
                                             fileInput(inputId = "bedFile", 
                                                       label = "Select PLINK .BED file",
                                                       multiple = FALSE,
                                                       accept = c(".bed")),
                                             
                                             # Input: Select .bim file ----
                                             fileInput(inputId = "bimFile", 
                                                       label = "Select PLINK .BIM file",
                                                       multiple = FALSE,
                                                       accept = c(".bim")),
                                             
                                             # Input: Select .fam file ----
                                             fileInput(inputId = "famFile", 
                                                       label = "Select PLINK .FAM file",
                                                       multiple = FALSE,
                                                       accept = c(".fam"))
                                             
                                         ), #END panel
                                         
                                         # Show only if file format is IMPUTE2 ----
                                         conditionalPanel(
                                             condition = "input.fileType == 'IMPUTE2'",
                                             
                                             # Input: Select .gens file ----
                                             fileInput(inputId = "imp2Imps1", 
                                                       label = "Select IMPUTE2 .gens file(s)",
                                                       multiple = TRUE,
                                                       accept = c(".gens")),
                                             
                                             # Input: Select .samples file ----
                                             fileInput(inputId = "imp2Imps2", 
                                                       label = "Select IMPUTE2 .samples file(s)",
                                                       multiple = TRUE,
                                                       accept = c(".samples"))
                                             
                                         ), #END panel
                                         
                                         # Show only if file format is BEAGLE ----
                                         conditionalPanel(
                                             condition = "input.fileType == 'BEAGLE'",
                                             
                                             # Input: Select .grobs or .dose file ----
                                             fileInput(inputId = "beagImps1", 
                                                       label = "Select BEAGLE .grobs or .dose file(s)",
                                                       multiple = TRUE,
                                                       accept = c(".grobs",
                                                                  ".dose")),
                                             
                                             # Input: Select .markers file ----
                                             fileInput(inputId = "beagImps2", 
                                                       label = "Select BEAGLE .markers file(s)",
                                                       multiple = TRUE,
                                                       accept = c(".markers"))
                                             
                                         ), #END panel
                                         
                                         # Show only if file format is MaCH ----
                                         conditionalPanel(
                                             condition = "input.fileType == 'MaCH'",
                                             
                                             # Input: Select .mlprob or .mldose file ----
                                             fileInput(inputId = "machImps1", 
                                                       label = "Select MaCH .mlprob/.mldose file(s)",
                                                       multiple = TRUE,
                                                       accept = c(".mlprob",
                                                                  ".mldose")),
                                             
                                             # Input: Select .mlinfo file ----
                                             fileInput(inputId = "machImps2", 
                                                       label = "Select MaCH .mlinfo file(s)",
                                                       multiple = TRUE,
                                                       accept = c(".mlinfo")),
                                             
                                             # Input: Select additional data file ----
                                             # NOT SURE WHAT THIS FILE IS??
                                             fileInput(inputId = "machImps3", 
                                                       label = "Select file(s) with columns 'SNP' and 'position'",
                                                       multiple = TRUE,
                                                       accept = c(".csv"))
                                             
                                         ), #END panel
                                         
                                ), #END tab
                                
                                
                                tabPanel(title = "Covariates", fluid = TRUE,
                                         
                                         br(),
                                         
                                         conditionalPanel(
                                             condition = "input.fileType == 'vcf'",
                                             
                                             tags$p(paste("VCF files must be accompanied by",
                                                          "a CSV file containing the fields:",
                                                          "sample ID, phenotype and sex. The",
                                                          "sample ID must correspond to the",
                                                          "ID used in the VCF file.")),
                                             br(),
                                             
                                             # Input: Select phenotype data file ----
                                             fileInput(inputId = "phenoFile", 
                                                       label = "Select sample data file (.csv)",
                                                       multiple = FALSE,
                                                       accept = c(".csv")),
                                             
                                             tags$hr(),
                                             
                                             tags$p(paste("Select the column corresponding",
                                                          "to each of the specified fields below.")),
                                             br(),
                                             
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
                                             
                                             tags$hr()
                                             
                                         ), #END panel
                                         
                                         tags$p(paste("Select any additional covariates",
                                                      "to be included in association analyses.",
                                                      "There must be an ID column selected",
                                                      "which corresponds to the ID",
                                                      "in the genotype data. May include phenotype",
                                                      "data, if so, select the relevant column below.")),
                                         br(),
                                         
                                         # Input: Select phenotype data file ----
                                         fileInput(inputId = "covarFile", 
                                                   label = "Select covariate data file (.csv)",
                                                   multiple = TRUE,
                                                   accept = c(".csv")),
                                         
                                         # Input: Select the id column ----
                                         selectInput(inputId = "id2",
                                                     label = "ID",
                                                     choices = character(0),
                                                     selectize = FALSE),
                                         
                                         # Input: Select the phenotype column ----
                                         selectInput(inputId = "pheno2",
                                                     label = "Phenotype",
                                                     choices = character(0),
                                                     selectize = FALSE),
                                         
                                         # Input: Select additional covariates ----
                                         selectInput(inputId = "moreCovars",
                                                     label = "Select additional covariates",
                                                     choices = c("N/A" = ""),
                                                     selectize = FALSE,
                                                     multiple = TRUE,
                                                     selected = ""),
                                         
                                         # Input: Submit button to upload additional data ----
                                         actionButton(inputId = "genoDo",
                                                      label = "Upload genome data"),
                                         
                                 ) #END tab
                            ), #END tabset
                            
                            tags$hr(),
                            
                            # Output: A busy status indicator
                            add_busy_spinner(spin = "fading-circle")
                            
                            ), #END column
                     
                     column(8,
                            
                            # Make correlation plot dynamic to size of window
                            # Code from https://stackoverflow.com/a/40539526
                            tags$head(
                            tags$script(
                            '$(document).on("shiny:connected", function(e) {
                            Shiny.onInputChange("innerWidth", window.innerWidth);
                            });
                            $(window).resize(function(e) {
                            Shiny.onInputChange("innerWidth", window.innerWidth);
                            });
                            ')),
                            
                            tags$hr(),
                            
                            tabsetPanel(
                                tabPanel(title = "Descriptive Statistics", fluid = TRUE,
                                         
                                         # Output: Summary statistics table
                                         htmlOutput(outputId = "phenoSummary")
                                         
                                         ), #END tab
                                    
                                tabPanel(title = "Pair-wise Comparisons", fluid = TRUE,
                                         
                                         br(), br(),
                                         
                                         # Output: Pair-wise comparisons
                                         plotOutput(outputId = "heatmap")
                                         
                                         ), #END tab
                                
                                tabPanel(title = "Files", fluid = TRUE,
                                         
                                         br(), br(),
                                         
                                         # Output: show the files that have been uploaded ----
                                         DT::dataTableOutput(outputId = "filesUploaded")
                                         
                                         )
                            ) #END tabset
                            ) #END column
                     ) #END row
        ), #END tab
        
        ################################################
        ## Tab 2 Configure Association Test           ##
        ################################################
        
        tabPanel("Association", fluid = TRUE,
                 
                 # Title ----
                 titlePanel("Association Test Settings"),
                 
                 fluidRow(
                     column(4,
                            tags$hr(),
                            
                            # Input: Select datatype of Phenotype ----
                            radioButtons(inputId = "phenoDatType",
                                         label = "Select phenotype data-type",
                                         choices = c("Binary" = "bin",
                                                     "Quantitative" = "quant",
                                                     "Categorical" = "cat"),
                                         selected = "bin"),
                            
                            tags$hr(),
                            
                            # INPUT: Select association model ----
                            radioButtons(inputId = "assocTest",
                                         label = "Select association model",
                                         choices = c("Logistic Regression" = "logistic", 
                                                     # "Logistic Mixed Model" = "logMM",
                                                     "Linear Regression" = "linear",
                                                     # "Linear Mixed Model" = "linMM",
                                                     #"Ordinal Logit" = "ordLog",
                                                     "(-_-)_/¯" = "idk")),
                            
                            tags$hr(),
                            
                            radioButtons(inputId = "multiCompars", 
                                         label = "Select correction method to calculate significance threshold",
                                         choices = c("Gemone-wide significance (5e-8)" = "wgs",
                                                     "Bonferroni at alpha = 0.05 (alpha/No. SNPs)" = "bonf"),
                                         selected = "bonf"),
                            
                            tags$hr(),
                            
                            # FUNCTIONALITY HIDDEN FOR NOW
                            # # Additional inputs/parameter config for linear/logistic mixed models ----
                            # conditionalPanel(
                            #     condition = "input.assocTest == 'linMM' || input.assocTest == 'logMM'",
                            #     
                            #     # Input: Select random effects (e.g. kinship matrix) ----
                            #     fileInput(inputId = "randEffect", 
                            #               label = "Select random effects file",
                            #               multiple = FALSE,
                            #               accept = c(".csv")) #Not sure what to accept here
                            #     
                            # ), #END panel
                            
                            tags$p(
                                tags$b(paste("Select a folder where the results",
                                             "of association analysis will be saved."))),
                            
                            shinyDirButton(id = "saveDir", 
                                           label = "Select a folder",
                                           title = "Select a folder"),
                            
                            br(), br(),
                            
                            verbatimTextOutput("directorypath"),
                            
                            tags$hr(),
                            
                            #FUNCTIONALITY REMOVED FOR NOW
                            # # Input: Submit button to run null model ----
                            # actionButton(inputId = "nullDo",
                            #              label = "Run Null Model"),
                            
                            # Input: Submit button to run association analysis ----
                            actionButton(inputId = "assocDo",
                                         label = "Run Association")
                            
                            ), #END column
                     
                     column(8,
                            
                            tags$hr(),
                            
                            # Output: A busy status indicator
                            add_busy_spinner(spin = "fading-circle"),
                            
                            #FUNCTIONALITY REMOVED FOR NOW
                            # # Output: show fixed effects table from null model ----
                            # DT::dataTableOutput(outputId = "nullSummary"),
                            
                            # Output: show summary table for gwas results ----
                            DT::dataTableOutput(outputId = "gwasSig"),
                            
                            # Output: show summary table for gwas results ----
                            DT::dataTableOutput(outputId = "gwasSummary"),
                            
                     ) #END column
                 ) #END row
        ),#END tab
        
        ################################################
        ## Tab 3 Display GWAS results and figures     ##
        ################################################
        
        tabPanel("Results", fluid = TRUE,
                 
                 # Title ----
                 titlePanel("Results and Figures"),
                 
                 fluidRow(
                     tags$hr(),
                     column(7,
                            
                            # Output: A busy status indicator
                            add_busy_spinner(spin = "fading-circle"),
                            
                            #Output: Manhattan plot ----
                            plotOutput(outputId = "manPlot")
                            
                     ), #END column
                     
                     column(5,
                            
                            #Output: QQ plot ----
                            plotOutput(outputId = "qqPlot")
                            
                     ) #END column
                 ), #END row
                 
        ) #END tab
     ) #END tabset
) #END ui



server <- function(input, output, session) {
    
    source("./shinyGWAS_fns.R")
    
    ################################################
    ## Tab 1 Uploading Data                       ##
    ################################################
    
    # --- UPLOAD GENOME DATA
    
    observeEvent(input$isImpute, {
        
        if (input$isImpute == 1){
            choices <- c("VCF" = "vcf",
                         "PLINK" = "plink",
                         "IMPUTE2" = "IMPUTE2",
                         "BEAGLE" = "BEAGLE",
                         "MaCH" = "MaCH")
            selected <- "vcf"
            
        } else {
            choices <- c("GDS" = "gds",
                         "VCF" = "vcf",
                         "PLINK" = "plink")
            selected <- "gds"
        }
        
        updateRadioButtons(session, inputId = "fileType",
                           choices = choices,
                           selected = selected)
    })
    
    
    # --- SAVE GDS FILE PATH
    # --- converts files if required
    # --- NON-IMPUTED DATA ONLY
    
    gdsPath <- eventReactive(input$genoDo, {
        req(input$fileType)
        req(input$isImpute != 1)
        
        showfile.gds(closeall = TRUE)
        
        if (input$fileType == "gds"){
            req(input$gdsFile)
            path <- input$gdsFile$datapath %>% as.character
            
        } else if (input$fileType == "vcf"){
            req(input$vcfFile)
            path <- input$vcfFile$datapath %>% as.character
            path <- convertToGds(path, "vcf")
            
        } else if (input$fileType == "plink"){
            req(input$bedFile)
            req(input$bimFile)
            req(input$famFile)
            paths <- list("bed" = input$bedFile$datapath,
                          "bim" = input$bimFile$datapath,
                          "fam" = input$famFile$datapath)
            path <- convertToGds(paths, "plink")
        }
        
        return(path)
    })
    
    
    # --- SAVE GDS FILEPATHS TO LIST
    # --- Converts files if required
    # --- IMPUTED DATA ONLY
    
    imputePaths <- eventReactive(input$genoDo, {
        req(input$fileType)
        req(input$isImpute == 1)
        req(input$isDosage)
        
        showfile.gds(closeall = TRUE)
        
        paths <- list()
        
        if (input$fileType == "vcf"){
            req(input$vcfImps)
            
            for (ii in 1:nrow(input$vcfImps)){
                inPath <- input$vcfImps$datapath[ii] %>% as.character
                
                paths[ii] <- imputeToGds(inPath, 
                                         input$fileType, 
                                         input$isDosage)
            }
            
        } else if (input$fileType == "IMPUTE2"){
            req(input$imp2Imps1)
            req(input$imp2Imps2)
            
            for (ii in 1:nrow(input$imp2Imps1)){
                inPath <- list()
                inPath[1] <- input$imp2Imps1$datapath[ii]
                inPath[2] <- input$imp2Imps2$datapath[ii]
                
                paths[ii] <- imputeToGds(inPath, 
                                         input$fileType, 
                                         input$isDosage,
                                         chr= ii)
            }
            
        } else if (input$fileType == "BEAGLE"){
            req(input$beagImps1)
            req(input$beagImps2)
            
            for (ii in 1:nrow(input$beagImps1)){
                inPath <- list()
                inPath[1] <- input$beagImps1$datapath[ii]
                inPath[2] <- input$beagImps2$datapath[ii]
                
                paths[ii] <- imputeToGds(inPath, 
                                         input$fileType, 
                                         input$isDosage,
                                         chr= ii)
            }
            
        } else if (input$fileType == "MaCH"){
            req(input$machImps1)
            req(input$machImps2)
            req(input$machImps3)
            
            for (ii in 1:nrow(input$machImps1)){
                inPath <- list()
                inPath[1] <- input$machImps1$datapath[ii]
                inPath[2] <- input$machImps2$datapath[ii]
                inPath[3] <- input$machImps3$datapath[ii]
                
                paths[ii] <- imputeToGds(inPath, 
                                         input$fileType, 
                                         input$isDosage,
                                         chr= ii)
            }}
        
        return(paths)
    })
    
    
    # --- SAMPLE DATA UPLOADED IN SEPARATE CSV FILE:
    
    sData <- reactive({
        req(input$phenoFile)
        
        path <- input$phenoFile$datapath %>% as.character
        read.csv(path)
    })
    
    observeEvent(sData(), {
        updateSelectInput(session, inputId = "id", 
                          choices = names(sData()))
    })
    
    observeEvent(sData(), {
        updateSelectInput(session, inputId = "pheno", 
                          choices = names(sData()))
    })
    
    observeEvent(sData(), {
        updateSelectInput(session, inputId = "sex", 
                          choices = names(sData()))
    })
    
    observeEvent(sData(), {
        updateSelectInput(session, inputId = "covars", 
                          choices = c("N/A" = "", names(sData())), 
                          selected = "")
    })
    
    
    # --- UPLOAD ADDITIONAL COVARIATE DATA:
    
    cData <- reactive({
        req(input$covarFile)
        
        path <- input$covarFile$datapath %>% as.character
        read.csv(path)
    })
    
    observeEvent(cData(), {
        updateSelectInput(session, inputId = "id2", 
                          choices = names(cData()))
    })
    
    observeEvent(cData(), {
        updateSelectInput(session, inputId = "pheno2", 
                          choices = c("N/A" = "", names(cData())),
                          selected = "")
    })
    
    observeEvent(cData(), {
        updateSelectInput(session, inputId = "moreCovars", 
                          choices = c("N/A" = "", names(cData())), 
                          selected = "")
    })
    
    
    # --- COLLATE AND SAVE SAMPLE DATA TO SCANANNOTATION OBJ
    
    scanAnnot <- eventReactive(input$genoDo, {
        
        if (input$moreCovars != "" | input$pheno2 != ""){
            id2 <- input$id2
            
            if (input$moreCovars != "" & input$pheno2 != ""){
                cols <- c(id2, input$pheno2, input$moreCovars)
                
            } else if (input$moreCovars != ""){
                cols <- c(id2, input$moreCovars)
                
            } else if (input$pheno2 != ""){
                cols <- c(id2, input$pheno2)
            }
            
            covars2 <- cData()[,cols]
            covars2[,id2] <- covars2[,id2] %>% as.factor
            
            if (input$pheno2 %in% colnames(covars2)){
                names(covars2)[names(covars2) == input$pheno2] <- 'XX_pheno_XX'
            }
        } else {
            covars2 <- NULL
        }

        if (input$fileType == "vcf"){
            req(sData())
            req(input$id)
            req(input$pheno)
            req(input$sex)
            
            id <- sData()[,input$id] %>% as.factor
            pheno <- sData()[,input$pheno]
            sex <- sData()[,input$sex]
            covars1 <- sData()[,input$covars]
            
            if (!is.null(covars2)){
                df <- buildScan(id, pheno, sex, covars1, id2, covars2)
            } else {
                df <- buildScan(id, pheno, sex, covars1)}
            
        } else {
            req(gdsPath())
            
            if (!is.null(covars2)){
                df <- extractScan(gdsPath(), id2, covars2)
            } else {
                df <- extractScan(gdsPath())}}
        
        if ("XX_pheno_XX" %in% colnames(df)){
            df$phenotype <- df$XX_pheno_XX
            df$XX_pheno_XX <- NULL
        }
        
        scan <- scanObj(df)
        scan %>% return
    })
    
    
    # --- MERGE GENOME AND SAMPLE DATA AS GENOTYPEDATA OBJ
    
    gData <- eventReactive(input$genoDo, {
        req(scanAnnot())
        req(gdsPath())
        
        openGDS(gdsPath(), scanAnnot())
    })
    
    
    # --- PRINT SUMMARY STATISTICS
    
    output$phenoSummary <- renderUI({
        req(scanAnnot())

        scanAnnot() %>% pData %>% dfSummary(display.labels = FALSE,
                                            graph.magnif = 0.65,
                                            headings = FALSE) %>% 
            print(method = 'render')
    })
    
    
    # --- PRINT PAIR WISE COMPARISONS
    # --- Code for dynamic re-sizing from 
    # --- https://stackoverflow.com/a/40539526
    
    output$heatmap <- renderPlot({
        req(scanAnnot())
        
        (scanAnnot() %>% pData)[, -1] %>% 
            ggpairs(upper = list(continuous = "cor",
                                 combo = "facethist",
                                 discrete = "facetbar"),
                    lower = list(continuous = "points",
                                 combo = "facethist",
                                 discrete = "blank"))
    }, height = reactive(ifelse(!is.null(input$innerWidth),
                                input$innerWidth*3/5,0)))
    
    # --- PRINT LIST OF UPLOADED FILES
    
    output$filesUploaded <- DT::renderDataTable({
        req(input$genoDo)
        req(input$fileType)
        
        if (input$fileType == "gds"){
            req(input$gdsFile)
            df <- input$gdsFile
            
        } else if (input$fileType == "vcf"){
            req(input$vcfFile)
            df <- input$vcfFile
            
        } else if (input$fileType == "plink"){
            req(input$bedFile)
            req(input$bimFile)
            req(input$famFile)
            df <- rbind(input$bedFile,
                        input$bimFile,
                        input$famFile)
            
        } else if (input$fileType == "IMPUTE2"){
            req(input$imp2Imps1)
            req(input$imp2Imps2)
            df <- rbind(input$imp2Imps1,
                        input$imp2Imps2)
            
        } else if (input$fileType == "BEAGLE"){
            req(input$beagImps1)
            req(input$beagImps2)
            df <- rbind(input$beagImps1,
                        input$beagImps2)
            
        } else if (input$fileType == "MaCH"){
            req(input$machImps1)
            req(input$machImps2)
            req(input$machImps3)
            df <- rbind(input$machImps1,
                        input$machImps2,
                        input$machImps3)
        }
        
        if (!is.null(input$phenoFile)){
            df <- rbind(df,
                        input$phenoFile)
        }
        
        if (!is.null(input$covarFile)){
            df <- rbind(df,
                        input$covarFile)
        }
        
        
        df <- df[,c(1,2)]
        colnames(df) <- c("Name", "Size_(MB)")
        df[,2] <- df[,2]/1000000
        df
    })
 
    
    ################################################
    ## Tab 2 Configure Association Test Settings  ##
    ################################################
    
    # --- UPDATE ASSOCIATION TEST OPTIONS
    
    observeEvent(input$phenoDatType, {
        
        if (input$phenoDatType == "bin"){
            choices <- c("Logistic Regression" = "logistic")
            selected <- "logistic"
            
        } else if (input$phenoDatType == "quant"){
            choices <- c("Linear Regression" = "linear")
            selected <- "linear"
            
        } else if (input$phenoDatType == "cat"){
            choices <- c("(-_-)_/¯" = "idk")
            selected <- "idk"
        }
        
        updateRadioButtons(session, inputId = "assocTest",
                           choices = choices,
                           selected = selected)
    })
    
    
    # --- USER SELECT A DIRECTORY TO SAVE ASSOCIATION TEST RESULTS
    # --- CODE REPLICATED FROM shinyFiles PACKAGE EXAMPLE APP
    # --- RUN APPLICATION BY PASSING shinyFilesExample() IN CONSOLE.
    # --- https://cran.r-project.org/web/packages/shinyFiles/shinyFiles.pdf
    
    volumes <- c(Home = fs::path_home(), 
                 "R Installation" = R.home(), 
                 getVolumes()())
    
    shinyDirChoose(input = input, 
                   id = "saveDir",
                   roots = volumes,
                   session = session, 
                   restrictions = system.file(package = "base"))
    
    output$directorypath <- renderPrint({
        if (is.integer(input$saveDir)) {
            cat("No directory selected.")
        } else {
            parseDirPath(roots = volumes, 
                         selection = input$saveDir)
        }
    })
    
    
    # # --- RUN NULL MODEL ON FIXED VARIABLES
    #
    # assocNull <- eventReactive(input$nullDo, {
    #     req(gData())
    #     gwasNull(gData(), input$phenoDatType)
    # })
    
    
    # --- RUN ASSOCIATION ANALYSIS AND SAVE RESULTS
    
    assocData <- eventReactive(input$assocDo, {
        req(gData())
        req(input$saveDir)
        
        assoc <- gwasLogLin(gData(), input$assocTest)

        # if (input$assocTest == "logistic" | input$assocTest == "linear"){
        #     assoc <- gwasLogLin(gData(), input$assocTest)
        # 
        # } else if (input$assocTest == "linMM" | input$assocTest == "logMM"){
        #     req(assocNull())
        #     gwasMixed(gData(), assocNull())
        # }
        
        assoc <- addPosition(gData(), assoc)
        
        sig <- Sys.time() %>% format("%Y%m%d-%H%M%S")
        savePath <- parseDirPath(roots = volumes, selection = input$saveDir) %>% 
            as.character %>% paste("/shinyGWAS.", sig, ".csv", sep = "")
        
        assoc %>% write.csv(savePath)
        assoc %>% return
    })
    
    
    # --- UPDATE SIGNIFICANCE VALUE
    
    sigValue <- reactive({
        
        if (input$multiCompars == "wgs"){
            thresh <- 5e-8
        } else if (input$multiCompars == "bonf"){
            thresh <- 0.05/(nrow(assocData()))
        }
        
        thresh %>% return
    })
    
    
    # # --- PRINT FIXED EFFECTS RESULTS CHI-SQUARE, BETA, P-VAL
    #
    # output$nullSummary <- DT::renderDataTable({
    #     req(assocNull())
    #     assocNull()$fixef
    #     })
    
    # --- PRINT SIGNIFICANT SNPS
    
    output$gwasSIG <- DT::renderDataTable({
        req(assocData())
        
        df <- assocData()[which(assocData()$Wald.pval < sigValue()),]
        df <- df[order(df$Wald.pval), ]
        
        datatable(assocData(), 
                  rownames = TRUE, 
                  options = list(scrollX = TRUE))
    })
    
    
    # --- PRINT RESULTS OF ASSOCIATION TEST
    
    output$gwasSummary <- DT::renderDataTable({
        req(assocData())
        
        datatable(assocData(), 
                  rownames = TRUE, 
                  options = list(scrollX = TRUE))
    })

    
    ################################################
    ## Tab 4 Display GWAS results and figures     ##
    ################################################
    
    output$manPlot <- renderCachedPlot({
        
        p <- assocData()$Wald.pval
        chr <- assocData()$chr
        sig <- sigValue()
        ylim <- range(0, abs(log10(length(p)) + 4), abs(log10(sig)))
        
        manhattanPlot(p = p,
                      chromosome = chr,
                      signif = sig, 
                      ylim = ylim)
        
    }, cacheKeyExpr = list(sigValue(), assocData()))
    
    
    output$qqPlot <- renderCachedPlot({
        
        par(pty = "s")
        qq(assocData()$Wald.pval)
    }, cacheKeyExpr = list(sigValue(), assocData()))


} #END server

# Run the application 
shinyApp(ui = ui, server = server)


