
library(shiny)
library(shinybusy)
library(magrittr)
# library(GENESIS)
library(SNPRelate)
# library(lme4)
library(GWASTools)
library(GWASdata)
library(qqman)
library(summarytools)
library(GGally)
library(shinyFiles)
library(DT)
library(plotly)
library(qqman)
library(dplyr)


# --- Set size of allowed file uploads.
gb <- 3 #number of GB
options(shiny.maxRequestSize = (gb*1000)*1024^2)

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
                            
                            # --- Side bar with two tabs for uploading data
                            tabsetPanel(
                                # --- Sidebar tab #1 where user uploads genome data
                                tabPanel(title = "Genome", fluid = TRUE,
                                         
                                         br(),
                                         
                                         # --- Input: User checks box if using imputed data.
                                         # --- Default selection is FALSE (not imputed)
                                         checkboxInput(inputId = "isImpute", 
                                                       label = "Check box if data has been imputed", 
                                                       value = FALSE),
                                         
                                         tags$hr(),
                                         
                                         # --- Input: User selects the file format of genome data.
                                         radioButtons(inputId = "fileType",
                                                      label = "Select file format",
                                                      choices = c("GDS" = "gds",
                                                                  "VCF" = "vcf",
                                                                  "PLINK" = "plink",
                                                                  "IMPUTE2" = "IMPUTE2",
                                                                  "BEAGLE" = "BEAGLE",
                                                                  "MaCH" = "MaCH")),
                                         
                                         tags$hr(),
                                         
                                         # --- If data is imputed, user is asked to specify whether imputed data
                                         # --- is recorded as dosages or genotype probabilities.
                                         # --- Default selection is dosages.
                                         conditionalPanel(
                                             condition = "input.isImpute == 1",
                                             
                                             radioButtons(inputId = "isDosage",
                                                          label = "Indicate how imputed data is represented",
                                                          choices = c("Dosage" = TRUE,
                                                                      "Genotype Probabilities" = FALSE),
                                                          selected = TRUE),
                                             
                                             tags$hr()
                                             
                                         ), #END panel
                                         
                                         # --- If genome data is GDS format, user is asked to upload GDS file.
                                         conditionalPanel(
                                             condition = "input.fileType == 'gds'",
                                             
                                             fileInput(inputId = "gdsFile", 
                                                       label = "Select GDS file",
                                                       multiple = FALSE,
                                                       accept = c(".gds"))
                                         ), #END panel
                                         
                                         # --- If genome data is VCF format, user is asked to upload VCF file
                                         # --- User can upload multiple VCF files if required.
                                         conditionalPanel(
                                             condition = "input.fileType == 'vcf'",
                                             
                                             fileInput(inputId = "vcfFile", 
                                                       label = "Select VCF file(s)",
                                                       multiple = TRUE,
                                                       accept = c(".vcf"))
                                         
                                          ), #END panel
                                         
                                         # --- If genome data is PLINK format, user is asked to upload
                                         # --- BED, BIM & FAM files.
                                         conditionalPanel(
                                             condition = "input.fileType == 'plink'",
                                             
                                             fileInput(inputId = "bedFile", 
                                                       label = "Select PLINK .BED file",
                                                       multiple = FALSE,
                                                       accept = c(".bed")),
                                             
                                             fileInput(inputId = "bimFile", 
                                                       label = "Select PLINK .BIM file",
                                                       multiple = FALSE,
                                                       accept = c(".bim")),
                                             
                                             fileInput(inputId = "famFile", 
                                                       label = "Select PLINK .FAM file",
                                                       multiple = FALSE,
                                                       accept = c(".fam"))
                                             
                                         ), #END panel
                                         
                                         # --- If genome data is IMPUTE2 format, user is asked to
                                         # --- upload GENS & SAMPLES files.
                                         # --- Multiple files can be selected.
                                         conditionalPanel(
                                             condition = "input.fileType == 'IMPUTE2'",
                                             
                                             fileInput(inputId = "imp2Imps1", 
                                                       label = "Select IMPUTE2 .gens file(s)",
                                                       multiple = TRUE,
                                                       accept = c(".gens")),
                                             
                                             fileInput(inputId = "imp2Imps2", 
                                                       label = "Select IMPUTE2 .samples file(s)",
                                                       multiple = TRUE,
                                                       accept = c(".samples"))
                                             
                                         ), #END panel
                                         
                                         # --- If genome data is BEAGLE format, user is asked to 
                                         # --- upload GROBS & MARKERS files.
                                         # --- Multiple files can be selected.
                                         conditionalPanel(
                                             condition = "input.fileType == 'BEAGLE'",
                                             
                                             fileInput(inputId = "beagImps1", 
                                                       label = "Select BEAGLE .grobs or .dose file(s)",
                                                       multiple = TRUE,
                                                       accept = c(".grobs",
                                                                  ".dose")),
                                             
                                             fileInput(inputId = "beagImps2", 
                                                       label = "Select BEAGLE .markers file(s)",
                                                       multiple = TRUE,
                                                       accept = c(".markers"))
                                             
                                         ), #END panel
                                         
                                         # --- If genome data is MaCH format, user is asked to 
                                         # --- upload MLPROB/MLDOSE, MLINFO & TXT/CSV files.
                                         # --- Multiple files can be selected.
                                         conditionalPanel(
                                             condition = "input.fileType == 'MaCH'",
                                             
                                             fileInput(inputId = "machImps1", 
                                                       label = "Select MaCH .mlprob/.mldose file(s)",
                                                       multiple = TRUE,
                                                       accept = c(".mlprob",
                                                                  ".mldose")),
                                             
                                             fileInput(inputId = "machImps2", 
                                                       label = "Select MaCH .mlinfo file(s)",
                                                       multiple = TRUE,
                                                       accept = c(".mlinfo")),
                                             
                                             fileInput(inputId = "machImps3", 
                                                       label = "Select file(s) with columns 'SNP' and 'position'",
                                                       multiple = TRUE,
                                                       accept = c(".csv", ".txt"))
                                             
                                         ), #END panel
                                ), #END tab
                                
                                # --- Sidebar tab #2 where user uploads additional 
                                # --- sample data if required and covariate data.
                                tabPanel(title = "Sample", fluid = TRUE,
                                         br(),
                                         
                                         # --- If genome data is VCF format, show instructions for
                                         # --- compulsary upload of additional sample information.
                                         conditionalPanel(
                                             condition = "input.fileType == 'vcf'",
                                             
                                             tags$p(paste("VCF files must be accompanied by",
                                                          "a text file containing at minimum, the fields:",
                                                          "sample ID, phenotype and sex.",
                                                          "An ID variable corresponding to the",
                                                          "genotype data must be selected from the ID field below.")),
                                             br(),
                                         ), #END panel
                                         
                                         # --- If genome data is NOT vcf format, show instructions for
                                         # --- optional upload of additional sample data.
                                         conditionalPanel(
                                             condition = "input.fileType != 'vcf'",
                                             
                                             tags$p(paste("Upload any additional sample data",
                                                          "for inclusion in association model.",
                                                          "An ID variable corresponding to the",
                                                          "genotype data must be selected from the ID field below.",
                                                          "If phenotype or sex variables are included",
                                                          ", these must be selected from the relevant",
                                                          "fields below.")),
                                             
                                             br(),
                                         ),
                        
                                         # --- User uploads additional sample data file in CSV or TXT format
                                         fileInput(inputId = "sampFile", 
                                                   label = "Select sample data file",
                                                   multiple = TRUE,
                                                   accept = c(".csv",
                                                              ".txt"),
                                                   placeholder = ".csv or .txt"),
                                         
                                         # --- User indicates which column contains the unique ID.
                                         # --- This field is compulsary if any sample data is uploaded.
                                         selectInput(inputId = "id",
                                                     label = "ID",
                                                     choices = character(0),
                                                     selectize = FALSE),

                                         # --- User indicates which column is the phenotype of interest.
                                         # --- This field is compulsary if genome data is VCF format,
                                         # --- otherwise this is optional.
                                         selectInput(inputId = "pheno",
                                                     label = "Phenotype",
                                                     choices = character(0),
                                                     selectize = FALSE),
                                         
                                         # --- User indicates which column contains sex information.
                                         # --- This field is compulsary if genome data is VCF format,
                                         # --- otherwise this is optional.
                                         selectInput(inputId = "sex",
                                                     label = "Sex",
                                                     choices = character(0),
                                                     selectize = FALSE),

                                         # --- User selects any additional covariates to be included in
                                         # --- the analyses. This field is always optional.
                                         # --- Multiple selections are allowed.
                                         selectInput(inputId = "covars",
                                                     label = "Other covariates",
                                                     choices = character(0),
                                                     selectize = FALSE,
                                                     multiple = TRUE),
                                         
                                         # --- User selects this button to upload data and see summary
                                         # --- information.
                                         actionButton(inputId = "genoDo",
                                                      label = "Upload genome data"),
                                         
                                 ) #END tab
                            ), #END tabset
                            
                            tags$hr(),
                            
                            # --- Places a spinning circle in the corner of the page to
                            # --- show when the app is busy computing.
                            add_busy_spinner(spin = "fading-circle")
                            
                            ), #END column
                     
                     # --- Right-hand panel which displays summary figures and information
                     # --- about the uploaded genome and sample data.
                     column(8,
                            
                            # --- The following code makes the pair-wise plot re-size to
                            # --- the size of window.
                            # --- Code from https://stackoverflow.com/a/40539526
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
                                # --- First tab shows descriptive information about the 
                                # --- data in a table format.
                                # --- htmlOutput() shows table as it would appear when 
                                # --- knitted in R, as opposed to the shiny table format.
                                tabPanel(title = "Descriptive Statistics", fluid = TRUE,
                                         br(), br(),
                                         
                                         htmlOutput(outputId = "phenoSummary")
                                         
                                         ), #END tab
                                
                                # --- Second tab shows a correlation plot with histograms.    
                                tabPanel(title = "Pair-wise Comparisons", fluid = TRUE,
                                         br(), br(),
                                         
                                         plotOutput(outputId = "heatmap")
                                         
                                         ), #END tab
                                
                                # --- Third tab shows the files that have been uploaded and their size.
                                # --- DT::dataTableOutput returns an interactive table.
                                tabPanel(title = "Files", fluid = TRUE,
                                         br(), br(),
                                         
                                         DT::dataTableOutput(outputId = "filesUploaded")
                                         
                                         ) #END tab
                            ) #END tabset
                            ) #END column
                     ) #END row
        ), #END tab
        
        ################################################
        ## Tab 2 Configure Association Test           ##
        ################################################
        
        tabPanel("Association", fluid = TRUE,
                 titlePanel("Association Test Settings"),
                 
                 fluidRow(
                     # --- Sidebar where association test options
                     # --- are displayed.
                     column(4,
                            tags$hr(),
                            
                            # --- User selects the data-type of the phenotype data.
                            # --- Currently only supports binary and quantitative data-types.
                            # --- Default selection is Binary.
                            radioButtons(inputId = "phenoDatType",
                                         label = "Select phenotype data-type",
                                         choices = c("Binary" = "bin",
                                                     "Quantitative" = "quant",
                                                     "Categorical" = "cat"),
                                         selected = "bin"),
                            
                            tags$hr(),
                            
                            # --- User selects the assoication model they wish to use.
                            # --- The available models will depend on the data-type selection
                            # --- made above. 
                            # --- Currently only linear and logistic regression models available.
                            radioButtons(inputId = "assocTest",
                                         label = "Select association model",
                                         choices = c("Logistic Regression" = "logistic", 
                                                     "Linear Regression" = "linear",
                                                     "(-_-)_/¯" = "idk")),
                            
                            tags$hr(),
                            
                            # --- User selects the type of multiple comparisons adjustment they want
                            # --- to calculate the significance threshold for association. 
                            # --- Currently only genome-wide significance and Bonferroni available.
                            # --- Default selection is Bonferroni.
                            radioButtons(inputId = "multiCompars", 
                                         label = "Select correction method to calculate significance threshold",
                                         choices = c("Gemone-wide significance (5e-8)" = "wgs",
                                                     "Bonferroni at alpha = 0.05 (alpha/No. SNPs)" = "bonf"),
                                         selected = "bonf"),
                            
                            tags$hr(),
                            
                            # --- User is asked to select a folder where the results of 
                            # --- association testing and the summary report will be
                            # --- automatically saved when association is complete.
                            tags$p(
                                tags$b(paste("Select a folder where the results",
                                             "of association analysis will be saved."))),
                            
                            shinyDirButton(id = "saveDir", 
                                           label = "Select a folder",
                                           title = "Select a folder"),
                            
                            br(), br(),
                            
                            # --- Shows the folder path that has been selected.
                            verbatimTextOutput("directorypath"),
                            
                            tags$hr(),
                            
                            # --- User clicks this button to confirm settings and 
                            # --- begin the association test.
                            actionButton(inputId = "assocDo",
                                         label = "Run Association")
                            
                            ), #END column
                     
                     # --- Right-hand panel where results of association test will be displayed.
                     column(8,
                            tags$hr(),
                            
                            # --- Places a spinning circle in the corner of the page to
                            # --- show when the app is busy computing.
                            add_busy_spinner(spin = "fading-circle"),
                            
                            # --- Shows the statistic results of the association test
                            # --- in an interactive table.
                            DT::dataTableOutput(outputId = "gwasSummary"),
                            
                     ) #END column
                 ) #END row
        ),#END tab
        
        ################################################
        ## Tab 3 Display GWAS results and figures     ##
        ################################################
        
        tabPanel("Results", fluid = TRUE,
                 titlePanel("Results and Figures"),
                 
                 fluidRow(
                     tags$hr(),
                     column(7,
                            
                            # --- Places a spinning circle in the corner of the page to
                            # --- show when the app is busy computing.
                            add_busy_spinner(spin = "fading-circle"),
                            
                            # --- Shows a manhattan plot of the p-values in order of
                            # --- chromosome and SNP position.
                            plotOutput(outputId = "manPlot")
                            
                     ), #END column
                     
                     column(5,
                            
                            # --- Shows a QQ plot
                            plotOutput(outputId = "qqPlot")
                            
                     ) #END column
                 ), #END row
        ), #END tab
        
        ################################################
        ## Tab 4 Interactive Manhattan Plot           ##
        ################################################
        
        tabPanel("Interactive Manhattan", fluid = TRUE,
                 titlePanel("Interactive Manhattan Plot"),
                 
                 fluidRow(
                     column(4,
                            tags$hr(),
                            
                            # --- Places a spinning circle in the corner of the page to
                            # --- show when the app is busy computing.
                            add_busy_spinner(spin = "fading-circle"),
                            
                            # --- User selects the chromosomes they wish to show on plot.
                            # --- Multiple chromosomes can be selected.
                            # --- Chromosome selection allows plot to render faster, especially
                            # --- if user is only interested in a few chromosomes.
                            selectInput(inputId = "chr",
                                        label = "Select Chromosome(s)",
                                        multiple = TRUE,
                                        choices = character(0),
                                        selectize = FALSE),
                            
                     ), #END column
                 ), #END row
                 
                 fluidRow(
                     tags$hr(),
                     column(11,
                            
                            # --- Show interactive Manhattan plot.
                            plotlyOutput(outputId = "manPlot2")
                            
                     ) #END column
                 ) #END row
            ) #END tab
     ) #END tabset
) #END ui



server <- function(input, output, session) {
    
    # --- This file must be in working directory.
    source("./shinyGWAS_fns.R")
    
    ################################################
    ## Tab 1 Uploading Data                       ##
    ################################################
    
    # --- UPLOAD GENOME DATA:
    # --- Options available will depend on whether
    # --- or not data has been imputed.
    
    observeEvent(input$isImpute, {
        
        if (input$isImpute == 1){ #If data is imputed...
            choices <- c("VCF" = "vcf",
                         "PLINK" = "plink",
                         "IMPUTE2" = "IMPUTE2",
                         "BEAGLE" = "BEAGLE",
                         "MaCH" = "MaCH")
            selected <- "vcf" #default option
        
        } else { #If data NOT imputed...
            choices <- c("GDS" = "gds",
                         "VCF" = "vcf",
                         "PLINK" = "plink")
            selected <- "gds" #default option
        }
        
        #Update the options
        updateRadioButtons(session, inputId = "fileType",
                           choices = choices,
                           selected = selected)
    })
    
    
    # --- SAVE GDS FILE PATH AS STRING.
    # --- Converts files to GDS if required.
    
    gdsPath <- eventReactive(input$genoDo, {
        req(input$fileType)
        
        showfile.gds(closeall = TRUE) #close any open gds files to avoid errors.
        
        outFile <- tempfile(pattern = "shinyGWAS", #temp file to save converted GDS file 
                            fileext = ".gds")
        
        if (input$fileType == "gds"){
            req(input$gdsFile)
            path <- input$gdsFile$datapath %>% as.character
            
        } else if (input$fileType == "vcf"){
            req(input$vcfFile)
            inPath <- input$vcfFile$datapath %>% as.character
            
            if (input$isImpute != 1){ #if data not imputed...
                path <- snpgdsVCF2GDS(vcf.fn = inPath, 
                                      out.fn = outFile)
                
            } else { #if data IS imputed...
                ##################################
                ########      WORKING       ######
                ##################################
                validate(
                    need(input$isDosage == "dosage", 
                         message = "Imputed data must be represented as dosages.")
                )
                
                seqVCF2GDS(vcf.fn = inPath, 
                           out.fn = outFile, 
                           fmt.import="DS")
                ##################################
                ##################################
                ##################################
            }
            
        } else if (input$fileType == "plink"){
            req(input$bedFile)
            req(input$bimFile)
            req(input$famFile)
            
            bed <- input$bedFile$datapath %>% as.character
            bim <- input$bimFile$datapath %>% as.character
            fam <- input$famFile$datapath %>% as.character
            
            path <- snpgdsBED2GDS(bed.fn = bed,
                                  bim.fn = bim,
                                  fam.fn = fam,
                                  out.gdsfn = outFile, 
                                  family = TRUE, 
                                  cvt.chr = "int")
        }
        return(path)
    })
    
    
    # --- SAVE GDS FILEPATHS TO LIST
    # --- Converts files if required
    # --- IMPUTED DATA ONLY
    
    # imputePaths <- eventReactive(input$genoDo, {
    #     req(input$fileType)
    #     req(input$isImpute == 1)
    #     req(input$isDosage)
    #     
    #     showfile.gds(closeall = TRUE)
    #     
    #     paths <- list()
    #     
    #     if (input$fileType == "vcf"){
    #         req(input$vcfImps)
    #         
    #         validate(
    #             need(input$isDosage == "dosage", 
    #                  message = "Imputed data must be represented as dosages.")
    #         )
    #         
    #         inPath <- input$vcfImps$datapath
    #         paths <- inputToGds(inPath,
    #                             input$fileType,
    #                             input$isDosage)
            # 
            # CODE FOR LOOPING THROUGH FILES AND MAKING MULTIPLE GDS
            # for (ii in 1:nrow(input$vcfImps)){
            #     inPath <- input$vcfImps$datapath[ii] %>% as.character
            #     
            #     paths[ii] <- imputeToGds(inPath, 
            #                              input$fileType, 
            #                              input$isDosage)
            # }
        # }
            # 
        # } else if (input$fileType == "IMPUTE2"){
        #     req(input$imp2Imps1)
        #     req(input$imp2Imps2)
        #     
        #     for (ii in 1:nrow(input$imp2Imps1)){
        #         inPath <- list()
        #         inPath[1] <- input$imp2Imps1$datapath[ii]
        #         inPath[2] <- input$imp2Imps2$datapath[ii]
        #         
        #         paths[ii] <- imputeToGds(inPath, 
        #                                  input$fileType, 
        #                                  input$isDosage,
        #                                  chr= ii)
        #     }
        #     
        # } else if (input$fileType == "BEAGLE"){
        #     req(input$beagImps1)
        #     req(input$beagImps2)
        #     
        #     for (ii in 1:nrow(input$beagImps1)){
        #         inPath <- list()
        #         inPath[1] <- input$beagImps1$datapath[ii]
        #         inPath[2] <- input$beagImps2$datapath[ii]
        #         
        #         paths[ii] <- imputeToGds(inPath, 
        #                                  input$fileType, 
        #                                  input$isDosage,
        #                                  chr= ii)
        #     }
        #     
        # } else if (input$fileType == "MaCH"){
        #     req(input$machImps1)
        #     req(input$machImps2)
        #     req(input$machImps3)
        #     
        #     for (ii in 1:nrow(input$machImps1)){
        #         inPath <- list()
        #         inPath[1] <- input$machImps1$datapath[ii]
        #         inPath[2] <- input$machImps2$datapath[ii]
        #         inPath[3] <- input$machImps3$datapath[ii]
        #         
        #         paths[ii] <- imputeToGds(inPath, 
        #                                  input$fileType, 
        #                                  input$isDosage,
        #                                  chr= ii)
        #     }}
    #     
    #     return(paths)
    # })
    
    
    # --- READ IN ADDITIONAL SAMPLE DATA
    
    sData <- reactive({
        req(input$sampFile)
        
        path <- input$sampFile$datapath %>% as.character
        read.csv(path)
    })
    
    # --- Update selection fields with column names from sample data file
    # --- for selecting the columns that correspond to ID, 
    # --- phenotype, sex and additional covariates.
    
    observeEvent(sData(), {
        updateSelectInput(session, inputId = "id", 
                          choices = names(sData()))
    })
    
    observeEvent(sData(), {
        if (input$fileType != "vcf"){ 
            updateSelectInput(session, inputId = "pheno", 
                              choices = c("N/A" = "", names(sData())), 
                              selected = "")
        } else { #if genome data is vcf format, this is compulsary
            updateSelectInput(session, inputId = "pheno", 
                              choices = names(sData()))
        }
    })
    
    observeEvent(sData(), {
        if (input$fileType != "vcf"){
            updateSelectInput(session, inputId = "sex", 
                              choices = c("N/A" = "", names(sData())), 
                              selected = "")
        } else { #if genome data is vcf format, this is compulsary
            updateSelectInput(session, inputId = "sex", 
                              choices = names(sData()))
        }
        
    })
    
    observeEvent(sData(), {
        updateSelectInput(session, inputId = "covars", 
                          choices = c("N/A" = "", names(sData())), 
                          selected = "")
    })

    
    # --- COLLATE AND SAVE SAMPLE DATA TO SCANANNOTATION OBJ
    # --- Depending on genome data format, will extract 
    # --- sample data from genome data file, if additional
    # --- sample data has been provided this will be used instead.
    
    scanAnnot <- eventReactive(input$genoDo, { #run this when file upload button is selected
        
        if (input$fileType == "vcf"){ #if genome data in VCF format...
            req(sData())
            req(input$id)
            req(input$pheno)
            req(input$sex)
            
            id <- sData()[,input$id] %>% as.factor
            pheno <- sData()[,input$pheno]
            sex <- sData()[,input$sex]
            covars <- sData()[,input$covars]
            
            df <- data.frame(scanID = id,
                             phenotype = pheno,
                             sex = sex)
            
            if (covars != ""){
                df <- data.frame(df,covars)
            }
            
        } else { #if genome data NOT in VCF format...
            req(gdsPath())
            
            if (!is.null(input$sampFile)){ #if additional sample data is provided...
                req(sData())
                req(input$id)
                
                validate( #at least one variable must be selected in addition to ID
                    need(input$sex != "" | input$pheno != "" | input$covars != "", 
                         message = "No variables selected.")
                )
                
                
                id <- input$id
                
                if (input$covars != ""){ 
                    if (input$pheno != ""){ 
                        if (input$sex != ""){ #covariates, phenotype and sex...
                            cols <- c(id, input$pheno, input$sex, input$covars)
                            
                        } else if (input$sex == ""){ #covariates and phenotype...
                            cols <- c(id, input$pheno, input$covars)}
                        
                    } else if (input$pheno == ""){
                        if (input$sex != ""){ #covariates and sex...
                            cols <- c(id, input$sex, input$covars)
                            
                        } else if (input$sex == ""){ #covariates only...
                            cols <- c(id, input$covars)}}
                    
                } else if (input$covars == ""){
                    if (input$pheno != ""){
                        if (input$sex != ""){ #phenotype and sex...
                            cols <- c(id, input$pheno, input$sex)
                            
                        } else if (input$sex == ""){ #phenotype only...
                            cols <- c(id, input$pheno)}
                        
                    } else if (input$pheno == ""){
                        if (input$sex != ""){ #sex only...
                            cols <- c(id, input$sex)
                            
                        } else if (input$sex == ""){ #no selections...
                            cols <- c(id)}} ### RAISE AN ERROR HERE
                }
                
                covars <- sData()[,cols] #extract the variables to add to scan object
                covars[,id] <- covars[,id] %>% as.factor
                
                # --- Rename columns if phenotype or sex have been selected
                # --- in a separate file in case they have the same names as the
                # --- columns extracted from the genome data.
                # --- Can then replace the columns simply.
                
                if (input$pheno %in% colnames(covars)){
                    names(covars)[names(covars) == input$pheno] <- 'XX_pheno_XX'
                }
                
                if (input$sex %in% colnames(covars)){
                    names(covars)[names(covars) == input$sex] <- 'XX_sex_XX'
                }
                
                df <- extractScan(gdsPath(), id, covars) #call function from shinyGWAS_fns.R
                
                if ("XX_pheno_XX" %in% colnames(df)){ #replace extracted phenotype with selected
                    df$phenotype <- df$XX_pheno_XX
                    df$XX_pheno_XX <- NULL
                }
                
                if ("XX_sex_XX" %in% colnames(df)){ #replace extracted sex with selected
                    df$sex <- df$XX_sex_XX
                    df$XX_sex_XX <- NULL
                }
                
            } else { #if no additional sample file has been selected...
                df <- extractScan(gdsPath())
            }
        }

        scan <- df %>% ScanAnnotationDataFrame
        scan$sex <- scan$sex %>% as.factor
        scan %>% return
    })
    
    
    # --- MERGE GENOME AND SAMPLE DATA AS GENOTYPEDATA OBJECT
    
    gData <- eventReactive(input$genoDo, {
        req(scanAnnot())
        req(gdsPath())
        
        GdsGenotypeReader(filename = gdsPath()) %>% 
            GenotypeData(scanAnnot = scanAnnot())
    })
    
    
    # --- PRINT SUMMARY STATISTICS
    
    output$phenoSummary <- renderUI({
        req(scanAnnot())

        scanAnnot() %>% pData %>% dfSummary(display.labels = FALSE,
                                            graph.magnif = 0.65,
                                            headings = FALSE) %>% 
            print(method = 'render') #will display as if knitted to HMTL in R
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
        
        if (!is.null(input$sampFile)){ #if an additional sample file uploaded, append to list
            df <- rbind(df,
                        input$sampFile)
        }
        
        
        df <- df[,c(1,2)] #display only file name and size
        colnames(df) <- c("Name", "Size_(MB)") #update column names
        df[,2] <- df[,2]/1000000 #display file size in MB, not bytes.
        df
    })
 
    
    ################################################
    ## Tab 2 Configure Association Test Settings  ##
    ################################################
    
    # --- UPDATE ASSOCIATION TEST OPTIONS
    # --- Currently only supports logistic model for binary
    # --- and linear model for quantitative data.
    
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
    
    volumes <- c(Home = fs::path_home(), #set initial path to Home
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
    
    
    # --- RUN ASSOCIATION ANALYSIS AND SAVE RESULTS
    
    assocData <- eventReactive(input$assocDo, { #run association when user clicks button
        req(gData())
        req(input$saveDir)
        
        covars <- (gData() %>% getScanVariableNames)[-c(1,2)]
        
        if (length(covars) == 0){
            covars <- NULL
        }
        
        chr <- gData() %>% getChromosome
        
        # --- If the y chromosome is in the genome data
        # --- loop over each chromosome one at a time and perform
        # --- association test, because cannot use "sex" as a covariate
        # --- for the Y chromosome.
        
        if (25 %in% chr){ 
            assoc.list <- unique(chr) %>% lapply(function(x) {
                
                #if sex is the only covariate, make covar = NULL for CHR25
                if (length(covars) == 1){
                    covar <- ifelse(x == 25, yes = c(""), no = covars)
                    
                    if (covar == ""){
                        covar <- NULL
                    }
                    
                #if more covariates just remove the sex variable for CHR25 
                } else if (length(covars) > 1){
                    y <- covars[-which(covars == "sex")]
                    covar <- ifelse(x == 25, yes = y, no = covars)
                }
                
                start <- which(chr == x)[1] #index where current chromosome starts
                end <- start + (which(chr == x) %>% length) -1 #index where it ends
                
                #run association
                assocRegression(genoData = gData(),
                                outcome = "phenotype",
                                model.type = input$assocTest,
                                covar = covar,
                                snpStart = start,
                                snpEnd = end)
            })
            
            assoc <- do.call(rbind, assoc.list) #bind results of all chromosomes together
            
        } else { #if the y chromosome is NOT in the data...
            assoc <- assocRegression(genoData = gData(),
                                     outcome = "phenotype",
                                     model.type = input$assocTest,
                                     covar = covars)
        }
        
        #add chromosome position to the results by extracting
        #and merging by snpID
        pos <- data.frame("snpID" = getSnpID(gData()),
                          "BP" = getPosition(gData()))
        
        assoc <- merge(assoc, pos, by = "snpID")
        assoc <- assoc[order(assoc[, "chr"], assoc[, "BP"]), ]
        
        #rename columns to standard PLINK format
        names(assoc)[names(assoc) == "chr"] <- "CHR"
        names(assoc)[names(assoc) == "snpID"] <- "SNP"
        names(assoc)[names(assoc) == "Wald.pval"] <- "P"
        
        #create a file name for saving the results
        sig <- Sys.time() %>% format("%Y%m%d-%H%M%S")
        savePath <- parseDirPath(roots = volumes, selection = input$saveDir) %>% 
            as.character %>% paste("/shinyGWAS.", sig, ".csv", sep = "")
        
        #save results
        assoc %>% write.csv(savePath)
        assoc %>% return
    })
    
    
    # --- UPDATE SIGNIFICANCE VALUE
    # --- according to user selection
    
    sigValue <- reactive({
        
        if (input$multiCompars == "wgs"){
            thresh <- 5e-8
        } else if (input$multiCompars == "bonf"){
            thresh <- 0.05/(nrow(assocData()))
        }
        thresh %>% return
    })
    
    
    # --- PRINT RESULTS OF ASSOCIATION TEST
    # --- as intereactive table
    
    output$gwasSummary <- DT::renderDataTable({
        req(assocData())
        
        datatable(assocData(), 
                  rownames = TRUE, 
                  options = list(scrollX = TRUE)) #allow table to scroll horizontally
    })

    
    ################################################
    ## Tab 3 Display GWAS results and figures     ##
    ################################################
    
    # --- PRINT STATIC MANHATTAN PLOT
    
    # output$manPlot <- renderCachedPlot({
    #     
    #     logs <- -log10(assocData()$P)
    # 
    #     manhattanPlot(p = assocData()$P,
    #                   chromosome = assocData()$CHR,
    #                   signif = sigValue(),
    #                   ylim = range(c(0, max(logs, na.rm = TRUE)+1)))
    # 
    # }, cacheKeyExpr = list(sigValue(), assocData()))
    
    output$manPlot <- renderPlot({
        req(assocData())
        req(sigValue())
        
        gg.manhattan(assocData(), sigValue()) #call function from shinyGWAS_fns.R
    })

    
    # --- PRINT STATIC QQ PLOT
    
    output$qqPlot <- renderCachedPlot({
        par(pty = "s")
        qq(assocData()$P)
        
    }, cacheKeyExpr = list(sigValue(), assocData()))


    ################################################
    ## Tab 4 INTERACTIVE MANHATTAN PLOT           ##
    ################################################
    
    # --- UPDATE CHROMOSOME CHOICES
    
    observeEvent(assocData(), {
        req(assocData())
        
        updateSelectInput(session, inputId = "chr", 
                          choices = c(" " = "", unique(assocData()$CHR)))
    })
    
    
    # --- SAVE SELECTED CHROMOSOME VALUE
    # --- If no selection made, make = NULL
    
    chr <- reactive({
        req(input$chr)
        
        if (input$chr != ""){
            chr <- input$chr
        } else {
            chr <- NULL
        }
        chr %>% return
    })
    
    
    # --- PRINT INTERACTIVE MANHATTAN PLOT
    
    output$manPlot2 <- renderPlotly({
        req(assocData())
        req(sigValue())
        req(input$chr != "")
        
        dat <- assocData()[which(assocData()$CHR == chr()), ]
        
        gg.manhattan(dat, sigValue()) %>% #call function from shinyGWAS_fns.R
            ggplotly(tooltip = "text", )
    })
    
} #END server

# Run the application 
shinyApp(ui = ui, server = server)


