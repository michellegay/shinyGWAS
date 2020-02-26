
# --- INSTALL PACKAGES

regPkgs = c('shiny','shinybusy','magrittr',
            'qqman', 'summarytools', 'GGally',
            'shinyFiles', 'DT', 'dplyr',
            'shinythemes')

# --- bioconductor packages
bioPkgs = c('SNPRelate', 'GWASTools', 'SeqArray')

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(version = "3.10")

for(p in regPkgs){
    if(!require(p, character.only = TRUE)){
        install.packages(p)
    } 
    library(p, character.only = TRUE)
}

for(p in bioPkgs){
    if(!require(p, character.only = TRUE)){
        BiocManager::install("GWASTools")
    } 
    library(p, character.only = TRUE)
}





# --- Set size of allowed file uploads.
gb <- 3 #number of GB
options(shiny.maxRequestSize = (gb*1000)*1024^2)

ui = fluidPage(
    navbarPage(
        theme = shinythemes::shinytheme("simplex"),
        "shinyGWAS",
        
        ################################################
        ## Tab 0 Instructions and specifications      ##
        ################################################
        
        tabPanel(title = "INFO", fluid = TRUE,
                 titlePanel("User Reference"),
                 
                 
                 
                 ), #END tab
        
        
        
        ################################################
        ## Tab 1 Uploading Phenotype & Covariate data ##
        ################################################
        
        tabPanel(title = "UPLOAD DATA", fluid = TRUE,
                 titlePanel("Upload Data"),
                 br(),
                 
                 # --- Places a spinning circle in the corner of the page to
                 # --- show when the app is busy computing.
                 add_busy_spinner(spin = "fading-circle"),
                 
                 fluidRow( 
                     sidebarPanel(
                            
                         # --- Side bar with two tabs for uploading data
                            tabsetPanel(
                                # --- Sidebar tab #1 where user uploads genome data
                                tabPanel(title = "Step one", fluid = TRUE,
                                         br(),
                                         h4("Upload genome data"),
                                         br(),
                                         
                                         # --- Input: User checks box if using imputed data.
                                         # --- Default selection is FALSE (not imputed)
                                         checkboxInput(inputId = "isImpute", 
                                                       label = "Imputed data", 
                                                       value = FALSE),
                                         
                                         # --- If data is imputed, user is asked to specify whether imputed data
                                         # --- is recorded as dosages or genotype probabilities.
                                         # --- Default selection is dosages.
                                         conditionalPanel(
                                             condition = "input.isImpute == 1",
                                             br(),
                                             
                                             radioButtons(inputId = "isDosage",
                                                          label = NULL,
                                                          choices = c("Dosage" = TRUE,
                                                                      "Genotype Probabilities" = FALSE),
                                                          selected = TRUE),
                                             
                                         ), #END panel
                                         
                                         # --- Input: User selects the file format of genome data.
                                         tags$hr(),
                                         radioButtons(inputId = "fileType",
                                                      label = "Select file format",
                                                      choices = c("GDS" = "gds",
                                                                  "VCF" = "vcf",
                                                                  "PLINK" = "plink",
                                                                  "IMPUTE2" = "IMPUTE2",
                                                                  "BEAGLE" = "BEAGLE",
                                                                  "MaCH" = "MaCH")),
                                         
                                        br(),
                                         
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
                                tabPanel(title = "Step two", fluid = TRUE,
                                         br(),
                                         h4("Upload sample data"),
                                         br(),
                                         
                                         # --- Input: User checks box if sample data has
                                         # --- no column names
                                         radioButtons(inputId = "hasHeader",
                                                      label = NULL,
                                                      choices = c("With column headers" = 1,
                                                                  "No column headers" = 0),
                                                      selected = 1),
                                         
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
                                         
                                         # --- Input: User checks to automatically re-format ID column
                                         # --- to sampleID123_sampleID123 format.
                                         checkboxInput(inputId = "formatId", 
                                                       label = "Re-format ID", 
                                                       value = FALSE),

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
                                         
                                         br(),
                                         
                                         # --- User selects this button to upload data and see summary
                                         # --- information.
                                         actionButton(inputId = "genoDo",
                                                      label = "Upload data"),
                                         
                                         br(),br(),
                                         textOutput("uploadComplete"),
                                         br(),
                                         
                                 ) #END tab
                            ), #END tabset
                            ), #END column
                     
                     # --- Right-hand panel which displays summary figures and information
                     # --- about the uploaded genome and sample data.
                     mainPanel(
                            
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
                            
                            tabsetPanel(
                                # --- tab #0 gives instructions for uploading data.
                                # --- Instructions will update depending on the 
                                # --- type of data uploaded
                                tabPanel(title = "Instructions", fluid = TRUE,
                                         
                                         column(6,
                                                br(),
                                                h4("Step One: Upload genome data"),
                                                tags$hr(),
                                                
                                                h5("1.  Imputed data:"),
                                                p(paste("If you are uploading imputed data,",
                                                        "check the 'Imputed data' box then specify",
                                                        "whether imputed data is represented as",
                                                        "dosages or genotype probabilities.")),
                                                
                                                br(),
                                                h5("2.  Data format:"),
                                                p(paste("Select the appropriate format that genome data",
                                                        "are stored in.")),
                                                
                                                br(),
                                                h5("3.  Upload genome data files:"),
                                                p(paste("If data are stored across",
                                                        "multiple file types (e.g. PLINK format files), upload",
                                                        "each file type in the correct upload field as directed.")),
                                                p(paste("If data are stored across multiple files of the same type",
                                                        "(e.g. data are saved by chromosome), hold 'Ctrl' to select",
                                                        "all relevant files.")),
                                                
                                                p(paste("NOTE: When uploading multiple files AND",
                                                        "multiple file types, the number of files",
                                                        "uploaded in each upload field must be the same",
                                                        "and correspond to the same data."))
                                                
                                                ), #END column
                                         
                                         column(6, 
                                                br(), 
                                                h4("Step Two: Upload sample data"),
                                                tags$hr(),
                                                
                                                h5("1.  Upload sample data file:"),
                                                
                                                # --- If genome data is VCF format, show instructions for
                                                # --- compulsary upload of additional sample information.
                                                conditionalPanel(
                                                    condition = "input.fileType == 'vcf'",
                                                    
                                                    p(paste("VCF files must be accompanied by",
                                                            "a text file in .txt or .csv format which contain",
                                                            "at minimum, the following fields:")),
                                                    p("- A unique sample ID, and"),
                                                    p("- Phenotype data"),
                                                    p(paste("NOTE: the sample ID must correspond to",
                                                            "the genotype data.")),
                                                    
                                                ), #END panel
                                                
                                                # --- If genome data is NOT vcf format, show instructions for
                                                # --- optional upload of additional sample data.
                                                conditionalPanel(
                                                    condition = "input.fileType != 'vcf'",
                                                    
                                                    p("-- OPTIONAL --"),
                                                    p(paste("Upload any additional sample data",
                                                            "for inclusion in association model.")),
                                                    p(paste("Data must be in .csv or .txt format and must",
                                                            "include a unique sample ID field which",
                                                            "corresponds to the genotype data.")),
                                                    
                                                ), #END panel
                                                
                                                p("Indicate whether the file contains column headers."),
                                                
                                                br(),
                                                h5("2.  Select variables:"),
                                                p("-- IGNORE IF NO SAMPLE DATA SELECTED --"),
                                                p(paste("From the 'ID' drop-down menu, select the column",
                                                        "that corresponds to the unique sample ID.")),
                                                p(paste("Check the 'Reformat ID' box to re-format the ID column to be in",
                                                        "'IDxyz_IDxyz' format.")),
                                                p(paste("If the sample data contains phenotype or sex data",
                                                        "select the corresponding columns from the respective",
                                                        "drop-down menus.")),
                                                p("NOTE: Sex data, if present, must be formatted as 'M' and 'F'"),
                                                p(paste("If the data contain any additional covariates (e.g.",
                                                        "principle components), that you wish you include in the model",
                                                        ", select these from the 'Other covariates' menu.",
                                                        "Click and drag or hold 'Ctrl' to select multiple covariates.")),
                                                
                                                br(),
                                                h5("3.  Upload data"),
                                                p("Click 'Upload data' button to confirm settings and upload data.")
                                                
                                                ) #END column
                                ), #END tab
                                
                                
                                # --- tab #1 shows descriptive information about the 
                                # --- data in a table format.
                                # --- htmlOutput() shows table as it would appear when 
                                # --- knitted in R, as opposed to the shiny table format.
                                tabPanel(title = "Data summary", fluid = TRUE,
                                         br(), br(),
                                         
                                         htmlOutput(outputId = "phenoSummary")
                                         
                                         ), #END tab
                                
                                # --- tab #2 shows a correlation plot with histograms.    
                                tabPanel(title = "Correlation", fluid = TRUE,
                                         br(), br(),
                                         
                                         plotOutput(outputId = "heatmap")
                                         
                                         ), #END tab
                                
                                # --- tab #3 shows the files that have been uploaded and their size.
                                # --- DT::dataTableOutput returns an interactive table.
                                tabPanel(title = "File summary", fluid = TRUE,
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
        
        tabPanel("RUN TEST", fluid = TRUE,
                 titlePanel("Select Model and Run Test"),
                 br(),
                 
                 # --- Places a spinning circle in the corner of the page to
                 # --- show when the app is busy computing.
                 add_busy_spinner(spin = "fading-circle"),
                 
                 fluidRow(
                     # --- Sidebar where association test options
                     # --- are displayed.
                     sidebarPanel(
                            
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
                                         label = "Select correction method",
                                         choices = c("Gemone-wide significance (5e-8)" = "wgs",
                                                     "Bonferroni at alpha = 0.05 (alpha/No. SNPs)" = "bonf"),
                                         selected = "bonf"),
                            
                            br(),
                            
                            # --- User is asked to select a folder where the results of 
                            # --- association testing and the summary report will be
                            # --- automatically saved when association is complete.
                            shinyDirButton(id = "saveDir", 
                                           label = "Select a folder",
                                           title = "Select a folder"),
                            
                            br(), br(),
                            
                            # --- Shows the folder path that has been selected.
                            verbatimTextOutput("directorypath"),
                            
                            br(), 
                            
                            # --- User clicks this button to confirm settings and 
                            # --- begin the association test.
                            actionButton(inputId = "assocDo",
                                         label = "Run test"),
                            
                            br(), br(),
                            textOutput("testComplete"),
                            br(),
                            
                            # --- User can download a summary report of the analyses
                            downloadButton("report", "Download report")
                            
                            ), #END column
                     
                     # --- Right-hand panel where results of association test will be displayed.
                     mainPanel(
                         tabsetPanel(
                             # --- tab #0 gives instructions for selecting model settings.
                             tabPanel(title = "Instructions", fluid = TRUE,
                                      br(),
                                      h4("Step Three: Select model settings"),
                                      tags$hr(),
                                      
                                      h5("1.  Select model:"),
                                      p("Indicate the type of phenotype data to be analysed."),
                                      p("From the available options, select the association model to be used."),
                                      
                                      br(),
                                      h5("2.  Select p-value correction:"),
                                      p(paste("Select the equation for adjusting the signifiance",
                                              "threshold to correct for multiple comparisons.")),
                                      p(paste("- Genome-wide significance: a standard threshold",
                                              "set to 5e-8 (0.00000005)")),
                                      p(paste("- Bonferroni: calculated by dividing alpha (set to 0.05)",
                                              "by the number of SNPs to be analysed. i.e. (0.05)/(no. SNPS).")),
                                      
                                      br(),
                                      h5("3.  Save results:"),
                                      p(paste("Select a folder where results of association test will be saved.")),
                                      
                                      br(),
                                      h5("4.  Run test:"),
                                      p("Click 'Run test' button to confirm settings and run association analyses."),
                                      
                                      br(),
                                      h5("5.  Download report:"),
                                      p(paste("Once test is complete you may choose to download a summary report of the",
                                        "test by clicking 'Download report'. The report details test settings and",
                                        "includes all summary tables and figures."))
                                      
                             ), #END tab
                             
                             tabPanel(title = "Significant SNPs", fluid = TRUE,
                                      
                                      br(),
                                      h4("Significant SNPs"),
                                      br(),
                                      
                                      # --- Shows the significant SNPs of the association test
                                      # --- in an interactive table.
                                      DT::dataTableOutput(outputId = "gwasSigs")
                                      
                            ) #END tab
                         ) #END tabset
                     ) #END main panel
                 ) #END row
        ),#END tab
        
        ################################################
        ## Tab 3 Display GWAS results and figures     ##
        ################################################
        
        tabPanel("FIGURES", fluid = TRUE,
                 titlePanel("Manhattan and QQ Plots"),
                 br(),
                 
                 # --- Places a spinning circle in the corner of the page to
                 # --- show when the app is busy computing.
                 add_busy_spinner(spin = "fading-circle"),
                 
                 fluidRow(
                     tags$hr(),
                     br(),br(),
                     column(7,
                            
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
        ## Tab 4 GWAS RESULTS                         ##
        ################################################
        
        tabPanel("RESULTS", fluid = TRUE,
                 titlePanel("GWAS Results"),
                 br(),
                 
                 fluidRow(
                     tags$hr(),
                     br(),
                     
                     column(1,
                            ), #END column
                     
                     column(10,
                            # --- Shows all results of the association test
                            # --- in an interactive table.
                            DT::dataTableOutput(outputId = "gwasSummary")
                            ), #END column
                     
                     column(1,
                            )
                     
                 ) #END row
        ) #END tab
        
        
        
        
        ################################################
        ## Tab 5 Interactive Manhattan Plot           ##
        ################################################
        # 
        # --- CURRENTLY DISABLED BECAUSE PLOTLY SEEMS TO SHOW ONLY
        # --- A SUBSET OF ALL POINTS AND CANNOT BE TRUSTED!
        # 
        # tabPanel("Interactive Manhattan", fluid = TRUE,
        #          titlePanel("Interactive Manhattan Plot"),
        #          
        #          fluidRow(
        #              column(4,
        #                     tags$hr(),
        #                     
        #                     # --- Places a spinning circle in the corner of the page to
        #                     # --- show when the app is busy computing.
        #                     add_busy_spinner(spin = "fading-circle"),
        #                     
        #                     # --- User selects the chromosomes they wish to show on plot.
        #                     # --- Multiple chromosomes can be selected.
        #                     # --- Chromosome selection allows plot to render faster, especially
        #                     # --- if user is only interested in a few chromosomes.
        #                     selectInput(inputId = "chr",
        #                                 label = "Select Chromosome(s)",
        #                                 multiple = TRUE,
        #                                 choices = character(0),
        #                                 selectize = FALSE),
        #                     
        #              ), #END column
        #          ), #END row
        #          
        #          fluidRow(
        #              tags$hr(),
        #              column(11,
        #                     
        #                     # --- Show interactive Manhattan plot.
        #                     plotlyOutput(outputId = "manPlot2")
            #                 
            #          ) #END column
            #      ) #END row
            # ) #END tab
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
                validate(
                    need(input$isDosage == TRUE, 
                         message = "Imputed data must be represented as dosages.")
                )
                
                temp <- tempfile()
                path <- seqVCF2GDS(vcf.fn = inPath, 
                                   out.fn = temp) %>% 
                    seqGDS2SNP(out.gdsfn = outFile, dosage = TRUE)
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
    #         paths <- imputeToGds(inPath,
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
    # --- Automatically detect .csv or .txt format
    # --- Return an error if neither of these formats
    # --- If data has no headers, creates header as 
    # --- column number + the value in the first row
    # --- e.g. 1_val1, 2_val2 etc.
    
    sData <- reactive({
        req(input$sampFile)
        
        path <- input$sampFile$datapath %>% as.character
        ext <- substr(x = path, 
                      start = (nchar(path)-2), 
                      stop = nchar(path)) %>% tolower()
        
        validate(
            need(ext == "txt" | ext == "csv", 
                 message = "File must be .csv or .txt format.")
        )
        
        if (input$hasHeader == 1){
            header <- TRUE
        } else if (input$hasHeader == 0){
            header <- FALSE
        }
        
        if (ext == "csv"){
            dat <- read.csv(path, header = header)
        } else {
            dat <- read.table(path, header = header)
        }
        
        if (input$hasHeader == 0){
            n <- c(1:ncol(dat)) %>% as.character()
            
            for (ii in 1:length(n)){
                n[ii] <- paste(n[ii], dat[1, ii], sep = "_")
            }
            
            names(dat) <- n
        }
        
        dat %>% return
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
            updateSelectInput(session, inputId = "sex",
                              choices = c("N/A" = "", names(sData())),
                              selected = "")
     })
    
    observeEvent(sData(), {
        updateSelectInput(session, inputId = "covars", 
                          choices = c("N/A" = "", names(sData())), 
                          selected = "")
    })
    

    # --- OPEN GENOME DATA IN READER OBJECT
    
    gReader <- eventReactive(input$genoDo, {
        req(gdsPath())
        GdsGenotypeReader(filename = gdsPath())
    })
    
    
    # --- COLLATE AND SAVE SAMPLE DATA TO SCANANNOTATION OBJ
    # --- Depending on genome data format, will extract 
    # --- sample data from genome data file, if additional
    # --- sample data has been provided this will be used instead.
    
    scanAnnot <- eventReactive(input$genoDo, { #run this when file upload button is selected
        req(gReader())
        
        if (input$fileType == "vcf"){ #if genome data in VCF format...
            req(sData())
            req(input$id)
            req(input$pheno)
            
            id <- sData()[,input$id] %>% as.factor
            pheno <- sData()[,input$pheno]
            
            df <- data.frame(scanID = id,
                             phenotype = pheno)
            
            if (input$sex != ""){
                sex <- sData()[,input$sex]
                df <- data.frame(df, sex = sex)
            }
            
            if (input$covars != ""){
                covars <- sData()[,input$covars]
                df <- data.frame(df, covars)
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
                
                df <- extractScan(gReader(), id, covars) #call function from shinyGWAS_fns.R
                
                if ("XX_pheno_XX" %in% colnames(df)){ #replace extracted phenotype with selected
                    df$phenotype <- df$XX_pheno_XX
                    df$XX_pheno_XX <- NULL
                }
                
                if ("XX_sex_XX" %in% colnames(df)){ #replace extracted sex with selected
                    df$sex <- df$XX_sex_XX
                    df$XX_sex_XX <- NULL
                }
                
            } else { #if no additional sample file has been selected...
                df <- extractScan(gReader())
            }
        }
        
        # --- Extract the sample IDs from genome data, the sample data must match these
        genomeIDs <- getScanID(gReader()) %>% as.factor() %>% as.data.frame()
        names(genomeIDs) <- "scanID"
        
        # --- Force ID column to match format of VCF files of UWA imputed data
        if (input$formatId == 1){
            df$scanID <- paste(df$scanID, "_", df$scanID, sep = "") %>% as.factor()
        }
        
        # --- MIGHT RESULT IN REMOVAL OF SAMPLES HERE
        # --- SHOULD THINK OF WAY TO DOCUMENT REMOVED
        # --- SAMPLES AND INCLUDE IN SUMMARY REPORT
        if (nrow(genomeIDs) > nrow(df)){
            df <- full_join(genomeIDs, df, by = "scanID")
        } else if (nrow(genomeIDs) < nrow(df)){
            df <- inner_join(genomeIDs, df, by = "scanID")
        }
        
        validate(
            need(nrow(df) == nrow(genomeIDs),
                 message = "ERROR: Sample ID's do not match genome data.")
        )
        
        scan <- df %>% ScanAnnotationDataFrame
        
        if ("sex" %in% names(scan)){
            scan$sex <- scan$sex %>% as.factor
        }
        
        scan %>% return
    })
    
    
    # --- MERGE GENOME AND SAMPLE DATA AS GENOTYPEDATA OBJECT
    
    gData <- eventReactive(input$genoDo, {
        req(scanAnnot())
        
        gReader() %>% 
            GenotypeData(scanAnnot = scanAnnot())
        
    })
    
    
    # --- INDICATE THAT GENOME DATA UPLOAD WAS SUCCESSFUL
    
    output$uploadComplete <- renderText({
        req(gData())
        "Upload complete."
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
    # --- Save to a separate parameter so can be
    # --- passed to summary report
    
    uploadedFiles <- eventReactive(input$genoDo, {
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
    
    output$filesUploaded <- DT::renderDataTable({
        req(input$genoDo)
        req(input$fileType)
        
        uploadedFiles()
        
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
    
    dirPath <- reactive({
        req(input$saveDir)
        
        if (is.integer(input$saveDir)) {
            cat("No directory selected.")
        } else {
            parseDirPath(roots = volumes, 
                         selection = input$saveDir)
        }
    })
    
    output$directorypath <- renderPrint({
        req(dirPath())
        
        dirPath()
    })
    
    
    # --- RUN ASSOCIATION ANALYSIS AND SAVE RESULTS
    
    runTest <- eventReactive(input$assocDo, { #run association when user clicks button
        req(gData())
        req(input$saveDir)
        
        startTime <- Sys.time()
        
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
        
        #make p-value column numeric (for some reason it
        #becomes a factor during imputation analysis...)
        #have to go to character first to avoid becoming integer sequence
        assoc$P <- assoc$P %>% as.character() %>% as.numeric()
        
        #create a file name for saving the results
        sig <- Sys.time() %>% format("%Y%m%d-%H%M%S")
        savePath <- parseDirPath(roots = volumes, selection = input$saveDir) %>% 
            as.character %>% paste("/shinyGWAS.", sig, ".csv", sep = "")
        
        #save results
        assoc %>% write.csv(savePath)
        
        endTime <- Sys.time()
        
        out <- list(results = assoc, start = startTime, end = endTime, path = savePath)
        out %>% return
    })
    
    
    # --- SAVE ASSOCIATION DATA TO OWN VARIABLE
    
    assocData <- reactive({
        req(runTest())
        
        runTest()$results %>% return()
    })
    
    
    # --- INDICATE THAT ASSOCIATION TEST WAS SUCCESSFUL
    
    output$testComplete <- renderText({
        req(assocData())
        "Test complete."
    })
    
    # --- GENERATE SUMMARY REPORT
    # --- Code from https://shiny.rstudio.com/articles/generating-reports.html
    
    output$report <- downloadHandler(
        
        filename = "report.html",
        content = function(file) {
            # 
            # Copy the report file to a temporary directory before processing it, in
            # case we don't have write permissions to the current working dir (which
            # can happen when deployed).
            tempReport <- file.path(tempdir(), "testhtml.Rmd")
            file.copy("testhtml.Rmd", tempReport, overwrite = TRUE)
            
            # Set up parameters to pass to Rmd document
            params <- list(imputed = input$isImpute,
                           dosages = input$isDosage,
                           file.format = input$fileType,
                           files = uploadedFiles(),
                           scan.annot = scanAnnot(),
                           dat.type = input$phenoDatType,
                           model.type = input$assocTest,
                           correction = input$multiCompars,
                           save.folder = dirPath(),
                           sig.snps = sigSNPs(),
                           sig.val = sigValue(),
                           results = assocData(),
                           start.time = runTest()$start,
                           end.time = runTest()$end,
                           results.file = runTest()$path)
            
            # Knit the document, passing in the `params` list, and eval it in a
            # child of the global environment (this isolates the code in the document
            # from the code in this app).
            rmarkdown::render(tempReport, output_file = file,
                              params = params,
                              envir = new.env(parent = globalenv())
            )
        }
    )
    
    
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
    
    
    # --- PRINT SIGNIFICANT RESULTS OF ASSOCIATION TEST
    # --- as intereactive table
    
    sigSNPs <- eventReactive(input$assocDo, {
        req(assocData())
        req(sigValue())
        
        cols <- c("SNP", "CHR", "P")
        dat <- assocData()[which(assocData()$P <= sigValue()), cols]
        
    })
    
    output$gwasSigs <- DT::renderDataTable({
        req(sigSNPs())
        
        datatable(sigSNPs(),
                  rownames = TRUE,
                  options = list(scrollX = TRUE))
    })
    
    
    ################################################
    ## Tab 3 Display GWAS figures                 ##
    ################################################
    
    # --- PRINT STATIC MANHATTAN PLOT
    
    output$manPlot <- renderCachedPlot({

        logs <- -log10(assocData()$P)

        manhattanPlot(p = assocData()$P,
                      chromosome = assocData()$CHR,
                      signif = sigValue(),
                      ylim = range(c(0, max(logs, na.rm = TRUE)+1)))

    }, cacheKeyExpr = list(sigValue(), assocData()))
    
    
    # --- PRINT STATIC QQ PLOT
    
    output$qqPlot <- renderCachedPlot({
        par(pty = "s")
        qq(assocData()$P)
        
    }, cacheKeyExpr = list(sigValue(), assocData()))
    
    
    ################################################
    ## Tab 4 ALL RESULTS TABLE                    ##
    ################################################
    
    # --- PRINT SIGNIFICANT RESULTS OF ASSOCIATION TEST
    # --- as intereactive table
    output$gwasSummary <- DT::renderDataTable({
        req(assocData())

        datatable(assocData(),
                  rownames = TRUE,
                  options = list(scrollX = TRUE)) #allow table to scroll horizontally
    })
    
    
    ###############################################
    ## Tab 5 INTERACTIVE MANHATTAN PLOT           ##
    ################################################
    #
    # --- CURRENTLY DISABLED BECAUSE PLOTLY SEEMS TO SHOW ONLY
    # --- A SUBSET OF ALL POINTS AND CANNOT BE TRUSTED!
    #
    # # --- UPDATE CHROMOSOME CHOICES
    # 
    # observeEvent(assocData(), {
    #     req(assocData())
    #     
    #     updateSelectInput(session, inputId = "chr", 
    #                       choices = c(" " = "", unique(assocData()$CHR)))
    # })
    # 
    # 
    # # --- SAVE SELECTED CHROMOSOME VALUE
    # # --- If no selection made, make = NULL
    # 
    # chr <- reactive({
    #     req(input$chr)
    #     
    #     if (input$chr != ""){
    #         chr <- input$chr
    #     } else {
    #         chr <- NULL
    #     }
    #     chr %>% return
    # })
    # 
    # 
    # # --- PRINT INTERACTIVE MANHATTAN PLOT
    # 
    # output$manPlot2 <- renderPlotly({
    #     req(assocData())
    #     req(sigValue())
    #     req(input$chr != "")
    #     
    #     dat <- assocData()[which(assocData()$CHR == chr()), ]
    #     
    #     gg.manhattan(dat, sigValue()) %>% #call function from shinyGWAS_fns.R
    #         ggplotly(tooltip = "text", )
    # })
    
} #END server

# Run the application 
shinyApp(ui = ui, server = server)


