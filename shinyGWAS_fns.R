library(magrittr)
library(GENESIS)
library(SNPRelate)
library(lme4)
library(GWASTools)
library(GWASdata)
library(qqman)

########################################################
## Read in and save phenotype and covariate data into ##
## a ScanAnnotationDataFrame object                   ##
########################################################

scanObj <- function(data, id, pheno, sex, covars){
  
  if (covars != ""){
    dat <- data.frame(scanID = data[,id],
                      phenotype = data[,pheno],
                      sex = data[,sex],
                      data[,covars])
  } else {
    dat <- data.frame(scanID = data[,id],
                      phenotype = data[,pheno],
                      sex = data[,sex])
  }
  
  scanAnnot <- ScanAnnotationDataFrame(dat)
  scanAnnot$sex <- as.factor(scanAnnot$sex)
  scanAnnot %>% return
}


########################################################
## Find file type of input files                      ##
########################################################

last4Char <- function(string){
  len <- nchar(string)
  substr(string, len-3, len) %>% return
}


########################################################
## Open GDS files                                     ##
########################################################

openGDS <- function(path){
  
  geno <- GdsGenotypeReader(filename = outPath)
  genoData <- GenotypeData(geno, scanAnnot = scan)
  
  return(genoData)
}


########################################################
## Close GDS files (and temp if required)             ##
########################################################

closeGDS <- function(genoData, temp = NULL){
  
  if (temp) {
    close(genoData)
    unlink(temp)
  } else {
    close(genoData)
  }
}


########################################################
## Convert plink format to GDS temp file              ##
########################################################

plinkToGds <- function(inPath) {
  
  outFile <- tempfile("shinyGWAS", fileext = ".gds")
  
  bed <- inPath[1]
  bim <- inPath[2]
  fam <- inPath[3]
  
  snpgdsBED2GDS(bed.fn = bed,
                bim.fn = bim,
                fam.fn = fam,
                out.gdsfn = outFile)
  
  return(outFile)
}




########################################################
## Read in and save genotype data into a              ##
## GenotypeData object                                ##
########################################################

genObj <- function(type, path, scan) {
  if (type == "gds") {
    NULL
    
  } else if (type == "PLINK") {
    bed <- path[1]
    bim <- path[2]
    fam <- path[3]
    
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
