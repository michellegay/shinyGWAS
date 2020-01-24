
library(magrittr)
library(GENESIS)
library(SNPRelate)
library(lme4)
library(GWASTools)
library(GWASdata)
library(qqman)

########################################################
## Read in and save phenotype and covariate data from ##
## a .CSV file into a ScanAnnotationDataFrame object  ##
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
## Open GDS files                                     ##
########################################################

openGDS <- function(path, scan = NULL){
  
  geno <- GdsGenotypeReader(filename = path)
  genoData <- GenotypeData(geno, scanAnnot = scan)
}

########################################################
## Convert plink format to GDS temp file              ##
########################################################

plinkToGds <- function(inPath) {
  
  outFile <- tempfile("shinyGWAS", fileext = ".gds")
  
  bed <- inPath[1] %>% as.character
  bim <- inPath[2] %>% as.character
  fam <- inPath[3] %>% as.character
  
  snpgdsBED2GDS(bed.fn = bed,
                bim.fn = bim,
                fam.fn = fam,
                out.gdsfn = outFile)
  
  return(outFile)
}

########################################################
## Perform association test by Linear Mixed Model     ##
########################################################

gwasNull <- function(genoData, fam){
  
  covars <- (genoData %>% getScanAnnotation %>% getVariableNames)[-c(1,2)]
  
  gaus <- c("ord", "cont")
  
  if (fam %in% gaus){
    null.model <- genoData %>% getScanAnnotation %>% 
      fitNullModel(outcome = "phenotype",
                   covars = covars,
                   family = gaussian)
    
  } else if (fam == "bin"){
    null.model <- genoData %>% getScanAnnotation %>% 
      fitNullModel(outcome = "phenotype",
                   covars = covars,
                   family = binomial)
  } else {
    null.model <- NULL
  }
  
  null.model %>% return
}

########################################################
## Perform association test by Linear Mixed Model     ##
########################################################

gwasMixed <- function(genoData, null.model) {
  
  iter <- GenotypeBlockIterator(genoData = genoData,
                                snpBlock = 10000)
  
  assoc <- assocTestSingle(iter, 
                           null.model = null.model)
  
  # Only returning some of the results from the null model
  assoc %>% return
}

########################################################
## Perform association test by Linear or Logistic     ##
## Regression                                         ##
########################################################

gwasLogLin <- function(genoData, model.type) {
  # model.type = "linear" or "logistic"
  
  covars <- (genoData %>% getScanAnnotation %>% getVariableNames)[-c(1,2)]
  
  chr <- genoData %>% getChromosome
  
  assoc.list <- unique(chr) %>% lapply(function(x) {
    ## copied from GWASTools vignette
    ## Y chromosome only includes males, cannot have sex as covariate
    covar <- ifelse(x == 25, yes = c(""), no = covars)
    
    if (covar == ""){
      covar <- NULL
    }
    
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



