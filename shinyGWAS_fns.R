
library(magrittr)
library(GENESIS)
library(SNPRelate)
library(lme4)
library(GWASTools)
library(GWASdata)
library(qqman)
#library(SeqArray) # should be imported by GENESIS

########################################################
## Read in and save phenotype and covariate data from ##
## a .CSV file into a df                              ##
########################################################

buildScan <- function(id, pheno, sex, covars1, id2 = NULL, covars2 = NULL){
  
  df <- data.frame(scanID = id,
                   phenotype = pheno,
                   sex = sex)

  if (covars1 != ""){
    if (!is.null(covars2)){
      df <- merge(x = data.frame(df,covars1), 
                  y = covars2,
                  by.x = "scanID",
                  by.y = id2)
    } else {
      df <- data.frame(df,covars1)
    }
  } else {
    if (!is.null(covars2)){
      df <- merge(x = df, 
                  y = covars2,
                  by.x = "scanID",
                  by.y = id2)
    }
  }
  
  df %>% return
}

########################################################
## extract scan df from gds file                      ##
########################################################

extractScan <- function(gdsFile, id2 = NULL, covars = NULL){
  
  gds <- GdsGenotypeReader(gdsFile)
  
  scanID <- getScanID(gds)
  sex <- getVariable(gds, "sample.annot/sex")
  sex[sex == ""] <- NA
  phenotype <- getVariable(gds, "sample.annot/phenotype")
  
  close(gds)
  
  df <- data.frame(scanID = scanID,
                   phenotype = phenotype,
                   sex = sex)
  
  if (!is.null(covars)){
    df <- left_join(x = df, y = covars, 
                    by = c("scanID" = id2))
  }
  
  df %>% return
  
}

########################################################
## Save phenotype and covariate data  into a          ##
## ScanAnnotationDataFrame object                     ##
########################################################

scanObj <- function(df){
  scanAnnot <- df %>% ScanAnnotationDataFrame
  scanAnnot$sex <- scanAnnot$sex %>% as.factor
  scanAnnot %>% return
}  


########################################################
## Open GDS files                                     ##
########################################################

openGDS <- function(path, scan = NULL){
  
  geno <- GdsGenotypeReader(filename = path)
  genoData <- geno %>% GenotypeData(scanAnnot = scan)
}

########################################################
## Convert genome files to GDS temp files             ##
## for non-imputed data                               ##
########################################################

convertToGds <- function(inPath, type) {
  
  # Create temp file
  outFile <- tempfile(pattern = "shinyGWAS", 
                      fileext = ".gds")
  
  if (type == "plink"){
    bed <- inPath[1] %>% as.character
    bim <- inPath[2] %>% as.character
    fam <- inPath[3] %>% as.character
    
    snpgdsBED2GDS(bed.fn = bed,
                  bim.fn = bim,
                  fam.fn = fam,
                  out.gdsfn = outFile, 
                  family = TRUE, 
                  cvt.chr = "int")
    
  } else if (type == "vcf"){
    snpgdsVCF2GDS(vcf.fn = inPath, 
                  out.fn = outFile)
  }

  outFile %>% return
}

########################################################
## Convert genome files to GDS temp files             ##
## for imputed data                                   ##
########################################################

imputeToGds <- function(inPath, input.type, dosage, chr){
  
  # Create temp file
  outFile <- tempfile(pattern = paste("shinyGWAS_", 
                                      as.character(chr)), 
                      fileext = ".gds")
  
  # ONLY NEEDS TO BE CALLED ONCE FOR VCF FILES
  if (input.type == "vcf"){
    seqVCF2GDS(vcf.fn = inPath, 
               out.fn = outFile, 
               fmt.import="DS")
    
  } else {
    
    if (dosage == TRUE){
      output.type <- "dosage"
    } else if (dosage == FALSE){
      output.type <- "genotype"
    }
    
    if (input.type == "imp2"){
      gen <- inPath[1] %>% as.character
      samps <- inPath[2] %>% as.character
      
      input.files <- c(gen, samps)
      
    } else if (input.type == "beagle"){
      grob.dose <- inPath[1] %>% as.character
      marker <- inPath[2] %>% as.character
      
      input.files <- c(grob.dose, marker)
      
    } else if (input.type == "mach"){
      prob.dose <- inPath[1] %>% as.character
      info <- inPath[2] %>% as.character
      snp.pos <- inPath[3] %>% as.character
      
      input.files <- c(prob.dose, info, snp.pos)
    }
    
    imputedDosageFile(input.files = input.files,
                      filename = outFile,
                      chromosome = chr,
                      input.type = input.type,
                      input.dosage = dosage,
                      file.type = "gds",
                      output.type = output.type)
  }
  
  outFile %>% return
}

########################################################
## Perform NULL association test                      ##
########################################################

# gwasNull <- function(genoData, fam){
#   
#   covars <- (genoData %>% getScanVariableNames)[-c(1,2)]
# 
#   if (length(covars) == 0){
#     covars <- NULL
#   }
#   
#   if (fam == "quant"){
#     family <- "gaussian"
#     
#   } else if (fam == "bin"){
#     family <- "binomial"
#   }
#   
#   null.model <- genoData %>%
#     fitNullModel(outcome = "phenotype",
#                  covars = covars,
#                  family = family)
#   
#   null.model %>% return
# }

########################################################
## Perform association test by Linear Mixed Model     ##
########################################################
# 
# gwasMixed <- function(genoData, null.model) {
#   
#   iter <- GenotypeBlockIterator(genoData = genoData,
#                                 snpBlock = 10000)
#   
#   assoc <- assocTestSingle(iter, 
#                            null.model = null.model)
#   
#   assoc %>% return
# }

########################################################
## Perform association test by Linear or Logistic     ##
## Regression                                         ##
########################################################

gwasLogLin <- function(genoData, model.type){
  # model.type = "linear" or "logistic"
  
  covars <- (genoData %>% getScanVariableNames)[-c(1,2)]
  
  if (length(covars) == 0){
    covars <- NULL
  }
  
  chr <- genoData %>% getChromosome
  
  if (25 %in% chr){
    assoc.list <- unique(chr) %>% lapply(function(x) {
      ## copied from GWASTools vignette
      ## Y chromosome only includes males, cannot have sex as covariate
      
      # If sex is the only covariate
      if (length(covars) == 1){
        covar <- ifelse(x == 25, yes = c(""), no = covars)
        
        if (covar == ""){
          covar <- NULL
        }
        
        # If more covariates  
      } else if (length(covars) > 1){
        y <- covars[-which(covars == "sex")]
        covar <- ifelse(x == 25, yes = y, no = covars)
        
      } else {
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
    
  } else {
    assoc <- assocRegression(genoData = genoData,
                             outcome = "phenotype",
                             model.type = model.type,
                             covar = covars)
  }
  
  assoc %>% return
}

#######################################################
# Add chromosome position to the GWAS results table  ##
# and sort by chromosome + position                  ##
#######################################################

addPosition <- function(genoData, assoc){
  
  pos <- data.frame("snpID" = getSnpID(genoData),
                    "pos" = getPosition(genoData))
  
  merged <- merge(assoc, pos, by = "snpID")
  merged <- merged[order(merged[, "chr"], merged[, "pos"]), ]
  merged %>% return
}




























