
library(magrittr)
library(SNPRelate)
library(GWASTools)
library(SeqArray)


########################################################
## Convert genome files to GDS temp files             ##
## for imputed data                                   ##
########################################################

# imputeToGds <- function(inPath, input.type, dosage, chr){
#   
#   # Create temp file
#   outFile <- tempfile(pattern = paste("shinyGWAS_", 
#                                       as.character(chr)), 
#                       fileext = ".gds")
#   
#   # ONLY NEEDS TO BE CALLED ONCE FOR VCF FILES
#   if (input.type == "vcf"){
#     
#     gds <- seqVCF2GDS(vcf.fn = inPath,
#                       out.fn = tempfile(fileext = ".gds"),
#                       storage.option = "ZIP_RA",
#                       parallel = 2L)
#     
#     seqGDS2SNP(gds, outFile,
#                dosage = TRUE)
#     
#     # # POSSIBLY WORK BUT NEED TO TEST
#     # seqVCF2GDS(vcf.fn = inPath, 
#     #            out.fn = outFile, 
#     #            fmt.import="DS")
#     
#   } else {
#     
#     if (dosage == TRUE){
#       output.type <- "dosage"
#     } else if (dosage == FALSE){
#       output.type <- "genotype"
#     }
#     
#     if (input.type == "imp2"){
#       gen <- inPath[1] %>% as.character
#       samps <- inPath[2] %>% as.character
#       
#       input.files <- c(gen, samps)
#       
#     } else if (input.type == "beagle"){
#       grob.dose <- inPath[1] %>% as.character
#       marker <- inPath[2] %>% as.character
#       
#       input.files <- c(grob.dose, marker)
#       
#     } else if (input.type == "mach"){
#       prob.dose <- inPath[1] %>% as.character
#       info <- inPath[2] %>% as.character
#       snp.pos <- inPath[3] %>% as.character
#       
#       input.files <- c(prob.dose, info, snp.pos)
#     }
#     
#     imputedDosageFile(input.files = input.files,
#                       filename = outFile,
#                       chromosome = chr,
#                       input.type = input.type,
#                       input.dosage = dosage,
#                       file.type = "gds",
#                       output.type = output.type)
#   }
#   
#   outFile %>% return
# }


########################################################
## extract scan df from gds file                      ##
########################################################

extractScan <- function(gds, id2 = NULL, covars = NULL){
  
  # gds <- GdsGenotypeReader(gdsFile)
  
  scanID <- getScanID(gds) %>% as.factor
  
  sex <- getVariable(gds, "sample.annot/sex")
  
  if (is.null(sex)){
    sex <- rep(0, length(scanID))
  }
  
  sex[sex == ""] <- NA
  
  phenotype <- getVariable(gds, "sample.annot/phenotype")
  
  if (is.null(phenotype)){
    phenotype <- rep(0, length(scanID))
  }
  
  # close(gds)
  
  
  df <- data.frame(scanID = scanID,
                   phenotype = phenotype,
                   sex = sex)
  
  if (!is.null(covars)){
    df <- left_join(x = df, y = covars, 
                    by = c("scanID" = id2))
  }
  
  df %>% return
  
}


#######################################################
# Manhattan plot function                            ##
#######################################################
#https://github.com/pcgoddard/Burchardlab_Tutorials/wiki/GGplot2-Manhattan-Plot-Function


# Libraries ====
library(readr)
library(ggrepel)
library(ggplot2)
library(dplyr)
library(RColorBrewer)

gg.manhattan <- function(df, threshold){
  # format df
  df.tmp <- df %>% 
    
    # Compute chromosome size
    group_by(CHR) %>% 
    summarise(chr_len=max(BP)) %>% 
    
    # Calculate cumulative position of each chromosome
    mutate(tot=cumsum(as.numeric(chr_len))-chr_len) %>%
    select(-chr_len) %>%
    
    # Add this info to the initial dataset
    left_join(df, ., by=c("CHR"="CHR")) %>%
    
    # Add a cumulative position of each SNP
    arrange(CHR, BP) %>%
    mutate(BPcum=BP+tot) %>%
    
    # Filter SNP to make the plot lighter
    filter(-log10(P)>0.5)
  
  # set colours
  col <- c("#88CCEE","#44AA99","#117733","#999933",
           "#DDCC77","#CC6677","#882255","#AA4499")
  
  # get chromosome center positions for x-axis
  axisdf <- df.tmp %>% group_by(CHR) %>% summarize(center=(max(BPcum)+min(BPcum))/2)
  
  # Prepare text description for each SNP:
  df.tmp$text <- paste("SNP: ", df.tmp$SNP, "\nPosition: ", df.tmp$BP,
                    "\nChromosome: ", df.tmp$CHR, "\n-log10(P-val):",
                    -log10(df.tmp$P) %>% round(2), sep="")
  
  # Get y limits
  ylims <- c(0, max(-log10(df.tmp$P), na.rm = TRUE)+1)
  
  # Make the plot
  p <- ggplot(df.tmp, aes(x=BPcum, y=-log10(P), text = text)) +
    
    # Show all points
    geom_point(aes(color=as.factor(CHR)), alpha=0.7, size=0.7) +
    scale_color_manual(values = rep(col, 22)) +
    
    # custom X axis:
    scale_x_continuous(label = axisdf$CHR, breaks=axisdf$center) +
    scale_y_continuous(expand = c(0, 0), limits = ylims) + # expand=c(0,0)removes space between plot area and x axis
    
    # add plot and axis titles
    ggtitle(paste0("")) +
    labs(x = "Chromosome") +
    
    # add genome-wide sig lines
    geom_hline(yintercept = -log10(threshold), color = "red", size = 0.3) +
    
    # Custom the theme:
    theme_bw(base_size = 10) +
    theme(
      axis.text = element_text(size = 7),
      plot.title = element_text(hjust = 0.5),
      legend.position="none",
      panel.border = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank()
    )
  
  return(p)
}


























