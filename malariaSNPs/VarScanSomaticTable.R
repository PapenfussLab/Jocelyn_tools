## Script to read multiple VarScan output files (in VarScan native format), 
## extract the somatic variants and make a table of all affected genes

## J Sietsma Penington January 2017
## Copied for GilsonT210617 files October 2017

source("Gilsonfile_paths.R")
library(dplyr)

list.files(file.path(varDir))
setwd(file.path(varDir, 'varscan02Oct'))
## read gff file and work with genomic intervals 
library("genomeIntervals")
pf_features <- readGff3(file.path(refDir, "PlasmoDB-29_Pfalciparum3D7.gff"), quiet = TRUE)


## destring function copied from : https://github.com/gsk3/taRifx/blob/master/R/Rfunctions.R#L1161
## Used to convert character percentages to numeric. 
## Only number-like characters are kept. Anything which does not contain numerals returns NA
destring <- function(x,keep="0-9.-") {
  return( as.numeric(gsub(paste("[^",keep,"]+",sep=""),"",x)) )
}

sigSomatic <- function(fname, pcrit=0.05) {
    sigsom <- 
      read.table(fname, header = TRUE, sep='\t') %>% 
      mutate(normal_var_freq = 
               as.numeric(sub('%', '', normal_var_freq))/100,
             tumor_var_freq = 
               as.numeric(sub('%', '', tumor_var_freq))/100) %>%
      filter(somatic_status=='Somatic' & somatic_p_value <= pcrit) 
    ## Not interested in strandedness, so discard all variable _plus or _minus
    sigsom <- sigsom[, !grepl('_plus|_minus', names(sigsom)) ]
    ## Discard a few more columns
    sigsom[, c('somatic_status', 'variant_p_value')] <- NULL
  return(sigsom)
}

## Make a GenomeIntervals object from the VarScan results.
## SNPs and Insertions have a point position. 
## Column 'var' starts with + for insertion, - for deletion.  
## Deletions end at position+lengthofdeletion which is nchar-1 because the 1st char is '-'. 
nativeVarscan.gi <- function(varscan_df){
  GenomeIntervals(chromosome=varscan_df$chrom, 
                             start=varscan_df$position,
                             end=ifelse(
                               substr(as.character(varscan_df$var),1,1) == "-",
                               varscan_df$position - 1 + 
                                 nchar(as.character(varscan_df$var)),
                               varscan_df$position) 
  )
}

filtByCDS <- function(var_df, feat=pf_features) {
  CDS_gff <- feat[annotation(feat)$type == "CDS", ] 
  vars_in_CDS <- interval_overlap(nativeVarscan.gi(var_df), CDS_gff)
  # This gives a list of same length as var_df, with a 0-length int vector if 
  # the variant is not in an CDS, and the relevant index or indexes in feat if it is.
  ## Code returns the first if there is more than one feature for one variant.
  filt_var <- var_df[sapply(vars_in_CDS, function(x) !length(x)==0), ]
  CDSofVars <- CDS_gff[na.omit(sapply(vars_in_CDS, function(x) x[1])), ] 
  varParents <-sub("[.].*", "", getGffAttribute(CDSofVars, "Parent") ) # remove any characters after a dot
  affectedGeneDesc <- sapply(varParents, function(vP) {
    getGffAttribute(feat[getGffAttribute(feat, "ID")==vP], # find records with ID matching a parent,
                    c("ID", "description"))})              # and extract 2 fields
  return(cbind(filt_var, t(affectedGeneDesc)) )
}

strain <- c('A', 'B', 'D', 'E')

## concatenate all 'somatic' variants for strain 
for (s in strain) {
  varfilenames <- list.files(pattern=paste0("mergeStoS4_", s))
  
  sigSom <- plyr::ldply(varfilenames, 
                        function(f) data.frame(sigSomatic(f), Sample = strsplit(f, "[[:punct:]]")[[1]][3])
  )
  
  # filter to only those in CDS regions, and add gene information
  sigSomCDS <- filtByCDS(sigSom)
  sigSomCDS <- sigSomCDS[order(sigSomCDS$ID, sigSomCDS$Sample),]
  # filter by tumor_var_freq
  sigSomCDS[sigSomCDS$tumor_var_freq>0.75,]
  
  write.table(sigSomCDS[sigSomCDS$tumor_var_freq>0.75,], sep='\t', quote=FALSE, 
              file=paste0(s,"_somatic_p05_tvf75_CDS.tsv"), row.names=FALSE)
  
}
