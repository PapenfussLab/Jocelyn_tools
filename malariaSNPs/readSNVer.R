## Reading vcf files output by SNVer, and filtered to subtract variants in 3D7 parent strains.  
## filtering by various parameters including 
## intersection with coding regions, and finding unions and intersections.
## VarScan was run with minimum depth 10 to call a variant, so I will filter by DP first

source("Gilsonfile_paths.R")
snverDir <- file.path(varDir, "snverP09Oct", "removedParentVariants")
library(dplyr)
library(GenomicRanges) 
library(VariantAnnotation)
#library(GenomicFeatures)    ## Needed for makeTxDbfromGFF
library("genomeIntervals")  ## Needed for readGff3, used for reading reference annotation

strain <- c('A', 'B', 'D', 'E')
EsnpVcf <- readVcf(file.path(snverDir, paste0(strain[4], "-parent.filter.vcf")))
View(info(header(EsnpVcf)))  # Table of descriptions of info codes
View(geno(header(EsnpVcf)))  # Table of descriptions of 'geno' codes
head(fixed(EsnpVcf))  # REF and ALT available  if wanted

pf_features <- readGff3(file.path(refDir, "PlasmoDB-29_Pfalciparum3D7.gff"), quiet = TRUE)

filt_vcf <- function(vcf) {
  rownames(colData(vcf)) <- basename(rownames(colData(vcf)))
  ## Filter by: all values of DP > 9
  vcf_depthfilt <- vcf[apply(geno(vcf)$DP>9, 1, all)]
  ## Filter by: at least 1 sample has alternate allele proportion > half
  vcf_AFfilt <- vcf_depthfilt[apply(geno(vcf_depthfilt)$AC/geno(vcf_depthfilt)$DP>0.5, 
                                          1, any)]
  return(vcf_AFfilt)
}

## Find filtered events that overlap features
filtSNPsE <- filt_vcf(EsnpVcf) # nrows 5
snpFeats <- findOverlaps(rowRanges(filtSNPsE), GRanges(pf_features)) # 1 hit

pf_features[subjectHits(snpFeats)]
info(filtSNPsE[unique(queryHits(snpFeats))])
geno(filtSNPsE[unique(queryHits(snpFeats))])$DP
geno(filtSNPsE[unique(queryHits(snpFeats))])$AC

#####################
snpVcf <- readVcf(file.path(snverDir, paste0(strain[1], "-parent.filter.vcf")))
filtSNPsA <- filt_vcf(snpVcf)  # nrows 2
writeVcf(filtSNPsA, filename=file.path(snverDir, "A_filt_DP10_AF5.vcf"))
snpFeatsA <- findOverlaps(rowRanges(filtSNPsA), GRanges(pf_features)) 

pf_features[subjectHits(snpFeatsA)]
info(filtSNPsA[unique(queryHits(snpFeatsA))])
geno(filtSNPsA[unique(queryHits(snpFeatsA))])$DP
geno(filtSNPsA[unique(queryHits(snpFeatsA))])$AC

snpVcf <- readVcf(file.path(snverDir, paste0(strain[2], "-parent.filter.vcf")))
filtSNPsB <- filt_vcf(snpVcf)
snpFeatsB <- findOverlaps(rowRanges(filtSNPsB), GRanges(pf_features)) 

pf_features[subjectHits(snpFeatsB)]  # nrows 3
info(filtSNPsB[unique(queryHits(snpFeatsB))])
geno(filtSNPsB[unique(queryHits(snpFeatsB))])$DP
geno(filtSNPsB[unique(queryHits(snpFeatsB))])$AC
#...
snpVcf <- readVcf(file.path(snverDir, paste0(strain[3], "-parent.filter.vcf")))
filtSNPsD <- filt_vcf(snpVcf)  # nrows 3
snpFeatsD <- findOverlaps(rowRanges(filtSNPsD), GRanges(pf_features)) 

pf_features[subjectHits(snpFeatsD)]
info(filtSNPsD[unique(queryHits(snpFeatsD))])
geno(filtSNPsD[unique(queryHits(snpFeatsD))])$DP
geno(filtSNPsD[unique(queryHits(snpFeatsD))])$AC
#...

#### indels ####
indelVcf <- readVcf(file.path(snverDir, paste0(strain[1], "-parent.indel.filter.vcf")))
indelFiltA <- filt_vcf(indelVcf)  ## nrows 1
indelFeatsA <- findOverlaps(rowRanges(indelFiltA), GRanges(pf_features)) # empty

indelFiltB <- filt_vcf(readVcf(file.path(snverDir, paste0(strain[2], "-parent.indel.filter.vcf"))))  ## empty

indelFiltD <- filt_vcf(readVcf(file.path(snverDir, paste0(strain[3], "-parent.indel.filter.vcf"))))  ## nrows 1
indelFeatsD <- findOverlaps(rowRanges(indelFiltD), GRanges(pf_features)) # empty

indelFiltE <- filt_vcf(readVcf(file.path(snverDir, paste0(strain[4], "-parent.indel.filter.vcf"))))  
## nrows 5, highest AF=0.4
indelFeatsE <- findOverlaps(rowRanges(indelFiltE), GRanges(pf_features)) ## 2 are in CDS features, both PfEMP1

pf_features[subjectHits(indelFeatsE)]
info(indelFiltE[unique(queryHits(indelFeatsE))])
geno(indelFiltE[unique(queryHits(indelFeatsE))])$DP
geno(indelFiltE[unique(queryHits(indelFeatsE))])$AC






