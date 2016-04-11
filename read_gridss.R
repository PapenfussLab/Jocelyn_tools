## This reads the vcf-format output of gridss
## Jocelyn Sietsma Penington  September 2015
## Re-running with updates, March 2016. First used 'git pull' at command line 
## in ~/src/gridss/ to get latest version of libgridss.R
setwd("~/src/gridss/src/test/r/")  
source("libgridss.R")

baseDir <- "/wehisan/bioinf/bioinf-data/Papenfuss_lab/projects/malaria/vaccine"
scriptDir <- file.path(baseDir, "analysis_tools/malvacR")
gridssDir <- file.path(baseDir, "strain_comparison/gridss-0.10.1")

library("VariantAnnotation") 
setwd(gridssDir)

vcffn <- file.path(gridssDir, "WTvsQTS.vcf")
WTvsQvcf <- readVcf(vcffn, "")
WTvsQdf <- gridss.vcftodf(WTvsQvcf, allColumns=TRUE)
WTvsQdf <- WTvsQdf[order(-WTvsQdf$QUAL),]  # Sort by QUAL value, descending

View(info(header(WTvsQvcf)))  # Table of descriptions of header codes
View(WTvsQdf)

## Only show events that are not also present in 3D7Q-TS, which is 'Normal'
## RP is ="Count of read pairs supporting breakpoint per category" , i.e. sum of both
## RP0 is count for 'IC=0' where IC is Input Category, which is 3D7Q-TS
## RP1 is count for 'IC=1' which is 3D7WT
## There are 3 main pieces of evidence: RP, SP=split reads, RSR=remote split reads
with(WTvsQdf, View(WTvsQdf[(RP+SR+RSR) > 1  & (RP0+SR0+RSR0) < 0.2 * (RP+SR+RSR) , ]))

## Save the highest quality
write.csv(with(WTvsQdf, WTvsQdf[((RP1+SR1+RSR1) < 0.2 * (RP+SR+RSR)) & (QUAL > 200) , ]),
          file="WTvsQTopQualfromR.csv", quote=FALSE)

## A different criterion, showing fewer columns for easier viewing:
with(WTvsQdf, View(WTvsQdf[(QUAL > 800)| (ASQ > 400) ,c('POS', 'QUAL', 'ASQ', 'BQ', 
                    'RP0', 'SR0', 'RP1', 'SR1', 
                    'HOMLEN',  'INSSEQ') ]))

