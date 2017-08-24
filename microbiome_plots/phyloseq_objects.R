## Import the biom data using phyloseq, and pre-process.
## This script reads a biom-format OTU table into a phyloseq object, then adds the 
## pynast-aligned tree and the sample metadata from the mapping file.
##
## Future work: parametise so files are inputs, instead of editing!
## Jocelyn Sietsma Penington

library("phyloseq")
library(ggplot2)
theme_set(theme_bw(base_size=18))
## file paths removed

phyloData = import_biom(biom_fp)  
# takes about three minutes. Makes phyloseq object. 

# Tree was made for otu_table filtered to OTUs with count >= 20
phylotree <- read_tree(tree_fp)

## Remove PCR wells with small counts
## Histograms of sample_sums gives guidance as to cut-off.
min(sample_sums(phyloData))
hist(sort(sample_sums(phyloData))[1:30], breaks = 20)
phyloData2 <- prune_samples(sample_sums(phyloData) >= 1000, phyloData)

## Suppress samples with 'No' or 'Very faint' PCR product.
psPCRyf <- prune_samples(sample_data(phyloData2)$PCRsuccess %in% c('Yes', 'Faint'),
                         phyloData2)
phyloData <- psPCRyf
## Suppress small OTUs
minOTUsize <- 20
psOTUge20 <- prune_taxa(taxa_sums(phyloData) >= minOTUsize, phyloData)
psWtree <- merge_phyloseq(psOTUge20, phylotree)

## Merge PCR replications, but keep sample replications. 

## Remove column names that don't apply to merged categories.
varsBySample <- c("SampleNumber", "SampleNumLvl2", "seqRunDate", "personID")
sample_data(otuBySample) <- sample_data(otuBySample)[, varsBySample]
## Restore sample data from the original object.
## For each merged sample, find the sample_data(phyloData) records
for (v in varsBySample) {
  sample_data(otuBySample)[[v]] <-
    sapply(sample_names(otuBySample), function(x) {
      unique(sample_data(psWtree)[
        which(sample_data(psWtree)$SampleNumLvl2== x)][[v]])
    } )
}

## Standardise counts to largest library size (=scaled proportions). Integers
libSizeMax <- max(sample_sums(otuBySample))
standf = function(x) round(libSizeMax * x / sum(x))
psOTUstd <- transform_sample_counts(otuBySample, standf)
## Make a data object of proportions, as well
psOTUprop <- transform_sample_counts(otuBySample, function(x) (x / sum(x)))

## Agglomerate counts with the same taxa to Rank 6, Genus-level
psL6prop <- tax_glom(psOTUprop, taxrank="Rank6", NArm=FALSE)
l6Std <- tax_glom(psOTUstd, taxrank="Rank6", NArm=FALSE)

## Remove the completely Unassigned OTU(s) from standardised genus counts
l6Std <- prune_taxa(as.vector(tax_table(l6Std)[,"Rank1"]!="Unassigned"), l6Std)

## Function to take a single row of a tax-table such as tax_table(l6Std) 
## and return the most detailed known Rank label, with NAs and '__'s appended.
## Added "uncultured" to the test, for Silva taxonomy format
fill_to_defined_level <- function(taxL6) {
  tax_char <- as.vector(taxL6, mode="character")
  if (tax_char[1]=="Unassigned") {
    filledRank <- "Unassigned"    }
  else {
    filledRank <- tax_char[6]    # Starting point is Rank6, genus-label
    for (rank in 6:3) {
      if (is.na(tax_char[rank]) | 
          length(unlist(strsplit(tax_char[rank], "__")))==1 |
          grepl("uncultured", unlist(strsplit(tax_char[rank], "__"))[2]) ) {
        filledRank <- paste(tax_char[rank-1], filledRank, sep=";") }
    }
  }
  filledRank
}
## Use this function to replace "Rank7" column with one filled to last defined label
tax_table(psL6prop)[,7] <- apply(tax_table(psL6prop), 1, fill_to_defined_level)
colnames(tax_table(psL6prop))[7] <- "lastDefRank"

tax_table(l6Std)[,7] <- apply(tax_table(l6Std), 1, fill_to_defined_level)
colnames(tax_table(l6Std))[7] <- "lastDefRank"
