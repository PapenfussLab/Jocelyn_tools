## Using phyloseq data structures to make abundance charts for Phylum, Order or Genus
## Box-and-whiskers plots.

library("phyloseq")
library(ggplot2)
library(reshape2)
library("RColorBrewer")
theme_set(theme_bw(base_size=8))

Ntop <- 10
palette_phyl <- c(brewer.pal(9, "Set1"), "#000000", brewer.pal(9, "Pastel1"))  
## extended for more than 9 phyla. 
## Swap some colours: Red, blue, green, brown, yellow, purple, orange
palette_phyl <- palette_phyl[c(1:3, 7, 6, 4, 5, 8:length(palette_phyl))]

## Consistent colours needed for 2 data sets, so hard-coding display order for phyla. 
phylOrder <- c( "Bacteroidetes",  
               "Firmicutes",     
               "Proteobacteria", 
               "Actinobacteria", 
               "Verrucomicrobia",
               "Tenericutes",    
               "Cyanobacteria",  
               "Euryarchaeota",  
               "Lentisphaerae",  
               "Fusobacteria",   
               "TM7",            
               "Deferribacteres")

## For same reason, hardcode orderOrder
## I want a mix of highest average proportion and highest maximum proportion
## The 12 that have a proportion in some sample greater than 2% makes a good list 
tax_table(psL4prop)[which(apply(otu_table(psL4prop), 2, max) >=0.02), c(2,4)]

ordOrder <- c( # 1 x bacteriodetes, 2 x firmicutes
  "Bacteroidales", "Clostridiales", "Erysipelotrichales", 
  # 4 x proteobacteria
  "Burkholderiales", "Pasteurellales", "RF32", "Enterobacteriales", 
  # 1 actinobacteria 
  "Bifidobacteriales", 
  # 1 verrucomicrobia 
  "Verrucomicrobiales", 
  # 2 x tenericutes 
  "Anaeroplasmatales", "RF39", 
  #1 cyanobacteria
  "YS2" )
palette_order <- c(palette_phyl[1], #  red for bacter,
                   brewer.pal(3, "Blues")[3:2],   # 2 x blue for firmi
                   brewer.pal(5, "Greens")[5:2],   # 4 shades of green, dark to pale
                   palette_phyl[4:5],  # brown for actino, yellow for verruco
                   brewer.pal(3, "Purples")[3:2], # 2 x purple for tenericutes
                   palette_phyl[7]     # orange for cyanobacteria
)

orderBoxplot <- function(phyl_L6counts,
                         phylOrd = phylOrder , 
                         palPhyl = palette_phyl) {
  psL4propNo0 <- transform_sample_counts(tax_glom(phyl_L6counts, taxrank="Rank4") , 
                                         function(x) ( (x+1)/sum(x) ))
  tax_table(psL4propNo0)[,"Rank2"] <- sapply(tax_table(psL4propNo0)[,"Rank2"],
                                             function(s){unlist(strsplit(s, '__'))[2]})
  tax_table(psL4propNo0)[,"Rank4"] <- sapply(tax_table(psL4propNo0)[,"Rank4"],
                                             function(s){unlist(strsplit(s, '__'))[2]})
  ## This uses the phyloseq 'plot_bar' function to set the ggplot data and aesthetics. 
  ## Two boxplots are added, 1st with border and point colours from Phylum, including 
  ## the default solid circle outliers. The 2nd has black borders, open circle outliers 
  ## and coloured fill. This has the effect of colouring outliers to match bars.
  ## The original barplot layer is then removed
  pbox <- plot_bar(psL4propNo0, x="Rank4") + scale_y_log10() + 
    geom_boxplot(aes(colour=Rank2)) + 
    geom_boxplot(aes(fill=Rank2), outlier.shape=21) + 
    labs(title="Proportion of bacterial Orders in each sample (log scale)", y="Proportion",
         x=element_blank()) +
    scale_fill_manual(values=palPhyl, name="Phylum") + 
    scale_colour_manual(values=palPhyl, name="Phylum") + 
    theme(axis.text.x = element_text(vjust=0.5))
  pbox$layers <- pbox$layers[-1]
  ## X-axis goes in order of decreasing total sum of proportions
  barOrder <- tax_table(psL4propNo0)[names(sort(taxa_sums(psL4propNo0), TRUE)),"Rank4"]
  pbox$data$Rank4 <- factor(pbox$data$Rank4, levels=barOrder) 
  pbox$data$Rank2 <- factor(pbox$data$Rank2, levels=phylOrd) 
  pbox  
}

#### Example usage ####
load(file.path(localRdata_fp, "Stdcounts_and_props.rdata"), verbose=TRUE) 

## l6Std was in above file. It is a phyloseq object containing the standardised counts, 
## combined using tax_glom(psOTUprop, taxrank="Rank6", NArm=FALSE), 
## and then "Rank7" replaced with "lastDefRank". Code is in phyloseq_objects.R

### Order (L4) box-and-whiskers, with log-scale ####
boxes <- orderBoxplot(phyl_L6counts = l6Std) 
print(boxes)
ggsave( filename= file.path(imageDir, "boxPropOrderColourPhylumLvl2.pdf"), 
        width=297, height=210, units="mm")
