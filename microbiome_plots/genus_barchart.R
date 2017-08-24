## Using phyloseq data structures to make bar charts for Genus
library("phyloseq")
library(ggplot2)
library(reshape2)
library("RColorBrewer")
theme_set(theme_bw(base_size=8))

## Consistent colours needed for 2 data sets, so hard-coding arrangement of phyla:
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

## List of genus, or higher label where genus undefined,  present in at least 10% in any sample. 
##  Compromise list for showing representative genera in WEHI and BCM data
genusLabelList <- c( # 4 o__Bacteroidales
  'g__Bacteroides', 'g__Prevotella' , 'f__S24-7', 'o__Bacteroidales'
  # 5 o__Clostridiales 
  , 'g__Lachnospira', 'g__Blautia', 'f__Lachnospiraceae', 'g__Faecalibacterium', 'o__Clostridiales' 
  # 3 p__Proteobacteria 
  , 'g__Sutterella', 'g__Haemophilus', 'Proteobacteria' 
  # 1 p__Actinobacteria
  , 'g__Bifidobacterium'
  # 1 p__Verrucomicrobia
  , 'g__Akkermansia'
  )
  
palette_genus <- c(brewer.pal(5, "Reds")[c(4,2,3,5)], #  4 x red for bacter,
                   brewer.pal(6, "Blues")[c(4,2,5,3,6)],   # 5 x blue for firmi
                   brewer.pal(4, "Greens")[c(3,2,4)],   # 3 shades of green for proteo
                   palette_phyl[4:5],  # brown for actino, yellow for verruco
                   palette_phyl[10]  # black for 'other'
)

applyNewGenusLabel <- function(tax){   ## Match the lowest possible taxa name to a name in genusLabelList
  tax_char <- as.vector(tax, mode="character")
  checkrank <- 7
  newRank <- tax_char[checkrank]    # 'lastDefRank' in column 7 is starting point 
  while (! newRank %in% genusLabelList & checkrank >2){ # if not in list, return phylum=Rank1
    checkrank <- checkrank - 1
    newRank <- tax_char[checkrank] 
  }
  newRank
}

genusBarplot <- function(phyl_L6props,
                         phylOrd = phylOrder,
                         Nphyla = 6, 
                         palGenus = palette_genus, 
                         genusOrder = genusLabelList
){
  tax_table(phyl_L6props)[,"Rank2"] <- sapply(tax_table(phyl_L6props)[,"Rank2"],
                                          function(s){if (grepl('__', s)) 
                                            unlist(strsplit(s, '__'))[2] else s})
  tax_table(phyl_L6props)[,'Rank2'] <- apply(tax_table(phyl_L6props), 1, applyNewGenusLabel)
  # 'Rank2' column now holds genus-label value instead of phylum
  psL6summary <- tax_glom(phyl_L6props, 'Rank2')
  summaryprop <- otu_table(psL6summary)
  colnames(summaryprop) <- as.data.frame(tax_table(psL6summary))[, "Rank2"]
  summarypropdf <- data.frame(summaryprop[, intersect(genusOrder, colnames(summaryprop))], 
                              # intersect drops items not in list, keeps in list order, and is OK if any missing
                              # if(taxa_are_rows(summaryprop)) this won't work, need to transpose
                              personID=sample_data(psL6summary)$personID,
                              DayXMethod=paste0(sample_data(psL6summary)$Stool,
                                                sample_data(psL6summary)$Method)  )
  summarypropdf$Otherbacteria <- apply(summaryprop[, setdiff(colnames(summaryprop), genusOrder)], 
                                       1-taxa_are_rows(summaryprop), sum)
  longProp <- melt(summarypropdf, value.name='propAbund',
                   id.vars=c("DayXMethod", "personID"),
                   variable.name='genuslabel')
  bars <- ggplot(longProp, aes(x=DayXMethod, y=propAbund, fill=genuslabel))
  bars + geom_bar(stat="identity") + 
    facet_wrap(~personID) + 
    labs(x="Day and method",
         y="Proportional abundance", 
         title="") +
    scale_fill_manual(values = palGenus, 
                      guide=guide_legend(nrow=4, title = "", keyheight = 0.5)) + 
    theme(axis.text.x = element_text(angle = -90, hjust=0, vjust=0.5)
          , strip.background=element_rect(fill=NA) 
          , legend.key = element_rect(size=NA)
          , legend.position = "bottom"
    )

}


#### Example usage ####
## Load the proportional counts for 72 samples: 
## (bar plots sum values by default - need to have same number of samples in each bar, in this case 1)
load(file.path(localRdata_fp, "Stdcounts_and_props.rdata"), verbose=TRUE) 

## psL6prop was in above file. It is a phyloseq object containing the proportions, 
## combined using tax_glom(psOTUprop, taxrank="Rank6", NArm=FALSE), and then "Rank7" replaced with "lastDefRank"

## Genus barplot (L6), list above of genus, or higher label where genus undefined,
## present in at least 10% in any sample. 
genuslabelbars <- genusBarplot(psL6prop) # other parameters as default
## Save to pdf
plot_name <- file.path(imageDir, "barPropGenus_phylaColours84w.pdf") 
genuslabelbars
ggsave( filename= plot_name,  
        width=84, height=110, units="mm")
