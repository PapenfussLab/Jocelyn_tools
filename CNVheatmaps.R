## Plot heatmaps of inferCNV output, incorporating
## cluster information from tSNE and Seurat
## Starting from Leon di Stefano's script redo_plots.R, 
## copied by Jocelyn SP September 2019

## TODO: add argument parsing, convert to functions of arguments
## instead of loop.
## TODO: optional extra: remove Group facet strip labels and instead use 
## geom_text() to place sample names on top of cluster bar. See:
## https://stackoverflow.com/questions/6455088/how-to-put-labels-over-geom-bar-in-r-with-ggplot2
## for tips
## TODO: optional extra: Move legends to the bottom to make main plot more square

require(tidyverse)
require(cowplot)

# if ( !requireNamespace("infercnv") ){
#   if (!requireNamespace("BiocManager"))
#     install.packages("BiocManager")
#   BiocManager::install( "infercnv" )
#   ## IT had to install 'rjags' for this to work
# }

source( 'filePaths.R' )

#### Cluster colour palette from Yunshun ####
clusterPal <- c( 
  "#1F77B4", "#FF7F0E", "#2CA02C", "#D62728", "#9467BD"
  , "#8C564B", "#E377C2", "#7F7F7F", "#BCBD22", "#17BECF"
  , "#AEC7E8","#FFBB78" ,"#98DF8A", "#FF9896", "#C5B0D5"
  , "#C49C94", "#F7B6D2", "#C7C7C7", "#DBDB8D", "#9EDAE5"
)

#### ggplot themes ####
theme_update(
  panel.background = element_blank(), 
  strip.background =element_blank() , # was  element_rect(fill = NA)
  strip.text = element_text(size = 10 )
  , strip.text.y = element_text(angle = 180)
  )

heatmaptheme <-   theme(
  # No axis marks; image up to the axis edge
  axis.line = element_blank(), 
  axis.ticks = element_blank(), 
  axis.text = element_blank(),
  panel.spacing = unit(0, "lines"),
  panel.border = element_rect( colour = "black", 
                               size = 0.5, linetype = 1,
                               fill = NA ), 
  plot.margin = unit( c(0,0,0,0), "cm" )
) 

## Function to abbreviate cell_group strip labels
justSample <- function(c) {
  return( str_split_fixed( as.character(c), "_", 2)[,1] %>%
            str_replace( ., "Total$", '') )
}

## Inputs and outputs
##### Start with MH0064 paired primary tumour /  lymph node metastasis 
samplegroup <- "MH0064"

for ( samplegroup in # list.files( infercnv_dir) ) {
      c( #"ERTotal",  "Her2"
        "MH0040", "MH0043", "MH0056", "MH0064", "MH0167", "MH0173" 
         ) ) {
  print( samplegroup )
  out_dir <- file.path( plotDir, samplegroup )
  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE)
  }
  
  ## Check for processed copy of data
  localfile <- file.path( localData, 
                          paste0( samplegroup, "_plot_data.Rdata" ) )
  if ( file.exists( localfile )
  ) { 
    print( "Loading local Rdata file" )
    load( localfile )
  } else {
    print( "Reading inferCNV object" )
    require( infercnv )
    infercnv_obj <-
      readRDS( file.path( infercnv_dir, samplegroup, "run.final.infercnv_obj" ) )
    
    ##### Put inferCNV data into tibble / data frames ####
    
    get_grouped_cell_indices <-
      function (x) {
        y <- unlist(map(x, length))
        tibble(
          cell_group = rep(names(y), times = y),
          cell_index      = unlist(x))
      }
    ## Extract names of reference samples. These samples are omitted from plot
    refsamples <- str_split_fixed( 
      names( infercnv_obj@reference_grouped_cell_indices ), 
      pattern = "_", n = 3 )[ , 2 ] %>% 
      unique( . )
    refpattern <- paste0( "^", 
                          paste0( refsamples, sep = "", collapse = "|^") )
    
    cell_tbl <- 
      get_grouped_cell_indices(
        infercnv_obj@observation_grouped_cell_indices
      ) %>%
      mutate(
        barcode = colnames(infercnv_obj@expr.data)[(.)$cell_index],
        cell_order = 1:n()) %>%
      mutate(
        sample_name = str_split_fixed(cell_group, "_", 3)[,2],
        cluster_name = str_split_fixed(cell_group, "_", 3)[,3] 
      ) %>% 
      mutate( 
        cluster_name = factor( cluster_name, 
                               levels = 1:max( as.numeric( cluster_name ) )
        )
      )
    cell_group_levels <- tidyr::expand( 
      cell_tbl, nesting( sample_name, cluster_name) ) %>% 
      arrange( cluster_name, sample_name) 
    cell_group_levels <- 
      with( cell_group_levels, 
            paste( sample_name, cluster_name, sep = "_")
      )
    ## Previously all tumor cell_group values started with "tumor_"
    ## As reference cells now omitted, this has been dropped
    cell_tbl <- mutate( cell_tbl,
                        cell_group = factor( 
                          str_replace( cell_group, "malignant_", "" ),
                          levels = cell_group_levels ) )
    ## Remove cells which are in sample (x) cluster groups with low counts
    ## That is, cells from a sample which is present in a cluster at a low level
    ## Replaced hard-coded minimum with mix of hard and relative
    # samXclustcrit <- 50
    ## Smaller of: 5% of largest sample (x) cluster
    ##                          100
    samXclustmax <- count( cell_tbl, sample_name, cluster_name ) %>%
      select( n ) %>% max( )
    samXclustcrit <- min( samXclustmax/20, 100 )
    cell_tbl <- cell_tbl %>% add_count( sample_name, cluster_name ) %>%
      filter( n > samXclustcrit ) %>%
      mutate( cell_order = 1:n(),
              n = NULL )
    
    gene_tbl <- 
      infercnv_obj@gene_order %>% 
      as_tibble() %>% 
      mutate(gene_name = rownames(infercnv_obj@expr.data)) %>%
      arrange(chr, start) %>%
      group_by(chr) %>%
      mutate(gene_order = 1:n()) %>%
      ungroup()
    
    expr_data <-
      infercnv_obj@expr.data %>% 
      as.data.frame() %>% 
      rownames_to_column("gene_name") %>% 
      ## Remove reference sample expressions
      select( -matches( refpattern ) ) %>% 
      gather(key = "barcode", value = "mod_expr", -gene_name) %>% 
      as.tbl()
    
    plot_data <-
      expr_data %>%
      left_join(gene_tbl) %>%
      right_join(cell_tbl) 
    
    lower_cutoff <- 1 - sd(plot_data$mod_expr) * 3
    upper_cutoff <- 1 + sd(plot_data$mod_expr) * 3
    ## Clip the expression so the isolated high values
    ## don't wash out all others in the heatmap palette
    plot_data$clip_expr <- pmax( pmin( upper_cutoff, plot_data$mod_expr ), 
                                 lower_cutoff )
    
    ## Save a local copy
    save( plot_data, cell_tbl, samXclustcrit, file = localfile )
  }
  print( table( cell_tbl$sample_name, cell_tbl$cluster_name ) )
  # }
  
  use_minidata <- FALSE
  if ( use_minidata ) {
    ## Cut-down data for experimentation
    ## Genes in 3 chromosomes x last 10% cells 
    usedata <- plot_data %>% filter(  chr %in% c("chr20", "chr21", "chr22" ) &
                                         cell_order > max( cell_order) * 0.90
    )
    hist( usedata$mod_expr )
    
    usecell <- cell_tbl  %>% filter( cell_order > max( cell_order) * 0.90 )
  } else {
    usedata <- plot_data
    usecell <- cell_tbl
  }

  #### Make and save plots ####
  ## Parameters for heatmap legend:  
  lowexpr <- ceiling( min( usedata$clip_expr ) * 100 ) / 100
  topexpr <- floor( max( usedata$clip_expr ) * 100 ) / 100
  
  ## Labeller function to add count of genes per chrom to facet strip label
  geneCounts <- usedata %>% group_by( chr ) %>% 
    summarise( n_genes = n_distinct( gene_name )  ) 
  geneCountLabs <- paste0(  geneCounts$chr , '\n',
                           geneCounts$n_genes, ' genes' )
  names( geneCountLabs ) <- geneCounts$chr
  
  p_rect <-
    usedata %>%
    ggplot(aes( )) + 
    ## because facetting by chrom has chr1 at the top to chr22 at the bottom
    ## we want genes arranged in reverse order up the y-axis
    geom_rect( aes(
      xmin = cell_order - 0.5,  xmax = cell_order + 0.5,
      ymin = -gene_order - 0.5, ymax = -gene_order + 0.5,
      fill = clip_expr )
    ) +
    facet_grid( chr ~ cell_group 
               , scale = "free"
               ## no change to 'space': all facets same size
               , switch = "y"  # move chrom facet labels to left
               , labeller = labeller( chr = geneCountLabs )
    ) +
    scale_fill_gradient2(low = "#00008b", high = "#8b0000", midpoint = 1,
                         name = "Expr", 
                         breaks = c( lowexpr, 1, topexpr ),
                         labels = c(  paste("<", lowexpr), "1.00",
                                      paste(">", topexpr) )
    ) +
    heatmaptheme +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0))
  if ( use_minidata ) { print( p_rect ) }
  legend_heatmap <- get_legend(p_rect + 
                                 theme(legend.position = "left") 
  )
  p_rect  <- p_rect + theme(legend.position='none') + 
    theme( strip.background.x = element_blank(),
           strip.text.x = element_blank() )
  
  #### Add coloured strip for clusters
  ## Named vector to add cell_group counts to labels
  cellgpcount <- summary( usecell$cell_group)[ which( summary( usecell$cell_group) > 0)]
  groupCountLabs <- paste0( justSample( names( cellgpcount ) ),
    '\n', cellgpcount, ' cells'
    )
  names( groupCountLabs ) <- names( cellgpcount )
    
  clusterbarp <- ggplot( usecell, aes( x = cell_order, y = 1, fill=cluster_name ) ) + 
    geom_tile() + scale_fill_manual(
      values = clusterPal[ 
        levels(usecell$cluster_name) %in% unique(usecell$cluster_name)
        ], 
      name = "Cluster" ) + 
    facet_grid( . ~ cell_group  
               , scale = "free"
               , labeller = labeller( cell_group = groupCountLabs )
    ) +
    heatmaptheme  + 
    scale_x_continuous(expand = c(0, 0) ) +
    scale_y_continuous(expand = c(0, 0), breaks = 1 ) + 
    labs( x = 'Cluster' )
  legend_clusterbar <- get_legend(clusterbarp)
  clusterbarp  <- clusterbarp + theme(legend.position='none',
                                      axis.title = element_blank() ) +
    theme( strip.background.y = element_rect( colour = NA )
            )
  
  ## Assemble with cowplot
  v <- 0
  main_plot  <- plot_grid( clusterbarp, p_rect, align = "v", axis = "lr",
                           rel_heights = c( 1, 12 ), ncol = 1 )
  legend_plot = plot_grid( legend_clusterbar, legend_heatmap, ncol = 1 )
  full_plot <- plot_grid( main_plot, legend_plot, rel_widths = c(9, 1), nrow = 1 )
  v <- v + 1
  plotfilename <- paste0("cnv_", "cluster_", 
                         format(Sys.Date(), "%d%b"),
                         "_v", v, ".jpg")
  # ggsave(plot = full_plot, filename =  file.path( out_dir, plotfilename ), 
  #        width=210, height=297, units="mm" )
  
  ### Simpler save after mysterious failure: 
  ### "Error in grid.newpage() : could not open file ... "
  jpeg( file = file.path( out_dir, plotfilename ), 
       width=210, height=297, units="mm", res = 300 )
  print( full_plot )
  dev.off()
  
}

