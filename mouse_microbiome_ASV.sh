## Pipeline for 16S microbiome analysis for mouse stool from Alan Yu
## Started August 2017 by Jocelyn Sietsma Penington

## Pipeline decisions: I'm not going to filter by alignment position using Mothur
## I will cluster using UParse
## CHANGE 30 Oct : I will not cluster, 
## I will 'denoise' and use Amplicon Sequence Variants instead of 97% similar clusters
## I will assign taxonomy using Silva

## Define some path shortcuts
export PAPDIR=/wehisan/bioinf/bioinf-data/Papenfuss_lab/projects
export BASEDIR=$PAPDIR/metagenomics/seth/AlanYu2017
export DATADIR=$BASEDIR/sequence_data
export OUT1DIR=$BASEDIR/01_processed_data
export TOOLDIR=$BASEDIR/analysis_tools
export SILVADIR=$PAPDIR/reference_genomes/bacterial16S/silva/SILVA123_QIIME_release
export GGDIR=/usr/local/bioinfsoftware/qiime/qiime_v2-1.8.0/qiime-deploy/qiime_software/

## Quality inspection of raw data:
mkdir $DATADIR/rawQC
fastqc -o $DATADIR/rawQC --noextract $DATADIR/*fastq.gz 

# Uncompress
gzip -cd $DATADIR/MISEQ2367_S1_L001_R1_001.fastq.gz > $DATADIR/MISEQ2367_R1.fastq
gzip -cd $DATADIR/MISEQ2367_S1_L001_R2_001.fastq.gz > $DATADIR/MISEQ2367_R2.fastq

# How many sequences contain Illumina universal forward primer ?
python $TOOLDIR/find_amplicon_primers.py $DATADIR/MISEQ2367_R1.fastq GTGACCTATGAACTCAGGAGTC 
## Reverse Illumina adapter :
python $TOOLDIR/find_amplicon_primers.py $DATADIR/MISEQ2367_R2.fastq CTGAGACTTGCACATCGCAGC -rc
## each present in majority of sequences

## Confirming amplicon primers
F341=CCTACGGGNGGCWGCAG
R805=GACTACHVGGGTATCTAATCC
python $TOOLDIR/find_amplicon_primers.py $DATADIR/MISEQ2367_R1.fastq $F341 -rc
python $TOOLDIR/find_amplicon_primers.py $DATADIR/MISEQ2367_R2.fastq $R805 -rc
## present in 14.6mil (R1) and 14.1mil (R2) sequences in normal orientation, almost none in rc.

## Merge overlapping paired-end reads with PEAR. 
## Options: -v is min overlap, -m is max assembled length, -n is min assembled length
## -q is quality threshold (not used), -j is number of threads
module load pear
pear -f $DATADIR/MISEQ2367_R1.fastq -r $DATADIR/MISEQ2367_R2.fastq  \
     -v 50 -m 600 -n 300 -j $(nproc) -o $DATADIR/MISEQ2367

cd $OUT1DIR
ln -s $DATADIR/MISEQ2367.assembled.fastq MISEQ2367.assembled.fastq

lamboot 

## Use Qiime extract_barcodes.py to remove the 8-mer barcodes that are at each end of the 
## merged sequences. 
## extract_barcodes.py -f $OUT1DIR/MISEQ2367.assembled.fastq   \
##  -o bar_exed -c barcode_paired_stitched -l 8 -L 8  \
##  -m mapping.txt 

## split_libraries_fastq: label sequences with sample ID based on index sequences. 
split_libraries_fastq.py   \
   -i $OUT1DIR/bar_exed/reads.fastq   \
   -b $OUT1DIR/bar_exed/barcodes.fastq -m $OUT1DIR/mapping.txt \
   --barcode_type 16 -q 29 -n 1 -o $OUT1DIR/labelled_hiqual 

## Remove universal primers
cd $OUT1DIR/labelled_hiqual/
python $TOOLDIR/trim_fasta_amplicons.py -i seqs.fna -o trimmed_seqs.fna

## Edit sequence headers to format suitable for usearch:
## QIIME format is <sample_id>_<seq_counter> , and we need just <sample_id> at start
## The seqid including underscore is added to the end.
perl -pe '$_ =~s />(.+?)(_\d+)(.*$)/>$1$3 seqid=$1$2/'   \
    $OUT1DIR/labelled_hiqual/trimmed_seqs.fna >  \
    $OUT1DIR/labelled_hiqual/trimmed_seqs_sampleID.fna

OUT2DIR=$BASEDIR/02_uSearch_OUT
mkdir $OUT2DIR
## As I am using uSearch pipeline, need to "dereplicate" i.e. extract unique sequence set
## The input sequences to denoise2 must be a set of unique sequences sorted in order 
## of decreasing abundance with size annotations in the labels. 
## Use all but 2 CPU cores
nohup vsearch --derep_full $OUT1DIR/labelled_hiqual/trimmed_seqs_sampleID.fna  \
   --output $OUT2DIR/unique_w_sizesV.fasta --sizeout --threads $(( $(nproc)-2)) \
    > nohup_vsearch.out &

## Pick Operational Taxonomic Units using uSearch  ####
## Parameters:  97% similarity, usearch - replaced 30 Oct 2017

# nohup usearch -cluster_otus $OUT2DIR/unique_w_sizesV.fasta -minsize 2  \
#    -otus $OUT2DIR/otu_rep_set.fasta -relabel OTU > nohup_cluster_otus.out &

#### Use usearch unoise2 (unoise3 not available in v9.2, which is what we have avail)
#### Accept defaults of -minampsize 4 ,  -unoise_alpha 2.0
nohup usearch -unoise2 $OUT2DIR/unique_w_sizesV.fasta  \
 -fastaout $OUT2DIR/denoised.fasta &

## Replace semicolon with space, after OTU ID in OTU header
perl -pe '$_ =~s /(>Otu\d+)(;)(.*)/$1 $3/' $OUT2DIR/denoised.fasta  \
   > $OUT2DIR/denoised_OTU_ID.fasta
   
## Assign sequences to OTUs, with 99% cut-off. Ties, with equal identity %,  
## are assigned to the 1st match, which will be the largest matching OTU as 
## unique_w_sizesV.fasta was size-sorted large to small

vsearch --usearch_global $OUT1DIR/labelled_hiqual/trimmed_seqs_sampleID.fna  \
   --db $OUT2DIR/denoised_OTU_ID.fasta --strand plus --id 0.99  \
   --otutabout $OUT2DIR/otutab99.txt \
   --biomout $OUT2DIR/otutab99.biom   

## Assign taxonomy using Silva 16S bacterial database (and QIIME)
OUT3DIR=$BASEDIR/03_OTUs_w_meta
## mkdir $OUT3DIR
TAXO_FP=$SILVADIR/taxonomy/16S_only/99/majority_taxonomy_7_levels.txt
REFSEQ=$SILVADIR/rep_set/rep_set_16S_only/99/99_otus_16S.fasta

## merge and delete step fails if there are too many processors - 
## set to 10 instead of nproc-2
nohup parallel_assign_taxonomy_uclust.py -i $OUT2DIR/denoised_OTU_ID.fasta \
   -o $OUT3DIR/denoised_Silva123 -O10 -t $TAXO_FP -r $REFSEQ  \
   -v  \
   > nohup_assign_silva.out &

## Taxonomy to OTU mapping file needs a header line
sed  -i '1 i\#OTU_ID\ttaxonomy\tconsensus_fraction\tnum_accepts'  \
   $OUT3DIR/denoised_Silva123/denoised_OTU_ID_tax_assignments.txt

## Add sample and taxonomy data to OTU table
biom add-metadata -i $OUT2DIR/otutab99.biom  \
   -o $OUT3DIR/denoised_allmeta.biom   \
   --sample-metadata-fp $OUT1DIR/mapping.txt  \
   --observation-metadata-fp $OUT3DIR/denoised_Silva123/denoised_OTU_ID_tax_assignments.txt  \
   --sc-separated taxonomy

## Make a phylogenetic tree. 
## As there are now  ~700 OTUs, I will not discard small first
## To make phylogenetic tree requires first aligning to ref using secondary structure
TEMPLATE=$SILVADIR/core_alignment/core_alignment_SILVA123.fasta 
nohup parallel_align_seqs_pynast.py -i $OUT2DIR/denoised_OTU_ID.fasta  \
  -t $TEMPLATE -O10 \
  -o $OUT3DIR/pynast_aligned_denoised/ \
  > nohup_align_silva.out &

nohup make_phylogeny.py -i $OUT3DIR/pynast_aligned_denoised/denoised_OTU_ID_aligned.fasta \
   -l $OUT3DIR/fasttree.log &

## Post-processing of OTU table ####
cd $OUT3DIR/
## Merge by sample (combining PCR wells/barcodes) i.e. mouse ID
collapse_samples.py -b $OUT3DIR/denoised_allmeta.biom  \
	-m $OUT1DIR/mapping.txt \
	--collapse_fields mouseID \
	--output_biom_fp denoised_byMouse.biom \
	--output_mapping_fp $OUT1DIR/mapping_byMouse.txt
## Inspect otu_table statistics
biom summarize-table -i denoised_allmeta_byMouse.biom -o denoised_byMouse.summary
more denoised_byMouse.summary 


