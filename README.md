# Jocelyn_tools
This will have various short pieces of code which might be useful to share.

Started on 11 April 2016 with 2 versions of the Python script trim_to_amplicon_primers.
Subsequently moved to folder amplicon_tools

read_gridss.R:  An example of reading GRIDSS output into R, and looking for changes in one genome that are not present in the other (very basic)

annotateSNP.py is based on a code fragment by Jan Schroeder. It is designed to read a vcf file of single nucleotide polymorphisms (1-1 changes which do not change the length of the sequence) and report whether they change a coded amino acid, or are non-coding or synonymous.
It reads a gff file, a reference file and a vcf file and outputs a tab-separated list with a single header line starting with #

Folder malariaSNPs created Nov 2017. Contains short shell scripts for running SNP callers, and R code for manipulating the outputs.

Folder microbiome_plots has code fragments in R for drawing plots with phyloseq-format microbiome data. As discussed at "Joining The Dots" August 2017.

mouse_microbiome_ASV.sh: mouse microbiome pipeline for making OTU table, as at Nov 2017. A blend of QIIME, Usearch and Vsearch tools

