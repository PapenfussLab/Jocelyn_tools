# Jocelyn_tools
This will have various short pieces of code which might be useful to share.

Started on 11 April 2016 with 2 versions of the Python script trim_to_amplicon_primers.

These python scripts were used for cleaning Illumina universal primers from short sequences of DNA for 16S rRNA sequencing for micobiome work. The original version reads a fastq file, the FNA version a fasta file.
Trim_to_amplicon_primers is designed to remove sequencing primers and other extraneous sequence that is outside the amplicon primers. The code has the DNA sequences for the forward and reverse primers hard-coded, so needs to be edited to change them.

trim_fasta_amplicons.py is a variation that removes the amplicon primers as well, leaving just the variable inner section of DNA bases.  
By default it looks for the 805R V4 amplicon and the reverse-complement of the 533F V4 amplicon, and assumes they are in that order. There is an option to use the more usual 5'-3' "forward" order.


read_gridss.R:  An example of reading GRIDSS output into R, and looking for changes in one genome that are not present in the other (very basic)

annotateSNP.py is based on a code fragment by Jan Schroeder. It is designed to read a vcf file of single nucleotide polymorphisms (1-1 changes which do not change the length of the sequence) and report whether they change a coded amino acid, or are non-coding or synonymous.
It reads a gff file, a reference file and a vcf file and outputs a tab-separated list with a single header line starting with #
