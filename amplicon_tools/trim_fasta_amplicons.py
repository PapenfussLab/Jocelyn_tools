from __future__ import print_function
"""
trim_amplicon_primers.py -i <inputfastaFilename> -o <outputfastailename>

Author: J Sietsma Penington
Original Date: 13 Nov 2014

Remove the 16S rRNA gene amplicon primer from each end of a fasta file of overlapped sequences 
assembled from paired-end reads. 
From the beginning of the line these sequences have 
possibly an 8-base barcode (which may have been removed, for example by QIIME's extract_barcodes.py),
followed by the 22-mer Illumina forward universal primer (which may not be complete),
then the V4 forward amplicon primer.
At the end there is the reverse-complement of the V4 reverse amplicon primer,
followed by the 21-mer reverse universal primer (which may not be complete),
optionally followed by an 8-base barcode
The program writes a fasta file with all bases before and including the V4 forward primer, and 
all bases from the start of the V4 reverse primer, removed.
Primer may be present with errors. Rather than implementing an inexact search, if a primer 
is not found the program removes the most common length of sequence found before or after primers 

A sequence is written for each input sequence. Sequences deemed too short for trimming are 
written unchanged.

To Do:
1. Import proper sequence handling from Bio.seq OR use reverse-complement function from Mungo
2. Import an inexact search algorithm?
3. Allow input or reading of amplicon primers

---------------------------------------
Changes for ENDIA background QC project  7 Dec 2015
---------------------------------------
Region is V4, amplicon primers as given by Steven Wilcox are:
805R GGACTACNVGGGTWTCTAAT
533F GTGYCAGCMGCCGCGGTAA
but direction is reversed: - reads contain 805R as supplied, followed further down by 
reverse-complement of 533F

10 Dec - changed to edit fasta, aka fna, instead of fastq
10 Aug 2015 - changed to remove amplicon primers as well as bases outside them
            - changed to add option to not reverse-complement amplicon primers, for use with 
            reverse-complemented sequences. That is, ...533F...seq...rc(805R)...
            - added framework, unimplemented yet, to specify region
To do: add fastq option

"""

import sys
import argparse
import re
from collections import Counter

def rc(seq): 
	complement= {'A':'T','T':'A','C':'G','G':'C', '[':']', ']':'['}
	c_seq = ''
	for b in seq:
		if b in complement: 
			b = complement[b]
		c_seq = c_seq + b
	return c_seq[::-1]

def trim_amplicons(in_name, out_name, **options_dict):
	if ('seqDirection' not in options_dict) or (options_dict['seqDirection']==None): 
		options_dict['seqDirection'] = 'Reverse'
	if ('regionV' not in options_dict) or (options_dict['regionV']==None): 
		options_dict['regionV'] = 'V4'
	f805R = "GGACTAC.[ACG]GGGT[AT]TCTAAT"
	f533F = "GTG[CT]CAGC[AC]GCCGCGGTAA"
#	f341F_V3 = "CCTACGGG.GGC[AT]GCAG" # copied from mouse work, which used V3-V4
	print "seqDirection:", options_dict['seqDirection'], " regionV:", options_dict['regionV']
	if (options_dict['regionV']=='V4'):
		if options_dict['seqDirection'][0] in ['R', 'r']:
			primer1 = re.compile(f805R, re.I)
			primer2 = re.compile(rc(f533F), re.I)
		else:
			primer1 = re.compile(f533F, re.I)
			primer2 = re.compile(rc(f805R), re.I)
	print 'primer1:', primer1.pattern, ', primer2:', primer2.pattern
	bothfound = 0
	pre_length = Counter(); pre_trim = 0
	post_length = Counter(); post_trim = 0
	out_fa = open(out_name, 'w')
	
	with open(in_name, 'r') as fasta:
	# Read through file twice, first counting primer positions:
		for line in fasta:
			if line[0] =='>':
				seqheader = line
				seq = fasta.next()
				
				pf_find = re.search(primer1, seq)
				if pf_find: 
					# Forward primer found: increment pre_length
					pre_length[pf_find.end()] += 1
				pr_find = re.search(primer2,seq)
				if pr_find:
					# Reverse primer found: increment post_length
					tail = len(seq) - pr_find.start() - 1
					post_length[tail] += 1
				if pf_find and pr_find: 
					bothfound += 1
		if len(pre_length) > 0:
			pre_trim = pre_length.most_common(1)[0][0]
		if len(post_length) > 0:
			post_trim = post_length.most_common(1)[0][0]
		print ("Forward primers found:", sum(pre_length.values()),"Counts of pre_length:", pre_length)
		print ("Number of bases trimmed from sequence start when forward primer not found:", pre_trim)
		print ("Reverse primers found:", sum(post_length.values()),"Counts of post_length", post_length)
		print ("Number of bases trimmed from sequence end when reverse primer not found:", post_trim)
		print ("Both primers found:", bothfound)
		
	# Read through file again, trimming sequences 
		fasta.seek(0)
		for line in fasta:
			if line[0] =='>':
				seqheader = line
				seq = fasta.next()
				
				pf_find = re.search(primer1, seq)
				if pf_find: 	# If forward primer is found, trim it and preceding bases
					start = pf_find.end()
				else:		# Otherwise, trim most common number of preceding bases
					if len(seq) > (pre_trim + post_trim):  # don't trim very short sequences
						start=pre_trim
					else: start=0
					
				pr_find = re.search(primer2,seq)
				if pr_find:		# If reverse primer is found, trim it and following bases
					end = pr_find.start()
				else:
					if len(seq) > (pre_trim + post_trim):
						end = len(seq) - post_trim - 1
					else: end = len(seq) - 1
				trimseq = seq[start:end]
				if trimseq[-1] != '\n': trimseq = trimseq + '\n'
				# write trimmed sequence and quality scores to the new file 
				out_fa.write(seqheader)
				out_fa.write(trimseq)
				
	out_fa.close()
	print ("Trimmed fasta records written to", out_name)					

if __name__=="__main__":
	parser = argparse.ArgumentParser(description=
	'Read a fasta file and write a new one with only bases inside the amplicon primers retained')
	parser.add_argument('-i', dest='in_fn', metavar='<infile>', 
					   help='the input fasta file of overlapped sequences. Default is seqs.fna')
	parser.add_argument('-o', dest='out_fn',  metavar='<outfile>',
					   help='filename for output fasta file. Default is trimmed<infile>')
	parser.add_argument('-d', dest='seqsRevDir',  metavar='<sequence_3_5>',
					   help='Is the sequence in \'Reverse\' direction, i.e. reverse primer first, or \'Forward\'? Default is Reverse')
	parser.add_argument('-v', dest='Vregion',  metavar='<V_region>',
					   help='V region. Default is V4: 533F to 805R. Other primers not yet available')

	args = parser.parse_args()
	if args.in_fn == None:
		args.in_fn = 'seqs.fna'
	if args.out_fn == None:
		args.out_fn = 'trimmed' + args.in_fn
		
	trim_amplicons(args.in_fn, args.out_fn, seqDirection=args.seqsRevDir, regionV=args.Vregion)
