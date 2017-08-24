from __future__ import print_function
"""
trim_to_amplicon_primers.py -i <inputfastqFilename> -o <outputfastqFilename>

Author: J Sietsma Penington
Date: 13 Nov 2014

Remove the 'universal primer' from each end of a fastq file of overlapped sequences 
assembled from paired-end reads. 
From the 5' end these sequences have 
possibly an 8-base barcode (which may have been removed, for example by QIIME's extract_barcodes.py),
followed by the 22-mer Illumina forward universal primer (which may not be complete),
then the V3 forward amplicon primer.
At the far (3') end there is the reverse-complement of the V4 reverse amplicon primer,
followed by the 21-mer reverse universal primer (which may not be complete),
optionally followed by an 8-base barcode
The program writes a fastq file with all bases before the V3 forward primer, and 
all bases after the V4 reverse primer, removed.
Primer may be present with errors. Rather than implementing an inexact search, if a primer 
is not found the program removes the most common length of sequence found before or after primers 

A sequence is written for each input sequence. Sequences deemed too short for trimming are 
written unchanged.

To Do:
1. Import proper sequence handling from Bio.seq OR write a proper reverse-complement function
2. Import an inexact search algorithm?
3. Allow input or reading of amplicon primers

---------------------------------------
Changes for ENDIA background QC project  7 Dec 2015
---------------------------------------
Region is V4, amplicon primers as given by Stephen Wilcox are:
805R GGACTACNVGGGTWTCTAAT
533F GTGYCAGCMGCCGCGGTAA
but direction is reversed: - reads contain 805R as supplied, followed further down by 
reverse-complement of 533F

"""

import sys
import argparse
import re
from collections import Counter

def rc(seq):
#this is lazy and ugly, but not worth fixing now
	complement= {'A':'T','T':'A','C':'G','G':'C', '[':']', ']':'['}
	c_seq = ''
	for b in seq:
		if b in complement: 
			b = complement[b]
		c_seq = c_seq + b
	return c_seq[::-1]

def trim_to_amplicons(in_name, out_name):
	f805R = re.compile("GGACTAC.[ACG]GGGT[AT]TCTAAT", re.I)
	f533F = "GTG[CT]CAGC[AC]GCCGCGGTAA"
	r533F = re.compile(rc(f533F), re.I)
	bothfound = 0
	pre_length = Counter(); pre_trim =0
	post_length = Counter(); post_trim = 0
	out_fq = open(out_name, 'w')
	
	with open(in_name, 'r') as fastq:
	# Read through file twice, first counting primer positions:
		for line in fastq:
			if line[0] =='@':
				seqheader = line
				seq = fastq.next()
				qualheader = fastq.next()
				qual = fastq.next()
				
				pf_find = re.search(f805R, seq)
				if pf_find: 
					# Forward primer found: increment pre_length
					pre_length[pf_find.start()] += 1
				pr_find = re.search(r533F,seq)
				if pr_find:
					# Reverse primer found: increment post_length
					tail = len(seq) - pr_find.end() - 1
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
		
	# Read through file again, trimming sequence and quality lines 
		fastq.seek(0)
		for line in fastq:
			if line[0] =='@':
				seqheader = line
				seq = fastq.next()
				qualheader = fastq.next()
				qual = fastq.next()
				
				pf_find = re.search(f805R, seq)
				if pf_find: 	# If forward primer is found, trim preceding bases
					start = pf_find.start()
				else:		# Otherwise, trim most common number of preceding bases
					if len(seq) > (pre_trim + post_trim):  # don't trim very short sequences
						start=pre_trim
					else: start=0
					
				pr_find = re.search(r533F,seq)
				if pr_find:		# If reverse primer is found, trim following bases
					end = pr_find.end()
				else:
					if len(seq) > (pre_trim + post_trim):
						end = len(seq) - post_trim - 1
					else: end = len(seq) - 1
				trimseq = seq[start:end]
				if trimseq[-1] != '\n': trimseq = trimseq + '\n'
				trimqual = qual[start:end]
				if trimqual[-1] != '\n': trimqual = trimqual + '\n'
				# write trimmed sequence and quality scores to the new file 
				out_fq.write(seqheader)
				out_fq.write(trimseq)
				out_fq.write(qualheader)
				out_fq.write(trimqual)
				
	out_fq.close()
	print ("Trimmed fastq records written to", out_name)					

if __name__=="__main__":
	parser = argparse.ArgumentParser(description=
	'Read a fastq file and write a new one with bases outside the amplicon primers removed')
	parser.add_argument('-i', dest='in_fn', metavar='<infile>', 
					   help='the input fastq file of overlapped sequences. Default is reads.fastq')
	parser.add_argument('-o', dest='out_fn',  metavar='<outfile>',
					   help='filename for output fastq file. Default is trimmed<infile>')

	args = parser.parse_args()
	if args.in_fn == None:
		args.in_fn = 'reads.fastq'
	if args.out_fn == None:
		args.out_fn = 'trimmed' + args.in_fn
		
	trim_to_amplicons(args.in_fn, args.out_fn)
