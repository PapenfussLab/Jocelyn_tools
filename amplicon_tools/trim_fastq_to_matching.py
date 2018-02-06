"""
trim_fastq_to_matching.py -f <inputfastqFilename> -m <fastqFiletoMatch> -o <outputfastqFilename>

Author: J Sietsma Penington
Date: 12 April 2016

Removing all records in fastq sequence file that are not present in matching fastq with 
same sequence ID. 
If the text after whitespace in the seq header differs, that from the match file is written.
The purpose is to take a pair of sequence and barcode files, and create a new pair with 
identical header lines, suitable for input to QIIME's split_libraries_fastq.py
"""

import sys
import argparse

def matching_headers(seq_name, match_name, out_name):
	matchcount = 0
	IDdict = {}
	seqlist = []
	with open(match_name, 'r') as matchfq:
		for line in matchfq:
			IDdict[ line.split()[0] ] = [line, matchfq.next(), matchfq.next(),matchfq.next()]
	with open(out_name, 'w') as out_fq:
		with open( match_name + 'trimmed', 'w') as out_match:
			with open(seq_name, 'r') as fastq:
				for line in fastq:
					seqID = line.split()[0]
					seq = fastq.next()
					qualheader = fastq.next()
					qual = fastq.next()
				
					if seqID in IDdict:
						matchcount+=1
						seqlist.append(seqID)
			# write header from match-list, sequence and quality scores from seq, to the new file 
						out_fq.write(IDdict[seqID][0])
						out_fq.write(seq)
						out_fq.write(qualheader)
						out_fq.write(qual)
			# write the match-file block of 4 lines to new match-file
						for i in range(4):
							out_match.write(IDdict[seqID][i])
			print "%i fastq records in %s that match sequence IDs in %s written to %s" \
				% (matchcount, seq_name, match_name, out_name)	
			print "%i fastq records from %s that are also in %s written to %s" \
				% (matchcount, match_name, seq_name, match_name + 'trimmed')


if __name__=="__main__":
	parser = argparse.ArgumentParser(description=
	("Read 2 fastq files and write a new one containing sequences from the first file "
	"which have sequence IDs in the second. \n Sequence header lines are from the 2nd file") )
	parser.add_argument('-f', dest='seq_fn', metavar='<seqfile>', 
						help=('the input fastq file of sequences. '
						'A subset will be written to the output file'),
						default='reads.fastq')
	parser.add_argument('-m', dest='match_fn', metavar='<matchfile>', required=True, 
						help=('the input fastq file of barcodes or other '
						'list of sequences to match. Sequences with matching IDs '
						'will be written to the output file'))
	parser.add_argument('-o', dest='out_fn',  metavar='<outfile>',
						help='filename for output fastq file. Default is trimmed<infile>')

	args = parser.parse_args()
	if args.out_fn == None:
		args.out_fn = args.seq_fn + 'trimmed'
# 	Note this is a misuse of the argparse package, and by convention these arguments 
# 	should be positional, as -option format is for options, not standard inputs.
# 	However it is good learning to use argparse
	matching_headers(args.seq_fn, args.match_fn, args.out_fn)

