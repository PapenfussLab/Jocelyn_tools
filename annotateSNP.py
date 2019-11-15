"""
annotate_SNP.py <inputGFFfile> <inputREFERENCEfastafile> <inputVCFfile>
For each SNV in vcf file, extracts the codon context from the fasta file and calculates
whether amino acid product is changed.
Outputs a list of annotations to stdout of either Synonymous or the aa change
Original author: Jan Schroeder
Modified: Jocelyn Sietsma Penington March 2016 to use on malaria vaccine project
Modified Nov 2019 to add position of AA in transcript, and correct omission of first 
'frame' bases of each CDS

My vcf files do not have strand information, my GFF does.
Original code assumes vcf reads have a start and end. I am not going to read indel files, 
as they are unable to be synonomous if in a coding region, so all changes are 1-1
Where there is more than one alternate allele in the vcf, all changes in coding sequences 
are reported. 

Features of type 'CDS' are used, rather than exons, because these include the frame.
Most exons have a matching CDS feature, but a few don't.
Now that I have changed to transcripts, 'frame' is not used and this could be changed.
Changes are now computed in transcript, and AA position in transcript is reported.

BUG: Each SNV is checked independently. If 2 fall in the same codon, it will be
reported as 2 (incorrect) AA changes from the reference
BUG: Code will only report one transcript per variant: last CDS in gff that contains 
variant will determine which transcript is used.
 
"""

#!/usr/bin/env python
from __future__ import print_function 
import sys
import HTSeq
'''
wehisan in July 2019 has 
	HTSeq.__version__ '0.7.2' if loaded python is 3.5.1
	HTSeq.__version__ '0.6.1p1' if loaded python is 2.7
'''
import itertools
import string
import math

def rc(dna):
 complements = str.maketrans('acgtrymkbdhvACGTRYMKBDHV', 'tgcayrkmvhdbTGCAYRKMVHDB')
 rcseq = dna.translate(complements)[::-1]
 return rcseq

amino = {"TTT":"F|Phe","TTC":"F|Phe","TTA":"L|Leu","TTG":"L|Leu","TCT":"S|Ser","TCC":"S|Ser",
  "TCA":"S|Ser","TCG":"S|Ser", "TAT":"Y|Tyr","TAC":"Y|Tyr","TAA":"*|Stp","TAG":"*|Stp",
  "TGT":"C|Cys","TGC":"C|Cys","TGA":"*|Stp","TGG":"W|Trp", "CTT":"L|Leu","CTC":"L|Leu",
  "CTA":"L|Leu","CTG":"L|Leu","CCT":"P|Pro","CCC":"P|Pro","CCA":"P|Pro","CCG":"P|Pro",
  "CAT":"H|His","CAC":"H|His","CAA":"Q|Gln","CAG":"Q|Gln","CGT":"R|Arg","CGC":"R|Arg",
  "CGA":"R|Arg","CGG":"R|Arg", "ATT":"I|Ile","ATC":"I|Ile","ATA":"I|Ile","ATG":"M|Met",
  "ACT":"T|Thr","ACC":"T|Thr","ACA":"T|Thr","ACG":"T|Thr", "AAT":"N|Asn","AAC":"N|Asn",
  "AAA":"K|Lys","AAG":"K|Lys","AGT":"S|Ser","AGC":"S|Ser","AGA":"R|Arg","AGG":"R|Arg",
  "GTT":"V|Val","GTC":"V|Val","GTA":"V|Val","GTG":"V|Val","GCT":"A|Ala","GCC":"A|Ala",
  "GCA":"A|Ala","GCG":"A|Ala", "GAT":"D|Asp","GAC":"D|Asp","GAA":"E|Glu",
  "GAG":"E|Glu","GGT":"G|Gly","GGC":"G|Gly","GGA":"G|Gly","GGG":"G|Gly"}

def amino_change(aa1, aa2):
  if len(aa1) != len(aa2): 
    print ("Error in amino_change: different length!")
    return
  change = []
  for i in range(len(aa1)):
    if aa1[i] != aa2[i]:
      change = [aa1[i], aa2[i]]
      return change

def protein(dna, frame):
  aaseq = ""
  x = frame
  while x < len(dna) - 2:
    try :
      aaseq += amino[dna[x:x+3]].split("|")[0]
    except (KeyError):
      aaseq += "X"
    x += 3
  return aaseq

sequences = dict( (s.name, s) for s in HTSeq.FastaReader(sys.argv[2]) )

gff_file = HTSeq.GFF_Reader(sys.argv[1], end_included=True)

CDSfeat = HTSeq.GenomicArray( "auto", stranded=False, typecode='O' )
# Strand information is available in array; 'False' means not required for indexing

transcript = {}

for feature in gff_file:
	if feature.type == "transcript":
		transcript[ feature.name ] = { 'iv' : feature.iv,  # .iv is GenomicInterval
		'CDSfeats' : [ ] }
	if feature.type == "CDS": 
		transcript[ feature.attr[ "Parent" ] ][ 'CDSfeats' ].append( feature.iv )
		## Future worry: do I need CDS.frame in transcript object?
		CDSfeat[ feature.iv ] = feature  
    
print ( "# Chrom\tPos\tPos in CDS\tBase change\tAA change\tAA pos in transcript\ttranscript ID" )
vcfr = HTSeq.VCF_Reader( sys.argv[3])

for vc in vcfr:
	vCDS = CDSfeat[ vc.pos ] 
	# vCDS.iv.start is base before 1st base of CDS
	if not vCDS==None and not vc.pos.start==vCDS.iv.start: 
		vTranscript = transcript[ vCDS.attr[ "Parent" ] ]
		refseq = str( HTSeq.Sequence( 
				 sequences[vCDS.iv.chrom].seq[vCDS.iv.start:vCDS.iv.end] ) )
		refseqT = ''.join( str( HTSeq.Sequence( 
					sequences[CDS.chrom].seq[CDS.start:CDS.end] ) )
							for CDS in vTranscript['CDSfeats'] )
		relpos = vc.pos.start - vCDS.iv.start 
		# if variant is 1st base of CDS, relpos=1
		if refseq[relpos-1] != vc.ref:
			print ("ERROR: Reference Base not according to SNP")
		for alt in vc.alt:        # vc.alt is a list of alternative base(s)
			alternateSeq = refseq[0:relpos-1] + alt + refseq[relpos:]  
			altSeqT = ''
			relposT = 0
			exonFound = False
			transLen = 0
			for exon in vTranscript['CDSfeats']: 
				transLen += exon.length
				if exon.contains( vc.pos ):
					altSeqT = ''.join( ( altSeqT, alternateSeq ) )
					relposT += vc.pos.start - exon.start 
					exonFound = True
				else:
					altSeqT = ''.join( ( altSeqT, str( HTSeq.Sequence( 
						sequences[exon.chrom].seq[exon.start:exon.end] ) )
						) )
					if not exonFound:
						relposT += exon.length
			if vTranscript['iv'].strand == "-":
				refseqT = rc(refseqT)
				altSeqT = rc(altSeqT)
				relpos = vCDS.iv.end + 1 - vc.pos.start 
				relposT = transLen + 1 - relposT
			aapos = int( math.ceil( relposT /3.0 ) )
			aa = protein(refseqT, 0)
			aa2 = protein(altSeqT, 0)
			if aa != aa2:
				change = amino_change(aa,aa2)
				ch="%s->%s" %( change[0], change[1])
			else:
				ch = 'Synonymous'
			snp = "%s->%s" %( vc.ref, alt)
			list = [vc.chrom, str(vc.pos.start), str(relpos), snp, ch, str(aapos), 
					vCDS.attr[ "Parent" ] ]
			print ( '\t'.join(list) )
	else:
		ch = 'X'
		gene_name= 'X'
		snp = "%s->%s" %( vc.ref, vc.alt[0]) # when not in exon, only 1st alternate is listed
		list = [vc.chrom, str(vc.pos.start), str(-1), snp, ch, str(-1), gene_name]
		print ( '\t'.join(list) )
