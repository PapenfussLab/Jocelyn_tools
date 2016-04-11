"""
annotate_SNP.py <inputGFFfile> <inputREFERENCEfastafile> <inputVCFfile>
For each SNV in vcf file, extracts the codon context from the fasta file and calculates
whether amino acid product is changed.
Outputs a list of annotations to stdout of either Synonymous or the aa change
Author: Jan Schroeder
Modified: Jocelyn Sietsma Penington March 2016 to use on malaria vaccine project

My vcf files do not have strand information, my GFF does
Original code assumes vcf reads have a start and end. I am not going to read indel files, 
as they are unable to be synonomous if in a coding region, so all changes are 1-1
Where there is more than one alternate allele in the vcf, all changes in coding sequences 
are reported. (Loop is inside test for exon)

'exons' extraction actually extracts features of type 'CDS' because these include the 
frame.
"""

#!/usr/bin/env python

import sys
import HTSeq
import itertools
import string

def rc(dna):
 complements = string.maketrans('acgtrymkbdhvACGTRYMKBDHV', 'tgcayrkmvhdbTGCAYRKMVHDB')
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
    print "Error in amino_change: different length!"
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

exons = HTSeq.GenomicArray( "auto", stranded=False, typecode='O' )
# Strand information is available in array; 'False' means not required for indexing

for feature in gff_file:
  if feature.type == "CDS": 
    exons[ feature.iv ] = feature
    
print "# Chrom\tPos\tPos in CDS\tBase change\tAA change\tgene ID"
vcfr = HTSeq.VCF_Reader( sys.argv[3])

for vc in vcfr:
	l = exons[ vc.pos ] 
	# l.iv.start is base before 1st base of CDS
	if not l==None and not vc.pos.start==l.iv.start: 
		refseq = sequences[l.iv.chrom].seq[l.iv.start:l.iv.end]
		if len(refseq) != l.iv.length:
			print "Error in sequence lengths"
		relpos = vc.pos.start - l.iv.start 
		# if variant is 1st base of exon, relpos=1
		if refseq[relpos-1] != vc.ref:
			print "Error: Reference Base not according to SNP!"
		for alt in vc.alt:        # vc.alt is a list of alternative base(s)
			alternateSeq = refseq[0:relpos-1] + alt + refseq[relpos:]    
			if len(refseq) != len(alternateSeq):
				print 'ALARM - length error when substituting alt for ref'
				print vc, relpos, l
			if l.iv.strand == "-":
				refseq = rc(refseq)
				alternateSeq = rc(alternateSeq)
			aa = protein(refseq, l.frame)
			aa2 = protein(alternateSeq, l.frame)
			if aa != aa2:
				change = amino_change(aa,aa2)
				ch="%s->%s" %( change[0], change[1])
			else:
				ch = 'Synonymous'
			snp = "%s->%s" %( vc.ref, alt)
			list = [vc.chrom, str(vc.pos.start), str(relpos), snp, ch, l.name]
			print '\t'.join(list)
	else:
		ch = 'X'
		gene_name= 'X'
		snp = "%s->%s" %( vc.ref, vc.alt[0]) # when not in exon, only 1st alternate is listed
		list = [vc.chrom, str(vc.pos.start), str(-1), snp, ch, gene_name]
		print '\t'.join(list)
