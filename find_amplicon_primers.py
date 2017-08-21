from __future__ import print_function
"""
find_amplicon_primers.py <inputfastqFilename> <primer> <rc_flag>

Author: J Sietsma Penington
Date: 7 Aug 2017

Based on code in trim_to_amplicon_primers.py
Count and report the position of a putative amplicon primer (or any short sequence) in a 
fastq file.
"""

import argparse
from collections import Counter
from Bio import SeqIO, SeqUtils
from Bio.Seq import Seq 
from Bio.Alphabet import IUPAC


def count_amplicons(in_name, fprimer, rc):
    Fprimer = Seq(fprimer, IUPAC.ambiguous_dna)
    pre_length = Counter()
    if rc:
        post_length = Counter()
        bothfound = 0
        Rprimer = Seq(fprimer, IUPAC.ambiguous_dna).reverse_complement()
        lenRprimer = len(Rprimer)   
    
    with open(in_name, 'r') as fastqF:
        for seqRecord in SeqIO.parse(fastqF, "fastq"):  
            Fpos = SeqUtils.nt_search(str(seqRecord.seq), str(Fprimer))            
            if len(Fpos) > 1:
                # SeqUtils.nt_search returns the pattern, followed by positions of any matches 
                # Forward primer found: increment pre_length
                pre_length[Fpos[1]] += 1
            if rc:
                RCpos = SeqUtils.nt_search(str(seqRecord.seq), str(Rprimer))
                if len(RCpos) > 1:
                    tail = len(seqRecord) - RCpos[-1] - lenRprimer 
                    post_length[tail] += 1
                    if len(Fpos) > 1:
                        bothfound += 1
    
    print ("Primers found:", sum(pre_length.values()))
    print ("Counts of pre_length:", pre_length)
    if rc:
        print ("Reverse primers found:", sum(post_length.values()))
        print("Counts of post_length", post_length)
        print ("Both primer and reverse_complement found:", bothfound)
        

if __name__=="__main__":
    parser = argparse.ArgumentParser(description=
    'Read a fastq file and count the occurences and positions of a short sequence ')
    parser.add_argument('in_fn', help='Input fastq file')
    parser.add_argument('primer', help='Ambiguous DNA code for amplicon')
    parser.add_argument('-rc', action = "store_true", 
                       help='Also search for reverse-complement of the sequence')

    args = parser.parse_args()
        
    count_amplicons(args.in_fn, args.primer, args.rc)
