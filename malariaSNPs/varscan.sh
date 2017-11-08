## VarScan2 with 'somatic' function to find differences between resistant strains 
## and parent-line clones.
## Only VarScan native format output for now. 
## Jocelyn Sietsma Penington, October 2017

## File-name format of sample bam files is <strainID>_S<sampleID>_nodup.bam

PAPDIR=/wehisan/bioinf/bioinf-data/Papenfuss_lab/projects
BASEDIR=$PAPDIR/malaria/cowman_lab/drug_resistance
CURDIR=$BASEDIR/currentProject
ALIGNDIR=$CURDIR/alignment
VARDIR=$CURDIR/variants
REFDIR=$PAPDIR/reference_genomes/plasmodium/PlasmoDB-29_Pfalciparum3D7

moday="$(date +"%d%b")"
VS2DIR=$VARDIR/varscan${moday}
mkdir $VS2DIR

## For each drug-resistant strain, make a mpileup file with combined parent strain,
## then run VarScan2. Parent first means it is read1='normal'

## Running with joint pileups of all samples for each strain, and a separate parent pileup
       samtools mpileup -f $REFDIR/PlasmoDB-29_Pfalciparum3D7_Genome.fasta  \
          -q 15 -B $ALIGNDIR/3D7-merge-B2_S1-F4_S4.bam   \
          --output $ALIGNDIR/mergeS1toS4.pileup

for s in A B D E
do
   samtools mpileup -f $REFDIR/PlasmoDB-29_Pfalciparum3D7_Genome.fasta  \
      -q 15 -B $ALIGNDIR/${s}_S*_nodup.bam  \
      --output $ALIGNDIR/${s}.mpileup

    VarScan somatic $ALIGNDIR/mergeS1toS4.pileup  \
      $ALIGNDIR/${s}.mpileup  $VS2DIR/mergeS1toS4_${s}  \
      --strandfilter 0 --min-coverage 10  --min-var-freq 0.2
done
