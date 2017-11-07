## Jocelyn Sietsma Penington
## October 2017

## Use SNVerPool to find single nucleotide variations in the joint Pf bam files
## for each set, SNVerIndividual for the merged parent bam.
## Use bedtools subtract to extract events in daughter strains that are not in parents
## Use bedtools intersect, and grep, to find and count events in coding regions.

## On 9 Oct 2017 SNVer version is 0.5.3

PAPDIR=/wehisan/bioinf/bioinf-data/Papenfuss_lab/projects
BASEDIR=$PAPDIR/malaria/cowman_lab/drug_resistance
CURDIR=$BASEDIR/currentProject
ALIGNDIR=$CURDIR/alignment
VARDIR=$CURDIR/variants
REFDIR=$PAPDIR/reference_genomes/plasmodium/PlasmoDB-29_Pfalciparum3D7
SNVERD=/usr/local/bioinfsoftware/SNVer/current

moday="$(date +"%d%b")"
SNVDIR=$VARDIR/snverP${moday}
mkdir $SNVDIR

for s in A B D E
do
   java -jar $SNVERD/SNVerPool.jar -i $ALIGNDIR -c $ALIGNDIR/${s}pool.txt \
        -r $REFDIR/PlasmoDB-29_Pfalciparum3D7_Genome.fasta \
        -o $OUTDIR/$s
done


java -jar $SNVERD/SNVerIndividual.jar -i $ALIGNDIR/3D7-merge-B2_S1-F4_S4.bam  \
     -r $REFDIR/PlasmoDB-29_Pfalciparum3D7_Genome.fasta -n 1  \
     -o $OUTDIR/mergedParent

## Take output of SNVerP for drug-resistant and parent 3D7 strains
## compared to reference, find the difference between resistant and parent strains. 
## Based on annotateIntersect.sh in malaria/vaccine/output/SNVerI and 
## intersect_snver_vcf.sh in drug_resistance/analysis_tools/
## (SNVer commands were in separate script, now inserted above.)

## Latest gff file from PlasmoDB (v29) does not need cleaning for this application.

## Check date in name of SNVDIR !
export OUTDIR=$SNVDIR/removedParentVariants
mkdir $OUTDIR

cd $SNVDIR

# Make vcf files of SNVs that are in resistant strain and not in parents. 
for s in A B D E
do
## First copy header/comment lines
   grep '^#' ${s}.filter.vcf > $OUTDIR/${s}-parent.filter.vcf
      
   bedtools subtract -a ${s}.filter.vcf -b mergedParent.filter.vcf  \
      >> $OUTDIR/${s}-parent.filter.vcf
   
   # Count non-comment lines to get number of SNVs
   echo Number of SNVs in $s:
   grep -vc '^#' $OUTDIR/${s}-parent.filter.vcf 
   
	# Intersect with reference annotations in both orders, to get SNVs that lie in features, 
	# and features that contain SNVs.  For the SNV-first version use -u so that results 
	# aren't repeated for each description/feature in gff. This filtering is repeated 
	# in maldrugGilsonR/readSNVer.R
    # Count events with grep -c
   grep '^#' $OUTDIR/${s}-parent.filter.vcf > $OUTDIR/${s}_SNPgff.vcf
   bedtools intersect -u -a $OUTDIR/${s}-parent.filter.vcf \
      -b $REFDIR/PlasmoDB-29_Pfalciparum3D7.gff >> $OUTDIR/${s}_SNPgff.vcf
   echo Number of SNVs in $s that are in Features:
   grep -vc '^#' $OUTDIR/${s}_SNPgff.vcf
    
   bedtools intersect -u -b $OUTDIR/${s}-parent.filter.vcf \
      -a $REFDIR/PlasmoDB-29_Pfalciparum3D7.gff > $OUTDIR/gffFilter_${s}.gff
    echo Number of genes in $s with SNVs\; number of exons: 
    grep -c $'\tgene\t' $OUTDIR/gffFilter_${s}.gff 
    grep -c $'\texon\t' $OUTDIR/gffFilter_${s}.gff 
   
   ## Repeat process for indels ##
   grep '^#' ${s}.indel.filter.vcf > $OUTDIR/${s}-parent.indel.filter.vcf
   
   bedtools subtract -a ${s}.indel.filter.vcf  \
      -b mergedParent.indel.filter.vcf  \
      >> $OUTDIR/${s}-parent.indel.filter.vcf
    
   # Count non-comment lines to get number of SNVs
   echo Number of indels in $s:
   grep -vc '^#' $OUTDIR/${s}-parent.indel.filter.vcf 
   
	# Intersect with reference annotations in both orders.
    # Count events with grep -c
   grep '^#' $OUTDIR/${s}-parent.indel.filter.vcf > $OUTDIR/${s}_indelGff.vcf
   bedtools intersect -u -a $OUTDIR/${s}-parent.indel.filter.vcf \
      -b $REFDIR/PlasmoDB-29_Pfalciparum3D7.gff >> $OUTDIR/${s}_indelGff.vcf
   echo Number of indels in $s that are in Features:
   grep -vc '^#'  $OUTDIR/${s}_indelGff.vcf 
   
   bedtools intersect -u -b $OUTDIR/${s}-parent.indel.filter.vcf \
      -a $REFDIR/PlasmoDB-29_Pfalciparum3D7.gff > $OUTDIR/gffindel.filter_${s}.gff
    echo Number of genes in $s with indels\; number of exons: 
    grep -c $'\tgene\t' $OUTDIR/gffindel.filter_${s}.gff 
    grep -c $'\texon\t' $OUTDIR/gffindel.filter_${s}.gff 

done 

###############################################
