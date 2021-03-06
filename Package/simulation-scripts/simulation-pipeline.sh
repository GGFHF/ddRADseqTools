#!/bin/bash

#-------------------------------------------------------------------------------

# This software has been developed by Forest Genetics and Physiology Research Group,
# Technical University of Madrid (UPM)

# Licence: GNU General Public Licence Version 3

#-------------------------------------------------------------------------------

# This script generates fragments from a genome; it simulates a double digest
# and generates their PE reads; it quantifies and removes PCR duplicates and
# calculates the loci number without data; it demultiplexes the individuals; it
# trims the adapters and other Illumina specific sequences; it aligns the reads
# and gets SAM, BAM, BED and VCF format files to study and visualize alignments;
# and it gets the distinct loci list with mutations.

#-------------------------------------------------------------------------------

# WARNING

# This script uses the following bioinformatics tools:
#    - BWA v0.7.13 (http://bio-bwa.sourceforge.net/)
#    - SAMtools & BCFtools v0.1.19 (http://www.htslib.org/)
#    - BEDtools v2.25.0 (http://bedtools.readthedocs.io/)
#    - VCFtools v0.1.14-14 (https://vcftools.github.io/)

#-------------------------------------------------------------------------------

# Control parameters

if [ -n "$*" ]; then echo 'This script has not parameters.'; exit 1; fi

#-------------------------------------------------------------------------------

# Set run environment

BWADIR=$APPS/BWA                  # BWA 0.7.13 programs directory
SAMTOOLSDIR=$APPS/SAMtools-0      # SAMTools 0.0.19 programs directory
BCFTOOLSDIR=$APPS/BCFtools-0      # BCFTools 0.1.19 programs directory
BEDTOOLSDIR=$APPS/BEDtools/bin    # BEDTools 2.25.0 programs directory
#VCFTOOLSDIR=$APPS/VCFtools        # VCFTools 0.1.14-14 programs directory

DDRADSEQTOOLSDIR=$TRABAJO/ddRADseqTools       # ddRADseqTools programs directory
GENOMESDIR=$TRABAJO/ddRADseqTools/genomes     # genomes file directory
FRAGSDIR=$TRABAJO/ddRADseqTools/fragments     # fragments directory
READSDIR=$TRABAJO/ddRADseqTools/reads         # reads directory
STATSDIR=$TRABAJO/ddRADseqTools/statistics    # statistics directory
ALIGNDIR=$TRABAJO/ddRADseqTools/alignments    # alignments directory

if [ ! -d "$FRAGSDIR" ]; then mkdir $FRAGSDIR; fi
if [ ! -d "$READSDIR" ]; then mkdir $READSDIR; else rm -f $READSDIR/*; fi
if [ ! -d "$STATSDIR" ]; then mkdir $STATSDIR; else rm -f $STATSDIR/*; fi
if [ ! -d "$ALIGNDIR" ]; then mkdir $ALIGNDIR; else rm -f $ALIGNDIR/*; fi

GENOME_SCEREVISIAE=GCF_000146045.2_R64_genomic.fna.gz                         # Saccharomyces cerevisiae genome
GENOME_CELEGANS=GCF_000002985.6_WBcel235_genomic.fna.gz                       # Caenorhabditis elegans genome
GENOME_DMELANOGASTER=GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.fna.gz    # Drosophila melanogaster genome
GENOME_WAUROPUNCTATA=GCF_000956235.1_wasmannia.A_1.0_genomic.fna.gz           # Wasmannia auropunctata genome
GENOME_VBERUS=GCA_000800605.1_Vber.be_1.0_genomic.fna.gz                      # Vipera berus genome
GENOME_HSAPIENS=GCF_000001405.29_GRCh38.p3_genomic.fna.gz                     # Homo sapiens genome
GENOME_BNAPUS=GCF_000686985.1_Brassica_napus_assembly_v1.0_genomic.fna.gz     # Brassica napus genome
GENOME_QROBUR=ena.fasta                                                       # Quercus robur genome
GENOME_PTAEDA=ptaeda.v1.01.scaffolds.fasta.gz                                 # Pinus taeda genome
GENOME=$GENOME_SCEREVISIAE                                                    # genome used in this run

ENZYME1=EcoRI
ENZYME2=MseI
TECHNIQUE=IND1_IND2_DBR
FORMAT=FASTQ
READTYPE=PE
INDEX1LEN=6
INDEX2LEN=6
DBRLEN=4
WEND=end91
CEND=end92
INDIVIDUALSFILE=individuals-8index1-6index2.txt

if [ `ulimit -n` -lt 1024 ]; then ulimit -n 1024; fi

#-------------------------------------------------------------------------------

# Generate genome fragments and get statistics

echo '**************************************************'
echo 'GENOME FRAGMENTS ARE BEING GENERATED FROM GENOME ...'

$DDRADSEQTOOLSDIR/rsitesearch.py \
    --genfile=$GENOMESDIR/$GENOME \
    --fragsfile=$FRAGSDIR/fragments-genome.fasta \
    --rsfile=$DDRADSEQTOOLSDIR/restrictionsites.txt \
    --enzyme1=$ENZYME1 \
    --enzyme2=$ENZYME2 \
    --minfragsize=101 \
    --maxfragsize=300 \
    --fragstfile=$STATSDIR/fragments-genome-stats.txt \
    --fragstinterval=25 \
    --plot=YES \
    --verbose=YES \
    --trace=NO
if [ $? -ne 0 ]; then echo 'Script ended with errors.'; exit 1; fi

#-------------------------------------------------------------------------------

# Generate ddRADseq simulated reads

echo '**************************************************'
echo 'DDRADSEQ SIMULATED READS ARE BEING GENERATED ...'

$DDRADSEQTOOLSDIR/simddradseq.py \
    --fragsfile=$FRAGSDIR/fragments-genome.fasta \
    --technique=$TECHNIQUE \
    --format=$FORMAT \
    --readsfile=$READSDIR/reads \
    --readtype=$READTYPE \
    --rsfile=$DDRADSEQTOOLSDIR/restrictionsites.txt \
    --enzyme1=$ENZYME1 \
    --enzyme2=$ENZYME2 \
    --endsfile=$DDRADSEQTOOLSDIR/ends.txt \
    --index1len=$INDEX1LEN \
    --index2len=$INDEX2LEN \
    --dbrlen=$DBRLEN \
    --wend=$WEND \
    --cend=$CEND \
    --individualsfile=$DDRADSEQTOOLSDIR/$INDIVIDUALSFILE \
    --locinum=3000 \
    --readsnum=300000 \
    --minreadvar=0.8 \
    --maxreadvar=1.2 \
    --insertlen=100 \
    --mutprob=0.2 \
    --locusmaxmut=1 \
    --indelprob=0.1 \
    --maxindelsize=10 \
    --dropout=0.05 \
    --pcrdupprob=0.2 \
    --pcrdistribution=MULTINOMIAL \
    --multiparam=0.167,0.152,0.136,0.121,0.106,0.091,0.076,0.061,0.045,0.030,0.015 \
    --poissonparam=1.0 \
    --gcfactor=0.2 \
    --verbose=YES \
    --trace=NO
if [ $? -ne 0 ]; then echo 'Script ended with errors.'; exit 1; fi

#-------------------------------------------------------------------------------

# Remove the PCR duplicates

echo '**************************************************'
echo 'THE PRC DUPLICATES ARE BEING REMOVED ...'

$DDRADSEQTOOLSDIR/pcrdupremoval.py \
    --format=$FORMAT \
    --readtype=$READTYPE \
    --readsfile1=$READSDIR/reads-1.fastq \
    --readsfile2=$READSDIR/reads-2.fastq \
    --clearfile=$READSDIR/reads-cleared \
    --dupstfile=$STATSDIR/pcrduplicates-stats.txt \
    --plot=YES \
    --verbose=YES \
    --trace=NO
if [ $? -ne 0 ]; then echo 'Script ended with errors.'; exit 1; fi

#-------------------------------------------------------------------------------

# Demultiplex the individual files

echo '**************************************************'
echo 'INDIVIDUAL FILES ARE BEING DEMULTIPLEXED ...'

$DDRADSEQTOOLSDIR/indsdemultiplexing.py \
    --technique=$TECHNIQUE \
    --format=$FORMAT \
    --readtype=$READTYPE \
    --endsfile=$DDRADSEQTOOLSDIR/ends.txt \
    --index1len=$INDEX1LEN \
    --index2len=$INDEX2LEN \
    --dbrlen=$DBRLEN \
    --wend=$WEND \
    --cend=$CEND \
    --individualsfile=$DDRADSEQTOOLSDIR/$INDIVIDUALSFILE \
    --readsfile1=$READSDIR/reads-cleared-1.fastq \
    --readsfile2=$READSDIR/reads-cleared-2.fastq \
    --verbose=YES \
    --trace=NO
if [ $? -ne 0 ]; then echo 'Script ended with errors.'; exit 1; fi

#-------------------------------------------------------------------------------

# Trim adapters and other Illumina-specific sequences from reads

echo '**************************************************'
echo 'ADAPTERS AND OTHER ILLUMINA SPECIFIC SEQUENCES ARE BEING TRIMMED IN INDIVIDUAL FILES ...'

ls $READSDIR/demultiplexed-ind*-1.fastq > $READSDIR/reads-files.txt

while read FILE_1; do

    FILE_2=`echo $FILE_1 | sed 's/-1.fastq/-2.fastq/g'`
    FILE_TRIMMED=`echo $FILE_1 | sed 's/-1.fastq/-trimmed/g'`

    $DDRADSEQTOOLSDIR/readstrim.py \
        --technique=$TECHNIQUE \
        --format=$FORMAT \
        --readtype=$READTYPE \
        --endsfile=$DDRADSEQTOOLSDIR/ends.txt \
        --index1len=$INDEX1LEN \
        --index2len=$INDEX2LEN \
        --dbrlen=$DBRLEN \
        --wend=$WEND \
        --cend=$CEND \
        --readsfile1=$FILE_1 \
        --readsfile2=$FILE_2 \
        --trimfile=$FILE_TRIMMED \
        --verbose=YES \
        --trace=NO
    if [ $? -ne 0 ]; then echo 'Script ended with errors.'; exit 1; fi

done < $READSDIR/reads-files.txt

#-------------------------------------------------------------------------------

# Index the genome

echo '**************************************************'
echo 'GENOME IS BEING INDEXED ...'

$BWADIR/bwa index -a bwtsw $GENOMESDIR/$GENOME
if [ $? -ne 0 ]; then echo 'Script ended with errors.'; exit 1; fi


#-------------------------------------------------------------------------------

# Align sequences of individual files in SAM format

echo '**************************************************'
echo 'SEQUENCES OF INDIVIDUAL FILES ARE BEING ALIGNED IN SAM FORMAT ...'

ls $READSDIR/demultiplexed-ind*-trimmed-1.fastq > $READSDIR/reads-trimmed-files.txt

while read FILE1_TRIM; do

    FILE2_TRIM=`echo $FILE1_TRIM | sed 's/-trimmed-1.fastq/-trimmed-2.fastq/g'`
    FILE_SAM=`echo $FILE1_TRIM | sed 's/-1.fastq/.sam/g' | sed "s|$READSDIR|$ALIGNDIR|g"`
    $BWADIR/bwa mem $GENOMESDIR/$GENOME $FILE1_TRIM $FILE2_TRIM > $FILE_SAM
    if [ $? -ne 0 ]; then echo 'Script ended with errors.'; exit 1; fi

done < $READSDIR/reads-trimmed-files.txt

#-------------------------------------------------------------------------------

# Convert SAM files to BED, BAM and VCF format

echo '**************************************************'
echo 'SAM FILES ARE BEING CONVERTED IN BAM, BED AND VCF FORMAT ...'

ls $ALIGNDIR/*.sam > $ALIGNDIR/sam-files.txt

while read FILE_SAM; do

    FILE_BAM=`echo $FILE_SAM | sed 's/.sam/.bam/g'`
    FILE_BAM_STATS=`echo $FILE_SAM | sed 's/.sam/-stats-bam.txt/g'`
    FILE_BED=`echo $FILE_SAM | sed 's/.sam/.bed/g'`
    FILE_SORTED=`echo $FILE_SAM | sed 's/.sam/.sorted/g'`
    FILE_SORTED_BAM=`echo $FILE_SAM | sed 's/.sam/.sorted.bam/g'`
    FILE_VCF=`echo $FILE_SAM | sed 's/.sam/.vcf/g'`
    FILE_VCF_STATS1=`echo $FILE_SAM | sed 's/.sam/-stats-vcf-1.txt/g'`
    FILE_VCF_STATS2=`echo $FILE_SAM | sed 's/.sam/-stats-vcf-2.txt/g'`

    # convert SAM file to BAM format
    $SAMTOOLSDIR/samtools  view -Sb -o $FILE_BAM $FILE_SAM
    if [ $? -ne 0 ]; then echo 'Script ended with errors.'; exit 1; fi
    
    # get the aligment statistic in a file
    $SAMTOOLSDIR/samtools flagstat $FILE_BAM >$FILE_BAM_STATS
    if [ $? -ne 0 ]; then echo 'Script ended with errors.'; exit 1; fi

    # convert the BAM file to BED format
    $BEDTOOLSDIR/bedtools bamtobed -i $FILE_BAM > $FILE_BED
    if [ $? -ne 0 ]; then echo 'Script ended with errors.'; exit 1; fi
    
    # sort the BAM file
    $SAMTOOLSDIR/samtools sort $FILE_BAM $FILE_SORTED
    if [ $? -ne 0 ]; then echo 'Script ended with errors.'; exit 1; fi
    
    # index the BAM file
    $SAMTOOLSDIR/samtools index $FILE_SORTED_BAM
    if [ $? -ne 0 ]; then echo 'Script ended with errors.'; exit 1; fi

    # convert the BAM file to VCF format
    $SAMTOOLSDIR/samtools mpileup -uf  $GENOMESDIR/$GENOME $FILE_SORTED_BAM | $BCFTOOLSDIR/bcftools view -vcg - > $FILE_VCF
    if [ $? -ne 0 ]; then echo 'Script ended with errors.'; exit 1; fi

    # get the variant statistic in a file with SAMtools
    $BCFTOOLSDIR/vcfutils.pl qstats $FILE_VCF > $FILE_VCF_STATS1
    if [ $? -ne 0 ]; then echo 'Script ended with errors.'; exit 1; fi

    # get the variant statistic in a file with VCFtools
    vcftools --vcf $FILE_VCF > $FILE_VCF_STATS2
    if [ $? -ne 0 ]; then echo 'Script ended with errors.'; exit 1; fi

done < $ALIGNDIR/sam-files.txt

#-------------------------------------------------------------------------------

# Get the distinct loci list with mutations

echo '**************************************************'
echo 'THE DISTINCT LOCI LIST WITH MUTATIONS ARE BEING GOT ...'

LOCI_FILE=$ALIGNDIR/loci.txt
cat >$LOCI_FILE <<EOF
EOF
DISTINCT_LOCI_FILE=$ALIGNDIR/distinct-loci.txt
DISTINCT_LOCI_TOTAL_FILE=$ALIGNDIR/distinct-loci-total.txt

ls $ALIGNDIR/*.vcf > $ALIGNDIR/vcf-files.txt

while read FILE_VCF; do

    while read -r LINE; do

        if [[ $LINE =~ ^#.*$ ]]; then continue; fi

        echo "$LINE" | cut -f 1,2 >> $LOCI_FILE

    done < $FILE_VCF

done < $ALIGNDIR/vcf-files.txt

sort -u $LOCI_FILE > $DISTINCT_LOCI_FILE

DISTINCT_LOCI_TOTAL=`wc -l $DISTINCT_LOCI_FILE | cut -f1 -d' '`
echo "The dictinct loci list contains $DISTINCT_LOCI_TOTAL loci."
echo "$DISTINCT_LOCI_TOTAL" > $DISTINCT_LOCI_TOTAL_FILE

#-------------------------------------------------------------------------------

# End
echo '**************************************************'
exit 0

#-------------------------------------------------------------------------------
