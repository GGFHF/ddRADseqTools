#!/bin/bash

#-------------------------------------------------------------------------------

# This script generates fragments from the Saccharomyces cerevisiae genome; it
# simulates a double digest and generates their PE reads; it quantifies and
# removes PCR duplicates and calculates the loci number without data; it
# calculates mutation statistics; it demultiplexes the individuals; it trims the
# adapters and other Illumina specific sequences; it aligns the reads and gets
# SAM, BAM, BED and VCF format files to study and visualize alignments; and it
# gets the distinct loci list with mutations.

#-------------------------------------------------------------------------------

# Control parameters

if [ -n "$*" ]; then echo 'This script has not parameters.'; exit 1; fi

#-------------------------------------------------------------------------------

# Set run environment

SAMTOOLSDIR=/usr/share/samtools    # SAMTools directory

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

GENOME=GCF_000146045.2_R64_genomic.fna.gz    # Saccharomyces cerevisiae genome
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

MUTPROB=(0.001 0.010 0.020 0.030 0.040 0.050 0.060 0.070 0.080 0.090 0.100)

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

for I in "${MUTPROB[@]}"
do

    # Generate ddRADseq simulated reads

    echo '**************************************************'
    echo "MUTPROB=$I"

    echo '--------------------------------------------------'
    echo 'DDRADSEQ SIMULATED READS ARE BEING GENERATED ...'

    $DDRADSEQTOOLSDIR/simddradseq.py \
        --fragsfile=$FRAGSDIR/fragments-genome.fasta \
        --technique=$TECHNIQUE \
        --format=$FORMAT \
        --readsfile=$READSDIR/'reads-'$I \
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
        --mutprob=$I \
        --locusmaxmut=1 \
        --indelprob=0.1 \
        --maxindelsize=10 \
        --dropout=0.0 \
        --pcrdupprob=0.2 \
        --pcrdistribution=MULTINOMIAL \
        --multiparam=0.167,0.152,0.136,0.121,0.106,0.091,0.076,0.061,0.045,0.030,0.015 \
        --poissonparam=1.0 \
        --gcfactor=0.2 \
        --trace=NO \
        --verbose=YES \
        --trace=NO
    if [ $? -ne 0 ]; then echo 'Script ended with errors.'; exit 1; fi

    # Remove the PCR duplicates

    echo '--------------------------------------------------'
    echo 'THE PRC DUPLICATES ARE BEING REMOVED ...'

    $DDRADSEQTOOLSDIR/pcrdupremoval.py \
        --format=$FORMAT \
        --readtype=$READTYPE \
        --readsfile1=$READSDIR/'reads-'$I'-1.fastq' \
        --readsfile2=$READSDIR/'reads-'$I'-2.fastq' \
        --clearfile=$READSDIR/'reads-cleared-'$I \
        --dupstfile=$STATSDIR/'pcrduplicates-stats-'$I'.txt' \
        --plot=YES \
        --verbose=YES \
        --trace=NO
    if [ $? -ne 0 ]; then echo 'Script ended with errors.'; exit 1; fi

    # Demultiplex the individual files

    echo '--------------------------------------------------'
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
        --readsfile1=$READSDIR/'reads-cleared-'$I'-1.fastq' \
        --readsfile2=$READSDIR/'reads-cleared-'$I'-2.fastq' \
        --verbose=YES \
        --trace=NO
    if [ $? -ne 0 ]; then echo 'Script ended with errors.'; exit 1; fi

    ls $READSDIR/demultiplexed-*-1.fastq > $READSDIR/reads-files.txt

    while read FILE_1; do

        FILE_1_R=`echo $FILE_1 | sed 's/demultiplexed/inddemultiplexed-'$I'/g'`
        FILE_2=`echo $FILE_1 | sed 's/-1.fastq/-2.fastq/g'`
        FILE_2_R=`echo $FILE_2 | sed 's/demultiplexed/inddemultiplexed-'$I'/g'`

        mv $FILE_1 $FILE_1_R
        mv $FILE_2 $FILE_2_R

    done < $READSDIR/reads-files.txt

    # Calcultes a statistics of mutated reads and no mutared reads of each individual

    echo '--------------------------------------------------'
    echo 'MUTATION STATISTICS ARE BEING CALCULATED ...'

    FILE_STATS=$STATSDIR/'mutations-stats-'$I'.csv'
cat >$FILE_STATS <<EOF
individual;no mutated reads;mutated reads;
EOF

    ls `echo $READSDIR/'inddemultiplexed-'$I'-*-1.fastq'` > $READSDIR/reads-files.txt

    while read FILE_1; do

        if [[ $FILE_1 =~ .*errors.* ]]; then continue; fi

        gawk 'function ltrim(s) { sub(/^[ \t\r\n]+/, "", s); return s }

              function rtrim(s) { sub(/[ \t\r\n]+$/, "", s); return s }

              function trim(s)  { return rtrim(ltrim(s)) }

              BEGIN { FS = "|"; mutated_reads = 0; no_mutated_reads = 0 }

              /@/ { if (FNR = 1) individual = trim(substr($6, 14, length($6))) }

              /mutated: True/ { ++mutated_reads }

              /mutated: False/ { ++no_mutated_reads }

              END { printf "%s;%d;%d;\n", individual, no_mutated_reads, mutated_reads }

             ' $FILE_1 >> $FILE_STATS
        if [ $? -ne 0 ]; then echo 'Script ended with errors.'; exit 1; fi

    done < $READSDIR/reads-files.txt

    # Trim adapters and other Illumina-specific sequences from reads

    echo '--------------------------------------------------'
    echo 'ADAPTERS AND OTHER ILLUMINA SPECIFIC SEQUENCES ARE BEING TRIMMED IN INDIVIDUAL FILES ...'

    ls `echo $READSDIR/'inddemultiplexed-'$I'-*-1.fastq'` > $READSDIR/reads-files.txt

    while read FILE_1; do

        if [[ $FILE_1 =~ .*errors.* ]]; then continue; fi

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

    # Align sequences of individual files in SAM format

    echo '--------------------------------------------------'
    echo 'SEQUENCES OF INDIVIDUAL FILES ARE BEING ALIGNED IN SAM FORMAT ...'

    ls `echo $READSDIR/'inddemultiplexed-'$I'-*-trimmed-1.fastq'` > $READSDIR/reads-trimmed-files.txt

    while read FILE1_TRIM; do

        FILE2_TRIM=`echo $FILE1_TRIM | sed 's/-trimmed-1.fastq/-trimmed-2.fastq/g'`
        FILE_SAM=`echo $FILE1_TRIM | sed 's/-1.fastq/.sam/g' | sed "s|$READSDIR|$ALIGNDIR|g"`
        bwa mem $GENOMESDIR/$GENOME $FILE1_TRIM $FILE2_TRIM > $FILE_SAM
        if [ $? -ne 0 ]; then echo 'Script ended with errors.'; exit 1; fi

    done < $READSDIR/reads-trimmed-files.txt

    #-------------------------------------------------------------------------------

    # Convert SAM files to BED, BAM and VCF format

    echo '--------------------------------------------------'
    echo 'SAM FILES ARE BEING CONVERTED IN BAM, BED AND VCF FORMAT ...'

    ls `echo $ALIGNDIR/'inddemultiplexed-'$I'-*.sam'` > $ALIGNDIR/sam-files.txt

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
        samtools view -bS $FILE_SAM >$FILE_BAM
        if [ $? -ne 0 ]; then echo 'Script ended with errors.'; exit 1; fi
        
        # get aligment statistics in a file
        samtools flagstat $FILE_BAM >$FILE_BAM_STATS
        if [ $? -ne 0 ]; then echo 'Script ended with errors.'; exit 1; fi

        # convert BAM file to BED format
        bedtools bamtobed -i $FILE_BAM > $FILE_BED
        if [ $? -ne 0 ]; then echo 'Script ended with errors.'; exit 1; fi
        
        # sort the BAM file
        samtools sort $FILE_BAM $FILE_SORTED
        if [ $? -ne 0 ]; then echo 'Script ended with errors.'; exit 1; fi
        
        # index the BAM file
        samtools index $FILE_SORTED_BAM
        if [ $? -ne 0 ]; then echo 'Script ended with errors.'; exit 1; fi

        # convert BAM file to VCF format
        samtools mpileup -uf  $GENOMESDIR/$GENOME $FILE_SORTED_BAM | bcftools view -vcg - > $FILE_VCF
        if [ $? -ne 0 ]; then echo 'Script ended with errors.'; exit 1; fi

        # get variant statistics in a file with SAMtools
        $SAMTOOLSDIR/vcfutils.pl qstats $FILE_VCF > $FILE_VCF_STATS1
        # -- if [ $? -ne 0 ]; then echo 'Script ended with errors.'; exit 1; fi

        # get variant statistics in a file with VCFtools
        vcftools --vcf $FILE_VCF > $FILE_VCF_STATS2
        if [ $? -ne 0 ]; then echo 'Script ended with errors.'; exit 1; fi

    done < $ALIGNDIR/sam-files.txt

    # Get the distinct loci list with mutations

    echo '--------------------------------------------------'
    echo 'THE DISTINCT LOCI LIST WITH MUTATIONS ARE BEING GOT ...'

    LOCI_FILE=$ALIGNDIR/loci-$I.txt
cat >$LOCI_FILE <<EOF
EOF
    DISTINCT_LOCI_FILE=$ALIGNDIR/distinct-loci-$I.txt
    DISTINCT_LOCI_TOTAL_FILE=$ALIGNDIR/distinct-loci-total-$I.txt

    ls `echo $ALIGNDIR/'inddemultiplexed-'$I'-*.vcf'` > $ALIGNDIR/vcf-files.txt

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

done

#-------------------------------------------------------------------------------

# End
echo '**************************************************'
exit 0

#-------------------------------------------------------------------------------
