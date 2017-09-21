#!/bin/bash

#-------------------------------------------------------------------------------

# This script generates fragments from a genome, simulates a double digest,
# generates their reads, removes PCR duplicates, demultiplexes the individuals
# calcultes a statistics of mutated reads and no mutared reads of each individual

#-------------------------------------------------------------------------------

# Control parameters

if [ -n "$*" ]; then echo 'This script has not parameters.'; exit 1; fi

#-------------------------------------------------------------------------------

# Set run environment

SAMTOOLSDIR=/usr/share/samtools    # SAMTools directory

DDRADSEQTOOLSDIR=$TRABAJO/ddRADseqTools/Exe    # ddRADseqTools programs directory
FRAGSDIR=$TRABAJO/ddRADseqTools/Fragments      # fragments directory
READSDIR=$TRABAJO/ddRADseqTools/Reads          # reads directory
STATSDIR=$TRABAJO/ddRADseqTools/Statistics     # statistics directory
ALIGNDIR=$TRABAJO/ddRADseqTools/Alignments     # alignments directory

if [ ! -d "$FRAGSDIR" ]; then mkdir $FRAGSDIR; fi
if [ ! -d "$READSDIR" ]; then mkdir $READSDIR; else rm -f $READSDIR/*; fi
if [ ! -d "$STATSDIR" ]; then mkdir $STATSDIR; else rm -f $STATSDIR/*; fi
if [ ! -d "$ALIGNDIR" ]; then mkdir $ALIGNDIR; else rm -f $ALIGNDIR/*; fi

SCEREVISIAE='Scerevisiae'                                           # Saccharomyces cerevisiae
SCEREVISIAE_FRAGS_FILE=$SCEREVISIAE'-fragments-EcoRI-MseI.fasta'    # Saccharomyces cerevisiae framents gotten by EcoI and MseI digest
MUTPROB=(0.0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9)

ENZYME1=EcoRI
ENZYME2=MseI
TECHNIQUE=IND1_IND2_DBR
FORMAT=FASTQ
READTYPE=PE
INDEX1LEN=6
INDEX2LEN=6
DBRLEN=4
WEND=end31
CEND=end32
INDIVIDUALSFILE=individuals-8index1-6index2.txt

if [ `ulimit -n` -lt 1024 ]; then ulimit -n 1024; fi

#-------------------------------------------------------------------------------

echo '**************************************************'
echo 'GENERATION OF READS AND STATISTICS'

for I in "${MUTPROB[@]}"
do

    # Generate ddRADseq simulated reads

    echo '--------------------------------------------------'
    echo "DDRADSEQ SIMULATED READS ARE BEING GENERATED WITH MUTPROB=$I ..."

    $DDRADSEQTOOLSDIR/simddradseq.py \
        --fragsfile=$FRAGSDIR/$SCEREVISIAE_FRAGS_FILE \
        --technique=$TECHNIQUE \
        --format=$FORMAT \
        --readsfile=$READSDIR/$SCEREVISIAE'-reads-'$I \
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
        --verbose=YES \
        --trace=NO
    if [ $? -ne 0 ]; then echo 'Script ended with errors.'; exit 1; fi

    # Remove the PCR duplicates

    echo '--------------------------------------------------'
    echo 'THE PRC DUPLICATES ARE BEING REMOVED ...'

    $DDRADSEQTOOLSDIR/pcrdupremoval.py \
        --format=$FORMAT \
        --readtype=$READTYPE \
        --readsfile1=$READSDIR/$SCEREVISIAE'-reads-'$I'-1.fastq' \
        --readsfile2=$READSDIR/$SCEREVISIAE'-reads-'$I'-2.fastq' \
        --clearfile=$READSDIR/$SCEREVISIAE'-reads-cleared-'$I \
        --dupstfile=$STATSDIR/$SCEREVISIAE'-pcrduplicates-stats-'$I'.txt' \
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
        --readsfile1=$READSDIR/$SCEREVISIAE'-reads-cleared-'$I'-1.fastq' \
        --readsfile2=$READSDIR/$SCEREVISIAE'-reads-cleared-'$I'-2.fastq' \
        --verbose=YES \
        --trace=NO
    if [ $? -ne 0 ]; then echo 'Script ended with errors.'; exit 1; fi

    ls $READSDIR/demultiplexed-*-1.fastq > $READSDIR/reads-files.txt

    while read FILE_1; do

        FILE_1_R=`echo $FILE_1 | sed 's/demultiplexed/'$SCEREVISIAE'-demultiplexed-'$I'/g'`
        FILE_2=`echo $FILE_1 | sed 's/-1.fastq/-2.fastq/g'`
        FILE_2_R=`echo $FILE_2 | sed 's/demultiplexed/'$SCEREVISIAE'-demultiplexed-'$I'/g'`

        mv $FILE_1 $FILE_1_R
        mv $FILE_2 $FILE_2_R

    done < $READSDIR/reads-files.txt

    # Calcultes a statistics of mutated reads and no mutared reads of each individual

    echo '--------------------------------------------------'
    echo 'MUTATION STATISTICS ARE BEING CALCULATED ...'

    FILE_STATS=$STATSDIR/$SCEREVISIAE'-mutations-stats-'$I'.csv'
cat >$FILE_STATS <<EOF
individual;no mutated reads;mutated reads;
EOF

    ls `echo $READSDIR/$SCEREVISIAE'-demultiplexed-'$I'-ind*-1.fastq'` > $READSDIR/reads-files.txt

    while read FILE_1; do

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

done

#-------------------------------------------------------------------------------

# End
echo '**************************************************'
exit 0

#-------------------------------------------------------------------------------
