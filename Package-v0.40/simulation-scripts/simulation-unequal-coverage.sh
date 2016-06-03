#!/bin/bash

#-------------------------------------------------------------------------------

# This software has been developed by Forest Genetics and Physiology Research Group,
# Technical University of Madrid (UPM)

# Licence: GNU General Public Licence Version 3

#-------------------------------------------------------------------------------

# This script analyses the coverage variation among individuals across loci.

#-------------------------------------------------------------------------------
 
# Control parameters

if [ -n "$*" ]; then echo 'This script has not parameters'; exit 1; fi

#-------------------------------------------------------------------------------

# Set run environment

DDRADSEQTOOLSDIR=$TRABAJO/ddRADseqTools/Exe    # ddRADseqTools programs directory
GENOMESDIR=$TRABAJO/ddRADseqTools/Genomes      # genomes file directory
FRAGSDIR=$TRABAJO/ddRADseqTools/Fragments      # fragments directory
READSDIR=$TRABAJO/ddRADseqTools/Reads          # reads directory
STATSDIR=$TRABAJO/ddRADseqTools/Statistics     # statistics directory

if [ ! -d "$READSDIR" ]; then mkdir $READSDIR; else rm -f $READSDIR/*; fi
if [ ! -d "$STATSDIR" ]; then mkdir $STATSDIR; else rm -f $STATSDIR/*; fi

GENOME=GCF_000146045.2_R64_genomic.fna.gz    # Saccharomyces cerevisiae genome
ENZYME1=EcoRI
ENZYME2=MseI
TECHNIQUE=IND1_IND2
FORMAT=FASTQ
READTYPE=PE
INDEX1LEN=6
INDEX2LEN=6
DBRLEN=0
WEND=end71
CEND=end72
INDIVIDUALSFILE=individuals-8index1-6index2.txt

READSNUM=(300000 600000 1200000 2400000)    # readsnum values
PCRDUPPROB=(0.0)                            # pcrdupprob values

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
    --fragstinterval=25
if [ $? -ne 0 ]; then echo 'Script ended with errors.'; exit 1; fi

#-------------------------------------------------------------------------------

for I in "${READSNUM[@]}"
do
    for J in "${PCRDUPPROB[@]}"
    do

        echo '**************************************************'
        echo "READSNUM=$I AND PCRDUPPROB=$J"

        # Generate ddRADseq simulated reads

        echo '--------------------------------------------------'
        echo 'SIMULATED READS ARE BEING GENERATED ...'

        $DDRADSEQTOOLSDIR/simddradseq.py \
            --fragsfile=$FRAGSDIR/fragments-genome.fasta \
            --technique=$TECHNIQUE \
            --format=$FORMAT \
            --readsfile=$READSDIR/'reads-'$I'-'$J \
            --readtype=$READTYPE \
            --rsfile=$DDRADSEQTOOLSDIR/restrictionsites.txt \
            --enzyme1=EcoRI \
            --enzyme2=MseI \
            --endsfile=$DDRADSEQTOOLSDIR/ends.txt \
            --index1len=$INDEX1LEN \
            --index2len=$INDEX2LEN \
            --dbrlen=$DBRLEN \
            --wend=$WEND \
            --cend=$CEND \
            --individualsfile=$DDRADSEQTOOLSDIR/$INDIVIDUALSFILE \
            --locinum=3000 \
            --readsnum=$I \
            --minreadvar=0.8 \
            --maxreadvar=1.2 \
            --insertlen=100 \
            --mutprob=0.2 \
            --locusmaxmut=1 \
            --indelprob=0.1 \
            --maxindelsize=10 \
            --dropout=0.0 \
            --pcrdupprob=$J \
            --pcrdistribution=MULTINOMIAL \
            --multiparam=0.167,0.152,0.136,0.121,0.106,0.091,0.076,0.061,0.045,0.030,0.015 \
            --poissonparam=1.0 \
            --gcfactor=0.2
        if [ $? -ne 0 ]; then echo 'Script ended with errors.'; exit 1; fi

        # Remove the PCR duplicates

        echo '--------------------------------------------------'
        echo 'REMOVING PCR DUPLICATES ...'

        $DDRADSEQTOOLSDIR/pcrdupremoval.py \
            --format=$FORMAT \
            --readtype=$READTYPE \
            --readsfile1=$READSDIR/'reads-'$I'-'$J'-1.fastq' \
            --readsfile2=$READSDIR/'reads-'$I'-'$J'-2.fastq' \
            --clearfile=$READSDIR/'reads-cleared-'$I'-'$J \
            --dupstfile=$STATSDIR/'pcrduplicates-stats-'$I'-'$J'.txt'
        if [ $? -ne 0 ]; then echo 'Script ended with errors.'; exit 1; fi

    done
done

#-------------------------------------------------------------------------------

# End
echo '**************************************************'
exit 0

#-------------------------------------------------------------------------------
