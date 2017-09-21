#!/bin/bash

#-------------------------------------------------------------------------------

# This software has been developed by Forest Genetics and Physiology Research Group,
# Technical University of Madrid (UPM)

# Licence: GNU General Public Licence Version 3

#-------------------------------------------------------------------------------

# This script analyzes how several enzymes pairs perform a double digest of the
# genome of Saccharomyces cerevisiae, Homo sapiens and Pinus taeda

#-------------------------------------------------------------------------------
 
# Control parameters

if [ -n "$*" ]; then echo 'This script has not parameters'; exit 1; fi

#-------------------------------------------------------------------------------

# Set run environment

DDRADSEQTOOLSDIR=$TRABAJO/ddRADseqTools       # ddRADseqTools programs directory
GENOMESDIR=$TRABAJO/ddRADseqTools/genomes     # genomes file directory
FRAGSDIR=$TRABAJO/ddRADseqTools/fragments     # fragments directory
STATSDIR=$TRABAJO/ddRADseqTools/statistics    # statistics directory

if [ ! -d "$FRAGSDIR" ]; then mkdir $FRAGSDIR; else rm -f $FRAGSDIR/*; fi
if [ ! -d "$STATSDIR" ]; then mkdir $STATSDIR; else rm -f $STATSDIR/*; fi

SCEREVISIAE='Scerevisiae'                                # Saccharomyces cerevisiae
SCEREVISIAE_GENOME=GCF_000146045.2_R64_genomic.fna.gz    # Saccharomyces cerevisiae genome

HSAPIENS='Hsapiens'                                          # Homo sapiens
HSAPIENS_GENOME=GCF_000001405.29_GRCh38.p3_genomic.fna.gz    # Homo sapiens genome

PTAEDA='Ptaeda'                                  # Pinus taeda
PTAEDA_GENOME=ptaeda.v1.01.scaffolds.fasta.gz    # Pinus taeda genome

ENZYME1=(EcoRI SbfI PstI)    # codes of 1st restriction enzyme to study
ENZYME2=(MseI)               # codes of 2nd restriction enzyme to study

#-------------------------------------------------------------------------------

# Generate Saccharomyces cerevisiae genome fragments

echo '**************************************************'
echo 'SACCHAROMYCES CEREVISIAE - GENOME FRAGMENTS ARE BEING GENERATED ...'

for I in "${ENZYME1[@]}"
do
    for J in "${ENZYME2[@]}"
    do

        echo '--------------------------------------------------'
        echo "SIMULATION WITH ENZYME1=$I AND ENZYME2=$J ..."

        $DDRADSEQTOOLSDIR/rsitesearch.py \
            --genfile=$GENOMESDIR/$SCEREVISIAE_GENOME \
            --fragsfile=$FRAGSDIR/$SCEREVISIAE'-fragments-'$I'-'$J'.fasta' \
            --rsfile=$DDRADSEQTOOLSDIR/restrictionsites.txt \
            --enzyme1=$I \
            --enzyme2=$J \
            --minfragsize=101 \
            --maxfragsize=300 \
            --fragstfile=$STATSDIR/$SCEREVISIAE'-fragments-'$I'-'$J'-stats.txt' \
            --fragstinterval=25 \
            --plot=YES \
            --verbose=YES \
            --trace=NO
        if [ $? -ne 0 ]; then echo 'Script ended with errors.'; exit 1; fi

    done
done

#-------------------------------------------------------------------------------

# Generate Homo sapiens genome fragments

echo '**************************************************'
echo 'HOMO SAPIENS - FRAGMENTS ARE BEING GENERATED ...'

for I in "${ENZYME1[@]}"
do
    for J in "${ENZYME2[@]}"
    do

        echo '--------------------------------------------------'
        echo "SIMULATION WITH ENZYME1=$I AND ENZYME2=$J ..."

        $DDRADSEQTOOLSDIR/rsitesearch.py \
            --genfile=$GENOMESDIR/$HSAPIENS_GENOME \
            --fragsfile=$FRAGSDIR/$HSAPIENS'-fragments-'$I'-'$J'.fasta' \
            --rsfile=$DDRADSEQTOOLSDIR/restrictionsites.txt \
            --enzyme1=$I \
            --enzyme2=$J \
            --minfragsize=201 \
            --maxfragsize=300 \
            --fragstfile=$STATSDIR/$HSAPIENS'-fragments-'$I'-'$J'-stats.txt' \
            --fragstinterval=25 \
            --plot=YES \
            --verbose=YES \
            --trace=NO
        if [ $? -ne 0 ]; then echo 'Script ended with errors.'; exit 1; fi

    done
done

#-------------------------------------------------------------------------------

# Generate Pinus taeda genome fragments

echo '**************************************************'
echo 'PINUS TAEDA - GENOME FRAGMENTS ARE BEING GENERATED ...'

for I in "${ENZYME1[@]}"
do
    for J in "${ENZYME2[@]}"
    do

        echo '--------------------------------------------------'
        echo "SIMULATION WITH ENZYME1=$I AND ENZYME2=$J ..."

        $DDRADSEQTOOLSDIR/rsitesearch.py \
            --genfile=$GENOMESDIR/$PTAEDA_GENOME \
            --fragsfile=$FRAGSDIR/$PTAEDA'-fragments-'$I'-'$J'.fasta' \
            --rsfile=$DDRADSEQTOOLSDIR/restrictionsites.txt \
            --enzyme1=$I \
            --enzyme2=$J \
            --minfragsize=201 \
            --maxfragsize=300 \
            --fragstfile=$STATSDIR/$PTAEDA'-fragments-'$I'-'$J'-stats.txt' \
            --fragstinterval=25 \
            --plot=YES \
            --verbose=YES \
            --trace=NO
        if [ $? -ne 0 ]; then echo 'Script ended with errors.'; exit 1; fi

    done
done

#-------------------------------------------------------------------------------

# End
echo '**************************************************'
exit 0

#-------------------------------------------------------------------------------
